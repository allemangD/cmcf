/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// CMCFlib Logic includes
#include "vtkSlicerCMCFlibLogic.h"

// MRML includes
#include <vtkMRMLModelNode.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLSequenceNode.h>

// VTK includes
#include <vtkCleanPolyData.h>
#include <vtkConnectivityFilter.h>
#include <vtkContourFilter.h>
#include <vtkCurvatures.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>

// STD includes
#include <cassert>
#include <iostream>

// IGL includes
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/principal_curvature.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_unreferenced.h>

/// Read vertex positions and face indices from polydata.
void from_polydata(vtkPolyData *pdata, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
  vtkNew<vtkTriangleFilter> triangulate;
  triangulate->SetInputData(pdata);
  triangulate->PassLinesOff();
  triangulate->PassVertsOff();
  triangulate->Update();

  pdata = triangulate->GetOutput();

  V.resize(pdata->GetNumberOfPoints(), 3);
  F.resize(pdata->GetNumberOfCells(), 3);

  for (int i = 0; i < pdata->GetNumberOfPoints(); ++i) {
    double v[3];
    pdata->GetPoint(i, v);
    V.row(i) << v[0], v[1], v[2];
  }

  for (int i = 0; i < pdata->GetNumberOfCells(); ++i) {
    auto *cell = pdata->GetCell(i);
    if (cell->GetNumberOfPoints() != 3) { throw std::runtime_error("Failed to triangulate input mesh."); }

    for (int j = 0; j < 3; ++j) { F(i, j) = static_cast<int>(cell->GetPointId(j)); }
  }
}

/// Write vertex positions into polydata.
void to_polydata(Eigen::MatrixXd const &V, vtkPolyData *pdata) {
  pdata->GetPoints()->SetNumberOfPoints(V.rows());
  for (int i = 0; i < V.rows(); ++i) {
    Eigen::Vector3d pt = V.row(i);
    pdata->GetPoints()->SetPoint(i, pt.data());
  }
  pdata->GetPoints()->Modified();
}

/// Rescale and recenter vertices. Else errors accumulate.
void renorm(Eigen::MatrixXd &V, Eigen::Matrix3Xi const &F) {
  Eigen::VectorXd double_area;
  igl::doublearea(V, F, double_area);
  double area = double_area.sum() / 2;

  Eigen::MatrixXd barycenter;
  igl::barycenter(V, F, barycenter);

  Eigen::RowVector3d centroid = Eigen::RowVector3d::Zero();
  for (int row = 0; row < barycenter.rows(); ++row) { centroid += double_area(row) * barycenter.row(row); }
  centroid /= area * 2;

  V.rowwise() -= centroid;
  V.array() /= sqrt(area);
}

/// Remove unreferenced vertices and update face indices. Else mass matrix will be singular.
void cleanup(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
  Eigen::VectorXi I, J;
  // Pass a copy of inputs to prevent aliasing.
  igl::remove_unreferenced(Eigen::MatrixXd(V), Eigen::MatrixXi(F), V, F, I, J);
}

/// Generate stages of conformalized mean curvature flow.
void vtkSlicerCMCFlibLogic::GenerateCMCFSequence(
  vtkMRMLModelNode *model,
  vtkMRMLSequenceNode *sequence,
  double rate,
  int stages
) {
  std::printf("Flowing %d stages at rate %f.\n", stages, rate);

  sequence->RemoveAllDataNodes();
  sequence->SetIndexType(vtkMRMLSequenceNode::NumericIndex);
  sequence->SetIndexName("stage");

  vtkNew<vtkMRMLModelNode> temp_model;
  temp_model->Copy(model);

  // Easier than computing normals for outputs. Let the renderer handle it or compute them downstream.
  temp_model->GetPolyData()->GetPointData()->SetNormals(nullptr);

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  from_polydata(temp_model->GetPolyData(), V, F);
  cleanup(V, F);
  renorm(V, F);

  to_polydata(V, temp_model->GetPolyData());
  sequence->SetDataNodeAtValue(temp_model, std::to_string(0));

  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;

  igl::cotmatrix(V, F, L);
  std::printf("Constructed cotmatrix with %ld nonzeros.", L.nonZeros());

  for (int stage = 1; stage <= stages; ++stage) {
    igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_BARYCENTRIC, M);
    if (stage == 1) std::printf("Constructed cotmatrix with %ld nonzeros.\n", L.nonZeros());

    solver.compute(M - rate * L);
    V = solver.solve(M * V).eval();
    renorm(V, F);

    to_polydata(V, temp_model->GetPolyData());
    sequence->SetDataNodeAtValue(temp_model, std::to_string(stage * rate));
  }
}

void vtkSlicerCMCFlibLogic::IdentifyParabolics(vtkMRMLSequenceNode *sequence, int skip) {
  // vtkNew<vtkCurvatures> H;
  // H->SetCurvatureTypeToMean();

  vtkNew<vtkCurvatures> K;
  K->SetCurvatureTypeToGaussian();
  // K->SetInputConnection(H->GetOutputPort());

  vtkNew<vtkContourFilter> contour;
  contour->SetInputConnection(K->GetOutputPort());
  contour->SetValue(0, 0);

  vtkNew<vtkCleanPolyData> clean;
  clean->SetInputConnection(contour->GetOutputPort());
  clean->SetTolerance(0.01);  // TODO parameterize; tune default.
  clean->ConvertLinesToPointsOn();
  clean->ConvertPolysToLinesOn();
  clean->ConvertStripsToPolysOn();
  clean->PointMergingOn();

  vtkNew<vtkConnectivityFilter> connect;
  connect->SetInputConnection(clean->GetOutputPort());
  connect->SetExtractionModeToAllRegions();

  for (int stage = skip; stage < sequence->GetNumberOfDataNodes(); ++stage) {
    auto model = dynamic_cast<vtkMRMLModelNode *>(sequence->GetNthDataNode(stage));

    // H->SetInputConnection(model->GetPolyDataConnection());
    K->SetInputConnection(model->GetPolyDataConnection());
    connect->Update();

    vtkPolyData *parabolic = contour->GetOutput();
    std::printf(
      "Stage %d parabolic curve has %ld points in %d regions.\n",
      stage,
      parabolic->GetNumberOfPoints(),
      connect->GetNumberOfExtractedRegions()
    );
  }
}

// region Slicer Module Logic (required boilerplate for Python wrapping)

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerCMCFlibLogic);

//----------------------------------------------------------------------------
vtkSlicerCMCFlibLogic::vtkSlicerCMCFlibLogic() {}

//----------------------------------------------------------------------------
vtkSlicerCMCFlibLogic::~vtkSlicerCMCFlibLogic() {}

//----------------------------------------------------------------------------
void vtkSlicerCMCFlibLogic::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerCMCFlibLogic::SetMRMLSceneInternal(vtkMRMLScene *newScene) {
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerCMCFlibLogic::RegisterNodes() { assert(this->GetMRMLScene() != 0); }

//---------------------------------------------------------------------------
void vtkSlicerCMCFlibLogic::UpdateFromMRMLScene() { assert(this->GetMRMLScene() != 0); }

//---------------------------------------------------------------------------
void vtkSlicerCMCFlibLogic::OnMRMLSceneNodeAdded(vtkMRMLNode *vtkNotUsed(node)) {}

//---------------------------------------------------------------------------
void vtkSlicerCMCFlibLogic::OnMRMLSceneNodeRemoved(vtkMRMLNode *vtkNotUsed(node)) {}

// endregion