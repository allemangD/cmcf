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
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkTriangleFilter.h>
#include <vtkVector.h>

// STD includes
#include <cassert>

// IGL includes
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <igl/cotmatrix.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/principal_curvature.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_unreferenced.h>

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

void from_polydata(vtkPolyData *pdata, Eigen::MatrixX3d &V, Eigen::MatrixX3i &F) {
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

void to_polydata(Eigen::MatrixX3d const &V, vtkPolyData *pdata) {
  pdata->GetPoints()->SetNumberOfPoints(V.rows());
  for (int i = 0; i < V.rows(); ++i) {
    Eigen::Vector3d pt = V.row(i);
    pdata->GetPoints()->SetPoint(i, pt.data());
  }
  pdata->GetPoints()->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerCMCFlibLogic::GenerateSequence(vtkMRMLModelNode *model, vtkMRMLSequenceNode *sequence) {
  sequence->RemoveAllDataNodes();

  auto rate = 1.3;
  // auto stages = 100;

  // TODO verify that the model is not empty.
  Eigen::MatrixX3d V;
  Eigen::MatrixX3i F;
  from_polydata(model->GetPolyData(), V, F);

  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> M;
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;

  igl::cotmatrix(V, F, L);
  std::cout << "cotmatrix: " << L.nonZeros() << " nonzeros." << std::endl;

  igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_BARYCENTRIC, M);
  std::cout << "massmatrix: " << M.nonZeros() << " nonzeros." << std::endl;

  solver.compute(M - rate * L);

  std::cout << "Solving..." << std::endl;

  V = solver.solve(M * V).eval();

  to_polydata(V, model->GetPolyData());

  std::cout << "Computed one update!" << std::endl;
}
