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
#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkConnectivityFilter.h>
#include <vtkContourFilter.h>
#include <vtkCurvatures.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkStaticPointLocator.h>
#include <vtkTriangleFilter.h>

// STD includes
#include <cassert>
#include <cstdio>

// IGL includes
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/gaussian_curvature.h>
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

vtkSmartPointer<vtkPolyData> vtkSlicerCMCFlibLogic::IdentifyParabolics(
  vtkMRMLSequenceNode *sequence,
  int skip,
  double tolerance
) {
  vtkSmartPointer<vtkConnectivityFilter> prev_connect;
  vtkSmartPointer<vtkPointSet> prev_points;

  vtkNew<vtkPolyData> output_curves;
  output_curves->SetPoints(vtkPoints::New());
  output_curves->Allocate();

  vtkNew<vtkIntArray> transition_ids;
  transition_ids->SetName("Transitions");
  int current_transition_id = 0;

  for (int stage = skip; stage < sequence->GetNumberOfDataNodes(); ++stage) {
    auto model = dynamic_cast<vtkMRMLModelNode *>(sequence->GetNthDataNode(stage));

    vtkNew<vtkCurvatures> K;
    K->SetInputConnection(model->GetPolyDataConnection());
    K->SetCurvatureTypeToGaussian();

    vtkNew<vtkContourFilter> contour;
    contour->SetInputConnection(K->GetOutputPort());
    contour->SetValue(0, 0);

    vtkNew<vtkConnectivityFilter> connect;
    connect->SetInputConnection(contour->GetOutputPort());
    connect->ColorRegionsOn();
    connect->SetExtractionModeToAllRegions();

    connect->Update();
    auto *points = dynamic_cast<vtkPolyData *>(connect->GetOutput());

    // todo Keep only VTK_LINES cells.

    vtkIdType start = output_curves->GetNumberOfPoints();
    for (vtkIdType idx = 0; idx < points->GetNumberOfPoints(); ++idx) {
      output_curves.GetPointer()->GetPoints()->InsertNextPoint(points->GetPoint(idx));
    }

    std::printf(
      "Stage %d parabolic curve has %lld points in %d regions.\n",
      stage,
      points->GetNumberOfPoints(),
      connect->GetNumberOfExtractedRegions()
    );
    if (prev_points && prev_connect) {
      vtkIdType prev_start = start - prev_points->GetNumberOfPoints();

      std::printf(
        "(prev had %lld in %d regions.)\n",
        prev_points->GetNumberOfPoints(),
        prev_connect->GetNumberOfExtractedRegions()
      );

      vtkNew<vtkStaticPointLocator> locator;
      locator->SetDataSet(prev_points);
      locator->BuildLocator();

      vtkDataArray *regions = points->GetPointData()->GetArray("RegionId");
      vtkDataArray *prev_regions = prev_points->GetPointData()->GetArray("RegionId");

      double scan_region = 0.0;
      double scan_coord[3];
      std::set<double> scan_src_regions;

      auto const log_scan = [&] {
        auto size = scan_src_regions.size();

        if (size == 0) {
          return;  // Lips event. Ignore.
        }

        if (size == 1) {
          return;  // Topology unchanged. Quiet.
        }

        // if (size != 2) { return; }  // Only show 2-curve mergers.

        std::printf("Merge.");

        {
          vtkDataArray *_regions = points->GetCellData()->GetArray("RegionId");

          for (vtkIdType idx = 0; idx < points->GetNumberOfCells(); ++idx) {
            if (_regions->GetComponent(idx, 0) == scan_region) {
              vtkNew<vtkIdList> pts;
              points->GetCellPoints(idx, pts);
              for (int i = 0; i < pts->GetNumberOfIds(); ++i) { pts->SetId(i, pts->GetId(i) + start); }
              output_curves->InsertNextCell(points->GetCellType(idx), pts);
              transition_ids->InsertNextValue(current_transition_id);
            }
          }
        }

        {
          vtkDataArray *_regions = prev_points->GetCellData()->GetArray("RegionId");

          for (double region: scan_src_regions) {
            for (vtkIdType idx = 0; idx < prev_points->GetNumberOfPoints(); ++idx) {
              if (_regions->GetComponent(idx, 0) == region) {
                vtkNew<vtkIdList> pts;
                prev_points->GetCellPoints(idx, pts);
                for (int i = 0; i < pts->GetNumberOfIds(); ++i) { pts->SetId(i, pts->GetId(i) + prev_start); }
                output_curves->InsertNextCell(prev_points->GetCellType(idx), pts);
                transition_ids->InsertNextValue(current_transition_id);
              }
            }
          }
        }

        current_transition_id++;

        std::printf(" %0.f <-", scan_region);
        for (auto const region: scan_src_regions) { std::printf(" %0.f", region); }
        std::printf("\n");
      };

      for (vtkIdType idx = 0; idx < points->GetNumberOfPoints(); ++idx) {
        auto const dst_region = regions->GetComponent(idx, 0);

        points->GetPoint(idx, scan_coord);

        double _dist2;
        vtkIdType const src_idx = locator->FindClosestPointWithinRadius(tolerance, scan_coord, _dist2);

        if (src_idx >= 0) {
          auto const src_region = prev_regions->GetComponent(src_idx, 0);

          if (scan_region != dst_region) {
            // the scan has entered a new dst_region.
            log_scan();
            scan_src_regions.clear();
            scan_region = dst_region;
          }
          scan_src_regions.insert(src_region);
        }
      }
      log_scan();
    }

    prev_connect = connect;
    prev_points = points;
  }

  output_curves->GetCellData()->AddArray(transition_ids);
  return output_curves;
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