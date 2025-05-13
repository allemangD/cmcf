#include <filesystem>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <imgui.h>

int main(int argc, char *argv[]) {
  if (argc <= 1) {
    std::cerr << "usage: [FILE]" << std::endl;
    return 1;
  }

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  auto path = std::filesystem::absolute(argv[1]);
  std::cout << "reading: " << path << std::endl;
  igl::read_triangle_mesh(path, V, F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);

  igl::opengl::glfw::imgui::ImGuiPlugin gui;
  viewer.plugins.push_back(&gui);

  Eigen::MatrixXd V_0 = V;
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  float rate = 0.001;
  int iter = 0;

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;  // about 90ms
  // Of all the other direct solvers and iterative solvers, and of different tolerances and float precision, this is the fastest.

  auto const &draw_menu = [&]() -> void {
    ImGui::Text("Rate: %.5f", rate);
    ImGui::Text("Iter: %d", iter);

    ImGui::SliderFloat("Rate", &rate, 0.001, 0.010);
  };

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  gui.widgets.push_back(&menu);
  menu.callback_draw_viewer_menu = draw_menu;

  auto const &key_down =
    [&solver, &L, &rate, &V_0, &iter](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod) -> bool {
    auto const &V = viewer.data().V;
    auto const &F = viewer.data().F;

    switch (key) {
      case 'R': {
        iter = 0;
        viewer.data().set_vertices(V_0);
        break;
      }
      case 'S': {
        Eigen::MatrixXd U = V;

        for (int k = 0; k < 10; ++k) {
          iter += 1;
          auto start_ = std::chrono::high_resolution_clock::now();

          Eigen::SparseMatrix<double> M;
          igl::massmatrix(U, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
          auto const &S = (M - rate * L);

          solver.compute(S);
          U = solver.solve(M * U).eval();

          Eigen::VectorXd double_area;
          igl::doublearea(U, F, double_area);
          double area = double_area.sum() / 2;

          Eigen::MatrixXd barycenter;
          igl::barycenter(U, F, barycenter);

          Eigen::RowVector3d centroid = Eigen::RowVector3d::Zero();
          for (int row = 0; row < barycenter.rows(); ++row) {
            centroid += double_area(row) * barycenter.row(row);
          }
          centroid /= area * 2;

          U.rowwise() -= centroid;
          U.array() /= sqrt(area);

          auto end_ = std::chrono::high_resolution_clock::now();
          std::cout << "took: "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_).count() << "ms"
                    << std::endl;
        }

        viewer.data().set_vertices(U);
        break;
      }
      case ' ': {
        iter += 1;
        auto start_ = std::chrono::high_resolution_clock::now();

        Eigen::SparseMatrix<double> M;
        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
        auto const &S = (M - rate * L);

        solver.compute(S);
        auto U = solver.solve(M * V).eval();

        Eigen::VectorXd double_area;
        igl::doublearea(U, F, double_area);
        double area = double_area.sum() / 2;

        Eigen::MatrixXd barycenter;
        igl::barycenter(U, F, barycenter);

        Eigen::RowVector3d centroid = Eigen::RowVector3d::Zero();
        for (int row = 0; row < barycenter.rows(); ++row) {
          centroid += double_area(row) * barycenter.row(row);
        }
        centroid /= area * 2;

        U.rowwise() -= centroid;
        U.array() /= sqrt(area);

        auto end_ = std::chrono::high_resolution_clock::now();
        std::cout << "took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_).count()
                  << "ms" << std::endl;

        viewer.data().set_vertices(U);
        break;
      }
      default:
        return false;
    }

    viewer.data().compute_normals();
    viewer.core().align_camera_center(V, F);
    return true;
  };

  viewer.callback_key_down = key_down;

  viewer.launch();
}
