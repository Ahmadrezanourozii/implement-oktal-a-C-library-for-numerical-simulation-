#include "oktal/io/VtkExport.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <vector>

namespace {
using namespace oktal;

double l2Norm(const std::vector<double> &vec) {
  double sum = 0.0;
  for (const auto val : vec) {
    sum += val * val;
  }
  return std::sqrt(sum);
}

void computeResidual(const CellGrid &grid, const std::vector<double> &u,
                     const std::vector<double> &f,
                     std::vector<double> &residual) {
  const auto &octree = grid.octree();
  residual.resize(grid.size(), 0.0);

  const std::array<Vec<std::ptrdiff_t, 3>, 6> offsets = {
      {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}}};

  for (size_t i = 0; i < grid.size(); ++i) {
    const auto cell_view = CellGrid::CellView(&grid, i);
    if (!cell_view.isValid()) {
      continue;
    }

    const auto octree_cell = cell_view.octreeCell();
    const size_t level = octree_cell.level();
    const double h = octree.geometry().dx(level);
    const double h2 = h * h;

    double neighbor_sum = 0.0;
    size_t neighbor_count = 0;
    for (const auto &offset : offsets) {
      try {
        const auto neighbor_indices = grid.neighborIndices(offset);
        if (neighbor_indices[i] != CellGrid::NO_NEIGHBOR) {
          const size_t neighbor_idx = neighbor_indices[i];
          neighbor_sum += u[neighbor_idx];
          ++neighbor_count;
        }
      } catch (const std::out_of_range &) { // NOLINT(bugprone-empty-catch)
      }
    }

    const double laplacian =
        (static_cast<double>(neighbor_count) * u[i] - neighbor_sum) / h2;

    residual[i] = f[i] + laplacian;
  }
}

void jacobiIteration(const CellGrid &grid, const std::vector<double> &u_old,
                     const std::vector<double> &f, std::vector<double> &u_new) {
  const auto &octree = grid.octree();

  u_new.resize(grid.size(), 0.0);

  const std::array<Vec<std::ptrdiff_t, 3>, 6> offsets = {
      {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}}};

  for (size_t i = 0; i < grid.size(); ++i) {
    const auto cell_view = CellGrid::CellView(&grid, i);
    if (!cell_view.isValid()) {
      u_new[i] = u_old[i];
      continue;
    }

    const auto octree_cell = cell_view.octreeCell();
    const size_t level = octree_cell.level();
    const double h = octree.geometry().dx(level);
    const double h2 = h * h;

    double neighbor_sum = 0.0;
    size_t neighbor_count = 0;
    for (const auto &offset : offsets) {
      try {
        const auto neighbor_indices = grid.neighborIndices(offset);
        if (neighbor_indices[i] != CellGrid::NO_NEIGHBOR) {
          const size_t neighbor_idx = neighbor_indices[i];
          neighbor_sum += u_old[neighbor_idx];
          ++neighbor_count;
        }
      } catch (const std::out_of_range &) { // NOLINT(bugprone-empty-catch)
      }
    }

    if (neighbor_count > 0) {
      u_new[i] =
          (f[i] * h2 + neighbor_sum) / static_cast<double>(neighbor_count);
    } else {
      u_new[i] = u_old[i];
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr
        << "Usage: "
        << argv[0] // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        << " <refinementLevel> <max-iterations> <epsilon> <output-file>"
        << '\n';
    return EXIT_FAILURE;
  }

  const size_t refinement_level = std::stoul(
      argv[1]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const size_t max_iterations = std::stoul(
      argv[2]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const double epsilon = std::stod(
      argv[3]); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const std::filesystem::path output_file =
      argv[4]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

  const auto octree = CellOctree::createUniformGrid(refinement_level);

  const std::array<Vec<std::ptrdiff_t, 3>, 6> offsets = {
      {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}}};

  const auto grid = CellGrid::create(octree)
                        .levels({refinement_level})
                        .neighborhood({offsets[0], offsets[1], offsets[2],
                                       offsets[3], offsets[4], offsets[5]})
                        .build();

  std::vector<double> u(grid.size(), 0.0);
  std::vector<double> f(grid.size(), 0.0);
  std::vector<double> residual(grid.size(), 0.0);

  std::ranges::fill(f, 1.0);

  size_t iteration = 0;
  double residual_norm = std::numeric_limits<double>::max();

  while (iteration < max_iterations && residual_norm > epsilon) {
    computeResidual(grid, u, f, residual);
    residual_norm = l2Norm(residual);

    std::vector<double> u_new;
    jacobiIteration(grid, u, f, u_new);
    u = std::move(u_new);

    ++iteration;
  }

  std::cout << "Iterations: " << iteration << '\n';
  std::cout << "Final residual L2-norm: " << residual_norm << '\n';

  io::vtk::exportCellGrid(grid, output_file)
      .writeGridVector<double>("u", u)
      .writeGridVector<double>("f", f)
      .writeGridVector<double>("residual", residual);

  return EXIT_SUCCESS;
}
