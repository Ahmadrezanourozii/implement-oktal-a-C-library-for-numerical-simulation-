#include "oktal/io/VtkExport.hpp"
#include "oktal/data/GridVector.hpp"
#include "oktal/octree/CellGrid.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
#include <experimental/mdspan>

#include <iostream>
#include <span>
#include <vector>

namespace oktal::io::vtk {

namespace {

template <typename T>
void write_dataset(HighFive::Group &group, const std::string &name,
                   const std::vector<T> &data) {
  if (!data.empty()) {
    group.createDataSet(name, data);
  } else {
    group.createDataSet<T>(name, HighFive::DataSpace({0}));
  }
}

void set_bit(std::vector<uint8_t> &buffer, size_t bitIndex) {
  if (buffer.empty()) {
    return;
  }
  const auto byteIndex = bitIndex / 8;
  const auto bitOffset = 7 - static_cast<int>(bitIndex % 8);
  buffer[byteIndex] |= static_cast<uint8_t>(1u << bitOffset);
}
void writeVtkAttributes(HighFive::Group &vtk_group) {
  const std::string type = "HyperTreeGrid";
  auto attr_type = vtk_group.createAttribute<std::string>(
      "Type", HighFive::DataSpace::From(type));
  attr_type.write(type);

  const std::array<int, 2> version = {1, 0};
  auto attr_version = vtk_group.createAttribute<int>(
      "Version", HighFive::DataSpace::From(version));
  attr_version.write(version);
}

void writeCoordinates(HighFive::Group &vtk_group, const Vec3D &origin,
                      double size) {
  const std::vector<double> x_coords = {origin[0], origin[0] + size};
  const std::vector<double> y_coords = {origin[1], origin[1] + size};
  const std::vector<double> z_coords = {origin[2], origin[2] + size};

  write_dataset(vtk_group, "XCoordinates", x_coords);
  write_dataset(vtk_group, "YCoordinates", y_coords);
  write_dataset(vtk_group, "ZCoordinates", z_coords);
}

void buildDescriptors(const CellOctree &ot, std::vector<uint8_t> &descriptors,
                      size_t &descriptorBitIndex) {
  const size_t levelCount = ot.numberOfLevels();
  for (size_t level = 0; level + 1 < levelCount; ++level) {
    for (const auto &node : ot.nodesStream(level)) {
      if (node.isRefined()) {
        set_bit(descriptors, descriptorBitIndex);
      }
      ++descriptorBitIndex;
    }
  }
}

void buildMask(const CellOctree &ot, std::vector<uint8_t> &mask,
               size_t &maskBitIndex) {
  const size_t levelCount = ot.numberOfLevels();
  for (size_t level = 0; level < levelCount; ++level) {
    for (const auto &node : ot.nodesStream(level)) {
      if (node.isPhantom()) {
        set_bit(mask, maskBitIndex);
      }
      ++maskBitIndex;
    }
  }
}

void fixRootMask(const CellOctree &ot, std::vector<uint8_t> &mask,
                 size_t totalNodes) {
  if (!mask.empty() && totalNodes > 0) {
    const bool rootIsPhantom = ot.nodesStream(0)[0].isPhantom();
    mask[0] = rootIsPhantom ? 1u : 0u;
  }
}

void buildLevels(const CellOctree &ot, std::vector<size_t> &levels) {
  const size_t levelCount = ot.numberOfLevels();
  for (size_t lvl = 0; lvl < levelCount; ++lvl) {
    const auto nodes = ot.nodesStream(lvl);
    for (size_t i = 0; i < nodes.size(); ++i) {
      levels.push_back(lvl);
    }
  }
}

} // namespace

void exportOctree(const CellOctree &ot, const std::filesystem::path &filename) {
  HighFive::File file(filename.string(), HighFive::File::Overwrite);
  auto vtk_group = file.createGroup("VTKHDF");

  writeVtkAttributes(vtk_group);

  const auto &geo = ot.geometry();
  const Vec3D origin = geo.origin();
  const double size = geo.sidelength();

  writeCoordinates(vtk_group, origin, size);

  const size_t levelCount = ot.numberOfLevels();

  size_t descriptorNodes = 0;
  for (size_t level = 0; level + 1 < levelCount; ++level) {
    descriptorNodes += ot.nodesStream(level).size();
  }

  std::vector<uint8_t> descriptors((descriptorNodes + 7) / 8, 0);
  size_t descriptorBitIndex = 0;
  buildDescriptors(ot, descriptors, descriptorBitIndex);

  const size_t totalNodes = ot.numberOfNodes();
  std::vector<uint8_t> mask((totalNodes + 7) / 8, 0);
  size_t maskBitIndex = 0;
  buildMask(ot, mask, maskBitIndex);
  fixRootMask(ot, mask, totalNodes);

  write_dataset(vtk_group, "Descriptors", descriptors);
  write_dataset(vtk_group, "Mask", mask);

  std::vector<size_t> levels;
  buildLevels(ot, levels);

  auto cell_data_group = vtk_group.createGroup("CellData");
  write_dataset(cell_data_group, "level", levels);
}

CellGridExporter::CellGridExporter(const CellGrid &grid,
                                   const std::filesystem::path &filename)
    : grid_(grid), filename_(filename),
      file_(filename.string(), HighFive::File::Overwrite),
      vtk_group_(file_.createGroup("VTKHDF")),
      cell_data_group_(vtk_group_.createGroup("CellData")) {
  writeVtkAttributes(vtk_group_);

  const auto &geo = grid_.octree().geometry();
  const Vec3D origin = geo.origin();
  const double size = geo.sidelength();
  writeCoordinates(vtk_group_, origin, size);
  const auto &octree = grid_.octree();
  const size_t levelCount = octree.numberOfLevels();

  size_t descriptorNodes = 0;
  for (size_t level = 0; level + 1 < levelCount; ++level) {
    descriptorNodes += octree.nodesStream(level).size();
  }

  std::vector<uint8_t> descriptors((descriptorNodes + 7) / 8, 0);
  size_t descriptorBitIndex = 0;
  buildDescriptors(octree, descriptors, descriptorBitIndex);

  const size_t totalNodes = octree.numberOfNodes();
  std::vector<uint8_t> mask((totalNodes + 7) / 8, 0);
  size_t maskBitIndex = 0;
  buildMask(octree, mask, maskBitIndex);
  fixRootMask(octree, mask, totalNodes);

  write_dataset(vtk_group_, "Descriptors", descriptors);
  write_dataset(vtk_group_, "Mask", mask);

  std::vector<size_t> levels;
  buildLevels(octree, levels);
  write_dataset(cell_data_group_, "level", levels);
}

template <typename T>
CellGridExporter &CellGridExporter::writeGridVector(const std::string &name,
                                                    std::span<const T> data) {
  const auto &octree = grid_.octree();
  const size_t totalNodes = octree.numberOfNodes();

  std::vector<T> vtk_data(totalNodes, T(0));

  for (size_t i = 0; i < data.size() && i < grid_.size(); ++i) {
    const auto cell_view = CellGrid::CellView(&grid_, i);
    if (cell_view.isValid()) {
      const auto octree_cell = cell_view.octreeCell();
      const size_t stream_idx = octree_cell.streamIndex();
      if (stream_idx < totalNodes) {
        vtk_data[stream_idx] = data[i];
      }
    }
  }

  write_dataset(cell_data_group_, name, vtk_data);
  return *this;
}
// Added for task 6
template <typename T, size_t Q>
CellGridExporter &
CellGridExporter::writeGridVector(const std::string &id,
                                  GridVectorView<const T, Q> vector) {
  const auto &octree = grid_.octree();
  const size_t totalNodes = octree.numberOfNodes();
  std::vector<T> vtk_data(totalNodes * Q, T(0));

  // VTK/HDF5 and MdCellDataView require Array-of-Structures (layout_right)
  namespace stdex = std::experimental;
  using AoSView =
      stdex::mdspan<T, stdex::extents<size_t, stdex::dynamic_extent, Q>,
                    stdex::layout_right>;
  const AoSView vtk_view(vtk_data.data(), totalNodes);

  // Transpose from SoA (GridVector) to AoS (VTK)
  for (size_t i = 0; i < grid_.size(); ++i) {
    const auto cell_view = CellGrid::CellView(&grid_, i);
    if (cell_view.isValid()) {
      const size_t stream_idx = cell_view.octreeCell().streamIndex();

      if constexpr (Q == 1) {
        vtk_view[stream_idx, 0] = vector[i];
      } else {
        for (size_t q = 0; q < Q; ++q) {
          vtk_view[stream_idx, q] = vector[i, q];
        }
      }
    }
  }

  if constexpr (Q == 1) {
    // For scalar fields, use 1D dataset
    auto dataset = cell_data_group_.createDataSet<T>(
        id, HighFive::DataSpace({totalNodes}));
    dataset.write(vtk_data);
  } else {
    // For vector fields, use 2D dataset and write_raw for proper
    // multi-dimensional writing
    auto dataset = cell_data_group_.createDataSet<T>(
        id, HighFive::DataSpace({totalNodes, Q}));
    dataset.write_raw(vtk_data.data());
  }
  return *this;
}

// Explicit template instantiations
template CellGridExporter &
CellGridExporter::writeGridVector<double>(const std::string &,
                                          std::span<const double>);
template CellGridExporter &
CellGridExporter::writeGridVector<float>(const std::string &,
                                         std::span<const float>);
template CellGridExporter &
CellGridExporter::writeGridVector<int>(const std::string &,
                                       std::span<const int>);
// Added for task 6
//  Explicit instantiations for GridVectorView
template CellGridExporter &
CellGridExporter::writeGridVector<double, 1>(const std::string &,
                                             GridVectorView<const double, 1>);
template CellGridExporter &
CellGridExporter::writeGridVector<double, 3>(const std::string &,
                                             GridVectorView<const double, 3>);

CellGridExporter exportCellGrid(const CellGrid &grid,
                                const std::filesystem::path &filename) {
  return {grid, filename};
}

} // namespace oktal::io::vtk
