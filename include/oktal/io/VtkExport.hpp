#pragma once

#include "oktal/data/GridVector.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"
#include <filesystem>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <span>
#include <string>
#include <vector>

namespace oktal::io::vtk {

/**
 * @brief Exports the CellOctree structure to a VTK HDF5 file format.
 * * @param ot The CellOctree to export.
 * @param filename The path to the output .vtkhdf file.
 */
void exportOctree(const CellOctree &ot, const std::filesystem::path &filename);

/// Exporter for CellGrid with grid vectors
class CellGridExporter {
private:
  const CellGrid
      &grid_; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
  std::filesystem::path filename_;
  HighFive::File file_;
  HighFive::Group vtk_group_;
  HighFive::Group cell_data_group_;

public:
  CellGridExporter(const CellGrid &grid, const std::filesystem::path &filename);

  template <typename T>
  CellGridExporter &writeGridVector(const std::string &name,
                                    std::span<const T> data);
  // Added this for Task 6
  template <typename T, size_t Q>
  CellGridExporter &writeGridVector(const std::string &id,
                                    GridVectorView<const T, Q> vector);
};

/**
 * @brief Exports the CellGrid structure to a VTK HDF5 file format.
 * @param grid The CellGrid to export.
 * @param filename The path to the output .vtkhdf file.
 * @return CellGridExporter object for chaining writeGridVector calls.
 */
CellGridExporter exportCellGrid(const CellGrid &grid,
                                const std::filesystem::path &filename);

} // namespace oktal::io::vtk