#include "oktal/octree/CellGrid.hpp"

#include <algorithm>
#include <optional>
#include <ranges>

namespace oktal {

std::optional<Vec<std::size_t, 3>>
NoPeriodicity::operator()(const Vec<std::size_t, 3> & /*coords*/) const {
  return std::nullopt;
}

std::optional<Vec<std::size_t, 3>>
Torus::operator()(const Vec<std::size_t, 3> &coords) const {
  Vec<std::size_t, 3> result = coords;

  for (std::size_t i = 0; i < 3; ++i) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
    if (periodic_directions_[i]) {
      result[i] = coords[i] % grid_resolution_;
    } else {
      if (coords[i] >= grid_resolution_) {
        return std::nullopt;
      }
      result[i] = coords[i];
    }
  }

  return result;
}

// CellGrid::CellView implementations

Vec3D CellGrid::CellView::center() const {
  if (!isValid()) {
    throw std::runtime_error("Invalid cell view");
  }
  const auto morton = grid_->morton_indices_[enumeration_index_];
  return grid_->octree_->geometry().cellCenter(morton);
}

Box3D CellGrid::CellView::boundingBox() const {
  if (!isValid()) {
    throw std::runtime_error("Invalid cell view");
  }
  const auto morton = grid_->morton_indices_[enumeration_index_];
  return grid_->octree_->geometry().cellBoundingBox(morton);
}

CellOctree::CellView CellGrid::CellView::octreeCell() const {
  if (!isValid()) {
    throw std::runtime_error("Invalid cell view");
  }
  const auto morton = grid_->morton_indices_[enumeration_index_];
  const auto cell_opt = grid_->octree_->getCell(morton);
  if (!cell_opt.has_value()) {
    throw std::runtime_error("Cell not found in octree");
  }
  return cell_opt.value();
}

CellGrid::CellView
CellGrid::CellView::neighbor(const Vec<std::ptrdiff_t, 3> &offset) const {
  if (!isValid()) {
    return {nullptr, NOT_ENUMERATED};
  }

  try {
    const auto neighbor_indices = grid_->neighborIndices(offset);
    const size_type neighbor_enum_idx = neighbor_indices[enumeration_index_];
    if (neighbor_enum_idx == NO_NEIGHBOR) {
      return {nullptr, NOT_ENUMERATED};
    }
    return {grid_, neighbor_enum_idx};
  } catch (const std::out_of_range &) {
    return {nullptr, NOT_ENUMERATED};
  }
}

MortonIndex CellGrid::CellView::mortonIndex() const {
  if (!isValid()) {
    throw std::runtime_error("Invalid cell view");
  }
  return grid_->morton_indices_[enumeration_index_];
}

CellGrid::size_type CellGrid::CellView::level() const {
  if (!isValid()) {
    throw std::runtime_error("Invalid cell view");
  }
  return mortonIndex().level();
}

// CellGrid implementations

CellGridBuilder CellGrid::create(std::shared_ptr<const CellOctree> octree) {
  return CellGridBuilder(std::move(octree));
}

CellGrid::size_type
CellGrid::getEnumerationIndex(size_type streamIndex) const noexcept {
  const auto it = stream_to_enumeration_.find(streamIndex);
  if (it == stream_to_enumeration_.end()) {
    return NOT_ENUMERATED;
  }
  return it->second;
}

CellGrid::size_type
CellGrid::getEnumerationIndex(const CellOctree::CellView &cell) const noexcept {
  return getEnumerationIndex(cell.streamIndex());
}

std::span<const CellGrid::size_type>
CellGrid::neighborIndices(const Vec<std::ptrdiff_t, 3> &offset) const {
  const auto it = adjacency_lists_.find(offset);
  if (it == adjacency_lists_.end()) {
    throw std::out_of_range("No adjacency list for given offset");
  }
  return it->second;
}

// CellGridBuilder implementation

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
CellGrid CellGridBuilder::build() const {
  using size_type = CellGrid::size_type;

  std::vector<MortonIndex> morton_indices;
  std::unordered_map<size_type, size_type> stream_to_enumeration;

  // Determine which levels to enumerate
  std::vector<size_type> levels_to_enumerate;
  if (levels_.empty()) {
    // Enumerate all levels
    for (size_type level = 0; level < octree_->numberOfLevels(); ++level) {
      levels_to_enumerate.push_back(level);
    }
  } else {
    levels_to_enumerate = levels_;
    // Sort levels to ensure top-to-bottom ordering
    std::ranges::sort(levels_to_enumerate);
  }

  // Enumerate cells in Z-order (pre-order depth-first)
  size_type enumeration_index = 0;
  for (const size_type level : levels_to_enumerate) {
    // Use horizontal range for each level to get cells in Z-order
    const auto range = octree_->horizontalRange(level);
    for (const auto &cell : range) {
      // Only enumerate non-phantom cells
      if (!octree_->nodesStream()[cell.streamIndex()].isPhantom()) {
        morton_indices.push_back(cell.mortonIndex());
        stream_to_enumeration[cell.streamIndex()] = enumeration_index;
        ++enumeration_index;
      }
    }
  }

  // Build adjacency lists
  std::map<Vec<std::ptrdiff_t, 3>, std::vector<size_type>> adjacency_lists;

  if (!neighborhood_.empty()) {
    const size_type grid_size = morton_indices.size();
    size_type max_level = 0;
    if (levels_to_enumerate.empty()) {
      max_level =
          (octree_->numberOfLevels() > 0) ? (octree_->numberOfLevels() - 1) : 0;
    } else {
      max_level = *std::ranges::max_element(levels_to_enumerate);
    }
    const std::size_t max_grid_resolution = static_cast<std::size_t>(1)
                                            << max_level;

    // Set grid resolution for Torus mapper if needed
    if (auto *torus = dynamic_cast<Torus *>(periodicity_mapper_.get())) {
      torus->setGridResolution(max_grid_resolution);
    }

    // Initialize adjacency lists
    for (const auto &offset : neighborhood_) {
      adjacency_lists[offset].resize(grid_size, CellGrid::NO_NEIGHBOR);
    }

    // For each cell, find neighbors
    for (size_type enum_idx = 0; enum_idx < grid_size; ++enum_idx) {
      const MortonIndex morton = morton_indices[enum_idx];
      const size_type level = morton.level();
      Vec<std::size_t, 3> coords = morton.gridCoordinates();
      // Grid resolution for this cell's level
      const std::size_t grid_resolution = static_cast<std::size_t>(1) << level;

      if (auto *torus = dynamic_cast<Torus *>(periodicity_mapper_.get())) {
        torus->setGridResolution(grid_resolution);
      }

      for (const auto &offset : neighborhood_) {
        Vec<std::ptrdiff_t, 3> neighbor_coords_ptrdiff;
        neighbor_coords_ptrdiff[0] =
            static_cast<std::ptrdiff_t>(coords[0]) + offset[0];
        neighbor_coords_ptrdiff[1] =
            static_cast<std::ptrdiff_t>(coords[1]) + offset[1];
        neighbor_coords_ptrdiff[2] =
            static_cast<std::ptrdiff_t>(coords[2]) + offset[2];

        Vec<std::size_t, 3> neighbor_coords;
        bool valid = true;

        for (std::size_t i = 0; i < 3; ++i) {
          if (neighbor_coords_ptrdiff[i] < 0 ||
              neighbor_coords_ptrdiff[i] >=
                  static_cast<std::ptrdiff_t>(grid_resolution)) {
            if (periodicity_mapper_) {
              Vec<std::size_t, 3> wrapped_coords;
              wrapped_coords[0] =
                  static_cast<std::size_t>(neighbor_coords_ptrdiff[0]);
              wrapped_coords[1] =
                  static_cast<std::size_t>(neighbor_coords_ptrdiff[1]);
              wrapped_coords[2] =
                  static_cast<std::size_t>(neighbor_coords_ptrdiff[2]);

              const auto mapped = (*periodicity_mapper_)(wrapped_coords);
              if (mapped.has_value()) {
                neighbor_coords = mapped.value();
              } else {
                valid = false;
                break;
              }
            } else {
              valid = false;
              break;
            }
          } else {
            neighbor_coords[i] =
                static_cast<std::size_t>(neighbor_coords_ptrdiff[i]);
          }
        }

        if (!valid) {
          continue;
        }

        for (std::size_t i = 0; i < 3; ++i) {
          if (neighbor_coords[i] >= grid_resolution) {
            valid = false;
            break;
          }
        }

        if (!valid) {
          continue;
        }

        const MortonIndex neighbor_morton =
            MortonIndex::fromGridCoordinates(level, neighbor_coords);

        const auto neighbor_cell_opt = octree_->getCell(neighbor_morton);
        if (neighbor_cell_opt.has_value()) {
          const size_type neighbor_stream_idx =
              neighbor_cell_opt.value().streamIndex();
          const auto enum_it = stream_to_enumeration.find(neighbor_stream_idx);
          if (enum_it != stream_to_enumeration.end()) {
            if (neighbor_cell_opt.value().level() == level) {
              adjacency_lists[offset][enum_idx] = enum_it->second;
            }
          }
        }
      }
    }
  }

  return {octree_, std::move(morton_indices), std::move(stream_to_enumeration),
          std::move(adjacency_lists)};
}

} // namespace oktal
