#ifndef OKTAL_OCTREE_CELLGRID_HPP
#define OKTAL_OCTREE_CELLGRID_HPP

#include "oktal/geometry/Box.hpp"
#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/CellOctree.hpp"
#include "oktal/octree/MortonIndex.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace oktal {

// Forward declaration
class CellGridBuilder;

/// Periodicity mapper interface
class PeriodicityMapper {
public:
  PeriodicityMapper() = default;
  PeriodicityMapper(const PeriodicityMapper &) = default;
  PeriodicityMapper &operator=(const PeriodicityMapper &) = default;
  PeriodicityMapper(PeriodicityMapper &&) noexcept = default;
  PeriodicityMapper &operator=(PeriodicityMapper &&) noexcept = default;
  virtual ~PeriodicityMapper() = default;
  [[nodiscard]] virtual std::optional<Vec<std::size_t, 3>>
  operator()(const Vec<std::size_t, 3> &coords) const = 0;
};

/// No periodicity mapper
class NoPeriodicity : public PeriodicityMapper {
public:
  NoPeriodicity() = default;
  NoPeriodicity(const NoPeriodicity &) = default;
  NoPeriodicity &operator=(const NoPeriodicity &) = default;
  NoPeriodicity(NoPeriodicity &&) noexcept = default;
  NoPeriodicity &operator=(NoPeriodicity &&) noexcept = default;
  ~NoPeriodicity() override = default;

  [[nodiscard]] std::optional<Vec<std::size_t, 3>>
  operator()(const Vec<std::size_t, 3> &coords) const override;
};

/// Torus periodicity mapper
class Torus : public PeriodicityMapper {
private:
  std::array<bool, 3> periodic_directions_;

public:
  explicit Torus(const std::array<bool, 3> &periodicDirections)
      : periodic_directions_(periodicDirections) {}

  [[nodiscard]] std::optional<Vec<std::size_t, 3>>
  operator()(const Vec<std::size_t, 3> &coords) const override;

  void setGridResolution(std::size_t resolution) {
    grid_resolution_ = resolution;
  }

private:
  mutable std::size_t grid_resolution_ = 0;
};

/// Cell grid data structure for numerical PDE solvers
class CellGrid {
public:
  using size_type = std::size_t;
  static constexpr size_type NOT_ENUMERATED =
      std::numeric_limits<size_type>::max();
  static constexpr size_type NO_NEIGHBOR = NOT_ENUMERATED;

  /// View of a single cell in the grid
  class CellView {
  private:
    const CellGrid *grid_;
    size_type enumeration_index_;

  public:
    CellView() : grid_(nullptr), enumeration_index_(NOT_ENUMERATED) {}

    CellView(const CellGrid *grid, size_type enumerationIndex)
        : grid_(grid), enumeration_index_(enumerationIndex) {}

    [[nodiscard]] bool isValid() const noexcept {
      return grid_ != nullptr && enumeration_index_ != NOT_ENUMERATED &&
             enumeration_index_ < grid_->size();
    }

    explicit operator bool() const noexcept { return isValid(); }

    operator size_type() const noexcept { return enumeration_index_; }

    [[nodiscard]] size_type enumerationIndex() const noexcept {
      return enumeration_index_;
    }

    [[nodiscard]] Vec3D center() const;
    [[nodiscard]] Box3D boundingBox() const;
    [[nodiscard]] CellOctree::CellView octreeCell() const;
    [[nodiscard]] CellView neighbor(const Vec<std::ptrdiff_t, 3> &offset) const;
    [[nodiscard]] MortonIndex mortonIndex() const;
    [[nodiscard]] size_type level() const;
  };

  /// Iterator for CellGrid
  class Iterator {
  private:
    const CellGrid *grid_;
    size_type index_;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = CellView;
    using difference_type = std::ptrdiff_t;
    using pointer = const value_type *;
    using reference = value_type;

    Iterator() : grid_(nullptr), index_(0) {}

    Iterator(const CellGrid *grid, size_type index)
        : grid_(grid), index_(index) {}

    reference operator*() const { return {grid_, index_}; }

    Iterator &operator++() {
      ++index_;
      return *this;
    }

    Iterator operator++(int) {
      Iterator temp = *this;
      ++(*this);
      return temp;
    }

    bool operator==(const Iterator &other) const {
      return grid_ == other.grid_ && index_ == other.index_;
    }

    bool operator!=(const Iterator &other) const { return !(*this == other); }
  };

  // Factory method for builder
  [[nodiscard]] static CellGridBuilder
  create(std::shared_ptr<const CellOctree> octree);

  // Accessors
  [[nodiscard]] const CellOctree &octree() const { return *octree_; }

  [[nodiscard]] std::span<const MortonIndex> mortonIndices() const noexcept {
    return morton_indices_;
  }

  [[nodiscard]] size_type size() const noexcept {
    return morton_indices_.size();
  }

  [[nodiscard]] size_type
  getEnumerationIndex(size_type streamIndex) const noexcept;

  [[nodiscard]] size_type
  getEnumerationIndex(const CellOctree::CellView &cell) const noexcept;

  [[nodiscard]] std::span<const size_type>
  neighborIndices(const Vec<std::ptrdiff_t, 3> &offset) const;

  [[nodiscard]] CellView operator[](size_type index) const {
    return {this, index};
  }

  // Range interface
  [[nodiscard]] Iterator begin() const noexcept { return {this, 0}; }

  [[nodiscard]] Iterator end() const noexcept { return {this, size()}; }

private:
  std::shared_ptr<const CellOctree> octree_;
  std::vector<MortonIndex> morton_indices_;
  std::unordered_map<size_type, size_type> stream_to_enumeration_;
  std::map<Vec<std::ptrdiff_t, 3>, std::vector<size_type>> adjacency_lists_;

  friend class CellGridBuilder;

  CellGrid(
      std::shared_ptr<const CellOctree> octree,
      std::vector<MortonIndex> mortonIndices,
      std::unordered_map<size_type, size_type> streamToEnumeration,
      std::map<Vec<std::ptrdiff_t, 3>, std::vector<size_type>> adjacencyLists)
      : octree_(std::move(octree)), morton_indices_(std::move(mortonIndices)),
        stream_to_enumeration_(std::move(streamToEnumeration)),
        adjacency_lists_(std::move(adjacencyLists)) {}
};

/// Builder for CellGrid
class CellGridBuilder {
private:
  std::shared_ptr<const CellOctree> octree_;
  std::vector<CellGrid::size_type> levels_;
  std::vector<Vec<std::ptrdiff_t, 3>> neighborhood_;
  mutable std::unique_ptr<PeriodicityMapper> periodicity_mapper_;

public:
  explicit CellGridBuilder(std::shared_ptr<const CellOctree> octree)
      : octree_(std::move(octree)),
        periodicity_mapper_(std::make_unique<NoPeriodicity>()) {}

  CellGridBuilder &levels(std::span<const CellGrid::size_type> levels) {
    levels_.assign(levels.begin(), levels.end());
    return *this;
  }

  CellGridBuilder &levels(std::initializer_list<CellGrid::size_type> levels) {
    levels_.assign(levels.begin(), levels.end());
    return *this;
  }

  CellGridBuilder &
  neighborhood(std::span<const Vec<std::ptrdiff_t, 3>> offsets) {
    neighborhood_.assign(offsets.begin(), offsets.end());
    return *this;
  }

  CellGridBuilder &
  neighborhood(std::initializer_list<Vec<std::ptrdiff_t, 3>> offsets) {
    neighborhood_.assign(offsets.begin(), offsets.end());
    return *this;
  }

  CellGridBuilder &
  periodicityMapper(std::unique_ptr<PeriodicityMapper> mapper) {
    periodicity_mapper_ = std::move(mapper);
    return *this;
  }

  template <typename MapperType>
  CellGridBuilder &periodicityMapper(const MapperType &mapper) {
    periodicity_mapper_ = std::make_unique<MapperType>(mapper);
    return *this;
  }

  CellGrid build() const;
};

} // namespace oktal

#endif // OKTAL_OCTREE_CELLGRID_HPP
