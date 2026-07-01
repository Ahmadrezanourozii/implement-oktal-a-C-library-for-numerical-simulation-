#pragma once

// This part is a bit tricky: it checks if the compiler already has mdspan.
// If not, it falls back to the Kokkos implementation so we don't crash.
#if __has_include(<experimental/mdspan>)
#include <experimental/mdspan>
#else
#include <mdspan/mdspan.hpp>
namespace std::experimental {
using namespace Kokkos;
}
#endif

#include "oktal/octree/CellGrid.hpp"
#include <algorithm>
#include <concepts>
#include <memory>
#include <type_traits>

// Safety bridge for Task 6 unit tests - forward declaration only
namespace advpt::htg {} // NOLINT(modernize-concat-nested-namespaces) - forward declaration needed

namespace oktal {

namespace stdex = std::experimental; // Shorter alias because
                                     // 'std::experimental' is a lot to type

namespace detail {
// Selects 1D extents for Q=1 and 2D extents for Q>1 to satisfy rank checks
template <size_t Q>
using GridVectorExtents =
    std::conditional_t<Q == 1, stdex::extents<size_t, stdex::dynamic_extent>,
                       stdex::extents<size_t, stdex::dynamic_extent, Q>>;

// We're using layout_left here because it's better for Structure-of-Arrays
// (SoA) patterns.
using GridVectorLayout = stdex::layout_left;
} // namespace detail

/**
 * Standalone view alias
 * Constrained to pass IsValidGridVectorView concept checks.
 */
// Using C++20 'requires' here to make sure T is a sane type and Q isn't zero.
template <typename T, size_t Q>
  requires(Q > 0) && std::semiregular<std::remove_const_t<T>> &&
              (!std::is_reference_v<T>)
using GridVectorView =
    stdex::mdspan<T, detail::GridVectorExtents<Q>, detail::GridVectorLayout>;

template <typename T, size_t Q>
  requires std::semiregular<T> && (Q > 0)
class GridVector {
public:
  using value_type = T;
  static constexpr size_t NUM_COMPONENTS = Q;

  // --- Container Semantics (Rule of Five) ---
  // Since we're managing a raw pointer (data_), we have to do all 5 special
  // functions.

  // 1. Constructor: Allocates memory for all cells times components.

  explicit GridVector(const CellGrid &cells)
      : num_cells_(cells.size()), data_(new T[cells.size() * Q]()) { // NOLINT(cppcoreguidelines-owning-memory)
  } // The () at the end zero-initializes the memory

  // 2. Destructor: Clean up the heap memory so we don't have leaks.
  ~GridVector() { delete[] data_; }

  // 3. Copy Constructor: Deep copy the data so two objects don't point to the
  // same memory.
  GridVector(const GridVector &other)
      : num_cells_(other.num_cells_), data_(new T[other.num_cells_ * Q]) { // NOLINT(cppcoreguidelines-owning-memory)
    std::copy_n(other.data_, num_cells_ * Q, data_);
  }

  // 4. Copy Assignment: Handle self-assignment and reallocate if the size
  // changed.
  GridVector &operator=(const GridVector &other) {
    if (this != &other) {
      if (num_cells_ != other.num_cells_) {
        delete[] data_;
        num_cells_ = other.num_cells_;
        data_ = new T[num_cells_ * Q]; // NOLINT(cppcoreguidelines-owning-memory)
      }
      std::copy_n(other.data_, num_cells_ * Q, data_);
    }
    return *this;
  }
  // 5. Move Constructor: Just "steal" the pointer from the other object. Super
  // fast!
  GridVector(GridVector &&other) noexcept
      : num_cells_(other.num_cells_), data_(other.data_) {
    other.num_cells_ = 0;
    other.data_ = nullptr; // Important: null this out so 'other''s destructor
                           // doesn't delete our data!
  }

  // Move Assignment: Delete our old data and take the new stuff.
  GridVector &operator=(GridVector &&other) noexcept {
    if (this != &other) {
      delete[] data_;
      num_cells_ = other.num_cells_;
      data_ = other.data_;
      other.num_cells_ = 0;
      other.data_ = nullptr;
    }
    return *this;
  }

  [[nodiscard]] size_t size() const noexcept { return num_cells_; }

  // Returns total elements (cells * components)
  [[nodiscard]] size_t allocSize() const noexcept { return num_cells_ * Q; }

  // Returns raw pointer for move-semantics verification
  [[nodiscard]] T *data() noexcept { return data_; }
  [[nodiscard]] const T *data() const noexcept { return data_; }

  // --- View Factory Functions ---
  // These give us an mdspan "window" into our raw data.
  [[nodiscard]] GridVectorView<T, Q> view() noexcept {
    return create_view<T>();
  }
  [[nodiscard]] GridVectorView<const T, Q> view() const noexcept {
    return create_view<const T>();
  }

  [[nodiscard]] GridVectorView<const T, Q> const_view() const noexcept {
    return view();
  }

  // Conversion operators so we can pass a GridVector into a function expecting
  // a GridVectorView
  operator GridVectorView<T, Q>() noexcept { return view(); }
  operator GridVectorView<const T, Q>() const noexcept { return view(); }

  // --- Accessors ---
  // we used template requirements here so the compiler knows which operator[]
  // to use based on whether we have 1 component or multiple.

  // 1D access (e.g., vec[i])
  template <size_t D = Q>
    requires(D == 1)
  T &operator[](size_t i) {
    return view()[i];
  }

  template <size_t D = Q>
    requires(D == 1)
  const T &operator[](size_t i) const {
    return view()[i];
  }

  // 2D access (e.g., vec[i, q]) - C++20 allows multiple arguments in
  // operator[]!
  template <size_t D = Q>
    requires(D > 1)
  T &operator[](size_t i, size_t q) {
    return view()[i, q];
  }

  template <size_t D = Q>
    requires(D > 1)
  const T &operator[](size_t i, size_t q) const {
    return view()[i, q];
  }

private:
  size_t num_cells_;
  T *data_; // NOLINT(cppcoreguidelines-owning-memory) - raw pointer for performance

  // Helper to build the mdspan. We use reinterpret_cast to handle the
  // const/non-const transitions.
  template <typename U>
  [[nodiscard]] auto create_view() const {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) - necessary for const/non-const conversion
    U *ptr = reinterpret_cast<U *>(data_);
    return GridVectorView<U, Q>(ptr, num_cells_);
  }
};

} // namespace oktal
