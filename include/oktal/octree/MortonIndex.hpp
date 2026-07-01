#ifndef OKTAL_OCTREE_MORTONINDEX_HPP
#define OKTAL_OCTREE_MORTONINDEX_HPP

#include "oktal/geometry/Vec.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace oktal {

using morton_bits_t = std::uint64_t;

class MortonIndex {
private:
  morton_bits_t bits_{1};

public:
  static constexpr std::size_t MAX_DEPTH = (sizeof(morton_bits_t) * 8 - 1) / 3;

  MortonIndex() = default;

  explicit MortonIndex(morton_bits_t bits) : bits_(bits) {}

  [[nodiscard]] morton_bits_t getBits() const { return bits_; }

  static MortonIndex fromPath(const std::vector<morton_bits_t> &path) {
    if (path.size() > MAX_DEPTH) {
      throw std::invalid_argument("Path too long for MortonIndex");
    }

    morton_bits_t b = 1;

    for (const morton_bits_t choice : path) {
      if (choice > 7) {
        throw std::invalid_argument("Invalid child index in path");
      }

      b = (b << 3) | choice;
    }

    return MortonIndex(b);
  }

  template <typename Range> static MortonIndex fromPath(const Range &path) {
    const auto path_length = static_cast<std::size_t>(
        std::ranges::distance(path.begin(), path.end()));

    if (path_length > MAX_DEPTH) {
      throw std::invalid_argument("Path too long for MortonIndex");
    }

    morton_bits_t b = 1;

    for (const auto &value : path) {
      const auto choice = static_cast<morton_bits_t>(value);

      if (choice > 7) {
        throw std::invalid_argument("Invalid child index in path");
      }

      b = (b << 3) | choice;
    }

    return MortonIndex(b);
  }

  [[nodiscard]] static std::vector<morton_bits_t> getPath(morton_bits_t bits) {
    std::vector<morton_bits_t> path_backwards;

    morton_bits_t b = bits;

    // NOLINTNEXTLINE(bugprone-infinite-loop)
    while (b > 1) {
      const morton_bits_t child_idx = b & 0b111u;
      path_backwards.push_back(child_idx);

      b = b >> 3;
    }

    std::ranges::reverse(path_backwards);

    return path_backwards;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  [[nodiscard]] std::vector<morton_bits_t> getPath() const {
    return getPath(bits_);
  }

  [[nodiscard]] static std::size_t level(morton_bits_t bits) {
    std::size_t depth = 0;
    morton_bits_t b = bits;

    // NOLINTNEXTLINE(bugprone-infinite-loop)
    while (b > 1) {
      b = b >> 3;
      depth = depth + 1;
    }

    return depth;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  [[nodiscard]] std::size_t level() const { return level(bits_); }

  [[nodiscard]] bool isRoot() const { return bits_ == 1; }

  [[nodiscard]] std::size_t siblingIndex() const {
    if (isRoot()) {
      return 0;
    }
    return static_cast<std::size_t>(bits_ & 0b111u);
  }

  [[nodiscard]] bool isFirstSibling() const {
    if (isRoot()) {
      return true;
    }
    return siblingIndex() == 0;
  }

  [[nodiscard]] bool isLastSibling() const {
    if (isRoot()) {
      return true;
    }
    return siblingIndex() == 7;
  }

  [[nodiscard]] static MortonIndex parent(morton_bits_t bits) {
    const morton_bits_t parent_bits = bits >> 3;
    return MortonIndex(parent_bits);
  }

  [[nodiscard]] MortonIndex parent() const { return parent(bits_); }

  [[nodiscard]] MortonIndex safeParent() const {
    if (isRoot()) {
      throw std::logic_error("Root has no parent");
    }
    return parent();
  }

  [[nodiscard]] static MortonIndex child(morton_bits_t bits,
                                         morton_bits_t idx) {
    const morton_bits_t child_bits = (bits << 3) | (idx & 0b111u);
    return MortonIndex(child_bits);
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  [[nodiscard]] MortonIndex child(morton_bits_t idx) const {
    return child(bits_, idx);
  }

  [[nodiscard]] static MortonIndex safeChild(morton_bits_t bits,
                                             morton_bits_t idx) {
    if (idx > 7) {
      throw std::logic_error("Invalid child index");
    }
    if (level(bits) >= MAX_DEPTH) {
      throw std::logic_error("Maximum depth exceeded");
    }
    const morton_bits_t child_bits = (bits << 3) | idx;
    return MortonIndex(child_bits);
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  [[nodiscard]] MortonIndex safeChild(morton_bits_t idx) const {
    return safeChild(bits_, idx);
  }

  friend bool operator==(MortonIndex a, MortonIndex b) {
    return a.bits_ == b.bits_;
  }

  friend bool operator!=(MortonIndex a, MortonIndex b) { return !(a == b); }

  friend bool operator>(MortonIndex a, MortonIndex b) {
    if (a.bits_ == b.bits_) {
      return false;
    }

    const std::size_t level_a = a.level();
    const std::size_t level_b = b.level();

    if (level_a >= level_b) {
      return false;
    }

    const std::size_t level_diff = level_b - level_a;
    const morton_bits_t b_ancestor = b.bits_ >> (3 * level_diff);

    return b_ancestor == a.bits_;
  }

  friend bool operator<(MortonIndex a, MortonIndex b) { return b > a; }

  friend bool operator>=(MortonIndex a, MortonIndex b) {
    return (a == b) || (a > b);
  }

  friend bool operator<=(MortonIndex a, MortonIndex b) {
    return (a == b) || (a < b);
  }

  [[nodiscard]] static Vec<std::size_t, 3> gridCoordinates(morton_bits_t bits) {
    Vec<std::size_t, 3> coords(0);

    const std::vector<morton_bits_t> path = getPath(bits);

    for (const morton_bits_t direction : path) {
      coords[0] = coords[0] << 1;
      coords[1] = coords[1] << 1;
      coords[2] = coords[2] << 1;

      coords[0] = coords[0] | ((direction >> 0) & 1u);
      coords[1] = coords[1] | ((direction >> 1) & 1u);
      coords[2] = coords[2] | ((direction >> 2) & 1u);
    }

    return coords;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  [[nodiscard]] Vec<std::size_t, 3> gridCoordinates() const {
    return gridCoordinates(bits_);
  }

  [[nodiscard]] static MortonIndex
  fromGridCoordinates(std::size_t level, const Vec<std::size_t, 3> &coords) {
    if (level > MAX_DEPTH) {
      throw std::invalid_argument("Level too deep for MortonIndex");
    }

    morton_bits_t b = 1;

    // NOLINTNEXTLINE(bugprone-infinite-loop)
    for (std::size_t i = level; i > 0; --i) {
      const std::size_t bit_pos = i - 1;

      morton_bits_t direction = 0;

      const morton_bits_t x_bit = (coords[0] >> bit_pos) & 1u;
      const morton_bits_t y_bit = (coords[1] >> bit_pos) & 1u;
      const morton_bits_t z_bit = (coords[2] >> bit_pos) & 1u;

      direction = direction | (x_bit << 0);
      direction = direction | (y_bit << 1);
      direction = direction | (z_bit << 2);

      b = (b << 3) | direction;
    }

    return MortonIndex(b);
  }
};

} // namespace oktal

#endif // OKTAL_OCTREE_MORTONINDEX_HPP