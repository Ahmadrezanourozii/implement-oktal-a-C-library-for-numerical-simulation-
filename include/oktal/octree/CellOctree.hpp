#ifndef OKTAL_OCTREE_CELLOCTREE_HPP
#define OKTAL_OCTREE_CELLOCTREE_HPP

#include "oktal/geometry/Box.hpp"
#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include "oktal/octree/OctreeGeometry.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace oktal {

using Vec3D = Vec<double, 3>;
using Box3D = Box<double>;

// Forward declarations
class OctreeCursor;

/// Concept for octree iterator policies
template <typename T>
concept OctreeIteratorPolicy =
    std::semiregular<T> && requires(const T t, OctreeCursor &c) {
      { t.advance(c) } -> std::same_as<void>;
    };

template <OctreeIteratorPolicy TPolicy> class OctreeCellsRange;

/// Main octree data structure using linear array representation
class CellOctree {
public:
  using size_type = std::size_t;

  static constexpr size_type OCTREE_CHILDREN_COUNT = 8;

  /// Compact node representation using bit packing
  struct Node {
  private:
    std::uint64_t data_{0};

    static constexpr size_type kRefinedBit = 0;
    static constexpr size_type kPhantomBit = 1;
    static constexpr size_type kIndexShift = 2;
    static constexpr std::uint64_t kIndexMask =
        (std::numeric_limits<std::uint64_t>::max() << kIndexShift);

  public:
    Node() = default;

    Node(bool isRefined, bool isPhantom, size_type childrenIndex) {
      setRefined(isRefined);
      setPhantom(isPhantom);
      setChildrenStartIndex(childrenIndex);
    }

    [[nodiscard]] bool isRefined() const noexcept {
      return (data_ & (1ULL << kRefinedBit)) != 0;
    }

    void setRefined(bool value) noexcept {
      if (value) {
        data_ |= (1ULL << kRefinedBit);
      } else {
        data_ &= ~(1ULL << kRefinedBit);
      }
    }

    [[nodiscard]] bool isPhantom() const noexcept {
      return (data_ & (1ULL << kPhantomBit)) != 0;
    }

    void setPhantom(bool value) noexcept {
      if (value) {
        data_ |= (1ULL << kPhantomBit);
      } else {
        data_ &= ~(1ULL << kPhantomBit);
      }
    }

    [[nodiscard]] size_type childrenStartIndex() const noexcept {
      return static_cast<size_type>((data_ & kIndexMask) >> kIndexShift);
    }

    void setChildrenStartIndex(size_type index) noexcept {
      data_ &= ~kIndexMask;
      data_ |= (static_cast<std::uint64_t>(index) << kIndexShift) & kIndexMask;
    }

    [[nodiscard]] size_type childIndex(size_type n) const noexcept {
      return childrenStartIndex() + n;
    }
  };

  /// View of a single octree cell
  class Cell {
  public:
    Cell(const CellOctree *octree, size_type streamIndex,
         MortonIndex mortonIndex)
        : octree_(octree), stream_index_(streamIndex),
          morton_index_(mortonIndex) {}

    [[nodiscard]] bool isRoot() const noexcept {
      return morton_index_.getBits() == 0b1uz;
    }

    [[nodiscard]] bool isRefined() const noexcept {
      return octree_->nodesStream()[stream_index_].isRefined();
    }

    [[nodiscard]] size_type level() const noexcept {
      return morton_index_.level();
    }

    [[nodiscard]] size_type streamIndex() const noexcept {
      return stream_index_;
    }

    [[nodiscard]] MortonIndex mortonIndex() const noexcept {
      return morton_index_;
    }

    [[nodiscard]] Vec3D center() const noexcept;
    [[nodiscard]] Box3D boundingBox() const noexcept;

  private:
    const CellOctree *octree_;
    size_type stream_index_;
    MortonIndex morton_index_;
  };

  using CellView = Cell;

  // Constructors
  CellOctree() : geometry_(OctreeGeometry(Vec3D(0.0), 1.0)) {
    initialize_root();
  }

  explicit CellOctree(OctreeGeometry geometry) : geometry_(geometry) {
    initialize_root();
  }

  // Factory methods
  [[nodiscard]] static CellOctree fromDescriptor(std::string_view descr);
  [[nodiscard]] static CellOctree fromDescriptor(std::string_view descr,
                                                 OctreeGeometry geometry);
  [[nodiscard]] static std::shared_ptr<const CellOctree>
  createUniformGrid(size_type level);
  [[nodiscard]] static std::shared_ptr<const CellOctree>
  createUniformGrid(OctreeGeometry geom, size_type level);

  // Cell access
  [[nodiscard]] std::optional<Cell> getRootCell() const noexcept;
  [[nodiscard]] std::optional<Cell> getCell(MortonIndex m) const noexcept;
  [[nodiscard]] bool cellExists(MortonIndex m) const noexcept;

  // Geometry and structure queries
  [[nodiscard]] const OctreeGeometry &geometry() const noexcept {
    return geometry_;
  }

  [[nodiscard]] size_type numberOfNodes() const noexcept {
    return nodes_.size();
  }

  [[nodiscard]] size_type numberOfLevels() const noexcept {
    return level_start_indices_.size();
  }

  [[nodiscard]] size_type numberOfNodes(size_type level) const noexcept;

  [[nodiscard]] std::span<const Node> nodesStream() const noexcept {
    return nodes_;
  }

  [[nodiscard]] std::span<const Node>
  nodesStream(size_type level) const noexcept;

  // Range factory methods
  [[nodiscard]] auto preOrderDepthFirstRange() const;
  [[nodiscard]] auto horizontalRange(size_t level) const;

private:
  std::vector<Node> nodes_;
  std::vector<size_type> level_start_indices_;
  OctreeGeometry geometry_;

  void initialize_root();
  std::span<Node> nodesStream(size_type level, bool mutable_access);

  [[nodiscard]] std::pair<size_type, size_type>
  find_node(MortonIndex m) const noexcept;
};

// -----------------------------------------------------------------------------
// Octree Cursor - Provides navigation through the octree structure
// -----------------------------------------------------------------------------

class OctreeCursor {
public:
  OctreeCursor() = default;

  explicit OctreeCursor(const CellOctree &octree)
      : octree_(&octree), path_{0} {}

  OctreeCursor(const CellOctree &octree, std::initializer_list<size_t> path)
      : octree_(&octree), path_(path) {}

  template <typename Range>
    requires std::ranges::input_range<Range> &&
                 std::convertible_to<std::ranges::range_value_t<Range>, size_t>
  OctreeCursor(const CellOctree &octree, const Range &path)
      : octree_(&octree), path_(std::begin(path), std::end(path)) {}

  // Accessors
  [[nodiscard]] const CellOctree *octree() const noexcept { return octree_; }

  [[nodiscard]] std::span<const size_t> path() const noexcept { return path_; }

  // State queries
  [[nodiscard]] bool empty() const noexcept { return octree_ == nullptr; }

  [[nodiscard]] bool end() const noexcept {
    return octree_ != nullptr && path_.empty();
  }

  [[nodiscard]] size_t currentLevel() const {
    return path_.empty() ? 0 : path_.size() - 1;
  }

  [[nodiscard]] size_t currentStreamIndex() const { return path_.back(); }

  [[nodiscard]] std::optional<CellOctree::CellView> currentCell() const {
    if (empty() || end()) {
      return std::nullopt;
    }

    const auto &node = getCurrentNode();
    if (node.isPhantom()) {
      return std::nullopt;
    }

    return CellOctree::CellView{octree_, currentStreamIndex(), mortonIndex()};
  }

  [[nodiscard]] bool firstSibling() const {
    if (path_.size() <= 1) {
      return true;
    }

    const auto parentIdx = path_[path_.size() - 2];
    const auto &parentNode = octree_->nodesStream()[parentIdx];
    return currentStreamIndex() == parentNode.childrenStartIndex();
  }

  [[nodiscard]] bool lastSibling() const {
    if (path_.size() <= 1) {
      return true;
    }

    const auto parentIdx = path_[path_.size() - 2];
    const auto &parentNode = octree_->nodesStream()[parentIdx];
    return currentStreamIndex() ==
           parentNode.childIndex(CellOctree::OCTREE_CHILDREN_COUNT - 1);
  }

  [[nodiscard]] MortonIndex mortonIndex() const {
    if (empty() || end()) {
      return {};
    }

    size_t mortonBits = 1;
    for (size_t i = 1; i < path_.size(); ++i) {
      const auto parentIdx = path_[i - 1];
      const auto currentIdx = path_[i];
      const auto &parentNode = octree_->nodesStream()[parentIdx];
      const auto childOctant = currentIdx - parentNode.childrenStartIndex();
      mortonBits = (mortonBits << 3) | childOctant;
    }
    return MortonIndex{mortonBits};
  }

  // Comparison operators
  bool operator==(const OctreeCursor &other) const {
    if (octree_ != other.octree_) {
      return false;
    }
    if (octree_ == nullptr) {
      return true;
    }
    return path_ == other.path_;
  }

  bool operator!=(const OctreeCursor &other) const { return !(*this == other); }

  // Navigation methods
  void ascend() {
    if (!path_.empty()) {
      path_.pop_back();
    }
  }

  void descend() {
    const auto &node = getCurrentNode();
    if (node.isRefined()) {
      path_.push_back(node.childrenStartIndex());
    }
  }

  void descend(size_t childIdx) {
    if (childIdx >= CellOctree::OCTREE_CHILDREN_COUNT) {
      throw std::out_of_range("Child index must be in range [0, 7]");
    }

    const auto &node = getCurrentNode();
    if (node.isRefined()) {
      path_.push_back(node.childIndex(childIdx));
    }
  }

  void previousSibling() {
    if (path_.size() > 1 && !firstSibling()) {
      --path_.back();
    }
  }

  void nextSibling() {
    if (path_.size() > 1 && !lastSibling()) {
      ++path_.back();
    }
  }

  void toSibling(size_t siblingIdx) {
    if (siblingIdx >= CellOctree::OCTREE_CHILDREN_COUNT) {
      throw std::out_of_range("Sibling index must be in range [0, 7]");
    }

    if (path_.size() == 1) {
      if (siblingIdx != 0) {
        throw std::out_of_range("Root node has no siblings");
      }
      return;
    }

    const auto parentIdx = path_[path_.size() - 2];
    const auto &parentNode = octree_->nodesStream()[parentIdx];
    path_.back() = parentNode.childIndex(siblingIdx);
  }

  void toEnd() { path_.clear(); }

private:
  const CellOctree *octree_{nullptr};
  std::vector<size_t> path_;

  [[nodiscard]] const CellOctree::Node &getCurrentNode() const {
    return octree_->nodesStream()[currentStreamIndex()];
  }
};

// -----------------------------------------------------------------------------
// Iterator and Range Types
// -----------------------------------------------------------------------------

template <OctreeIteratorPolicy TPolicy> class OctreeCellsIterator {
public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = CellOctree::CellView;
  using difference_type = std::ptrdiff_t;
  using pointer = const value_type *;
  using reference = value_type;

  OctreeCellsIterator() = default;

  OctreeCellsIterator(OctreeCursor cursor, TPolicy policy)
      : cursor_(std::move(cursor)), policy_(std::move(policy)) {}

  reference operator*() const {
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    return cursor_.currentCell().value();
  }

  OctreeCellsIterator &operator++() {
    policy_.advance(cursor_);
    return *this;
  }

  OctreeCellsIterator operator++(int) {
    auto temp = *this;
    ++(*this);
    return temp;
  }

  bool operator==(const OctreeCellsIterator &other) const {
    return cursor_ == other.cursor_;
  }

  bool operator!=(const OctreeCellsIterator &other) const {
    return !(*this == other);
  }

private:
  OctreeCursor cursor_;
  TPolicy policy_;
};

// Compatibility alias
template <typename T> using OctreeIterator = OctreeCellsIterator<T>;

template <OctreeIteratorPolicy TPolicy> class OctreeCellsRange {
public:
  OctreeCellsRange(OctreeCursor start, OctreeCursor end, TPolicy policy)
      : start_(std::move(start)), end_(std::move(end)),
        policy_(std::move(policy)) {}

  [[nodiscard]] OctreeCellsIterator<TPolicy> begin() const {
    return {start_, policy_};
  }

  [[nodiscard]] OctreeCellsIterator<TPolicy> end() const {
    return {end_, policy_};
  }

private:
  OctreeCursor start_;
  OctreeCursor end_;
  TPolicy policy_;
};

// -----------------------------------------------------------------------------
// Traversal Policies
// -----------------------------------------------------------------------------

/// Depth-first pre-order traversal policy
struct PreOrderDepthFirstPolicy {
  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  void advance(OctreeCursor &cursor) const {
    do {
      const size_t previousLevel = cursor.currentLevel();
      cursor.descend();

      if (cursor.currentLevel() > previousLevel) {
        // Successfully descended to child node
        continue;
      }

      // Cannot descend (leaf node), move to next sibling or backtrack
      while (!cursor.end() && cursor.lastSibling()) {
        cursor.ascend();
      }

      if (!cursor.end()) {
        cursor.nextSibling();
      }

    } while (!cursor.end() && !cursor.currentCell().has_value());
  }
};

/// Level-order (breadth-first) traversal policy for a specific level
class HorizontalPolicy {
public:
  HorizontalPolicy() = default;
  explicit HorizontalPolicy(size_t targetLevel) : target_level_(targetLevel) {}

  void advance(OctreeCursor &cursor) const {
    do {
      advanceInternal(cursor);
    } while (!cursor.end() && (!cursor.currentCell().has_value() ||
                               cursor.currentLevel() != target_level_));
  }

private:
  size_t target_level_{0};

  void advanceInternal(OctreeCursor &cursor) const {
    // Ascend if too deep or finished with current sibling group
    while (!cursor.end() &&
           (cursor.currentLevel() > target_level_ || cursor.lastSibling())) {
      if (cursor.currentLevel() == 0) {
        cursor.toEnd();
        return;
      }
      cursor.ascend();
    }

    if (cursor.end()) {
      return;
    }

    // Move to next sibling
    cursor.nextSibling();

    // Descend to target level
    while (!cursor.end() && cursor.currentLevel() < target_level_) {
      const size_t previousLevel = cursor.currentLevel();
      cursor.descend();

      if (cursor.currentLevel() == previousLevel) {
        // Could not descend (not refined) - dead end
        advanceInternal(cursor);
        return;
      }
    }
  }
};

// -----------------------------------------------------------------------------
// CellOctree Range Factory Method Implementations
// -----------------------------------------------------------------------------

inline auto CellOctree::preOrderDepthFirstRange() const {
  const PreOrderDepthFirstPolicy policy;
  OctreeCursor start(*this);
  OctreeCursor end(*this);
  end.toEnd();

  // Skip phantom root if necessary
  if (!start.empty() && !start.currentCell().has_value()) {
    policy.advance(start);
  }

  return OctreeCellsRange<PreOrderDepthFirstPolicy>(std::move(start),
                                                    std::move(end), policy);
}

inline auto CellOctree::horizontalRange(size_t level) const {
  const HorizontalPolicy policy(level);
  OctreeCursor start(*this);
  OctreeCursor end(*this);
  end.toEnd();

  if (!start.empty()) {
    // Descend to target level
    while (!start.end() && start.currentLevel() < level) {
      const size_t previousLevel = start.currentLevel();
      start.descend();

      if (start.currentLevel() == previousLevel) {
        // Cannot reach target level from this path
        policy.advance(start);
        break;
      }
    }

    // Validate current position (skip phantom nodes)
    if (!start.end() && start.currentLevel() == level) {
      if (!start.currentCell().has_value()) {
        policy.advance(start);
      }
    }
  }

  return OctreeCellsRange<HorizontalPolicy>(std::move(start), std::move(end),
                                            policy);
}

} // namespace oktal

#endif // OKTAL_OCTREE_CELLOCTREE_HPP