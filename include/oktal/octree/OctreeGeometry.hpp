#ifndef OKTAL_OCTREE_OCTREEGEOMETRY_HPP
#define OKTAL_OCTREE_OCTREEGEOMETRY_HPP

#include "oktal/geometry/Box.hpp"
#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include <cstddef>
#include <cstdint>

namespace oktal {

class OctreeGeometry {
private:
  Vec3D origin_;
  double sidelength_;

public:
  constexpr OctreeGeometry() noexcept : origin_(0.0), sidelength_(1.0) {}

  constexpr OctreeGeometry(const Vec3D &origin, double sidelength) noexcept
      : origin_(origin), sidelength_(sidelength) {}

  [[nodiscard]] constexpr Vec3D origin() const noexcept { return origin_; }

  [[nodiscard]] constexpr double sidelength() const noexcept {
    return sidelength_;
  }

  [[nodiscard]] constexpr double dx(std::size_t level) const noexcept {
    const std::uint64_t divisions = std::uint64_t{1} << level;
    const double cell_size = sidelength_ / static_cast<double>(divisions);
    return cell_size;
  }

  [[nodiscard]] static constexpr Vec3D cellExtents(double sidelength,
                                                   std::size_t level) noexcept {
    // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
    const double h =
        sidelength / static_cast<double>(std::uint64_t{1} << level);
    Vec3D extents(h);
    return extents;
  }

  [[nodiscard]] constexpr Vec3D cellExtents(std::size_t level) const noexcept {
    return cellExtents(sidelength_, level);
  }

  [[nodiscard]] static Vec3D cellMinCorner(const Vec3D &origin,
                                           double sidelength,
                                           const MortonIndex &m) noexcept {
    const std::size_t level = m.level();
    // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
    const double h =
        sidelength / static_cast<double>(std::uint64_t{1} << level);
    const Vec<std::size_t, 3> coords = m.gridCoordinates();
    Vec3D min_corner = origin;
    min_corner[0] = min_corner[0] + static_cast<double>(coords[0]) * h;
    min_corner[1] = min_corner[1] + static_cast<double>(coords[1]) * h;
    min_corner[2] = min_corner[2] + static_cast<double>(coords[2]) * h;
    return min_corner;
  }

  [[nodiscard]] Vec3D cellMinCorner(const MortonIndex &m) const noexcept {
    return cellMinCorner(origin_, sidelength_, m);
  }

  [[nodiscard]] static Vec3D cellMaxCorner(const Vec3D &origin,
                                           double sidelength,
                                           const MortonIndex &m) noexcept {
    const Vec3D min_corner = cellMinCorner(origin, sidelength, m);
    const Vec3D extents = cellExtents(sidelength, m.level());
    Vec3D max_corner = min_corner + extents;
    return max_corner;
  }

  [[nodiscard]] Vec3D cellMaxCorner(const MortonIndex &m) const noexcept {
    return cellMaxCorner(origin_, sidelength_, m);
  }

  [[nodiscard]] static Box<double>
  cellBoundingBox(const Vec3D &origin, double sidelength,
                  const MortonIndex &m) noexcept {
    const Vec3D min_corner = cellMinCorner(origin, sidelength, m);
    const Vec3D max_corner = cellMaxCorner(origin, sidelength, m);
    Box<double> bbox(min_corner, max_corner);
    return bbox;
  }

  [[nodiscard]] Box<double>
  cellBoundingBox(const MortonIndex &m) const noexcept {
    return cellBoundingBox(origin_, sidelength_, m);
  }

  [[nodiscard]] static Vec3D cellCenter(const Vec3D &origin, double sidelength,
                                        const MortonIndex &m) noexcept {
    const Vec3D min_corner = cellMinCorner(origin, sidelength, m);
    const Vec3D extents = cellExtents(sidelength, m.level());
    const Vec3D half_extents = extents / 2.0;
    Vec3D center = min_corner + half_extents;
    return center;
  }

  [[nodiscard]] Vec3D cellCenter(const MortonIndex &m) const noexcept {
    return cellCenter(origin_, sidelength_, m);
  }
};

} // namespace oktal

#endif // OKTAL_OCTREE_OCTREEGEOMETRY_HPP