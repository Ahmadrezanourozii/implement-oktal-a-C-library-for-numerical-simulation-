#pragma once

#include <array>
#include <cmath>
#include <cstdint>

namespace oktal {

/**
 * @brief Partially periodic cuboid volume exposing periodic point mapping and
 * distance computations.
 */
class PeriodicBox {
public:
  /**
   * @brief Construct a periodic box
   *
   * @param minCorner The box's minimum (left-front-bottom) corner
   * @param maxCorner The box's maximum (right-back-top) corner
   * @param periodicity The box's periodicity in x-, y-, and z-direction
   */
  PeriodicBox(std::array<double, 3> minCorner, std::array<double, 3> maxCorner,
              std::array<bool, 3> periodicity)
      : minCorner_{minCorner}, maxCorner_{maxCorner},
        periodicity_{periodicity} {}

  /**
   * @brief The box's minimum (left-front-bottom) corner
   */
  [[nodiscard]] const std::array<double, 3> &minCorner() const {
    return minCorner_;
  }

  /**
   * @brief The box's maximum (right-back-top) corner
   */
  [[nodiscard]] const std::array<double, 3> &maxCorner() const {
    return maxCorner_;
  }

  /**
   * @brief The box's periodicity in x-, y-, and z-direction
   */
  [[nodiscard]] const std::array<bool, 3> &periodicity() const {
    return periodicity_;
  }

  /**
   * @brief Map a point in 3D space to its image with respect to this periodic
   * box
   */
  [[nodiscard]] std::array<double, 3>
  mapIntoBox(std::array<double, 3> point) const;

  /**
   * @brief Compute the square euclidean distance of two points with respect to
   * the periodic box
   *
   * @warning Results are undefined for points outside of the box.
   * Transform them into the box using `mapIntoBox` first.
   */
  [[nodiscard]] double sqrDistance(std::array<double, 3> pointA,
                                   std::array<double, 3> pointB) const;

  /**
   * @brief Compute the euclidean distance of two points with respect to the
   * periodic box
   *
   * @warning Results are undefined for points outside of the box.
   * Transform them into the box using `mapIntoBox` first.
   */
  [[nodiscard]] double distance(std::array<double, 3> pointA,
                                std::array<double, 3> pointB) const {
    return std::sqrt(sqrDistance(pointA, pointB));
  }

private:
  std::array<double, 3> minCorner_;
  std::array<double, 3> maxCorner_;
  std::array<bool, 3> periodicity_;
};

} // namespace oktal
