#ifndef OKTAL_GEOMETRY_BOX_HPP
#define OKTAL_GEOMETRY_BOX_HPP

#include "oktal/geometry/Vec.hpp"
#include <concepts>
#include <type_traits>

namespace oktal {

template <typename T = double>
  requires std::is_arithmetic_v<T> && std::is_scalar_v<T>
class Box {
public:
  using vector_type = Vec<T, 3>;

private:
  vector_type minCorner_;
  vector_type maxCorner_;

public:
  Box() : minCorner_(T(0)), maxCorner_(T(0)) {}

  Box(const vector_type &min, const vector_type &max)
      : minCorner_(min), maxCorner_(max) {}

  vector_type &minCorner() { return minCorner_; }

  [[nodiscard]] const vector_type &minCorner() const { return minCorner_; }

  vector_type &maxCorner() { return maxCorner_; }

  [[nodiscard]] const vector_type &maxCorner() const { return maxCorner_; }

  [[nodiscard]] vector_type center() const
    requires std::is_floating_point_v<T>
  {
    vector_type sum = minCorner_ + maxCorner_;
    vector_type result;

    for (std::size_t i = 0; i < 3; i++) {
      result[i] = sum[i] / T(2);
    }

    return result;
  }

  [[nodiscard]] T volume() const {
    vector_type ext = extents();
    T vol = ext[0] * ext[1] * ext[2];
    return vol;
  }

  [[nodiscard]] vector_type extents() const {
    vector_type result;

    for (std::size_t i = 0; i < 3; i++) {
      result[i] = maxCorner_[i] - minCorner_[i];
    }

    return result;
  }
};

} // namespace oktal

#endif // OKTAL_GEOMETRY_BOX_HPP