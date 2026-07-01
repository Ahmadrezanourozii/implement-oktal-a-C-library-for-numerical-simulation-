#ifndef OKTAL_GEOMETRY_VEC_HPP
#define OKTAL_GEOMETRY_VEC_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <initializer_list>

namespace oktal {

template <typename T>
concept Arithmetic = std::integral<T> || std::floating_point<T>;

template <Arithmetic T, std::size_t DIM> class Vec {
private:
  std::array<T, DIM> v_{};

public:
  constexpr Vec() noexcept = default;

  constexpr explicit Vec(T value) noexcept { v_.fill(value); }

  constexpr Vec(std::initializer_list<T> values) noexcept {
    const auto copy_count = std::min(values.size(), DIM);
    std::copy_n(values.begin(), copy_count, v_.begin());
    if (copy_count < DIM) {
      std::fill(v_.begin() + static_cast<std::ptrdiff_t>(copy_count), v_.end(),
                T());
    }
  }

  template <typename S>
    requires std::convertible_to<S, T>
  constexpr Vec(const Vec<S, DIM> &other) noexcept {
    for (std::size_t i = 0; i < DIM; i++) {
      v_.at(i) = static_cast<T>(other[i]);
    }
  }

  [[nodiscard]] constexpr std::size_t size() const noexcept { return DIM; }

  [[nodiscard]] T *data() noexcept { return v_.data(); }

  [[nodiscard]] const T *data() const noexcept { return v_.data(); }

  T &operator[](std::size_t index) { return v_.at(index); }

  const T &operator[](std::size_t index) const { return v_.at(index); }

  [[nodiscard]] auto begin() noexcept { return v_.begin(); }

  [[nodiscard]] auto begin() const noexcept { return v_.begin(); }

  [[nodiscard]] auto end() noexcept { return v_.end(); }

  [[nodiscard]] auto end() const noexcept { return v_.end(); }

  bool operator==(const Vec &other) const {
    for (std::size_t i = 0; i < DIM; i++) {
      if (v_.at(i) != other[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const Vec &other) const { return !(*this == other); }

  bool operator<(const Vec &other) const {
    for (std::size_t i = 0; i < DIM; i++) {
      if (v_.at(i) < other[i]) {
        return true;
      }
      if (v_.at(i) > other[i]) {
        return false;
      }
    }
    return false; // Equal
  }

  template <typename U = T>
    requires std::signed_integral<U> || std::floating_point<U>
  Vec operator-() const {
    Vec result;
    for (std::size_t i = 0; i < DIM; i++) {
      result[i] = -v_.at(i);
    }
    return result;
  }

  Vec &operator+=(const Vec &other) {
    for (std::size_t i = 0; i < DIM; i++) {
      v_.at(i) = v_.at(i) + other[i];
    }
    return *this;
  }

  Vec &operator-=(const Vec &other) {
    for (std::size_t i = 0; i < DIM; i++) {
      v_.at(i) = v_.at(i) - other[i];
    }
    return *this;
  }

  Vec &operator*=(T scalar) {
    for (std::size_t i = 0; i < DIM; i++) {
      v_.at(i) = v_.at(i) * scalar;
    }
    return *this;
  }

  Vec &operator/=(T scalar) {
    for (std::size_t i = 0; i < DIM; i++) {
      v_.at(i) = v_.at(i) / scalar;
    }
    return *this;
  }

  Vec operator+(const Vec &other) const {
    Vec result;
    for (std::size_t i = 0; i < DIM; i++) {
      result[i] = v_.at(i) + other[i];
    }
    return result;
  }

  Vec operator-(const Vec &other) const {
    Vec result;
    for (std::size_t i = 0; i < DIM; i++) {
      result[i] = v_.at(i) - other[i];
    }
    return result;
  }

  Vec operator*(T scalar) const {
    Vec result;
    for (std::size_t i = 0; i < DIM; i++) {
      result[i] = v_.at(i) * scalar;
    }
    return result;
  }

  Vec operator/(T scalar) const {
    Vec result;
    for (std::size_t i = 0; i < DIM; i++) {
      result[i] = v_.at(i) / scalar;
    }
    return result;
  }

  friend Vec operator*(T scalar, const Vec &vec) { return vec * scalar; }

  [[nodiscard]] T sqrMagnitude() const {
    T sum = T();
    for (std::size_t i = 0; i < DIM; i++) {
      const auto value = v_.at(i);
      sum = sum + (value * value);
    }
    return sum;
  }

  template <typename U = T>
    requires std::floating_point<U>
  [[nodiscard]] T magnitude() const {
    const auto sqrMag = static_cast<double>(sqrMagnitude());
    const auto mag = std::sqrt(sqrMag);
    return static_cast<T>(mag);
  }

  template <typename U = T>
    requires std::floating_point<U>
  [[nodiscard]] T sqrtMagnitude() const {
    return magnitude();
  }
};

using Vec3F = Vec<float, 3>;
using Vec3D = Vec<double, 3>;

} // namespace oktal

#endif // OKTAL_GEOMETRY_VEC_HPP