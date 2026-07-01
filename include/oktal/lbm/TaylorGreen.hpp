#ifndef OKTAL_LBM_TAYLORGREEN_HPP
#define OKTAL_LBM_TAYLORGREEN_HPP

#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/OctreeGeometry.hpp"
#include <cstdint>
#include <stdexcept>
#include <cmath>
#include <numbers>

namespace oktal::lbm {

class TaylorGreen {
public:
  constexpr TaylorGreen(size_t refinementLevel,
                        double latticeViscosity = 1. / 6.)
      : refinementLevel_{refinementLevel}, latticeViscosity_{latticeViscosity} {
    if (refinementLevel < 5) {
      throw std::invalid_argument{"Refinement level must be >= 5"};
    }
  }

  //  NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  [[nodiscard]] OctreeGeometry geometry() const {
    return {Vec3D(0.), 2. * std::numbers::pi};
  }

  [[nodiscard]] constexpr double rho(Vec3D x, double t) const {
    const double uMax = latticeMaxVelocity();
    const double uMaxSq = uMax * uMax;
    const double decay = std::exp(-4. * physicalViscosity() * t);
    const double _3o4 = 3. / 4.;

    return 1. -
           _3o4 * uMaxSq * (std::cos(2. * x[0]) + std::cos(2. * x[1])) * decay;
  }

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  [[nodiscard]] constexpr Vec3D u(Vec3D x, double t) const {
    const double decay = std::exp(-2. * physicalViscosity() * t);
    return latticeMaxVelocity() * decay * Vec3D{
                                      std::sin(x[0]) * std::cos(x[1]),
                                      -std::cos(x[0]) * std::sin(x[1]),
                                      0.0,
                                  };
  }

  [[nodiscard]] constexpr double latticeViscosity() const {
    return latticeViscosity_;
  }

  [[nodiscard]] constexpr double omega() const {
    return 1. / (3. * latticeViscosity_ + 0.5);
  }

  [[nodiscard]] constexpr double dx() const {
    return 2.0 * std::numbers::pi / double(1 << refinementLevel_);
  }

  [[nodiscard]] constexpr double latticeMaxVelocity() const {
    return 0.04 / double(1 << (refinementLevel_ - 5));
  }

  [[nodiscard]] constexpr double dt() const {
    return latticeMaxVelocity() * dx();
  }

  [[nodiscard]] constexpr double dNu() const { return dx() * dx() / dt(); }

  [[nodiscard]] constexpr double physicalViscosity() const {
    return latticeViscosity_ * dNu();
  }

  [[nodiscard]] constexpr double decayTime() const {
    return 1. / (latticeViscosity_ * 2. * dx() * dx());
  }

  [[nodiscard]] constexpr size_t numberOfTimesteps() const {
    return size_t(1.93 * decayTime());
  }

private:
  size_t refinementLevel_;
  double latticeViscosity_;
};
} // namespace oktal::lbm

#endif // OKTAL_LBM_TAYLORGREEN_HPP
