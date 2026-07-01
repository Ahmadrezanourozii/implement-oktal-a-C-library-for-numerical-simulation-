#ifndef OKTAL_LBM_LBMKERNELS_HPP
#define OKTAL_LBM_LBMKERNELS_HPP

#include "D3Q19.hpp"
#include "../data/GridVector.hpp"
#include "../octree/CellGrid.hpp"
#include <cmath>

namespace oktal::lbm {

namespace detail {
[[nodiscard]] inline double computeEquilibrium(size_t q, double rho, const Vec3D &u) {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  const double w = D3Q19::W[q];
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  const auto &c = D3Q19::CS[q]; // Vec<int64_t, 3>

  const double cu = double(c[0]) * u[0] + double(c[1]) * u[1] + double(c[2]) * u[2];
  const double uSq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

  const double A = D3Q19::C_S_SQ_RECIP * cu;
  const double B = 0.5 * D3Q19::C_S_SQ_RECIP * D3Q19::C_S_SQ_RECIP * cu * cu;
  const double C = -0.5 * D3Q19::C_S_SQ_RECIP * uSq;

  return w * rho * (1.0 + A + B + C);
}
} // namespace detail

class InitializePdfs {
public:
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  void operator()(CellGrid::CellView cell, D3Q19LatticeView pdfs,
                  GridVectorView<const double, 1> rho,
                  GridVectorView<const double, 3> u) const {
    const size_t cellIdx{cell};
    const double localRho = rho[cellIdx];
    const Vec3D localU{u[cellIdx, 0], u[cellIdx, 1], u[cellIdx, 2]};

    for (size_t q = 0; q < D3Q19::Q; ++q) {
      pdfs[cellIdx, q] = detail::computeEquilibrium(q, localRho, localU);
    }
  }
};

class ComputeMacroscopicQuantities {
public:
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  void operator()(CellGrid::CellView cell, D3Q19LatticeConstView pdfs,
                  GridVectorView<double, 1> rho, GridVectorView<double, 3> u) const {
    const size_t cellIdx{cell};
    double localRho = 0.0;
    Vec3D localU{0.0, 0.0, 0.0};

    for (size_t q = 0; q < D3Q19::Q; ++q) {
      localRho += pdfs[cellIdx, q];
    }
    rho[cellIdx] = localRho;

    if (localRho > 1e-12) {
      for (size_t q = 0; q < D3Q19::Q; ++q) {
        const double f = pdfs[cellIdx, q];
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
        const auto &c = D3Q19::CS[q];
        localU[0] += f * double(c[0]);
        localU[1] += f * double(c[1]);
        localU[2] += f * double(c[2]);
      }
      localU = localU / localRho;
    }

    u[cellIdx, 0] = localU[0];
    u[cellIdx, 1] = localU[1];
    u[cellIdx, 2] = localU[2];
  }
};

class Collide {
public:
  explicit Collide(double omega) : omega_(omega) {}

  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  void operator()(CellGrid::CellView cell, D3Q19LatticeView pdfs,
                  GridVectorView<const double, 1> rho,
                  GridVectorView<const double, 3> u) const {
    const size_t cellIdx{cell};
    const double localRho = rho[cellIdx];
    const Vec3D localU{u[cellIdx, 0], u[cellIdx, 1], u[cellIdx, 2]};
    const double omega = omega_;

    for (size_t q = 0; q < D3Q19::Q; ++q) {
      const double eq = detail::computeEquilibrium(q, localRho, localU);
      pdfs[cellIdx, q] = (1.0 - omega) * pdfs[cellIdx, q] + omega * eq;
    }
  }

private:
  double omega_;
};

class Stream {
public:
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  void operator()(CellGrid::CellView cell, D3Q19LatticeView pdfsDst,
                  D3Q19LatticeConstView pdfsSrc) const {
    const size_t cellIdx{cell};

    for (size_t q = 0; q < D3Q19::Q; ++q) {
      const size_t oppQ = D3Q19::opposite(q);
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
      const auto& offsetVec = D3Q19::CS[oppQ];
      Vec<std::ptrdiff_t, 3> offset{ std::ptrdiff_t(offsetVec[0]), std::ptrdiff_t(offsetVec[1]), std::ptrdiff_t(offsetVec[2]) };
      
      if (offset[0] == 0 && offset[1] == 0 && offset[2] == 0) {
          pdfsDst[cellIdx, q] = pdfsSrc[cellIdx, q];
          continue;
      }

      const CellGrid::CellView neighbor = cell.neighbor(offset);
      
      if (!neighbor.isValid()) {
           pdfsDst[cellIdx, q] = 0.0; 
           continue; 
      }

      const size_t neighborIdx{neighbor};
      pdfsDst[cellIdx, q] = pdfsSrc[neighborIdx, q];
    }
  }
};

} // namespace oktal::lbm

#endif // OKTAL_LBM_LBMKERNELS_HPP
