#ifndef OKTAL_LBM_D3Q19_HPP
#define OKTAL_LBM_D3Q19_HPP

#include <array>
#include <cstdint>
#include <cstddef> // for ptrdiff_t

#include "../data/GridVector.hpp"
#include "../geometry/Vec.hpp"

namespace oktal::lbm {

struct D3Q19 {
public:
  /* Dimensionality*/
  static constexpr size_t D = 3;

  /* Number of velocities */
  static constexpr size_t Q = 19;

  /* Square speed of sound */
  static constexpr double C_S_SQ = 1.0 / 3.0;

  /* Reciprocal value of the square speed of sound (`1 / C_S_SQ`) */
  static constexpr double C_S_SQ_RECIP = 3.0;

  /* Lattice velocities */
  // Use ptrdiff_t to match CellGrid API
  static constexpr std::array<Vec<std::ptrdiff_t, D>, Q> CS = {{
      {0L, 0L, 0L}, // C,

      {1L, 0L, 0L},  // W
      {-1L, 0L, 0L}, // E

      {0L, 1L, 0L},  // N
      {0L, -1L, 0L}, // S

      {0L, 0L, 1L},  // T
      {0L, 0L, -1L}, // B

      {1L, 1L, 0L},   // NW
      {-1L, -1L, 0L}, // SE

      {1L, -1L, 0L}, // SW
      {-1L, 1L, 0L}, // NE

      {1L, 0L, 1L},   // TW
      {-1L, 0L, -1L}, // BE

      {1L, 0L, -1L}, // BW
      {-1L, 0L, 1L}, // TE

      {0L, 1L, 1L},   // TN
      {0L, -1L, -1L}, // BS

      {0L, 1L, -1L}, // BN
      {0L, -1L, 1L}, // TS
  }};

  /* Lattice velocities without center */
  static constexpr std::array<Vec<std::ptrdiff_t, D>, Q - 1> CS_NO_CENTER = {{
      {1L, 0L, 0L},  // W
      {-1L, 0L, 0L}, // E

      {0L, 1L, 0L},  // N
      {0L, -1L, 0L}, // S

      {0L, 0L, 1L},  // T
      {0L, 0L, -1L}, // B

      {1L, 1L, 0L},   // NW
      {-1L, -1L, 0L}, // SE

      {1L, -1L, 0L}, // SW
      {-1L, 1L, 0L}, // NE

      {1L, 0L, 1L},   // TW
      {-1L, 0L, -1L}, // BE

      {1L, 0L, -1L}, // BW
      {-1L, 0L, 1L}, // TE

      {0L, 1L, 1L},   // TN
      {0L, -1L, -1L}, // BS

      {0L, 1L, -1L}, // BN
      {0L, -1L, 1L}, // TS
  }};

  /* Lattice weights */
  static constexpr std::array<double, Q> W = {{
      //  Stationary
      1. / 3.,
      // Cardinals
      1. / 18.,
      1. / 18.,
      1. / 18.,
      1. / 18.,
      1. / 18.,
      1. / 18.,
      // Diagonals
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
      1. / 36.,
  }};

  /* Get the index of the opposite velocity */
  [[nodiscard]] static constexpr size_t opposite(size_t i) {
    if (i == 0uz) {
      return 0uz;
    }

    //  Invert least significant bit
    const size_t tmp{i - 1uz};
    return 1uz + ((tmp & 0xffffffffffffffeuz) | (tmp ^ 0x1));
  }
};

/** D3Q19 lattice PDF field */
using D3Q19Lattice = GridVector<double, D3Q19::Q>;

/** D3Q19 lattice PDF field view */
using D3Q19LatticeView = GridVectorView<double, D3Q19::Q>;

/** D3Q19 lattice PDF field const view */
using D3Q19LatticeConstView = GridVectorView<const double, D3Q19::Q>;

} // namespace oktal::lbm

#endif // OKTAL_LBM_D3Q19_HPP
