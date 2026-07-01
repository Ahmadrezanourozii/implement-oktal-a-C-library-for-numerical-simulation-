#pragma once

#include <cstddef>
#include <functional>
#include <numeric>
#include <vector>

namespace HighFive {

inline size_t compute_total_size(const std::vector<size_t> &dims) {
  return std::accumulate(dims.begin(), dims.end(), size_t{1u},
                         std::multiplies<size_t>());
}

} // namespace HighFive
