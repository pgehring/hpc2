#ifndef LPC_MATVEC_DOT_HPP
#define LPC_MATVEC_DOT_HPP

#include <cstddef>
#include <type_traits>

namespace lpc {
namespace matvec {

/**
 * computes the dot product <x, y>
 */
template <typename TX, typename TY, typename T>
T dot(std::size_t n, TX *x, std::ptrdiff_t incX, TY *y, std::ptrdiff_t incY) {
    if (n == 0) return T(0);
    T val{0};
    for (std::size_t i = 0; i < n; i++) {
        val += x[i * incX] * y[i * incY];
    }
    return val;
}

}  // namespace matvec
}  // namespace lpc

#endif
