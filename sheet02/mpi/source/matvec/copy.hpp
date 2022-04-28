#ifndef LPC_MATVEC_COPY_HPP
#define LPC_MATVEC_COPY_HPP

#include <cstddef>

namespace lpc {
namespace matvec {

/*
 * copies y <- x
 */
template <typename TX, typename TY>
void copy(std::size_t n,               //
          TX *x, std::ptrdiff_t incX,  //
          TY *y, std::ptrdiff_t incY) {
    for (std::size_t i = 0; i < n; i++) {
        y[i * incY] = x[i * incX];
    }
}

}  // namespace matvec
}  // namespace lpc

#endif
