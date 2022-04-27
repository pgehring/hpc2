#ifndef MV_HPP
#define MV_HPP

#include <cstddef>

namespace lpc {

/**
 * matrix vector product y <- alpha*A*x + beta*y
 * @param A matrix stored in row-major
 */
template <typename ALPHA, typename TA, typename TX, typename BETA, typename TY>
void mv(std::size_t m, std::size_t n, ALPHA alpha, TA *A,
	    std::ptrdiff_t incRow, TX *x, std::ptrdiff_t incX, BETA beta, TY *y,
	    std::ptrdiff_t incY) {
	if (m == 0) return;
	if (beta == BETA(0)) {
		for (std::size_t i = 0; i < m; i++) {
			y[i*incY] = TY(0);
		}
	} else if (beta != BETA(1)) {
		for (std::size_t i = 0; i < m; i++) {
			y[i*incY] *= beta;
		}
	}
	if (n == 0 || alpha == ALPHA(0)) return;
	for (std::size_t i = 0; i < m; i++) {
		for (std::size_t j = 0; j < n; j++) {
			y[i*incY] += A[i*incRow + j] * x[j*incX];
		}
	}
}

} // namespace lpc

#endif 
