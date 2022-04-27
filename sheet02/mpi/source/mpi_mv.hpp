#ifndef MV_HPP
#define MV_HPP

#include <mpi.h>

#include <cassert>
#include <cstddef>

#include "fundamental.hpp"
#include "mpi_mv_div_p.hpp"

/**
 * matrix vector product y <- alpha*A*x + beta*y
 * @param A matrix stored in row-major
 */
template <typename ALPHA, typename TA, typename TX, typename BETA, typename TY>
void mpi_mv(std::size_t m, std::size_t n, ALPHA alpha, TA *A,
	    std::ptrdiff_t incRow, TX *x, std::ptrdiff_t incX, BETA beta, TY *y,
	    std::ptrdiff_t incY, std::size_t p) {
	if (m % p == 0) {
		/**FIXME mpi_mv_div_p*/
	} else {
		/**FIXME mpi_mv_not_div_p*/
	}
}

#endif
