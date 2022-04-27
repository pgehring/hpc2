#ifndef PRINT_HPP
#define PRINT_HPP

#include <cstddef>

#include "printf.hpp"

namespace lpc {

template<typename T>
void print_matrix(std::size_t m, std::size_t n, T *A, std::ptrdiff_t incRow) {
	for (std::size_t i = 0; i < m; i++) {
		for (std::size_t j = 0; j < n; j++) {
			fmt::printf(" %3.1lf", A[i*incRow + j]);
		}
		fmt::printf("\n");
	}	
}

}

#endif
