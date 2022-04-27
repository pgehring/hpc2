#ifndef PRINT_HPP
#define PRINT_HPP

#include <cstddef>

#include <fmt/printf.hpp>

namespace lpc {
namespace tools {

/**
 * prints a vector
 */
template<typename T>
void print_vector(std::size_t n, T* x, std::ptrdiff_t inc) {
    for (std::size_t i = 0; i < n; i++) {
        fmt::printf(" %7.3lf", x[i * inc]);
    }
    fmt::printf("\n");
}

/**
 * prints a matrix
 * @param A matrix stored in row-major
 */
template<typename T>
void print_matrix(std::size_t m, std::size_t n, T *A, std::ptrdiff_t incRow) {
	for (std::size_t i = 0; i < m; i++) {
		for (std::size_t j = 0; j < n; j++) {
			fmt::printf(" %7.3lf", A[i*incRow + j]);
		}
		fmt::printf("\n");
	}	
}

} // namespace tools
} // namespace lpc

#endif
