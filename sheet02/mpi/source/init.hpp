#ifndef INIT_HPP
#define INIT_HPP

#include <cstddef>

namespace lpc {

template <typename T>
void incremental_init(std::size_t n, T *x, std::ptrdiff_t inc) {
	for (std::size_t i = 0; i < n; i++) {
		x[i * inc] = i;
	}
}

}  // namespace lpc

#endif
