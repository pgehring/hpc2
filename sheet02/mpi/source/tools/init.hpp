#ifndef INIT_HPP
#define INIT_HPP

#include <cstddef>
#include <random>

namespace lpc {
namespace tools {

/**
 * initializes a vector x with 0, 1, 2, ...
 */
template <typename T>
void incremental_init(std::size_t n, T* x, std::ptrdiff_t inc) {
    for (std::size_t i = 0; i < n; i++) {
        x[i * inc] = i;
    }
}

std::mt19937 engine{std::random_device()()};
std::uniform_real_distribution<double> uniform(-1, 1);

/**
 * initializes a vector x with random values in [-1, 1]
 */
template <typename T>
void random_init(std::size_t n, T* x, std::ptrdiff_t inc) {
    for (std::size_t i = 0; i < n; i++) {
        x[i * inc] = uniform(engine);
    }
}

}  // namespace tools
}  // namespace lpc

#endif
