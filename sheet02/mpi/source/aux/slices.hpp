#ifndef HPC_AUX_SLICE_HPP
#define HPC_AUX_SLICE_HPP

#include <cassert>
#include <type_traits>

/*
   UniformSlices and AlignedSlices support partitioning where
    - N, i.e. the problem size (e.g. a dimension of a matrix or vector), and
    - P, i.e. the number of partitions
   are fixed. Both classes deliver then
    - Oi, i.e. the offset for the i-th partition, and
    - Li, i.e. the length of the i-th partition
   for i = 0, ..., P-1 where
    - sum Li = P
    - Oi = sum Lj for j = 0, ..., i-1, and
    - Li <= Lj if i <= j
*/

namespace hpc {
namespace aux {

std::size_t align(std::size_t size, std::size_t alignment) {
    std::size_t remainder = size % alignment;
    if (remainder == 0) {
        return size;
    } else {
        return size + alignment - remainder;
    }
}

template <typename T = std::size_t>
struct UniformSlices {
    UniformSlices(T nof_partitions, T problem_size)
        : nof_partitions((assert(nof_partitions > 0), nof_partitions)),
          problem_size(problem_size),
          remainder(problem_size % nof_partitions),
          slice_size(problem_size / nof_partitions) {}
    T offset(T index) const {
        assert(index < nof_partitions);
        if (index < remainder) {
            return index * (slice_size + 1);
        } else {
            return remainder * (slice_size + 1) +
                   (index - remainder) * slice_size;
        }
    }
    T size(T index) const {
        assert(index < nof_partitions);
        if (index < remainder) {
            return slice_size + 1;
        } else {
            return slice_size;
        }
    }
    T nof_partitions;
    T problem_size;
    T remainder;
    T slice_size;
};

struct AlignedSlices {
    AlignedSlices(std::size_t nof_partitions, std::size_t problem_size,
                  std::size_t granularity = 1)
        : nof_partitions((assert(nof_partitions > 0), nof_partitions)),
          problem_size(problem_size),
          slice_size{
              problem_size / granularity >= nof_partitions
                  ? align((problem_size + nof_partitions - 1) / nof_partitions,
                          granularity)
                  : granularity},
          full_slices(problem_size / slice_size) {}
    std::size_t offset(std::size_t index) const {
        assert(index < nof_partitions);
        if (slice_size > 0) {
            return index * slice_size;
        } else {
            return index;
        }
    }
    std::size_t size(std::size_t index) const {
        assert(index < nof_partitions);
        if (index < full_slices) {
            return slice_size;
        } else if (index == full_slices) {
            return problem_size - full_slices * slice_size;
        } else {
            return 0;
        }
    }
    std::size_t nof_partitions;
    std::size_t problem_size;
    std::size_t slice_size;
    std::size_t full_slices;
};

template <typename Slices, typename Body>
void foreach_slice(std::size_t nof_partitions, std::size_t problem_size,
                   Body body) {
    Slices slices{nof_partitions, problem_size};
    for (std::size_t index = 0; index < nof_partitions; ++index) {
        std::size_t size = slices.size(index);
        if (size > 0) {
            body(slices.offset(index), size);
        }
    }
}

template <typename Slices, typename Body>
void foreach_slice(std::size_t nof_partitions, std::size_t problem_size,
                   std::size_t granularity, Body body) {
    Slices slices{nof_partitions, problem_size, granularity};
    for (std::size_t index = 0; index < nof_partitions; ++index) {
        std::size_t size = slices.size(index);
        if (size > 0) {
            body(slices.offset(index), size);
        }
    }
}

}  // namespace aux
}  // namespace hpc

#endif  // HPC_AUX_SLICE_H

