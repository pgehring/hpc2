#ifndef HPC_MPI_FUNDAMENTAL_HPP
#define HPC_MPI_FUNDAMENTAL_HPP 1

#include <mpi.h>

#include <complex>

namespace hpc {
namespace mpi {

template <typename T>
struct fundamental_type {
    static constexpr bool defined = false;
};

template <>
struct fundamental_type<char> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_CHAR; }
};
template <>
struct fundamental_type<signed char> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_SIGNED_CHAR; }
};
template <>
struct fundamental_type<unsigned char> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_UNSIGNED_CHAR; }
};
template <>
struct fundamental_type<short> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_SHORT; }
};
template <>
struct fundamental_type<unsigned short> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_UNSIGNED_SHORT; }
};
template <>
struct fundamental_type<int> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_INT; }
};
template <>
struct fundamental_type<unsigned> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_UNSIGNED; }
};
template <>
struct fundamental_type<long> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_LONG; }
};
template <>
struct fundamental_type<unsigned long> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_UNSIGNED_LONG; }
};
template <>
struct fundamental_type<long long> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_LONG_LONG; }
};
template <>
struct fundamental_type<unsigned long long> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_UNSIGNED_LONG_LONG; }
};
template <>
struct fundamental_type<float> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_FLOAT; }
};
template <>
struct fundamental_type<double> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_DOUBLE; }
};
template <>
struct fundamental_type<long double> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_LONG_DOUBLE; }
};

/* the following types are not supported by older
   MPI implementations; hence we make it dependent
   on the corresponding preprocessor symbols */
#ifdef MPI_CXX_FLOAT_COMPLEX
template <>
struct fundamental_type<std::complex<float>> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_CXX_FLOAT_COMPLEX; }
};
#endif
#ifdef MPI_CXX_DOUBLE_COMPLEX
template <>
struct fundamental_type<std::complex<double>> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_CXX_DOUBLE_COMPLEX; }
};
#endif
#ifdef MPI_CXX_LONG_DOUBLE_COMPLEX
template <>
struct fundamental_type<std::complex<long double>> {
    static constexpr bool defined = true;
    MPI_Datatype get() const { return MPI_CXX_LONG_DOUBLE_COMPLEX; }
};
#endif

template <typename T>
typename std::enable_if<fundamental_type<T>::defined, MPI_Datatype>::type
get_type(const T& value) {
    return fundamental_type<T>().get();
}

}  // namespace mpi
}  // namespace hpc

#endif
