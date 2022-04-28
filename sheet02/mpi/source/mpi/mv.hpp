#ifndef LPC_MPI_MV_HPP
#define LPC_MPI_MV_HPP

#include <mpi.h>

#include <aux/slices.hpp>
#include <cassert>
#include <cstddef>
#include <mpi/fundamental.hpp>

namespace lpc {
namespace mpi {

constexpr int root = 0;

/**
 * matrix vector product y <- alpha*A*x + beta*y
 * @param A if node is root, A is a matrix stored in row-major;
 *          otherwise it should be a nullptr
 * @param x vector of length n
 * @param y vector of length m
 */
template <typename ALPHA, typename TA, typename TX, typename BETA, typename TY>
void mv(std::size_t m, std::size_t n, ALPHA alpha, TA *A, std::ptrdiff_t incRow,
        TX *x, std::ptrdiff_t incX, BETA beta, TY *y, std::ptrdiff_t incY) {
    using namespace hpc;
    using namespace hpc::aux;
    using namespace hpc::mpi;

    // gather mpi infos
    int rank;
    int p;  // = nof_processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // compute variables for scatter
    UniformSlices<std::size_t> slices(p, m);
    std::size_t m_ = slices.size(rank);
    std::size_t offset = slices.offset(rank);
    TA A_[m_ * n];
    int sendcounts[p];
    int displs[p];
    for (int i = 0; i < p; i++) {
        sendcounts[i] = slices.size(i) * n;
        displs[i] = slices.offset(i) * n;
    }
    int recvcount = m_ * n;

    // scatter A, fully send x and y
    MPI_Scatterv(A, sendcounts, displs, get_type(TA{0}), /* send info */
                 A_, recvcount, get_type(TA{0}),         /* recv info */
                 root, MPI_COMM_WORLD);
    MPI_Bcast(x, n, get_type(TX{0}), root, MPI_COMM_WORLD);
    MPI_Bcast(y, m, get_type(TY{0}), root, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, get_type(alpha), root, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, get_type(beta), root, MPI_COMM_WORLD);

    // compute this nodes part of y
    lpc::matvec::mv(m_, n,      /**/
                    alpha,      /**/
                    A_, incRow, /**/
                    x, incX,    /**/
                    beta,       /**/
                    &y[offset], incY);

    // compute variables for gather
    int recvcounts[p];
    for (int i = 0; i < p; i++) {
        recvcounts[i] = slices.size(i);
        displs[i] = slices.offset(i);
    }

    // gather y
    MPI_Gatherv(rank ? &y[offset] : MPI_IN_PLACE, m_,
                get_type(TY{0}),                        /* send info */
                y, recvcounts, displs, get_type(TY{0}), /* recv info */
                root, MPI_COMM_WORLD);
}

/**
 * matrix vector product y <- alpha*A*x + beta*y
 * @param A if node is root, A is a matrix stored in row-major;
 *          otherwise it should be a nullptr
 * @param x vector of length n
 * @param y vector of length m
 */
template <typename ALPHA, typename TA, typename TX, typename BETA, typename TY>
void mv_div_p(std::size_t m, std::size_t n, ALPHA alpha, TA *A,
              std::ptrdiff_t incRow, TX *x, std::ptrdiff_t incX, BETA beta,
              TY *y, std::ptrdiff_t incY) {
    using namespace hpc;
    using namespace hpc::mpi;

    // gather mpi infos
    int root = 0;
    int rank;
    int p;  // = nof_processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    assert(m % p == 0);

    // compute variables for scatter
    std::size_t m_ = m / p;
    TA A_[m_ * n];
    std::size_t sendcount = m_ * n;

    // scatter A, fully send x and y
    MPI_Scatter(A, sendcount, get_type(TA{0}),  /* send info */
                A_, sendcount, get_type(TA{0}), /* recv info */
                root, MPI_COMM_WORLD);
    MPI_Bcast(x, n, get_type(TX{0}), root, MPI_COMM_WORLD);
    MPI_Bcast(y, m, get_type(TY{0}), root, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, get_type(alpha), root, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, get_type(beta), root, MPI_COMM_WORLD);

    // compute this nodes part of y
    lpc::matvec::mv(m_, n,      /**/
                    alpha,      /**/
                    A_, incRow, /**/
                    x, incX,    /**/
                    beta,       /**/
                    &y[m_ * rank], incY);

    // gather y
    MPI_Gather(rank ? &y[m_ * rank] : MPI_IN_PLACE, m_,
               get_type(TY{0}),        /* send info */
               y, m_, get_type(TY{0}), /* recv info */
               root, MPI_COMM_WORLD);
}

}  // namespace mpi
}  // namespace lpc

#endif
