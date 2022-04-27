#ifndef MPI_MV_DIV_P
#define MPI_MV_DIV_P

#include <mpi.h>

#include <cassert>
#include <cstddef>
#include <mpi/fundamental.hpp>
//#include <tools/print.hpp>

namespace lpc {
namespace mpi {

/**
 * matrix vector product y <- alpha*A*x + beta*y
 * @param A matrix stored in row-major
 */
template <typename ALPHA, typename TA, typename TX, typename BETA, typename TY>
void mv_div_p(std::size_t m, std::size_t n, ALPHA alpha, TA *A,
              std::ptrdiff_t incRow, TX *x, std::ptrdiff_t incX, BETA beta,
              TY *y, std::ptrdiff_t incY) {
    using namespace hpc;
    using namespace hpc::mpi;

    int root = 0;
    int rank;
    int p;  // = nof_processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    assert(m % p == 0);

    std::size_t m_ = m / p;
    TA *A_ = new TA[m_ * n];
    std::size_t sendcount = m_ * n;

    MPI_Scatter(A, sendcount,
                A ? get_type(A[0]) : MPI_DATATYPE_NULL, /* send info */
                A_, sendcount, get_type(A_[0]),         /* recv info */
                root, MPI_COMM_WORLD);
    MPI_Bcast(x, n, get_type(x[0]), root, MPI_COMM_WORLD);
    MPI_Bcast(y, m, get_type(y[0]), root, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, get_type(alpha), root, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, get_type(beta), root, MPI_COMM_WORLD);

    /*
    if (rank) {
        fmt::printf("A_ =\n");
        lpc::tools::print_matrix(m_, n, A_, incRow);
        fmt::printf("x = ");
        lpc::tools::print_vector(n, x, 1);
        fmt::printf("y = ");
        lpc::tools::print_vector(m, y, 1);
    }
    */

    lpc::matvec::mv(m_, n,      /**/
                    alpha,      /**/
                    A_, incRow, /**/
                    x, incX,    /**/
                    beta,       /**/
                    &y[m_ * rank], incY);

    /*
    if (rank) {
        fmt::printf("y = ");
        lpc::tools::print_vector(m, y, 1);
    }
    */

    MPI_Gather(rank ? &y[m_ * rank] : MPI_IN_PLACE, m_,
               get_type(y[0]),        /* send info */
               y, m_, get_type(y[0]), /* recv info */
               root, MPI_COMM_WORLD);

    delete[] A_;
}

}  // namespace mpi
}  // namespace lpc

#endif
