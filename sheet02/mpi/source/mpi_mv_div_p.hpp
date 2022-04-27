#ifndef MPI_MV_DIV_P
#define MPI_MV_DIV_P

#include <mpi.h>

#include <cassert>
#include <cstddef>

#include "fundamental.hpp"
#include "print.hpp"

namespace lpc {

template <typename ALPHA, typename TA, typename TX, typename BETA, typename TY>
void mpi_mv_div_p(std::size_t m, std::size_t n, ALPHA alpha, TA *A,
		  std::ptrdiff_t incRow, TX *x, std::ptrdiff_t incX, BETA beta,
		  TY *y, std::ptrdiff_t incY) {
	using namespace hpc;
	using namespace hpc::mpi;
	
	int rank;
	int nof_processes;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
	assert(m % nof_processes == 0);

	std::size_t m_ = m / nof_processes;
	TA *A_ = new TA[m_ * n];
	std::size_t sendcount = m_ * n;

	if (rank == 0) {
		MPI_Scatter(A, sendcount, get_type(A[0]), A_, sendcount,
			    get_type(A_[0]), 0, MPI_COMM_WORLD);
	} else {
		MPI_Scatter(nullptr, 0, MPI_DATATYPE_NULL, A_, sendcount,
			    get_type(A_[0]), 0, MPI_COMM_WORLD);
	}

	delete[] A_;
}

} // namespace lpc
#endif
