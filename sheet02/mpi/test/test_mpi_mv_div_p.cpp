#include <mpi.h>

#include <init.hpp>
#include <mpi_mv_div_p.hpp>
#include <print.hpp>

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	std::size_t m = 10;
	std::size_t n = 10;
	std::ptrdiff_t incRow = n;
	double A[m * n];
	double x[n];
	double y[m];

	lpc::incremental_init(m * n, A, 1);
	lpc::mpi_mv_div_p(m, n, 1, A, incRow, x, 1, 1, y, 1);

	MPI_Finalize();
}
