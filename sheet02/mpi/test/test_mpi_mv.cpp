#include <mpi.h>

#include <cassert>
#include <fmt/printf.hpp>
#include <matvec/copy.hpp>
#include <matvec/mv.hpp>
#include <mpi/mv.hpp>
#include <string>
#include <tools/init.hpp>
#include <tools/print.hpp>

#ifndef DIM_M
#define DIM_M 9
#endif

#ifndef DIM_N
#define DIM_N 9
#endif

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    int nof_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    std::size_t m = DIM_M;
    std::size_t n = DIM_N;
    std::ptrdiff_t incRow = n;
    double* A = nullptr;
    double x[n];
    double y_ref[m], y_tst[m];
    double ab[2];  // alpha and beta

    if (rank == 0) {
        A = new double[m * n];
        lpc::tools::random_init(m * n, A, 1);
        lpc::tools::random_init(n, x, 1);
        lpc::tools::random_init(m, y_ref, 1);
        lpc::matvec::copy(m, y_ref, 1, y_tst, 1);
        lpc::tools::random_init(2, ab, 1);
    }

    lpc::mpi::mv(m, n, ab[0], A, incRow, x, 1, ab[1], y_tst, 1);

    if (rank == 0) {
        lpc::matvec::mv(m, n, ab[0], A, incRow, x, 1, ab[1], y_ref, 1);

        std::string line(32, '-');
        fmt::printf("%s\n", line);
        fmt::printf("[i] m=%zu, n=%zu, p=%zu\n", m, n, nof_processes);
        fmt::printf("[i] Testing equality of y_ref and y_tst...\n");

        fmt::printf("[i] (Showing <= 10 entries)\n");
        fmt::printf("[i] y_ref = ");
        lpc::tools::print_vector(m <= 10 ? m : 10, y_ref, 1);
        fmt::printf("[i] y_tst = ");
        lpc::tools::print_vector(m <= 10 ? m : 10, y_tst, 1);

        fmt::printf("[i] (No output = no error)\n");
        for (std::size_t i = 0; i < m; i++) {
            if (y_ref[i] != y_tst[i]) {
                fmt::printf("[e] y_ref[%zu] != y_tst[%zu]\n", i, i);
            }
        }

        fmt::printf("[i] Done.\n");
        fmt::printf("%s\n", line);

        delete[] A;
    }

    MPI_Finalize();
}
