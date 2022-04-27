CXX := mpic++ -std=c++17
FLAGS := -O3 -Wall -I ../source -xc++ -lstdc++ -shared-libgcc

all: test_mpi_mv_div_p

test_mpi_mv_div_p: test_mpi_mv_div_p.cpp
	- $(CXX) $(FLAGS) -o ../build/$@ $@.cpp
	- mpirun -np 2 ../build/$@
