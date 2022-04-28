CXX := mpic++ -std=c++17
FLAGS := -O3 -Wall -I ../source

.PHONY: tests

tests: test_mpi_mv_div_p test_mpi_mv

test_mpi_mv_div_p: test_mpi_mv_div_p.cpp
	$(foreach DIM, 10 100 1000, \
		$(CXX) $(FLAGS) -DDIM_M=$(DIM) -DDIM_N=$(DIM) -o ../build/$@ $@.cpp \
		&& mpirun -np 2 ../build/$@;)
	
test_mpi_mv: test_mpi_mv.cpp
	$(foreach DIM, 9 99 999, \
		$(CXX) $(FLAGS) -DDIM_M=$(DIM) -DDIM_N=$(DIM) -o ../build/$@ $@.cpp \
		&& mpirun -np 4 ../build/$@;)
