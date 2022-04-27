CXX := mpicc -std=c++17
FLAGS := -O3 -Wall -I .
HEADER := *.hpp

all: $(HEADER)

%.hpp: Makefile
	$(CXX) $(FLAGS) -c $@ -o ../build/$*.o 
