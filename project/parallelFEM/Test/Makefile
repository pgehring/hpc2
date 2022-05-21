
CC = mpicc

ifdef OPENMP
  CFLAGS += -g -fopenmp 
  LDFLAGS += -g -fopenmp
endif

ifdef OPT
  CFLAGS += -g -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

DEBUG =

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(DEBUG) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib test_mesh_split test_vec_accumulate test_mesh_transfer

lib:
	( cd ../Lib ; $(MAKE) )

test_mesh_split: lib test_mesh_split.c Makefile
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)
	- ./$@

test_vec_accumulate: lib test_vec_accumulate.c Makefile
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)
	- mpirun -np 4 ./$@

test_mesh_transfer: lib test_mesh_transfer.c Makefile
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r test_mesh_split \
    *.a *.dSYM *.obj *.dll
