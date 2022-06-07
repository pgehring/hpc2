
CC = mpicc

ifdef OPENMP
  CFLAGS += -g -wall -fopenmp 
  LDFLAGS += -g -fopenmp
endif

ifdef OPT
  CFLAGS += -g -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

ifdef _DEBUG
DEBUG = -D_DEBUG
else
DEBUG=
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(DEBUG) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

HPC_TARGETS = \
	test_mesh_split \
	test_vec_accumulate \
	test_mesh_transfer \
	test_build_sm \
	test_sed_spmv \
	test_dot_dist \
	test_solver

all: lib $(HPC_TARGETS)

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
	- mpirun -np 8 ./$@

test_build_sm: lib test_build_sm.c Makefile
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)

test_sed_spmv: lib test_sed_spmv.c Makefile
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)
	- ./$@

test_dot_dist: lib test_dot_dist.c Makefile
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)
	- mpirun -np 4 ./$@

test_solver: lib test_solver.c Makefile
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)
	- mpirun -np 8 ./$@

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r $(HPC_TARGETS) \
    *.a *.dSYM *.obj *.dll
