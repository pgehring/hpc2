
CC = mpicc

ifdef OPENMP
  CFLAGS += -g -wall -fopenmp 
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

all: lib bench_ref_mg

lib:
	( cd ../Lib ; $(MAKE) )

bench_ref_mg: lib bench_ref_mg.c Makefile
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r bench_ref_MG\
    *.a *.dSYM *.obj *.dll
