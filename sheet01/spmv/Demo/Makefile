CC = gcc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib hpc_demo bench
	- ./hpc_demo < ../Matrix/t1
	- ./hpc_demo < ../Matrix/t2
	- ./hpc_bench

lib:
	( cd ../Lib ; $(MAKE) )

bench: lib hpc_bench.c Makefile
	$(CC) $(CF) $(I) -o hpc_bench hpc_bench.c $(CS)

hpc_demo: lib hpc_demo.c Makefile
	$(CC) $(CF) $(I) -o hpc_demo hpc_demo.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo hpc_bench \
    *.a *.dSYM *.obj *.dll
