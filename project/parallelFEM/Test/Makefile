
CC = gcc

ifdef OPENMP
  CFLAGS += -g -fopenmp 
  LDFLAGS += -g -fopenmp
endif

ifdef OPT
  CFLAGS += -g -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib test_mesh_split
	- ./test_mesh_split

lib:
	( cd ../Lib ; $(MAKE) )

test_mesh_split: lib test_mesh_split.c Makefile
	$(CC) $(CF) $(I) -o $@ $@.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo1 hpc_demo2 hpc_demo3 hpc_demo4 hpc_demo5 \
    *.a *.dSYM *.obj *.dll
