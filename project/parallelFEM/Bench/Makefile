
CC = mpicc

ifdef OPENMP
  CFLAGS += -g -wall -fopenmp 
  LDFLAGS += -g -fopenmp
endif

ifdef OPT
  CFLAGS += -g -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

DEBUG =

RESULT_DIR=bench_results

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(DEBUG) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib bench_ref_mg bench_cg_sep

lib:
	( cd ../Lib ; $(MAKE) )

bench_ref_mg: lib bench_ref_mg.c Makefile $(RESULT_DIR)
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)

bench_cg_sep:   lib bench_cg_sep.c Makefile $(RESULT_DIR)
	 - $(CC) $(CF) $(I) -o $@ $@.c $(CS)

$(RESULT_DIR):
	mkdir $@

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r bench_ref_MG bench_cg_sep
	- $(RM) *.a *.dSYM *.obj *.dll
	- $(RM) *.o* *.e* *.p* core*
	- $(RM) bench_results/*.csv
	     

