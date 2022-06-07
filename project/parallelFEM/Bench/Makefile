

CC = mpicc
CFLAGS = -g
LDFLAGS = -g
ifdef OPENMP
  CFLAGS += -wall -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

ifdef _DEBUG
DEBUG = -D_DEBUG
else
DEBUG=
endif

RESULT_DIR=bench_results

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(DEBUG) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib bench_ref_mg bench_cg_sep bench_cg_full bench_jcb

lib:
	( cd ../Lib ; $(MAKE) _DEBUG=1 )

bench_ref_mg: lib bench_ref_mg.c Makefile $(RESULT_DIR)
	- $(CC) $(CF) $(I) -o $@ $@.c $(CS)

bench_cg_sep:   lib bench_cg_sep.c Makefile $(RESULT_DIR)
	 - $(CC) $(CF) $(I) -o $@ $@.c $(CS)

bench_cg_full:	lib bench_cg_full.c Makefile $(RESULT_DIR)
	 - $(CC) $(CF) $(I) -o $@ $@.c $(CS)

bench_jcb: lib bench_jcb.c Makefile $(RESULT_DIR)
	 - $(CC) $(CF) $(I) -o $@ $@.c $(CS)

$(RESULT_DIR):
	mkdir $@

clean:
	- $(RM) *.o
	- $(RM) *.a *.dSYM *.obj *.dll
	- $(RM) *.o* *.e* *.p* core* .wrapper*
	- (cd jobs_bench_cg_sep; $(RM) *.o* *.e* *.p* core*, .wrapper*)
	- (cd jobs_bench_cg_full; $(RM) *.o* *.e* *.p* core*, .wrapper*)

purge: clean
	- $(RM) -r bench_ref_mg bench_cg_sep
	- $(RM) bench_results/*.csv
	     

