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

all: test_cs_spmv test_jds_spmv test_jds_spmv_diag test_sky_spmv_diag
	- ./test_cs_spmv
	- ./test_jds_spmv
	- ./test_jds_spmv_diag
	- ./test_sky_spmv_diag

test_cs_spmv: test_cs_spmv.c
	$(CC) $(CF) $(I) -o $@ $@.c utils.c $(CS)

test_jds_spmv: test_jds_spmv.c
	$(CC) $(CF) $(I) -o $@ $@.c utils.c $(CS)

test_jds_spmv_diag: test_jds_spmv_diag.c
	$(CC) $(CF) $(I) -o $@ $@.c utils.c $(CS)

test_sky_spmv_diag: test_sky_spmv_diag.c
	$(CC) $(CF) $(I) -o $@ $@.c utils.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r ./test_cs_spmv \
    *.a *.dSYM *.obj *.dll
	- $(RM) -r ./test_jds_spmv \
    *.a *.dSYM *.obj *.dll
	- $(RM) -r ./test_jds_spmv_diag \
    *.a *.dSYM *.obj *.dll
	- $(RM) -r ./test_sky_spmv_diag \
    *.a *.dSYM *.obj *.dll