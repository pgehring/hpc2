CC = gcc

LIBRARY = libhpc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG -Wpedantic  # -mfpmath=sse #-march=native #
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fPIC
I = -I../Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

HPC = \
 cs_compress.o \
 cs_entry.o \
 cs_spmv.o \
 cs_load.o \
 cs_print.o \
 cs_util.o \
 gem_compress.o \
 gem_gauss.o \
 gem_gausssol.o \
 gem_print.o \
 gem_util.o \
 gem_spmv.o \
 hpc_cumsum.o \
 hpc_malloc.o \
 hpc_util.o \
 jds_compress.o \
 jds_spmv.o \
 jds_print.o \
 jds_util.o  \
 sed_compress.o \
 sed_spmv.o \
 sed_print.o \
 sed_util.o  \
 sky_compress.o \
 sky_print.o  \
 sky_util.o \
 sky_load.o \
 sky_spmv.o

$(HPC): ../Include/hpc.h Makefile

%.o: ../Source/%.c ../Include/hpc.h
	$(CC) $(CF) $(I) -c $<

static: $(AR_TARGET)

$(AR_TARGET): $(HPC)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install archive in this directory
install: $(AR_TARGET)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

