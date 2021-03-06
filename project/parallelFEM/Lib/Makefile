# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
# resp.  
#       make VERBOSE=1
#       make CC=gcc OPENMP=1
#       make CC=/opt/local/bin/clang-mp-3.7 OPT=1 OPENMP=1
#
# To run a demo using the library 
#       cd ../Demo ; make
#

CC = mpicc
CFLAGS += -g
LIBRARY = libhpc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG -Wpedantic  # -mfpmath=sse #-march=native #
endif

ifdef DEBUG
DEBUG = -D_DEBUG
else
DEBUG =
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(DEBUG) -fPIC
I = -I../Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

HPC = \
 sed_buildS.o \
 sed_dupl.o \
 sed_gaxpy.o \
 sed_print.o \
 sed_gs_constr.o \
 sed_util.o  \
 sed_sm_build.o \
 sed_spmv_sym.o \
 cs_print.o \
 cs_util.o \
 cs_spmv.o \
 cs_entry.o \
 hpc_cumsum.o \
 hpc_malloc.o \
 hpc_fmg.o  \
 hpc_mg.o  \
 hpc_mg_cycle.o  \
 hpc_solve_cg.o \
 hpc_solve_util.o \
 hpc_jacobi.o \
 mesh_buildRhs.o \
 mesh_load.o  \
 hpc_prol.o  \
 hpc_prol_quad.o  \
 mesh_getEdge2no.o \
 hpc_rest.o  \
 mesh_print.o  \
 mesh_refine.o  \
 mesh_getFixed.o  \
 mesh_util.o \
 mesh_split.o \
 mesh_initRefinement.o \
 mesh_build_rhs.o \
 slices.o \
 mesh_transfer.o \
 vec_accumulate.o \
 vec_util.o \
 vec_dot_dist.o \
 blas_level1.o


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

