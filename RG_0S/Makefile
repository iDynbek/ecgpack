# This is a makefile to compile an ECG code. Examples of compilations with this
# makefile:
#    make release COMPILER_TYPE=gnu MACHINE=ubuntu-generic PREC=8 OPTLPKBLS=yes EXEFILE=ecg
#    make debug COMPILER_TYPE=gnu MACHINE=shabyt PREC=10 OPTLPKBLS=no EXEFILE=ecg
#    make release COMPILER_TYPE=intel MACHINE=muon PREC=8 OPTLPKBLS=yes EXEFILE=ecg


# Set default compiler type (intel, gnu, sun, etc..)
COMPILER_TYPE=gnu

# Set default machine where the compilation takes place (mylinuxbox, ice, marin, etc.)
MACHINE=ubuntu-generic

# Set the default precision (kind of real type, 8, 10, or 16)
# Value PREC=8 corresponds to standard double precision
# Value PREC=16 corresponds to quadruple precision
PREC=8

# Set whether selected linear algebra routines should be taken from optimized
# LAPACK/BLAS library (OPTLPKBLS=yes) such as Intel Math Kernel library or Sun Performance
# Library. Otherwise (OPTLPKBLS=no or any other value) the supplied nonoptimized versions
# of these routines will be used.
OPTLPKBLS=yes

# CUDA backend (USE_CUDA=yes requires nvcc, PREC=8, and an NVIDIA GPU).
# Set CUDA_ARCH to match your GPU's compute capability (sm_86 = RTX 3060 Ti / A100).
# Set CUDA_HOME to the CUDA toolkit installation prefix.
USE_CUDA=no
NVCC=nvcc
CUDA_ARCH=sm_86
CUDA_HOME=/usr/local/cuda

# OpenACC backend (USE_OPENACC=yes requires the NVIDIA HPC SDK / nvfortran and an
# NVIDIA GPU). Unlike CUDA, the ENTIRE program is compiled with nvfortran -acc; the
# device routine lives in matelem_acc.f90. GPU/PREC=8 only. Build with OPTLPKBLS=no so
# the bundled reference BLAS/LAPACK are compiled by nvfortran (avoids linking a system
# optimized library against nvfortran). Set OPENACC_ARCH to your GPU's compute
# capability (cc70 = V100, cc80 = A100, cc86 = RTX 3060 Ti). Requires the nvhpc module
# to be loaded so mpif90 resolves to the HPC SDK wrapper around nvfortran.
USE_OPENACC=no
OPENACC_FC=mpif90
OPENACC_ARCH=cc70

# Set the default name for the generated executable file
EXEFILE=ecg

# Directory where all source files are located
SRC_DIR=./src

# Set compiler flags for each machine and compiler type (for both release and debug builds).

# Generic Ubuntu machine, generic Linux machine, shabyt, muon
ifeq ($(shell [ "$(MACHINE)" = "ubuntu-generic" ] || [ "$(MACHINE)" = "linux-generic" ] || [ "$(MACHINE)" = "shabyt" ] || [ "$(MACHINE)" = "muon" ] && echo true),true)
  # gfortran and Open MPI (system defaults)
  ifeq ($(COMPILER_TYPE), gnu)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-J$(CONFIG) -cpp
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -fbounds-check -fbacktrace -finit-real=nan -fcheck=all -fimplicit-none -ffpe-trap=invalid,zero,overflow
      export LINKFLAGS+=-C -O0 -g -fbounds-check -fbacktrace -finit-real=nan -fcheck=all -fimplicit-none -ffpe-trap=invalid,zero,overflow
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -march=native
      export LINKFLAGS+=-O3 -march=native
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with Lapack and BLAS
        export LINKFLAGS+=-lblas -llapack
      endif
    endif
    ifeq ($(PREC), 10)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
    ifeq ($(PREC), 16)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # Intel Fortran and Intel MPI
  ifeq ($(COMPILER_TYPE), intel)
    FC=mpiifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -check bounds -check all -check uninit -traceback -ftrapuv -debug full
      export LINKFLAGS+=-C -O0 -g -check bounds -check all -check uninit -traceback -ftrapuv -debug full
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ip -xhost -fp-model precise
      export LINKFLAGS+=-O3 -ip -xhost -fp-model precise
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export FCFLAGS+=-I$(MKLROOT)/include
        export LINKFLAGS+=-L$(MKLROOT)
        export LINKFLAGS+=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
      endif
    endif
    ifeq ($(PREC), 16)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif  
endif

# puma, ocelote, elgato
ifeq ($(shell [ "$(MACHINE)" = "puma" ] || [ "$(MACHINE)" = "ocelote" ] || [ "$(MACHINE)" = "elgato" ] && echo true),true)
  # gfortran and Open MPI
  ifeq ($(COMPILER_TYPE), gnu)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-J$(CONFIG) -cpp
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -fbounds-check -fbacktrace -finit-real=nan -fcheck=all -fimplicit-none -ffpe-trap=invalid,zero,overflow
      export LINKFLAGS+=-C -O0 -g -fbounds-check -fbacktrace -finit-real=nan -fcheck=all -fimplicit-none -ffpe-trap=invalid,zero,overflow
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -march=native
      export LINKFLAGS+=-O3 -march=native
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with Lapack and BLAS
        export LINKFLAGS+=-lblas -llapack
      endif
    endif
    ifeq ($(PREC), 10)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
    ifeq ($(PREC), 16)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # Intel Fortran and Intel MPI
  ifeq ($(COMPILER_TYPE), intel)
    FC=mpiifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -check bounds -check all -check uninit -traceback -ftrapuv -debug full
      export LINKFLAGS+=-C -O0 -g -check bounds -check all -check uninit -traceback -ftrapuv -debug full
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ip -xhost -fp-model precise
      export LINKFLAGS+=-O3 -ip -xhost -fp-model precise
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export FCFLAGS+=-I$(MKLROOT)/include
        export LINKFLAGS+=-L$(MKLROOT)
        export LINKFLAGS+=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
      endif
    endif
    ifeq ($(PREC), 16)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif  
endif

# Creating lists of objects for convenience

ifeq ($(OPTLPKBLS), yes)
  export BARE_OBJS_LPKBLS=
else
  export BARE_OBJS_LPKBLS=BLAS.o LAPACK.o
endif
export OBJS_LPKBLS=$(addprefix $(CONFIG)/,$(BARE_OBJS_LPKBLS))

export BARE_OBJS=wp_def_$(PREC).o X1MACH.o dmng.o globvars.o misc.o linalg.o \
              matelem.o matform.o workproc.o main.o spin.o
export OBJS=$(addprefix $(CONFIG)/,$(BARE_OBJS))
export OBJS+=$(OBJS_LPKBLS)

# CUDA object files and flags (only when USE_CUDA=yes and PREC=8)
export CUDA_OBJS=
GPU_DEP=
ifeq ($(USE_CUDA), yes)
  ifeq ($(PREC), 8)
    export FCFLAGS+=-DUSE_CUDA
    export LINKFLAGS+=-L$(CUDA_HOME)/lib64 -lcudart -lstdc++
    export CUDA_OBJS=$(CONFIG)/matelem_cuda.o
    export BARE_OBJS+=matelem_gpu.o
    export OBJS=$(addprefix $(CONFIG)/,$(BARE_OBJS)) $(OBJS_LPKBLS)
    GPU_DEP=$(CONFIG)/matelem_gpu.o
  else
    $(error USE_CUDA=yes is only supported with PREC=8)
  endif
endif

# OpenACC backend (only when USE_OPENACC=yes and PREC=8). The whole program is
# rebuilt with nvfortran -acc, so FC and the compile/link flags are OVERRIDDEN here
# (replacing the gfortran/-march=native settings from the machine block above).
ACC_DEP=
ifeq ($(USE_OPENACC), yes)
  ifeq ($(PREC), 8)
    FC=$(OPENACC_FC)
    export FCFLAGS=-module $(CONFIG) -Mpreprocess -DUSE_OPENACC -acc -gpu=$(OPENACC_ARCH) -Kieee
    export LINKFLAGS=-module $(CONFIG) -acc -gpu=$(OPENACC_ARCH) -Kieee
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3
      export LINKFLAGS+=-O3
    endif
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-O0 -g -Mbounds -traceback
      export LINKFLAGS+=-O0 -g -Mbounds -traceback
    endif
    export BARE_OBJS+=matelem_acc.o
    export OBJS=$(addprefix $(CONFIG)/,$(BARE_OBJS)) $(OBJS_LPKBLS)
    ACC_DEP=$(CONFIG)/matelem_acc.o
  else
    $(error USE_OPENACC=yes is only supported with PREC=8)
  endif
endif


.PHONY: release debug all clean cleaner cleanest cleandebug cleanrelease ptarg

release:
ifndef CONFIG
	-@mkdir -p $(@)
	@$(MAKE) $(MAKEFILE) --no-print-directory $(@) CONFIG=$(@)
else
	@$(MAKE) $(MAKEFILE) --no-print-directory $(CONFIG)/$(EXEFILE)
endif

debug:
ifndef CONFIG
	-@mkdir -p $(@)
	@$(MAKE) $(MAKEFILE) --no-print-directory $(@) CONFIG=$(@)
else
	@$(MAKE) $(MAKEFILE) --no-print-directory $(CONFIG)/$(EXEFILE)
endif

all: release debug

cleanrelease:
	-rm -rf release/*.o release/*.mod
cleandebug:
	-rm -rf debug/*.o debug/*.mod
clean: cleanrelease cleandebug
cleaner:
	-rm -rf debug/*.o debug/*.mod debug/$(EXEFILE)
	-rm -rf release/*.o release/*.mod release/$(EXEFILE)
cleanest:
	-rm -rf debug release

# This phony target is just used to output some info and error check

ptarg:
	@echo CONFIG=$(CONFIG)
ifeq ($(COMPILER_SUPPORTED), yes)
	@echo Using COMPILER_TYPE=$(COMPILER_TYPE) and settings for MACHINE=$(MACHINE) PREC=$(PREC) OPTLPKBLS=$(OPTLPKBLS) ...
else
	@echo Error - there is no support for COMPILER_TYPE=$(COMPILER_TYPE) and MACHINE=$(MACHINE) PREC=$(PREC) OPTLPKBLS=$(OPTLPKBLS)
	@false
endif

# The real targets

$(CONFIG)/$(EXEFILE): ptarg $(OBJS) $(CUDA_OBJS)
	$(FC) $(OBJS) $(CUDA_OBJS) -o $(CONFIG)/$(EXEFILE) $(LINKFLAGS) $(MPILINKFLAGS)

$(CONFIG)/wp_def_$(PREC).o: $(SRC_DIR)/wp_def_$(PREC).f90
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/wp_def_$(PREC).f90 -o $(CONFIG)/wp_def_$(PREC).o

$(CONFIG)/BLAS.o: $(SRC_DIR)/BLAS.f $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/BLAS.f -o $(CONFIG)/BLAS.o

$(CONFIG)/LAPACK.o: $(SRC_DIR)/LAPACK.f $(CONFIG)/BLAS.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/LAPACK.f -o $(CONFIG)/LAPACK.o

$(CONFIG)/X1MACH.o: $(SRC_DIR)/X1MACH.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/X1MACH.f90 -o $(CONFIG)/X1MACH.o

$(CONFIG)/dmng.o: $(SRC_DIR)/dmng.f $(CONFIG)/X1MACH.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/dmng.f -o $(CONFIG)/dmng.o

$(CONFIG)/globvars.o: $(SRC_DIR)/globvars.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/globvars.f90 -o $(CONFIG)/globvars.o

$(CONFIG)/linalg.o: $(SRC_DIR)/linalg.f90 $(CONFIG)/globvars.o $(OBJ_LPKBLS)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/linalg.f90 -o $(CONFIG)/linalg.o

$(CONFIG)/spin.o: $(SRC_DIR)/spin.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/spin.f90 -o $(CONFIG)/spin.o

$(CONFIG)/misc.o: $(SRC_DIR)/misc.f90 $(CONFIG)/globvars.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/misc.f90 -o $(CONFIG)/misc.o

$(CONFIG)/matelem.o: $(SRC_DIR)/matelem.f90 $(CONFIG)/globvars.o $(CONFIG)/spin.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matelem.f90 -o $(CONFIG)/matelem.o

$(CONFIG)/matelem_gpu.o: $(SRC_DIR)/matelem_gpu.f90
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matelem_gpu.f90 -o $(CONFIG)/matelem_gpu.o

$(CONFIG)/matelem_cuda.o: $(SRC_DIR)/matelem_cuda.cu
	$(NVCC) -c $(SRC_DIR)/matelem_cuda.cu -o $(CONFIG)/matelem_cuda.o \
	  -O3 -arch=$(CUDA_ARCH) --compiler-options '-fPIC'

$(CONFIG)/matelem_acc.o: $(SRC_DIR)/matelem_acc.f90 $(CONFIG)/globvars.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matelem_acc.f90 -o $(CONFIG)/matelem_acc.o

$(CONFIG)/matform.o: $(SRC_DIR)/matform.f90 $(CONFIG)/matelem.o $(GPU_DEP) $(ACC_DEP)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matform.f90 -o $(CONFIG)/matform.o

$(CONFIG)/workproc.o: $(SRC_DIR)/workproc.f90 $(CONFIG)/misc.o $(CONFIG)/linalg.o $(CONFIG)/matform.o $(CONFIG)/matelem.o $(CONFIG)/dmng.o $(OBJS_LPKBLS)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/workproc.f90 -o $(CONFIG)/workproc.o

$(CONFIG)/main.o: $(SRC_DIR)/main.f90 $(CONFIG)/workproc.o $(GPU_DEP)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/main.f90 -o $(CONFIG)/main.o

