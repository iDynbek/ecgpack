# This is a makefile to compile an ECG code. Examples of compilations with this
# makefile:
#    make release COMPILER_TYPE=intel MACHINE=mylinuxbox PREC=8 OPTLPKBLS=no EXEFILE=ecg  
#    make debug   COMPILER_TYPE=intel MACHINE=mylinuxbox PREC=8 OPTLPKBLS=no EXEFILE=ecg  
#    make all     COMPILER_TYPE=intel MACHINE=mylinuxbox PREC=8 OPTLPKBLS=no EXEFILE=ecg  

# Set default compiler type (intel, gnu, sun, etc..)
COMPILER_TYPE=intel

# Set default machine where the compilation takes place (mylinuxbox, ice, marin, etc.)
MACHINE=mylinuxbox

# Set the default precision (kind of real type, 8, 10, or 16)
# Value PREC=8 corresponds to standard double precision
# Value PREC=16 corresponds to quadruple precision
PREC=8

# Set whether selected linear algebra routines should be taken from optimized
# LAPACK/BLAS library (OPTLPKBLS=yes) such as Intel Math Kernel library or Sun Performance
# Library. Otherwise (OPTLPKBLS=no or any other value) the supplied nonoptimized versions 
# of these routines will be used.
OPTLPKBLS=no

# Set the default name for the generated executable file
EXEFILE=ecg

# Directory where all source files are located
SRC_DIR=./src

# Compiler flags for each machine and compiler type (for both release and debug builds). 

# My work desktop computer 
ifeq ($(MACHINE), mylinuxbox) 
  # Intel Fortran and Open MPI (mylinuxbox)
  ifeq ($(COMPILER_TYPE), intel)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -traceback
      export LINKFLAGS+=-C -O0 -g -traceback
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ipo -xhost -fp-model precise
      export LINKFLAGS+=-O3 -ipo -xhost -fp-model precise
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
  # gfortran and Open MPI (mylinuxbox)
  ifeq ($(COMPILER_TYPE), gnu)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-J$(CONFIG)
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3
      export LINKFLAGS+=-O3
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export LINKFLAGS+=-lblas -llapack
      endif
    endif
    ifeq ($(PREC), 10)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # Sun Fortran and Sun MPI (mylinuxbox)
  ifeq ($(COMPILER_TYPE), sun)
    FC=/opt/sun/sunstudio12.1/bin/f90
    export MPICOMPFLAGS=$(shell /opt/SUNWhpc/HPC8.2.1/sun/bin/mpif90 --showme:compile)
    export MPILINKFLAGS=$(shell /opt/SUNWhpc/HPC8.2.1/sun/bin/mpif90 --showme:link)
    export FCFLAGS=-M$(CONFIG) -moddir=$(CONFIG)
    export LINKFLAGS=-M$(CONFIG) -moddir=$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-fast
      export LINKFLAGS+=-fast
    endif
    ifeq ($(PREC), 8)
      export COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with Sun Performance Library
        export LINKFLAGS+=-xlic_lib=sunperf
      endif
    endif
  endif
endif

ifeq ($(MACHINE), bluehive) 
  # Intel Fortran and Intel MPI (bluehive)
  ifeq ($(COMPILER_TYPE), intel)
    FC=mpiifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ipo -xhost -fp-model precise
      export LINKFLAGS+=-O3 -ipo -xhost -fp-model precise
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        #export LINKFLAGS+=-L$(MKLROOT) 
        export LINKFLAGS+=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
        export LINKFLAGS+=-L/usr/local/intel/mkl/10.2.0.013/lib/em64t
      endif
    endif
    ifeq ($(PREC), 16)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # gfortran and Open MPI (bluehive)
  ifeq ($(COMPILER_TYPE), gnu)
    FC=gfortran
    export MPICOMPFLAGS=$(shell mpif90 --showme:compile)
    export MPILINKFLAGS=$(shell mpif90 --showme:link)
    export FCFLAGS=-J$(CONFIG)
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -march=native -mtune=native
      export LINKFLAGS+=-O3 -march=native -mtune=native
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export LINKFLAGS+=-lblas -llapack
      endif
    endif
    ifeq ($(PREC), 10)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
endif

ifeq ($(MACHINE), ocelote) 
  # Intel Fortran and Intel MPI (ocelote)
  ifeq ($(COMPILER_TYPE), intel)
    FC=mpiifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -traceback
      export LINKFLAGS+=-C -O0 -g -traceback
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ipo -xhost -fp-model precise
      export LINKFLAGS+=-O3 -ipo -xhost -fp-model precise
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
  # gfortran and Open MPI (ocelote, same as mylinuxbox)
  ifeq ($(COMPILER_TYPE), gnu)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-J$(CONFIG)
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3
      export LINKFLAGS+=-O3
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export LINKFLAGS+=-lblas -llapack
      endif
    endif
    ifeq ($(PREC), 10)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
endif

ifeq ($(MACHINE), elgato) 
  # Intel Fortran and Open MPI (elgato)
  ifeq ($(COMPILER_TYPE), intel)
    FC=mpiifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -traceback
      export LINKFLAGS+=-C -O0 -g -traceback
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ipo -xhost -fp-model precise
      export LINKFLAGS+=-O3 -ipo -xhost -fp-model precise
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
  # gfortran and Open MPI (elgato, same as mylinuxbox)
  ifeq ($(COMPILER_TYPE), gnu)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-J$(CONFIG)
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O2
      export LINKFLAGS+=-O2
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export LINKFLAGS+=-lblas -llapack
      endif
    endif
    ifeq ($(PREC), 10)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
endif

ifeq ($(MACHINE), ice_8400) 
  # Intel Fortran and Open MPI (ice_8400)
  ifeq ($(COMPILER_TYPE), intel)
    FC=mpiifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g -traceback
      export LINKFLAGS+=-C -O0 -g -traceback
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ipo -xhost -fp-model precise
      export LINKFLAGS+=-O3 -ipo -xhost -fp-model precise
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        #export LINKFLAGS+=-L$(MKLROOT) 
        export LINKFLAGS+=-lmkl_solver_lp64_sequential -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
        export LINKFLAGS+=-L/uaopt/intel/2012.0.032/mkl/lib/intel64 
      endif
    endif
    ifeq ($(PREC), 16)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # gfortran and Open MPI (ice_8400)
  ifeq ($(COMPILER_TYPE), gnu)
    FC=gfortran
    export MPICOMPFLAGS=$(shell mpif90 --showme:compile)
    export MPILINKFLAGS=$(shell mpif90 --showme:link)
    export FCFLAGS=-J$(CONFIG)
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -march=native -mtune=native
      export LINKFLAGS+=-O3 -march=native -mtune=native
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export LINKFLAGS+=-L/uaopt/lapack/3.4.0/lib -L/uaopt/blas/lib -lblas -llapack
      endif
    endif
    ifeq ($(PREC), 10)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
endif


# ICE cluster at U of A
ifeq ($(MACHINE), ice) 
  # Intel Fortran and Intel MPI (ice)
  # User must load modules "intel" and "intel-mpi"
  # prior to compilation
  ifeq ($(COMPILER_TYPE), intel)
    FC=mpiifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ipo
      export LINKFLAGS+=-O3 -ipo
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export LINKFLAGS+=-L/uaopt/intel/fc/11.1/064/mkl/lib/em64t/        
        export LINKFLAGS+=-Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
      endif
    endif
    ifeq ($(PREC), 16)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif
  # gfortran and Intel MPI (ice)
  # User must load modules "gcc" and "intel-mpi"
  # prior to compilation
  ifeq ($(COMPILER_TYPE), gnu)
    FC=mpif90
    export MPICOMPFLAGS=
    export MPILINKFLAGS=
    export FCFLAGS=-J$(CONFIG)
    export LINKFLAGS=-J$(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3
      export LINKFLAGS+=-O3
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        # Link with MKL
        export LINKFLAGS+=-L/uaopt/intel/fc/11.1/064/mkl/lib/em64t/
        export LINKFLAGS+=-Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
      endif
    endif
    ifeq ($(PREC), 10)
      ifneq ($(OPTLPKBLS),yes)
        COMPILER_SUPPORTED=yes
      endif
    endif
  endif  
endif

# MARIN supercomputer at U of A
ifeq ($(MACHINE), marin) 
  # Intel Fortran and Intel MPI (marin)
  # User must load module "intel-fc-10.0.025"
  # prior to compilation
  ifeq ($(COMPILER_TYPE), intel)
    FC=ifort
    export MPICOMPFLAGS=
    export MPILINKFLAGS=-lmpi
    export FCFLAGS=-module $(CONFIG)
    export LINKFLAGS=-module $(CONFIG)
    ifeq ($(CONFIG), debug)
      export FCFLAGS+=-C -O0 -g
      export LINKFLAGS+=-C -O0 -g
    endif
    ifeq ($(CONFIG), release)
      export FCFLAGS+=-O3 -ipo -ftz -fno-alias -fno-fnalias -align all
      export LINKFLAGS+=-O3 -ipo -ftz -fno-alias -fno-fnalias -align all
    endif
    ifeq ($(PREC), 8)
      COMPILER_SUPPORTED=yes
      ifeq ($(OPTLPKBLS),yes)
        export MKLROOT=/opt/intel/mkl/10.0.2.018
        export LINKFLAGS+=-L$(MKLROOT)/lib/64 -lmkl_solver_lp64_sequential -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
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
              matelem.o matform.o workproc.o main.o
export OBJS=$(addprefix $(CONFIG)/,$(BARE_OBJS))
export OBJS+=$(OBJS_LPKBLS)


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

$(CONFIG)/$(EXEFILE): ptarg $(OBJS)
	$(FC) $(OBJS) -o $(CONFIG)/$(EXEFILE) $(LINKFLAGS) $(MPILINKFLAGS)

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

$(CONFIG)/misc.o: $(SRC_DIR)/misc.f90 $(CONFIG)/globvars.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/misc.f90 -o $(CONFIG)/misc.o

$(CONFIG)/matelem.o: $(SRC_DIR)/matelem.f90 $(CONFIG)/globvars.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matelem.f90 -o $(CONFIG)/matelem.o

$(CONFIG)/matform.o: $(SRC_DIR)/matform.f90 $(CONFIG)/matelem.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matform.f90 -o $(CONFIG)/matform.o

$(CONFIG)/workproc.o: $(SRC_DIR)/workproc.f90 $(CONFIG)/misc.o $(CONFIG)/linalg.o $(CONFIG)/matform.o $(CONFIG)/matelem.o $(CONFIG)/dmng.o $(OBJS_LPKBLS)
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/workproc.f90 -o $(CONFIG)/workproc.o

$(CONFIG)/main.o: $(SRC_DIR)/main.f90 $(CONFIG)/workproc.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/main.f90 -o $(CONFIG)/main.o

