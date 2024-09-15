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
    export FCFLAGS=-J$(CONFIG)
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
    export FCFLAGS=-J$(CONFIG)
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

export BARE_OBJS=wp_def_$(PREC).o globvars.o matelem.o workproc.o main.o spin.o
export OBJS=$(addprefix $(CONFIG)/,$(BARE_OBJS))


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
	@echo Using COMPILER_TYPE=$(COMPILER_TYPE) and settings for MACHINE=$(MACHINE) PREC=$(PREC) ...
else
	@echo Error - there is no support for COMPILER_TYPE=$(COMPILER_TYPE) and MACHINE=$(MACHINE) PREC=$(PREC)
	@false
endif

# The real targets

$(CONFIG)/$(EXEFILE): ptarg $(OBJS)
	$(FC) $(OBJS) -o $(CONFIG)/$(EXEFILE) $(LINKFLAGS) $(MPILINKFLAGS)

$(CONFIG)/wp_def_$(PREC).o: $(SRC_DIR)/wp_def_$(PREC).f90
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/wp_def_$(PREC).f90 -o $(CONFIG)/wp_def_$(PREC).o

$(CONFIG)/globvars.o: $(SRC_DIR)/globvars.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/globvars.f90 -o $(CONFIG)/globvars.o

$(CONFIG)/spin.o: $(SRC_DIR)/spin.f90 $(CONFIG)/wp_def_$(PREC).o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/spin.f90 -o $(CONFIG)/spin.o

$(CONFIG)/matelem.o: $(SRC_DIR)/matelem.f90 $(CONFIG)/globvars.o $(CONFIG)/spin.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/matelem.f90 -o $(CONFIG)/matelem.o

$(CONFIG)/workproc.o: $(SRC_DIR)/workproc.f90 $(CONFIG)/matelem.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/workproc.f90 -o $(CONFIG)/workproc.o

$(CONFIG)/main.o: $(SRC_DIR)/main.f90 $(CONFIG)/workproc.o
	$(FC) $(FCFLAGS) $(MPICOMPFLAGS) -c $(SRC_DIR)/main.f90 -o $(CONFIG)/main.o

