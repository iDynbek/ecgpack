#!/bin/csh
# Compilation script to compile ECG codes. Typical usages:
#   ./compile_XXXXX shabyt gcc9 8 10   <--  gcc9 toolchain, precision 8 and 10
#   ./compile_XXXXX wsl_ubuntu20 gcc 8 10 16   <-- system default gcc, precisions 8, 10, 16
#   ./compile_XXXXX worker26 gcc12 8 10 16   <-- system default gcc, precisions 8, 10, 16
# This script must be executed in the directory one level above
# directories RGL0, RGL1, RGL2P, RGL2D

# Set the default machine (if the script executed with no argument)
# Possible choices are: shabyt, mylinuxbox, wsl_ubuntu20
if ($#argv < 1) then
  set MACHINE = wsl_ubuntu20
else
  set MACHINE = $argv[1]
endif

#Set the default compiler toolchain (if the script is executed with no more than 1 argument)
#Possible choices: gcc (system default), gcc9
if ($#argv < 2) then
  set TOOLCHAIN=gcc9
else
  set TOOLCHAIN = $argv[2]
endif

#Read which precisions should be used (beginning with the third argument)
set PRECVALS=()

if ($#argv > 2) then
  set i=3
  while ($i <= $#argv)
    if ($argv[$i] == 8 || $argv[$i] == 10 || $argv[$i] == 16) then
      set PRECVALS=($PRECVALS $argv[$i])
    else
      echo "Wrong precision specified: "  $argv[$i]
      exit 1
    endif
    @ i++
  end
else
  set PRECVALS=($PRECVALS 8)
  set PRECVALS=($PRECVALS 10)
  set PRECVALS=($PRECVALS 16)
endif

# Set configuration (debug or release). This is not read as a command line argument,
#but it can be set manually here. Generally, one should always use release, unless 
#there is a need to debug codes
set CONFIG = release

set BINSUBDIR = ${TOOLCHAIN}

#Initialize Lmod (or environment modules)
if (${MACHINE} == mylinuxbox) then
  source /usr/share/modules/init/csh
else if (${MACHINE} == shabyt) then
  source /usr/share/lmod/lmod/init/csh
  module use /shared/modulefiles
else if (${MACHINE} == wsl_ubuntu20) then
  source /usr/share/lmod/lmod/init/csh
  module use /shared/modulefiles
else if (${MACHINE} == worker26) then
  source /usr/share/modules/init/csh
endif

#Load moduli corresponding to each specific compiler and machine
if (${MACHINE} == shabyt || ${MACHINE} == wsl_ubuntu20) then
  if (${TOOLCHAIN} == gcc9) then
    module load gcc/9.5.0
    module load openmpi/gcc9/4.1.5
  else if (${TOOLCHAIN} == gcc12) then
    module load gcc/12.2.0
    module load openmpi/gcc12/4.1.5
  else if (${TOOLCHAIN} == gcc) then
    #Nothing to load for wsl_ubuntu20
    #System default gcc/openmpi should be accessible without loading any moduli
    if (${MACHINE} == shabyt) module load mpi/openmpi-x86_64
  else
    echo "Error: Toolchain " ${TOOLCHAIN} " is not supported on machine " ${MACHINE}
    exit 1
  endif
else if (${MACHINE} == worker26 || ${MACHINE} == mylinuxbox) then
  if (${TOOLCHAIN} == gcc9) then
    module load gcc/9.5.0
    module load openmpi-gcc9/4.1.5
  else if (${TOOLCHAIN} == gcc12) then
    module load gcc/12.2.0
    module load openmpi-gcc12/4.1.5
  else if (${TOOLCHAIN} == gcc) then
    #Nothing to load for worker26 or mylinuxbox
    #System default gcc/openmpi should be accessible without loading any moduli
  else
    echo "Error: Toolchain " ${TOOLCHAIN} " is not supported on machine " ${MACHINE}
    exit 1
  endif
else
    echo "Machine " ${MACHINE} " is not supported"
    exit 1
endif


# Create directory called "bin" (and subdirectories) in current directory if it does not yet exist
if (! -e bin/) mkdir bin
if (! -e bin/${BINSUBDIR}) mkdir bin/${BINSUBDIR}
if (! -e bin/${BINSUBDIR}/${CONFIG}) mkdir bin/${BINSUBDIR}/${CONFIG}

foreach PREC ($PRECVALS)

  #Set parameters COMPILER, REALX_NAME depending on the machine/compiler used
  if (${TOOLCHAIN} == gcc || ${TOOLCHAIN} == gcc9 || ${TOOLCHAIN} == gcc12) then
    set COMPILER=gnu
    if (${MACHINE} == shabyt || ${MACHINE} == wsl_ubuntu20 || ${MACHINE} == worker26 || ${MACHINE} == mylinuxbox) then
      if (${PREC} == 8) then
        set REALX_NAME=MPI_DOUBLE_PRECISION
      else if (${PREC} == 10) then
        set REALX_NAME=MPI_REAL16
      else if (${PREC} == 16) then
        set REALX_NAME=MPI_REAL16
      endif
    else
      echo "Toolchain/Machine commbination not supported: " ${TOOLCHAIN}/${MACHINE}
      exit 1
    endif
  endif

  #This specifies if optimized BLAS/LAPACK should be invoked
  if (${PREC} == 8) then
    set OPTLPKBLS="yes"
  else
    set OPTLPKBLS="no"
  endif

  foreach L (0 1 2P 2D)
    cd RGL${L}
    foreach NPARTICLES (3 4 5 6 7 8)
      echo " "
      echo " ==========================================================="
      echo " Compiling code with the following options: "
      echo " L=${L} PREC=${PREC} NPARTICLES=${NPARTICLES} CONFIG=${CONFIG}"
      echo " COMPILER=${COMPILER} TOOLCHAIN=${TOOLCHAIN} MACHINE=${MACHINE}"
      echo " -----------------------------------------------------------"
      echo " mpif90 wrapper check:"
      mpif90 --version | head -1
      echo " -----------------------------------------------------------"

      # Replace "Glob_MaxAllowedNumOfParticles=..." with "Glob_MaxAllowedNumOfParticles=${NPARTICLES}" in file src/wp_def_${PREC}.f90 
      sed -e "s/Glob_MaxAllowedNumOfParticles=[0-9][0-9]*/Glob_MaxAllowedNumOfParticles=${NPARTICLES}/" < src/wp_def_${PREC}.f90 > src/temp.f90
      mv -f src/temp.f90 src/wp_def_${PREC}.f90
      # Replace in file "MPI_DPREC=..." with "MPI_DPREC=${REALX_NAME}" src/wp_def_${PREC}.f90
      sed -e "s/MPI_DPREC=[^ ][^ ]*/MPI_DPREC=${REALX_NAME}/" < src/wp_def_${PREC}.f90 > src/temp.f90
      mv -f src/temp.f90 src/wp_def_${PREC}.f90

      make clean
      make ${CONFIG} COMPILER_TYPE=${COMPILER} MACHINE=${MACHINE} PREC=${PREC} OPTLPKBLS=${OPTLPKBLS} EXEFILE=ecg
      mv ${CONFIG}/ecg ../bin/${BINSUBDIR}/${CONFIG}/RGL${L}_N${NPARTICLES}_${PREC}

      #Replace back the above two values to default
      sed -e "s/Glob_MaxAllowedNumOfParticles=[0-9][0-9]*/Glob_MaxAllowedNumOfParticles=5/" < src/wp_def_${PREC}.f90 > src/temp.f90
      mv -f src/temp.f90 src/wp_def_${PREC}.f90
      sed -e "s/MPI_DPREC=[^ ][^ ]*/MPI_DPREC=MPI_DOUBLE_PRECISION/" < src/wp_def_${PREC}.f90 > src/temp.f90
      mv -f src/temp.f90 src/wp_def_${PREC}.f90
      echo " ==========================================================="
    end
    cd ../
  end

end



exit 0

