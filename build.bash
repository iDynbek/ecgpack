#!/bin/bash
# This is a script that builds all ECGPACK codes in different environments and on different machines
# For usage/help run it without any arguments

# Function to display usage information
usage_print() {
  echo "Missing arguments or invalid arguments."
  echo ""  
  echo "PROPER USAGE:"
  echo "$0 machine=<machinename> toolchain=<toolchainnames> config=<confignames> code=<codenames> nparticles=<nparticles> precision=<precisions> use_optimized_lapack=<yesnoflag>"
  echo ""
  echo "NOTE:"
  echo "All arguments are optional except nparticles. If multiple values are specified for an argument, they must be separated by a comma."
  echo "The resulting binary files are stored in the directory bin/<toolchainname>/<configname>/. This script must be executed in the directory where it is located (ecgpack root directory)."
  echo ""
  echo "DESCRIPTION OF ARGUMENTS:"
  echo "<machinename> is the name of the machine. Only a single value may be specified. It could be linux-generic (default), ubuntu-generic, irgetas, shabyt, etc."
  echo "<toolchainnames> are the names of the toolchains." 
  echo "<confignames> is the list of configurations that need to be built. Currently these could be debug or release. The default includes all configurations."
  echo "<codenames> is the list of codes that need to be built. Currently these could be RG_0S,RG_1P,RG_2P,RG_2D,RG_0S-1P,RG_0S2P,RG_0S-2D,RG_1P-1P,RG_2P-2P,RG_2P-2D. The default includes all codes."
  echo "<nparticles> defines for how many particles each code must be build for. There is no default value. This argument must be present."
  echo "<precisions> is the kind parameter for real type. 8 corresponds to double precision (fp64), 10 corresponds to extended precision (fp80), 16 corresponds to quadruple precision. Different compilers/toolchain support different kinds. For example, Intel compilers supports only 8 and 16, while modern GNU compilers support 8, 10, and 16. The default value is 8."
  echo "<yesnoflag> is a flag that specifies whether to use optimized LAPACK libraries (yes) or use nonoptimized LAPACK from source (no). The default value is yes for precision=8. For precision=10 and precision=16 the value is always no, regardless how the flag is set because optimized LAPACK for these two cases is unavailable."
  echo "" 
  echo "Supported toolchains on different machines are listed below."
  echo ""   
  echo "Machine: linux-generic (default), ubuntu-generic, irgetas, shabyt, muon"
  echo "Supported toolchain names :"
  echo "    systemdefault (system default compiler and MPI installed on the system, if any)"
  echo "    foss-2025b (requires Easybuild module foss/2025b)"
  echo "    foss-2025a (requires Easybuild module foss/2025a)"
  echo "    foss-2024a (requires Easybuild module foss/2024a)"
  echo "    foss-2023b (requires Easybuild module foss/2023b)"
  echo "    intel-2025b (requires Easybuild module intel/2025b)"
  echo "    intel-2025a (requires Easybuild module intel/2025a)"
  echo "    intel-2024a (requires Easybuild module intel/2024a)"  
  echo "    intel-2023b (requires Easybuild module intel/2023b)"
  echo "    nvhpc-25.9 (requires Easybuild module NVHPC/25.9-CUDA-12.9.1)"
  echo "    nvhpc-25.3 (requires Easybuild module NVHPC/25.3-CUDA-12.8.0)"  
  echo ""
   echo "Machine: puma, ocelote, elgato"
  echo "Supported toolchain names :"
  echo "    systemdefault (system default compiler and MPI installed on the system)"  
  echo "    gnu8-8.3.0 (requires module gnu8/8.3.0)"
  echo "    intel-2020.4 (requires module intel/2020.4)"

  echo "EXECUTION EXAMPLES:"
  echo ""  
  echo "    $0 machine=linux-generic toolchain=foss-2025b config=release code=RG_0S,RG_1P,RG_2D,RG_2P nparticles=3,4,5,6,7,8 precision=8,10,16 use_optimized_lapack=no"
  echo ""  
  echo "    $0 machine=linux-generic toolchain=foss-2025a config=release code=RG_0S,RG_1P,RG_2D,RG_2P,RG_0S-1P,RG_1P-2D,RG_1P-2P,RG_0S-2D,RG_0S-2P,RG_1P-1P,RG_2D-2D,RG_2P-2D,RG_2P-2P nparticles=3,4,5,6,7,8 precision=8 use_optimized_lapack=yes"
  echo ""  
  echo "    $0 machine=irgetas toolchain=foss-2023b config=release code=RG_0S,RG_1P nparticles=3,4,5,6 precision=8,10,16 use_optimized_lapack=no"
  echo ""  
  echo "    $0 machine=muon toolchain=foss-2023b,intel-2023b config=debug,release code=RG_0S nparticles=4,5 precision=8,16 use_optimized_lapack=yes,no"
  echo "" 
  echo "    $0 machine=puma toolchain=intel-2020.4 config=debug,release code=RG_0S nparticles=5 precision=8,16 use_optimized_lapack=yes,no"  
  echo ""  
  echo "    $0 machine=shabyt nparticles=4,5 precision=10"
  echo ""  
  echo "    $0 nparticles=6"
  exit 1
}

# Set default values of arguments
machine="linux-generic"
toolchain="systemdefault"
config="debug,release"
# code="CG_0S, RG_0S, RG_1P, RG_2D, RG_2P, RG_0S-1P, RG_1P-2D, RG_1P-2P, RG_0S-2D, RG_0S-2P, RG_1P-1P, RG_2D-2D, RG_2P-2D, RG_2P-2P"
code="RG_0S, RG_1P, RG_2D, RG_2P, RG_0S-1P, RG_1P-2D, RG_1P-2P, RG_0S-2D, RG_0S-2P, RG_1P-1P, RG_2D-2D, RG_2P-2D, RG_2P-2P"
nparticles=""
precision="8"
use_optimized_lapack="yes"

# Parse the arguments
for arg in "$@"; do
  case $arg in
    machine=*) machine="${arg#*=}" ;;
    toolchain=*) toolchain="${arg#*=}" ;;
    config=*) config="${arg#*=}" ;;
    code=*) code="${arg#*=}" ;;
    nparticles=*) nparticles="${arg#*=}" ;;
    precision=*) precision="${arg#*=}" ;;
    use_optimized_lapack=*) use_optimized_lapack="${arg#*=}" ;;
    *) echo "ERROR, INVALID ARGUMENT: $arg" ; echo "" ; usage_print ;;
  esac
done

# Check if required argument 'nparticles' is set
if [[ -z "$nparticles" ]]; then
  echo "ERROR, MISSING REQUIRED ARGUMENT: nparticles"
  usage_print
  exit 1
fi

# Parse all arguments with possible multiple values
IFS=',' read -ra machine_list <<< $machine
IFS=',' read -ra toolchain_list <<< $toolchain
IFS=',' read -ra config_list <<< $config
IFS=',' read -ra code_list <<< $code
IFS=',' read -ra nparticles_list <<< $nparticles
IFS=',' read -ra precision_list <<< $precision
IFS=',' read -ra use_optimized_lapack_list <<< $use_optimized_lapack

# Check if machine is set properly. Only a single machine can be used as an argument
if [[ " linux-generic ubuntu-generic shabyt muon puma ocelote elgato " != *" $machine "* ]]; then
  echo "ERROR, INVALID MACHINE: $machine"
  usage_print
  exit 1
fi

# Chech if config is set properly
for config_value in ${config_list[@]}; do
  if [[ " debug release " != *" $config_value "* ]]; then
    echo "ERROR, WRONG VALUE(S) OF ARGUMENT: config"
    usage_print
    exit 1
  fi
done

# Check if code is set properly
for code_value in ${code_list[@]}; do
  if [[ " CG_0S RG_0S RG_1P RG_2D RG_2P RG_0S-1P RG_0S-2D RG_0S-2P RG_1P-1P RG_1P-2D RG_2D-2D RG_2P-2D RG_2P-2P " != *" $code_value "* ]]; then
    echo "ERROR, WRONG VALUE(S) OF ARGUMENT: code"
    usage_print
    exit 1
  fi
done

# Check if nparticles is set properly
for nparticles_value in ${nparticles_list[@]}; do
  if ! [[ $nparticles_value =~ ^[0-9]+$ ]] ; then
    echo "ERROR, NONINTEGER VALUE(S) OF ARGUMENT: nparticles"
    usage_print
    exit 1
  fi
  if [ $nparticles_value -le 1 ] || [ $nparticles_value -ge 20 ]; then
    echo "ERROR, OUT-OF-RANGE VALUE(S) OF ARGUMENT: nparticles"
    usage_print
    exit 1
  fi
done

# Check if precision is set properly
for precision_value in ${precision_list[@]}; do
  if [[ " 8 10 16 " != *" $precision_value "* ]]; then
    echo "ERROR, WRONG VALUE(S) OF ARGUMENT: precision"
    usage_print
    exit 1
  fi
done

# Check if use_optimized_lapack is set properly
for use_optimized_lapack_value in ${use_optimized_lapack_list[@]}; do
  if [[ " yes no " != *" $use_optimized_lapack_value "* ]]; then
    echo "ERROR, WRONG VALUE(S) OF ARGUMENT: use_optimized_lapack"
    usage_print
    exit 1
  fi
done

# Set the name of the directory where all binaries will be stored
bindirname="bin"

# Set the counters for counting builds
counter_attempted_builds=0
counter_successful_builds=0
counter_failed_builds=0

# Helper that prepares a toolchain: it loads the corresponding Lmod module(s)
# (one or more, space-separated; skipped when <modules> is empty, e.g.
# systemdefault), sets the global compiler, and verifies that the Fortran
# compiler and MPI wrapper resolve to their expected installation paths. On any
# failure it prints a diagnostic and returns 1, so a caller can write
# `load_toolchain ... || continue`.
# Usage: load_toolchain "<modules>" <compiler> <fc_bin> <fc_path> <mpi_bin> <mpi_path>
load_toolchain() {
  module_for_toolchain="$1"
  compiler="$2"
  local fc_bin="$3" fc_path="$4" mpi_bin="$5" mpi_path="$6" hint
  if [ -n "$module_for_toolchain" ]; then
    hint="Check that module $module_for_toolchain is properly configured."
    module load $module_for_toolchain
    if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
      echo "Module $module_for_toolchain for toolchain $toolchain_value could not be loaded"
      echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
      echo "Skipping this toolchain."
      return 1
    fi
  else
    hint="Check that it is installed on your system."
  fi
  if ! which "$fc_bin" 2>&1 | grep -qi "$fc_path"; then
    echo "$fc_bin compiler matching toolchain $toolchain_value is not accessible."
    echo "$hint"
    echo "Skipping this toolchain."
    return 1
  fi
  if ! which "$mpi_bin" 2>&1 | grep -qi "$mpi_path"; then
    echo "$mpi_bin compiler matching toolchain $toolchain_value is not accessible."
    echo "$hint"
    echo "Skipping this toolchain."
    return 1
  fi
  return 0
}


for toolchain_value in ${toolchain_list[@]}; do
  # Drop any modules that may have been left loaded and suppress any output of this command
  module reset > /dev/null 2>&1
  # Check if the toolchain is valid for the specified machine and load the corresponding module.
  # If the module is not loaded or/and the proper Fortran compiler and MPI are inaccessible,
  # print an error message.
  if [[ " linux-generic ubuntu-generic irgetas shabyt muon " == *" $machine "* ]]; then 
    # Possibly insert $machine in the path given by $bindirname (this can be left empty)
    machinedirname=""
    if [ "$toolchain_value" = "foss-2025b" ]; then
      load_toolchain "foss/2025b" gfortran gfortran "software/GCCcore/14.3.0/bin/gfortran" mpif90 "software/OpenMPI/5.0.8-GCC-14.3.0/bin/mpif90" || continue
    elif [ "$toolchain_value" = "foss-2025a" ]; then
      load_toolchain "foss/2025a" gfortran gfortran "software/GCCcore/14.2.0/bin/gfortran" mpif90 "software/OpenMPI/5.0.7-GCC-14.2.0/bin/mpif90" || continue
    elif [ "$toolchain_value" = "foss-2024a" ]; then
      load_toolchain "foss/2024a" gfortran gfortran "software/GCCcore/13.3.0/bin/gfortran" mpif90 "software/OpenMPI/5.0.3-GCC-13.3.0/bin/mpif90" || continue
    elif [ "$toolchain_value" = "foss-2023b" ]; then
      load_toolchain "foss/2023b" gfortran gfortran "software/GCCcore/13.2.0/bin/gfortran" mpif90 "software/OpenMPI/4.1.6-GCC-13.2.0/bin/mpif90" || continue
    elif [ "$toolchain_value" = "intel-2025b" ]; then
      load_toolchain "intel/2025b" ifx ifx "software/intel-compilers/2025.2.0/compiler/2025.2/bin/ifx" mpiifx "software/impi/2021.16.1-intel-compilers-2025.2.0/mpi/2021.16/bin/mpiifx" || continue
    elif [ "$toolchain_value" = "intel-2025a" ]; then
      load_toolchain "intel/2025a" ifx ifx "software/intel-compilers/2025.1.1/compiler/2025.1/bin/ifx" mpiifx "software/impi/2021.15.0-intel-compilers-2025.1.1/mpi/2021.15/bin/mpiifx" || continue
    elif [ "$toolchain_value" = "intel-2024a" ]; then
      load_toolchain "intel/2024a" ifx ifx "software/intel-compilers/2024.2.0/compiler/2024.2/bin/ifx" mpiifx "software/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/bin/mpiifx" || continue
    elif [ "$toolchain_value" = "intel-2023b" ]; then
      load_toolchain "intel/2023b" ifort ifort "software/intel-compilers/2023.2.1/compiler/2023.2.1/linux/bin/intel64/ifort" mpiifort "software/impi/2021.10.0-intel-compilers-2023.2.1/mpi/2021.10.0/bin/mpiifort" || continue
    elif [ "$toolchain_value" = "nvhpc-25.3" ]; then
      load_toolchain "NVHPC/25.3-CUDA-12.8.0" nvfortran nvfortran "software/nvidia-compilers/25.3-CUDA-12.8.0/Linux_x86_64/25.3/compilers/bin/nvfortran" mpif90 "software/NVHPC/25.3-CUDA-12.8.0/Linux_x86_64/25.3/comm_libs/hpcx/bin/mpif90" || continue
    elif [ "$toolchain_value" = "nvhpc-25.9" ]; then
      load_toolchain "NVHPC/25.9-CUDA-12.9.1" nvfortran nvfortran "software/nvidia-compilers/25.9-CUDA-12.9.1/Linux_x86_64/25.9/compilers/bin/nvfortran" mpif90 "software/NVHPC/25.9-CUDA-12.9.1/Linux_x86_64/25.9/comm_libs/hpcx/bin/mpif90" || continue
    elif [ "$toolchain_value" = "systemdefault" ]; then
      load_toolchain "" gfortran gfortran "/usr/bin/gfortran" mpif90 "/usr/bin/mpif90" || continue
    else
      echo "ERROR, INVALID TOOLCHAIN $toolchain FOR MACHINE $machine"
      usage_print
      exit 1
    fi
  fi
  if [[ " puma ocelote elgato " == *" $machine "* ]]; then
    # Possibly insert $machine in the path given by $bindirname (this can be left empty). Be sure to add a trailing slash.
    machinedirname=${machine}"/"
    if [ "$toolchain_value" = "gnu8-8.3.0" ]; then
      load_toolchain "gnu8/8.3.0" gfortran gfortran "/opt/ohpc/pub/compiler/gcc/8.3.0/bin/gfortran" mpif90 "/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.4/bin/mpif90" || continue
    elif [ "$toolchain_value" = "intel-2020.4" ]; then
      load_toolchain "intel/2020.4" ifort ifort "/opt/ohpc/pub/compiler/intel_2020_u4/compilers_and_libraries/linux/bin/intel64/ifort" mpiifort "/opt/ohpc/pub/compiler/intel_2020_u4/compilers_and_libraries/linux/mpi/intel64/bin/mpiifort" || continue
    elif [ "$toolchain_value" = "systemdefault" ]; then
      load_toolchain "" gfortran gfortran "/opt/ohpc/pub/compiler/gcc/8.3.0/bin/gfortran" mpif90 "/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.4/bin/mpif90" || continue
    else
      echo "ERROR, INVALID TOOLCHAIN $toolchain FOR MACHINE $machine"
      usage_print
      exit 1
    fi
  fi
  # Loop over all configurations
  for config_value in ${config_list[@]}; do
    # Loop over all codes
    for code_value in ${code_list[@]}; do
      # Loop over all nparticles
      for nparticles_value in ${nparticles_list[@]}; do
        # Loop over all precisions
        for precision_value in ${precision_list[@]}; do
          # Set the name of the real type for MPI that will be later inserted in the proper place of the Fortran source
          if [[ "$precision_value" = "8" ]]; then
            MPI_REALX_name=MPI_DOUBLE_PRECISION
          elif [[ "$precision_value" = "10" ]]; then
            MPI_REALX_name=MPI_REAL16
          elif [[ "$precision_value" = "16" ]]; then
            MPI_REALX_name=MPI_REAL16
          fi
          # Loop over all values of use_optimized_lapack, but for precision=10,16 only the `no` is allowed
          for use_optimized_lapack_value in ${use_optimized_lapack_list[@]}; do
            binsubdirname=${bindirname}/${machinedirname}${toolchain_value}/${config_value}
            binaryfilename=${code_value}_N${nparticles_value}_P${precision_value}
            if [[ "$use_optimized_lapack_value" = "yes" ]]; then
              if [[ "$precision_value" = "10" ]]; then
                continue
              elif [[ "$precision_value" = "16" ]]; then
                continue
              else
                # Add suffix to the binary file name if optimized LAPACK is used
                binaryfilename=${binaryfilename}_optlapack
              fi
            fi
            echo ""
            echo "════════════════════════ Starting a new build ═════════════════════════"
            echo "machine="$machine "   toolchain="$toolchain_value "   config="$config_value
            echo "code="$code_value "   nparticles="$nparticles_value "   precision="$precision_value "   use_optimized_lapack="$use_optimized_lapack_value
            echo "───────────────────────────── make output ─────────────────────────────"
            # Check if file ${code_value}/src/wp_def_${precision_value}.f90 exists. This way
            # we also automtically test if the directory ${code_value} for this specific code exists
            if ! [ -f "${code_value}/src/wp_def_${precision_value}.f90" ]; then
              echo "ERROR, FILE ${code_value}/src/wp_def_${precision_value}.f90 DOES NOT EXIST"
              echo "Skipping this build."
              continue
            fi
            cd ${code_value}
            counter_attempted_builds=$((counter_attempted_builds+1))
            # Make a copy of the original file src/wp_def_${precision_value}.f90
            cp -f src/wp_def_${precision_value}.f90 src/wp_def_temporary.f90
            # Replace "Glob_MaxAllowedNumOfParticles=..." with "Glob_MaxAllowedNumOfParticles=${nparticles_value}" in file src/wp_def_${precision_value}.f90
            sed -i "s/Glob_MaxAllowedNumOfParticles=[0-9]\+/Glob_MaxAllowedNumOfParticles=${nparticles_value}/g" src/wp_def_${precision_value}.f90
            # Replace "MPI_DPREC=..." with "MPI_DPREC=${MPI_REALX_name}" in file src/wp_def_${precision_value}.f90
            sed -i "s/MPI_DPREC=[^ ][^ ]*/MPI_DPREC=${MPI_REALX_name}/g" src/wp_def_${precision_value}.f90
            # Build the code
            make clean > /dev/null 2>&1
            make ${config_value} COMPILER=${compiler} MACHINE=${machine} PREC=${precision_value} OPTLPKBLS=${use_optimized_lapack_value} EXEFILE=ecg
            # Check if the build was successful
            if [ $? -eq 0 ]; then
              echo "═════════════════════ Build finished succesfully ══════════════════════"
              # Copy the code to the bin directory
              mkdir -p ../${binsubdirname}
              mv ${config_value}/ecg ../${binsubdirname}/${binaryfilename}
              counter_successful_builds=$((counter_successful_builds+1))
            else
              echo "════════════════════════════ Build failed ═════════════════════════════"
              counter_failed_builds=$((counter_failed_builds+1))
            fi
            # Restore the original file src/wp_def_${precision_value}.f90
            mv -f src/wp_def_temporary.f90 src/wp_def_${precision_value}.f90
            # Go back to upper level ecg directory, where the script is located
            cd ../
          done
        done
      done
    done
  done
done

echo ""
echo "════════════════════════ Summary of all builds ════════════════════════"
echo "Total number of attempted builds:  " $counter_attempted_builds
echo "Total number of successful builds: " $counter_successful_builds
echo "Total number of failed builds:     " $counter_failed_builds

exit 0
