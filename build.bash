#!/bin/bash
# This is a script that builds all ECG codes in environments and machines used
# in Bubin Research Group. For usage/help run it without any arguments

# Function to display usage information
usage_print() {
  echo ""
  echo "PROPER USAGE:"
  echo "$0 machine=<machinename> toolchain=<toolchainnames> config=<confignames> code=<codenames> nparticles=<nparticles> precision=<precisions> use_optimized_lapack=<yesnoflag>"
  echo ""
  echo "NOTE:"
  echo "All arguments are optional except nparticles. If multiple values are specified for an argument, they must be separated by a comma."
  echo "The resulting binary files are stored in the directory bin/<toolchainname>/<configname>/. This script must be executed in the directory where it is located (ecg root directory)."
  echo ""
  echo "DESCRIPTION OF ARGUMENTS:"
  echo "<machinename> is the name of the machine. Only a single value may be specified. It could be ubuntu-generic (default), linux-generic, shabyt, muon, puma, ocelote, elgato."
  echo "<toolchainnames> are the names of the toolchains. Supported values for different machines are"
  echo "  ubuntu-generic  ::  foss-2025b, foss-2025a, foss-2024a, foss-2023b, intel-2025b, intel-2025a, intel-2023b, systemdefault (default)"
  echo "  linux-generic   ::  foss-2025b, foss-2025a, foss-2024a, foss-2023b, intel-2025b, intel-2025a, intel-2023b, systemdefault (default)"
  echo "  irgetas         ::  foss-2025b, foss-2025a, foss-2024a, foss-2023b, intel-2025b, intel-2025a, intel-2023b, systemdefault (default)"
  echo "  shabyt          ::  foss-2025b, foss-2025a, foss-2024a, foss-2023b, intel-2025b, intel-2025a, intel-2023b, systemdefault (default)"
  echo "  muon            ::  foss-2025b, foss-2025a, foss-2024a, foss-2023b, intel-2025b, intel-2025a, intel-2023b, systemdefault (default)"
  echo "  puma            ::  gnu8-8.3.0, intel-2020.4, systemdefault (default)"
  echo "  ocelote         ::  gnu8-8.3.0, intel-2020.4, systemdefault (default)"
  echo "  elgato          ::  gnu8-8.3.0, intel-2020.4, systemdefault (default)"
  echo "<confignames> is the list of configurations that need to be built. Currently these could be debug or release. The default includes all configurations."
  echo "<codenames> is the list of codes that need to be built. Currently these could be RG_0S,RG_1P,RG_2P,RG_2D,RG_0S-1P,RG_0S2P,RG_0S-2D,RG_1P-1P,RG_2P-2P,RG_2P-2D. The default includes all codes."
  echo "<nparticles> defines for how many particles each code must be build for. There is no default value. This argument must be present."
  echo "<precisions> is the kind parameter for real type. 8 corresponds to double precision (fp64), 10 corresponds to extended precision (fp80), 16 corresponds to quadruple precision. Different compilers/toolchain support different kinds. For example, Intel compilers supports only 8 and 16, while modern GNU compilers support 8, 10, and 16. The default value is 8."
  echo "<yesnoflag> is a flag that specifies whether to use optimized LAPACK libraries (yes) or use nonoptimized LAPACK from source (no). The default value is yes for precision=8. For precision=10 and precision=16 the value is always no, regardless how the flag is set because optimized LAPACK for these two cases is unavailable."
  echo ""
  echo "EXECUTION EXAMPLES:"
  echo "$0 machine=ubuntu-generic toolchain=foss-2025b config=release code=RG_0S,RG_1P,RG_2D,RG_2P nparticles=3,4,5,6,7,8 precision=8,10,16 use_optimized_lapack=no"
  echo "$0 machine=ubuntu-generic toolchain=foss-2025a config=release code=RG_0S,RG_1P,RG_2D,RG_2P nparticles=3,4,5,6,7,8 precision=8 use_optimized_lapack=yes"
  echo "$0 machine=ubuntu-generic toolchain=foss-2023b config=release code=RG_0S,RG_1P nparticles=3,4,5,6 precision=8,10,16 use_optimized_lapack=no"
  echo "$0 machine=muon toolchain=foss-2023b,intel-2023b config=debug,release code=RG_0S nparticles=4,5 precision=8,16 use_optimized_lapack=yes,no"
  echo "$0 machine=shabyt nparticles=4,5 precision=10"
  echo "$0 nparticles=7"
  exit 1
}

# Set default values for arguments
machine="ubuntu-generic"
toolchain="systemdefault"
config="debug,release"
# code="CG_0S, RG_0S, RG_1P, RG_2D, RG_2P, RG_0S-1P, RG_1P-2D, RG_1P-2P, RG_0S-2D, RG_0S-2P, RG_1P-1P, RG_2D-2D, RG_2P-2D, RG_2P-2P"
code="CG_0S, RG_0S, RG_1P, RG_2D, RG_2P, RG_0S-1P, RG_1P-2D, RG_1P-2P, RG_0S-2D, RG_0S-2P, RG_1P-1P, RG_2D-2D, RG_2P-2D, RG_2P-2P"
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
if [[ ( "$machine" != "ubuntu-generic" ) && ( "$machine" != "linux-generic" ) && ( "$machine" != "shabyt" ) && ( "$machine" != "muon" ) && ( "$machine" != "puma" ) && ( "$machine" != "ocelote" ) && ( "$machine" != "elgato" ) ]]; then
  echo "ERROR, INVALID MACHINE: $machine"
  usage_print
  exit 1
fi

# Chech if config is set properly
for config_value in ${config_list[@]}; do
  if [ "$config_value" != "debug" ] && [ "$config_value" != "release" ]; then
    echo "ERROR, WRONG VALUE(S) OF ARGUMENT: config"
    usage_print
    exit 1
  fi
done

# Check if code is set properly
for code_value in ${code_list[@]}; do
  if [[ ( "$code_value" != "CG_0S" ) && ( "$code_value" != "RG_0S" ) && ( "$code_value" != "RG_1P" ) && ( "$code_value" != "RG_2D" ) && ( "$code_value" != "RG_2P" ) && ( "$code_value" != "RG_0S-1P" ) && ( "$code_value" != "RG_0S-2D" ) && ( "$code_value" != "RG_0S-2P" ) && ( "$code_value" != "RG_1P-1P" ) && ( "$code_value" != "RG_1P-2D" ) && ( "$code_value" != "RG_2D-2D" )  && ( "$code_value" != "RG_2P-2D" ) && ( "$code_value" != "RG_2P-2P" ) ]]; then
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
  if [ "$precision_value" != "8" ] && [ "$precision_value" != "10" ] && [ "$precision_value" != "16" ]; then
    echo "ERROR, WRONG VALUE(S) OF ARGUMENT: precision"
    usage_print
    exit 1
  fi
done

# Check if use_optimized_lapack is set properly
for use_optimized_lapack_value in ${use_optimized_lapack_list[@]}; do
  if [ "$use_optimized_lapack_value" != "yes" ] && [ "$use_optimized_lapack_value" != "no" ]; then
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

for toolchain_value in ${toolchain_list[@]}; do
  # Drop any modules that may have been left loaded and suppress any output of this command
  module reset > /dev/null 2>&1
  # Check if the toolchain is valid for the specified machine and load the corresponding module.
  # If the module is not loaded or/and the proper Fortran compiler and MPI are inaccessible,
  # print an error message.
  if [ "$machine" = "ubuntu-generic" ] || [ "$machine" = "linux-generic" ] || [ "$machine" = "irgetas" ] || [ "$machine" = "shabyt" ] || [ "$machine" = "muon" ]; then
    # Possibly insert $machine in the path given by $bindirname (this can be left empty)
    machinedirname=""
    if [ "$toolchain_value" = "foss-2025b" ]; then
      module_for_toolchain="foss/2025b"
      module load $module_for_toolchain
      compiler_type=gnu
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which gfortran 2>&1 | grep -qi "software/GCCcore/14.3.0/bin/gfortran"; then
        echo "gfortran compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpif90 2>&1 | grep -qi "software/OpenMPI/5.0.8-GCC-14.3.0/bin/mpif90"; then
        echo "mpif90 compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "foss-2025a" ]; then
      module_for_toolchain="foss/2025a"
      module load $module_for_toolchain
      compiler_type=gnu
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which gfortran 2>&1 | grep -qi "software/GCCcore/14.2.0/bin/gfortran"; then
        echo "gfortran compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpif90 2>&1 | grep -qi "software/OpenMPI/5.0.7-GCC-14.2.0/bin/mpif90"; then
        echo "mpif90 compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "foss-2024a" ]; then
      module_for_toolchain="foss/2024a"
      module load $module_for_toolchain
      compiler_type=gnu
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which gfortran 2>&1 | grep -qi "software/GCCcore/13.3.0/bin/gfortran"; then
        echo "gfortran compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpif90 2>&1 | grep -qi "software/OpenMPI/5.0.3-GCC-13.3.0/bin/mpif90"; then
        echo "mpif90 compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "foss-2023b" ]; then
      module_for_toolchain="foss/2023b"
      module load $module_for_toolchain
      compiler_type=gnu
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which gfortran 2>&1 | grep -qi "software/GCCcore/13.2.0/bin/gfortran"; then
        echo "gfortran compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpif90 2>&1 | grep -qi "software/OpenMPI/4.1.6-GCC-13.2.0/bin/mpif90"; then
        echo "mpif90 compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "intel-2025b" ]; then
      module_for_toolchain="intel/2025b"
      compiler_type=intel
      module load $module_for_toolchain
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which ifx 2>&1 | grep -qi "software/intel-compilers/2025.2.0/compiler/2025.2/bin/ifx"; then
        echo "ifort compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpiifort 2>&1 | grep -qi "software/impi/2021.16.1-intel-compilers-2025.2.0/mpi/2021.16/bin/mpiifort"; then
        echo "mpiifort compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "intel-2025a" ]; then
      module_for_toolchain="intel/2025a"
      compiler_type=intel
      module load $module_for_toolchain
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which ifx 2>&1 | grep -qi "software/intel-compilers/2025.1.1/compiler/2025.1/bin/ifx"; then
        echo "ifort compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpiifort 2>&1 | grep -qi "software/impi/2021.15.0-intel-compilers-2025.1.1/mpi/2021.15/bin/mpiifort"; then
        echo "mpiifort compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "intel-2023b" ]; then
      module_for_toolchain="intel/2023b"
      compiler_type=intel
      module load $module_for_toolchain
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which ifort 2>&1 | grep -qi "software/intel-compilers/2023.2.1/compiler/2023.2.1/linux/bin/intel64/ifort"; then
        echo "ifort compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpiifort 2>&1 | grep -qi "software/impi/2021.10.0-intel-compilers-2023.2.1/mpi/2021.10.0/bin/mpiifort"; then
        echo "mpiifort compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "systemdefault" ]; then
      compiler_type=gnu
      if ! which gfortran 2>&1 | grep -qi "/usr/bin/gfortran"; then
        echo "gfortran compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that it is installed on your system."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpif90 2>&1 | grep -qi "/usr/bin/mpif90"; then
        echo "mpif90 compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that it is installed on your system."
        echo "Skipping this toolchain."
        continue
      fi
    else
      echo "ERROR, INVALID TOOLCHAIN $toolchain FOR MACHINE $machine"
      usage_print
      exit 1
    fi
  fi
  if [ "$machine" = "puma" ] || [ "$machine" = "ocelote" ] || [ "$machine" = "elgato" ]; then
    # Possibly insert $machine in the path given by $bindirname (this can be left empty). Be sure to add a trailing slash.
    machinedirname=${machine}"/"
    if [ "$toolchain_value" = "gnu8-8.3.0" ]; then
      module_for_toolchain="gnu8/8.3.0"
      module load $module_for_toolchain
      compiler_type=gnu
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which gfortran 2>&1 | grep -qi "/opt/ohpc/pub/compiler/gcc/8.3.0/bin/gfortran"; then
        echo "gfortran compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpif90 2>&1 | grep -qi "/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.4/bin/mpif90"; then
        echo "mpif90 compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "intel-2020.4" ]; then
      module_for_toolchain="intel/2020.4"
      compiler_type=intel
      module load $module_for_toolchain
      if module load $module_for_toolchain 2>&1 | grep -qi "error"; then
        echo "Module" $module_for_toolchain " for toolchain" $toolchain_value "could not be loaded"
        echo "Either the machine argument you use is incorrect or the lmod configuration has been changed."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which ifort 2>&1 | grep -qi "/opt/ohpc/pub/compiler/intel_2020_u4/compilers_and_libraries/linux/bin/intel64/ifort"; then
        echo "ifort compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpiifort 2>&1 | grep -qi "/opt/ohpc/pub/compiler/intel_2020_u4/compilers_and_libraries/linux/mpi/intel64/bin/mpiifort"; then
        echo "mpiifort compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that module" $module_for_toolchain " is properly configured."
        echo "Skipping this toolchain."
        continue
      fi
    elif [ "$toolchain_value" = "systemdefault" ]; then
      compiler_type=gnu
      if ! which gfortran 2>&1 | grep -qi "/opt/ohpc/pub/compiler/gcc/8.3.0/bin/gfortran"; then
        echo "gfortran compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that it is installed on your system."
        echo "Skipping this toolchain."
        continue
      fi
      if ! which mpif90 2>&1 | grep -qi "/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.4/bin/mpif90"; then
        echo "mpif90 compiler matching toolchain" $toolchain_value" is not accessible."
        echo "Check that it is installed on your system."
        echo "Skipping this toolchain."
        continue
      fi
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
            echo "========================== Starting new build =========================="
            echo "machine="$machine "   toolchain="$toolchain_value "   config="$config_value
            echo "code="$code_value "   nparticles="$nparticles_value "   precision="$precision_value "   use_optimized_lapack="$use_optimized_lapack_value
            echo "-------------------------- Command make output -------------------------"
            # Check if file ${code_value}/src/wp_def_${precision_value}.f90 exists. This way
            # we aslo automtically test if the directory ${code_value} for this specific code exists
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
            sed -i "s/Glob_MaxAllowedNumOfParticles=[0-9]\+/Glob_MaxAllowedNumOfParticles=${precision_value}/g" src/wp_def_${precision_value}.f90
            # Replace "MPI_DPREC=..." with "MPI_DPREC=${MPI_REALX_name}" in file src/wp_def_${precision_value}.f90
            sed -i "s/MPI_DPREC=[^ ][^ ]*/MPI_DPREC=${MPI_REALX_name}/g" src/wp_def_${precision_value}.f90
            # Build the code
            make clean > /dev/null 2>&1
            make ${config_value} COMPILER_TYPE=${compiler_type} MACHINE=${machine} PREC=${precision_value} OPTLPKBLS=${use_optimized_lapack_value} EXEFILE=ecg
            # Check if the build was successful
            if [ $? -eq 0 ]; then
              echo "===================== Build finished succesfully ======================="
              # Copy the code to the bin directory
              mkdir -p ../${binsubdirname}
              mv ${config_value}/ecg ../${binsubdirname}/${binaryfilename}
              counter_successful_builds=$((counter_successful_builds+1))
            else
              echo "============================= Build failed ============================="
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
echo "========================== Summary of all builds ========================="
echo "Total number of attempted builds:  " $counter_attempted_builds
echo "Total number of successful builds: " $counter_successful_builds
echo "Total number of failed builds:     " $counter_failed_builds

exit 0
