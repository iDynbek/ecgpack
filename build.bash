#!/bin/bash
# This is a script that builds all ECG codes in environments and machines used 
# in Bubin Research Group. For usage/help run it without any arguments

# Function to display usage information
usage_print() {
    echo ""
    echo "PROPER USAGE:"
    echo "$0 machine=<machinename> toolchain=<toolchainname> config=<confignames> code=<codenames> nparticles=<nparticles> precision=<precisions> use_optimized_lapack=<yesnoflag>"
    echo ""   
    echo "NOTE:" 
    echo "All arguments optional except nparticles. If multiple values are specified for an argument, they must be separated by a comma."
    echo "" 
    echo "DESCRIPTION:"
    echo "<machinename> is the name of the machine or OS distro currently supported. Only a single value may be specified. It could be wsl-ubuntu-22.04 (default), shabyt, muon."
    echo "<toolchain> are the names of the toolchains. Supported values for different machines are"
    echo "    wsl-ubuntu-22.04: foss-2023b (default), intel-2023b"
    echo "    shabyt: foss-2023b (default), intel-2023b"
    echo "    muon: foss-2023b (default), intel-2023b"
    echo "<confignames> is the list of configurations that need to be built. Currently these could be debug or release. The default includes all configurations."
    echo "<codenames> is the list of codes that need to be built. Currently these could be RGL0,RGL1,RGL2P,RGL2D,RGL01,RGL02P,RGL02D,RGL11,RGL2P2P,RGL2P2D. The default includes is all codes."
    echo "<nparticles> defines for how many particles each code must be build for. It can include a list separated by a comma, e.g. 3,4,5 or a similar expression. There is no default value. This argument must be present."
    echo "<precisions> is the kind parameter for real type. 8 corresponds to double precision (fp64), 10 corresponds to extended precision (fp80), 16 corresponds to quadruple precision. Different compilers/toolchain support different kinds. For example, Intel compilers supports only 8 and 16, while modern GNU compilers support 8, 10, and 16. It is possible to specify multiple choices separated by a comma. The default value is 8."
    echo "<yesnoflag> is a flag that specifies whether to use optimized LAPACK libraries (yes) or use nonoptimized LAPACK from source (no). The default value is yes for precision=8. For precision=10 and precision=16 the value is always no, regardless how the flag is set because optimized LAPACK for these two cases is unavailable."
    echo ""
    echo "EXECUTION EXAMPLES:"
    echo "$0 machine=wsl-ubuntu-22.04 toolchain=foss-2023b config=release code=RGL0,RGL1,RGL2P,RGL2D nparticles=3,4,5,6 precision=8,10,16 use_optimized_lapack=no"
    echo "$0 machine=muon toolchain=foss-2023b,intel-2023b config=debug,release code=RGL0 nparticles=4,5 precision=8,16 use_optimized_lapack=yes,no"
    echo "$0 machine=shabyt nparticles=4,5 precision=10"
    echo "$0 nparticles=7"
    exit 1
}

# Set default values for arguments
machine="wsl-ubuntu-22.04"
toolchain="foss-2023b"
config="debug,release"
code="RGL0,RGL1,RGL2P,RGL2D,RGL01,RGL02P,RGL02D,RGL11,RGL2P2P,RGL2P2D"
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
    if [[ ( "$code_value" != "RGL0" ) && ( "$code_value" != "RGL1" ) && ( "$code_value" != "RGL2P" ) && ( "$code_value" != "RGL2D" ) && ( "$code_value" != "RGL01" ) && ( "$code_value" != "RGL02P" ) && ( "$code_value" != "RGL02D" ) && ( "$code_value" != "RGL11" ) && ( "$code_value" != "RGL2P2P" ) && ( "$code_value" != "RGL2P2D" ) ]]; then
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

# Drop any modules that may have been left loaded and suppress any output of this command
module reset > /dev/null 2>&1

# Execute a branch corresponding to a specific machine
if [ "$machine" = "wsl-ubuntu-22.04" ]; then
    # Loop over all toolchains
    for toolchain_value in ${toolchain_list[@]}; do
        # Check if the toolchain is valid for wsl-ubuntu-22.04 and load the corresponding module
        if [ "$toolchain" = "foss-2023b" ]; then
            module load foss/2023b
        elif [ "$toolchain" = "intel-2023b" ]; then
            module load intel/2023b
        else
            echo "ERROR, INVALID TOOLCHAIN $toolchain FOR MACHINE $machine"
            usage_print
            exit 1
        fi
        # Loop over all configurations
        for config_value in ${config_list[@]}; do
            # Loop over all codes
            for code_value in ${code_list[@]}; do
                # Loop over all nparticles
                for nparticles_value in ${nparticles_list[@]}; do
                    # Loop over all precisions
                    for precision_value in ${precision_list[@]}; do
                        # Loop over all values of use_optimized_lapack, but for precision=10,16 only `no` is allowed
                        for use_optimized_lapack_value in ${use_optimized_lapack_list[@]}; do
                            if [[ "$use_optimized_lapack_value" = "yes" && "$precision_value" = "10"  ]]; then
                                continue
                            fi
                            if [[ "$use_optimized_lapack_value" = "yes" && "$precision_value" = "16"  ]]; then
                                continue
                            fi                            
                            echo "========================== New build started ==========================="
                            echo "machine="$machine "   toolchain="$toolchain_value "   config="$config_value
                            echo "code="$code_value "   nparticles="$nparticles_value "   precision="$precision_value "   use_optimized_lapack="$use_optimized_lapack_value
                            echo "------------------------------------------------------------------------"
                            binsubdirname=${bindirname}/${toolchain_value}/${config_value}
                            #mkdir -p ${binsubdirname}
                            binaryfilename=${code_value}_N${nparticles_value}_P${precision_value}
                            if [[ "$use_optimized_lapack_value" = "yes" && precision_value = "8" ]]; then
                                    binaryfilename=${binaryfilename}_optlapack
                            fi
                            #cd ${code_value}
                            #make clean
                            #make clean
                            #cd ../
                            echo "============================ build finished ============================"
                        done
                    done
                done
            done
        done
        if [ "$toolchain" = "foss-2023b" ]; then
            module unload foss/2023b
        elif [ "$toolchain" = "intel-2023b" ]; then
            module unload intel/2023b
        else  
            echo "WARNING: TOOLCHAIN" $toolchain_value "IS NOT SUPPORTED FOR" $machine
            echo "SKIPPING THIS TOOLCHAIN"
        fi      
    done
elif [ "$machine" = "shabyt" ]; then
    echo "This is work in progress"
elif [ "$machine" = "muon" ]; then
    echo "This is work in progress"
else
    echo "ERROR, INVALID MACHINE: $machine"
    usage_print
    exit 1
fi
