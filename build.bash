#!/bin/bash
# This is a script that builds all ECG codes in environments and machines used 
# in Bubin Research Group. General usage:
#   ./build.bash <machine> <toolchain> 
#
# Examples of invoking this script with arguments:
#   ./build.bash 
#   ./build.bash shabyt gcc9 8 10   <--  gcc9 toolchain, precision 8 and 10
#   ./compile.csh wsl_ubuntu20 gcc 8 10 16   <-- system default gcc, precisions 8, 10, 16
#   ./compile.csh worker26 gcc12 8 10 16   <-- system default gcc, precisions 8, 10, 16
# This script must be executed in the directory one level above
# directories RGL0, RGL1, RGL2P, RGL2D

# ================ This is currently work in progress =======================

# Function to display usage information
usage() {
    echo "Usage:"
    echo "$0 machine=<machinename> toolchain=<toolchainname> code=<codenames> nparticles=<nparticles> precision=<precisions>"
    echo "Description:"
    echo "<machinename> is the name of the machine or OS distro currently supported. It could take the values wsl_ubuntu_22.04 (default), shabyt, muon."
    echo "<toolchain> is the name of the toolchain. Supported values for different machines are"
    echo "    wsl_ubuntu_22.04: foss-2023b (default), intel-2023b"
    echo "    shabyt: foss-2023b (default), intel-2023b"
    echo "    muon: foss-2023b (default), intel-2023b"
    echo "<codenames> is the list of codes (if more than one then separated by a comma) that need to be built. Currently these could be RGL0,RGL1,RGL2P,RGL2D,RGL01,RGL11,RGL02P,RGL2P2P,RGL2P2D. The default includes is all codes."
    echo "<nparticles> defines for how many particles each code must be build for. It can include a list separated by a comma, e.g. 3,4,5 or a similar expression. There is no default value. This argument must be present."
    echo "<precisions> is the kind parameter for real type. 8 corresponds to double precision (fp64), 10 corresponds to extended precision (fp80), 16 corresponds to quadruple precision. Different compilers/toolchain support different kinds. For example, Intel compilers supports only 8 and 16, while modern GNU compilers support 8, 10, and 16. It is possible to specify multiple choices separated by a comma. The default value is 8."
    echo "Examples:"
    echo "$0 machine=wsl_ubuntu_22.04 toolchain=foss-2023b codenames=RGL0,RGL1,RGL2P,RGL2D nparticles=3,4,5,6 precision=8,10,16"
    echo "$0 machine=muon toolchain=intel-2023b codenames=RGL0 nparticles=4,5 precision=8,16"
    echo "$0 nparticles=7"
    exit 1
}

# Parse arguments
for arg in "$@"; do
    case $arg in
        argument1=*)
            IFS=',' read -r -a array1 <<< "${arg#*=}"
            ;;
        argument2=*)
            IFS=',' read -r -a array2 <<< "${arg#*=}"
            ;;
        *)
            echo "Invalid argument: $arg"
            usage
            ;;
    esac
done