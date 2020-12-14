#!/bin/bash

#######################################################################
##                       BASH ENVIRONMENT SETUP                      ##
#######################################################################

# set DYNAMON variable and include executables folder in $PATH
export DYNAMON=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
export PATH=$PATH:${DYNAMON}/exe

# handy function to compile a DYNAMON executable
function cdynamon() {
  # check no argument
  if [ "$1" == "" ]; then
    echo "ERROR: No input arguments. Use -h for help."
    return
  fi
  # read EXE name
  exe_name=$1; shift
  # help
  if [ "$exe_name" == "-h" ] || [ "$exe_name" == "--help" ]; then
    echo "##  Compilation wrapper for DYNAMON"
    echo ""
    echo "USAGE:  cdynamon <exe_name>"
    echo ""
    echo "OPTIONS:"
    echo "  EXEPATH=<absolute-path>"
    echo "  DYNPATH=<absolute-path>"
    echo ""
    echo "Current folder must contain 'nofix_qm.f90' file"
    return
  fi
  # compilation
  make dynamon -s -C $DYNAMON/src EXE="$exe_name" "$@" && echo "Succesful compilation of '$exe_name'"
}
export -f cdynamon
