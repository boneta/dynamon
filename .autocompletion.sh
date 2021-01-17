#!/bin/bash

#######################################################################
##                        BASH AUTOCOMPLETION                        ##
#######################################################################

_dynamon_compl() {
  local n
  local IFS=$'\n'
  local c=${COMP_WORDS[COMP_CWORD]}
  for ((n=1;n<COMP_CWORD;++n)) ; do
    [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break
  done
  local p=${COMP_WORDS[COMP_CWORD-n]}
  COMPREPLY=()
  if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then
    COMPREPLY=( $(compgen -S ' ' -X '!*.dynn' -f -- $c ; compgen -S '/' -d -- $c ; compgen -S ' ' -W $'-h\n--help\n--MODE\n--NAME\n--SYS\n--BIN\n--SELE\n--COORD\n--FF\n--SEQ\n--CORES\n--MEMORY\n--CHARGE\n--MULTI\n--SEMIEMP\n--GAUSS\n--FUNC\n--BASIS\n--CG_TOLERANCE\n--LBFGSB_TOLERANCE\n--TEMP\n--EQUI\n--PROD\n--VEL\n--LOC_STEPS\n--TS\n--IRC_DIR\n--KIE_SKIP\n--KIE_HESS\n--INT_DCD\n--N\n--DIST' -- $c) )
    return 0
  fi
  case "$p" in
    --MODE) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -W $'BIN\nSP\nMINI\nLOCATE\nIRC\nMD\nINTERACTION\nKIE' -- $c));;
    --SYS) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -W "$(ls $DYNAMON/user/ | awk -F'.' '{print $1}' )" -- $c));;
    --BIN) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -W "$(ls $DYNAMON/user/ | grep .bin )" -- $c; compgen -S ' ' -X '!*.bin' -f -- $c));;
    --SELE) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -W "$(ls $DYNAMON/user/ | grep .dynn )" -- $c; compgen -S ' ' -X '!*.dynn' -f -- $c));;
    --COORD) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.crd' -f -- $c ; compgen -S '/' -d $c));;
    --FF) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ff' -f -- $c ; compgen -S '/' -d $c));;
    --SEQ) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.seq' -f -- $c ; compgen -S '/' -d $c));;

    --SEMIEMP) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -W $'AM1\nRM1\nPM3\nMNDO\nPDDG' -- $c));;
    --GAUSS) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -W $'T\nF' -- $c));;
    --VEL) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.vel' -f -- $c ; compgen -S '/' -d $c));;
    --TS) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -W $'T\nF' -- $c));;
    --IRC_DIR) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -W $'-1\n1' -- $c));;
    --INT_DCD) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.dcd' -f -- $c ; compgen -S '/' -d $c));;
  esac }

complete -o nospace -F _dynamon_compl dynamon
