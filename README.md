# DYNAMO<sup>N</sup>

![Platform](https://img.shields.io/badge/platform-linux-lightgrey.svg)
[![License: GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


*A general-purpose script for common calculations with fDynamo*

## Installation
No installation needed.\
For a better usability, *source* the file `.dynamon.sh` in your *.bashrc*.
The *$DYNAMON* variable will be set to its directory.

## Usage
```
cdynamon <exe>
<exe> <options.dynn>
```

An executable is compiled only once for every system to be worked on,
defined by the fixed residues and atoms included in the QM region.
This definition is stated in a file named `nofix_qm.f90`.
*fDynamo* will be compiled the first time a DYNAMO<sup>N</sup> compilation is invoked.

The wrapper `cdynamon` for a friendly compilation is available after sourcing.
The executable generated is stored by default at *DYNAMON* and included
in the execution path.

Any calculation is configured at runtime by an options file (.dynn) passed
as first argument.

*Remark:* Only calculations starting from a single structure can be performed.
Composed calculation such as SCAN/PES/PMF have to be launched individually.

#### Manual compilation
`make -C $DYNAMON/src EXE=<exe_name> [EXEPATH=<absolute-path>]`

## Options
Mandatory: MODE | SYSBIN | COORD

Available modes:
  - SP
  - MINIMIZATION
  - LOCATE
  - IRC
  - DYNAMIC
  - INTERACTION
  - KIE

## Example
*Single point calculation in human beta-defensin 1*
> cd $DYNAMON/example \
> cdynamon hbdef1 \
> hbdef1 options.dynn > hbdef1-sp.log
