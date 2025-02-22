<img width="300" height="150" src="./docs/dynamon_logo.svg" align="right" />

# DYNAMO<sup>N</sup>

![Platform](https://img.shields.io/badge/platform-linux-lightgrey.svg)
[![License: GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) \
*A general-purpose script for common calculations with fDynamo*

## Installation

1. Clone the repository in your preferred location. Change directory to the new folder.
2. Compile DYNAMON with a Fortran compiler. Tested with *gfortran*.
3. For a better usability, *source* the file `dynamon.rc` in your *.bashrc*.
   The *$DYNAMON* variable will be set to the installation directory and the *bin*
   folder appended to your *$PATH*

```
git clone https://github.com/boneta/dynamon
cd dynamon
make -C src
source dynamon.rc
```

## Usage
Only compiled once and configured at runtime on every execution! \
For more information, see the [documentation](./docs/README.md).
```
dynamon [.dynn] [[--option arg] ...]
```

Two non-exclusive ways of pass specific calculation options:
- With an [input file](./docs/options_file.md) (.dynn)
- Through [command line arguments](./docs/options_cli.md)

For every system to work with a binary file (.bin) and a [selection file](./docs/sele_file.md) (.dynn) with QM and NOFIX atoms are needed. Every calculation starts from a structure file (.crd).

Additionally, a [PyMOL plug-in](./plugin/README.md) is available to helpful with the system setup.
