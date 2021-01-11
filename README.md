<img width="200" height="200" src="./docs/dynamon_logo.svg" align="right" />

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
For more information, see the [documentation](./docs/README.md).
```
dynamon [.dynn] [[--option arg] ...]
```
Only compiled once and configured at runtime on every execution. \
There are two complementary ways of pass specific calculation options: an [input file](./docs/options_file.md) (.dynn)
and through [command line arguments](./docs/options_cli.md). Additionally, for every system to work with a binary file
(.bin) and a [selection file](./docs/sele_file.md) (.dynn) with QM and NOFIX atoms are needed.
*Remark:* Only calculations starting from a single structure can be performed.
