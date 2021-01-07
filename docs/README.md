# DYNAMO<sup>N</sup>

*A general-purpose script for common calculations with fDynamo*

## Modes
Several main calculation modes are available to be chosen at runtime. It must be
specified in upper cases.

| MODE NAME | ALT NAMES | DESCRIPTION |
| :-------- | :-------- | ----------- |
| SP          | CORR      | Single point calculation |
| MINI        | SCAN, PES | Geometrical optimization |
| LOCATE      | -         | Find a stationary point or TS |
| IRC         | -         | Follow the internal reaction coordinate from a TS |
| MD          | PMF       | BO-Molecular dynamics |
| INTERACTION | -         | Electrostatic and VdW interactions between ligand-protein |
| KIE         | -         | Kinetic isotope effect |

## Calculation options
All the options for every calculation are set at runtime.
There are two ways of doing so: an input file that is passed as first argument or
directly with specific arguments but with limited capabilities. Both methods can
be independently employed or together, in which case the arguments overwrite the file options.
 - [Input file (.dynn)](./options_file.md)
 - [Command line arguments](./options_cli.md)
