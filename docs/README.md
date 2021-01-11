# DYNAMO<sup>N</sup>

*A general-purpose script for common calculations with fDynamo*

## Modes
Several main calculation modes are available to be chosen at runtime. Must be
specified in upper cases. Composed calculation such as SCAN/PES/PMF have to be launched individually.

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

#### Selection file
A [selection file (.dynn)](./sele_file.md) (.dynn) is also needed on every calculation
to provide the QM atoms and NOFIX atoms. Usually prepared once for a system.

#### Binary file
A binary file (.bin) which related the coordinates and topology and created with
the fDynamo libraries is necessary for every system.
