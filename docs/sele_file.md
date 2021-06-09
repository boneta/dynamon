# Selection file
A plain-text file (.dynn), following the [input file](./options_file.md#format-specifications)
format, that is used to specify a selection, such as the QM atoms or NOFIX atoms.

#### Selection blocks
A selection is specified in a delimited block, beginning with the keyword **SELECTION** and the name of
the selection (*QM* or *NOFIX*). It must be ended with another **SELECTION**.

Inside, the available hierarchical options are: *subsystem*, *residue number* and *atom name*. \
The are selected with the corresponding option and followed by its argument. \
If a high-level stratum of the hierarchy is selected but not inferior levels are specified, the whole is taken.

|   OPTION NAME    | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| S                | *char* | Subsystem name (chain) |
| R                | *int*  | Residue number |
| A                | *char* | Atom name |

#### Calculation options
The same options as in the [input file](./options_file.md) can be specified as defaults
for the system (i.e. CHARGE, MULTI). But must be keep in mind that this selection file
is the last to be read and therefore will always overwrite the proper input file or
command line arguments.
