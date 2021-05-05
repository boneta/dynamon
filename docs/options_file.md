# Input file
A plain-text file (.dynn) containing the necessary calculation options
that can be passed as first argument. Additional [arguments](options_cli.md)
can also be passed to add or modify them.

### Format specifications
Options must be specified in upper cases, but the arguments not necessarily. \
A line preceded with *'#'* or *'!'* will be treated as a comment and skipped. \
Only one option can be passed per line and all must present an argument, separated
by spaces or tabs. \
The order of the options is irrelevant. \
String arguments that contains slashes or commas should be between quotes.

#### Basics
|   OPTION NAME    | TYPE   | DEFAULT | DESCRIPTION |
| :--------------- | :--:   | :-----: | ----------- |
| MODE             | *char* | -       | Type of calculation |
| | | |
| NAME             | *char* | *mode*  | Base name for files |
| OUT              | *char* | *mode*.out  | File with distances and energies |
| | | |
| SYS              | *char* | -       | System base name (for searching BIN/SELE) |
| BIN              | *char* | -       | Binary file of system (.bin) |
| SELE             | *char* | -       | Selection file for QM and NOFIX (.dynn) |
| COORD            | *char* | -       | Coordinates file |
| FF               | *char* | -       | Force field file (.ff) |
| SEQ              | *char* | -       | Sequence file (.seq) |

#### Computational Settings
|   OPTION NAME    | TYPE   | DEFAULT | DESCRIPTION |
| :--------------- | :--:   | :-----: | ----------- |
| CORES            | *int*  | 1       | Total number of cores (for Gaussian) |
| MEMORY           | *char* | 3000MB  | Total RAM memory (for Gaussian) |

#### QM Region
|   OPTION NAME    | TYPE   | DEFAULT | DESCRIPTION |
| :--------------- | :--:   | :-----: | ----------- |
| CHARGE           | *int*  | 0       | QM-region charge |
| MULTI            | *int*  | 1       | QM-region multiplicity |
| FORCE_UHF        | *bool* | F       | Force an unrestricted calculation |

#### Method
|   OPTION NAME    | TYPE   | DEFAULT     | DESCRIPTION |
| :--------------- | :--:   | :-----:     | ----------- |
| SEMIEMP          | *char* | AM1         | Semi-empirical hamiltonian |
| GAUSS            | *bool* | F           | Use Gaussian software for the qm-region |
| FUNC             | *char* | M062X       | DFT functional (for Gaussian) |
| BASIS            | *char* | 6-31+G(d,p) | DFT basis set (for Gaussian) |

#### Minimization
|   OPTION NAME    | TYPE   | DEFAULT | DESCRIPTION |
| :--------------- | :--:   | :-----: | ----------- |
| CG_STEPS         | *int*  | 10000   | Conjugate-Gradient maximum number of steps |
| CG_TOLERANCE     | *real* | 0.2     | Conjugate-Gradient convergence criteria |
| LBFGSB_STEPS     | *int*  | 1000    | L-BFGSB maximum number of steps |
| LBFGSB_TOLERANCE | *real* | 0.1     | L-BFGSB convergence criteria |

#### Molecular Dynamics
|   OPTION NAME    | TYPE   | DEFAULT | DESCRIPTION |
| :--------------- | :--:   | :-----: | ----------- |
| TEMP             | *real* | 298.0   | Temperature [k] |
| MD_STEP          | *real* | 0.001   | MD time-step [ps] |
| EQUI             | *int*  | 0       | Number of steps of equilibration |
| PROD             | *int*  | 1000    | Number of steps of production |
| DCD_FREQ         | *int*  | 100     | Frequency to save structures to the trajectory file |
| VEL              | *char* | -       | Velocities file to read instead of generate random (for continuations) |

#### PBC
|   OPTION NAME    | TYPE   | DEFAULT | DESCRIPTION |
| :--------------- | :--:   | :-----: | ----------- |
| PBC              | *bool* | F       | Scheme for non-bonding interactions (minimum_image) |

#### Locate
|   OPTION NAME    | TYPE   | DEFAULT | DESCRIPTION |
| :--------------- | :--:   | :-----: | ----------- |
| LOC_STEPS        | *int*  | 100     | Baker search maximum number of steps |
| LOC_TOLERANCE    | *real* | 1.0     | Convergence criteria for location |
| TS               | *bool* | F       | Search for a saddle point (transition state) |

#### IRC
|   OPTION NAME    | TYPE   | DEFAULT | DESCRIPTION |
| :--------------- | :--:   | :-----: | ----------- |
| IRC_DIR          | *int*  | -       | Initial direction to follow {-1,1} |
| IRC_STEPS        | *int*  | 400     | IRC maximum number of steps |
| IRC_DSP          | *real* | 0.01    | Displacement on every step |

#### Interaction
|   OPTION NAME    | TYPE   | DEFAULT      | DESCRIPTION |
| :--------------- | :--:   | :-----:      | ----------- |
| DCD_STRIDE       | *int*  | 1            | Read only every n-th frame of the trajectory |
| INT_DCD          | *char* | -            | Trajectory file along which calculate interactions |
| INT_WBOX         | *char* | BOX, SOL     | Subsystem name to be considered water box (added to def) |
| INT_IONS         | *char* | IONS, NA, CL | Subsystem name to be considered ions (added to def) |

#### KIE
|   OPTION NAME    | TYPE                | DEFAULT     | DESCRIPTION |
| :--------------- | :--:                | :-----:     | ----------- |
| KIE_ATOM         | *char* *int* *char* | -           | Atom to calculate KIE (subsystem, residue number, atom name) |
| KIE_SKIP         | *int*               | 0           | Number of frequencies to skip |
| KIE_MASS         | *real*              | 2.01410     | New mass for the atom |
| KIE_HESS         | *char*              | update.dump | Hessian file |

#### Constraints
An arbitrary number of constraints can be applied on almost every kind of calculation. \
The atoms involved must be independently defined and every constraint placed in
its own block with specific options.

An atom is defined with the option **ATOM** or **A** followed by three arguments: chain letter,
number of residue and atom name. Each atom is labeled with a number in order of appearance, starting from 1.

Every constraint block is delimited by the options **CONSTR** or **C**. The first must be followed
by an argument with the type of constraint but the last does not need argument.
Between them the options for each constraint are specified in the same way as previously explained.

The distances for the constraint are taken from the coordinates file (.crd) if DCRD is set *True*.
If not, it would be directly set through DIST argument and if absent, it is calculated taking into
account N, DINIT and STEP options.

##### Constraint types
| TYPE NAME         | ALT NAME |
| :--------         | :------- |
| ANGLE             | -        |
| DIHEDRAL          | -        |
| DISTANCE          | D, d     |
| >DISTANCE         | >D, >d   |
| <DISTANCE         | <d, <d   |
| MULTIPLE_DISTANCE | M, m     |

##### Constraint options
|   OPTION NAME    | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| SYMM             | *int*  | Symmetry for multiple distance {-1,1} |
| N                | *int*  | Index number |
| ATOMS            | *int*  | Atoms numbers corresponding to definition order |
| FORCE            | *real* | Force to apply [J·A] |
| DCRD             | *bool* | Read distance from coordinates file |
| DINIT            | *real* | Initial distance to consider |
| DIST             | *real* | Distance to constraint directly |
| STEP             | *real* | Step distance to consider |
| DFILE            | *char* | File name to write distance evolution (dat_) |
