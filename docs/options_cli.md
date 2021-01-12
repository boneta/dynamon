# Command line arguments
Some options can be passed as arguments while invoking the program.
In case of conflict with the [input file](./options_file.md), it will be overwritten.
Upper cases must be used for the option names but not necessary in
the following argument.

#### Basics
|   ARGUMENT NAME  | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| --MODE           | *char* | Type of calculation |
| --NAME           | *char* | Base name for files |
| --SYS            | *char* | System base name (for searching BIN/SELE) |
| --BIN            | *char* | Binary file of system (.bin) |
| --SELE           | *char* | Selection file for QM and NOFIX (.dynn) |
| --COORD          | *char* | Coordinates file |
| --FF             | *char* | Force field file (.ff) |
| --SEQ            | *char* | Sequence file (.seq) |

#### Computational settings
|   ARGUMENT NAME  | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| --CORES          | *int*  | Total number of cores (for Gaussian) |
| --MEMORY         | *char* | Total RAM memory (for Gaussian) |

#### QM region
|   ARGUMENT NAME  | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| --CHARGE         | *int*  | QM-region charge |
| --MULTI          | *int*  | QM-region multiplicity |

#### Method
|   ARGUMENT NAME  | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| --SEMIEMP        | *char* | Semi-empirical hamiltonian |
| --GAUSS          | *bool* | Use Gaussian software for the qm-region |
| --FUNC           | *char* | DFT functional (for Gaussian) |
| --BASIS          | *char* | DFT basis set (for Gaussian) |

#### Mode-specific
|   ARGUMENT NAME    | TYPE   | DESCRIPTION |
| :--------------    | :--:   | ----------- |
| --CG_TOLERANCE     | *real* | Conjugate-Gradient convergence criteria |
| --LBFGSB_TOLERANCE | *real* | L-BFGSB convergence criteria |
| | |
| --TEMP             | *real* | Temperature [k] |
| --EQUI             | *int*  | Number of steps of equilibration |
| --PROD             | *int*  | Number of steps of production |
| --VEL              | *char* | Velocities file to read instead of generate random (for continuations) |
| | |
| --LOC_STEPS        | *int*  | Baker search maximum number of steps |
| --TS               | *bool* | Search for a saddle point (transition state) |
| | |
| --IRC_DIR          | *int*  | Initial direction to follow {-1,1} |
| | |
| --KIE_SKIP         | *int*  | Number of frequencies to skip |
| --KIE_HESS         | *char* | Hessian file |
| | |
| --INT_DCD          | *char* | Trajectory file along which calculate interactions |

#### Constraints
The definition of constraints is necessarily made through an input file, although
some option can be added or modified with command line arguments.

The same number of argument in an option must be used as number of constraints defined,
assigning to them in corresponding order.

|   ARGUMENT NAME  | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| --N              | *int*  | Index number |
| --DIST           | *real* | Distance to constraint directly |
