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
| --BIN            | *char* | Binary file of system |
| --COORD          | *char* | Coordinates file |

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
|   ARGUMENT NAME  | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| --TEMP           | *real* | Temperature [k] |
| --EQUI           | *int*  | Number of steps of equilibration |
| --PROD           | *int*  | Number of steps of production |
| --VEL            | *char* | Velocities file to read instead of generate random (for continuations) |
| | |
| --TS             | *bool* | Search for a saddle point (transition state) |
| | |
| --IRC_DIR        | *int*  | Initial direction to follow {-1,1} |

#### Constraints
The definition of constraints is necessarily made through an input file, although
some option can be added or modified with command line arguments.

The same number of argument in an option must be used as number of constraints defined,
assigning to them in corresponding order.

|   ARGUMENT NAME  | TYPE   | DESCRIPTION |
| :--------------- | :--:   | ----------- |
| --N              | *int*  | Index number |
| --DIST           | *real* | Distance to constraint directly |
