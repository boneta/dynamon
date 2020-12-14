!===============================================================================
!
!        fDynamo v2.2 - a program for performing molecular simulations.
!                    Copyright (C) 2005-2007 Martin J. Field
!
!===============================================================================
!
!       This program is free software; you can redistribute it and/or
!       modify it under the terms of the GNU General Public License
!       as published by the Free Software Foundation; either version 2
!       of the License, or (at your option) any later version.
!
!       This program is distributed in the hope that it will be useful,
!       but WITHOUT ANY WARRANTY; without even the implied warranty of
!       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!       GNU General Public License for more details.
!
!       You should have received a copy of the GNU General Public License
!       along with this program; if not, write to the Free Software
!       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!       MA  02110-1301, USA.
!
!===============================================================================
!
!  Email: martin.field@ibs.fr
!  WWWs:  http://www.ibs.fr and http://www.pdynamo.org
!
!===============================================================================
!                               The DYNAMO Module
!===============================================================================
!
! . Subroutines:
!
!   DYNAMO_FOOTER                   Clean up at the end of a job.
!   DYNAMO_HEADER                   Start up a job.
!
!   This module contains a list of all the DYNAMO program modules. Inclusion of
!   the statement "USE DYNAMO" at the top of a program will allow access to all
!   modules in the DYNAMO library.
!
!===============================================================================
MODULE DYNAMO

! . Module declarations (in order of function and dependency).
USE DEFINITIONS

USE CONSTANTS
USE ELEMENTS
USE IO_UNITS
USE PRINTING
USE STRING
USE TIME

USE FILES
USE PARSING

USE CONJUGATE_GRADIENT
USE DIAGONALIZATION
USE LINEAR_ALGEBRA
USE LINEAR_EQUATIONS
USE RANDOM_NUMBERS
USE SORT
USE SPECIAL_FUNCTIONS
USE STATISTICS

USE BAKER_OPTIMIZATION

USE ATOMS
USE SEQUENCE
USE SYMMETRY

USE ATOM_MANIPULATION

USE CONNECTIVITY
USE GEOMETRY
USE SOLVENT_ACCESSIBLE_SURFACE
USE TRANSFORMATION

USE SUPERIMPOSE

USE COORDINATE_IO
USE PDB_IO
USE XYZ_IO

USE ZMATRIX
USE ZMATRIX_IO

USE CONSTRAINT
#ifdef	HACK_SPLN
USE SPLINE
#endif
USE WHAM

USE MM_FILE_DATA
USE MM_FILE_IO

USE MM_TERMS

USE MM_SYSTEM
USE MM_SYSTEM_EDIT
USE MM_SYSTEM_IO

USE BUILD_COORDINATES

USE ENERGY_COVALENT
USE ENERGY_NON_BONDING

USE MOPAC_DATA
USE MOPAC_PARAMETERS

USE HESSIAN_UPDATE

USE GAUSSIAN_BASIS

USE MOPAC_DENSITY
USE MOPAC_FOCK_MATRIX
USE MOPAC_HAMILTONIAN
USE MOPAC_INTEGRALS
USE MOPAC_SCF

USE MOPAC_ANALYSIS
USE MOPAC_GRADIENTS

USE QUANTUM_ENERGY
USE QUANTUM_PROPERTIES

USE ABINITIO

#ifdef	HACK_GMS
USE GAMESS
#endif

#ifdef	HACK_HTSC
USE HTSC
#endif

#ifdef	HACK_URC
USE URC
#endif

#ifdef	HACK_VRC
USE VRC
#endif

#ifdef	HACK_ZBND
USE ZBND
#endif

USE POTENTIAL_ENERGY
USE NUMERICAL_DERIVATIVES

USE DCD_IO
USE MULTIPOLES
USE NORMAL_MODE_UTILITIES

USE OPTIMIZE_BFGS
USE OPTIMIZE_COORDINATES
USE REACTION_PATH
USE SELF_AVOIDING_WALK
USE NEB_CLIMB
USE NEB_SPLINE

USE DCD_ANALYSIS
USE NORMAL_MODE
USE THERMODYNAMICS_RRHO

USE VELOCITY

USE DYNAMICS_UTILITIES

USE DYNAMICS_LANGEVIN_VERLET
USE DYNAMICS_LEAPFROG_VERLET
USE DYNAMICS_VELOCITY_VERLET

USE MONTE_CARLO_ENERGY
USE MONTE_CARLO_SIMULATION

USE QUICKDOCK

USE AMBERFILE_IO

USE LBFGSB_WRAPPER

#ifdef	HACK_KIE
USE KIE
#endif

! . Remaining statements.
IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------
   SUBROUTINE DYNAMO_FOOTER
   !-----------------------

#ifdef F95
   ! . Print out the CPU time.
   CALL TIME_CPU_PRINT
#endif

   ! . Print out the time.
   CALL TIME_PRINT

   ! . Finish printing.
   CALL PRINT_STOP

   END SUBROUTINE DYNAMO_FOOTER

   !---------------------------------
   SUBROUTINE DYNAMO_HEADER ( TITLE )
   !---------------------------------

   ! . Arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: TITLE

   ! . Start printing.
   IF ( PRESENT ( TITLE ) ) THEN
      CALL PRINT_START ( TITLE )
   ELSE
      CALL PRINT_START ( "Dynamo Output" )
   END IF

   ! . Write out the print line.
   WRITE ( PRINT_LINE, "(A,F3.1,A)" ) "DYNAMO Module Library (Version ", VERSION, ")"

   ! . Print the line as a header.
   CALL PRINT_HEADING ( "olive", "H1", 2, 50, 80, .FALSE. )

   ! . Print out the time.
   CALL TIME_PRINT

#ifdef F95
   ! . Initialize the CPU time counter.
   CALL TIME_CPU_INITIALIZE
#endif

   END SUBROUTINE DYNAMO_HEADER

   !---------------------------
   SUBROUTINE DYNAMO_INITIALIZE
   !---------------------------

   ! . Initialize all data structures to do with system definitions.
   CALL    ATOMS_INITIALIZE
   CALL  MM_FILE_INITIALIZE
   CALL MM_TERMS_INITIALIZE
   CALL SEQUENCE_INITIALIZE
   CALL SYMMETRY_INITIALIZE

   ! . Energy definitions.
   CALL ENERGY_INITIALIZE
   CALL ENERGY_NON_BONDING_OPTIONS ! . Should be replaced by ENB_INITIALIZE.

   END SUBROUTINE DYNAMO_INITIALIZE

END MODULE DYNAMO
