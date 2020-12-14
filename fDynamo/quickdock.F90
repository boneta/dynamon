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
!                              The QuickDock Module
!===============================================================================
!
! . Procedures (public):
!
!   QUICKDOCK_DISTANCE_MATRIX     Calculate a quickdock distance matrix for a
!                                 set of structures.
!   QUICKDOCK_GRID_SETUP          Set up the interaction grids.
!   QUICKDOCK_OPTIMIZE            Perform a series of quickdock optimizations.
!   QUICKDOCK_OPTIONS             Define the quickdock options.
!   QUICKDOCK_RANDOMIZE_REFERENCE Randomize the reference structure.
!   QUICKDOCK_SETUP               Setup the quickdock calculation.
!   QUICKDOCK_TRAJECTORY          Generate a trajectory of optimized structures.
!
! . Procedures (private):
!
!   QUICKDOCK_COORDINATES_LIGAND  Generate the ligand coordinates.
!   QUICKDOCK_ENERGY              Calculate the energy of a configuration.
!   QUICKDOCK_ENERGY_DIHEDRAL     Calculate the dihedral energy.
!   QUICKDOCK_INTERACTIONS        Calculate the ligand non-bonding energy.
!
!   QUICKDOCK_GRID_GENERATE       Generate the interaction grid.
!   QUICKDOCK_GRID_READ           Read in the interaction grids.
!   QUICKDOCK_GRID_WRITE          Write out the interaction grids.
!
!   QUICKDOCK_INITIALIZE_GRID     Initialize the quickdock grid data.
!   QUICKDOCK_INITIALIZE_OPTIONS  Initialize the quickdock options.
!   QUICKDOCK_INITIALIZE_SYSTEM   Initialize the quickdock system data.
!
!   STRUCTURES_ORDER              Order structures by energy (lowest first).
!   STRUCTURES_READ               Read in a set of structures from a file.
!   STRUCTURES_REDUCE             Remove structures within each other's
!                                 excluded range.
!   STRUCTURES_WRITE              Write out a set of structures to a file.
!
!   VARIABLES_COMPARE             Compare sets of variables.
!   VARIABLES_DISTANCE            Get the distance between two structures.
!   VARIABLES_RANGE_INITIALIZE    Initialize the variable search range.
!
! . Notes:
!
!   This module implements a quick (!) docking strategy that is based upon an
!   MM energy function.
!
!   The search ranges for dihedrals, rotations and translations are symmetric.
!
!===============================================================================
MODULE QUICKDOCK

! . Module declarations.
USE CONSTANTS,       ONLY : ELECT_CONST, GAS_CONSTANT => R, PI, TO_RADIANS
USE DEFINITIONS,     ONLY : DP
USE FILES,           ONLY : NEXT_UNIT
USE IO_UNITS,        ONLY : OUTPUT
USE LINEAR_ALGEBRA,  ONLY : NORMALIZE, ROTATION_MATRIX_AXIS, ROTATION_MATRIX_XYZ
USE PRINTING,        ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH, PRINT_PARAGRAPH_START, PRINT_PARAGRAPH_STOP, &
                            PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, &
			    PRINT_TEXT
USE RANDOM_NUMBERS,  ONLY : RANDOM, RANDOM_VECTOR
USE SORT,            ONLY : SORT_REAL_RANK

USE ATOMS,           ONLY : ATMCRD, NATOMS
USE CONNECTIVITY,    ONLY : CONNECTIVITY_ANGLES, CONNECTIVITY_DIHEDRALS, ORDER_BONDS
USE DCD_IO,          ONLY : DCD_ACTIVATE_WRITE, DCD_DEACTIVATE, DCD_INITIALIZE, DCD_TYPE, DCD_WRITE
USE GEOMETRY,        ONLY : GEOMETRY_ANGLES, GEOMETRY_DIHEDRALS, GEOMETRY_DISTANCES
USE MM_TERMS,        ONLY : ATMCHG, ATMCHG14, ATMEPS, ATMEPS14, ATMEXCI, ATMEXCJ, ATME14I, ATME14J, ATMSIG, &
                            BONDS, DIHEDRAL_TERM, SYSDIH => DIHEDRALS, NBONDS, NSYSDIH => NDIHEDRALS

IMPLICIT NONE
PRIVATE
PUBLIC :: QUICKDOCK_GRID_SETUP, QUICKDOCK_OPTIMIZE, QUICKDOCK_OPTIONS, QUICKDOCK_RANDOMIZE_REFERENCE, &
          QUICKDOCK_SETUP, QUICKDOCK_TRAJECTORY
#ifndef PGPC
SAVE
#endif

!-------------------------------------------------------------------------------
! . Quickdock options.
!-------------------------------------------------------------------------------
! . Scalars.
LOGICAL            :: QDEBUGM = .FALSE., QDEBUGQ  = .FALSE., QEXCLUDE = .FALSE., QREDUCE = .FALSE.
REAL ( KIND = DP ) :: EOFFGRID = 60.0E+6_DP, EPSILON = 4.0_DP, ERECMAX = 10.0_DP, EXCLUDED_DIHEDRAL = 15.0_DP * TO_RADIANS, &
                      EXCLUDED_ROTATION = 7.5_DP * TO_RADIANS, EXCLUDED_TRANSLATION = 0.5_DP, MINDIST2 = 0.01_DP**2,        &
		      RECOMBINATIONFACTOR = 0.1_DP, REDUCEFACTOR = 3.0_DP

!-------------------------------------------------------------------------------
! . System variables.
!-------------------------------------------------------------------------------
! . Scalars.
INTEGER            :: NDATOMS = 0, NDBONDS = 0, NxDIHEDRALS = 0, NLIGNB = 0, NLIGNB14 = 0, NRATOMS = 0, &
                      NROTATABLE = 0, NVARIABLES = 0, REFATM = 0
LOGICAL            :: QROTATION = .FALSE., QTRANSLATION = .FALSE.
REAL ( KIND = DP ) :: TRANSLATION_RANGE = 0.0_DP

! . Allocatable arrays.
INTEGER,             ALLOCATABLE, DIMENSION(:)     :: DATMIND, LISTI, LISTI14, LISTJ, LISTJ14
INTEGER,             ALLOCATABLE, DIMENSION(:,:)   :: DBONDS, ROTATABLE
LOGICAL,             ALLOCATABLE, DIMENSION(:)     :: QDFIXED, QRFIXED, QROTATABLE, QTRANSLATABLE
LOGICAL,             ALLOCATABLE, DIMENSION(:,:)   :: QDIHEDRAL
REAL ( KIND = DP ),  ALLOCATABLE, DIMENSION(:)     :: DATMCHG, DATMCHG14, DATMREP, DATMREP14, DATMVDW, DATMVDW14
REAL ( KIND = DP ),  ALLOCATABLE, DIMENSION(:,:)   :: DATMCRD, DREFCRD, REFCRD
TYPE(DIHEDRAL_TERM), ALLOCATABLE, DIMENSION(:)     :: xDIHEDRALS

!-------------------------------------------------------------------------------
! . System variables (for debugging).
!-------------------------------------------------------------------------------
! . Scalars.
INTEGER :: NDANGLES = 0, NDDIHEDRALS = 0

! . Allocatable arrays.
INTEGER,                          DIMENSION(:,:), POINTER :: DANGLES => NULL ( ), DDIHEDRALS => NULL ( )
REAL ( KIND = DP ),  ALLOCATABLE, DIMENSION(:)            :: REFANG, REFBND, REFDIH

!-------------------------------------------------------------------------------
! . Grid variables.
!-------------------------------------------------------------------------------
! . Scalars.
LOGICAL            :: QGRID       = .FALSE.
REAL ( KIND = DP ) :: GRIDSPACING = 0.2_DP

! . Allocatable arrays.
REAL ( KIND = DP ),  ALLOCATABLE, DIMENSION(:,:,:) :: GRIDQ, GRIDR, GRIDV

! . Fixed arrays.
INTEGER,            DIMENSION(1:3) :: NGRID
REAL ( KIND = DP ), DIMENSION(1:3) :: GRIDEXTENT, GRIDORIGIN

!==============================================================================
CONTAINS
!==============================================================================

   !----------------------------------------------------------------------------
   SUBROUTINE QUICKDOCK_GRID_SETUP ( NPOINTS, CENTER, SPACING, INFILE, OUTFILE )
   !----------------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: INFILE, OUTFILE
   REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: SPACING

   ! . Array arguments.
   INTEGER,            DIMENSION(1:3), OPTIONAL :: NPOINTS
   REAL ( KIND = DP ), DIMENSION(1:3), OPTIONAL :: CENTER

   ! . Initialization.
   CALL QUICKDOCK_INITIALIZE_GRID

   ! . Check that a system exists.
   IF ( NRATOMS <= 0 ) CALL PRINT_ERROR ( "QUICKDOCK_GRID_SETUP", "The reference system has not been defined." )

   ! . Read in a previously generated grid.
   IF ( PRESENT ( INFILE ) ) THEN

      ! . Check that no other options are present.
      IF ( PRESENT ( NPOINTS ) .OR. PRESENT ( CENTER ) .OR. PRESENT ( SPACING ) .OR. PRESENT ( OUTFILE ) ) THEN
         CALL PRINT_ERROR ( "QUICKDOCK_GRID_SETUP", "Other arguments present with INFILE." )
      END IF

      ! . Read the grid.
      CALL QUICKDOCK_GRID_READ ( INFILE )

   ! . Generate a new grid.
   ELSE

      ! . Check for CENTER and NPOINTS which must be present.
      IF ( .NOT. PRESENT ( CENTER  ) ) CALL PRINT_ERROR ( "QUICKDOCK_GRID_SETUP", "CENTER argument missing."  )
      IF ( .NOT. PRESENT ( NPOINTS ) ) CALL PRINT_ERROR ( "QUICKDOCK_GRID_SETUP", "NPOINTS argument missing." )

      ! . Get the grid-spacing.
      IF ( PRESENT ( SPACING ) ) GRIDSPACING = SPACING

      ! . Set NGRID.
      NGRID = NPOINTS

      ! . Set GRIDEXTENT.
      GRIDEXTENT = GRIDSPACING * REAL ( NGRID - 1, DP )

      ! . Set GRIDORIGIN.
      GRIDORIGIN = CENTER - 0.5_DP * GRIDEXTENT

      ! . Allocate the grid arrays.
      ALLOCATE ( GRIDQ(1:NGRID(1),1:NGRID(2),1:NGRID(3)), &
                 GRIDR(1:NGRID(1),1:NGRID(2),1:NGRID(3)), &
		 GRIDV(1:NGRID(1),1:NGRID(2),1:NGRID(3))  )

      ! . Generate the grid.
      CALL QUICKDOCK_GRID_GENERATE

      ! . Do some printing.
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "FF0000" )
      CALL PRINT_SUMMARY_START ( "Quick-Dock Grid Setup" )
      WRITE ( PRINT_LINE, "(I14)"   ) PRODUCT ( NGRID ) ; CALL PRINT_SUMMARY_ELEMENT ( "Grid Points"  )
      WRITE ( PRINT_LINE, "(G14.8)" ) GRIDSPACING       ; CALL PRINT_SUMMARY_ELEMENT ( "Grid Spacing" )
      CALL PRINT_SUMMARY_STOP

      ! . Write out the grid if necessary.
      IF ( PRESENT ( OUTFILE ) ) CALL QUICKDOCK_GRID_WRITE ( OUTFILE )

   END IF

   END SUBROUTINE QUICKDOCK_GRID_SETUP

   !------------------------------------------------------------------------------------------------------
   SUBROUTINE QUICKDOCK_OPTIMIZE ( NOPTIMIZATIONS, NSEARCH, NHUNT, NFINE, NEW_STRUCTURES, OLD_STRUCTURES )
   !------------------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: NEW_STRUCTURES
   INTEGER,               INTENT(IN) :: NFINE, NHUNT, NOPTIMIZATIONS, NSEARCH

   ! . Optional arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: OLD_STRUCTURES

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: VARTOL = 1.0E-2_DP

   ! . Local scalars.
   INTEGER            :: I, IFINE, IHUNT, IOPT, ISEARCH, NDONE, NSTROLD, NSTRUCTURES
   LOGICAL            :: QSKIP
   REAL ( KIND = DP ) :: E, EBEST, EREFERENCE

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES) :: VARIABLES, VARIABLES_BEST, VARIABLES_DEBUG, VARIABLES_RANGE, &
                                                  VARIABLES_RANGE_EXCLUDED
   REAL ( KIND = DP ), DIMENSION(:),   POINTER :: ENERGY_STORE
   REAL ( KIND = DP ), DIMENSION(:,:), POINTER :: VARIABLES_STORE

   ! . Check that there are variables.
   IF ( NVARIABLES <= 0 ) RETURN

   ! . Check that there is a grid.
   IF ( .NOT. QGRID ) RETURN

   ! . There are structures on an old file.
   IF ( PRESENT ( OLD_STRUCTURES ) ) THEN
      CALL STRUCTURES_READ ( OLD_STRUCTURES, NOPTIMIZATIONS, NVARIABLES, NSTRUCTURES, ENERGY_STORE, VARIABLES_STORE )

   ! . There are no old structures.
   ELSE
      NSTRUCTURES = 0
      ALLOCATE ( ENERGY_STORE(1:NOPTIMIZATIONS), VARIABLES_STORE(1:NVARIABLES,1:NOPTIMIZATIONS) )
   END IF

   ! . Save NSTRUCTURES.
   NSTROLD = NSTRUCTURES

   ! . Do some printing.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "FF0000" )
   CALL PRINT_SUMMARY_START ( "Quick-Dock Optimization" )
   WRITE ( PRINT_LINE, "(I14)"   ) NOPTIMIZATIONS ; CALL PRINT_SUMMARY_ELEMENT ( "Optimizations"       )
   WRITE ( PRINT_LINE, "(I14)"   ) NSEARCH        ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Searches"  )
   WRITE ( PRINT_LINE, "(I14)"   ) NHUNT          ; CALL PRINT_SUMMARY_ELEMENT ( "Hunt Configurations" )
   WRITE ( PRINT_LINE, "(I14)"   ) NFINE          ; CALL PRINT_SUMMARY_ELEMENT ( "Fine Configurations" )
   WRITE ( PRINT_LINE, "(I14)"   ) NSTRUCTURES    ; CALL PRINT_SUMMARY_ELEMENT ( "Previous Structures" )
   WRITE ( PRINT_LINE, "(I14)"   ) NVARIABLES     ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Variables" )
   CALL PRINT_SUMMARY_STOP

   ! . Initialize the excluded search range.
   VARIABLES_RANGE_EXCLUDED = VARIABLES_RANGE_INITIALIZE ( .TRUE. )

   ! . Get the energy of the starting structure as a reference.
   VARIABLES = 0.0_DP
   CALL QUICKDOCK_COORDINATES_LIGAND ( VARIABLES )
   EREFERENCE = QUICKDOCK_ENERGY ( )

   ! . Write out some data.
   WRITE ( PRINT_LINE, "(A,F20.5,A)" ) "Energy of reference structure = ", EREFERENCE, " kJ mol^-1."
   CALL PRINT_PARAGRAPH_START
   CALL PRINT_TEXT ( PRINT_LINE )

   ! . Initialize the debug array.
   VARIABLES_DEBUG = 0.0_DP

   ! . Loop over the optimizations.
   DO IOPT = 1,NOPTIMIZATIONS

      ! . Initialize the starting structure as the reference structure.
      EBEST = EREFERENCE ; VARIABLES_BEST = 0.0_DP

      ! . Loop over the searches.
      DO ISEARCH = 1,NSEARCH

         ! . Set the search range.
         VARIABLES_RANGE = VARIABLES_RANGE_INITIALIZE ( .FALSE. )

         ! . Initialization.
	 NDONE = 0

	 ! . Loop over the configurations for the search.
	 DO IHUNT = 1,NHUNT

            ! . Generate a new set of variables.
	    CALL VARIABLES_NEW

            ! . Perform a recombination.
	    CALL VARIABLES_RECOMBINE ( NDONE, NHUNT )

            ! . Reject the configuration if all variables are in the excluded zone.
	    IF ( VARIABLES_COMPARE ( VARIABLES, VARIABLES_BEST, VARIABLES_RANGE_EXCLUDED ) ) CYCLE

            ! . Reject the configuration if it is in the exclusion range of any old structure.
	    IF ( QEXCLUDE ) THEN
	       QSKIP = .FALSE.
               DO I = 1,NSTROLD
	          IF ( VARIABLES_COMPARE ( VARIABLES, VARIABLES_STORE(:,I), VARIABLES_RANGE_EXCLUDED ) ) THEN
                     QSKIP = .TRUE. ; EXIT
                  END IF
	       END DO
	       IF ( QSKIP ) CYCLE
	    END IF

            ! . The configuration has been accepted so increment NDONE.
	    NDONE = NDONE + 1

            ! . Generate the coordinates.
	    CALL QUICKDOCK_COORDINATES_LIGAND ( VARIABLES )

            ! . Calculate the energy.
	    E = QUICKDOCK_ENERGY ( )

            ! . Save the configuration if it is of lower energy.
	    IF ( E < EBEST ) THEN
	       EBEST = E ; VARIABLES_BEST = VARIABLES
	    END IF

            ! . Reduce the search range.
	    CALL VARIABLES_RANGE_REDUCE ( NDONE, NHUNT )

            ! . Do some debug printing.
	    IF ( QDEBUGQ ) THEN
	       WRITE ( OUTPUT, "(A,4I10,G20.8)" ) "Hunt Optimization = ", IOPT, ISEARCH, IHUNT, NDONE, E
	    END IF

            ! . Do not continue if the variable range is smaller than the exclusion zone.
	    IF ( ALL ( ( VARIABLES_RANGE - VARIABLES_RANGE_EXCLUDED ) < 0.0_DP ) ) EXIT

	 END DO

      END DO

      ! . Set the search range to the excluded zone.
      VARIABLES_RANGE = VARIABLES_RANGE_EXCLUDED

      ! . Loop over the configurations for the final refinement.
      DO IFINE = 1,NFINE

         ! . Generate a new set of variables.
	 CALL VARIABLES_NEW

         ! . Generate the coordinates.
	 CALL QUICKDOCK_COORDINATES_LIGAND ( VARIABLES )

         ! . Calculate the energy.
	 E = QUICKDOCK_ENERGY ( )

         ! . Save the configuration if it is of lower energy.
	 IF ( E < EBEST ) THEN
	    EBEST = E ; VARIABLES_BEST = VARIABLES
	 END IF

         ! . Reduce the search range.
	 CALL VARIABLES_RANGE_REDUCE ( IFINE, NFINE )

         ! . Do some debug printing.
	 IF ( QDEBUGQ ) THEN
	    WRITE ( OUTPUT, "(A,2I10,G20.8)" ) "Fine Optimization = ", IOPT, IFINE, E
	 END IF

         ! . Do not continue if the variable range is too small.
	 IF ( ALL ( ABS ( VARIABLES_RANGE ) < VARTOL ) ) EXIT

      END DO

      ! . Increment the number of structures.
      NSTRUCTURES = NSTRUCTURES + 1

      ! . Save the best energy and variables.
      ENERGY_STORE   (  NSTRUCTURES) = EBEST
      VARIABLES_STORE(:,NSTRUCTURES) = VARIABLES_BEST

   END DO

   ! . Write out some data.
   WRITE ( PRINT_LINE, "(A,F20.5,A)" ) "Energy of lowest    structure = ", MINVAL ( ENERGY_STORE ), " kJ mol^-1."
   CALL PRINT_TEXT ( PRINT_LINE )
   CALL PRINT_PARAGRAPH_STOP

   ! . Order the structures.
   CALL STRUCTURES_ORDER ( NSTRUCTURES, NVARIABLES, ENERGY_STORE, VARIABLES_STORE )

   ! . Remove redundant structures.
   IF ( QREDUCE ) CALL STRUCTURES_REDUCE ( NSTRUCTURES, NVARIABLES, VARIABLES_RANGE_EXCLUDED, ENERGY_STORE, VARIABLES_STORE )

   ! . Write out the structures to the new file.
   CALL STRUCTURES_WRITE ( NEW_STRUCTURES, ENERGY_STORE(1:NSTRUCTURES), VARIABLES_STORE(:,1:NSTRUCTURES) )

   ! . Deallocate space.
   DEALLOCATE ( ENERGY_STORE, VARIABLES_STORE )

   ! . Write out some debug information.
   IF ( QDEBUGQ ) THEN
      WRITE ( OUTPUT, "(/A)" ) "Maximum Variable Displacements:"
      WRITE ( OUTPUT, "(6G20.10)" ) VARIABLES_DEBUG
   END IF

   !============================================================================
   CONTAINS
   !============================================================================

      !-----------------------
      SUBROUTINE VARIABLES_NEW
      !-----------------------

      ! . Local scalars.
      INTEGER :: I, N

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:NVARIABLES) :: NUMBERS

      ! . Get a vector of random numbers in the range -1/2,+1/2.
      NUMBERS = RANDOM_VECTOR ( NVARIABLES ) - 0.5_DP

      ! . Generate the new set of variables.
      VARIABLES = VARIABLES_BEST + NUMBERS * VARIABLES_RANGE

      ! . Make sure the new variables are within range.
      ! . Rotations (-PI/+PI) and translations (-TR/+TR).
      ! . Initialization.
      N = 0

      ! . Rotatable angles.
      IF ( NROTATABLE > 0 ) THEN
         VARIABLES(N+1:N+NROTATABLE) = VARIABLES(N+1:N+NROTATABLE) - ( 2.0_DP * PI ) * &
	                                                    ANINT ( VARIABLES(N+1:N+NROTATABLE) / ( 2.0_DP * PI ), DP )
         N = N + NROTATABLE
      END IF

      ! . Translations.
      IF ( QTRANSLATION ) THEN
         DO I = 1,3
	    IF ( VARIABLES(I+N) > TRANSLATION_RANGE ) THEN
	       VARIABLES(I+N) =  TRANSLATION_RANGE
	    ELSE IF ( VARIABLES(I+N) < -TRANSLATION_RANGE ) THEN
	       VARIABLES(I+N) = -TRANSLATION_RANGE
	    END IF
	 END DO
	 N = N + 3
      END IF

      ! . Rotations.
      IF ( QROTATION ) THEN
         VARIABLES(N+1:N+3) = VARIABLES(N+1:N+3) - ( 2.0_DP * PI ) * ANINT ( VARIABLES(N+1:N+3) / ( 2.0_DP * PI ), DP )
         N = N + 3
      END IF

      ! . Check the variable range when debugging.
      IF ( QDEBUGQ ) VARIABLES_DEBUG = MAX ( VARIABLES_DEBUG, ABS ( NUMBERS * VARIABLES_RANGE ) )

      END SUBROUTINE VARIABLES_NEW

      !--------------------------------------------------
      SUBROUTINE VARIABLES_RANGE_REDUCE ( IDONE, ITOTAL )
      !--------------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(IN) :: IDONE, ITOTAL

      ! . Local scalars.
      REAL ( KIND = DP ) :: FREDUCE

      ! . Return if there is no reduction to be done.
      IF ( REDUCEFACTOR == 0.0_DP ) RETURN

      ! . Calculate the factor by which to reduce the search range.
      FREDUCE = EXP ( - REDUCEFACTOR * REAL ( IDONE, DP ) / REAL ( ITOTAL, DP ) )

      ! . Reduce the search range.
      VARIABLES_RANGE = FREDUCE * VARIABLES_RANGE

      END SUBROUTINE VARIABLES_RANGE_REDUCE

      !-----------------------------------------------
      SUBROUTINE VARIABLES_RECOMBINE ( IDONE, ITOTAL )
      !-----------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(IN) :: IDONE, ITOTAL

      ! . Local scalars.
      INTEGER            :: ISTRUCTURE
      REAL ( KIND = DP ) :: PRECOMBINATION

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:NVARIABLES) :: PROBABILITIES

      ! . Return if there are no previous structures.
      IF ( NSTRUCTURES <= 0 ) RETURN

      ! . There is no recombination if IDONE is 0.
      IF ( IDONE <= 0 ) RETURN

      ! . Return if there is no recombination to be done.
      IF ( ( RECOMBINATIONFACTOR == 0.0_DP ) .OR. ( REDUCEFACTOR == 0.0_DP ) ) RETURN

      ! . Choose a previous structure at random.
      IF ( NSTRUCTURES == 1 ) THEN
         ISTRUCTURE = 1
      ELSE
         ISTRUCTURE = INT ( RANDOM ( ) * REAL ( NSTRUCTURES, DP ) ) + 1
      END IF

      ! . Calculate the recombination probability.
      PRECOMBINATION = 1.0_DP - EXP ( - RECOMBINATIONFACTOR * REDUCEFACTOR * REAL ( IDONE, DP ) / REAL ( ITOTAL, DP ) )

      ! . Get a vector of probabilities.
      PROBABILITIES = RANDOM_VECTOR ( NVARIABLES )

      ! . Swap the selected variables.
      WHERE ( PROBABILITIES < PRECOMBINATION ) VARIABLES = VARIABLES_STORE(:,ISTRUCTURE)

      END SUBROUTINE VARIABLES_RECOMBINE

   END SUBROUTINE QUICKDOCK_OPTIMIZE

   !-------------------------------------------------------------------------------------------------------------
   SUBROUTINE QUICKDOCK_OPTIONS ( DIELECTRIC, LIGAND_RECEPTOR_CUTOFF_ENERGY, OFFGRID_ENERGY,                    &
                                  DIHEDRAL_EXCLUSION_ZONE, ROTATION_EXCLUSION_ZONE, TRANSLATION_EXCLUSION_ZONE, &
				  INTERACTION_CUTOFF, RECOMBINATION_FACTOR, REDUCTION_FACTOR )
   !-------------------------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: DIELECTRIC, LIGAND_RECEPTOR_CUTOFF_ENERGY, OFFGRID_ENERGY,                    &
                                               DIHEDRAL_EXCLUSION_ZONE, ROTATION_EXCLUSION_ZONE, TRANSLATION_EXCLUSION_ZONE, &
				               INTERACTION_CUTOFF, RECOMBINATION_FACTOR, REDUCTION_FACTOR

   ! . Initialize the options.
   CALL QUICKDOCK_INITIALIZE_OPTIONS

   ! . Check the input arguments.
   IF ( PRESENT ( DIELECTRIC                    ) ) EPSILON              = DIELECTRIC
   IF ( PRESENT ( LIGAND_RECEPTOR_CUTOFF_ENERGY ) ) ERECMAX              = LIGAND_RECEPTOR_CUTOFF_ENERGY
   IF ( PRESENT ( OFFGRID_ENERGY                ) ) EOFFGRID             = OFFGRID_ENERGY
   IF ( PRESENT ( DIHEDRAL_EXCLUSION_ZONE       ) ) EXCLUDED_DIHEDRAL    = DIHEDRAL_EXCLUSION_ZONE * TO_RADIANS
   IF ( PRESENT ( ROTATION_EXCLUSION_ZONE       ) ) EXCLUDED_ROTATION    = ROTATION_EXCLUSION_ZONE * TO_RADIANS
   IF ( PRESENT ( TRANSLATION_EXCLUSION_ZONE    ) ) EXCLUDED_TRANSLATION = TRANSLATION_EXCLUSION_ZONE
   IF ( PRESENT ( RECOMBINATION_FACTOR          ) ) RECOMBINATIONFACTOR  = RECOMBINATION_FACTOR
   IF ( PRESENT ( REDUCTION_FACTOR              ) ) REDUCEFACTOR         = REDUCTION_FACTOR
   IF ( PRESENT ( INTERACTION_CUTOFF            ) ) MINDIST2             = INTERACTION_CUTOFF**2

   ! . Do some printing.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "FF0000" )
   CALL PRINT_SUMMARY_START ( "Quick-Dock Options" )
   WRITE ( PRINT_LINE, "(G14.8)" ) EOFFGRID             ; CALL PRINT_SUMMARY_ELEMENT ( "Off-Grid Energy"       )
   WRITE ( PRINT_LINE, "(G14.8)" ) EPSILON              ; CALL PRINT_SUMMARY_ELEMENT ( "Dielectric" 	       )
   WRITE ( PRINT_LINE, "(G14.8)" ) ERECMAX              ; CALL PRINT_SUMMARY_ELEMENT ( "Max. LR Energy"        )
   WRITE ( PRINT_LINE, "(G14.8)" ) EXCLUDED_DIHEDRAL    ; CALL PRINT_SUMMARY_ELEMENT ( "Dihedral Exclusion"    )
   WRITE ( PRINT_LINE, "(G14.8)" ) EXCLUDED_ROTATION    ; CALL PRINT_SUMMARY_ELEMENT ( "Rotation Exclusion"    )
   WRITE ( PRINT_LINE, "(G14.8)" ) EXCLUDED_TRANSLATION ; CALL PRINT_SUMMARY_ELEMENT ( "Translation Exclusion" )
   WRITE ( PRINT_LINE, "(G14.8)" ) RECOMBINATIONFACTOR  ; CALL PRINT_SUMMARY_ELEMENT ( "Recombination Factor"  )
   WRITE ( PRINT_LINE, "(G14.8)" ) REDUCEFACTOR         ; CALL PRINT_SUMMARY_ELEMENT ( "Reduction Factor"      )
   WRITE ( PRINT_LINE, "(G14.8)" ) SQRT ( MINDIST2 )    ; CALL PRINT_SUMMARY_ELEMENT ( "Interaction Cutoff"    )
   CALL PRINT_SUMMARY_STOP

   END SUBROUTINE QUICKDOCK_OPTIONS

   !---------------------------------------
   SUBROUTINE QUICKDOCK_RANDOMIZE_REFERENCE
   !---------------------------------------

   ! . Local scalars.
   INTEGER            :: I
   REAL ( KIND = DP ) :: E, EREFERENCE

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES) :: VARIABLES, VARIABLES_RANGE

   ! . Check that there are variables.
   IF ( NVARIABLES <= 0 ) RETURN

   ! . Check that there is a grid.
   IF ( .NOT. QGRID ) RETURN

   ! . Get the energy of the starting structure as a reference.
   VARIABLES = 0.0_DP
   CALL QUICKDOCK_COORDINATES_LIGAND ( VARIABLES )
   EREFERENCE = QUICKDOCK_ENERGY ( )

   ! . Set the range within which the variables can vary.
   VARIABLES_RANGE = VARIABLES_RANGE_INITIALIZE ( .FALSE. )

   ! . Generate the new set of variables.
   VARIABLES = VARIABLES_RANGE * ( RANDOM_VECTOR ( NVARIABLES ) - 0.5_DP )

   ! . Generate the coordinates and its new energy.
   CALL QUICKDOCK_COORDINATES_LIGAND ( VARIABLES )
   E = QUICKDOCK_ENERGY ( )

   ! . Write out some data.
   CALL PRINT_PARAGRAPH_START
   WRITE ( PRINT_LINE, "(A,F20.5,A)" ) "Energy of starting   structure = ", EREFERENCE, " kJ mol^-1."
   CALL PRINT_TEXT ( PRINT_LINE )
   WRITE ( PRINT_LINE, "(A,F20.5,A)" ) "Energy of randomized structure = ", E,          " kJ mol^-1."
   CALL PRINT_TEXT ( PRINT_LINE )
   CALL PRINT_PARAGRAPH_STOP

   ! . Put the coordinates into ATMCRD.
   ATMCRD = REFCRD
   DO I = 1,NDATOMS
      ATMCRD(1:3,DATMIND(I)) = DATMCRD(1:3,I)
   END DO

   ! . Write out a summary.
   CALL PRINT_PARAGRAPH ( "Reference structure randomized." )

   END SUBROUTINE QUICKDOCK_RANDOMIZE_REFERENCE

   !----------------------------------------------------------------------------------------------------------------------
   SUBROUTINE QUICKDOCK_SETUP ( LIGAND_ATOMS_BY_GROUPS, ROTATABLE_BONDS, ROTATABLE_FLAGS, ROTATABLE_TRANSLATABLE_GROUPS, &
                                REFERENCE_ATOM, MAXIMUM_TRANSLATION )
   !----------------------------------------------------------------------------------------------------------------------

   ! . Array arguments.
   INTEGER, DIMENSION(1:NATOMS), INTENT(IN) :: LIGAND_ATOMS_BY_GROUPS

   ! . Optional arguments.
   INTEGER,                 INTENT(IN), OPTIONAL :: REFERENCE_ATOM
   REAL ( KIND = DP ),      INTENT(IN), OPTIONAL :: MAXIMUM_TRANSLATION
   INTEGER, DIMENSION(:),   INTENT(IN), OPTIONAL :: ROTATABLE_TRANSLATABLE_GROUPS
   INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: ROTATABLE_BONDS
   LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: ROTATABLE_FLAGS

   ! . Local scalars.
   INTEGER :: IATOM, IROT, JATOM, NLATOMS

   ! . Local allocatable arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DATMGRP, RATMIND

   ! . Local fixed arrays.
   LOGICAL, DIMENSION(1:NATOMS) :: LIGAND_ATOMS

   ! . Initialize the quickdock grid and system data structures.
   CALL QUICKDOCK_INITIALIZE_GRID
   CALL QUICKDOCK_INITIALIZE_SYSTEM

   ! . Check to see if there are atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Fill LIGAND_ATOMS.
   LIGAND_ATOMS = ( LIGAND_ATOMS_BY_GROUPS > 0 )

   ! . Check the number of ligand atoms.
   NLATOMS = COUNT ( LIGAND_ATOMS )
   IF ( ( NLATOMS <= 0 ) .OR. ( NLATOMS >= NATOMS ) ) THEN
      CALL PRINT_ERROR ( "QUICKDOCK_SETUP", "There are too few or too many ligand atoms." )
   END IF

   ! . Set NRATOMS.
   NRATOMS = NATOMS

   ! . Save the reference coordinates.
   ALLOCATE ( REFCRD(1:3,1:NRATOMS) ) ; REFCRD = ATMCRD

   ! . Check for rotatable bonds.
   IF ( PRESENT ( ROTATABLE_BONDS ) ) THEN

      ! . Get the number of bonds.
      NROTATABLE = SIZE ( ROTATABLE_BONDS, 2 )

      ! . Save the bonds.
      ALLOCATE ( ROTATABLE(1:2,1:NROTATABLE) ) ; ROTATABLE = ROTATABLE_BONDS

      ! . ROTATABLE_FLAGS is present.
      IF ( PRESENT ( ROTATABLE_FLAGS ) ) THEN

         ! . Check the dimensions.
	 IF ( ( SIZE ( ROTATABLE_FLAGS, 1 ) /= NLATOMS    ) .OR. &
	      ( SIZE ( ROTATABLE_FLAGS, 2 ) /= NROTATABLE ) ) THEN
            CALL PRINT_ERROR ( "QUICKDOCK_SETUP", "ROTATABLE_FLAGS is of the wrong dimension." )
	 END IF

      ! . ROTATABLE_FLAGS is absent.
      ELSE
         CALL PRINT_ERROR ( "QUICKDOCK_SETUP", "ROTATABLE_FLAGS is missing." )
      END IF

   ! . There are no rotatable bonds.
   ELSE
      IF ( PRESENT ( ROTATABLE_FLAGS ) ) THEN
         CALL PRINT_ERROR ( "QUICKDOCK_SETUP", "ROTATABLE_FLAGS is present without ROTATABLE_BONDS." )
      END IF
   END IF

   ! . Set up the dockable atoms and bonds.
   CALL SETUP_DOCKABLE_ATOMS
   CALL SETUP_DOCKABLE_BONDS

   ! . Check the bond indices.
   IF ( NROTATABLE > 0 ) THEN

      ! . Loop over the bonds.
      DO IROT = 1,NROTATABLE

         ! . Get the bond indices.
	 IATOM = RATMIND(ROTATABLE(1,IROT))
	 JATOM = RATMIND(ROTATABLE(2,IROT))

         ! . Check to see if they are out of range.
	 IF ( ( IATOM <= 0 ) .OR. ( JATOM <= 0 ) ) THEN
	    CALL PRINT_ERROR ( "QUICKDOCK_SETUP", "There is a bond involving a non-dockable atom." )
	 END IF

         ! . Save IATOM and JATOM.
	 ROTATABLE(1,IROT) = IATOM
	 ROTATABLE(2,IROT) = JATOM

      END DO

      ! . Allocate and set QDIHEDRAL.
      ALLOCATE ( QDIHEDRAL(1:NDATOMS,1:NROTATABLE) ) ; QDIHEDRAL = .FALSE. ; QDIHEDRAL(1:NLATOMS,1:NROTATABLE) = ROTATABLE_FLAGS

   END IF

   ! . Set up the rotatable and translatable atoms.
   CALL SETUP_ROTTRANS_ATOMS

   ! . Set the number of variables.
   NVARIABLES = NROTATABLE
   IF ( QTRANSLATION ) NVARIABLES = NVARIABLES + 3
   IF ( QROTATION    ) NVARIABLES = NVARIABLES + 3

   ! . Check NVARIABLES.
   IF ( NVARIABLES <= 0 ) CALL PRINT_ERROR ( "QUICKDOCK_SETUP", "No varaibles have been specified." )

   ! . Set up the ligand arrays.
   CALL SETUP_LIGAND_ARRAYS

   ! . Deallocate space.
   DEALLOCATE ( DATMGRP, RATMIND )

   ! . Do some printing.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "FF0000" )
   CALL PRINT_SUMMARY_START ( "Quick-Dock Setup" )
   WRITE ( PRINT_LINE, "(I14)"   ) NDATOMS           ; CALL PRINT_SUMMARY_ELEMENT ( "Dockable Atoms"    )
   WRITE ( PRINT_LINE, "(I14)"   ) NRATOMS           ; CALL PRINT_SUMMARY_ELEMENT ( "Reference Atoms"   )
   WRITE ( PRINT_LINE, "(I14)"   ) NROTATABLE        ; CALL PRINT_SUMMARY_ELEMENT ( "Rotatable Bonds"   )
   WRITE ( PRINT_LINE, "(I14)"   ) NVARIABLES        ; CALL PRINT_SUMMARY_ELEMENT ( "Variables"         )
   WRITE ( PRINT_LINE, "(I14)"   ) NLIGNB            ; CALL PRINT_SUMMARY_ELEMENT ( "NB Interactions"   )
   WRITE ( PRINT_LINE, "(I14)"   ) NLIGNB14          ; CALL PRINT_SUMMARY_ELEMENT ( "1-4 Interactions"  )
   WRITE ( PRINT_LINE, "(I14)"   ) NDBONDS           ; CALL PRINT_SUMMARY_ELEMENT ( "Ref. Bonds"        )
   WRITE ( PRINT_LINE, "(I14)"   ) NxDIHEDRALS        ; CALL PRINT_SUMMARY_ELEMENT ( "Energy Dihedrals"  )
   WRITE ( PRINT_LINE, "(I14)"   ) REFATM            ; CALL PRINT_SUMMARY_ELEMENT ( "Reference Atom"    )
   WRITE ( PRINT_LINE, "(G14.8)" ) TRANSLATION_RANGE ; CALL PRINT_SUMMARY_ELEMENT ( "Translation Range" )
   IF ( QDEBUGQ ) THEN
      WRITE ( PRINT_LINE, "(I14)" ) NDANGLES	; CALL PRINT_SUMMARY_ELEMENT ( "Ref. Angles"	)
      WRITE ( PRINT_LINE, "(I14)" ) NDDIHEDRALS ; CALL PRINT_SUMMARY_ELEMENT ( "Ref. Dihedrals" )
   END IF
   CALL PRINT_SUMMARY_STOP

   !============================================================================
   CONTAINS
   !============================================================================

      !------------------------------
      SUBROUTINE SETUP_DOCKABLE_ATOMS
      !------------------------------

      ! . Local scalars.
      INTEGER :: I, IBOND, IROT, J, JBOND, K, LIGATM, NGROUP, NROT, RECATM
      LOGICAL :: QLIGREC, QROT

      ! . Local arrays.
      INTEGER, DIMENSION(1:NATOMS) :: TMPGRP, TMPIND

      ! . Initialize NDATOMS.
      NDATOMS = 0

      ! . Initialize TMPGRP and TMPIND.
      TMPIND = 0
      DO I = 1,NRATOMS
         IF ( LIGAND_ATOMS(I) ) THEN
	    NDATOMS = NDATOMS + 1
	    TMPGRP(NDATOMS) = LIGAND_ATOMS_BY_GROUPS(I)
	    TMPIND(NDATOMS) = I
	 END IF
      END DO

      ! . Get the maximum group index.
      NGROUP = MAXVAL ( LIGAND_ATOMS_BY_GROUPS )

      ! . Loop over the bonds defined for the system.
      DO IBOND = 1,NBONDS

         ! . Get the atom indices.
	 I = BONDS(IBOND)%I
	 J = BONDS(IBOND)%J

         ! . Check to see if this bond is between ligand and receptor.
	 QLIGREC = ( LIGAND_ATOMS(I) .AND. ( .NOT. LIGAND_ATOMS(J) ) ) .OR. &
	           ( LIGAND_ATOMS(J) .AND. ( .NOT. LIGAND_ATOMS(I) ) )

         ! . There is a ligand-receptor bond.
	 IF ( QLIGREC ) THEN

            ! . Initialize NROT and QROT.
	    NROT = 0 ; QROT = .FALSE.

            ! . Verify that this bond is in the rotatable set.
            IF ( NROTATABLE > 0 ) THEN
               DO IROT = 1,NROTATABLE
                  IF ( ( ( ROTATABLE(1,IROT) == I ) .AND. ( ROTATABLE(2,IROT) == J ) ) .OR. &
		       ( ( ROTATABLE(2,IROT) == I ) .AND. ( ROTATABLE(1,IROT) == J ) ) ) THEN
		     NROT = IROT ; QROT = .TRUE. ; EXIT
		  END IF
               END DO
	    END IF

            ! . The bond wasn't found.
	    IF ( .NOT. QROT ) CALL PRINT_ERROR ( "SETUP_DOCKABLE_ATOMS", "A ligand/receptor bond was not specified as rotatable." )

            ! . Find the ligand and receptor atoms.
	    IF ( LIGAND_ATOMS(I) ) THEN
	       LIGATM = I ; RECATM = J 
	    ELSE
	       LIGATM = J ; RECATM = I
	    END IF

            ! . Make the TMPIND entry for LIGATM negative.
	    DO K = 1,NLATOMS
	       IF ( ABS ( TMPIND(K) ) == LIGATM ) THEN
	          TMPIND(K) = - LIGATM
	       END IF
	    END DO

            ! . Increment NGROUP.
            NGROUP = NGROUP + 1

            ! . Include RECATM in the dockable set.
	    NDATOMS = NDATOMS + 1
	    TMPIND(NDATOMS) = - RECATM
	    TMPGRP(NDATOMS) = NGROUP

            ! . Find all non-ligand atoms bound to RECATM and store them.
	    DO JBOND = 1,NBONDS
               IF ( BONDS(JBOND)%I == RECATM ) THEN
	          IF ( .NOT. LIGAND_ATOMS(BONDS(JBOND)%J) ) THEN
		     IF ( .NOT. ANY ( TMPIND(NLATOMS+1:NDATOMS) == - BONDS(JBOND)%J ) ) THEN
  		        NDATOMS = NDATOMS + 1
		        TMPIND(NDATOMS) = - BONDS(JBOND)%J
	                TMPGRP(NDATOMS) = NGROUP
	             END IF
		  END IF
	       ELSE IF ( BONDS(JBOND)%J == RECATM ) THEN
	          IF ( .NOT. LIGAND_ATOMS(BONDS(JBOND)%I) ) THEN
		     IF ( .NOT. ANY ( TMPIND(NLATOMS+1:NDATOMS) == - BONDS(JBOND)%I ) ) THEN
  		        NDATOMS = NDATOMS + 1
		        TMPIND(NDATOMS) = - BONDS(JBOND)%I
	                TMPGRP(NDATOMS) = NGROUP
	             END IF
		  END IF
	       END IF
	    END DO
	 END IF
      END DO

      ! . Allocate some arrays.
      ALLOCATE ( DATMGRP(1:NDATOMS), DATMIND(1:NDATOMS), RATMIND(1:NRATOMS), QDFIXED(1:NDATOMS), QRFIXED(1:NRATOMS) )

      ! . Save DATMIND.
      DATMIND = ABS ( TMPIND(1:NDATOMS) )

      ! . Save DATMGRP.
      DATMGRP = TMPGRP(1:NDATOMS)

      ! . Save QDFIXED.
      QDFIXED = ( TMPIND(1:NDATOMS) < 0 )

      ! . Initialize QRFIXED and RATMIND.
      QRFIXED = .TRUE. ; RATMIND = 0

      ! . Fill QRFIXED and RATMIND.
      DO I = 1,NDATOMS
         IATOM = DATMIND(I)
	 QRFIXED(IATOM) = .FALSE.
	 RATMIND(IATOM) = I
      END DO

      ! . Do some debug printing.
      IF ( QDEBUGQ ) THEN
         WRITE ( OUTPUT, "(/A,I6)" ) "Number of dockable atoms = ", NDATOMS
	 WRITE ( OUTPUT, "(4(2I6,L6,6X))" ) ( DATMIND(I), DATMGRP(I), QDFIXED(I), I = 1,NDATOMS )
      END IF

      END SUBROUTINE SETUP_DOCKABLE_ATOMS

      !------------------------------
      SUBROUTINE SETUP_DOCKABLE_BONDS
      !------------------------------

      ! . Local scalars.
      INTEGER :: IATOM, IBOND, JATOM, N

      ! . Local arrays.
      INTEGER, DIMENSION(1:NBONDS) :: TMPIND

      ! . Initialization.
      N = 0

      ! . Loop over the bonds.
      DO IBOND = 1,NBONDS

         ! . Get the dockable atom indices.
	 IATOM = RATMIND(BONDS(IBOND)%I)
	 JATOM = RATMIND(BONDS(IBOND)%J)

         ! . Save the bond index if both atoms are dockable.
	 IF ( ( IATOM > 0 ) .AND. ( JATOM > 0 ) ) THEN
	    N = N + 1
	    TMPIND(N) = IBOND
	 END IF

      END DO

      ! . Save the number of bonds.
      NDBONDS = N

      ! . Allocate DBONDS.
      ALLOCATE ( DBONDS(1:2,1:NDBONDS) )

      ! . Save the bonds.
      DO IBOND = 1,NDBONDS
         DBONDS(1,IBOND) = RATMIND(BONDS(TMPIND(IBOND))%I)
	 DBONDS(2,IBOND) = RATMIND(BONDS(TMPIND(IBOND))%J)
      END DO

      ! . Make sure the bonds are in the correct order.
      CALL ORDER_BONDS ( DBONDS )

      END SUBROUTINE SETUP_DOCKABLE_BONDS

      !-----------------------------
      SUBROUTINE SETUP_LIGAND_ARRAYS
      !-----------------------------

      ! . Local scalars.
      INTEGER :: I, IATOM, IDIH, IINT, IROT, J, JATOM, N
      LOGICAL :: QEXCL, QMATCH

      ! . Local arrays.
      INTEGER, ALLOCATABLE, DIMENSION(:) :: DIHIND, TMPJ

      ! . Allocate the ligand arrays.
      ALLOCATE ( DATMCHG(1:NDATOMS), DATMCHG14(1:NDATOMS), DATMREP(1:NDATOMS), DATMREP14(1:NDATOMS),      &
                 DATMVDW(1:NDATOMS), DATMVDW14(1:NDATOMS), DATMCRD(1:3,1:NDATOMS), DREFCRD(1:3,1:NDATOMS) )

      ! . Loop over the ligand atoms.
      DO IATOM = 1,NDATOMS

         ! . Get the atom index in the full array.
	 I = DATMIND(IATOM)

         ! . Set the coordinates.
	 DATMCRD(1:3,IATOM) = 0.0_DP
	 DREFCRD(1:3,IATOM) = REFCRD(1:3,I)

         ! . Set up the parameters.
	 DATMCHG  (IATOM) = ATMCHG  (I)
	 DATMCHG14(IATOM) = ATMCHG14(I)
	 DATMREP  (IATOM) = ATMEPS  (I) * ATMSIG(I)**12
	 DATMREP14(IATOM) = ATMEPS14(I) * ATMSIG(I)**12
	 DATMVDW  (IATOM) = ATMEPS  (I) * ATMSIG(I)**6
	 DATMVDW14(IATOM) = ATMEPS14(I) * ATMSIG(I)**6

      END DO

      ! . If there are rotatable bonds get the dihedral energy terms.
      IF ( NROTATABLE > 0 ) THEN

         ! . Allocate space.
	 ALLOCATE ( DIHIND(1:NSYSDIH) ) ; DIHIND = 0

         ! . Initialization.
	 NxDIHEDRALS = 0

         ! . Loop over the system dihedrals.
	 DO I = 1,NSYSDIH

            ! . Initialize for a match.
	    QMATCH = .FALSE.

            ! . Loop over the rotatable bonds for a match (central JK atoms only).
	    DO IROT = 1,NROTATABLE
	       IF ( ( ( SYSDIH(I)%J == DATMIND(ROTATABLE(1,IROT)) ) .AND. ( SYSDIH(I)%K == DATMIND(ROTATABLE(2,IROT)) ) ) .OR. &
	            ( ( SYSDIH(I)%J == DATMIND(ROTATABLE(2,IROT)) ) .AND. ( SYSDIH(I)%K == DATMIND(ROTATABLE(1,IROT)) ) ) ) THEN
	          QMATCH = .TRUE. ; EXIT
	       END IF
	    END DO

            ! . There has been a match so save the dihedral index.
	    IF ( QMATCH ) THEN
	       NxDIHEDRALS = NxDIHEDRALS + 1
	       DIHIND(NxDIHEDRALS) = I
	    END IF

	 END DO

         ! . Allocate space.
	 ALLOCATE ( xDIHEDRALS(1:NxDIHEDRALS) )

         ! . Fill DIHEDRALS.
	 DO I = 1,NxDIHEDRALS
	    xDIHEDRALS(I)   = SYSDIH(DIHIND(I))
	    xDIHEDRALS(I)%I = RATMIND(xDIHEDRALS(I)%I)
	    xDIHEDRALS(I)%J = RATMIND(xDIHEDRALS(I)%J)
	    xDIHEDRALS(I)%K = RATMIND(xDIHEDRALS(I)%K)
	    xDIHEDRALS(I)%L = RATMIND(xDIHEDRALS(I)%L)
	 END DO

         ! . Deallocate space.
	 DEALLOCATE ( DIHIND )

      END IF

      ! . If there are rotatable bonds get the internal non-bonding energy lists.
      IF ( NROTATABLE > 0 ) THEN

         ! . Allocate and initialize the LISTI arrays.
         ALLOCATE ( LISTI(1:NDATOMS+1), LISTI14(1:NDATOMS+1) ) ; LISTI = 0 ; LISTI14 = 0

         ! . 1-4 interactions.
	 IF ( SIZE ( ATME14J ) > 0 ) THEN

            ! . Allocate some temporary lists.
	    ALLOCATE ( TMPJ(1:SIZE(ATME14J)) )

            ! . Initialization.
	    N = 0

            ! . Loop over atoms.
            DO IATOM = 1,NDATOMS

               ! . Set LISTI14 for the atom.
	       LISTI14(IATOM) = N

               ! . Get the atom index.
	       I = DATMIND(IATOM)

               ! . Loop over the interactions for the atom.
	       DO IINT = ATME14I(I)+1,ATME14I(I+1)

                  ! . Get the second atom of the interaction.
		  J = ATME14J(IINT)

                  ! . Get the ligand index of this atom.
		  JATOM = RATMIND(J)

                  ! . Skip if the atom is not a ligand atom.
		  IF ( JATOM <= 0 ) CYCLE

                  ! . Exclude the interaction if both atoms are fixed or both are in the same group.
		  IF ( ( DATMGRP(IATOM) == DATMGRP(JATOM) ) .OR. ( QDFIXED(IATOM) .AND. QDFIXED(JATOM) ) ) CYCLE

                  ! . Include the interaction.
		  N = N + 1
		  TMPJ(N) = JATOM

	       END DO

	    END DO

            ! . Set the last element of LISTI14.
	    LISTI14(NDATOMS+1) = N

            ! . Save the list.
	    NLIGNB14 = N ; ALLOCATE ( LISTJ14(1:N) ) ; LISTJ14 = TMPJ(1:N) ; DEALLOCATE ( TMPJ )

	 END IF

         ! . Allocate temporary space for the normal interactions.
	 ALLOCATE ( TMPJ(1:(NDATOMS*(NDATOMS-1))/2) )

         ! . Initialization.
	 N = 0

         ! . Outer loop over atoms.
	 DO IATOM = 1,NDATOMS

            ! . Set LISTI for the atom.
	    LISTI(IATOM) = N

            ! . Get the atom index.
	    I = DATMIND(IATOM)

            ! . Inner loop over atoms.
	    DO JATOM = (IATOM+1),NDATOMS

               ! . Get the atom index.
	       J = DATMIND(JATOM)

               ! . Exclude the interaction if both atoms are fixed or both are in the same group.
	       IF ( ( DATMGRP(IATOM) == DATMGRP(JATOM) ) .OR. ( QDFIXED(IATOM) .AND. QDFIXED(JATOM) ) ) CYCLE

               ! . Check to see if the atoms are excluded.
	       IF ( J > I ) THEN
	          QEXCL = ANY ( ATMEXCJ(ATMEXCI(I)+1:ATMEXCI(I+1)) == J )
	       ELSE
	          QEXCL = ANY ( ATMEXCJ(ATMEXCI(J)+1:ATMEXCI(J+1)) == I )
	       END IF

               ! . Include non-excluded interactions.
	       IF ( .NOT. QEXCL ) THEN
	          N = N + 1
		  TMPJ(N) = JATOM
	       END IF

	    END DO

	 END DO

         ! . Set the last element of LISTI.
	 LISTI(NDATOMS+1) = N

         ! . Save the lists.
	 IF ( N > 0 ) THEN
	    NLIGNB = N ; ALLOCATE ( LISTJ(1:N) ) ; LISTJ = TMPJ(1:N)
	 END IF

         ! . Deallocate space.
	 DEALLOCATE ( TMPJ )

      END IF

      ! . Set up the bond, angle and dihedral lists if the debug flag is on.
      IF ( QDEBUGQ ) THEN

         ! . Check the number of bonds.
	 IF ( NDBONDS <= 0 ) RETURN

         ! . Allocate and calculate the reference bond distances.
	 ALLOCATE ( REFBND(1:NDBONDS) ) ; REFBND = GEOMETRY_DISTANCES ( DREFCRD, DBONDS )

         ! . Calculate the angle lists from the bond list.
         CALL CONNECTIVITY_ANGLES ( DANGLES, DBONDS )

         ! . Set NDANGLES.
         NDANGLES = SIZE ( DANGLES, 2 )

         ! . Check NDANGLES.
	 IF ( NDANGLES <= 0 ) RETURN

         ! . Allocate and calculate the reference bond angles.
	 ALLOCATE ( REFANG(1:NDANGLES) ) ; REFANG = GEOMETRY_ANGLES ( DREFCRD, DANGLES )

         ! . Calculate the dihderal lists from the bond and angle lists.
         CALL CONNECTIVITY_DIHEDRALS ( DDIHEDRALS, DBONDS, DANGLES )

         ! . Set NDDIHEDRALS.
         NDDIHEDRALS = SIZE ( DDIHEDRALS, 2 )

         ! . Initialization.
	 N = 0

         ! . Loop over the dihedrals to remove rotatable bonds.
	 DO IDIH = 1,NDDIHEDRALS

            ! . Get the central atoms of the dihedral.
	    I = DDIHEDRALS(2,IDIH)
	    J = DDIHEDRALS(3,IDIH)

            ! . Initialize QMATCH.
	    QMATCH = .FALSE.

            ! . Loop over the rotatable bonds for a match.
	    DO IROT = 1,NROTATABLE
	       IF ( ( ( ROTATABLE(1,IROT) == I ) .AND. ( ROTATABLE(2,IROT) == J ) ) .OR. &
	            ( ( ROTATABLE(2,IROT) == I ) .AND. ( ROTATABLE(1,IROT) == J ) ) ) THEN
	          QMATCH = .TRUE. ; EXIT
	       END IF
	    END DO

            ! . Include the dihedral if there is no match.
	    IF ( .NOT. QMATCH ) THEN
	       N = N + 1
	       DDIHEDRALS(:,N) = DDIHEDRALS(:,IDIH)
	    END IF

	 END DO

         ! . Set NDDIHEDRALS.
         NDDIHEDRALS = N

         ! . Check NDDIHEDRALS.
	 IF ( NDDIHEDRALS <= 0 ) RETURN

         ! . Allocate and calculate the reference dihedrals.
	 ALLOCATE ( REFDIH(1:NDDIHEDRALS) ) ; REFDIH = GEOMETRY_DIHEDRALS ( DREFCRD, DDIHEDRALS(:,1:NDDIHEDRALS) )

      END IF

      END SUBROUTINE SETUP_LIGAND_ARRAYS

      !------------------------------
      SUBROUTINE SETUP_ROTTRANS_ATOMS
      !------------------------------

      ! . Local scalars.
      INTEGER :: IATOM, IBOND

      ! . Local arrays.
      LOGICAL, DIMENSION(1:NDATOMS) :: QFLAG

      ! . A set of rotatable and translatable atoms have been specified.
      IF ( PRESENT ( ROTATABLE_TRANSLATABLE_GROUPS ) ) THEN

         ! . Check for the presence of REFERENCE_ATOM and MAXIMUM_TRANSLATION.
	 IF ( .NOT. PRESENT ( REFERENCE_ATOM ) .OR. .NOT. PRESENT ( MAXIMUM_TRANSLATION ) ) THEN
            CALL PRINT_ERROR ( "SETUP_ROTTRANS_ATOMS", "REFERENCE_ATOM or MAXIMUM_TRANSLATION argument missing." )
	 END IF

         ! . Initialize QFLAG.
	 QFLAG = .FALSE.

         ! . Loop over the dockable atoms.
	 DO IATOM = 1,NDATOMS

            ! . Set the flag for this atom.
	    QFLAG(IATOM) = ANY ( DATMGRP(IATOM) == ROTATABLE_TRANSLATABLE_GROUPS )

	 END DO

         ! . Check QFLAG.
	 IF ( COUNT ( QFLAG ) <= 0 ) CALL PRINT_ERROR ( "SETUP_ROTTRANS_ATOMS", "No rotatable or translatable atoms exist." )

         ! . Check to see that there are no bonds between movable and fixed atoms.
	 DO IBOND = 1,NDBONDS
            IF ( ( QFLAG(DBONDS(1,IBOND)) .AND. .NOT. QFLAG(DBONDS(2,IBOND)) ) .OR. &
	         ( QFLAG(DBONDS(2,IBOND)) .AND. .NOT. QFLAG(DBONDS(1,IBOND)) ) ) THEN
	       CALL PRINT_ERROR ( "SETUP_ROTTRANS_ATOMS", "There are bonds between movable and fixed atoms." )
            END IF
	 END DO

         ! . Set the reference atom.
	 REFATM = RATMIND(REFERENCE_ATOM)

         ! . Check to see if the reference atom is in the rotatable/translatable set.
         IF ( REFATM <= 0 ) THEN
	    CALL PRINT_ERROR ( "SETUP_ROTTRANS_ATOMS", "Reference atom is not a dockable atom." )
	 ELSE IF ( .NOT. QFLAG(REFATM) ) THEN
	    CALL PRINT_ERROR ( "SETUP_ROTTRANS_ATOMS", "Reference atom is not among the rotatable or translatable atoms." )
	 END IF

         ! . Set the translation range.
	 TRANSLATION_RANGE = ABS ( MAXIMUM_TRANSLATION )

         ! . Set the rotation and translation flags.
	 QROTATION    = .TRUE.
	 QTRANSLATION = .TRUE.

         ! . Set the QROTATABLE and QTRANSLATABLE arrays.
	 ALLOCATE ( QROTATABLE(1:NDATOMS), QTRANSLATABLE(1:NDATOMS) )
	 QROTATABLE = QFLAG ; QTRANSLATABLE = QFLAG

      ! . No rotatable and translatable atoms have been specified.
      ELSE

         ! . Check for the presence of REFERENCE_ATOM and MAXIMUM_TRANSLATION.
	 IF ( PRESENT ( REFERENCE_ATOM ) .OR. PRESENT ( MAXIMUM_TRANSLATION ) ) THEN
            CALL PRINT_ERROR ( "SETUP_ROTTRANS_ATOMS", "REFERENCE_ATOM or MAXIMUM_TRANSLATION argument "//&
	                                               " present without ROTATABLE_TRANSLATABLE_GROUP." )
	 END IF

      END IF

      END SUBROUTINE SETUP_ROTTRANS_ATOMS

   END SUBROUTINE QUICKDOCK_SETUP

   !--------------------------------------------------
   SUBROUTINE QUICKDOCK_TRAJECTORY ( INFILE, OUTFILE )
   !--------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: INFILE, OUTFILE

   ! . Local scalars.
   INTEGER        :: I, IFRAME, NRFIXED, NSTRUCTURES
   TYPE(DCD_TYPE) :: TRAJECTORY

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(:),   POINTER :: ENERGY_STORE
   REAL ( KIND = DP ), DIMENSION(:,:), POINTER :: VARIABLES_STORE

   ! . Check that there are variables.
   IF ( NVARIABLES <= 0 ) RETURN

   ! . Read the structures.
   CALL STRUCTURES_READ ( INFILE, 0, NVARIABLES, NSTRUCTURES, ENERGY_STORE, VARIABLES_STORE )

   ! . Check NSTRUCTURES.
   IF ( NSTRUCTURES <= 0 ) RETURN

   ! . Get the number of fixed atoms.
   NRFIXED = COUNT ( QRFIXED )

   ! . Activate the trajectory.
   CALL DCD_INITIALIZE ( TRAJECTORY )
   CALL DCD_ACTIVATE_WRITE ( OUTFILE, TRAJECTORY, "CORD", NRATOMS, NRFIXED, ( NSTRUCTURES + 1 ), QFIX = QRFIXED )
   CALL DCD_WRITE ( TRAJECTORY, REFCRD )

   ! . Loop over the structures.
   DO IFRAME = 1,NSTRUCTURES

      ! . Convert the structure to coordinates.
      CALL QUICKDOCK_COORDINATES_LIGAND ( VARIABLES_STORE(:,IFRAME) )

      ! . Put the coordinates into ATMCRD.
      ATMCRD = REFCRD
      DO I = 1,NDATOMS
         ATMCRD(1:3,DATMIND(I)) = DATMCRD(1:3,I)
      END DO

      ! . Write out the structure.
      CALL DCD_WRITE ( TRAJECTORY, ATMCRD )

   END DO

   ! . Deactivate the trajectory.
   CALL DCD_DEACTIVATE ( TRAJECTORY )

   ! . Deallocate space.
   DEALLOCATE ( ENERGY_STORE, VARIABLES_STORE )

   END SUBROUTINE QUICKDOCK_TRAJECTORY

!-------------------------------------------------------------------------------
! . Private coordinate and energy procedures.
!-------------------------------------------------------------------------------

   !-----------------------------------
   SUBROUTINE QUICKDOCK_CHECK_STRUCTURE
   !-----------------------------------

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: ANGTOL = 1.0E-2_DP, BNDTOL = 1.0E-5_DP

   ! . Local scalars.
   REAL ( KIND = DP ) :: MAXDIF

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NDBONDS)     :: STRBND
   REAL ( KIND = DP ), DIMENSION(1:NDANGLES)    :: STRANG
   REAL ( KIND = DP ), DIMENSION(1:NDDIHEDRALS) :: STRDIH

   ! . Check the number of bonds.
   IF ( NDBONDS <= 0 ) RETURN

   ! . Calculate the distances.
   STRBND = GEOMETRY_DISTANCES ( DATMCRD, DBONDS )

   ! . Determine the largest deviation.
   MAXDIF = MAXVAL ( ABS ( REFBND - STRBND ) )

   ! . Write out the largest deviation if it is too large.
   IF ( MAXDIF > BNDTOL ) THEN
      WRITE ( OUTPUT, "(/A,F20.5)" ) "Largest bond distance deviation = ", MAXDIF
   END IF

   ! . Check the number of angles.
   IF ( NDANGLES <= 0 ) RETURN

   ! . Calculate the angles.
   STRANG = GEOMETRY_ANGLES ( DATMCRD, DANGLES )

   ! . Determine the largest deviation.
   MAXDIF = MAXVAL ( ABS ( REFANG - STRANG ) )

   ! . Write out the largest deviation if it is too large.
   IF ( MAXDIF > ANGTOL ) THEN
      WRITE ( OUTPUT, "(/A,F20.5)" ) "Largest bond angle deviation = ", MAXDIF
   END IF

   ! . Check the number of dihedrals.
   IF ( NDDIHEDRALS <= 0 ) RETURN

   ! . Calculate the angles.
   STRDIH = GEOMETRY_DIHEDRALS ( DATMCRD, DDIHEDRALS(:,1:NDDIHEDRALS) )

   ! . Determine the largest deviation.
   MAXDIF = MAXVAL ( ABS ( REFDIH - STRDIH ) )

   ! . Write out the largest deviation if it is too large.
   IF ( MAXDIF > ANGTOL ) THEN
      WRITE ( OUTPUT, "(/A,F20.5)" ) "Largest dihedral angle deviation = ", MAXDIF
   END IF

   END SUBROUTINE QUICKDOCK_CHECK_STRUCTURE

   !----------------------------------------------------
   SUBROUTINE QUICKDOCK_COORDINATES_LIGAND ( VARIABLES )
   !----------------------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: VARIABLES

   ! . Local scalars.
   INTEGER :: I, IROT, N

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3)     :: T
   REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: R

   ! . Copy the reference to the new coordinates.
   DATMCRD = DREFCRD

   ! . Initialization.
   N = 0

   ! . Loop over the rotatable bonds.
   DO IROT = 1,NROTATABLE

      ! . Increment N.
      N = N + 1

      ! . Get the axis for the rotation.
      T = NORMALIZE ( DATMCRD(1:3,ROTATABLE(1,IROT)) - DATMCRD(1:3,ROTATABLE(2,IROT)) )

      ! . Create the rotation matrix (about an axis).
      R = ROTATION_MATRIX_AXIS ( VARIABLES(N), T )

      ! . Rotate the necessary atoms.
      DO I = 1,NDATOMS
         IF ( QDIHEDRAL(I,IROT) ) THEN
	    T = DATMCRD(1:3,I) - DATMCRD(1:3,ROTATABLE(1,IROT))
	    DATMCRD(1:3,I) = MATMUL ( R, T ) + DATMCRD(1:3,ROTATABLE(1,IROT))
	 END IF
      END DO

   END DO

   ! . Translation.
   IF ( QTRANSLATION ) THEN

      ! . Get the translation consisting of a shift of the reference atom to the reference point followed by the variable translation.
      T = DREFCRD(1:3,REFATM) - DATMCRD(1:3,REFATM) + VARIABLES(N+1:N+3)

      ! . Translate the necessary atoms.
      DO I = 1,NDATOMS
         IF ( QTRANSLATABLE(I) ) DATMCRD(1:3,I) = DATMCRD(1:3,I) + T
      END DO

      ! . Increment N.
      N = N + 3

   END IF

   ! . Rotation.
   IF ( QROTATION ) THEN

      ! . Create the rotation matrix (from XYZ angles).
      R = ROTATION_MATRIX_XYZ ( VARIABLES(N+1:N+3) )

      ! . Rotate the necessary atoms.
      DO I = 1,NDATOMS
         IF ( QROTATABLE(I) ) THEN
	    T = DATMCRD(1:3,I) - DATMCRD(1:3,REFATM)
	    DATMCRD(1:3,I) = MATMUL ( R, T ) + DATMCRD(1:3,REFATM)
	 END IF
      END DO

   END IF

   ! . Debug the structure.
   IF ( QDEBUGQ ) CALL QUICKDOCK_CHECK_STRUCTURE

   END SUBROUTINE QUICKDOCK_COORDINATES_LIGAND

   !----------------------------
   FUNCTION QUICKDOCK_ENERGY ( )
   !----------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: QUICKDOCK_ENERGY

   ! . Local scalars.
   INTEGER            :: I, IATOM, J, K
   LOGICAL            :: QOFF
   REAL ( KIND = DP ) :: AM, AP, BM, BP, CM, CP, EDIH, ELIG, ELIG14, EREC, POTQ, POTR, POTV

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: R

   ! . Initialization.
   EDIH   = 0.0_DP
   ELIG   = 0.0_DP
   ELIG14 = 0.0_DP
   EREC   = 0.0_DP
   QOFF   = .FALSE.

   ! . Find the dihedral energy.
   IF ( NxDIHEDRALS > 0 ) EDIH = QUICKDOCK_ENERGY_DIHEDRAL ( )

   ! . Find the internal (or explicitly-calculated) non-bonding energy.
   IF ( NLIGNB   > 0 ) ELIG   = QUICKDOCK_INTERACTIONS ( LISTI,   LISTJ,   DATMCHG,   DATMREP,   DATMVDW   )
   IF ( NLIGNB14 > 0 ) ELIG14 = QUICKDOCK_INTERACTIONS ( LISTI14, LISTJ14, DATMCHG14, DATMREP14, DATMVDW14 )

   ! . Find the external (or grid-calculated) energy.
   ! . Loop over the atoms to be docked.
   DO IATOM = 1,NDATOMS

      ! . Skip fixed atoms.
      IF ( QDFIXED(IATOM) ) CYCLE

      ! . Get the atom coordinates relative to the grid origin.
      R = DATMCRD(1:3,IATOM) - GRIDORIGIN

      ! . The atom lies outside the grid.
      IF ( ANY ( R <= 0.0_DP ) .OR. ANY ( R >= GRIDEXTENT ) ) THEN

         ! . Set the off-grid flag and exit.
         QOFF = .TRUE. ; EXIT

      ! . The atom lies inside the grid.
      ELSE

         ! . Get the indices of the grid point within which the point lies.
         I = INT ( R(1) / GRIDSPACING ) + 1
         J = INT ( R(2) / GRIDSPACING ) + 1
         K = INT ( R(3) / GRIDSPACING ) + 1

         ! . Get the fractional distances within the grid point.
	 AP = ( R(1) - REAL ( I - 1, DP ) * GRIDSPACING ) / GRIDSPACING
	 BP = ( R(2) - REAL ( J - 1, DP ) * GRIDSPACING ) / GRIDSPACING
	 CP = ( R(3) - REAL ( K - 1, DP ) * GRIDSPACING ) / GRIDSPACING

         ! . Get some distance factors.
	 AM = 1.0_DP - AP ; BM = 1.0_DP - BP ; CM = 1.0_DP - CP

         ! . Get the various potentials by trilinear interpolation.
         POTQ = AM * BM * CM * GRIDQ(I,  J,  K  ) + &
	        AP * BM * CM * GRIDQ(I+1,J,  K  ) + &
	        AM * BP * CM * GRIDQ(I,  J+1,K  ) + &
	        AP * BP * CM * GRIDQ(I+1,J+1,K  ) + &
	        AM * BM * CP * GRIDQ(I,  J,  K+1) + &
	        AP * BM * CP * GRIDQ(I+1,J,  K+1) + &
	        AM * BP * CP * GRIDQ(I,  J+1,K+1) + &
	        AP * BP * CP * GRIDQ(I+1,J+1,K+1)
         POTR = AM * BM * CM * GRIDR(I,  J,  K  ) + &
	        AP * BM * CM * GRIDR(I+1,J,  K  ) + &
	        AM * BP * CM * GRIDR(I,  J+1,K  ) + &
	        AP * BP * CM * GRIDR(I+1,J+1,K  ) + &
	        AM * BM * CP * GRIDR(I,  J,  K+1) + &
	        AP * BM * CP * GRIDR(I+1,J,  K+1) + &
	        AM * BP * CP * GRIDR(I,  J+1,K+1) + &
	        AP * BP * CP * GRIDR(I+1,J+1,K+1)
         POTV = AM * BM * CM * GRIDV(I,  J,  K  ) + &
	        AP * BM * CM * GRIDV(I+1,J,  K  ) + &
	        AM * BP * CM * GRIDV(I,  J+1,K  ) + &
	        AP * BP * CM * GRIDV(I+1,J+1,K  ) + &
	        AM * BM * CP * GRIDV(I,  J,  K+1) + &
	        AP * BM * CP * GRIDV(I+1,J,  K+1) + &
	        AM * BP * CP * GRIDV(I,  J+1,K+1) + &
	        AP * BP * CP * GRIDV(I+1,J+1,K+1)

         ! . Sum the energy into the total energy.
         EREC = EREC + DATMCHG(IATOM) * POTQ + DATMREP(IATOM) * POTR - DATMVDW(IATOM) * POTV      

      END IF

   END DO

   ! . Check for an off-grid energy.
   IF ( QOFF ) THEN
      EREC = EOFFGRID
   ! . Check that the ligand-receptor energy is not too large.
   ELSE
      EREC = MIN ( ERECMAX, EREC )
   END IF

   ! . Set the energy.
   QUICKDOCK_ENERGY = EDIH + ELIG + ELIG14 + EREC

   ! . Print out the energies if debugging is on.
   IF ( QDEBUGQ ) WRITE ( OUTPUT, "(A,5G20.10)" ) "Energies = ", QUICKDOCK_ENERGY, EDIH, ELIG, ELIG14, EREC

   END FUNCTION QUICKDOCK_ENERGY

   !-------------------------------------
   FUNCTION QUICKDOCK_ENERGY_DIHEDRAL ( )
   !-------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: QUICKDOCK_ENERGY_DIHEDRAL

   ! . Local scalars.
   INTEGER            :: IDIHE, I, J, K, L
   REAL ( KIND = DP ) :: COSNPHI, COSPHI, COSPHI2, DOTIJ, DOTLK, EDIHEDRAL, R, RKJ, S

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DR, DRIJ, DRKJ, DRLK, DS

   ! . Parameter declarations.
   REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP

   ! . Initialization.
   QUICKDOCK_ENERGY_DIHEDRAL = 0.0_DP

   ! . Return if there are no dihedrals.
   IF ( NxDIHEDRALS <= 0 ) RETURN

   ! . Initialization.
   EDIHEDRAL = 0.0_DP

   ! . Loop over the dihedrals.
   DO IDIHE = 1,NxDIHEDRALS

      ! . Obtain the atoms for the dihedral.
      I = xDIHEDRALS(IDIHE)%I
      J = xDIHEDRALS(IDIHE)%J
      K = xDIHEDRALS(IDIHE)%K
      L = xDIHEDRALS(IDIHE)%L

      ! . Calculate some displacement vectors.
      DRIJ = DATMCRD(1:3,I) - DATMCRD(1:3,J)
      DRKJ = DATMCRD(1:3,K) - DATMCRD(1:3,J)
      DRLK = DATMCRD(1:3,L) - DATMCRD(1:3,K)

      ! . Calculate the size of the RKJ vector.
      RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

      ! . Normalize the vector.
      DRKJ = DRKJ / RKJ

      ! . Calculate some intermediate dot products.
      DOTIJ = DOT_PRODUCT ( DRIJ, DRKJ )
      DOTLK = DOT_PRODUCT ( DRLK, DRKJ )

      ! . Calculate the DR and DS vectors.
      DR = DRIJ - DOTIJ * DRKJ
      DS = DRLK - DOTLK * DRKJ

      ! . Calculate the magnitudes of DR and DS.
      R = SQRT ( DOT_PRODUCT ( DR, DR ) )
      S = SQRT ( DOT_PRODUCT ( DS, DS ) )

      ! . Normalize DR and DS.
      DR = DR / R
      DS = DS / S

      ! . Calculate the dot product of the vectors.
      COSPHI = DOT_PRODUCT ( DR, DS )

      ! . Ensure DOTFAC is between -1 and 1.
      COSPHI = SIGN ( MIN ( ABS ( COSPHI ), DOT_LIMIT ), COSPHI )

      ! . Calculate the square.
      COSPHI2 = COSPHI * COSPHI

      ! . Calculate Cos(n phi) and its derivative for the appropriate periodicity.
      SELECT CASE ( xDIHEDRALS(IDIHE)%PERIOD )
      CASE ( 0 ) ; COSNPHI = 1.0_DP
      CASE ( 1 ) ; COSNPHI = COSPHI
      CASE ( 2 ) ; COSNPHI = 2.0_DP * COSPHI2 - 1.0_DP
      CASE ( 3 ) ; COSNPHI = ( 4.0_DP * COSPHI2 - 3.0_DP ) * COSPHI
      CASE ( 4 ) ; COSNPHI = 8.0_DP * ( COSPHI2 - 1.0_DP ) * COSPHI2 + 1.0_DP
      CASE ( 5 ) ; COSNPHI = ( ( 16.0_DP * COSPHI2 - 20.0_DP ) * COSPHI2 + 5.0_DP ) * COSPHI
      CASE ( 6 ) ; COSNPHI = ( ( 32.0_DP * COSPHI2 - 48.0_DP ) * COSPHI2 + 18.0_DP ) * COSPHI2 - 1.0_DP
      END SELECT

      ! . Calculate the dihedral's contribution to the energy.
      EDIHEDRAL = EDIHEDRAL + xDIHEDRALS(IDIHE)%FC * ( 1.0_DP + xDIHEDRALS(IDIHE)%PHASE * COSNPHI )

   END DO

   ! . Set the function value.
   QUICKDOCK_ENERGY_DIHEDRAL = EDIHEDRAL

   END FUNCTION QUICKDOCK_ENERGY_DIHEDRAL

   !---------------------------------------------------------------------
   FUNCTION QUICKDOCK_INTERACTIONS ( LISTI, LISTJ, CHARGE, LJREP, LJVDW )
   !---------------------------------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: QUICKDOCK_INTERACTIONS

   ! . Array arguments.
   INTEGER,            DIMENSION(1:NDATOMS+1), INTENT(IN) :: LISTI
   INTEGER,            DIMENSION(:),           INTENT(IN) :: LISTJ
   REAL ( KIND = DP ), DIMENSION(1:NDATOMS),   INTENT(IN) :: CHARGE, LJREP, LJVDW

   ! . Local scalars.
   INTEGER            :: IATOM, IINT, JATOM
   REAL ( KIND = DP ) :: AI, AIJ, BI, BIJ, ENB, EPSFAC, QI, QIJ, SIJ2, SIJ6, SIJ12

   ! . Initialize the energies.
   ENB = 0.0_DP

   ! . Calculate the conversion factor for the electrostatic interactions.
   EPSFAC = ELECT_CONST / EPSILON

   ! . Outer loop over atoms.
   DO IATOM = 1,NDATOMS

      ! . Get some information for the atom.
      AI = LJREP(IATOM)
      BI = LJVDW(IATOM)
      QI = EPSFAC * CHARGE(IATOM)

      ! . Loop over the interactions for the atom.
      DO IINT = LISTI(IATOM)+1,LISTI(IATOM+1)

         ! . Get the second atom of the interaction.
         JATOM = LISTJ(IINT)

         ! . Get the inverse distance squared between atoms.
         SIJ2 = 1.0_DP / MAX ( MINDIST2, SUM ( ( DATMCRD(1:3,IATOM) - DATMCRD(1:3,JATOM) )**2 ) )

         ! . Get some data for the IATOM/JATOM interaction.
         AIJ = AI * LJREP(JATOM)
         BIJ = BI * LJVDW(JATOM)
         QIJ = QI * CHARGE(JATOM)

         ! . Calculate some distance factors.
         SIJ6 = SIJ2 * SIJ2 * SIJ2 ; SIJ12 = SIJ6 * SIJ6

         ! . Calculate the electrostatic and Lennard-Jones interactions.
         ENB = ENB + QIJ * SIJ2 + AIJ * SIJ12 - BIJ * SIJ6

      END DO
   END DO

   ! .Set the function value.
   QUICKDOCK_INTERACTIONS = ENB

   END FUNCTION QUICKDOCK_INTERACTIONS

!-------------------------------------------------------------------------------
! . Private grid procedures.
!-------------------------------------------------------------------------------

   !---------------------------------
   SUBROUTINE QUICKDOCK_GRID_GENERATE
   !---------------------------------

   ! . Local scalars.
   INTEGER            :: IATOM, IX, IY, IZ
   REAL ( KIND = DP ) :: AI, BI, EPSFAC, QI, SIJ2, SIJ6, SIJ12

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: R

   ! . Initialize the grid potentials.
   GRIDQ = 0.0_DP ; GRIDR = 0.0_DP ; GRIDV = 0.0_DP

   ! . Get the dielectric factor.
   EPSFAC = ELECT_CONST / EPSILON

   ! . Loop over the fixed atoms.
   DO IATOM = 1,NRATOMS

      ! . Skip moving atoms.
      IF ( .NOT. QRFIXED(IATOM) ) CYCLE

      ! . Set some factors for the atom.
      AI = ATMEPS(IATOM) * ATMSIG(IATOM)**12
      BI = ATMEPS(IATOM) * ATMSIG(IATOM)**6
      QI = EPSFAC * ATMCHG(IATOM)

      ! . Loop over the grid points and get their coordinates.
      DO IZ = 1,NGRID(3)
         R(3) = GRIDORIGIN(3) + REAL ( IZ - 1, DP ) * GRIDSPACING
         DO IY = 1,NGRID(2)
            R(2) = GRIDORIGIN(2) + REAL ( IY - 1, DP ) * GRIDSPACING
	    DO IX = 1,NGRID(1)
	       R(1) = GRIDORIGIN(1) + REAL ( IX - 1, DP ) * GRIDSPACING

               ! . Get the inverse distance squared between the point and the atom.
	       SIJ2 = 1.0_DP / MAX ( MINDIST2, SUM ( ( REFCRD(1:3,IATOM) - R )**2 ) )

               ! . Get some distance factors.
	       SIJ6 = SIJ2 * SIJ2 * SIJ2 ; SIJ12 = SIJ6 * SIJ6

               ! . Set the potentials.
               GRIDQ(IX,IY,IZ) = GRIDQ(IX,IY,IZ) + QI * SIJ2
               GRIDR(IX,IY,IZ) = GRIDR(IX,IY,IZ) + AI * SIJ12
               GRIDV(IX,IY,IZ) = GRIDV(IX,IY,IZ) + BI * SIJ6

	    END DO
	 END DO
      END DO

   END DO

   ! . Set the QGRID flag.
   QGRID = .TRUE.

   END SUBROUTINE QUICKDOCK_GRID_GENERATE

   !--------------------------------------
   SUBROUTINE QUICKDOCK_GRID_READ ( FILE )
   !--------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Local scalars.
   INTEGER :: UNIT

   ! . Get the next file number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, STATUS = "UNKNOWN" )

   ! . READ out the header data.
   READ ( UNIT, "(F25.15)"  ) GRIDSPACING
   READ ( UNIT, "(3I25)"    ) NGRID
   READ ( UNIT, "(3F25.15)" ) GRIDEXTENT
   READ ( UNIT, "(3F25.15)" ) GRIDORIGIN

   ! . Allocate the grid arrays.
   ALLOCATE ( GRIDQ(1:NGRID(1),1:NGRID(2),1:NGRID(3)), &
              GRIDR(1:NGRID(1),1:NGRID(2),1:NGRID(3)), &
              GRIDV(1:NGRID(1),1:NGRID(2),1:NGRID(3))  )

   ! . READ out the grid data.
   READ ( UNIT, "(5G25.15)" ) GRIDQ
   READ ( UNIT, "(5G25.15)" ) GRIDR
   READ ( UNIT, "(5G25.15)" ) GRIDV

   ! . Close the file.
   CLOSE ( UNIT )

   ! . Set the QGRID flag.
   QGRID = .TRUE.

   ! . Write out some data.
   WRITE ( PRINT_LINE, "(A)" ) "Quick-dock grid read from " // TRIM ( FILE ) //"."
   CALL PRINT_PARAGRAPH

   END SUBROUTINE QUICKDOCK_GRID_READ

   !---------------------------------------
   SUBROUTINE QUICKDOCK_GRID_WRITE ( FILE )
   !---------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Local scalars.
   INTEGER :: UNIT

   ! . Check QGRID.
   IF ( .NOT. QGRID ) RETURN

   ! . Get the next file number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, STATUS = "UNKNOWN" )

   ! . Write out the header data.
   WRITE ( UNIT, "(F25.15)"  ) GRIDSPACING
   WRITE ( UNIT, "(3I25)"    ) NGRID
   WRITE ( UNIT, "(3F25.15)" ) GRIDEXTENT
   WRITE ( UNIT, "(3F25.15)" ) GRIDORIGIN

   ! . Write out the grid data.
   WRITE ( UNIT, "(5G25.15)" ) GRIDQ
   WRITE ( UNIT, "(5G25.15)" ) GRIDR
   WRITE ( UNIT, "(5G25.15)" ) GRIDV

   ! . Close the file.
   CLOSE ( UNIT )

   ! . Write out some data.
   WRITE ( PRINT_LINE, "(A)" ) "Quick-dock grid written to " // TRIM ( FILE ) //"."
   CALL PRINT_PARAGRAPH

   END SUBROUTINE QUICKDOCK_GRID_WRITE

!-------------------------------------------------------------------------------
! . Private initialization procedures.
!-------------------------------------------------------------------------------

   !-----------------------------------
   SUBROUTINE QUICKDOCK_INITIALIZE_GRID
   !-----------------------------------

   ! . Scalars.
   QGRID       = .FALSE.
   GRIDSPACING = 0.2_DP

   ! . Allocatable arrays.
   IF ( ALLOCATED ( GRIDQ ) ) DEALLOCATE ( GRIDQ )
   IF ( ALLOCATED ( GRIDR ) ) DEALLOCATE ( GRIDR )
   IF ( ALLOCATED ( GRIDV ) ) DEALLOCATE ( GRIDV )

   ! . Fixed arrays.
   NGRID      = 0
   GRIDEXTENT = 0.0_DP
   GRIDORIGIN = 0.0_DP

   END SUBROUTINE QUICKDOCK_INITIALIZE_GRID

   !--------------------------------------
   SUBROUTINE QUICKDOCK_INITIALIZE_OPTIONS
   !--------------------------------------

   ! . Initialize the quickdock options.
   EOFFGRID             = 60.0E+6_DP
   EPSILON              = 4.0_DP
   ERECMAX              = 10.0_DP
   EXCLUDED_DIHEDRAL    = 15.0_DP * TO_RADIANS
   EXCLUDED_ROTATION    = 7.5_DP * TO_RADIANS
   EXCLUDED_TRANSLATION = 0.5_DP
   MINDIST2             = 0.01_DP**2
   RECOMBINATIONFACTOR  = 0.1_DP
   REDUCEFACTOR         = 3.0_DP

   END SUBROUTINE QUICKDOCK_INITIALIZE_OPTIONS

   !-------------------------------------
   SUBROUTINE QUICKDOCK_INITIALIZE_SYSTEM
   !-------------------------------------

   ! . Scalars.
   NDATOMS           = 0
   NDBONDS           = 0
   NxDIHEDRALS       = 0
   NLIGNB            = 0
   NLIGNB14          = 0
   NRATOMS           = 0
   NROTATABLE        = 0
   NVARIABLES        = 0
   REFATM            = 0
   QROTATION         = .FALSE.
   QTRANSLATION      = .FALSE.
   TRANSLATION_RANGE = 0.0_DP

   ! . Arrays.
   IF ( ALLOCATED ( DATMIND       ) ) DEALLOCATE ( DATMIND       )
   IF ( ALLOCATED ( LISTI         ) ) DEALLOCATE ( LISTI         )
   IF ( ALLOCATED ( LISTI14       ) ) DEALLOCATE ( LISTI14       )
   IF ( ALLOCATED ( LISTJ         ) ) DEALLOCATE ( LISTJ         )
   IF ( ALLOCATED ( LISTJ14       ) ) DEALLOCATE ( LISTJ14       )
   IF ( ALLOCATED ( DBONDS        ) ) DEALLOCATE ( DBONDS        )
   IF ( ALLOCATED ( ROTATABLE     ) ) DEALLOCATE ( ROTATABLE     )
   IF ( ALLOCATED ( QDFIXED       ) ) DEALLOCATE ( QDFIXED       )
   IF ( ALLOCATED ( QRFIXED       ) ) DEALLOCATE ( QRFIXED       )
   IF ( ALLOCATED ( QROTATABLE    ) ) DEALLOCATE ( QROTATABLE    )
   IF ( ALLOCATED ( QTRANSLATABLE ) ) DEALLOCATE ( QTRANSLATABLE )
   IF ( ALLOCATED ( QDIHEDRAL     ) ) DEALLOCATE ( QDIHEDRAL     )
   IF ( ALLOCATED ( DATMCHG       ) ) DEALLOCATE ( DATMCHG       )
   IF ( ALLOCATED ( DATMCHG14     ) ) DEALLOCATE ( DATMCHG14     )
   IF ( ALLOCATED ( DATMREP       ) ) DEALLOCATE ( DATMREP       )
   IF ( ALLOCATED ( DATMREP14     ) ) DEALLOCATE ( DATMREP14     )
   IF ( ALLOCATED ( DATMVDW       ) ) DEALLOCATE ( DATMVDW       )
   IF ( ALLOCATED ( DATMVDW14     ) ) DEALLOCATE ( DATMVDW14     )
   IF ( ALLOCATED ( DATMCRD       ) ) DEALLOCATE ( DATMCRD       )
   IF ( ALLOCATED ( DREFCRD       ) ) DEALLOCATE ( DREFCRD       )
   IF ( ALLOCATED ( REFCRD        ) ) DEALLOCATE ( REFCRD        )
   IF ( ALLOCATED ( xDIHEDRALS     ) ) DEALLOCATE ( xDIHEDRALS     )

   ! . Debug variables.
   NDANGLES    = 0
   NDDIHEDRALS = 0

   ! . Arrays.
   IF ( ASSOCIATED ( DANGLES    ) ) DEALLOCATE ( DANGLES    )
   IF ( ASSOCIATED ( DDIHEDRALS ) ) DEALLOCATE ( DDIHEDRALS )
   IF ( ALLOCATED  ( REFANG     ) ) DEALLOCATE ( REFANG     )
   IF ( ALLOCATED  ( REFBND     ) ) DEALLOCATE ( REFBND     )
   IF ( ALLOCATED  ( REFDIH     ) ) DEALLOCATE ( REFDIH     )

   END SUBROUTINE QUICKDOCK_INITIALIZE_SYSTEM

!-------------------------------------------------------------------------------
! . Private structure procedures.
!-------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   SUBROUTINE STRUCTURES_ORDER ( NSTRUCTURES, NVARIABLES, ENERGIES, VARIABLES )
   !---------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: NSTRUCTURES, NVARIABLES

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(INOUT) :: ENERGIES
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: VARIABLES

   ! . Local scalars.
   INTEGER :: I

   ! . Local arrays.
   INTEGER,            DIMENSION(1:NSTRUCTURES)              :: INDICES
   REAL ( KIND = DP ), DIMENSION(1:NSTRUCTURES)              :: ETEMP
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES,1:NSTRUCTURES) :: VTEMP

   ! . Rank the structures.
   CALL SORT_REAL_RANK ( ENERGIES, INDICES )

   ! . Reorder the energies and variables.
   ETEMP = ENERGIES
   VTEMP = VARIABLES
   DO I = 1,NSTRUCTURES
      ENERGIES (  I) = ETEMP(  INDICES(I))
      VARIABLES(:,I) = VTEMP(:,INDICES(I))
   END DO

   END SUBROUTINE STRUCTURES_ORDER

   !--------------------------------------------------------------------------------
   SUBROUTINE STRUCTURES_READ ( FILE, NTOADD, NVARIN, NSTROUT, ENERGIES, VARIABLES )
   !--------------------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)  :: FILE
   INTEGER,               INTENT(IN)  :: NTOADD, NVARIN
   INTEGER,               INTENT(OUT) :: NSTROUT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   POINTER :: ENERGIES
   REAL ( KIND = DP ), DIMENSION(:,:), POINTER :: VARIABLES

   ! . Local scalars.
   INTEGER :: I, NVAR, UNIT

   ! . Get the next file number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, STATUS = "UNKNOWN" )

   ! . Read in the structure and variable counters.
   READ ( UNIT, "(2I10)" ) NSTROUT, NVAR

   ! . Check NVAR.
   IF ( NVAR /= NVARIN ) CALL PRINT_ERROR ( "STRUCTURES_READ", "Variable number mismatch." )

   ! . Allocate the arrays.
   ALLOCATE ( ENERGIES(1:NSTROUT+NTOADD), VARIABLES(1:NVAR,1:NSTROUT+NTOADD) )

   ! . Read in the data for each structure.
   DO I = 1,NSTROUT
      READ ( UNIT, "(G20.10)" ) ENERGIES(I)
      READ ( UNIT, "(20X,5G20.10)" ) VARIABLES(1:NVAR,I)
   END DO

   ! . Close the file.
   CLOSE ( UNIT )

   END SUBROUTINE STRUCTURES_READ

   !------------------------------------------------------------------------------------------------------
   SUBROUTINE STRUCTURES_REDUCE ( NSTRUCTURES, NVARIABLES, VARIABLES_RANGE_EXCLUDED, ENERGIES, VARIABLES )
   !------------------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN)    :: NVARIABLES
   INTEGER, INTENT(INOUT) :: NSTRUCTURES

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES), INTENT(IN)    :: VARIABLES_RANGE_EXCLUDED
   REAL ( KIND = DP ), DIMENSION(:),            INTENT(INOUT) :: ENERGIES
   REAL ( KIND = DP ), DIMENSION(:,:),          INTENT(INOUT) :: VARIABLES

   ! . Local scalars.
   INTEGER :: I, J, N
   LOGICAL :: QNOT

   ! . The first structure automatically goes on the list.
   N = 1

   ! . Outer loop over structures (the ones that need to be checked).
   DO I = 2,NSTRUCTURES

      ! . Initialization.
      QNOT = .TRUE.

      ! . Inner loop over structures (the ones that are unique).
      DO J = 1,N
         IF ( VARIABLES_COMPARE ( VARIABLES(:,I), VARIABLES(:,J), VARIABLES_RANGE_EXCLUDED ) ) THEN
	    QNOT = .FALSE. ; EXIT
	 END IF
      END DO

      ! . Check for no match.
      IF ( QNOT ) THEN
         N = N + 1
	 ENERGIES (  N) = ENERGIES (  I)
	 VARIABLES(:,N) = VARIABLES(:,I)
      END IF

   END DO

   ! . Reset NSTRUCTURES.
   NSTRUCTURES = N

   END SUBROUTINE STRUCTURES_REDUCE

   !--------------------------------------------------------
   SUBROUTINE STRUCTURES_WRITE ( FILE, ENERGIES, VARIABLES )
   !--------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(IN) :: ENERGIES
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: VARIABLES

   ! . Local scalars.
   INTEGER :: I, NSTR, NVAR, UNIT

   ! . Get the number of structures and variables.
   NSTR = SIZE ( VARIABLES, 2 )
   NVAR = SIZE ( VARIABLES, 1 )

   ! . Get the next file number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, STATUS = "UNKNOWN" )

   ! . Write out the structure and variable counters.
   WRITE ( UNIT, "(2I10)" ) NSTR, NVAR

   ! . Write out the data for each structure.
   DO I = 1,NSTR
      WRITE ( UNIT, "(G20.10)" ) ENERGIES(I)
      WRITE ( UNIT, "(20X,5G20.10)" ) VARIABLES(1:NVAR,I)
   END DO

   ! . Close the file.
   CLOSE ( UNIT )

   END SUBROUTINE STRUCTURES_WRITE

!-------------------------------------------------------------------------------
! . Private variable procedures.
!-------------------------------------------------------------------------------

   !--------------------------------------------------------------------
   FUNCTION VARIABLES_COMPARE ( VSET1, VSET2, VARIABLES_RANGE_EXCLUDED )
   !--------------------------------------------------------------------

   ! . Function declarations.
   LOGICAL :: VARIABLES_COMPARE

   !. Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES), INTENT(IN) :: VSET1, VSET2, VARIABLES_RANGE_EXCLUDED

   ! . Local scalars.
   INTEGER :: N

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES) :: VDIFF

   ! . Get the difference in the variables.
   VDIFF = VSET1 - VSET2

   ! . Adjust the differences for rotations to make sure that they fall within the minimum image distance.
   ! . Initialization.
   N = 0

   ! . Rotatable bonds.
   IF ( NROTATABLE > 0 ) THEN
      VDIFF(N+1:N+NROTATABLE) = VDIFF(N+1:N+NROTATABLE) - ( 2.0_DP * PI ) * ANINT ( VDIFF(N+1:N+NROTATABLE) / ( 2.0_DP * PI ), DP )
      N = N + NROTATABLE
   END IF

   ! . Translations.
   IF ( QTRANSLATION ) N = N + 3

   ! . Rotations.
   IF ( QROTATION ) THEN
      VDIFF(N+1:N+3) = VDIFF(N+1:N+3) - ( 2.0_DP * PI ) * ANINT ( VDIFF(N+1:N+3) / ( 2.0_DP * PI ), DP )
      N = N + 3
   END IF

   ! . Check to see if the differences are all less than the excluded range.
   VARIABLES_COMPARE = ALL ( ABS ( VDIFF ) < 0.5_DP * VARIABLES_RANGE_EXCLUDED )

   END FUNCTION VARIABLES_COMPARE

   !-------------------------------------------
   FUNCTION VARIABLES_DISTANCE ( VSET1, VSET2 )
   !-------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: VARIABLES_DISTANCE

   !. Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES), INTENT(IN) :: VSET1, VSET2

   ! . Local scalars.
   INTEGER :: N

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES) :: VDIFF

   ! . Get the difference in the variables.
   VDIFF = VSET1 - VSET2

   ! . Adjust the differences for rotations to make sure that they fall within the minimum image distance.
   ! . Initialization.
   N = 0

   ! . Rotatable bonds.
   IF ( NROTATABLE > 0 ) THEN
      VDIFF(N+1:N+NROTATABLE) = VDIFF(N+1:N+NROTATABLE) - ( 2.0_DP * PI ) * ANINT ( VDIFF(N+1:N+NROTATABLE) / ( 2.0_DP * PI ), DP )
      N = N + NROTATABLE
   END IF

   ! . Translations.
   IF ( QTRANSLATION ) N = N + 3

   ! . Rotations.
   IF ( QROTATION ) THEN
      VDIFF(N+1:N+3) = VDIFF(N+1:N+3) - ( 2.0_DP * PI ) * ANINT ( VDIFF(N+1:N+3) / ( 2.0_DP * PI ), DP )
      N = N + 3
   END IF

   ! . Get the distance between structures.
   VARIABLES_DISTANCE = SQRT ( DOT_PRODUCT ( VDIFF, VDIFF ) )

   END FUNCTION VARIABLES_DISTANCE

   !------------------------------------------------
   FUNCTION VARIABLES_RANGE_INITIALIZE ( QEXCLUDED )
   !------------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:NVARIABLES) :: VARIABLES_RANGE_INITIALIZE

   ! . Scalar arguments.
   LOGICAL, INTENT(IN) :: QEXCLUDED

   ! . Local scalars.
   INTEGER :: N

   ! . Initialization.
   N = 0

   ! . Set the excluded search range.
   IF ( QEXCLUDED ) THEN

      ! . Rotatable angles.
      IF ( NROTATABLE > 0 ) THEN
         VARIABLES_RANGE_INITIALIZE(N+1:N+NROTATABLE) = EXCLUDED_DIHEDRAL
         N = N + NROTATABLE
      END IF

      ! . Translations.
      IF ( QTRANSLATION ) THEN
         VARIABLES_RANGE_INITIALIZE(N+1:N+3) = EXCLUDED_TRANSLATION
         N = N + 3
      END IF

      ! . Rotations.
      IF ( QROTATION ) THEN
         VARIABLES_RANGE_INITIALIZE(N+1:N+3) = EXCLUDED_ROTATION
         N = N + 3
      END IF

   ! . Set the full search range.
   ELSE

      ! . Rotatable angles.
      IF ( NROTATABLE > 0 ) THEN
         VARIABLES_RANGE_INITIALIZE(N+1:N+NROTATABLE) = 2.0_DP * PI
         N = N + NROTATABLE
      END IF

      ! . Translations.
      IF ( QTRANSLATION ) THEN
         VARIABLES_RANGE_INITIALIZE(N+1:N+3) = 2.0_DP * TRANSLATION_RANGE
         N = N + 3
      END IF

      ! . Rotations.
      IF ( QROTATION ) THEN
         VARIABLES_RANGE_INITIALIZE(N+1:N+3) = 2.0_DP * PI
         N = N + 3
      END IF

   END IF

   END FUNCTION VARIABLES_RANGE_INITIALIZE

!-------------------------------------------------------------------------------
! . Public Monte Carlo procedures.
!-------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE QUICKDOCK_MONTE_CARLO ( OLD_STRUCTURES, NEW_STRUCTURES, NOPTIMIZATIONS, NEQUILIBRATION, NANNEALING, NCONFORMATIONS, &
                                      TEMPERATURE_INITIAL, ACCEPTANCE, ADJUSTMENT_FREQUENCY, MOVES, TEMPERATURE_FINAL,            &
				      TEMPERATURE_SCALING )
   !-------------------------------------------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: OLD_STRUCTURES, NEW_STRUCTURES
   INTEGER,               INTENT(IN) :: NOPTIMIZATIONS, NEQUILIBRATION, NANNEALING, NCONFORMATIONS
   REAL ( KIND = DP ),    INTENT(IN) :: TEMPERATURE_INITIAL

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: TEMPERATURE_SCALING
   INTEGER,               INTENT(IN), OPTIONAL :: ADJUSTMENT_FREQUENCY
   REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: ACCEPTANCE, TEMPERATURE_FINAL

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT), OPTIONAL :: MOVES

   ! . Local scalars.
   CHARACTER ( LEN = 16 ) :: TOPTION
   INTEGER                :: ADJFRQ, I, IBLOCK, ICONF, IMOVE, IOPTIMIZATION, ISTRUCTURE, NBLOCK, NPARAMETERS, NREJECTB, &
                             NSTROLD, NSTRUCTURES
   REAL ( KIND = DP )     :: ACCRAT, FBEST, FCURRENT, TEMPERATURE, TSCALE, TSTART, TSTOP

   ! . Local arrays.
   INTEGER,            ALLOCATABLE, DIMENSION(:)            :: NREJECT, NREJECTT, NTRY, NTRYT
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)            :: MOVE_SIZES, PARAMETERS, VARIABLES_BEST
   REAL ( KIND = DP ),              DIMENSION(:),   POINTER :: ENERGY_NEW, ENERGY_OLD
   REAL ( KIND = DP ),              DIMENSION(:,:), POINTER :: VARIABLES_NEW, VARIABLES_OLD

   !---------------------------------------------------------------------------
   ! . Check the input options.
   !---------------------------------------------------------------------------
   ! . Check that there are variables.
   IF ( NVARIABLES <= 0 ) RETURN

   ! . Check that there is a grid.
   IF ( .NOT. QGRID ) RETURN

   ! . Get the old structures.
   CALL STRUCTURES_READ ( OLD_STRUCTURES, 0, NVARIABLES, NSTROLD, ENERGY_OLD, VARIABLES_OLD )

   ! . Get and check the total number of structures to be stored.
   NSTRUCTURES = ( NOPTIMIZATIONS + 1 ) * NSTROLD
   IF ( NSTRUCTURES <= NSTROLD ) RETURN

   ! . Allocate space for the new structures.
   ALLOCATE ( ENERGY_NEW(1:NSTRUCTURES), VARIABLES_NEW(1:NVARIABLES,1:NSTRUCTURES) )

   ! . Save the old structures.
   ENERGY_NEW   (  1:NSTROLD) = ENERGY_OLD   (  1:NSTROLD) / GAS_CONSTANT
   VARIABLES_NEW(:,1:NSTROLD) = VARIABLES_OLD(:,1:NSTROLD)

   ! . Reset NSTRUCTURES.
   NSTRUCTURES = NSTROLD

   ! . Get the number of parameters.
   NPARAMETERS = NVARIABLES

   ! . Allocate some temporary arrays.
   ALLOCATE ( MOVE_SIZES(1:NPARAMETERS), NREJECT(1:NPARAMETERS), NREJECTT(1:NPARAMETERS), &
              NTRY(1:NPARAMETERS), NTRYT(1:NPARAMETERS),PARAMETERS(1:NPARAMETERS),        &
              VARIABLES_BEST(1:NPARAMETERS) )

   ! . Initialize the rejection and attempts arrays.
   NREJECT  = 0 ; NTRY  = 0
   NREJECTT = 0 ; NTRYT = 0

   ! . Set the number of blocks.
   NBLOCK = NEQUILIBRATION + NANNEALING

   ! . Check the essential input arguments.
   IF ( ( NBLOCK * NCONFORMATIONS ) <= 0 ) RETURN

   ! . Initialize the temperature.
   TSTART = TEMPERATURE_INITIAL

   ! . Initialize the optional arguments (move parameters).
   ACCRAT = 0.4_DP
   ADJFRQ = 10 * NPARAMETERS

   ! . Initialize the move sizes.
   MOVE_SIZES = 0.1_DP
   CALL QUICKDOCK_MC_ADJUST_MOVE_SIZES ( MOVE_SIZES )

   ! . Initialize the optional arguments (temperature).
   TOPTION = "CONSTANT"
   TSTOP   = TSTART / 1000.0_DP

   ! . Check for optional arguments.
   IF ( PRESENT ( ACCEPTANCE           ) ) ACCRAT     = ACCEPTANCE
   IF ( PRESENT ( ADJUSTMENT_FREQUENCY ) ) ADJFRQ     = ADJUSTMENT_FREQUENCY
   IF ( PRESENT ( TEMPERATURE_FINAL    ) ) TSTOP      = TEMPERATURE_FINAL
   IF ( PRESENT ( TEMPERATURE_SCALING  ) ) TOPTION    = TEMPERATURE_SCALING
   IF ( PRESENT ( MOVES                ) ) MOVE_SIZES = MOVES

   ! . Check the temperature options.
   CALL CHECK_TEMPERATURE_OPTIONS

   ! . Do some printing.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "FF0000" )
   CALL PRINT_SUMMARY_START ( "Quick-Dock Monte Carlo Refinement" )
   WRITE ( PRINT_LINE, "(I14)"   ) NSTROLD             ; CALL PRINT_SUMMARY_ELEMENT ( "Input Structures"     )
   WRITE ( PRINT_LINE, "(I14)"   ) NOPTIMIZATIONS      ; CALL PRINT_SUMMARY_ELEMENT ( "Optimizations"	     )
   WRITE ( PRINT_LINE, "(I14)"   ) NEQUILIBRATION      ; CALL PRINT_SUMMARY_ELEMENT ( "Equilibration Blocks" )
   WRITE ( PRINT_LINE, "(I14)"   ) NANNEALING          ; CALL PRINT_SUMMARY_ELEMENT ( "Annealing Blocks"     )
   WRITE ( PRINT_LINE, "(I14)"   ) NCONFORMATIONS      ; CALL PRINT_SUMMARY_ELEMENT ( "Conformations/Block"  )
   WRITE ( PRINT_LINE, "(G14.4)" ) TSTART              ; CALL PRINT_SUMMARY_ELEMENT ( "Initial Temperature"  )
   WRITE ( PRINT_LINE, "(G14.4)" ) ACCRAT              ; CALL PRINT_SUMMARY_ELEMENT ( "Acceptance Ratio"     )
   WRITE ( PRINT_LINE, "(I14)"   ) ADJFRQ              ; CALL PRINT_SUMMARY_ELEMENT ( "Adjustment Frequency" )
   IF ( NANNEALING > 0 ) THEN
      WRITE ( PRINT_LINE, "(A14)"   ) ADJUSTR ( TOPTION ) ; CALL PRINT_SUMMARY_ELEMENT ( "Temperature Scaling" )
      WRITE ( PRINT_LINE, "(G14.4)" ) TSTOP               ; CALL PRINT_SUMMARY_ELEMENT ( "Final Temperature"   )
   END IF
   CALL PRINT_SUMMARY_STOP

   !----------------------------------------------------------------------------
   ! . Do the simulations.
   !----------------------------------------------------------------------------
   ! . Loop over structures.
   DO ISTRUCTURE = 1,NSTROLD

      ! . Loop over optimizations.
      DO IOPTIMIZATION = 1,NOPTIMIZATIONS

         ! . Set the parameter values.
	 PARAMETERS = VARIABLES_OLD(:,ISTRUCTURE)

         ! . Determine the current function value.
         FCURRENT = QUICKDOCK_MC_ENERGY ( PARAMETERS )

         ! . Initialize the best parameter variables.
	 FBEST          = FCURRENT
	 VARIABLES_BEST = PARAMETERS

	 ! . Initialize the move frequency.
	 IMOVE = 0

         ! . Initialize the temperature.
	 TEMPERATURE = TSTART

	 ! . Loop over the number of blocks.
	 DO IBLOCK = 1,NBLOCK

	    ! . Adjust the temperature.
	    CALL ADJUST_TEMPERATURE

	    ! . Loop over the configurations.
	    DO ICONF = 1,NCONFORMATIONS

               ! . Increment the move number.
               IMOVE = IMOVE + 1

               ! . Move a parameters.
               CALL MOVE_PARAMETER

               ! . Reset the best parameters.
	       IF ( FCURRENT < FBEST ) THEN
	          FBEST          = FCURRENT
	          VARIABLES_BEST = PARAMETERS
	       END IF

               ! . Check to see if the move sizes need to be adjusted.
               IF ( ADJFRQ > 0 ) THEN
        	  IF ( MOD ( IMOVE, ADJFRQ ) == 0 ) CALL ADJUST_MOVE_SIZES ( NREJECT, NTRY )
               END IF

	    END DO

	 END DO

         ! . Increment the number of structures.
	 NSTRUCTURES = NSTRUCTURES + 1

         ! . Save the best energy and variables.
         ENERGY_NEW   (  NSTRUCTURES) = FBEST
         VARIABLES_NEW(:,NSTRUCTURES) = VARIABLES_BEST

      END DO
   END DO

   !----------------------------------------------------------------------------
   ! . Finish up.
   !----------------------------------------------------------------------------
   ! . Calculate the acceptance ratios.
   DO I = 1,NPARAMETERS
      IF ( NTRYT(I) > 0 ) THEN
         PARAMETERS(I) = REAL ( ( NTRYT(I) - NREJECTT(I) ), DP ) / REAL ( NTRYT(I), DP )
      ELSE
         PARAMETERS(I) = 0.0_DP
      END IF
   END DO

   ! . Write out the final move sizes.
   WRITE ( OUTPUT, "(/25('-'),A,24('-'))" ) " Adjusted Move Sizes "
   WRITE ( OUTPUT, "(3I10,2G20.10)" ) ( I, NTRYT(I), NREJECTT(I), PARAMETERS(I), MOVE_SIZES(I), I = 1,NPARAMETERS )
   WRITE ( OUTPUT, "(70('-'))" )

   ! . Save the final move sizes if required.
   IF ( PRESENT ( MOVES ) ) MOVES = MOVE_SIZES

   ! . Reset ENERGY_NEW.
   ENERGY_NEW = ENERGY_NEW * GAS_CONSTANT

   ! . Order the structures.
   CALL STRUCTURES_ORDER ( NSTRUCTURES, NVARIABLES, ENERGY_NEW, VARIABLES_NEW )

   ! . Write out the structures to the new file.
   CALL STRUCTURES_WRITE ( NEW_STRUCTURES, ENERGY_NEW(1:NSTRUCTURES), VARIABLES_NEW(:,1:NSTRUCTURES) )

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( ENERGY_NEW, ENERGY_OLD, MOVE_SIZES, NREJECT, NREJECTT, NTRY, NTRYT, VARIABLES_BEST, VARIABLES_NEW, VARIABLES_OLD )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !---------------------------------------------
      SUBROUTINE ADJUST_MOVE_SIZES ( NREJECT, NTRY )
      !---------------------------------------------

      ! . Array arguments.
      INTEGER, DIMENSION(:), INTENT(INOUT) :: NREJECT, NTRY

      ! . Local parameters.
      REAL ( KIND = DP ), PARAMETER :: DOWN = 0.95_DP, UP = 1.05_DP

      ! . Local scalars.
      INTEGER :: I

      ! . Loop over the parameters in turn.
      DO I = 1,NPARAMETERS
         IF ( NTRY(I) > 0 ) THEN
            IF ( ( REAL ( ( NTRY(I) - NREJECT(I) ), DP ) / REAL ( NTRY(I), DP ) ) > ACCRAT ) THEN
               MOVE_SIZES(I) =   UP * MOVE_SIZES(I)
            ELSE
               MOVE_SIZES(I) = DOWN * MOVE_SIZES(I)
            END IF
         END IF
      END DO

      ! . Reinitialize the adjustment arrays.
      NREJECT = 0 ; NTRY = 0

      ! . Check the move sizes.
      CALL QUICKDOCK_MC_ADJUST_MOVE_SIZES ( MOVE_SIZES )

      END SUBROUTINE ADJUST_MOVE_SIZES

      !----------------------------
      SUBROUTINE ADJUST_TEMPERATURE
      !----------------------------

      ! . Local scalars.
      INTEGER :: ABLOCK

      ! . Get the annealing block number.
      ABLOCK = IBLOCK - NEQUILIBRATION

      ! . Check the block number.
      IF ( ABLOCK <= 0 ) RETURN

      ! . Branch on TOPTION.
      SELECT CASE ( TOPTION )
      CASE ( "CONSTANT   " ) ; TEMPERATURE = TSTART
      CASE ( "EXPONENTIAL" ) ; TEMPERATURE = TSTART * EXP ( TSCALE * REAL ( ABLOCK - 1, DP ) )
      CASE ( "LINEAR     " ) ; TEMPERATURE = TSCALE * REAL ( ABLOCK - 1, DP ) + TSTART 
      CASE DEFAULT ; RETURN
      END SELECT

      END SUBROUTINE ADJUST_TEMPERATURE

      !-----------------------------------
      SUBROUTINE CHECK_TEMPERATURE_OPTIONS
      !-----------------------------------

      ! . Make sure the temperatures are positive.
      TSTART = ABS ( TSTART ) ; TSTOP = ABS ( TSTOP )

      ! . Branch on TOPTION.
      SELECT CASE ( TOPTION )
      CASE ( "CONSTANT   " ) ; TSCALE = 0.0_DP ; TSTOP = TSTART
      CASE ( "EXPONENTIAL" ) ; TSCALE = LOG ( TSTOP / TSTART ) / REAL ( NANNEALING - 1, DP )
      CASE ( "LINEAR     " ) ; TSCALE = ( TSTOP - TSTART )     / REAL ( NANNEALING - 1, DP )
      CASE DEFAULT ; CALL PRINT_ERROR ( "MONTE_CARLO", "Unknown temperature scaling option." )
      END SELECT

      END SUBROUTINE CHECK_TEMPERATURE_OPTIONS

      !------------------------
      SUBROUTINE MOVE_PARAMETER
      !------------------------

      ! . Local scalars.
      INTEGER            :: CHOSEN
      REAL ( KIND = DP ) :: OLDFUNC, OLDPARM

      ! . Choose a parameter to move.
      CHOSEN = INT ( REAL ( NPARAMETERS, DP ) * RANDOM ( ) ) + 1

      ! . Increment NTRY.
      NTRY (CHOSEN) = NTRY (CHOSEN) + 1
      NTRYT(CHOSEN) = NTRYT(CHOSEN) + 1

      ! . Save the current function value.
      OLDFUNC = FCURRENT

      ! . Save the old parameter.
      OLDPARM = PARAMETERS(CHOSEN)

      ! . Calculate the new parameter.
      PARAMETERS(CHOSEN) = 2.0_DP * MOVE_SIZES(CHOSEN) * ( RANDOM ( ) - 0.5_DP )

      ! . Check the new parameter.
      CALL QUICKDOCK_MC_CHECK_VARIABLE ( CHOSEN, PARAMETERS(CHOSEN) )

      ! . Calculate the new function value.
      FCURRENT =  QUICKDOCK_MC_ENERGY ( PARAMETERS )

      ! . Do some debug printing.
      IF ( QDEBUGM ) THEN
         WRITE ( OUTPUT, "(A,2I10,4G20.10)" ) "MC Move>", IMOVE, CHOSEN, MOVE_SIZES(CHOSEN), FCURRENT, OLDFUNC, TEMPERATURE
      END IF

      ! . Check to see if the move is rejected.
      IF ( QUICKDOCK_MC_QREJECT ( ( FCURRENT - OLDFUNC ) / TEMPERATURE ) ) THEN

         ! . Increment NREJECT.
         NREJECT(CHOSEN)  = NREJECT (CHOSEN) + 1
         NREJECTT(CHOSEN) = NREJECTT(CHOSEN) + 1
         NREJECTB         = NREJECTB         + 1

         ! . Reactivate the old configuration.
         FCURRENT           = OLDFUNC
         PARAMETERS(CHOSEN) = OLDPARM

      END IF

      END SUBROUTINE MOVE_PARAMETER

   END SUBROUTINE QUICKDOCK_MONTE_CARLO

!-------------------------------------------------------------------------------
! . Private MC procedures.
!-------------------------------------------------------------------------------

   !----------------------------------------------------
   SUBROUTINE QUICKDOCK_MC_ADJUST_MOVE_SIZES ( MCMOVES )
   !----------------------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: MCMOVES

   ! . Local scalars.
   INTEGER :: N

   ! . Set a lower limit on the move size.
   WHERE ( ABS ( MCMOVES ) < 1.0E-4_DP ) MCMOVES = 1.0E-4_DP

   ! . Make sure the new moves are within range.
   ! . Rotations (-PI/+PI) and translations (-TR/+TR).
   ! . Initialization.
   N = 0

   ! . Rotatable angles.
   IF ( NROTATABLE > 0 ) THEN
      WHERE ( MCMOVES(N+1:N+NROTATABLE) > PI ) MCMOVES(N+1:N+NROTATABLE) = PI
      N = N + NROTATABLE
   END IF

   ! . Translations.
   IF ( QTRANSLATION ) THEN
      WHERE ( MCMOVES(N+1:N+3) > TRANSLATION_RANGE ) MCMOVES(N+1:N+3) = TRANSLATION_RANGE
      N = N + 3
   END IF

   ! . Rotations.
   IF ( QROTATION ) THEN
      WHERE ( MCMOVES(N+1:N+3) > PI ) MCMOVES(N+1:N+3) = PI
      N = N + 3
   END IF

   END SUBROUTINE QUICKDOCK_MC_ADJUST_MOVE_SIZES

   !----------------------------------------------------------
   SUBROUTINE QUICKDOCK_MC_CHECK_VARIABLE ( CHOSEN, VARIABLE )
   !----------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN)    :: CHOSEN
   REAL ( KIND = DP ), INTENT(INOUT) :: VARIABLE

   ! . Local scalars.
   INTEGER :: N

   ! . Make sure the new moves are within range.
   ! . Rotations (-PI/+PI) and translations (-TR/+TR).
   ! . Initialization.
   N = 0

   ! . Rotatable angles.
   IF ( NROTATABLE > 0 ) THEN
      IF ( ( CHOSEN >= ( N + 1 ) ) .AND. ( CHOSEN <= ( N + NROTATABLE ) ) ) THEN
         VARIABLE = VARIABLE - ( 2.0_DP * PI ) * ANINT ( VARIABLE / ( 2.0_DP * PI ), DP ) ; RETURN
      END IF
      N = N + NROTATABLE
   END IF

   ! . Translations.
   IF ( QTRANSLATION ) THEN
      IF ( ( CHOSEN >= ( N + 1 ) ) .AND. ( CHOSEN <= ( N + 3 ) ) ) THEN
	 IF ( VARIABLE > TRANSLATION_RANGE ) THEN
	    VARIABLE =  TRANSLATION_RANGE
	 ELSE IF ( VARIABLE < -TRANSLATION_RANGE ) THEN
	    VARIABLE = -TRANSLATION_RANGE
	 END IF
	 RETURN
      END IF
      N = N + 3
   END IF

   ! . Rotations.
   IF ( QROTATION ) THEN
      IF ( ( CHOSEN >= ( N + 1 ) ) .AND. ( CHOSEN <= ( N + 3 ) ) ) THEN
         VARIABLE = VARIABLE - ( 2.0_DP * PI ) * ANINT ( VARIABLE / ( 2.0_DP * PI ), DP ) ; RETURN
      END IF
      N = N + 3
   END IF

   END SUBROUTINE QUICKDOCK_MC_CHECK_VARIABLE

   !-------------------------------------------
   FUNCTION QUICKDOCK_MC_ENERGY ( MCVARIABLES )
   !-------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: QUICKDOCK_MC_ENERGY

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: MCVARIABLES

   ! . Generate the coordinates of the ligand.
   CALL QUICKDOCK_COORDINATES_LIGAND ( MCVARIABLES )

   ! . Get the energy (in temperature units).
   QUICKDOCK_MC_ENERGY = QUICKDOCK_ENERGY ( ) / GAS_CONSTANT ! . Equivalent to R.

   END FUNCTION QUICKDOCK_MC_ENERGY

   !----------------------------------------
   FUNCTION QUICKDOCK_MC_QREJECT ( DELTAEB )
   !----------------------------------------

   ! . Function declarations.
   LOGICAL :: QUICKDOCK_MC_QREJECT

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: DELTAEB

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: UNDERFLOW = 75.0_DP

   ! . Local scalars.
   LOGICAL :: QACCEPT

   ! . Initialize QACCEPT.
   QACCEPT = .FALSE.

   ! . Check to see if the move has been accepted.
   IF ( DELTAEB < UNDERFLOW ) THEN
      IF ( DELTAEB <= 0.0_DP ) THEN
         QACCEPT = .TRUE.
      ELSE IF ( EXP ( - DELTAEB ) > RANDOM ( ) ) THEN
         QACCEPT = .TRUE.
      END IF
   END IF

   ! . Set QREJECT.
   QUICKDOCK_MC_QREJECT = .NOT. QACCEPT

   END FUNCTION QUICKDOCK_MC_QREJECT

END MODULE QUICKDOCK
