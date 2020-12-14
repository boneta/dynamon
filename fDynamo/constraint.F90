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
!
! . Boundary terms added by Ignacio Galvan.
!
!===============================================================================
!                             The Constraint Module
!===============================================================================
!
! . General subroutines:
!
!   ENERGY_CONSTRAINT                  Calculate the constraint energy and
!                                      derivatives.
!
! . Constraint data:
!
!   NCONSTRAINT                        The number of constraints.
!   CONSTRAINT_FIRST, CONSTRAINT_LAST  The constraint definitions.
!
! . Constraint subroutines:
!
!   CONSTRAINT_DEFINE                  Define a constraint.
!   CONSTRAINT_INITIALIZE              Initialize the constraint data.
!   CONSTRAINT_PRINT                   Print a list of constraints.
!   CONSTRAINT_WRITING_START           Start writing constraint data.
!   CONSTRAINT_WRITING_STOP            Stop writing constraint data.
!   ENERGY_SOFT                        Calculate the soft constraint energy and
!                                      derivatives.
!
! . Constraint point data:
!
!   NPOINTS                            The current number of defined points.
!   POINT_FIRST, POINT_LAST            The point definitions.
!
! . Constraint point subroutines:
!
!   CONSTRAINT_POINT_DEFINE            Define a constraint point.
!   CONSTRAINT_POINT_INITIALIZE        Initialize the constraint point data.
!
! . Tether or harmonic position constraint data:
!
!   QTETHER                            The tether constraint flag.
!   TETHER_FCS                         The tether force constants.
!   TETHER_REFCRD                      The tether reference coordinates.
!
! . Tether subroutines:
!
!   TETHER_DEFINE                      Define the tether constraints.
!   TETHER_INITIALIZE                  Initialize the tether constraints.
!   ENERGY_TETHER                      Calculate the tether energy.
!
! . Notes:
!
!   The constraints are of two general types.
!
!   (a) Tether constraints in which atoms are harmonically bound to a reference
!       position in space.
!
!   (b) Constraints between points defined as weighted sums of atom positions.
!       The allowed types of these constraints are:
!
!       DISTANCE, ANGLE, DIHEDRAL, MULTIPLE_DISTANCE, TRANSFER
!
!       Multiple distance constraints are of the form:
!
!                            Sum_i a_i | R1_i - R2_i |
!
!       whereas transfer constraints are of the form:
!
!                 ( R1 - ( R2 + R3 ) / 2 ) . ( R2 - R3 ) / | R2 - R3 |
!
!===============================================================================
MODULE CONSTRAINT

! .Module declarations.
USE CONSTANTS,      ONLY : TO_DEGREES
USE DEFINITIONS,    ONLY : DP
USE FILES,          ONLY : NEXT_UNIT
USE LINEAR_ALGEBRA, ONLY : CROSS_PRODUCT
USE PRINTING,       ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, &
                           PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, PRINT_SUMMARY_STOP,  &
                           PRINT_TABLE_ELEMENT, PRINT_TABLE_OPTIONS, PRINT_TABLE_START,     &
			   PRINT_TABLE_STOP

USE ATOMS, 	    ONLY : ATMCRD, NATOMS, NFREE

IMPLICIT NONE
PRIVATE
PUBLIC :: BOUNDARY_DEFINE, &
          CONSTRAINT_DEFINE, CONSTRAINT_INITIALIZE, CONSTRAINT_POINT_DEFINE, CONSTRAINT_POINT_INITIALIZE,  &
          CONSTRAINT_PRINT, CONSTRAINT_WRITING_START, CONSTRAINT_WRITING_STOP, ENERGY_CONSTRAINT, QTETHER, & 
	  TETHER_DEFINE, TETHER_INITIALIZE
#ifndef PGPC
SAVE
#endif

! . The point type definition.
TYPE POINT_TYPE
   INTEGER                                   :: NATOMS
   INTEGER,            DIMENSION(:), POINTER :: INDICES
   REAL ( KIND = DP ), DIMENSION(:), POINTER :: WEIGHTS
END TYPE POINT_TYPE

! . The linked point type definition.
TYPE LINKED_POINT_TYPE
   TYPE(POINT_TYPE),        POINTER :: POINT
   TYPE(LINKED_POINT_TYPE), POINTER :: NEXT_POINT
END TYPE LINKED_POINT_TYPE

! . The constraint type definition.
TYPE CONSTRAINT_TYPE
   CHARACTER ( LEN = 32 )                    :: TYPE
   INTEGER                                   :: NPOINTS, UNIT
   LOGICAL                                   :: QWRITE
   REAL ( KIND = DP )                        :: EQ, FC
   REAL ( KIND = DP ), DIMENSION(:), POINTER :: WEIGHTS
   TYPE(POINT_TYPE),   DIMENSION(:), POINTER :: POINTS
   TYPE(CONSTRAINT_TYPE),            POINTER :: NEXT_CONSTRAINT
END TYPE CONSTRAINT_TYPE

! . Point data.
INTEGER :: NPOINTS = 0

! . Type data.
TYPE(LINKED_POINT_TYPE), POINTER :: POINT_FIRST, POINT_LAST

! . Constraint data.
INTEGER :: NCONSTRAINT = 0

! . Type data.
TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT_FIRST, CONSTRAINT_LAST

! . Tether constraint data.
LOGICAL                                         :: QTETHER = .FALSE.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: TETHER_FCS
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: TETHER_REFCRD

! . Boundary.
LOGICAL                            :: QBOUNDARY = .FALSE.
LOGICAL, ALLOCATABLE, DIMENSION(:) :: BOUNDARY_SEL

!===============================================================================
CONTAINS
!===============================================================================

   !----------------------------------------------------------------------
   SUBROUTINE ENERGY_CONSTRAINT ( ECONSTRAINT, VIRIAL, GRADIENT, HESSIAN )
   !----------------------------------------------------------------------
   
   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: ECONSTRAINT
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Local scalars.
   REAL ( KIND = DP ) :: EBOUND, ESOFT, ETETHER

   ! . Initialization.
   ECONSTRAINT = 0.0_DP

   ! . Check to see if any energy terms need to be calculated.
   IF ( ( NCONSTRAINT <= 0 ) .AND. ( .NOT. QBOUNDARY ) .AND. ( .NOT. QTETHER ) ) RETURN

   ! . Check the consistency of the derivative options.
   IF ( PRESENT ( HESSIAN ) ) CALL PRINT_ERROR ( "ENERGY_CONSTRAINT", "Second derivatives unavailable for constraints." )

   ! . Call the various energy subroutines.
   CALL ENERGY_BOUNDARY ( EBOUND,  VIRIAL, GRADIENT )
   CALL ENERGY_SOFT     ( ESOFT,   VIRIAL, GRADIENT )
   CALL ENERGY_TETHER   ( ETETHER, VIRIAL, GRADIENT )

   ! . Sum up the total constraint energy.
   ECONSTRAINT = EBOUND + ESOFT + ETETHER

   END SUBROUTINE ENERGY_CONSTRAINT

!===============================================================================
! . Soft Constraint Subroutines.
!===============================================================================

   !------------------------------------------------------------------
   SUBROUTINE CONSTRAINT_DEFINE ( TYPE, FC, EQ, FILE, PRINT, WEIGHTS )
   !------------------------------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: TYPE
   REAL ( KIND = DP ),    INTENT(IN) :: EQ, FC

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE
   LOGICAL,               INTENT(IN), OPTIONAL :: PRINT

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   INTEGER                          :: I, IOSTAT, NDIST
   LOGICAL                          :: QOK, QPRINT
   TYPE(CONSTRAINT_TYPE),   POINTER :: CLOCAL
   TYPE(LINKED_POINT_TYPE), POINTER :: PCURRENT, PNEXT

   ! . Check to see if there are any atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Initialization.
!   NDIST = 0
   NDIST = 1
   QOK   = .TRUE.

   ! . Branch on the type of constraint.
   SELECT CASE ( TYPE )
   CASE ( "ANGLE"             ) ; QOK = ( NPOINTS == 3 )
   CASE ( "DIHEDRAL"          ) ; QOK = ( NPOINTS == 4 )
   CASE ( "DISTANCE"          ) ; QOK = ( NPOINTS == 2 )
   CASE ( ">DISTANCE"         ) ; QOK = ( NPOINTS == 2 )
   CASE ( "<DISTANCE"         ) ; QOK = ( NPOINTS == 2 )
   CASE ( "MULTIPLE_DISTANCE" ) ; NDIST = NPOINTS / 2 ; QOK = ( NPOINTS >  0 ) .AND. ( 2 * NDIST == NPOINTS )
   CASE ( "TRANSFER"          ) ; QOK = ( NPOINTS == 3 )
   CASE DEFAULT ; CALL PRINT_ERROR ( "CONSTRAINT_DEFINE", "Unknown constraint type." )
   END SELECT

   ! . Check that the number of constraint points is valid.
   IF ( .NOT. QOK ) CALL PRINT_ERROR ( "CONSTRAINT_DEFINE", "The constraint has an invalid number of points." )

   ! . Allocate the new constraint.
   ALLOCATE ( CLOCAL )

   ! . Fill the local type variables.
   CLOCAL%TYPE    = TYPE
   CLOCAL%NPOINTS = NPOINTS
   CLOCAL%UNIT    = -1
   CLOCAL%QWRITE  = .FALSE.
   CLOCAL%EQ      = EQ
   CLOCAL%FC      = FC

   ! . The file argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      CLOCAL%UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( CLOCAL%UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "UNKNOWN", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "CONSTRAINT_DEFINE", "I/O Error on constraint data file.", IOSTAT )

   END IF

   ! . Allocate the weights array.
   ALLOCATE ( CLOCAL%WEIGHTS(1:NDIST) )

   ! . Check for a MULTIPLE_DISTANCE constraint.
   IF ( TYPE == "MULTIPLE_DISTANCE" ) THEN

      ! . The weights argument is present.
      IF ( PRESENT ( WEIGHTS ) ) THEN

         ! . Check the dimension of WEIGHTS.
	 IF ( SIZE ( WEIGHTS ) /= NDIST ) THEN
	    CALL PRINT_ERROR ( "CONSTRAINT_DEFINE", "The WEIGHTS array has the wrong size." )
	 ELSE
	    CLOCAL%WEIGHTS = WEIGHTS
	 END IF

      ! . No weights argument is present.
      ELSE

         ! . Set all the weights to 1.
	 CLOCAL%WEIGHTS = 1.0_DP

      END IF

   ! . Other constraints.
   ELSE

      ! . The weights argument is present.
      IF ( PRESENT ( WEIGHTS ) ) THEN
	    CALL PRINT_ERROR ( "CONSTRAINT_DEFINE", "A WEIGHTS array was present with a non-MULTIPLE_DISTANCE constraint." )
      END IF

   END IF

   ! . Allocate the points in the constraint.
   ALLOCATE ( CLOCAL%POINTS(1:NPOINTS) )

   ! . Loop over the points.
   DO I = 1,NPOINTS

      ! . Get the next constraint in the list.
      IF ( I == 1 ) THEN
         PCURRENT => POINT_FIRST
         PNEXT    => POINT_FIRST%NEXT_POINT
      ELSE
         PCURRENT => PNEXT
         PNEXT    => PNEXT%NEXT_POINT
      END IF

      ! . Set the number of atoms in the point.
      CLOCAL%POINTS(I)%NATOMS = PCURRENT%POINT%NATOMS

      ! . Set the point arrays.
      CLOCAL%POINTS(I)%INDICES => PCURRENT%POINT%INDICES
      CLOCAL%POINTS(I)%WEIGHTS => PCURRENT%POINT%WEIGHTS

      ! . Nullify the constraint point arrays.
      NULLIFY ( PCURRENT%POINT%INDICES, PCURRENT%POINT%WEIGHTS )

      ! . Deallocate the point.
      DEALLOCATE ( PCURRENT%POINT )

      ! . Nullify the pointer in PCURRENT.
      NULLIFY ( PCURRENT%NEXT_POINT )

      ! . Deallocate the constraint.
      DEALLOCATE ( PCURRENT )

  END DO

   ! . Reset the number of constraint points.
   NPOINTS = 0

   ! . Nullify the first and last constraint point pointers.
   NULLIFY ( POINT_FIRST, POINT_LAST )

   ! . Nullify the pointer in CLOCAL.
   NULLIFY ( CLOCAL%NEXT_CONSTRAINT )

   ! . Increment the number of constraints.
   NCONSTRAINT = NCONSTRAINT + 1

   ! . This is the first constraint on the list.
   IF ( NCONSTRAINT == 1 ) THEN

      ! . Define the first constraint pointer.
      CONSTRAINT_FIRST => CLOCAL

   ! . There are other constraints in the list.
   ELSE

      ! . Set NEXT_CONSTRAINT of the old last constraint.
      CONSTRAINT_LAST%NEXT_CONSTRAINT => CLOCAL

   END IF

   ! . Update the last constraint pointer.
   CONSTRAINT_LAST => CLOCAL

   ! . Get the print flag.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Write out some information about the constraint.
   IF ( QPRINT ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "0000FF", VARIABLEWIDTH = 18 )
      CALL PRINT_SUMMARY_START ( "Constraint Data" )
      CALL PRINT_SUMMARY_ELEMENT ( "Type", TEXT =  CLOCAL%TYPE(1:18) )
      WRITE ( PRINT_LINE, "(I18)"   ) CLOCAL%UNIT    ; CALL PRINT_SUMMARY_ELEMENT ( "Unit Number"      )
      WRITE ( PRINT_LINE, "(G18.8)" ) CLOCAL%FC      ; CALL PRINT_SUMMARY_ELEMENT ( "Force Constant"   )
      WRITE ( PRINT_LINE, "(G18.8)" ) CLOCAL%EQ      ; CALL PRINT_SUMMARY_ELEMENT ( "Reference Value"  )
      WRITE ( PRINT_LINE, "(I18)"   ) CLOCAL%NPOINTS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Points" )
      WRITE ( PRINT_LINE, "(I18)"   ) SIZE ( CLOCAL%WEIGHTS )
      CALL PRINT_SUMMARY_ELEMENT ( "Distances" )
      CALL PRINT_SUMMARY_STOP
   END IF

   END SUBROUTINE CONSTRAINT_DEFINE

   !-------------------------------
   SUBROUTINE CONSTRAINT_INITIALIZE
   !-------------------------------

   ! . Local scalars.
   INTEGER                        :: I, P
   TYPE(CONSTRAINT_TYPE), POINTER :: CURRENT, NEXT

   ! . Check the number of constraints.
   IF ( NCONSTRAINT <= 0 ) RETURN

   ! . Loop over the constraints.
   DO I = 1,NCONSTRAINT

      ! . Get the next constraint in the list.
      IF ( I == 1 ) THEN
         CURRENT => CONSTRAINT_FIRST
         NEXT    => CONSTRAINT_FIRST%NEXT_CONSTRAINT
      ELSE
         CURRENT => NEXT
         NEXT    => NEXT%NEXT_CONSTRAINT
      END IF

      ! . Loop over the points.
      DO P = 1,CURRENT%NPOINTS

         ! . Free the point arrays.
	 DEALLOCATE ( CURRENT%POINTS(P)%INDICES )
	 DEALLOCATE ( CURRENT%POINTS(P)%WEIGHTS )

      END DO

      ! . Free the constraint arrays.
      DEALLOCATE ( CURRENT%POINTS  )
      DEALLOCATE ( CURRENT%WEIGHTS )

      ! . Nullify the pointer in CURRENT.
      NULLIFY ( CURRENT%NEXT_CONSTRAINT )

      ! . Deallocate the constraint.
      DEALLOCATE ( CURRENT )

  END DO

   ! . Initialize the number of constraints.
   NCONSTRAINT = 0

   ! . Nullify the first and last constraint pointers.
   NULLIFY ( CONSTRAINT_FIRST, CONSTRAINT_LAST )

   ! . Write out some information.
   CALL PRINT_PARAGRAPH ( TEXT = "All constraints have been removed." )

   END SUBROUTINE CONSTRAINT_INITIALIZE

   !--------------------------
   SUBROUTINE CONSTRAINT_PRINT
   !--------------------------

   ! . Local scalars.
   INTEGER                        :: I
   TYPE(CONSTRAINT_TYPE), POINTER :: CURRENT

   ! . There are constraints.
   IF ( NCONSTRAINT > 0 ) THEN

      ! . Write out the header.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#0000FF", VARIABLEWIDTHS = (/ 20, 10, 20, 20, 10 /) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Summary of Defined Constraints", COLSPAN = 5, HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( COLOR = "#00AAFF", TEXT = "Type",            HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( COLOR = "#00AAFF", TEXT = "Points",          HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( COLOR = "#00AAFF", TEXT = "Reference Value", HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( COLOR = "#00AAFF", TEXT = "Force Constant",  HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( COLOR = "#00AAFF", TEXT = "File",            HEADER = .TRUE. )
 
      ! . Loop over the constraints.
      DO I = 1,NCONSTRAINT

	 ! . Get the next constraint on the list.
	 IF ( I == 1 ) THEN
            CURRENT => CONSTRAINT_FIRST
	 ELSE
            CURRENT => CURRENT%NEXT_CONSTRAINT
	 END IF

	 ! . Write out information about the constraint.
         CALL PRINT_TABLE_ELEMENT ( TEXT = CURRENT%TYPE(1:20) )
	 WRITE ( PRINT_LINE, "(I10)"   ) CURRENT%NPOINTS      ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F20.5)" ) CURRENT%EQ           ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F20.5)" ) CURRENT%FC           ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(L10)"   ) ( CURRENT%UNIT > 0 ) ; CALL PRINT_TABLE_ELEMENT

      END DO

      ! . Write out the terminator.
      CALL PRINT_TABLE_STOP

   ! . There are no constraints.
   ELSE
      CALL PRINT_PARAGRAPH ( TEXT = "No constraints are defined." )
   END IF

   END SUBROUTINE CONSTRAINT_PRINT

   !----------------------------------
   SUBROUTINE CONSTRAINT_WRITING_START
   !----------------------------------

   ! . Local scalars.
   INTEGER                        :: I
   TYPE(CONSTRAINT_TYPE), POINTER :: CURRENT

   ! . Return if there are no constraints.
   IF ( NCONSTRAINT <= 0 ) RETURN

   ! . Set the WRITE flags for all those constraints for which files have been defined.
   DO I = 1,NCONSTRAINT

      ! . Get the next constraint on the list.
      IF ( I == 1 ) THEN
         CURRENT => CONSTRAINT_FIRST
      ELSE
         CURRENT => CURRENT%NEXT_CONSTRAINT
      END IF

      ! . If a file has been defined activate it.
      IF ( CURRENT%UNIT > 0 ) THEN
         CURRENT%QWRITE = .TRUE.
         WRITE ( CURRENT%UNIT, "(2F20.10)" ) CURRENT%FC, CURRENT%EQ
      END IF
   END DO

   ! . Write out some information.
   CALL PRINT_PARAGRAPH ( TEXT = "Writing of constraint data started." )

   END SUBROUTINE CONSTRAINT_WRITING_START

   !-------------------------------------------
   SUBROUTINE CONSTRAINT_WRITING_STOP ( CLOSE )
   !-------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: CLOSE

   ! . Local scalars.
   INTEGER                        :: I
   LOGICAL                        :: QCLOSE
   TYPE(CONSTRAINT_TYPE), POINTER :: CURRENT

   ! . Return if there are no constraints.
   IF ( NCONSTRAINT <= 0 ) RETURN

   ! . Set the QCLOSE flag.
   IF ( PRESENT ( CLOSE ) ) THEN
      QCLOSE = CLOSE
   ELSE
      QCLOSE = .TRUE.
   END IF

   ! . Loop over the constraints.
   DO I = 1,NCONSTRAINT

      ! . Get the next constraint on the list.
      IF ( I == 1 ) THEN
         CURRENT => CONSTRAINT_FIRST
      ELSE
         CURRENT => CURRENT%NEXT_CONSTRAINT
      END IF

      ! . A file has been defined.
      IF ( CURRENT%UNIT > 0 ) THEN

         ! . Turn off printing.
         CURRENT%QWRITE = .FALSE.

         ! . Close the file.
         IF ( QCLOSE ) THEN
            CLOSE ( CURRENT%UNIT )
            CURRENT%UNIT = -1
         END IF

      END IF

   END DO

   ! . Write out some information.
   IF ( QCLOSE ) THEN
      CALL PRINT_PARAGRAPH ( TEXT = "Writing of constraint data stopped and files closed." )
   ELSE
      CALL PRINT_PARAGRAPH ( TEXT = "Writing of constraint data stopped but files remain open." )
   END IF

   END SUBROUTINE CONSTRAINT_WRITING_STOP

   !-------------------------------------------------
   SUBROUTINE ENERGY_SOFT ( ESOFT, VIRIAL, GRADIENT )
   !-------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: ESOFT
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL
   
   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: GRADIENT

   ! . Local scalars.
   INTEGER                        :: I, IC, N, P
   LOGICAL                        :: QGRADIENT
   REAL ( KIND = DP )             :: VALUE, W
   TYPE(CONSTRAINT_TYPE), POINTER :: CURRENT

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: CENTERS, DEDC

   ! . Initialize the constraint energy.
   ESOFT = 0.0_DP

   ! . Check whether there are constraints.
   IF ( NCONSTRAINT <= 0 ) RETURN

   ! . Set the gradient flag.
   QGRADIENT = PRESENT ( GRADIENT )

   ! . Loop over the constraints.
   DO IC = 1,NCONSTRAINT

      ! . Get the next constraint on the list.
      IF ( IC == 1 ) THEN
         CURRENT => CONSTRAINT_FIRST
      ELSE
         CURRENT => CURRENT%NEXT_CONSTRAINT
      END IF

      ! . Allocate space for the points and their derivatives.
      ALLOCATE ( CENTERS(1:3,1:CURRENT%NPOINTS) )
      IF ( QGRADIENT ) ALLOCATE ( DEDC(1:3,1:CURRENT%NPOINTS) )

      ! . Calculate the positions of the points.
      CENTERS = 0.0_DP
      DO P = 1,CURRENT%NPOINTS
         DO I = 1,CURRENT%POINTS(P)%NATOMS
            N              = CURRENT%POINTS(P)%INDICES(I)
            CENTERS(1:3,P) = CENTERS(1:3,P) + CURRENT%POINTS(P)%WEIGHTS(I) * ATMCRD(1:3,N)
         END DO
      END DO

      ! . Branch on the type of constraint.
      SELECT CASE ( CURRENT%TYPE )
      CASE ( "ANGLE"             ) ; CALL ENERGY_ANGLE
      CASE ( "DIHEDRAL"          ) ; CALL ENERGY_DIHEDRAL
      CASE ( "DISTANCE"          ) ; CALL ENERGY_DISTANCE
      CASE ( "MULTIPLE_DISTANCE" ) ; CALL ENERGY_MULTIPLE
      CASE ( "TRANSFER"          ) ; CALL ENERGY_TRANSFER
      CASE ( ">DISTANCE"         ) ; CALL ENERGY_BTDISTANCE
      CASE ( "<DISTANCE"         ) ; CALL ENERGY_LTDISTANCE
      END SELECT

      ! . Write out the coordinate value if necessary.
      IF ( CURRENT%QWRITE ) WRITE ( CURRENT%UNIT, "(F20.10)" ) VALUE

      ! . Deallocate the center array.
      DEALLOCATE ( CENTERS )

      ! . Convert the gradients with respect to points to those with respect to atoms.
      IF ( QGRADIENT ) THEN

         ! . Loop over the points.
	 DO P = 1,CURRENT%NPOINTS
            DO I = 1,CURRENT%POINTS(P)%NATOMS
               N               = CURRENT%POINTS(P)%INDICES(I)
	       W               = CURRENT%POINTS(P)%WEIGHTS(I)
               GRADIENT(1:3,N) = GRADIENT(1:3,N) + W * DEDC(1:3,P)

               ! . Determine the virial. Note that this expression is OK
	       ! . as the result is invariant to an arbitrary translation
	       ! . and the minimum image convention has not been applied.
	       VIRIAL = VIRIAL + W * DOT_PRODUCT ( ATMCRD(1:3,N), DEDC(1:3,P) )
	    END DO
	 END DO

         ! . Deallocate the derivative array.
	 DEALLOCATE ( DEDC )

      END IF
   END DO

   !============================================================================
   CONTAINS
   !============================================================================

      !----------------------
      SUBROUTINE ENERGY_ANGLE
      !----------------------

      ! . Parameter declarations.
      REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP - 1.0E-6_DP

      ! . Local scalars.
      REAL ( KIND = DP ) :: DF, DIFF, DOTFAC, RIJ, RKJ

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ, DRKJ, DTI, DTJ, DTK

      ! . Get the interatomic vectors.
      DRIJ = CENTERS(1:3,1) - CENTERS(1:3,2)
      DRKJ = CENTERS(1:3,3) - CENTERS(1:3,2)

      ! . Calculate their sizes.
      RIJ = SQRT ( DOT_PRODUCT ( DRIJ, DRIJ ) )
      RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

      ! . Normalize the vectors.
      DRIJ = DRIJ / RIJ
      DRKJ = DRKJ / RKJ

      ! . Calculate the dot product between the two vectors.
      DOTFAC = DOT_PRODUCT ( DRIJ, DRKJ )

      ! . Ensure DOTFAC is between -1 and 1.
      DOTFAC = SIGN ( MIN ( ABS ( DOTFAC ), DOT_LIMIT ), DOTFAC )

      ! . Get the angle.
      VALUE = ACOS ( DOTFAC ) * TO_DEGREES

      ! . Calculate some intermediate factors.
      DIFF = VALUE - CURRENT%EQ
      DF   = CURRENT%FC * DIFF

      ! . Calculate the contribution to the energy.
      ESOFT = ESOFT + 0.5_DP * DF * DIFF

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) RETURN

      ! . Calculate DF for the first derivative calculation.
      DF = - TO_DEGREES * DF / SQRT ( 1.0_DP - DOTFAC * DOTFAC )

      ! . Calculate the components of the derivatives.
      DTI  = ( DRKJ - DOTFAC * DRIJ ) / RIJ
      DTK  = ( DRIJ - DOTFAC * DRKJ ) / RKJ
      DTJ  = - ( DTI + DTK )

      ! . Calculate the gradients.
      DEDC(1:3,1) = DF * DTI
      DEDC(1:3,2) = DF * DTJ
      DEDC(1:3,3) = DF * DTK

      END SUBROUTINE ENERGY_ANGLE

      !-------------------------
      SUBROUTINE ENERGY_DIHEDRAL
      !-------------------------

      ! . Parameter declarations.
      REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP !- 1.0E-6_DP

      ! . Local scalars.
      REAL ( KIND = DP ) :: DF, DIFF, DOTFAC, FACTIJ, FACTKL, MSIZE, NSIZE, RKJ, SGNFAC

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ, DRKJ, DRKL, DTI, DTJ, DTK, DTL, M, N

      ! . Get the interatomic vectors.
      DRIJ = CENTERS(1:3,1) - CENTERS(1:3,2)
      DRKJ = CENTERS(1:3,3) - CENTERS(1:3,2)
      DRKL = CENTERS(1:3,3) - CENTERS(1:3,4)

      ! . Get M and N.
      M = CROSS_PRODUCT ( DRIJ, DRKJ )
      N = CROSS_PRODUCT ( DRKJ, DRKL )

      ! . Get the size of M and N.
      MSIZE = SQRT ( DOT_PRODUCT ( M, M ) )
      NSIZE = SQRT ( DOT_PRODUCT ( N, N ) )

      ! . Normalize M and N.
      M = M / MSIZE
      N = N / NSIZE

      ! . Ensure DOTFAC is between -1 and 1.
      DOTFAC = DOT_PRODUCT ( M, N )
      DOTFAC = SIGN ( MIN ( ABS ( DOTFAC ), DOT_LIMIT ), DOTFAC )

      ! . Get the sign of the angle.
      SGNFAC = 1.0_DP
      IF ( DOT_PRODUCT ( DRIJ, N ) < 0.0_DP ) SGNFAC = -1.0_DP

      ! . Determine the dihedral.
      VALUE = SGNFAC * TO_DEGREES * ACOS ( DOTFAC )

      ! . Calculate some intermediate factors shifting DIFF so that its value is a minimum wrt translations of 2 Pi.
      DIFF = VALUE - CURRENT%EQ
! -- orig
!      DIFF = DIFF  - 360.0_DP * ANINT ( DIFF / 360.0_DP )
! -- xexo
      IF ( DIFF < -180._DP ) THEN
          VALUE = VALUE + 360._DP
          DIFF  = DIFF  + 360._DP
      END IF
      IF ( DIFF > 180._DP ) THEN
          VALUE = VALUE - 360._DP
          DIFF  = DIFF  - 360._DP
      END IF
! ---------------------------------------------------------
      DF   = CURRENT%FC * DIFF

      ! . Calculate the contribution to the energy.
      ESOFT = ESOFT + 0.5_DP * DF * DIFF

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) RETURN

      ! . Scale DF.
      DF = TO_DEGREES * DF

      ! . Calculate RKJ.
      RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

      ! . Calculate DTI and DTL.
      DTI =   DF * RKJ * M / MSIZE
      DTL = - DF * RKJ * N / NSIZE

      ! . Calculate some additional factors.
      FACTIJ = DOT_PRODUCT ( DRIJ, DRKJ ) / ( RKJ * RKJ )
      FACTKL = DOT_PRODUCT ( DRKL, DRKJ ) / ( RKJ * RKJ )

      ! . Calculate DTJ and DTK.
      DTJ = DTI * ( FACTIJ - 1.0_DP ) - FACTKL * DTL
      DTK = DTL * ( FACTKL - 1.0_DP ) - FACTIJ * DTI

      ! . Calculate the gradients.
      DEDC(1:3,1) = DTI
      DEDC(1:3,2) = DTJ
      DEDC(1:3,3) = DTK
      DEDC(1:3,4) = DTL

      END SUBROUTINE ENERGY_DIHEDRAL

      !-------------------------
      SUBROUTINE ENERGY_DISTANCE
      !-------------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: DF, DIFF

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      ! . Calculate the distance between the centers.
      DR    = CENTERS(1:3,1) - CENTERS(1:3,2)
      VALUE = SQRT ( DOT_PRODUCT ( DR, DR ) )

      ! . Calculate some intermediate factors.
      DIFF  = VALUE - CURRENT%EQ
      DF    = CURRENT%FC * DIFF

      ! . Calculate the contribution to the energy.
      ESOFT = ESOFT + 0.5_DP * DF * DIFF

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) RETURN

      ! . Calculate the gradients.
      DF = DF / VALUE
      DEDC(1:3,1) =   DF * DR
      DEDC(1:3,2) = - DF * DR

      END SUBROUTINE ENERGY_DISTANCE

      !----------------------------
      SUBROUTINE ENERGY_BTDISTANCE
      !----------------------------
      REAL ( KIND = DP ) :: DF, DIFF
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      DEDC  = .0_DP
      DR    = CENTERS(1:3,1) - CENTERS(1:3,2)
      VALUE = SQRT ( DOT_PRODUCT ( DR, DR ) )
      IF( VALUE <= CURRENT%EQ ) RETURN
      DIFF  = VALUE - CURRENT%EQ
      DF    = CURRENT%FC * DIFF
      ESOFT = ESOFT + 0.5_DP * DF * DIFF
      IF ( .NOT. QGRADIENT ) RETURN
      DF = DF / VALUE
      DEDC(1:3,1) =   DF * DR
      DEDC(1:3,2) = - DF * DR
      END SUBROUTINE ENERGY_BTDISTANCE

      !----------------------------
      SUBROUTINE ENERGY_LTDISTANCE
      !----------------------------
      REAL ( KIND = DP ) :: DF, DIFF
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      DEDC  = .0_DP
      DR    = CENTERS(1:3,1) - CENTERS(1:3,2)
      VALUE = SQRT ( DOT_PRODUCT ( DR, DR ) )
      IF( VALUE >= CURRENT%EQ ) RETURN
      DIFF  = VALUE - CURRENT%EQ
      DF    = CURRENT%FC * DIFF
      ESOFT = ESOFT + 0.5_DP * DF * DIFF
      IF ( .NOT. QGRADIENT ) RETURN
      DF = DF / VALUE
      DEDC(1:3,1) =   DF * DR
      DEDC(1:3,2) = - DF * DR
      END SUBROUTINE ENERGY_LTDISTANCE


      !-------------------------
      SUBROUTINE ENERGY_MULTIPLE
      !-------------------------

      ! . Local scalars.
      INTEGER            :: I, II, NDIST
      REAL ( KIND = DP ) :: DF, DIFF

      ! . Local arrays.
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: R
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: DR

      ! . Get the number of distances involved in the constraint.
      NDIST = SIZE ( CURRENT%WEIGHTS )

      ! . Allocate space for the distances and intercenter vectors.
      ALLOCATE ( DR(1:3,NDIST), R(1:NDIST) )

      ! . Calculate the distances between the centers.
      DO I = 1,NDIST
         II        = 2 * ( I - 1 )
         DR(1:3,I) = CENTERS(1:3,II+1) - CENTERS(1:3,II+2)
	 R(I)      = SQRT ( DOT_PRODUCT ( DR(1:3,I), DR(1:3,I) ) )
      END DO

      ! . Calculate the value of the constraint.
      VALUE = DOT_PRODUCT ( CURRENT%WEIGHTS, R )

      ! . Calculate some intermediate factors.
      DIFF  = VALUE - CURRENT%EQ
      DF    = CURRENT%FC * DIFF

      ! . Calculate the contribution to the energy.
      ESOFT = ESOFT + 0.5_DP * DF * DIFF

      ! . Check for a gradient calculation.
      IF ( QGRADIENT ) THEN

	 ! . Calculate the derivatives.
	 DO I = 1,NDIST
            II             = 2 * ( I - 1 )
	    DEDC(1:3,II+1) =   CURRENT%WEIGHTS(I) * DF * DR(1:3,I) / R(I)
	    DEDC(1:3,II+2) = - CURRENT%WEIGHTS(I) * DF * DR(1:3,I) / R(I)
	 END DO
 
      END IF
 
      ! . Deallocate any temporary space.
      DEALLOCATE ( DR, R )

      END SUBROUTINE ENERGY_MULTIPLE

      !-------------------------
      SUBROUTINE ENERGY_TRANSFER
      !-------------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: DF, DIFF, R23

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR1M, DR23, DV

      ! . Calculate the vector between the transfering atom and the mid-point between the transfering centers.
      DR1M = CENTERS(1:3,1) - 0.5_DP * ( CENTERS(1:3,2) + CENTERS(1:3,3) )

      ! . Calculate the distance between the transfering centers.
      DR23 = CENTERS(1:3,2) - CENTERS(1:3,3)
      R23  = SQRT ( DOT_PRODUCT ( DR23, DR23 ) )

      ! . Normalize DR23.
      DR23 = DR23 / R23

      ! . Calculate the constraint.
      VALUE = DOT_PRODUCT ( DR1M, DR23 )

      ! . Calculate some intermediate factors.
      DIFF  = VALUE - CURRENT%EQ
      DF    = CURRENT%FC * DIFF

      ! . Calculate the contribution to the energy.
      ESOFT = ESOFT + 0.5_DP * DF * DIFF

      ! . The gradients are to be calculated.
      IF ( .NOT. QGRADIENT ) RETURN

      ! . Calculate some derivative factors.
      DV = ( DR1M - VALUE * DR23 ) / R23

      ! . Calculate the derivatives.
      DEDC(1:3,1) =   DF * DR23
      DEDC(1:3,2) =   DF * ( DV - 0.5_DP * DR23 )
      DEDC(1:3,3) = - DF * ( DV + 0.5_DP * DR23 )

      END SUBROUTINE ENERGY_TRANSFER

   END SUBROUTINE ENERGY_SOFT

!===============================================================================
! . Constraint Point Subroutines.
!===============================================================================

   !---------------------------------------------------------------
   SUBROUTINE CONSTRAINT_POINT_DEFINE ( SELECTION, WEIGHTS, PRINT )
   !---------------------------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Array arguments.
   LOGICAL,            DIMENSION(:), INTENT(IN)           :: SELECTION
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   INTEGER                          :: I, N, NSELECTED
   LOGICAL                          :: QPRINT, QWEIGHTS
   TYPE(LINKED_POINT_TYPE), POINTER :: PLOCAL

   ! . Check to see if there are any atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Check the selection array size.
   IF ( SIZE ( SELECTION ) /= NATOMS ) CALL PRINT_ERROR ( "CONSTRAINT_POINT_DEFINE", "Invalid selection array size." )

   ! . Set the weights array flag.
   QWEIGHTS = PRESENT ( WEIGHTS )

   ! . Check the weights array size.
   IF ( QWEIGHTS ) THEN
      IF ( SIZE ( WEIGHTS ) /= NATOMS ) CALL PRINT_ERROR ( "CONSTRAINT_POINT_DEFINE", "Invalid weight array size." )
   END IF

   ! . Get the number of atoms in each selection.
   NSELECTED = COUNT ( SELECTION )

   ! . Check the number of atoms.
   IF ( NSELECTED <= 0 ) CALL PRINT_ERROR ( "CONSTRAINT_POINT_DEFINE", "No atoms have been selected." )

   ! . Allocate the new constraint point.
   ALLOCATE ( PLOCAL )

   ! . Allocate the point type.
   ALLOCATE ( PLOCAL%POINT )

   ! . Set the number of atoms in the point.
   PLOCAL%POINT%NATOMS = NSELECTED

   ! . Allocate the constraint point arrays.
   ALLOCATE ( PLOCAL%POINT%INDICES(1:NSELECTED), PLOCAL%POINT%WEIGHTS(1:NSELECTED) )

   ! . Fill the point index array.
   N = 0
   DO I = 1,NATOMS
      IF ( SELECTION(I) ) THEN
         N = N + 1
         PLOCAL%POINT%INDICES(N) = I
      END IF
   END DO

   ! . A weights array is present.
   IF ( QWEIGHTS ) THEN

      ! . Fill the weights array.
      N = 0
      DO I = 1,NATOMS
         IF ( SELECTION(I) ) THEN
            N = N + 1
            PLOCAL%POINT%WEIGHTS(N) = WEIGHTS(I)
         END IF
      END DO

   ! . A weights array is not present.
   ELSE

      ! . All the weights are one.
      PLOCAL%POINT%WEIGHTS = 1.0_DP

   END IF

   ! . Normalize the weights array.
   PLOCAL%POINT%WEIGHTS = PLOCAL%POINT%WEIGHTS / SUM ( PLOCAL%POINT%WEIGHTS )   

   ! . Nullify the pointer in WLOCAL.
   NULLIFY ( PLOCAL%NEXT_POINT )

   ! . Increment the number of constraint points.
   NPOINTS = NPOINTS + 1

   ! . This is the first point on the list.
   IF ( NPOINTS == 1 ) THEN

      ! . Define the first constraint point pointer.
      POINT_FIRST => PLOCAL

   ! . There are other constraint points on the list.
   ELSE

      ! . Set NEXT_POINT of the old last constraint point.
      POINT_LAST%NEXT_POINT => PLOCAL

   END IF

   ! . Update the last constraint point pointer.
   POINT_LAST => PLOCAL

   ! . Get the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Write out some information about the point.
   IF ( QPRINT ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "0000FF", VARIABLEWIDTH = 12 )
      CALL PRINT_SUMMARY_START ( "Constraint Point Data" )
      WRITE ( PRINT_LINE, "(I12)" ) NPOINTS   ; CALL PRINT_SUMMARY_ELEMENT ( "Current Point Number"     )
      WRITE ( PRINT_LINE, "(I12)" ) NSELECTED ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Atoms in Point" )
      CALL PRINT_SUMMARY_STOP
   END IF

   END SUBROUTINE CONSTRAINT_POINT_DEFINE

   !-------------------------------------
   SUBROUTINE CONSTRAINT_POINT_INITIALIZE
   !-------------------------------------

   ! . Local scalars.
   INTEGER                          :: I
   TYPE(LINKED_POINT_TYPE), POINTER :: CURRENT, NEXT

   ! . Check the number of constraint points.
   IF ( NPOINTS <= 0 ) RETURN

   ! . Loop over the constraint points.
   DO I = 1,NPOINTS

      ! . Get the next constraint point on the list.
      IF ( I == 1 ) THEN
         CURRENT => POINT_FIRST
         NEXT    => POINT_FIRST%NEXT_POINT
      ELSE
         CURRENT => NEXT
         NEXT    => NEXT%NEXT_POINT
      END IF

      ! . Free the point arrays.
      DEALLOCATE ( CURRENT%POINT%INDICES )
      DEALLOCATE ( CURRENT%POINT%WEIGHTS )

      ! . Deallocate the point.
      DEALLOCATE ( CURRENT%POINT )

      ! . Nullify the pointer in CURRENT.
      NULLIFY ( CURRENT%NEXT_POINT )

      ! . Deallocate the constraint point.
      DEALLOCATE ( CURRENT )

  END DO

   ! . Initialize the number of constraint points.
   NPOINTS = 0

   ! . Nullify the first and last constraint point pointers.
   NULLIFY ( POINT_FIRST, POINT_LAST )

   ! . Write out some information.
   CALL PRINT_PARAGRAPH ( TEXT = "All constraint points have been removed." )

   END SUBROUTINE CONSTRAINT_POINT_INITIALIZE

!===============================================================================
! . Boundary Constraint Subroutines.
!===============================================================================

   !---------------------------------------
   SUBROUTINE BOUNDARY_DEFINE ( SELECTION )
   !---------------------------------------

   ! . Array arguments.
   LOGICAL, DIMENSION(:), INTENT(IN) :: SELECTION

   ! . Local scalars.
   INTEGER :: NUM

   ! . Check the selection array size.
   IF ( SIZE ( SELECTION ) /= NATOMS ) CALL PRINT_ERROR ( "BOUNDARY_DEFINE", "Invalid selection array size." )

   ! . Set the boundary selection
   IF ( ALLOCATED ( BOUNDARY_SEL ) ) DEALLOCATE ( BOUNDARY_SEL )
   ALLOCATE ( BOUNDARY_SEL ( NATOMS ) )
   BOUNDARY_SEL = SELECTION
   NUM = COUNT ( SELECTION )
   QBOUNDARY = ( NUM > 0 )

   IF ( QBOUNDARY ) THEN
     WRITE ( PRINT_LINE, '(I10)' ) NUM
     CALL PRINT_PARAGRAPH ( TEXT = "Boundary constraints added for " // TRIM ( ADJUSTL ( PRINT_LINE ) ) // " atoms." )
   ELSE
     CALL PRINT_PARAGRAPH ( TEXT = "Boundary constraints removed." )
   END IF

   END SUBROUTINE BOUNDARY_DEFINE

   !-----------------------------------------------------------
   SUBROUTINE ENERGY_BOUNDARY ( ECONSTRAINT, VIRIAL, GRADIENT )
   !-----------------------------------------------------------

   ! . External user-defined function for the boundary
   INTERFACE
     FUNCTION BOUND_FUNC(COORD)
       DOUBLE PRECISION, DIMENSION(4) :: BOUND_FUNC
       DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: COORD
     END FUNCTION
   END INTERFACE

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: ECONSTRAINT
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: GRADIENT

   ! . Local scalars.
   REAL ( KIND = DP ), DIMENSION(4) :: VALUE
   LOGICAL :: QGRADIENT
   INTEGER :: I

   ECONSTRAINT = 0.0_DP
   VIRIAL = 0.0_DP

   IF ( .NOT. QBOUNDARY ) RETURN

   QGRADIENT = PRESENT ( GRADIENT )

   DO I = 1, NATOMS
     IF ( .NOT. BOUNDARY_SEL(I) ) CYCLE
     VALUE = BOUND_FUNC ( ATMCRD(1:3,I) )
     ECONSTRAINT = ECONSTRAINT + VALUE(1)
     IF ( .NOT. QGRADIENT ) CYCLE
     VIRIAL = VIRIAL + DOT_PRODUCT ( ATMCRD(1:3,I), VALUE(2:4) )
     GRADIENT(1:3,I) = GRADIENT(1:3,I) + VALUE(2:4)
   END DO

   END SUBROUTINE ENERGY_BOUNDARY

!===============================================================================
! . Tether Constraint Subroutines.
!===============================================================================

   !--------------------------------------------------------------------------
   SUBROUTINE TETHER_DEFINE ( FORCE_CONSTANT, REFERENCE_STRUCTURE, SELECTION )
   !--------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: FORCE_CONSTANT

   ! . Array arguments.
   LOGICAL,            DIMENSION(1:NATOMS),     INTENT(IN), OPTIONAL :: SELECTION
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(IN), OPTIONAL :: REFERENCE_STRUCTURE

   ! . Initialize the tether constraint data structures.
   CALL TETHER_INITIALIZE

   ! . The force constant is positive so define the constraints.
   IF ( FORCE_CONSTANT > 0.0_DP ) THEN

      ! . Allocate the constraint arrays.
      ALLOCATE ( TETHER_FCS(1:NATOMS), TETHER_REFCRD(1:3,1:NATOMS) )

      ! . The reference structure is present.
      IF ( PRESENT ( REFERENCE_STRUCTURE ) ) THEN
         TETHER_REFCRD = REFERENCE_STRUCTURE
      ! . The reference structure is absent.
      ELSE
         CALL PRINT_ERROR( "TETHER_DEFINE", "The tether constraint reference structure is missing." )
      END IF

      ! . Initialize the force constant array.
      TETHER_FCS = FORCE_CONSTANT

      ! . A selection array is present.
      IF ( PRESENT ( SELECTION ) ) THEN

         ! . Apply the tether constraints only to specified atoms.
         WHERE ( .NOT. SELECTION ) TETHER_FCS = 0.0_DP

      END IF

      ! . Set the tether flag.
      QTETHER = .TRUE.

   END IF

   ! . Write out a summary of the data.
   IF ( QTETHER ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "0000FF" )
      CALL PRINT_SUMMARY_START ( "Summary of Tether Constraint Data" )
      WRITE ( PRINT_LINE, "(I14)"   ) COUNT ( TETHER_FCS > 0.0_DP )
      CALL PRINT_SUMMARY_ELEMENT ( "Number of Constraints" )
      WRITE ( PRINT_LINE, "(G14.8)" ) FORCE_CONSTANT
      CALL PRINT_SUMMARY_ELEMENT ( "Force Constant" )
      CALL PRINT_SUMMARY_STOP
   ELSE
      CALL PRINT_PARAGRAPH ( TEXT = "All harmonic tether constraints on the atoms have been removed." )
   END IF

   END SUBROUTINE TETHER_DEFINE

   !---------------------------
   SUBROUTINE TETHER_INITIALIZE
   !---------------------------

   ! . Initialize the tether flag.
   QTETHER = .FALSE.

   ! . Deallocate the tether arrays.
   IF ( ALLOCATED ( TETHER_FCS    ) ) DEALLOCATE ( TETHER_FCS    )
   IF ( ALLOCATED ( TETHER_REFCRD ) ) DEALLOCATE ( TETHER_REFCRD )

   END SUBROUTINE TETHER_INITIALIZE

   !---------------------------------------------------------
   SUBROUTINE ENERGY_TETHER ( ECONSTRAINT, VIRIAL, GRADIENT )
   !---------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: ECONSTRAINT
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL 
   
   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: GRADIENT

   ! . Local scalars.
   INTEGER            :: I
   LOGICAL            :: QGRADIENT
   REAL ( KIND = DP ) :: K

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DR

   ! . Initialize the constraint energy.
   ECONSTRAINT = 0.0_DP

   ! . Check whether there are constraints.
   IF ( .NOT. QTETHER ) RETURN

   ! . Check the size of the arrays.
   IF ( SIZE ( TETHER_FCS ) /= NATOMS ) CALL PRINT_ERROR ( "ENERGY_TETHER", "Invalid CONSTRAINT arrays." )

   ! . Set the gradient flag.
   QGRADIENT = PRESENT ( GRADIENT )

   ! . Loop over the atoms.
   DO I = 1,NATOMS

      ! . Get the force constant.
      K = TETHER_FCS(I)

      ! . There is a constraint.
      IF ( K > 0.0_DP ) THEN

         ! . Get the displacement vector.
         DR = ATMCRD(1:3,I) - TETHER_REFCRD(1:3,I)

         ! . Calculate the energy.
         ECONSTRAINT = ECONSTRAINT + 0.5_DP * K * DOT_PRODUCT ( DR, DR )

         ! . Check for a gradient calculation.
	 IF ( .NOT. QGRADIENT ) CYCLE

         ! . Calculate the virial.
	 VIRIAL = VIRIAL + K * DOT_PRODUCT ( ATMCRD(1:3,I), DR )

         ! . Calculate the gradient.
         IF ( QGRADIENT ) GRADIENT(1:3,I) = GRADIENT(1:3,I) + K * DR

      END IF
   END DO

   END SUBROUTINE ENERGY_TETHER

END MODULE CONSTRAINT

!===============================================================================
! . Boundary Constraint Function.
!===============================================================================

! . I have to put a default BOUND_FUNC, even if it's not used
! . Pass --allow-multiple-definition to the linker if you provide your own one.

FUNCTION BOUND_FUNC ( COORD )
  USE DEFINITIONS
  REAL ( KIND = DP ), DIMENSION(4)             :: BOUND_FUNC
  REAL ( KIND = DP ), DIMENSION(3), INTENT(IN) :: COORD
  WRITE ( 77, * )
  BOUND_FUNC = 0.0_DP
END FUNCTION BOUND_FUNC
