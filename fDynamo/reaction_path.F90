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
!                            The Reaction Path Module
!===============================================================================
!
! . Subroutines:
!
!   REACTION_PATH_TRACE             Trace half a path from a saddle point.
!
!===============================================================================
MODULE REACTION_PATH

! . Module declarations.
USE DEFINITIONS,           ONLY : DP
USE PRINTING,              ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, &
                                  PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT, PRINT_TABLE_OPTIONS, PRINT_TABLE_START,            &
				  PRINT_TABLE_STOP
USE DIAGONALIZATION,       ONLY : SYMMETRIC_UPPER

USE ATOMS,                 ONLY : ATMCRD, ATMFIX, NATOMS, NFIXED, NFREE
USE DCD_IO
USE NORMAL_MODE_UTILITIES
USE POTENTIAL_ENERGY,      ONLY : ATMDER, ATMHES, ETOTAL, GRADIENT, HESSIAN

IMPLICIT NONE
PRIVATE
PUBLIC :: REACTION_PATH_TRACE

!==============================================================================
CONTAINS
!==============================================================================

   !-------------------------------------------------------------------------------------------------------------------
   SUBROUTINE REACTION_PATH_TRACE ( DIRECTION, MAXSTP, PRINT_FREQUENCY, FROM_SADDLE, USE_MASS_WEIGHTING, ENERGY_STEP, &
                                                       PATH_STEP, SAVE_FREQUENCY, FILE )
   !-------------------------------------------------------------------------------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = 1 ), INTENT(IN) :: DIRECTION
   INTEGER,               INTENT(IN) :: MAXSTP

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE
   INTEGER,               INTENT(IN), OPTIONAL :: PRINT_FREQUENCY, SAVE_FREQUENCY
   LOGICAL,               INTENT(IN), OPTIONAL :: FROM_SADDLE, USE_MASS_WEIGHTING
   REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: ENERGY_STEP, PATH_STEP

   ! . Local scalars.
   INTEGER            :: I, NPRINT, NSAVE
   LOGICAL            :: QMASS, QSADDLE
   REAL ( KIND = DP ) :: DELTA, DIRFAC, ELOWER

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: WEIGHTS

   ! . Local types.
   TYPE(DCD_TYPE) :: TRAJECTORY

   ! . Check the number of free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Check the direction.
   SELECT CASE ( DIRECTION )
   CASE ( "+" ) ; DIRFAC = + 1.0_DP
   CASE ( "-" ) ; DIRFAC = - 1.0_DP
   CASE DEFAULT ; CALL PRINT_ERROR ( "REACTION_PATH_TRACE", "Unknown direction." )
   END SELECT

   ! . Initialize some options.
   NPRINT  = 1
   NSAVE   = 0
   QMASS   = .FALSE.
   QSADDLE = .TRUE.
   DELTA   = 0.1_DP
   ELOWER  = 1.0_DP

   ! . Process the input options.
   IF ( PRESENT ( ENERGY_STEP        ) ) ELOWER  = ENERGY_STEP
   IF ( PRESENT ( FROM_SADDLE        ) ) QSADDLE = FROM_SADDLE
   IF ( PRESENT ( PATH_STEP          ) ) DELTA   = PATH_STEP
   IF ( PRESENT ( PRINT_FREQUENCY    ) ) NPRINT  = PRINT_FREQUENCY
   IF ( PRESENT ( SAVE_FREQUENCY     ) ) NSAVE   = SAVE_FREQUENCY
   IF ( PRESENT ( USE_MASS_WEIGHTING ) ) QMASS   = USE_MASS_WEIGHTING

   ! . Check NPRINT and NSAVE.
   IF ( ( NPRINT <= 0 ) .OR. ( NPRINT > MAXSTP ) ) NPRINT = 0
   IF ( ( NSAVE  <= 0 ) .OR. ( NSAVE  > MAXSTP ) ) NSAVE  = 0

   ! . Check for printing.
   IF ( NPRINT > 0 ) THEN

      ! . Write out the options.
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#008888", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Steepest Descent Reaction Path Tracing" )
      WRITE ( PRINT_LINE, "(I16)"   ) MAXSTP  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Steps"     )
      WRITE ( PRINT_LINE, "(I16)"   ) NPRINT  ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency"     )
      WRITE ( PRINT_LINE, "(L16)"   ) QSADDLE ; CALL PRINT_SUMMARY_ELEMENT ( "From Saddle"         )
      WRITE ( PRINT_LINE, "(L16)"   ) QMASS   ; CALL PRINT_SUMMARY_ELEMENT ( "Mass Weighting"      )
      WRITE ( PRINT_LINE, "(G16.4)" ) ELOWER  ; CALL PRINT_SUMMARY_ELEMENT ( "Initial Energy Step" )
      WRITE ( PRINT_LINE, "(G16.4)" ) DELTA   ; CALL PRINT_SUMMARY_ELEMENT ( "Path Step"           )
      WRITE ( PRINT_LINE, "(G16.4)" ) DIRFAC  ; CALL PRINT_SUMMARY_ELEMENT ( "Direction Factor"    )
      WRITE ( PRINT_LINE, "(I16)"   ) NSAVE   ; CALL PRINT_SUMMARY_ELEMENT ( "Save Frequency"      )
      CALL PRINT_SUMMARY_STOP

   END IF

   ! . Activate the trajectory if necessary.
   IF ( NSAVE > 0 ) THEN

      ! . Check the FILE argument.
      IF ( .NOT. PRESENT ( FILE ) ) CALL PRINT_ERROR ( "REACTION_PATH_TRACE", "Path file name absent." )

      ! . Open the trajectory.
      CALL DCD_INITIALIZE ( TRAJECTORY )
      CALL DCD_ACTIVATE_WRITE ( FILE, TRAJECTORY, "CORD", NATOMS, NFIXED, ( MAXSTP / NSAVE ) + 1, QFIX = ATMFIX )

   END IF

   ! . Fill the WEIGHTS array if necessary.
   IF ( QMASS ) THEN
      DO I = 1,NATOMS
         IF ( ATMFIX(I) ) THEN
	    WEIGHTS(1:3,I) = 0.0_DP
	 ELSE
            WEIGHTS(1:3,I) = 1.0_DP / SQRT ( ATMMAS(I) )
	 END IF
      END DO
   END IF

   ! . Check the saddle point structure and get the initial displaced structure.
   IF ( QSADDLE ) CALL CHECK_SADDLE_POINT

   ! . Calculate the path.
   CALL STEEPEST_DESCENT_PATH

   ! . Deactivate the trajectory.
   IF ( NSAVE > 0 ) CALL DCD_DEACTIVATE ( TRAJECTORY )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !----------------------------
      SUBROUTINE CHECK_SADDLE_POINT
      !----------------------------

      ! . Local scalars.
      INTEGER            :: I, IFREE, II, N
      REAL ( KIND = DP ) :: EIGVAL1, GRMS, STPINI

      ! . Local arrays.
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: EIGVAL
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: DXSTEP, EIGVEC

      ! . Get the number of modes.
      N = 3 * NFREE

      ! . Allocate some temporary arrays.
      ALLOCATE ( DXSTEP(1:3,1:NATOMS), EIGVAL(1:N), EIGVEC(1:N,1:N) )

      ! . Calculate the energy, gradient and Hessian at the saddle point.
      CALL HESSIAN ( PRINT = .FALSE. )

      ! . Calculate the RMS gradient.
      GRMS = SQRT ( SUM ( ATMDER * ATMDER ) / REAL ( N, DP ) )

      ! . Mass-weight the Hessian if necessary.
      IF ( QMASS ) CALL MASS_WEIGHT_HESSIAN ( ATMHES, PACK ( WEIGHTS, MASK = ( WEIGHTS /= 0.0_DP ) ) )

      ! . Remove the rotational and translational modes.
      CALL RAISE_ROTATION_TRANSLATION ( ATMHES, QMASS )

      ! . Diagonalize the Hessian.
      CALL SYMMETRIC_UPPER ( ATMHES, EIGVAL, EIGVEC )

      ! . Check the number of modes.
      IF ( COUNT ( EIGVAL < 0.0_DP ) /= 1 ) THEN
         CALL PRINT_ERROR ( "TRACE_PATH", "The starting structure is not a saddle point." )
      END IF

      ! . Save the negative mode.
      EIGVAL1 = EIGVAL(1)
      IFREE   = 0
      DO I = 1,NATOMS
         IF ( ATMFIX(I) ) THEN
            DXSTEP(1:3,I) = 0.0_DP
	 ELSE
	    IFREE = IFREE + 1
            II    = 3 * ( IFREE - 1 )
            DXSTEP(1:3,I) = EIGVEC(II+1:II+3,1)
	 END IF
      END DO

      ! . Convert the mode to Cartesian coordinates if necessary.
      IF ( QMASS ) DXSTEP = DXSTEP * WEIGHTS

      ! . Calculate the initial step.
      STPINI = DIRFAC * SQRT ( - 2.0_DP * ELOWER / EIGVAL1 )

      ! . Scale the displacement vector by the initial step size.
      DXSTEP = STPINI * DXSTEP

      ! . Calculate the coordinates of the initial point.
      ATMCRD = ATMCRD + DXSTEP

      ! . Deallocate the temporary space.
      DEALLOCATE ( DXSTEP, EIGVAL, EIGVEC )

      ! . Write out some information about the saddle point.
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#008888", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Saddle Point Information" )
      WRITE ( PRINT_LINE, "(F16.4)" ) ETOTAL ; CALL PRINT_SUMMARY_ELEMENT ( "Potential Energy" )
      WRITE ( PRINT_LINE, "(F16.4)" ) GRMS   ; CALL PRINT_SUMMARY_ELEMENT ( "RMS Gradient"     )
      CALL PRINT_SUMMARY_STOP

      END SUBROUTINE CHECK_SADDLE_POINT

      !-------------------------------
      SUBROUTINE STEEPEST_DESCENT_PATH
      !-------------------------------

      ! . Local scalars.
      INTEGER            :: ISTEP
      REAL ( KIND = DP ) :: GNORM

      ! . Calculate the energy and derivatives at the first point.
      CALL GRADIENT ( PRINT = .FALSE. )

      ! . Write out the header and the data about the first path point.
      IF ( NPRINT > 0 ) THEN
	 CALL PRINT_TABLE_OPTIONS ( COLUMNS = 4, HEADER_COLOR = "#008888", PAGEWIDTH = 50, VARIABLEWIDTHS = (/ 5, 15, 15, 15 /) )
	 CALL PRINT_TABLE_START
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",         HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Distance",     HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Energy",       HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS Gradient", HEADER = .TRUE. )
         WRITE ( PRINT_LINE, "(I5)"    ) 0      ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F15.4)" ) 0.0_DP ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F15.4)" ) ETOTAL ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F15.4)" ) SQRT ( SUM ( ATMDER * ATMDER ) / REAL ( 3 * NFREE, DP ) )
	 CALL PRINT_TABLE_ELEMENT
      END IF

      ! . Write out the initial structure to the trajectory.
      IF ( NSAVE > 0 ) CALL DCD_WRITE ( TRAJECTORY, ATMCRD )

      ! . Loop over the number of steps.
      DO ISTEP = 1,MAXSTP

         ! . Mass weight the gradients if necessary.
         IF ( QMASS ) ATMDER = ATMDER * WEIGHTS

         ! . Calculate the gradient norm at the point.
         GNORM = SQRT ( SUM ( ATMDER * ATMDER ) )

         ! . Normalize the gradient vector.
         ATMDER = ATMDER / GNORM

         ! . Calculate the coordinates of the new point.
         IF ( QMASS ) THEN
            ATMCRD = ATMCRD - DELTA * ATMDER * WEIGHTS
         ELSE
            ATMCRD = ATMCRD - DELTA * ATMDER
         END IF

         ! . Calculate the energy and derivatives at the current point.
         CALL GRADIENT ( PRINT = .FALSE. )

         ! . Write out the data about the path.
         IF ( NPRINT > 0 ) THEN
            IF ( MOD ( ISTEP, NPRINT ) == 0 ) THEN
               WRITE ( PRINT_LINE, "(I5)"    ) ISTEP                      ; CALL PRINT_TABLE_ELEMENT
               WRITE ( PRINT_LINE, "(F15.4)" ) DELTA * REAL ( ISTEP, DP ) ; CALL PRINT_TABLE_ELEMENT
               WRITE ( PRINT_LINE, "(F15.4)" ) ETOTAL                     ; CALL PRINT_TABLE_ELEMENT
               WRITE ( PRINT_LINE, "(F15.4)" ) SQRT ( SUM ( ATMDER * ATMDER ) / REAL ( 3 * NFREE, DP ) )
	       CALL PRINT_TABLE_ELEMENT
            END IF
         END IF

         ! . Write out the structure to the trajectory.
         IF ( NSAVE > 0 ) THEN
            IF ( MOD ( ISTEP, NSAVE ) == 0 ) CALL DCD_WRITE ( TRAJECTORY, ATMCRD )
         END IF

      END DO

      ! . Write out the terminator.
      IF ( NPRINT > 0 ) CALL PRINT_TABLE_STOP

      END SUBROUTINE STEEPEST_DESCENT_PATH

   END SUBROUTINE REACTION_PATH_TRACE

END MODULE REACTION_PATH
