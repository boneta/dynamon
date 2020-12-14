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
!                       The Coordinate Optimization Module
!===============================================================================
!
! . Subroutines:
!
!   OPTIMIZE_BAKER                  Find a minimum or a saddle point.
!   OPTIMIZE_CONJUGATE_GRADIENT     Find a minimum using a CG method.
!
!===============================================================================
MODULE OPTIMIZE_COORDINATES

! . Utility data structure declarations.
USE DEFINITIONS,           ONLY : DP
USE PRINTING,              ONLY : PRINT_ERROR

USE BAKER_OPTIMIZATION,    ONLY : BAKER_SEARCH
USE CONJUGATE_GRADIENT,    ONLY : CONJUGATE_GRADIENT_MINIMIZE

USE ATOMS,                 ONLY : ATMCRD, ATMFIX, NATOMS, NFREE
USE NORMAL_MODE_UTILITIES, ONLY : RAISE_ROTATION_TRANSLATION
USE POTENTIAL_ENERGY,      ONLY : ATMDER, ATMHES, ETOTAL, GRADIENT, HESSIAN

IMPLICIT NONE
PRIVATE
PUBLIC :: OPTIMIZE_BAKER, OPTIMIZE_CONJUGATE_GRADIENT

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------------------------------------------------------------
   SUBROUTINE OPTIMIZE_BAKER ( POINT, STATUS, FOLLOW_MODE, FOLLOW_VARIABLE, PRINT_FREQUENCY, &
                                              STEP_NUMBER, USE_NR_STEP, GRADIENT_TOLERANCE,  &
                                        MAXIMUM_EIGENVALUE, MAXIMUM_STEP, MINIMUM_EIGENVALUE )
   !------------------------------------------------------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: POINT

   ! . Optional arguments.
   INTEGER,            INTENT(IN),  OPTIONAL :: FOLLOW_MODE, FOLLOW_VARIABLE, PRINT_FREQUENCY, STEP_NUMBER
   INTEGER,            INTENT(OUT), OPTIONAL :: STATUS
   LOGICAL,            INTENT(IN),  OPTIONAL :: USE_NR_STEP
   REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: GRADIENT_TOLERANCE, MAXIMUM_EIGENVALUE, MAXIMUM_STEP, MINIMUM_EIGENVALUE

   ! . Local scalars.
   LOGICAL :: QSADDLE

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: X

   ! . Initialize the STATUS counter if necessary.
   IF ( PRESENT ( STATUS ) ) STATUS = 0

   ! . Check the number of free atoms.
   IF ( NFREE < 0 ) RETURN

   ! . Determine the type of point to optimize.
   IF ( POINT == "MINIMUM" ) THEN
      QSADDLE = .FALSE.
   ELSE IF ( POINT == "SADDLE" ) THEN
      QSADDLE = .TRUE.
   ELSE
      CALL PRINT_ERROR ( "OPTIMIZE_BAKER", "Unrecognized stationary point type." )
   END IF

   ! . Allocate space for the variables.
   ALLOCATE ( X(1:3*NFREE) )

   ! . Move the coordinates to X.
   CALL VARIABLES_FILL ( X, ATMCRD )

   ! . Perform the minimization.
   CALL BAKER_SEARCH ( EGHCALC, X, STATUS, FOLLOW_MODE, FOLLOW_VARIABLE, PRINT_FREQUENCY, STEP_NUMBER, &
                                   QSADDLE, USE_NR_STEP, GRADIENT_TOLERANCE, MAXIMUM_EIGENVALUE,       &
                                   MAXIMUM_STEP, MINIMUM_EIGENVALUE )

   ! . Deallocate the temporary array.
   CALL VARIABLES_EMPTY ( X, ATMCRD )

   ! . Deallocate the temporary array.
   DEALLOCATE ( X )

   END SUBROUTINE OPTIMIZE_BAKER

   !-------------------------------------------------------------------------------------------------------------
   SUBROUTINE OPTIMIZE_CONJUGATE_GRADIENT ( STATUS, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, GRADIENT_TOLERANCE )
   !-------------------------------------------------------------------------------------------------------------

   ! . Optional scalar arguments.
   INTEGER,            INTENT(IN),  OPTIONAL :: PRINT_FREQUENCY, STEP_NUMBER
   INTEGER,            INTENT(OUT), OPTIONAL :: STATUS
   REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: GRADIENT_TOLERANCE, STEP_SIZE

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: X

   ! . Initialize the STATUS counter if necessary.
   IF ( PRESENT ( STATUS ) ) STATUS = 0

   ! . Check the number of free atoms.
   IF ( NFREE < 0 ) RETURN

   ! . Allocate space for the variables.
   ALLOCATE ( X(1:3*NFREE) )

   ! . Move the coordinates to X.
   CALL VARIABLES_FILL ( X, ATMCRD )

   ! . Perform the minimization.
   CALL CONJUGATE_GRADIENT_MINIMIZE ( EGCALC, X, STATUS, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, GRADIENT_TOLERANCE )

   ! . Save the optimized coordinates.
   CALL VARIABLES_EMPTY ( X, ATMCRD )

   ! . Deallocate the temporary array.
   DEALLOCATE ( X )

   END SUBROUTINE OPTIMIZE_CONJUGATE_GRADIENT

   !---------------------------------------
   SUBROUTINE VARIABLES_EMPTY ( X, ATMDAT )
   !---------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(IN)  :: X
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: ATMDAT

   ! . Local scalars.
   INTEGER :: IATOM, II

   ! . Loop over the free atoms.
   II = -3
   DO IATOM = 1,NATOMS
      IF ( .NOT. ATMFIX(IATOM) ) THEN
         II = II + 3
         ATMDAT(1:3,IATOM) = X(II+1:II+3)
      END IF
   END DO

   END SUBROUTINE VARIABLES_EMPTY

   !--------------------------------------
   SUBROUTINE VARIABLES_FILL ( X, ATMDAT )
   !--------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT) :: X
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)  :: ATMDAT

   ! . Local scalars.
   INTEGER :: IATOM, II

   ! . Loop over the free atoms.
   II = -3
   DO IATOM = 1,NATOMS
      IF ( .NOT. ATMFIX(IATOM) ) THEN
         II = II + 3
         X(II+1:II+3) = ATMDAT(1:3,IATOM)
      END IF
   END DO

   END SUBROUTINE VARIABLES_FILL

   !----------------------------
   SUBROUTINE EGCALC ( X, E, G )
   !----------------------------

   ! . Argument declarations.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)  :: X
   REAL ( KIND = DP ),               INTENT(OUT) :: E
   REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: G

   ! . Move X to the coordinates.
   CALL VARIABLES_EMPTY ( X, ATMCRD )

   ! . Calculate the energy and its first derivatives.
   CALL GRADIENT ( PRINT = .FALSE. ) ; E = ETOTAL ; CALL VARIABLES_FILL ( G, ATMDER )

   END SUBROUTINE EGCALC

   !--------------------------------
   SUBROUTINE EGHCALC ( X, E, G, H )
   !--------------------------------

   ! . Argument declarations.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(IN)  :: X
   REAL ( KIND = DP ),                 INTENT(OUT) :: E
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT) :: G, H

   ! . Move X to the coordinates.
   CALL VARIABLES_EMPTY ( X, ATMCRD )

   ! . Calculate the energy and its derivatives.
   CALL HESSIAN ( PRINT = .FALSE. ) ; E = ETOTAL ; CALL VARIABLES_FILL ( G, ATMDER ) ; H = ATMHES
! --------------------------------------------------------------------------------------------------
!write( 99 ) etotal
!write( 99 ) x
!write( 99 ) g
!write( 99 ) h
! --------------------------------------------------------------------------------------------------
   ! . Change the eigenvalues of the rotational and translational modes.
   CALL RAISE_ROTATION_TRANSLATION ( H, .FALSE. )

   END SUBROUTINE EGHCALC

END MODULE OPTIMIZE_COORDINATES
