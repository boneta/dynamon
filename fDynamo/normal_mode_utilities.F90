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
!                         The Normal Mode Utility Module
!===============================================================================
!
! . Subroutines:
!
!   GET_ROTATION_TRANSLATION        Get the rotation and translation modes.
!   MASS_WEIGHT_HESSIAN             Mass weight a Hessian matrix.
!   PROJECT_ROTATION_TRANSLATION    Project out the zero modes eigenvectors.
!   RAISE_ROTATION_TRANSLATION      Increase the eigenvalues of the zero modes.
!
!===============================================================================
MODULE NORMAL_MODE_UTILITIES

! . Module declarations.
USE DEFINITIONS,    ONLY : DP
USE ATOMS,          ONLY : ATMCRD, ATMMAS, NATOMS, ATMFIX, NFREE

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !--------------------------------------------------------------
   SUBROUTINE GET_ROTATION_TRANSLATION ( NTRROT, TRROT, QWEIGHTS )
   !--------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(OUT) :: NTRROT
   LOGICAL, INTENT(IN)  :: QWEIGHTS

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: TRROT

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: TOL = 1.0E-10_DP

   ! . Local scalars.
   INTEGER            :: I, II, J, N
   REAL ( KIND = DP ) :: MI, SUM

   ! . Other local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: RCG

   ! . Determine the matrix dimension.
   N = 3 * NFREE

   ! . Initialize the modes.
   TRROT = 0.0_DP

   ! . Check the mode dimensions.
   IF ( ( SIZE ( TRROT, 1 ) < N ) .OR. ( SIZE ( TRROT, 2 ) < 6 ) ) RETURN

   ! . Mass weighting is required.
   IF ( QWEIGHTS ) THEN

      ! . Calculate the center of mass of the system.
      SUM = .0_DP
      RCG = .0_DP
      DO I = 1, NATOMS
         IF( .NOT. ATMFIX(I) ) THEN
             SUM = SUM + ATMMAS(I)
             RCG = RCG + ATMCRD(1:3,I) * ATMMAS(I)
         END IF
      END DO
      RCG = RCG / SUM

      ! . Calculate the translational and rotational modes.
      J = 1
      DO I = 1,NATOMS
         IF( .NOT. ATMFIX(I) ) THEN
             II = 3 * (J-1)
             MI = SQRT ( ATMMAS(I) )
             TRROT(II+1,1) = MI
             TRROT(II+2,2) = MI
             TRROT(II+3,3) = MI
             TRROT(II+3,4) =   MI * ( ATMCRD(2,I) - RCG(2) )
             TRROT(II+2,4) = - MI * ( ATMCRD(3,I) - RCG(3) )
             TRROT(II+1,5) =   MI * ( ATMCRD(3,I) - RCG(3) )
             TRROT(II+3,5) = - MI * ( ATMCRD(1,I) - RCG(1) )
             TRROT(II+2,6) =   MI * ( ATMCRD(1,I) - RCG(1) )
             TRROT(II+1,6) = - MI * ( ATMCRD(2,I) - RCG(2) )
             J = J + 1
         END IF
      END DO

   ! . No mass weighting is needed.
   ELSE

      ! . Calculate the center of geometry of the system.
      RCG = .0_DP
      DO I = 1, NATOMS
         IF( .NOT. ATMFIX(I) ) THEN
             RCG = RCG + ATMCRD(1:3,I)
         END IF
      END DO
      RCG = RCG / REAL( NFREE, DP )

      ! . Calculate the translational modes.
      J = 1
      DO I = 1,NATOMS
         IF( .NOT. ATMFIX(I) ) THEN
            II = 3 * (J-1)
            TRROT(II+1,1) = 1.0_DP
            TRROT(II+2,2) = 1.0_DP
            TRROT(II+3,3) = 1.0_DP
            TRROT(II+3,4) =   ( ATMCRD(2,I) - RCG(2) )
            TRROT(II+2,4) = - ( ATMCRD(3,I) - RCG(3) )
            TRROT(II+1,5) =   ( ATMCRD(3,I) - RCG(3) )
            TRROT(II+3,5) = - ( ATMCRD(1,I) - RCG(1) )
            TRROT(II+2,6) =   ( ATMCRD(1,I) - RCG(1) )
            TRROT(II+1,6) = - ( ATMCRD(2,I) - RCG(2) )
             J = J + 1
         END IF
      END DO

   END IF

   ! . Create a linearly independent set of modes.
   NTRROT = 0
   DO I = 1,6
      DO J = 1,NTRROT
         TRROT(1:N,I) = TRROT(1:N,I) - DOT_PRODUCT ( TRROT(1:N,I), TRROT(1:N,J) ) * TRROT(1:N,J)
      END DO
      SUM = DOT_PRODUCT ( TRROT(1:N,I), TRROT(1:N,I) )
      IF ( SUM > TOL ) THEN
         NTRROT = NTRROT + 1
         TRROT(1:N,NTRROT) = TRROT(1:N,I) / SQRT ( SUM )
      END IF
   END DO

   END SUBROUTINE GET_ROTATION_TRANSLATION

   !--------------------------------------------------
   SUBROUTINE MASS_WEIGHT_HESSIAN ( HESSIAN, WEIGHTS )
   !--------------------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)    :: WEIGHTS
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: HESSIAN

   ! . Local scalars.
   INTEGER :: I, II, J, N
   
   ! . Get the dimension of the problem.
   N = SIZE ( WEIGHTS )

   ! . Loop over the elements in the matrix.
   II = 0
   DO I = 1,N
      DO J = 1,I
         II = II + 1
         HESSIAN(II) = HESSIAN(II) * WEIGHTS(I) * WEIGHTS(J)
      END DO
   END DO

   END SUBROUTINE MASS_WEIGHT_HESSIAN

   !------------------------------------------------------
   SUBROUTINE PROJECT_ROTATION_TRANSLATION ( H, QWEIGHTS )
   !------------------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN) :: QWEIGHTS

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: H

   ! . Local scalars.
   INTEGER            :: I, II, ITR, J, JTR, N, NTRROT
   REAL ( KIND = DP ) :: SUM

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: HV, TRROT

   ! . Other local arrays.
   REAL ( KIND = DP ), DIMENSION(1:6,1:6) :: VHV

   ! . Determine the matrix dimension.
   N = 3 * NFREE

   ! . Allocate space for the modes.
   ALLOCATE ( HV(1:N,1:6), TRROT(1:N,1:6) )

   ! . Calculate the rotation and translation modes.
   CALL GET_ROTATION_TRANSLATION ( NTRROT, TRROT, QWEIGHTS )

   ! . Remove the modes from the second derivative matrix.
   IF ( NTRROT > 0 ) THEN
      
      ! . Loop over the modes.
      DO ITR = 1,NTRROT

         ! . Calculate the dot product of the Hessian with the mode.
         DO I = 1,N

            II = ( I * ( I - 1 ) ) / 2
            SUM = 0.0_DP
            DO J = 1,I
               II = II + 1
               SUM = SUM + H(II) * TRROT(J,ITR)
            END DO
            DO J = (I+1),N
               II = II + J - 1
               SUM = SUM + H(II) * TRROT(J,ITR)
            END DO

            HV(I,ITR) = SUM

         END DO

         ! . Calculate VHV.
         DO JTR = 1,ITR
            VHV(JTR,ITR) = DOT_PRODUCT ( TRROT(1:N,ITR), HV(1:N,JTR) )
            VHV(ITR,JTR) = VHV(JTR,ITR)
         END DO

      END DO

      ! . Modify the Hessian matrix.
      II = 0
      DO I = 1,N
         DO J = 1,I

            ! . Calculate the projection terms.
            SUM = 0.0_DP
            DO ITR = 1,NTRROT
               SUM = SUM - HV(I,ITR) * TRROT(J,ITR) - HV(J,ITR) * TRROT(I,ITR)
               DO JTR = 1,NTRROT
                  SUM = SUM + TRROT(J,JTR) * VHV(JTR,ITR) * TRROT(I,ITR)
               END DO
            END DO

            ! . Modify the matrix.
            II = II + 1
            H(II) = H(II) + SUM

         END DO
      END DO

   END IF

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( HV, TRROT )

   END SUBROUTINE PROJECT_ROTATION_TRANSLATION

   !----------------------------------------------------
   SUBROUTINE RAISE_ROTATION_TRANSLATION ( H, QWEIGHTS )
   !----------------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN) :: QWEIGHTS

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: H

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: LARGE = 5.0E+4_DP

   ! . Local scalars.
   INTEGER :: I, II, ITR, J, N, NTRROT

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: TRROT

   ! . Determine the matrix dimension.
   N = 3 * NFREE

   ! . Allocate space for the modes.
   ALLOCATE ( TRROT(1:N,1:6) )

   ! . Calculate the rotation and translation modes.
   CALL GET_ROTATION_TRANSLATION ( NTRROT, TRROT, QWEIGHTS )

   ! . Remove the modes from the second derivative matrix.
   IF ( NTRROT > 0 ) THEN
      II = 0
      DO I = 1,N
         DO J = 1,I
            II = II + 1
            DO ITR = 1,NTRROT
               H(II) = H(II) + LARGE * TRROT(I,ITR) * TRROT(J,ITR)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( TRROT )

   END SUBROUTINE RAISE_ROTATION_TRANSLATION

END MODULE NORMAL_MODE_UTILITIES
