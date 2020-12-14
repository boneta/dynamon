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
!                           The Linear Algebra Module
!===============================================================================
!
! . Functions:
!
!   CROSS_PRODUCT                Calculate the cross-product of two 3-vectors.
!   DETERMINANT                  Calculate the determinant of a 3x3 matrix.
!   NORM                         The norm of a vector.
!   NORMALIZE                    Normalize a vector.
!   ROTATION_MATRIX_AXIS         Construct a rotation matrix from an axis.
!   ROTATION_MATRIX_EULER        Construct a rotation matrix from Euler angles.
!   ROTATION_MATRIX_XYZ          Construct a rotation matrix from XYZ angles.
!
! . Subroutines:
!
!   PROJECT_OUT                  Make a vector orthogonal to a set of other ones.
!   SCHMIDT_ORTHOGONALIZE        Orthogonalize a set of vectors.  
!
!===============================================================================
MODULE LINEAR_ALGEBRA

! . Module declarations.
USE DEFINITIONS, ONLY : DP

IMPLICIT NONE
PUBLIC
PRIVATE :: NORMALIZATION_TOLERANCE

! . Module parameters.
REAL ( KIND = DP ), PARAMETER :: NORMALIZATION_TOLERANCE = 1.0E-16_DP

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------
   FUNCTION CROSS_PRODUCT ( A, B )
   !------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3) :: CROSS_PRODUCT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN) :: A, B

   ! . Calculate the cross product A x B.
   CROSS_PRODUCT(1) = A(2) * B(3) - A(3) * B(2)
   CROSS_PRODUCT(2) = A(3) * B(1) - A(1) * B(3)
   CROSS_PRODUCT(3) = A(1) * B(2) - A(2) * B(1)

   END FUNCTION CROSS_PRODUCT

   !-------------------------
   FUNCTION DETERMINANT ( A )
   !-------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: DETERMINANT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:3), INTENT(IN) :: A

   ! . Calculate the determinant of A.
   DETERMINANT = A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) - &
                 A(1,2) * ( A(2,1) * A(3,3) - A(3,1) * A(2,3) ) + &
                 A(1,3) * ( A(2,1) * A(3,2) - A(3,1) * A(2,2) )

   END FUNCTION DETERMINANT

   !-----------------------
   FUNCTION NORM ( VECTOR )
   !-----------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: VECTOR

   ! . Function declarations.
   REAL ( KIND = DP ) :: NORM

   ! . Local scalars.
   INTEGER :: N

   ! . Initialization.
   NORM = 0.0_DP

   ! . Determine and check the size of the vector.
   N = SIZE ( VECTOR )
   IF ( N <= 0 ) RETURN

   ! . Get the size of the vector.
   NORM = SQRT ( DOT_PRODUCT ( VECTOR(1:N), VECTOR(1:N) ) )

   ! . If NORM is too small make it zero.
   IF ( NORM < NORMALIZATION_TOLERANCE ) NORM = 0.0_DP

   END FUNCTION NORM

   !----------------------------
   FUNCTION NORMALIZE ( VECTOR )
   !----------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: VECTOR

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:SIZE(VECTOR)) :: NORMALIZE

   ! . Local scalars.
   INTEGER            :: N
   REAL ( KIND = DP ) :: FACT

   ! . Determine and check the size of the vector.
   N = SIZE ( VECTOR )
   IF ( N <= 0 ) RETURN

   ! . Get the NORM of the vector.
   FACT = NORM ( VECTOR )

   ! . Normalize the vector.
   IF ( FACT > 0.0_DP ) NORMALIZE = VECTOR / FACT

   END FUNCTION NORMALIZE

   !--------------------------------------------
   FUNCTION ROTATION_MATRIX_AXIS ( ANGLE, AXIS )
   !--------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: ROTATION_MATRIX_AXIS

   ! . Arguments.
   REAL ( KIND = DP ),                 INTENT(IN) :: ANGLE
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN) :: AXIS

   ! . Local scalars.
   REAL ( KIND = DP ) :: E0, E1, E2, E3

   ! . Get the matrix.
   E0 =           COS ( ANGLE / 2.0_DP )
   E1 = AXIS(1) * SIN ( ANGLE / 2.0_DP ) 
   E2 = AXIS(2) * SIN ( ANGLE / 2.0_DP )
   E3 = AXIS(3) * SIN ( ANGLE / 2.0_DP )
   ROTATION_MATRIX_AXIS(1,1) = E0*E0 + E1*E1 - E2*E2 - E3*E3
   ROTATION_MATRIX_AXIS(1,2) = 2.0_DP * ( E1*E2 + E0*E3 )
   ROTATION_MATRIX_AXIS(1,3) = 2.0_DP * ( E1*E3 - E0*E2 )
   ROTATION_MATRIX_AXIS(2,1) = 2.0_DP * ( E1*E2 - E0*E3 )
   ROTATION_MATRIX_AXIS(2,2) = E0*E0 - E1*E1 + E2*E2 - E3*E3
   ROTATION_MATRIX_AXIS(2,3) = 2.0_DP * ( E2*E3 + E0*E1 )
   ROTATION_MATRIX_AXIS(3,1) = 2.0_DP * ( E1*E3 + E0*E2 )
   ROTATION_MATRIX_AXIS(3,2) = 2.0_DP * ( E2*E3 - E0*E1 )
   ROTATION_MATRIX_AXIS(3,3) = E0*E0 - E1*E1 - E2*E2 + E3*E3

   END FUNCTION ROTATION_MATRIX_AXIS

   !----------------------------------------
   FUNCTION ROTATION_MATRIX_EULER ( ANGLES )
   !----------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: ROTATION_MATRIX_EULER

   ! . Arguments.
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN) :: ANGLES

   ! . Local scalars.
   REAL ( KIND = DP ) :: E0, E1, E2, E3, PHI, PSI, THETA

   ! . Get the half angles.
   PHI = 0.5_DP * ANGLES(1) ; PSI = 0.5_DP * ANGLES(2) ; THETA = 0.5_DP * ANGLES(3)

   ! . Get the matrix.
   E0 = COS ( PHI + PSI ) * COS ( THETA )
   E1 = COS ( PHI - PSI ) * SIN ( THETA )
   E2 = SIN ( PHI - PSI ) * SIN ( THETA )
   E3 = SIN ( PHI + PSI ) * COS ( THETA )
   ROTATION_MATRIX_EULER(1,1) = E0*E0 + E1*E1 - E2*E2 - E3*E3
   ROTATION_MATRIX_EULER(1,2) = 2.0_DP * ( E1*E2 + E0*E3 )
   ROTATION_MATRIX_EULER(1,3) = 2.0_DP * ( E1*E3 - E0*E2 )
   ROTATION_MATRIX_EULER(2,1) = 2.0_DP * ( E1*E2 - E0*E3 )
   ROTATION_MATRIX_EULER(2,2) = E0*E0 - E1*E1 + E2*E2 - E3*E3
   ROTATION_MATRIX_EULER(2,3) = 2.0_DP * ( E2*E3 + E0*E1 )
   ROTATION_MATRIX_EULER(3,1) = 2.0_DP * ( E1*E3 + E0*E2 )
   ROTATION_MATRIX_EULER(3,2) = 2.0_DP * ( E2*E3 - E0*E1 )
   ROTATION_MATRIX_EULER(3,3) = E0*E0 - E1*E1 - E2*E2 + E3*E3

   END FUNCTION ROTATION_MATRIX_EULER

   !--------------------------------------
   FUNCTION ROTATION_MATRIX_XYZ ( ANGLES )
   !--------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: ROTATION_MATRIX_XYZ

   ! . Arguments.
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN) :: ANGLES

   ! . Local scalars.
   INTEGER            :: I, IANGLE
   REAL ( KIND = DP ) :: COSA, SINA

   ! . Local arrays.
   INTEGER,            DIMENSION(1:2)     :: COMPONENTS
   REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: R

   ! . Initialization.
   ROTATION_MATRIX_XYZ = 0.0_DP
   DO I = 1,3
      ROTATION_MATRIX_XYZ(I,I) = 1.0_DP
   END DO

   ! . Loop over the matices.
   DO IANGLE = 1,3

      ! . Fill components.
      SELECT CASE ( IANGLE )
      CASE ( 1 ) ; COMPONENTS(1) = 2 ; COMPONENTS(2) = 3
      CASE ( 2 ) ; COMPONENTS(1) = 3 ; COMPONENTS(2) = 1
      CASE ( 3 ) ; COMPONENTS(1) = 1 ; COMPONENTS(2) = 2
      END SELECT

      ! . Get COSA and SINA.
      COSA = COS ( ANGLES(IANGLE) ) ; SINA = SIN ( ANGLES(IANGLE) )

      ! . Initialize the local matrix.
      R = 0.0_DP

      ! . Diagonal terms.
      R(IANGLE,IANGLE)               =   1.0_DP
      R(COMPONENTS(1),COMPONENTS(1)) =   COSA
      R(COMPONENTS(2),COMPONENTS(2)) =   COSA

      ! . Off-diagonal terms.
      R(COMPONENTS(1),COMPONENTS(2)) =   SINA
      R(COMPONENTS(2),COMPONENTS(1)) = - SINA

      ! . Sum in the total matrix.
      ROTATION_MATRIX_XYZ = MATMUL ( R, ROTATION_MATRIX_XYZ )

   END DO

   END FUNCTION ROTATION_MATRIX_XYZ

   !-----------------------------------------
   SUBROUTINE PROJECT_OUT ( VECTOR, VECTORS )
   !-----------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(INOUT) :: VECTOR
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)    :: VECTORS

   ! . Scalar arguments.
   INTEGER            :: IVEC, N, NVEC
   REAL ( KIND = DP ) :: DOTFAC

   ! . Get the dimension of the vectors.
   N    = SIZE ( VECTOR,  1 )
   NVEC = SIZE ( VECTORS, 2 )

   ! . Return if there is dimension mismatch or the second set of vectors is null.
   IF ( ( N /= SIZE ( VECTORS, 1 ) ) .OR. ( NVEC <= 0 ) ) RETURN

   ! . Project out VECTORS from VECTOR.
   DO IVEC = 1,NVEC
      DOTFAC = DOT_PRODUCT ( VECTOR(1:N), VECTORS(1:N,IVEC) )
      VECTOR(1:N) = VECTOR(1:N) - DOTFAC * VECTORS(1:N,IVEC)
   END DO

   END SUBROUTINE PROJECT_OUT

   !---------------------------------------------------------
   SUBROUTINE SCHMIDT_ORTHOGONALIZE ( VECTORS, NINDEPENDENT )
   !---------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(OUT) :: NINDEPENDENT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: VECTORS

   ! . Scalar arguments.
   INTEGER            :: IVEC, JVEC, N, NVEC
   REAL ( KIND = DP ) :: DOTFAC

   ! . Initialization.
   NINDEPENDENT = 0

   ! . Get the dimensions of the vectors.
   N    = SIZE ( VECTORS, 1 )
   NVEC = SIZE ( VECTORS, 2 )

   ! . Return if there is any null dimension.
   IF ( ( N <= 0 ) .OR. ( NVEC <= 0 ) ) RETURN

   ! . Loop over the vectors.
   DO IVEC = 1,NVEC

      ! . Take out the contributions from vectors of lower index.
      DO JVEC = 1,NINDEPENDENT
         DOTFAC = DOT_PRODUCT ( VECTORS(1:N,IVEC), VECTORS(1:N,JVEC) )
         VECTORS(1:N,IVEC) = VECTORS(1:N,IVEC) - DOTFAC * VECTORS(1:N,JVEC)
      END DO

      ! . Get the norm of the vector.
      DOTFAC = NORM ( VECTORS(1:N,IVEC) )

      ! . Add the vector to the set.
      IF ( DOTFAC > 0.0_DP ) THEN

         ! . Increment the number of independent vectors.
         NINDEPENDENT = NINDEPENDENT + 1

         ! . Store the vector.
         VECTORS(1:N,NINDEPENDENT) = VECTORS(1:N,IVEC) / DOTFAC

      END IF

   END DO

   ! . Zero out the remaining vector components.
   VECTORS(1:N,NINDEPENDENT+1:NVEC) = 0.0_DP

   END SUBROUTINE SCHMIDT_ORTHOGONALIZE

END MODULE LINEAR_ALGEBRA

