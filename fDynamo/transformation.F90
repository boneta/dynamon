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
!                      The Coordinate Transformation Module
!===============================================================================
!
! . Functions:
!
!   CENTER                         Calculate the geometric center of a set of
!                                  coordinates (with or without weighting).
!   INERTIA_MATRIX                 Calculate the inertia_matrix.
!
! . Subroutines:
!
!   MOMENTS_OF_INERTIA             Calculate the moments and axes of inertia.
!   MOMENTS_OF_INERTIA_AT_CENTER   Calculate the moments of inertia at a center.
!   ROTATE                         Rotate a set of coordinates.
!   TRANSLATE                      Translate a set of coordinates.
!   TRANSLATE_TO_CENTER            Translate a set of coordinates to its center.
!   TO_PRINCIPAL_AXES              Transform a set of coordinates to its
!                                  principal axes.
!
! . Notes:
!
!   The procedures in this module perform a number of useful transformations
!   on a set of coordinates. Most of the procedures take as arguments an
!   array of coordinates and an optional array of weights which will often
!   contain the atom masses. The weights are used in the calculation of the
!   center of the coordinate set or in the calculation of the inertia matrix.
!
!===============================================================================
MODULE TRANSFORMATION

! . Module declarations.
USE DEFINITIONS,     ONLY : DP
USE DIAGONALIZATION, ONLY : SYMMETRIC_UPPER
USE LINEAR_ALGEBRA,  ONLY : DETERMINANT
USE PRINTING,        ONLY : PRINT_PARAGRAPH

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !---------------------------------------
   FUNCTION CENTER ( COORDINATES, WEIGHTS )
   !---------------------------------------

   ! . Function declaration.
   REAL ( KIND = DP ), DIMENSION(1:3) :: CENTER

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   INTEGER :: I, NATOMS

   ! . Initialize CENTER.
   CENTER = 0.0_DP

   ! . Determine the number of atoms.
   NATOMS = SIZE ( COORDINATES, 2 )

   ! . Check the dimensions of COORDINATES.
   IF ( ( SIZE ( COORDINATES, 1 ) /= 3 ) .OR. ( NATOMS <  1 ) ) RETURN

   ! . The WEIGHTS array is present.
   IF ( PRESENT ( WEIGHTS ) ) THEN

      ! . Check the dimension of the WEIGHTS array.
      IF ( SIZE ( WEIGHTS ) /= NATOMS ) RETURN

      ! . Calculate the center.
      DO I = 1,NATOMS
         CENTER(1:3) = CENTER(1:3) + WEIGHTS(I) * COORDINATES(1:3,I)
      END DO

      ! . Normalize the center.
      CENTER(1:3) = CENTER(1:3) / SUM ( WEIGHTS(1:NATOMS) )

   ! . The WEIGHTS array is absent.
   ELSE

      ! . Calculate the center.
      DO I = 1,NATOMS
         CENTER(1:3) = CENTER(1:3) + COORDINATES(1:3,I)
      END DO

      ! . Normalize the center.
      CENTER(1:3) = CENTER(1:3) / REAL ( NATOMS, DP )

   END IF

   END FUNCTION CENTER

   !-----------------------------------------------
   FUNCTION INERTIA_MATRIX ( COORDINATES, WEIGHTS )
   !-----------------------------------------------

   ! . Function declaration.
   REAL ( KIND = DP ), DIMENSION(1:6) :: INERTIA_MATRIX

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   INTEGER            :: I, NATOMS
   REAL ( KIND = DP ) :: W, X, Y, Z, XX, XY, XZ, YY, YZ, ZZ

   ! . Determine the number of atoms.
   NATOMS = SIZE ( COORDINATES, 2 )

   ! . Check the dimensions of COORDINATES.
   IF ( ( SIZE ( COORDINATES, 1 ) /= 3 ) .OR. ( NATOMS <  1 ) ) RETURN

   ! . Initialize the elements of the inertia matrix.
   XX = 0.0_DP
   XY = 0.0_DP
   XZ = 0.0_DP
   YY = 0.0_DP
   YZ = 0.0_DP
   ZZ = 0.0_DP

   ! . The WEIGHTS array is present.
   IF ( PRESENT ( WEIGHTS ) ) THEN

      ! . Check the dimension of the WEIGHTS array.
      IF ( SIZE ( WEIGHTS ) /= NATOMS ) RETURN

      ! . Calculate the inertia matrix.
      DO I = 1,NATOMS
         W = WEIGHTS(I)
         X = COORDINATES(1,I)
         Y = COORDINATES(2,I)
         Z = COORDINATES(3,I)
         XX = XX + W * X * X
         XY = XY + W * X * Y
         XZ = XZ + W * X * Z
         YY = YY + W * Y * Y
         YZ = YZ + W * Y * Z
         ZZ = ZZ + W * Z * Z
      END DO

   ! . The WEIGHTS array is absent.
   ELSE

      ! . Calculate the inertia matrix.
      DO I = 1,NATOMS
         X = COORDINATES(1,I)
         Y = COORDINATES(2,I)
         Z = COORDINATES(3,I)
         XX = XX + X * X
         XY = XY + X * Y
         XZ = XZ + X * Z
         YY = YY + Y * Y
         YZ = YZ + Y * Z
         ZZ = ZZ + Z * Z         
      END DO

   END IF

   ! . Set the elements of the inertia matrix (upper triangle format).
   INERTIA_MATRIX(1) = YY + ZZ ! XX
   INERTIA_MATRIX(2) = - XY    ! YX
   INERTIA_MATRIX(3) = XX + ZZ ! YY
   INERTIA_MATRIX(4) = - XZ    ! ZX
   INERTIA_MATRIX(5) = - YZ    ! ZY
   INERTIA_MATRIX(6) = XX + YY ! ZZ

   END FUNCTION INERTIA_MATRIX

   !--------------------------------------------------------------------
   SUBROUTINE MOMENTS_OF_INERTIA ( COORDINATES, MOMENTS, AXES, WEIGHTS )
   !--------------------------------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)  :: COORDINATES
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(OUT) :: MOMENTS

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:3), INTENT(OUT), OPTIONAL :: AXES
   REAL ( KIND = DP ), DIMENSION(:),       INTENT(IN),  OPTIONAL :: WEIGHTS

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: TOL = 1.0E-6_DP

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:6) :: INERTIA

   ! . Initialize the moments array.
   MOMENTS = HUGE ( 0.0_DP )

   ! . Calculate the inertia matrix.
   INERTIA = INERTIA_MATRIX ( COORDINATES, WEIGHTS )

   ! . Check the matrix.
   IF ( MAXVAL ( ABS ( INERTIA ) ) < TOL ) RETURN

   ! . The axes are required.
   IF ( PRESENT ( AXES ) ) THEN

      ! . Diagonalize the inertia matrix.
      CALL SYMMETRIC_UPPER ( INERTIA, EIGENVALUES = MOMENTS, EIGENVECTORS = AXES )

      ! . Change the sign of the first eigenvector so that the determinant is positive.
      IF ( DETERMINANT ( AXES ) < 0.0_DP ) AXES(1:3,1) = - AXES(1:3,1)

   ! . The axes are not required.
   ELSE

      ! . Diagonalize the inertia matrix.
      CALL SYMMETRIC_UPPER ( INERTIA, EIGENVALUES = MOMENTS )

   END IF

   ! . Zero out the small moments of inertia.
   WHERE ( ABS ( MOMENTS ) < TOL ) MOMENTS = 0.0_DP

   END SUBROUTINE MOMENTS_OF_INERTIA

   !------------------------------------------------------------------------------
   SUBROUTINE MOMENTS_OF_INERTIA_AT_CENTER ( COORDINATES, MOMENTS, AXES, WEIGHTS )
   !------------------------------------------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)  :: COORDINATES
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(OUT) :: MOMENTS

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:3), INTENT(OUT), OPTIONAL :: AXES
   REAL ( KIND = DP ), DIMENSION(:),       INTENT(IN),  OPTIONAL :: WEIGHTS

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:SIZE(COORDINATES,1),1:SIZE(COORDINATES,2)) :: XYZ

   ! . Save the input coordinates in XYZ.
   XYZ = COORDINATES

   ! . Move the system to its center.
   CALL TRANSLATE_TO_CENTER ( XYZ, WEIGHTS )

   ! . Calculate the moments of inertia.
   CALL MOMENTS_OF_INERTIA ( XYZ, MOMENTS, AXES, WEIGHTS )

   END SUBROUTINE MOMENTS_OF_INERTIA_AT_CENTER

   !------------------------------------------
   SUBROUTINE ROTATE ( COORDINATES, ROTATION )
   !------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:),     INTENT(INOUT) :: COORDINATES
   REAL ( KIND = DP ), DIMENSION(1:3,1:3), INTENT(IN)    :: ROTATION

   ! . Local scalars.
   INTEGER            :: I, NATOMS
   REAL ( KIND = DP ) :: X, Y, Z

   ! . Determine the number of atoms.
   NATOMS = SIZE ( COORDINATES, 2 )

   ! . Check the dimensions of COORDINATES.
   IF ( ( SIZE ( COORDINATES, 1 ) /= 3 ) .OR. ( NATOMS <  1 ) ) RETURN

   ! . Rotate the coordinates.
   DO I = 1,NATOMS
      X = COORDINATES(1,I)
      Y = COORDINATES(2,I)
      Z = COORDINATES(3,I)
      COORDINATES(1,I) = ROTATION(1,1)*X + ROTATION(1,2)*Y + ROTATION(1,3)*Z
      COORDINATES(2,I) = ROTATION(2,1)*X + ROTATION(2,2)*Y + ROTATION(2,3)*Z
      COORDINATES(3,I) = ROTATION(3,1)*X + ROTATION(3,2)*Y + ROTATION(3,3)*Z
   END DO

   END SUBROUTINE ROTATE

   !------------------------------------------------
   SUBROUTINE TRANSLATE ( COORDINATES, TRANSLATION )
   !------------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: COORDINATES
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN)    :: TRANSLATION

   ! . Local scalars.
   INTEGER :: I, NATOMS

   ! . Determine the number of atoms.
   NATOMS = SIZE ( COORDINATES, 2 )

   ! . Check the dimensions of COORDINATES.
   IF ( ( SIZE ( COORDINATES, 1 ) /= 3 ) .OR. ( NATOMS <  1 ) ) RETURN

   ! . Translate the coordinates.
   DO I = 1,NATOMS
      COORDINATES(1:3,I) = COORDINATES(1:3,I) + TRANSLATION(1:3)
   END DO

   END SUBROUTINE TRANSLATE

   !------------------------------------------------------
   SUBROUTINE TRANSLATE_TO_CENTER ( COORDINATES, WEIGHTS )
   !------------------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: COORDINATES

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: XYZ

   ! . Calculate the center of the system.
   XYZ = - CENTER ( COORDINATES, WEIGHTS )

   ! . Translate the coordinates to the center.
   CALL TRANSLATE ( COORDINATES, XYZ )

   END SUBROUTINE TRANSLATE_TO_CENTER

   !-----------------------------------------------------------
   SUBROUTINE TO_PRINCIPAL_AXES ( COORDINATES, WEIGHTS, PRINT )
   !-----------------------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: COORDINATES

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   LOGICAL :: QPRINT

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3)     :: MOMENTS
   REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: AXES

   ! . Check the dimensions of COORDINATES.
   IF ( ( SIZE ( COORDINATES, 1 ) /= 3 ) .OR. &
        ( SIZE ( COORDINATES, 2 ) <  1 ) ) RETURN

   ! . Move the system to its center of mass.
   CALL TRANSLATE_TO_CENTER ( COORDINATES, WEIGHTS )

   ! . Calculate the moments and axes of inertia.
   CALL MOMENTS_OF_INERTIA ( COORDINATES, MOMENTS, AXES, WEIGHTS )

   ! . Transpose the AXES matrix.
   AXES = TRANSPOSE ( AXES )

   ! . Rotate the coordinates.
   CALL ROTATE ( COORDINATES, AXES )

   ! . Set the print option.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Do some printing if necessary.
   IF ( QPRINT ) CALL PRINT_PARAGRAPH ( TEXT = "Coordinate set transformed to principal axes." )

   END SUBROUTINE TO_PRINCIPAL_AXES

END MODULE TRANSFORMATION
