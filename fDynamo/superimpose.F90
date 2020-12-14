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
!   The Module for the Comparison and the Superimposition of Coordinates Sets
!===============================================================================
!
! . Subroutines:
!
!   RMS_DEVIATION                  Calculate the RMS coordinate deviation.
!   SUPERIMPOSE_KABSCH             Superimpose two structures (Kabsch's method).
!   SUPERIMPOSE_QUATERNION         Superimpose two structures (quaternions).
!
! . Notes:
!
!   RMS_DEVIATION calculates the RMS coordinate deviation between two sets of
!   coordinate sets.
!
!   SUPERIMPOSE_KABSCH superimposes two coordinate sets. The algorithm used is
!   taken exactly from the papers:
!
!   Wolfgang Kabsch, `A solution for the best rotation to relate two sets of
!                     vectors', Acta Cryst. A32, 922--923, 1976.
!
!   Wolfgang Kabsch, `A discussion of the solution for the best rotation to
!                     relate two sets of vectors', Acta Cryst. A34, 827--828,
!                     1978.
!
!   SUPERIMPOSE_QUATERNION superimposes two coordinate sets using quaternions.
!   The reference for the method is:
!
!   Gerald Kneller, `Superposition of Molecular Structures using Quaternions',
!                    Mol. Sim. 7, 113--119, 1991.
!
!   All the subroutines take as input two sets of coordinates and an optional
!   array that can contain weights for use in the calculations.
!
!===============================================================================
MODULE SUPERIMPOSE

! . Module declarations.
USE DEFINITIONS,     ONLY : DP
USE DIAGONALIZATION, ONLY : SYMMETRIC_UPPER
USE IO_UNITS,        ONLY : OUTPUT
USE LINEAR_ALGEBRA,  ONLY : CROSS_PRODUCT, DETERMINANT, NORMALIZE
USE PRINTING,        ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH
USE TRANSFORMATION,  ONLY : CENTER, ROTATE, TRANSLATE, TRANSLATE_TO_CENTER

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------------------------------------
   SUBROUTINE RMS_DEVIATION ( SET1, SET2, WEIGHTS, RMS, PRINT )
   !-----------------------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: SET1
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: SET2

   ! . Optional scalar arguments.
   LOGICAL,            INTENT(IN),  OPTIONAL :: PRINT
   REAL ( KIND = DP ), INTENT(OUT), OPTIONAL :: RMS

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   INTEGER            :: I, NATOMS
   LOGICAL            :: QPRINT
   REAL ( KIND = DP ) :: DX, DY, DZ, RMST

   ! . Check the dimensions of the coordinate arrays.
   NATOMS = SIZE ( SET1, 2 )
   IF ( ( SIZE ( SET1, 1 ) /= 3 ) .OR. &
        ( SIZE ( SET2, 1 ) /= 3 ) .OR. &
        ( SIZE ( SET2, 2 ) /= NATOMS ) ) RETURN

   ! . Initialize the RMS counter.
   RMST = 0.0_DP

   ! . The WEIGHTS array is present.
   IF ( PRESENT ( WEIGHTS ) ) THEN

      ! . Check the dimensions of the WEIGHTS array.
      IF ( SIZE ( WEIGHTS ) /= NATOMS ) RETURN

      ! . Loop over the atoms.
      DO I = 1,NATOMS
         DX   = SET1(1,I) - SET2(1,I)
         DY   = SET1(2,I) - SET2(2,I)
         DZ   = SET1(3,I) - SET2(3,I)
         RMST = RMST + WEIGHTS(I) * ( DX*DX + DY*DY + DZ*DZ )
      END DO

      ! . Calculate the RMS.
      RMST = SQRT ( RMST / SUM ( WEIGHTS ) )

   ! . The WEIGHTS array is absent.
   ELSE

      ! . Loop over the atoms.
      DO I = 1,NATOMS
         DX   = SET1(1,I) - SET2(1,I)
         DY   = SET1(2,I) - SET2(2,I)
         DZ   = SET1(3,I) - SET2(3,I)
         RMST = RMST + ( DX*DX + DY*DY + DZ*DZ )
      END DO

      ! . Calculate the RMS.
      RMST = SQRT ( RMST / REAL ( NATOMS, DP ) )

   END IF

   ! . Check the RMS option.
   IF ( PRESENT ( RMS ) ) RMS = RMST

   ! . Set the print option.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Do some printing if necessary.
   IF ( QPRINT ) THEN
      WRITE ( PRINT_LINE, "(A,F20.5,A)" ) "RMS Difference = ", RMST, "."
      CALL PRINT_PARAGRAPH
   END IF

   END SUBROUTINE RMS_DEVIATION

   !-----------------------------------------------------------
   SUBROUTINE SUPERIMPOSE_KABSCH ( SET1, SET2, WEIGHTS, PRINT )
   !-----------------------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)    :: SET1
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: SET2

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   INTEGER :: NATOMS
   LOGICAL :: QPRINT

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: TMP1

   ! . Other local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3)     :: COM1
   REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: R, U

   ! . Check the dimensions of the coordinate arrays.
   NATOMS = SIZE ( SET1, 2 )
   IF ( ( SIZE ( SET1, 1 ) /= 3 ) .OR. &
        ( SIZE ( SET2, 1 ) /= 3 ) .OR. &
        ( SIZE ( SET2, 2 ) /= NATOMS ) ) RETURN

   ! . Check the dimensions of the WEIGHTS array.
   IF ( PRESENT ( WEIGHTS ) ) THEN
      IF ( SIZE ( WEIGHTS ) /= NATOMS ) RETURN
   END IF

   ! . Allocate the temporary array.
   ALLOCATE ( TMP1(1:3,1:NATOMS) )

   ! . Save the coordinates of SET1 in TMP1.
   TMP1 = SET1

   ! . Calculate the center of the first set.
   COM1 = CENTER ( TMP1, WEIGHTS )

   ! . Translate the sets to their centers of mass.
   CALL TRANSLATE           ( TMP1,  - COM1 )
   CALL TRANSLATE_TO_CENTER ( SET2, WEIGHTS )

   ! . Compute the Kabsch matrix for the two sets.
   CALL KABSCH_MATRIX

   ! . Deallocate the temporary array.
   DEALLOCATE ( TMP1 )

   ! . Find the rotation matrix.
   CALL KABSCH_ROTATION

   ! . Rotate the coordinates of the second set.
   CALL ROTATE ( SET2, U )

   ! . Translate the coordinates of the second set to the center of mass of the first set.
   CALL TRANSLATE ( SET2, COM1 )

   ! . Set the print option.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Do some printing if necessary.
   IF ( QPRINT ) CALL PRINT_PARAGRAPH ( TEXT = "Coordinate set 2 superimposed onto coordinate set 1." )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !-----------------------
      SUBROUTINE KABSCH_MATRIX
      !-----------------------

      ! . Local scalars.
      INTEGER            :: I, J, P
      REAL ( KIND = DP ) :: SUM

      ! . The WEIGHTS array is present.
      IF ( PRESENT ( WEIGHTS ) ) THEN

         DO I = 1,3
            DO J = 1,3
               SUM = 0.0_DP
               DO P = 1,NATOMS
                  SUM = SUM + WEIGHTS(P) * TMP1(J,P) * SET2(I,P)
               END DO
               R(J,I) = SUM
            END DO
         END DO

      ! . The WEIGHTS array is absent.
      ELSE

         DO I = 1,3
            DO J = 1,3
               SUM = 0.0_DP
               DO P = 1,NATOMS
                  SUM = SUM + TMP1(J,P) * SET2(I,P)
               END DO
               R(J,I) = SUM
            END DO
         END DO

      END IF

      END SUBROUTINE KABSCH_MATRIX

      !-------------------------
      SUBROUTINE KABSCH_ROTATION
      !-------------------------

      ! . Local scalars.
      INTEGER            :: I, IJ, J
      LOGICAL            :: QSECOND
      REAL ( KIND = DP ) :: DET

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3)     :: MU
      REAL ( KIND = DP ), DIMENSION(1:6)     :: RTR
      REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: A, B, RA

      ! . Construct a positive definite symmetric matrix by multiplying R by its transpose.
      A = MATMUL ( TRANSPOSE ( R ), R )

      ! . Extract the upper triangle.
      IJ = 0
      DO I = 1,3
         DO J = 1,I
            IJ = IJ + 1
            RTR(IJ) = A(J,I)
         END DO
      END DO

      ! . Diagonalize RTR.
      CALL SYMMETRIC_UPPER ( RTR, EIGENVALUES = MU, EIGENVECTORS = B )

      ! . Put the eigenvectors in order of descending eigenvalue.
      A(1:3,1) = B(1:3,3) ; A(1:3,2) = B(1:3,2) ; A(1:3,3) = B(1:3,1)
      DET      = MU(1)    ; MU(1)    = MU(3)    ; MU(3)    = DET

      ! . Reset the eigenvector with the lowest eigenvalue.
      A(1:3,3) = CROSS_PRODUCT ( A(1:3,1), A(1:3,2) )

      ! . Determine RA = R * A.
      RA = MATMUL ( R, A )

      ! . Normalize B1 and B2.
      DO I = 1,2
         B(1:3,I) = NORMALIZE ( RA(1:3,I) )
      END DO

      ! . Calculate B3.
      B(1:3,3) = CROSS_PRODUCT ( B(1:3,1), B(1:3,2) )

      ! . If necessary change the sign of the lowest-eigenvalue eigenvector.
      IF ( DOT_PRODUCT ( B(1:3,3), RA(1:3,3) ) < 0.0_DP ) B(1:3,3) = - B(1:3,3)

      ! . Initialize QSECOND.
      QSECOND = .FALSE.

      ! . Top of loop over U.
      10 CONTINUE

      ! . Form U.
      U = MATMUL ( B, TRANSPOSE ( A ) )

      ! . Check the determinant of U.
      DET = DETERMINANT ( U )
      IF ( ( ABS ( DET - 1.0_DP ) > 1.0E-6_DP ) .OR. &
           ( ABS ( SUM ( MATMUL ( U, TRANSPOSE ( U ) ) ) - 3.0_DP ) > 1.0E-6_DP ) ) THEN
         ! . This is the second time so signal an error.
	 IF ( QSECOND ) THEN
            WRITE ( OUTPUT, "(/A,G10.4)" ) "Invalid rotation matrix determinant = ", DET
  	    CALL PRINT_ERROR ( "SUPERIMPOSE_KABSCH", "Invalid rotation matrix." )
         ! . Redo U by resetting B and flagging QSECOND.
         ELSE
	    QSECOND  = .TRUE.
	    B(1:3,3) = - B(1:3,3)
	    GO TO 10
         END IF
      END IF

      END SUBROUTINE KABSCH_ROTATION
  
   END SUBROUTINE SUPERIMPOSE_KABSCH

   !---------------------------------------------------------------
   SUBROUTINE SUPERIMPOSE_QUATERNION ( SET1, SET2, WEIGHTS, PRINT )
   !---------------------------------------------------------------

   ! . Essential array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)    :: SET1
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: SET2

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   INTEGER :: NATOMS
   LOGICAL :: QPRINT

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: TMP1

   ! . Other local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3)     :: COM1
   REAL ( KIND = DP ), DIMENSION(1:10)    :: M
   REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: U

   ! . Check the dimensions of the coordinate arrays.
   NATOMS = SIZE ( SET1, 2 )
   IF ( ( SIZE ( SET1, 1 ) /= 3 ) .OR. &
        ( SIZE ( SET2, 1 ) /= 3 ) .OR. &
        ( SIZE ( SET2, 2 ) /= NATOMS ) ) RETURN

   ! . Check the dimensions of the WEIGHTS array.
   IF ( PRESENT ( WEIGHTS ) ) THEN
      IF ( SIZE ( WEIGHTS ) /= NATOMS ) RETURN
   END IF

   ! . Allocate the temporary array.
   ALLOCATE ( TMP1(1:3,1:NATOMS) )

   ! . Save the coordinates of SET1 in TMP1.
   TMP1 = SET1

   ! . Calculate the center of the first set.
   COM1 = CENTER ( TMP1, WEIGHTS )

   ! . Translate the sets to their centers of mass.
   CALL TRANSLATE           ( TMP1,  - COM1 )
   CALL TRANSLATE_TO_CENTER ( SET2, WEIGHTS )

   ! . Compute the matrix to be diagonalized.
   CALL CONSTRUCT_MATRIX

   ! . Deallocate the temporary array.
   DEALLOCATE ( TMP1 )

   ! . Find the rotation matrix.
   CALL CONSTRUCT_ROTATION

   ! . Rotate the coordinates of the second set.
   CALL ROTATE ( SET2, U )

   ! . Translate the coordinates of the second set to the center of mass of the first set.
   CALL TRANSLATE ( SET2, COM1 )

   ! . Set the print option.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Do some printing if necessary.
   IF ( QPRINT ) CALL PRINT_PARAGRAPH ( TEXT = "Coordinate set 2 superimposed onto coordinate set 1." )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !--------------------------
      SUBROUTINE CONSTRUCT_MATRIX
      !--------------------------

      ! . Local scalars.
      INTEGER            :: P
      REAL ( KIND = DP ) :: FACT

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3)  :: C
      REAL ( KIND = DP ), DIMENSION(1:10) :: R

      ! . Initialize M.
      M = 0.0_DP

      ! . The WEIGHTS array is present.
      IF ( PRESENT ( WEIGHTS ) ) THEN

         ! . Loop over the atoms.
         DO P = 1,NATOMS

            ! . Calculate some intermediate factors.
            C    = 2.0_DP * CROSS_PRODUCT ( TMP1(1:3,P), SET2(1:3,P) )
            FACT = SUM ( ( TMP1(1:3,P) + SET2(1:3,P) )**2 )

            ! . Calculate the elements of the matrix.
            R(1)  = SUM ( ( TMP1(1:3,P) - SET2(1:3,P) )**2 )
            R(2)  = C(1)
            R(3)  = FACT - 4.0_DP * TMP1(1,P) * SET2(1,P)
            R(4)  = C(2)
            R(5)  = - 2.0_DP * ( TMP1(1,P) * SET2(2,P) + TMP1(2,P) * SET2(1,P) )
            R(6)  = FACT - 4.0_DP * TMP1(2,P) * SET2(2,P)
            R(7)  = C(3)
            R(8)  = - 2.0_DP * ( TMP1(1,P) * SET2(3,P) + TMP1(3,P) * SET2(1,P) )
            R(9)  = - 2.0_DP * ( TMP1(2,P) * SET2(3,P) + TMP1(3,P) * SET2(2,P) )
            R(10) = FACT - 4.0_DP * TMP1(3,P) * SET2(3,P)
            
            ! . Add in the elements to the full matrix.
            M = M + WEIGHTS(P) * R

         END DO

      ! . The WEIGHTS array is absent.
      ELSE

         ! . Loop over the atoms.
         DO P = 1,NATOMS

            ! . Calculate some intermediate factors.
            C    = 2.0_DP * CROSS_PRODUCT ( TMP1(1:3,P), SET2(1:3,P) )
            FACT = SUM ( ( TMP1(1:3,P) + SET2(1:3,P) )**2 )

            ! . Calculate the elements of the matrix.
            R(1)  = SUM ( ( TMP1(1:3,P) - SET2(1:3,P) )**2 )
            R(2)  = C(1)
            R(3)  = FACT - 4.0_DP * TMP1(1,P) * SET2(1,P)
            R(4)  = C(2)
            R(5)  = - 2.0_DP * ( TMP1(1,P) * SET2(2,P) + TMP1(2,P) * SET2(1,P) )
            R(6)  = FACT - 4.0_DP * TMP1(2,P) * SET2(2,P)
            R(7)  = C(3)
            R(8)  = - 2.0_DP * ( TMP1(1,P) * SET2(3,P) + TMP1(3,P) * SET2(1,P) )
            R(9)  = - 2.0_DP * ( TMP1(2,P) * SET2(3,P) + TMP1(3,P) * SET2(2,P) )
            R(10) = FACT - 4.0_DP * TMP1(3,P) * SET2(3,P)
            
            ! . Add in the elements to the full matrix.
            M = M + R

         END DO

      END IF

      END SUBROUTINE CONSTRUCT_MATRIX

      !----------------------------
      SUBROUTINE CONSTRUCT_ROTATION
      !----------------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: DET, Q0, Q1, Q2, Q3

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:4)     :: MU
      REAL ( KIND = DP ), DIMENSION(1:4,1:4) :: A

      ! . Diagonalize M.
      CALL SYMMETRIC_UPPER ( M, EIGENVALUES = MU, EIGENVECTORS = A )

      ! . Construct the rotation matrix using the values in the eigenvector with the lowest eigenvalue.
      Q0 = A(1,1) ; Q1 = A(2,1) ; Q2 = A(3,1) ; Q3 = A(4,1)
      U(1,1) = Q0**2 + Q1**2 - Q2**2 - Q3**2
      U(2,1) = 2.0_DP * (   Q0*Q3 + Q1*Q2 )
      U(3,1) = 2.0_DP * ( - Q0*Q2 + Q1*Q3 )
      U(1,2) = 2.0_DP * ( - Q0*Q3 + Q1*Q2 )
      U(2,2) = Q0**2 - Q1**2 + Q2**2 - Q3**2
      U(3,2) = 2.0_DP * (   Q0*Q1 + Q2*Q3 )
      U(1,3) = 2.0_DP * (   Q0*Q2 + Q1*Q3 )
      U(2,3) = 2.0_DP * ( - Q0*Q1 + Q2*Q3 )
      U(3,3) = Q0**2 - Q1**2 - Q2**2 + Q3**2

      ! . Check the determinant of U.
      DET = DETERMINANT ( U )
      IF ( ( ABS ( DET - 1.0_DP ) > 1.0E-6_DP ) .OR. &
           ( ABS ( SUM ( MATMUL ( U, TRANSPOSE ( U ) ) ) - 3.0_DP ) > 1.0E-6_DP ) ) THEN
         WRITE ( OUTPUT, "(/A,G10.4)" ) "Invalid rotation matrix determinant = ", DET
	 CALL PRINT_ERROR ( "SUPERIMPOSE_QUATERNION", "Invalid rotation matrix." )
      END IF

      END SUBROUTINE CONSTRUCT_ROTATION
  
   END SUBROUTINE SUPERIMPOSE_QUATERNION

END MODULE SUPERIMPOSE
