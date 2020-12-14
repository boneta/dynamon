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
!                           The Diagonalization Module
!===============================================================================
!
! . Subroutines:
!
!   SYMMETRIC_UPPER              Diagonalize a symmetric matrix (upper triange).
!
! . Notes:
!
!   This module makes use of the LAPACK subroutine DSPEV.
!
!===============================================================================
MODULE DIAGONALIZATION

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_ERROR

IMPLICIT NONE
PRIVATE
PUBLIC :: SYMMETRIC_UPPER

# define LAPACK_DSPEV

!===============================================================================
CONTAINS
!===============================================================================

#ifdef LAPACK_DSPEV
   !---------------------------------------------------------------
   SUBROUTINE SYMMETRIC_UPPER ( MATRIX, EIGENVALUES, EIGENVECTORS )
   !---------------------------------------------------------------

   ! . Essential array argument declarations.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: MATRIX
   REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT)   :: EIGENVALUES

   ! . Optional array argument declarations.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT), OPTIONAL :: EIGENVECTORS

   ! . Local scalars.
   INTEGER :: IFAIL, N
   LOGICAL :: QOK

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: WORK
   REAL ( KIND = DP ),              DIMENSION(1:1) :: EVTEMP

   ! . Find the dimension of the matrix.
   N = SIZE ( EIGENVALUES )

   ! . Check the dimensions of the matrix and the other input arguments.
   QOK = ( ( N * ( N + 1 ) ) / 2 ) == SIZE ( MATRIX )
   IF ( PRESENT ( EIGENVECTORS ) ) THEN
      QOK = QOK .AND. ( SIZE ( EIGENVECTORS, 1 ) == N ) &
                .AND. ( SIZE ( EIGENVECTORS, 2 ) == N )
   END IF
   IF ( .NOT. QOK ) CALL PRINT_ERROR ( "SYMMETRIC_UPPER", "Array dimension error." )

   ! . Allocate the temporary arrays.
   ALLOCATE ( WORK(1:3*N) )

   ! . Initialize IFAIL.
   IFAIL = 0

   ! . The eigenvectors are to be found.
   IF ( PRESENT ( EIGENVECTORS ) ) THEN

      ! . Perform the diagonalization.
      CALL DSPEV ( "V", "U", N, MATRIX, EIGENVALUES, EIGENVECTORS, N, WORK, IFAIL )
     
   ! . Only the eigenvalues are to be found.
   ELSE

      ! . Perform the diagonalization.
      CALL DSPEV ( "N", "U", N, MATRIX, EIGENVALUES, EVTEMP,       N, WORK, IFAIL )

   END IF

   ! . Check for errors.
   IF ( IFAIL /= 0 ) CALL PRINT_ERROR ( "SYMMETRIC_UPPER", "Diagonalization error.", CODE = IFAIL )

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( WORK )

   END SUBROUTINE SYMMETRIC_UPPER
#else
   !---------------------------------------------------------------
   SUBROUTINE SYMMETRIC_UPPER ( MATRIX, EIGENVALUES, EIGENVECTORS )
   !---------------------------------------------------------------

   ! . Essential array argument declarations.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: MATRIX
   REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT)   :: EIGENVALUES

   ! . Optional array argument declarations.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT), OPTIONAL :: EIGENVECTORS

   ! . Local scalars.
   INTEGER :: IFAIL, N
   INTEGER :: I,J,LWORK,LIWORK
   LOGICAL :: QOK

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: WORK
   REAL ( KIND = DP ),              DIMENSION(1:1) :: EVTEMP
   INTEGER,            ALLOCATABLE, DIMENSION(:)   :: IWORK
   REAL ( KIND = DP ),              DIMENSION(1:10):: WORK_TMP
   INTEGER,                         DIMENSION(1:10):: IWORK_TMP

   ! . Find the dimension of the matrix.
   N = SIZE ( EIGENVALUES )
   ! . Check the dimensions of the matrix and the other input arguments.
   QOK = ( ( N * ( N + 1 ) ) / 2 ) == SIZE ( MATRIX )
   IF ( PRESENT ( EIGENVECTORS ) ) THEN
      QOK = QOK .AND. ( SIZE ( EIGENVECTORS, 1 ) == N ) &
                .AND. ( SIZE ( EIGENVECTORS, 2 ) == N )
   END IF
   IF ( .NOT. QOK ) CALL PRINT_ERROR ( "SYMMETRIC_UPPER", "Array dimension error." )

   ! . Initialize IFAIL.
   IFAIL = 0

   ! . The eigenvectors are to be found.
   IF ( PRESENT ( EIGENVECTORS ) ) THEN
      ! . Copy over the EigenVector
      do i=1,N
         do j=1,i
	    EIGENVECTORS(j,i) = MATRIX((i*(i-1))/2+j)
         enddo 
      enddo

      ! . Find the proper work variables
      LWORK  = -1
      LIWORK = -1
      CALL DSYEVD("V","U",N,EIGENVECTORS,N,EIGENVALUES,WORK_TMP, &
                  LWORK,IWORK_TMP,LIWORK,IFAIL)
      LWORK  = WORK_TMP(1)
      LIWORK = IWORK_TMP(1)
      
      ! . Allocate work space
      ALLOCATE(WORK(1:LWORK))
      ALLOCATE(IWORK(1:LIWORK)) 

      ! . Perform the diagonalization.
      CALL DSYEVD("V","U",N,EIGENVECTORS,N,EIGENVALUES,WORK, &
                  LWORK,IWORK,LIWORK,IFAIL)

      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
 
   ! . Only the eigenvalues are to be found.
   ELSE
      ! . Find the proper work variables
      LWORK  = -1
      LIWORK = -1
      CALL DSYEVD("N","U",N,EVTEMP,N,EIGENVALUES,WORK_TMP, &
                  LWORK,IWORK_TMP,LIWORK,IFAIL)
      LWORK  = WORK_TMP(1)
      LIWORK = IWORK_TMP(1)

      ! . Allocate work space
      ALLOCATE(WORK(1:LWORK))
      ALLOCATE(IWORK(1:LIWORK))

      ! . Perform the diagonalization.
      CALL DSYEVD("N","U",N,EVTEMP,N,EIGENVALUES,WORK, &
                  LWORK,IWORK,LIWORK,IFAIL)
   
      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)

   END IF

   ! . Check for errors.
   IF ( IFAIL /= 0 ) CALL PRINT_ERROR ( "SYMMETRIC_UPPER", "Diagonalization error.", CODE = IFAIL )

   ! . Deallocate the temporary arrays.
   END SUBROUTINE SYMMETRIC_UPPER
#endif
END MODULE DIAGONALIZATION
