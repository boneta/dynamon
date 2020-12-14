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
!                    The Finite Difference Derivatives Module
!===============================================================================
!
! . Subroutines:
!
!   NUMERICAL_GRADIENT             Calculate a finite difference gradient.
!   NUMERICAL_HESSIAN              Calculate a finite difference Hessian.
!
! . Notes:
!
!   The subroutines use a simple central difference formula to calculate the
!   first and second derivatives of the potential energy with respect to the
!   atomic coordinates.
!
!===============================================================================
MODULE NUMERICAL_DERIVATIVES

! . Module declarations.
USE DEFINITIONS,      ONLY : DP

USE ATOMS,            ONLY : ATMCRD, ATMFIX, NATOMS, NFIXED, NFREE
USE POTENTIAL_ENERGY, ONLY : ATMDER, ENERGY, ETOTAL, GRADIENT

IMPLICIT NONE
PRIVATE
PUBLIC :: NUMERICAL_GRADIENT, NUMERICAL_HESSIAN

!===============================================================================
CONTAINS
!===============================================================================

   !-----------------------------------------------------------
   SUBROUTINE NUMERICAL_GRADIENT ( DELTA, GRADIENT, SELECTION )
   !-----------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: DELTA

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(OUT) :: GRADIENT

   ! . Optional array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN), OPTIONAL :: SELECTION

   ! . Local scalars.
   INTEGER            :: I, IATOM
   LOGICAL            :: QSELECT
   REAL ( KIND = DP ) :: EMINUS, EPLUS, TEMP

   ! . Check for SELECTION.
   QSELECT = PRESENT ( SELECTION )

   ! . Initialize GRADIENT.
   GRADIENT = 0.0_DP

   ! . Loop over the atoms.
   DO IATOM = 1,NATOMS

      ! . Check for fixed atoms.
      IF ( ATMFIX(IATOM) ) CYCLE

      ! . Check for a selected atom.
      IF ( QSELECT ) THEN
         IF ( .NOT. SELECTION(IATOM) ) CYCLE
      END IF

      ! . Loop over the Cartesian components.
      DO I = 1,3

         ! . Save the coordinate.
         TEMP = ATMCRD(I,IATOM)

         ! . Increment the coordinate.
         ATMCRD(I,IATOM) = TEMP + DELTA

         ! . Calculate and save the energy.
         CALL ENERGY ( PRINT = .FALSE. ) ; EPLUS = ETOTAL

         ! . Decrement the coordinate.
         ATMCRD(I,IATOM) = TEMP - DELTA

         ! . Calculate and save the energy.
         CALL ENERGY ( PRINT = .FALSE. ) ; EMINUS = ETOTAL

         ! . Restore the coordinate.
         ATMCRD(I,IATOM) = TEMP

         ! . Calculate the component of the gradient.
         GRADIENT(I,IATOM) = EPLUS - EMINUS

      END DO
   END DO

   ! . Scale the gradient array.
   GRADIENT = GRADIENT / ( 2.0_DP * DELTA )

   END SUBROUTINE NUMERICAL_GRADIENT

   !----------------------------------------------
   SUBROUTINE NUMERICAL_HESSIAN ( DELTA, HESSIAN )
   !----------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: DELTA

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(OUT) :: HESSIAN

   ! . Local scalars.
   INTEGER            :: I, IATOM, IFREE, INDEX, J
   REAL ( KIND = DP ) :: TEMP

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: GMINUS, SD_MATRIX

   ! . Allocate and initialize some temporary arrays.
   ALLOCATE ( GMINUS(1:3,1:NFREE), SD_MATRIX(1:3*NFREE,1:3*NFREE) ) ; SD_MATRIX = 0.0_DP

   ! . Initialization.
   IFREE = 0

   ! . Loop over the atoms.
   DO IATOM = 1,NATOMS

      ! . Skip fixed atoms.
      IF ( ATMFIX(IATOM) ) CYCLE

      ! . Increment the free atom counter.
      IFREE = IFREE + 1

      ! . Loop over the Cartesian components.
      DO I = 1,3

         ! . Save the coordinate.
         TEMP = ATMCRD(I,IATOM)

         ! . Decrement the coordinate.
         ATMCRD(I,IATOM) = TEMP - DELTA

         ! . Calculate and save the gradient.
         CALL GRADIENT ( PRINT = .FALSE. ) ; GMINUS = CONTRACT ( )

         ! . Increment the coordinate.
         ATMCRD(I,IATOM) = TEMP + DELTA

         ! . Calculate and save the gradient.
         CALL GRADIENT ( PRINT = .FALSE. )

         ! . Restore the coordinate.
         ATMCRD(I,IATOM) = TEMP

         ! . Calculate the components of the second derivative matrix.
         INDEX = 3 * ( IFREE - 1 ) + I
         SD_MATRIX(1:3*NFREE,INDEX) = PACK ( ( CONTRACT ( ) - GMINUS ), .TRUE. )

      END DO
   END DO

   ! . Copy the elements of SD_MATRIX to HESSIAN.
   INDEX = 0
   DO I = 1,3*NFREE
      DO J = 1,I
         INDEX = INDEX + 1
         HESSIAN(INDEX) = 0.5_DP * ( SD_MATRIX(I,J) + SD_MATRIX(J,I) )
      END DO
   END DO

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( GMINUS, SD_MATRIX )

   ! . Scale the hessian.
   HESSIAN = HESSIAN / ( 2.0_DP * DELTA )

   !============================================================================
   CONTAINS
   !============================================================================

      !--------------------
      FUNCTION CONTRACT ( )
      !--------------------

      ! . Function declarations.
      REAL ( KIND = DP ), DIMENSION(1:3,1:NFREE) :: CONTRACT

      ! . Local scalars.
      INTEGER :: IATOM, IFREE

      ! . There are no fixed atoms.
      IF ( NFIXED == 0 ) THEN
         CONTRACT = ATMDER
      ! . There are fixed atoms.
      ELSE

         ! . Loop over the atoms.
         IFREE = 0
         DO IATOM = 1,NATOMS
            IF ( ATMFIX(IATOM) ) CYCLE
            IFREE = IFREE + 1
            CONTRACT(1:3,IFREE) = ATMDER(1:3,IATOM)
         END DO

      END IF

      END FUNCTION CONTRACT

   END SUBROUTINE NUMERICAL_HESSIAN

END MODULE NUMERICAL_DERIVATIVES
