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
!                             The Multipoles Module
!===============================================================================
!
! . Subroutines:
!
!   ATOM_CHARGES                    Calculate the atomic charges.
!   DIPOLE_DERIVATIVES              Calculate the dipole moment derivatives.
!   DIPOLE_MOMENT                   Calculate the dipole moment.
!
!===============================================================================
MODULE MULTIPOLES

! . Module declarations.
USE CONSTANTS,          ONLY : ANGSTROMS_TO_BOHRS, AU_TO_DB
USE DEFINITIONS,        ONLY : DP
USE ELEMENTS,           ONLY : SYMBOL
USE PRINTING,           ONLY : PRINT_LINE, PRINT_TABLE_ELEMENT, PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP

USE ATOMS,              ONLY : ATMCRD, ATMFIX, ATMMAS, ATMNUM, NATOMS, NATOMSMM
USE MM_TERMS,           ONLY : ATMCHG
USE QUANTUM_PROPERTIES, ONLY : QUANTUM_CHARGES, QUANTUM_DIPOLE, QUANTUM_DIPOLE_DERIVATIVES
USE TRANSFORMATION,     ONLY : CENTER_OF_MASS => CENTER

IMPLICIT NONE
PRIVATE
PUBLIC :: ATOM_CHARGES, DIPOLE_DERIVATIVES, DIPOLE_MOMENT

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------------------
   SUBROUTINE ATOM_CHARGES ( CHARGES, PRINT )
   !-----------------------------------------

   ! . Optional scalar argument declarations.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Optional array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:NATOMS), INTENT(OUT), OPTIONAL :: CHARGES

   ! . Local scalars.
   INTEGER :: I
   LOGICAL :: QPRINT

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: ATOMQ

   ! . Check that there are atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Check the PRINT argument.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Allocate and initialize the charge array.
   ALLOCATE ( ATOMQ(1:NATOMS) ) ; ATOMQ = 0.0_DP

   ! . Calculate the charges on the quantum atoms if any.
   CALL QUANTUM_CHARGES ( ATOMQ )

   ! . Add in any MM charge contributions.
   IF ( NATOMSMM > 0 ) ATOMQ = ATOMQ + ATMCHG

   ! . Print the results.
   IF ( QPRINT ) THEN
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 4, HEADER_COLOR = "#AAAA00" )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Atom Charges", COLSPAN = 4, HEADER = .TRUE. )
      DO I = 1,NATOMS
         WRITE ( PRINT_LINE, "(I4,2X,A,I3,2X,F7.3)" ) I, SYMBOL(ATMNUM(I)), ATMNUM(I), ATOMQ(I) ; CALL PRINT_TABLE_ELEMENT
      END DO
      CALL PRINT_TABLE_STOP
   END IF

   ! . Check the CHARGES argument.
   IF ( PRESENT ( CHARGES ) ) CHARGES = ATOMQ

   ! . Deallocate the charge arrays.
   DEALLOCATE ( ATOMQ )

   END SUBROUTINE ATOM_CHARGES

   !--------------------------------------
   SUBROUTINE DIPOLE_DERIVATIVES ( DMUDR )
   !--------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3*NATOMS,1:3),INTENT(OUT) :: DMUDR

   ! . Local scalars.
   INTEGER :: IATOM, II
   LOGICAL :: QOK

   ! . Check that there are atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Check the size of DMUDR and the validity of ATMCHG.
   QOK = ( SIZE ( DMUDR, 1 ) == 3 * NATOMS ) .AND. ( SIZE ( DMUDR, 2 ) == 3 ) .AND. &
         ( ALLOCATED ( ATMCHG ) ) .AND. ( SIZE ( ATMCHG ) == NATOMS )

   ! . Initialize DMUDR.
   DMUDR = 0.0_DP

   ! . Calculate the quantum dipole derivatives.
   CALL QUANTUM_DIPOLE_DERIVATIVES ( DMUDR )

   ! . Calculate the MM dipole derivatives.
   IF ( NATOMSMM > 0 ) THEN

      ! . Loop over the atoms.
      DO IATOM = 1,NATOMS

         ! . The atom is free.
         IF ( .NOT. ATMFIX(IATOM) ) THEN

            ! . Calculate the index into the array.
            II = 3 * ( IATOM - 1 )

            ! . The dipole derivatives are diagonal for each atom.
            DMUDR(II+1,1) = DMUDR(II+1,1) + ATMCHG(IATOM)
            DMUDR(II+2,2) = DMUDR(II+2,2) + ATMCHG(IATOM)
            DMUDR(II+3,3) = DMUDR(II+3,3) + ATMCHG(IATOM)

         END IF
      END DO
   END IF

   END SUBROUTINE DIPOLE_DERIVATIVES

   !-----------------------------------------
   SUBROUTINE DIPOLE_MOMENT ( DIPOLE, PRINT )
   !-----------------------------------------

   ! . Optional scalar argument declarations.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Optional array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(OUT), OPTIONAL :: DIPOLE

   ! . Local parameters.
   CHARACTER ( LEN = 1 ), DIMENSION(1:3), PARAMETER :: LABEL = (/ "X", "Y", "Z" /)

   ! . Local scalars.
   INTEGER            :: I, IATOM
   LOGICAL            :: QPRINT
   REAL ( KIND = DP ) :: DIPTOT, DTOTMM, DTOTQM

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: CENMAS, DIPMM, DIPMOM, DIPQM

   ! . Check that there are atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Check the PRINT argument.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Initialization.
   DIPMM = 0.0_DP
   DIPQM = 0.0_DP

   ! . Calculate the centre of mass of the system.
   CENMAS = CENTER_OF_MASS ( ATMCRD, ATMMAS )

   ! . Calculate the dipole moment due to the static charges.
   IF ( NATOMSMM > 0 ) THEN
      DO IATOM = 1,NATOMS
         DIPMM = DIPMM + ATMCHG(IATOM) * ( ATMCRD(1:3,IATOM) - CENMAS )
      END DO
   END IF

   ! . Convert the dipole moment to Debyes.
   DIPMM = ( ANGSTROMS_TO_BOHRS * AU_TO_DB ) * DIPMM

   ! . Calculate the total MM dipole moment.
   DTOTMM = SQRT ( DOT_PRODUCT ( DIPMM, DIPMM ) )

   ! . Calculate the quantum dipole.
   CALL QUANTUM_DIPOLE ( DIPQM, CENTER = CENMAS )

   ! . Calculate the total QM dipole moment.
   DTOTQM = SQRT ( DOT_PRODUCT ( DIPQM, DIPQM ) )

   ! . Calculate the total dipole moment for the system.
   DIPMOM = DIPMM + DIPQM
   DIPTOT = SQRT ( DOT_PRODUCT ( DIPMOM, DIPMOM ) )

   ! . Print the dipole moments.
   IF ( QPRINT ) THEN

      ! . Write out header.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#AAAA00", VARIABLEWIDTHS = (/ 10, 18, 18, 18, 18 /) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Dipole Moments (Debyes)", COLSPAN = 5, HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Component", HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Center",    HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Total",     HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "QM",        HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "MM",        HEADER = .TRUE. )

      ! . Write out the total dipole moments.
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Total", ALIGN = "LEFT", COLSPAN = 2 )
      WRITE ( PRINT_LINE, "(F18.4)" ) DIPTOT ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F18.4)" ) DTOTQM ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F18.4)" ) DTOTMM ; CALL PRINT_TABLE_ELEMENT

      ! . Write out the components.
      DO I = 1,3
         CALL PRINT_TABLE_ELEMENT ( TEXT = LABEL(I), ALIGN = "LEFT" )
	 WRITE ( PRINT_LINE, "(F18.4)" ) CENMAS(I) ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F18.4)" ) DIPMOM(I) ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F18.4)" ) DIPQM(I)  ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F18.4)" ) DIPMM(I)  ; CALL PRINT_TABLE_ELEMENT
      END DO
      
      ! . Write out the terminator.
      CALL PRINT_TABLE_STOP

   END IF

   ! . Check the DIPOLE argument.
   IF ( PRESENT ( DIPOLE ) ) DIPOLE = DIPMOM

   END SUBROUTINE DIPOLE_MOMENT

END MODULE MULTIPOLES
