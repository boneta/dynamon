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
!                                The Atoms Module
!===============================================================================
!
! . Parameter data:
!
!   ATOM_NAME_LENGTH                The maximum length of an atom name.
!
! . Scalar data:
!
!   NATOMS                          The total number of atoms.
!   NFIXED                          The number of fixed atoms.
!   NFREE                           The number of free atoms (NATOMS - NFIXED).
!
!   NATOMSMM                        The number of MM atoms.
!   NATOMSQM                        The number of QM atoms.
!
! . Array data:
!
!   ATMCRD                          The coordinates of the atoms.
!   ATMFIX                          The array to indicate which atoms are fixed.
!   ATMIND                          The free atom index array.
!   ATMMAS                          The atom masses.
!   ATMNAM                          The atom names.
!   ATMNUM                          The atom numbers.
!   ATMQMI                          The indices of the QM atoms.
!
! . Subroutines:
!
!   ATOMS_ALLOCATE                  Allocate the atoms variables.
!   ATOMS_FIX                       Fill the ATMFIX array.
!   ATOMS_FORMULA                   Write out the formula for the structure.
!   ATOMS_INITIALIZE                Initialize the atoms variables.
!   ATOMS_SUMMARY                   Print a summary of the atom data.
!
! . Notes:
!
!   The atoms module stores the data for the atoms in the system. In addition
!   to the number of atoms in the system, for each atom there are its mass,
!   its name and its coordinates. All data is public and, hence, accessible.
!
!   All atoms are supposed to be MM by default.
!
!===============================================================================
MODULE ATOMS

! . Module declarations.
USE CONSTANTS,   ONLY : UNDEFINED
USE DEFINITIONS, ONLY : DP
USE ELEMENTS,    ONLY : NELEMENTS, SYMBOL
USE PRINTING,    ONLY : PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS,         &
                        PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT, PRINT_TABLE_ENDLINE, &
                        PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP

IMPLICIT NONE
PUBLIC
#ifndef PGPC
SAVE
#endif

! . Parameter declarations.
INTEGER, PARAMETER :: ATOM_NAME_LENGTH = 8

! . Scalar data.
INTEGER :: NATOMS = 0, NATOMSMM = 0, NATOMSQM = 0, NFIXED = 0, NFREE = 0

! . Character array data.
CHARACTER ( LEN = ATOM_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: ATMNAM

! . Integer array data.
INTEGER, ALLOCATABLE, DIMENSION(:) :: ATMIND, ATMNUM, ATMQMI

! . Logical array data.
LOGICAL, ALLOCATABLE, DIMENSION(:) :: ATMFIX

! . Real array data.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: ATMMAS
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: ATMCRD

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------
   SUBROUTINE ATOMS_ALLOCATE ( N )
   !------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: N

   ! . Local scalars.
   INTEGER :: I

   ! . Initialize the scalar variables.
   NATOMS   = N
   NATOMSMM = N
   NATOMSQM = 0
   NFIXED   = 0
   NFREE    = N

   ! . Allocate all the arrays.
   ALLOCATE ( ATMCRD(1:3,1:NATOMS), ATMFIX(1:NATOMS), ATMIND(1:NATOMS), ATMMAS(1:NATOMS), ATMNAM(1:NATOMS), &
              ATMNUM(1:NATOMS),     ATMQMI(1:NATOMS) )

   ! . Initialize all the atom arrays except ATMIND.
   ATMCRD = UNDEFINED
   ATMFIX = .FALSE.
   ATMMAS = 0.0_DP
   ATMNAM = "        "
   ATMNUM = 0
   ATMQMI = 0

   ! . Initialize ATMIND.
   DO I = 1,NATOMS
      ATMIND(I) = I
   END DO

   END SUBROUTINE ATOMS_ALLOCATE

   !----------------------------
   SUBROUTINE ATOMS_FIX ( QFIX )
   !----------------------------

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: QFIX

   ! . Local scalars.
   INTEGER :: I, IFREE

   ! . There are atoms.
   IF ( NATOMS > 0 ) THEN

      ! . Copy the QFIX array to ATMFIX.
      ATMFIX = QFIX

      ! . Set the number of fixed atoms.
      NFIXED = COUNT ( ATMFIX )

      ! . Set the number of free atoms
      NFREE = NATOMS - NFIXED

      ! . Fill the ATMIND array.
      IFREE = 0
      DO I = 1,NATOMS
         IF ( ATMFIX(I) ) THEN
	    ATMIND(I) = 0
	 ELSE
	    IFREE = IFREE + 1
	    ATMIND(I) = IFREE
	 END IF
      END DO

   ! . There are no atoms.
   ELSE
      NFIXED = 0
   END IF

   ! . Write out some data.
   WRITE ( PRINT_LINE, "(A,I6,A)" ) "The coordinates of ", NFIXED, " atoms have been fixed."
   CALL PRINT_PARAGRAPH

   END SUBROUTINE ATOMS_FIX

   !-----------------------
   SUBROUTINE ATOMS_FORMULA
   !-----------------------

   ! . Local scalars.
   INTEGER :: I, N, NOTHER

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:) :: FREQUENCY, INDEX

   ! . Check that there are atoms.
   IF ( NATOMS > 0 ) THEN

      ! . Allocate and initialize some temporary arrays.
      ALLOCATE ( FREQUENCY(1:NELEMENTS), INDEX(1:NELEMENTS) )
      FREQUENCY = 0

      ! . Loop over the atoms.
      NOTHER = 0
      DO I = 1,NATOMS
         IF ( ( ATMNUM(I) <= 0 ) .OR. ( ATMNUM(I) > NELEMENTS ) ) THEN
            NOTHER = NOTHER + 1
         ELSE
            FREQUENCY(ATMNUM(I)) = FREQUENCY(ATMNUM(I)) + 1
         END IF
      END DO

      ! . Contract the frequency array and create the index array.
      N = 0
      DO I = 1,NELEMENTS
         IF ( FREQUENCY(I) > 0 ) THEN
            N = N + 1
            FREQUENCY(N) = FREQUENCY(I)
            INDEX(N)     = I
         END IF
      END DO

      ! . Write out the formula.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 8 )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( COLOR = "#FF0000", TEXT = "Atom Formula", COLSPAN = 8, HEADER = .TRUE. )
      DO I = 1,N
         WRITE ( PRINT_LINE, "(2X,A2,I6)" ) SYMBOL(INDEX(I)), FREQUENCY(I)
	 CALL PRINT_TABLE_ELEMENT
      END DO
      CALL PRINT_TABLE_ENDLINE
      IF ( NOTHER > 0 ) THEN
         WRITE ( PRINT_LINE, "(A,I12,A)" ) "Number of Unknown Atoms = ", NOTHER, "."
	 CALL PRINT_TABLE_ELEMENT ( ALIGN = "LEFT", COLSPAN = 8 )
      END IF
      CALL PRINT_TABLE_STOP

      ! . Deallocate the temporary arrays.
      DEALLOCATE ( FREQUENCY, INDEX )

   END IF

   END SUBROUTINE ATOMS_FORMULA

   !--------------------------
   SUBROUTINE ATOMS_INITIALIZE
   !--------------------------

   ! . Initialize the scalar variables.
   NATOMS   = 0
   NATOMSMM = 0
   NATOMSQM = 0
   NFIXED   = 0
   NFREE    = 0

   ! . Deallocate all the arrays.
   IF ( ALLOCATED ( ATMCRD ) ) DEALLOCATE ( ATMCRD )
   IF ( ALLOCATED ( ATMFIX ) ) DEALLOCATE ( ATMFIX )
   IF ( ALLOCATED ( ATMIND ) ) DEALLOCATE ( ATMIND )
   IF ( ALLOCATED ( ATMMAS ) ) DEALLOCATE ( ATMMAS )
   IF ( ALLOCATED ( ATMNAM ) ) DEALLOCATE ( ATMNAM )
   IF ( ALLOCATED ( ATMNUM ) ) DEALLOCATE ( ATMNUM )
   IF ( ALLOCATED ( ATMQMI ) ) DEALLOCATE ( ATMQMI )

   END SUBROUTINE ATOMS_INITIALIZE

   !-----------------------
   SUBROUTINE ATOMS_SUMMARY
   !-----------------------

   ! . Local scalars.
   INTEGER            :: NHEAVY, NHYDRO, NOTHER
   REAL ( KIND = DP ) :: TOTMAS

   ! . Check that there are some atoms.
   IF ( NATOMS > 0 ) THEN

      ! . Determine the number of atoms in each class and the total mass.
      NHYDRO = COUNT ( ATMNUM(1:NATOMS) == 1 )
      NOTHER = COUNT ( ATMNUM(1:NATOMS) <= 0 ) + COUNT ( ATMNUM(1:NATOMS) > NELEMENTS )
      NHEAVY = NATOMS - NHYDRO - NOTHER
      TOTMAS = SUM ( ATMMAS(1:NATOMS) )

      ! . Write out the data.
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "FF0000" )
      CALL PRINT_SUMMARY_START ( "Summary of Atom Data" )
      WRITE ( PRINT_LINE, "(I14)"   ) NATOMS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Atoms"       )
      WRITE ( PRINT_LINE, "(I14)"   ) NHYDRO ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Hydrogens"   )
      WRITE ( PRINT_LINE, "(I14)"   ) NHEAVY ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Heavy Atoms" )
      WRITE ( PRINT_LINE, "(I14)"   ) NOTHER ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Unknowns"    )
      WRITE ( PRINT_LINE, "(I14)"   ) NFIXED ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Fixed Atoms" )
      WRITE ( PRINT_LINE, "(F14.1)" ) TOTMAS ; CALL PRINT_SUMMARY_ELEMENT ( "Total Mass (a.m.u.)"   )
      CALL PRINT_SUMMARY_STOP

   END IF

   END SUBROUTINE ATOMS_SUMMARY

END MODULE ATOMS
