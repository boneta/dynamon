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
!                The Molecular Mechanics System Definition Module
!===============================================================================
!
! . Parameter data:
!
!   MM_IO_BATCH                    The maximum number of MM terms in a binary
!                                  file read or write operation.
!
! . Atom data:
!
!   ATMTYP                         The atom types (for reference only).
!
! . Bond data:
!
!   NBONDS                         The number of bonds.
!   BONDS                          The bond terms.
!
! . Angle data:
!
!   NANGLES                        The number of angles.
!   ANGLES                         The angle terms.
!
! . Dihedral data:
!
!   NDIHEDRALS                     The number of dihedrals.
!   DIHEDRALS                      The dihedral terms.
!
! . Improper data:
!
!   NIMPROPERS                     The number of impropers
!   IMPROPERS                      The improper terms.
!
! . Non-bonding data:
!
!   ATMEXCI, ATMEXCJ               The atom 1-2, 1-3 and 1-4 exclusions.
!   ATME14I, ATME14J               The atom 1-4 interactions.
!   ATMCHG                         The atom charges.
!   ATMCHG14                       ATMCHG for 1-4 interactions.
!   ATMEPS                         The atom LJ epsilon (twice the square roots).
!   ATMEPS14                       ATMEPS for 1-4 interactions.
!   ATMSIG                         The atom LJ sigma values (the square roots).
!   SCALE_EL14                     The electrostatics 1-4 scale factor.
!   SCALE_LJ14                     The LJ 1-4 scale factor.
!
! . Routines:
!
!   MM_TERMS_ALLOCATE              Allocate the MM terms data structure.
!   MM_TERMS_INITIALIZE            Initialize the data structure.
!   MM_TERMS_SUMMARY               Write a summary of the data structure.
!
! . Notes:
!
! 1. MM_TERMS stores all the data necessary for the calculation of the MM
!    energy. There are covalent energy terms as well as parameters for the
!    calculation of the non-bond energy.
!
! 2. The force constants stored for the internal energy terms are used as
!    is. There is no need to multiply them by a half or any other factor.
!
!===============================================================================
MODULE MM_TERMS

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, &
                        PRINT_SUMMARY_STOP

USE ATOMS,       ONLY : ATOM_NAME_LENGTH, NATOMS

IMPLICIT NONE
PUBLIC
#ifndef PGPC
SAVE
#endif

! . Angle term type definition.
TYPE ANGLE_TERM
   INTEGER            :: I, J, K
   REAL ( KIND = DP ) :: EQ, FC
END TYPE ANGLE_TERM

! . Bond term type definition.
TYPE BOND_TERM
   INTEGER            :: I, J
   REAL ( KIND = DP ) :: EQ, FC
END TYPE BOND_TERM

! . Dihedral term type definition.
TYPE DIHEDRAL_TERM
   INTEGER            :: I, J, K, L, PERIOD
   REAL ( KIND = DP ) :: FC, PHASE
END TYPE DIHEDRAL_TERM

! . Module parameters.
INTEGER, PARAMETER :: MM_IO_BATCH = 1024

! . The atom data.
CHARACTER ( LEN = ATOM_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: ATMTYP

! . The bond data.
INTEGER                                    :: NBONDS
TYPE(BOND_TERM), ALLOCATABLE, DIMENSION(:) ::  BONDS

! . The angle data.
INTEGER                                     :: NANGLES
TYPE(ANGLE_TERM), ALLOCATABLE, DIMENSION(:) ::  ANGLES

! . The dihedral data.
INTEGER                                        :: NDIHEDRALS
TYPE(DIHEDRAL_TERM), ALLOCATABLE, DIMENSION(:) ::  DIHEDRALS

! . The improper data.
INTEGER                                        :: NIMPROPERS
TYPE(DIHEDRAL_TERM), ALLOCATABLE, DIMENSION(:) ::  IMPROPERS

! . Scalar non-bonding data.
REAL ( KIND = DP ) :: SCALE_EL14 = 1.0_DP, SCALE_LJ14 = 1.0_DP

! . Array non-bonding data.
INTEGER,            ALLOCATABLE, DIMENSION(:) :: ATMEXCI, ATMEXCJ, ATME14I, ATME14J
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: ATMCHG, ATMCHG14, ATMEPS, ATMEPS14, ATMSIG

!===============================================================================
CONTAINS
!===============================================================================

   !------------------------------------------------------
   SUBROUTINE MM_TERMS_ALLOCATE ( NBND, NANG, NDIH, NIMP )
   !------------------------------------------------------

   ! . Scalar argument declarations.
   INTEGER, INTENT(IN) :: NBND, NANG, NDIH, NIMP

   ! . Assign the counters.
   NANGLES    = NANG
   NBONDS     = NBND
   NDIHEDRALS = NDIH
   NIMPROPERS = NIMP

   ! . Allocate the atom arrays.
   ALLOCATE ( ATMCHG(1:NATOMS),   ATMCHG14(1:NATOMS), ATMEPS(1:NATOMS), &
              ATMEPS14(1:NATOMS), ATMSIG(1:NATOMS),   ATMTYP(1:NATOMS)  )

   ! . Allocate the atom exclusion arrays.
   ALLOCATE ( ATMEXCI(1:NATOMS+1), ATME14I(1:NATOMS+1) )

   ! . Allocate the MM terms arrays.
   ALLOCATE ( ANGLES(1:NANGLES), BONDS(1:NBONDS), DIHEDRALS(1:NDIHEDRALS), IMPROPERS(1:NIMPROPERS) )

   END SUBROUTINE MM_TERMS_ALLOCATE

   !-----------------------------
   SUBROUTINE MM_TERMS_INITIALIZE
   !-----------------------------

   ! . Initialize the MM term counters.
   NANGLES    = 0
   NBONDS     = 0
   NDIHEDRALS = 0
   NIMPROPERS = 0

   ! . Initialize the non-bond options.
   SCALE_EL14 = 1.0_DP
   SCALE_LJ14 = 1.0_DP

   ! . Initialize the atom arrays.
   IF ( ALLOCATED ( ATMCHG   ) ) DEALLOCATE ( ATMCHG   )
   IF ( ALLOCATED ( ATMCHG14 ) ) DEALLOCATE ( ATMCHG14 )
   IF ( ALLOCATED ( ATMEPS   ) ) DEALLOCATE ( ATMEPS   )
   IF ( ALLOCATED ( ATMEPS14 ) ) DEALLOCATE ( ATMEPS14 )
   IF ( ALLOCATED ( ATMSIG   ) ) DEALLOCATE ( ATMSIG   )
   IF ( ALLOCATED ( ATMTYP   ) ) DEALLOCATE ( ATMTYP   )

   ! . Initialize the atom exclusion arrays.
   IF ( ALLOCATED ( ATMEXCI ) ) DEALLOCATE ( ATMEXCI )
   IF ( ALLOCATED ( ATMEXCJ ) ) DEALLOCATE ( ATMEXCJ )
   IF ( ALLOCATED ( ATME14I ) ) DEALLOCATE ( ATME14I )
   IF ( ALLOCATED ( ATME14J ) ) DEALLOCATE ( ATME14J )

   ! . Initialize the MM term arrays.
   IF ( ALLOCATED ( ANGLES    ) ) DEALLOCATE ( ANGLES    )
   IF ( ALLOCATED ( BONDS     ) ) DEALLOCATE ( BONDS     )
   IF ( ALLOCATED ( DIHEDRALS ) ) DEALLOCATE ( DIHEDRALS )
   IF ( ALLOCATED ( IMPROPERS ) ) DEALLOCATE ( IMPROPERS )

   END SUBROUTINE MM_TERMS_INITIALIZE

   !--------------------------
   SUBROUTINE MM_TERMS_SUMMARY
   !--------------------------

   ! . Write out a summary of the MM terms.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#AA00AA" )
   CALL PRINT_SUMMARY_START ( "Summary of MM Terms Data" )
   WRITE ( PRINT_LINE, "(I14)" ) NBONDS     ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Bonds"	 )
   WRITE ( PRINT_LINE, "(I14)" ) NANGLES    ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Angles"	 )
   WRITE ( PRINT_LINE, "(I14)" ) NDIHEDRALS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Dihedrals" )
   WRITE ( PRINT_LINE, "(I14)" ) NIMPROPERS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Impropers" )
   WRITE ( PRINT_LINE, "(I14)" ) SIZE ( ATMEXCJ ) - SIZE ( ATME14J )
   CALL PRINT_SUMMARY_ELEMENT ( "1-2/1-3 Exclusions" )
   WRITE ( PRINT_LINE, "(I14)" ) SIZE ( ATME14J )
   CALL PRINT_SUMMARY_ELEMENT ( "1-4 Exclusions"     )
   WRITE ( PRINT_LINE, "(F14.4)" ) SCALE_EL14 ; CALL PRINT_SUMMARY_ELEMENT ( "1-4 Scale Elect." )
   WRITE ( PRINT_LINE, "(F14.4)" ) SCALE_LJ14 ; CALL PRINT_SUMMARY_ELEMENT ( "1-4 Scale LJ"     )
   WRITE ( PRINT_LINE, "(F14.4)" ) SUM ( ATMCHG(1:NATOMS) )
   CALL PRINT_SUMMARY_ELEMENT ( "Total System Charge" )
   CALL PRINT_SUMMARY_STOP

   END SUBROUTINE MM_TERMS_SUMMARY

END MODULE MM_TERMS
