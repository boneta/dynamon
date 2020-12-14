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
!                            The MM I/O System Module
!===============================================================================
!
! . Subroutines:
!
!   MM_SYSTEM_READ                 Read a system file.
!   MM_SYSTEM_WRITE                Write a system file.
!
! . Notes:
!
!   MM_SYSTEM_READ and MM_SYSTEM_WRITE read and write an MM system file which
!   contains information about the ATOMS, MM_TERMS and SEQUENCE. SYMMETRY
!   information is not contained although MM_SYSTEM_READ will initialize the
!   symmetry module.
!
!===============================================================================
MODULE MM_SYSTEM_IO

! . Utility data structure declarations.
USE DEFINITIONS, ONLY : DP, FORCE_FIELD, MAX_RECORD_LENGTH, VERSION
USE FILES,       ONLY : NEXT_UNIT
USE PRINTING,    ONLY : PRINT_ERROR, PRINT_PARAGRAPH

USE ATOMS
USE MM_TERMS
USE SEQUENCE
USE SYMMETRY

IMPLICIT NONE
PUBLIC

!===============================================================================
CONTAINS
!===============================================================================

   !----------------------------------------
   SUBROUTINE MM_SYSTEM_READ ( FILE, PRINT )
   !----------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Local I/O scalars.
   CHARACTER ( LEN = LEN ( FORCE_FIELD ) ) :: HEADER
   INTEGER                                 :: IOSTAT, UNIT
   LOGICAL                                 :: QPRINT
   REAL ( KIND = DP )                      :: VERNUM

   ! . Local counters.
   INTEGER :: I, NANG, NATM, NBND, NDIH, NEXCL, NEX14, NIMP, NRES, NSUB

   ! . Initialize all the data structures.
   CALL    ATOMS_INITIALIZE
   CALL MM_TERMS_INITIALIZE
   CALL SEQUENCE_INITIALIZE
   CALL SYMMETRY_INITIALIZE

   ! . Get the next unit number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, ACTION = "READ", FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "OLD", IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "MM_SYSTEM_READ", "I/O Error.", IOSTAT )

   ! . Header.
   READ ( UNIT ) HEADER, VERNUM
   IF ( ( HEADER /= FORCE_FIELD ) .OR. ( VERNUM > VERSION ) ) CALL PRINT_ERROR ( "MM_SYSTEM_READ", "Invalid system file header." )

   ! . Residue and subsystem data.
   READ ( UNIT ) NRES, NSUB
   CALL SEQUENCE_ALLOCATE ( NRES, NSUB )
   READ ( UNIT ) RESIND
   READ ( UNIT ) RESNAM
   READ ( UNIT ) SUBIND
   READ ( UNIT ) SUBNAM

   ! . Atom data.
   READ ( UNIT ) NATM
   CALL ATOMS_ALLOCATE ( NATM )
   READ ( UNIT ) ATMMAS
   READ ( UNIT ) ATMNAM
   READ ( UNIT ) ATMNUM

   ! . MM terms data.
   READ ( UNIT ) NANG, NBND, NDIH, NIMP, NEXCL, NEX14, SCALE_EL14, SCALE_LJ14
   CALL MM_TERMS_ALLOCATE ( NBND, NANG, NDIH, NIMP )
   ALLOCATE ( ATMEXCJ(1:NEXCL), ATME14J(1:NEX14) )
   READ ( UNIT ) ATMCHG
   READ ( UNIT ) ATMCHG14
   READ ( UNIT ) ATMEPS
   READ ( UNIT ) ATMEPS14
   READ ( UNIT ) ATMSIG
   READ ( UNIT ) ATMTYP

   READ ( UNIT ) ATMEXCI
   READ ( UNIT ) ATMEXCJ
   READ ( UNIT ) ATME14I
   READ ( UNIT ) ATME14J

   DO I = 1,NANGLES,    MM_IO_BATCH ; READ ( UNIT )    ANGLES(I:MIN(I+MM_IO_BATCH,NANGLES))    ; ENDDO
   DO I = 1,NBONDS,     MM_IO_BATCH ; READ ( UNIT )     BONDS(I:MIN(I+MM_IO_BATCH,NBONDS))     ; ENDDO
   DO I = 1,NDIHEDRALS, MM_IO_BATCH ; READ ( UNIT ) DIHEDRALS(I:MIN(I+MM_IO_BATCH,NDIHEDRALS)) ; ENDDO
   DO I = 1,NIMPROPERS, MM_IO_BATCH ; READ ( UNIT ) IMPROPERS(I:MIN(I+MM_IO_BATCH,NIMPROPERS)) ; ENDDO

   ! . Close the file.
   CLOSE ( UNIT )

   ! . Set the print flag.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Print out the system data if required.
   IF ( QPRINT ) THEN

      ! . Write out some general information.
      CALL PRINT_PARAGRAPH ( TEXT = "System data read from " // FILE )

      ! . Write out some specific information.
      CALL ATOMS_SUMMARY
      CALL MM_TERMS_SUMMARY
      CALL SEQUENCE_SUMMARY
      CALL SYMMETRY_SUMMARY

   END IF

   END SUBROUTINE MM_SYSTEM_READ

   !----------------------------------
   SUBROUTINE MM_SYSTEM_WRITE ( FILE )
   !----------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Local I/O scalars.
   INTEGER :: I, IOSTAT, UNIT

   ! . Local counters.
   INTEGER :: NEXCL, NEX14

   ! . Get the next unit number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "UNKNOWN", IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "MM_SYSTEM_WRITE", "I/O Error.", IOSTAT )

   ! . Header.
   WRITE ( UNIT ) FORCE_FIELD, VERSION

   ! . Residue and subsystem data.
   WRITE ( UNIT ) NRESID, NSUBSYS
   WRITE ( UNIT ) RESIND
   WRITE ( UNIT ) RESNAM
   WRITE ( UNIT ) SUBIND
   WRITE ( UNIT ) SUBNAM

   ! . Atom data.
   WRITE ( UNIT ) NATOMS
   WRITE ( UNIT ) ATMMAS
   WRITE ( UNIT ) ATMNAM
   WRITE ( UNIT ) ATMNUM

   ! . Determine the size of the exclusion arrays.
   NEXCL = SIZE ( ATMEXCJ )
   NEX14 = SIZE ( ATME14J )

   ! . MM terms data.
   WRITE ( UNIT ) NANGLES, NBONDS, NDIHEDRALS, NIMPROPERS, NEXCL, NEX14, SCALE_EL14, SCALE_LJ14
   WRITE ( UNIT ) ATMCHG
   WRITE ( UNIT ) ATMCHG14
   WRITE ( UNIT ) ATMEPS
   WRITE ( UNIT ) ATMEPS14
   WRITE ( UNIT ) ATMSIG
   WRITE ( UNIT ) ATMTYP

   WRITE ( UNIT ) ATMEXCI
   WRITE ( UNIT ) ATMEXCJ
   WRITE ( UNIT ) ATME14I
   WRITE ( UNIT ) ATME14J

   DO I = 1,NANGLES,    MM_IO_BATCH ; WRITE ( UNIT )    ANGLES(I:MIN(I+MM_IO_BATCH,NANGLES))    ; ENDDO
   DO I = 1,NBONDS,     MM_IO_BATCH ; WRITE ( UNIT )     BONDS(I:MIN(I+MM_IO_BATCH,NBONDS))     ; ENDDO
   DO I = 1,NDIHEDRALS, MM_IO_BATCH ; WRITE ( UNIT ) DIHEDRALS(I:MIN(I+MM_IO_BATCH,NDIHEDRALS)) ; ENDDO
   DO I = 1,NIMPROPERS, MM_IO_BATCH ; WRITE ( UNIT ) IMPROPERS(I:MIN(I+MM_IO_BATCH,NIMPROPERS)) ; ENDDO

   ! . Close the file.
   CLOSE ( UNIT )

   ! . Write out some information.
   CALL PRINT_PARAGRAPH ( TEXT = "System data written to " // FILE )

   END SUBROUTINE MM_SYSTEM_WRITE

END MODULE MM_SYSTEM_IO
