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
!                              The Sequence Module
!===============================================================================
!
! . Parameter data:
!
!   RESIDUE_NAME_LENGTH            The maximum length of a residue name.
!   SUBSYSTEM_NAME_LENGTH          The maximum length of a subsystem name.
!
! . Scalar data:
!
!   NRESID                         The number of residues.
!   NSUBSYS                        The number of subsystems.
!
! . Array data:
!
!   RESIND                         The residue index array (to atoms).
!   RESNAM                         The residue names.
!
!   SUBIND                         The subsystem index array (to residues).
!   SUBNAM                         The subsystem names.
!
! . Subroutines:
!
!   SEQUENCE_ALLOCATE              Create the sequence data structure.
!   SEQUENCE_INITIALIZE            Initialize the sequence data structure.
!   SEQUENCE_PRINT                 Print the sequence.
!   SEQUENCE_SUMMARY               Print a summary of the sequence.
!   SEQUENCE_UNKNOWN               Create unknown sequence.
!   SEQUENCE_WRITE                 Write out a sequence file.
!
! . Notes:
!
!   The sequence module stores data about the residues and subsystems into
!   which a system is divided. In addition to the number and names of the
!   residues and subsystems there are index arrays.
!
!   The index arrays are structured as follows. For residue i, the element
!   RESIND(i)+1 gives the first atom of the residue and RESIND(i+1) gives
!   the last element. The total number of atoms in the residue is, thus,
!   RESIND(i+1) - RESIND(i). Note that the total length of the array is
!   NRESID + 1, not NRESID. The SUBIND array has the same structure, except
!   that it gives an index into the residue array for each subsystem.
!
!===============================================================================
MODULE SEQUENCE

! . Module declarations.
USE DEFINITIONS, ONLY : LINE_LENGTH
USE FILES,       ONLY : NEXT_UNIT
USE IO_UNITS,    ONLY : OUTPUT
USE PRINTING,    ONLY : PRINT_BLANKLINE, PRINT_ERROR, PRINT_LINE, PRINT_NOFORMAT_START, PRINT_NOFORMAT_STOP, PRINT_PARAGRAPH, &
                        PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, PRINT_SUMMARY_STOP,                &
			PRINT_TABLE_ELEMENT, PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP

USE ATOMS, ONLY : NATOMS

IMPLICIT NONE
PUBLIC
#ifndef PGPC
SAVE
#endif

! . Parameter declarations.
INTEGER, PARAMETER ::   RESIDUE_NAME_LENGTH = 32, &
                      SUBSYSTEM_NAME_LENGTH = 32

! . Scalar data.
INTEGER :: NRESID = 0, NSUBSYS = 0

! . Character array data.
CHARACTER ( LEN =   RESIDUE_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: RESNAM
CHARACTER ( LEN = SUBSYSTEM_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: SUBNAM

! . Integer array data.
INTEGER, ALLOCATABLE, DIMENSION(:) :: RESIND, SUBIND

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------------
   SUBROUTINE SEQUENCE_ALLOCATE ( NRES, NSUB )
   !------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: NRES, NSUB

   ! . Initialize the scalar variables.
   NRESID  = NRES   
   NSUBSYS = NSUB

   ! . Allocate all the arrays.
   ALLOCATE ( RESIND(1:NRESID+1), RESNAM(1:NRESID), SUBIND(1:NSUBSYS+1), SUBNAM(1:NSUBSYS) )

   END SUBROUTINE SEQUENCE_ALLOCATE

   !-----------------------------
   SUBROUTINE SEQUENCE_INITIALIZE
   !-----------------------------

   ! . Initialize the scalar variables.
   NRESID  = 0   
   NSUBSYS = 0

   ! . Deallocate all the arrays.
   IF ( ALLOCATED ( RESIND ) ) DEALLOCATE ( RESIND )
   IF ( ALLOCATED ( RESNAM ) ) DEALLOCATE ( RESNAM )
   IF ( ALLOCATED ( SUBIND ) ) DEALLOCATE ( SUBIND )
   IF ( ALLOCATED ( SUBNAM ) ) DEALLOCATE ( SUBNAM )

   END SUBROUTINE SEQUENCE_INITIALIZE

   !------------------------
   SUBROUTINE SEQUENCE_PRINT
   !------------------------

   ! . Local scalars.
   INTEGER :: I, ISUB, LENGTH

   ! . Check that a sequence is present.
   IF ( ( NRESID > 0 ) .AND. ( NSUBSYS > 0 ) ) THEN

      ! . Write out the header.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 2, HEADER_COLOR = "#FF0000", VARIABLEWIDTHS = (/ 40, 40 /) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "System Sequence", COLSPAN = 2, HEADER = .TRUE. )

      ! . Loop over the subsystems.
      DO ISUB = 1,NSUBSYS

         ! . Write out the subsystem name.
         LENGTH = LEN_TRIM ( SUBNAM(ISUB) )
         WRITE ( PRINT_LINE, "(A,I4,2X,A)" ) "Subsystem ", ISUB, SUBNAM(ISUB)(1:LENGTH)//":"
         CALL PRINT_TABLE_ELEMENT ( ALIGN = "LEFT", COLSPAN = 2, ITALIC = .TRUE. )

         ! . Write out the residue names.
	 DO I = (SUBIND(ISUB)+1),SUBIND(ISUB+1)
	    CALL PRINT_TABLE_ELEMENT ( ALIGN = "LEFT", TEXT = RESNAM(I) )
	 END DO
      END DO

      ! . Write out the terminator.
      CALL PRINT_TABLE_STOP

   END IF

   END SUBROUTINE SEQUENCE_PRINT

   !--------------------------
   SUBROUTINE SEQUENCE_SUMMARY
   !--------------------------

   ! . Check that a sequence is present.
   IF ( ( NRESID > 0 ) .AND. ( NSUBSYS > 0 ) ) THEN

      ! . Write out the data.
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF0000" )
      CALL PRINT_SUMMARY_START ( "Summary of Sequence Data" )
      WRITE ( PRINT_LINE, "(I14)" ) NSUBSYS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Subsystems" )
      WRITE ( PRINT_LINE, "(I14)" ) NRESID  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Residues"   )
      CALL PRINT_SUMMARY_STOP

   END IF

   END SUBROUTINE SEQUENCE_SUMMARY

   !--------------------------
   SUBROUTINE SEQUENCE_UNKNOWN
   !--------------------------

   IF ( ( NRESID <= 0 ) .OR. ( NSUBSYS <= 0 ) ) THEN
       CALL SEQUENCE_INITIALIZE
       CALL SEQUENCE_ALLOCATE ( 1, 1 )
       RESIND(1) = 0 ; RESIND(2) = NATOMS ; RESNAM(1) = "UNK"
       SUBIND(1) = 0 ; SUBIND(2) = 1      ; SUBNAM(1) = "UNKNOWN"
   END IF

   END SUBROUTINE SEQUENCE_UNKNOWN

   !---------------------------------
   SUBROUTINE SEQUENCE_WRITE ( FILE )
   !---------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Local scalars.
   CHARACTER ( LEN = LINE_LENGTH ) :: LINE
   INTEGER                         :: II, IOSTAT, IRES, ISUB, LENGTH, N, START, STOP, UNIT

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "SEQUENCE_WRITE", "I/O Error.", IOSTAT )

      ! . Do some printing.
      CALL PRINT_PARAGRAPH ( TEXT = "Sequence written to "//FILE )

   ! . The unit for writing is the output stream.
   ELSE

      ! . Assign the unit number.
      UNIT = OUTPUT

      ! . Do some printing.
      CALL PRINT_PARAGRAPH ( TEXT = "Sequence written to the output stream." )
      CALL PRINT_BLANKLINE
      CALL PRINT_NOFORMAT_START

   END IF

   ! . Write out the sequence file header.
   WRITE ( UNIT, "(A)" ) "Sequence"

   ! . Write out the number of subsystems.
   WRITE ( UNIT, "(I6)" ) NSUBSYS

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Write out the subsystem header.
      WRITE ( UNIT, "(/A)" ) "Subsystem "//TRIM ( SUBNAM(ISUB) )

      ! . Write out the number of residues in the subsystem.
      WRITE ( UNIT, "(I6)" ) ( SUBIND(ISUB+1) - SUBIND(ISUB) )

      ! . Find the longest residue name.
      LENGTH = 0
      DO IRES = (SUBIND(ISUB)+1),SUBIND(ISUB+1)
         LENGTH = MAX ( LENGTH, LEN_TRIM ( RESNAM(IRES) ) )
      END DO

      ! . Calculate the number of residue names that can be fit onto a line.
      N = ( LINE_LENGTH + 3 ) / ( LENGTH + 3 )

      ! . Top of the loop over residues.
      START = SUBIND(ISUB) + 1
      10 CONTINUE

         ! . Calculate the stopping position.
	 STOP = MIN ( ( START + N - 1 ), SUBIND(ISUB+1) )

         ! . Blank out the line.
	 LINE = REPEAT ( " ", LINE_LENGTH )

         ! . Fill the line.
	 DO IRES = START,STOP
	    II = ( IRES - START ) * ( LENGTH + 3 )
	    LINE(II+1:II+LENGTH) = RESNAM(IRES)(1:LENGTH)
	    IF ( IRES /= STOP ) LINE(II+LENGTH+1:II+LENGTH+3) = " ; "
	 END DO

         ! . Write out the lines.
         WRITE ( UNIT, "(A)" ) TRIM ( LINE )

         ! . Reset START.
	 START = STOP + 1

      ! . Check for termination.
      IF ( START <= SUBIND(ISUB+1) ) GO TO 10

      ! . Write out the subsystem terminator.
      WRITE ( UNIT, "(A)" ) "End"

   END DO

   ! . Write out the sequence terminator.
   WRITE ( UNIT, "(/A)" ) "End"

   ! . Finish up.
   IF ( UNIT == OUTPUT ) THEN
      CALL PRINT_NOFORMAT_STOP
   ELSE
      CLOSE ( UNIT )
   END IF

   END SUBROUTINE SEQUENCE_WRITE

END MODULE SEQUENCE
