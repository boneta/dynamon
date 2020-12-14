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
!                              The Z-Matrix Module
!===============================================================================
!
! . Subroutines:
!
!   ZMATRIX_DEFINE                  Define a system from a Z-matrix file.
!   ZMATRIX_READ                    Read a Z-matrix file.
!   ZMATRIX_WRITE                   Write a Z-matrix file.
!
! . Notes:
!
!   The three ZMATRIX I/O subroutines work in exactly the same way as the
!   equivalent subroutines in the module COORDINATE_IO except that there
!   are no optional arguments for the subroutines ZMATRIX_READ and
!   ZMATRIX_WRITE as there are in the corresponding COORDINATE_READ and
!   COORDINATE_WRITE cases.
!
!===============================================================================
MODULE ZMATRIX_IO

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE ELEMENTS,    ONLY : MASS
USE FILES,       ONLY : NEXT_UNIT
USE IO_UNITS,    ONLY : INPUT, OUTPUT
USE PARSING
USE PRINTING,    ONLY : PRINT_BLANKLINE, PRINT_ERROR, PRINT_NOFORMAT_START, PRINT_NOFORMAT_STOP, &
                        PRINT_PARAGRAPH

USE ATOMS
USE SEQUENCE
USE SYMMETRY
USE ZMATRIX

IMPLICIT NONE
PRIVATE
PUBLIC :: ZMATRIX_DEFINE, ZMATRIX_READ, ZMATRIX_WRITE

!==============================================================================
CONTAINS
!==============================================================================

   !---------------------------------
   SUBROUTINE ZMATRIX_DEFINE ( FILE )
   !---------------------------------

   ! . Optional scalar arugments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Local scalars.
   INTEGER :: I, IATOM, IOSTAT, IRES, ISUB, JATOM, JSUB, NATM, NATMIN, NRES, NRESIN, NUMATM, NUMRES, NSUBIN, UNIT

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "ZMATRIX_DEFINE", "I/O Error.", IOSTAT )

   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Initialize the atom, sequence, symmetry and Z-matrix data structures.
   CALL ATOMS_INITIALIZE
   CALL SEQUENCE_INITIALIZE
   CALL SYMMETRY_INITIALIZE
   CALL ZMATRIX_INITIALIZE

   ! . Read in the first line of the file.
   CALL GET_LINE
   CALL GET_INTEGER ( NATMIN )
   CALL GET_INTEGER ( NRESIN )
   CALL GET_INTEGER ( NSUBIN )

   ! . Check the counters.
   IF ( ( NATMIN <= 0 ) .OR. ( NRESIN <= 0 ) .OR. ( NSUBIN <= 0 ) ) THEN
      CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "Invalid atoms, residue or subsystem counter." )
   END IF

   ! . Allocate the atom, sequence and Z-matrix data arrays.
   CALL ATOMS_ALLOCATE ( NATMIN )
   CALL SEQUENCE_ALLOCATE ( NRESIN, NSUBIN )
   CALL ZMATRIX_ALLOCATE

   ! . Read in any symmetry definitions.
   CALL SYMMETRY_RECORD_READ

   ! . Initialization.
   NUMATM = 0
   NUMRES = 0

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Get the first line of the subsystem.
      IF ( ISUB > 1 ) THEN
         CALL GET_LINE
         CALL GET_WORD
      END IF
      IF ( WORD(1:WRDLEN) /= "SUBSYSTEM" ) THEN
         CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "SUBSYSTEM label invalid." )
      END IF
      CALL GET_INTEGER ( I )
      CALL GET_WORD
      IF ( ( WRDLEN <= 0 ) .OR. ( WRDLEN > SUBSYSTEM_NAME_LENGTH ) ) THEN
         CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "SUBSYSTEM name length invalid." )
      END IF

      ! . Check the subsystem names.
      DO JSUB = 1,(ISUB-1)
         ! . Compare the names.
         IF ( SUBNAM(JSUB) == WORD(1:WRDLEN) ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "Two SUBSYSTEM names are the same." )
         END IF
      END DO

      ! . Save the subsystem name.
      SUBNAM(ISUB) = WORD(1:WRDLEN)

      ! . Fill the subsystem arrays.
      SUBIND(ISUB) = NUMRES

      ! . Get the number of residues in the subsystem.
      CALL GET_LINE
      CALL GET_INTEGER ( NRES )
      IF ( ( NRES <= 0 ) .OR. ( ( NUMRES+NRES ) > NRESID ) ) THEN
         CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "Invalid number of residues in a subsystem." )
      END IF

      ! . Loop over the residues in the subsystem.
      DO IRES = 1,NRES

         ! . Increment the NUMRES counter.
         NUMRES = NUMRES + 1

         ! . Get the first line of the residue.
         CALL GET_LINE
         CALL GET_WORD
         IF ( WORD(1:WRDLEN) /= "RESIDUE" ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "RESIDUE label invalid." )
         END IF
         CALL GET_INTEGER ( I )
         IF ( I /= IRES ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "RESIDUE number invalid." )
         END IF
         CALL GET_WORD
         IF ( ( WRDLEN <= 0 ) .OR. ( WRDLEN > RESIDUE_NAME_LENGTH ) ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "RESIDUE name length invalid." )
         END IF

         ! . Save the residue name.
         RESNAM(NUMRES) = WORD(1:WRDLEN)

         ! . Fill the residue arrays.
         RESIND(NUMRES) = NUMATM

         ! . Get the number of atoms in the residue.
         CALL GET_LINE
         CALL GET_INTEGER ( NATM )
         IF ( ( NATM <= 0 ).OR.( (NUMATM+NATM) > NATOMS ) ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "Invalid number of atoms in a residue." )
         END IF

         ! . Loop over the atoms/cards in the residue.
         DO IATOM = 1,NATM

            ! . Increment the atom counter.
            NUMATM = NUMATM + 1

            ! . Parse an atom line ( the order, name and number ).
            CALL GET_LINE
            CALL GET_INTEGER ( I )
            CALL GET_WORD
            IF ( WRDLEN <= ATOM_NAME_LENGTH ) THEN
               ATMNAM(NUMATM) = WORD(1:WRDLEN)
            ELSE
               CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "ATOM name length invalid." )
            END IF

            ! . Check the atom names.
            DO JATOM = (RESIND(NUMRES)+1),(NUMATM-1)
               ! . Compare the names.
               IF ( ATMNAM(JATOM) == WORD(1:WRDLEN) ) THEN
                  CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "Two ATOM names are the same." )
               END IF
            END DO

            ! . Read in the remaining data.
            CALL GET_INTEGER ( ATMNUM(NUMATM)          )
            CALL GET_REAL    ( ZMCARDS(NUMATM)%BOND     )
            CALL GET_REAL    ( ZMCARDS(NUMATM)%ANGLE    )
            CALL GET_REAL    ( ZMCARDS(NUMATM)%DIHEDRAL )
            CALL GET_INTEGER ( ZMCARDS(NUMATM)%J        )
            CALL GET_INTEGER ( ZMCARDS(NUMATM)%K        )
            CALL GET_INTEGER ( ZMCARDS(NUMATM)%L        )

            ! . Save the order of the card.
            ZMCARDS(NUMATM)%IORDER = NUMATM

         END DO
      END DO
   END DO

   ! . Check NUMATM and NUMRES.
   IF ( ( NUMATM /= NATOMS ) .OR. ( NUMRES /= NRESID ) ) THEN
      CALL PARSE_ERROR ( "ZMATRIX_DEFINE", "ATOM or RESIDUE number mismatch." )
   END IF

   ! . Fill the last RESIND and SUBIND elements.
   RESIND(NRESID+1)  = NUMATM
   SUBIND(NSUBSYS+1) = NUMRES

   ! . Fill the atom mass array.
   DO IATOM = 1,NATOMS
      IF ( ATMNUM(IATOM) > 0 ) THEN
         ATMMAS(IATOM) = MASS(ATMNUM(IATOM))
      ELSE
         ATMMAS(IATOM) = 0.0_DP
      END IF
   END DO

   ! . Close the input file if necessary.
   IF ( UNIT /= INPUT ) CLOSE ( UNIT )

   ! . Reset the parsing unit.
   CALL POP_UNIT

   ! . Check the Z-matrix syntax.
   CALL ZMATRIX_CHECK

   ! . Do some printing.
   IF ( PRESENT ( FILE ) ) THEN
      CALL PRINT_PARAGRAPH ( TEXT = "System read from " // FILE )
   ELSE
      CALL PRINT_PARAGRAPH ( TEXT = "System read from the input stream." )
   END IF

   END SUBROUTINE ZMATRIX_DEFINE

   !-------------------------------
   SUBROUTINE ZMATRIX_READ ( FILE )
   !-------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Local scalars.
   INTEGER :: I, IATOM, IOSTAT, IRES, ISUB, JATOM, NATM, NATMIN, NRES, NRESIN, NSUBIN, NUMATM, UNIT

   ! . Return if no residues have been defined.
   IF ( NRESID <= 0 ) RETURN

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . Check for an error.
      IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "ZMATRIX_READ", "I/O Error.", IOSTAT )

   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Initialize the Z-matrix data structure.
   CALL ZMATRIX_INITIALIZE

   ! . Read in the first line of the file.
   CALL GET_LINE
   CALL GET_INTEGER ( NATMIN )
   CALL GET_INTEGER ( NRESIN )
   CALL GET_INTEGER ( NSUBIN )

   ! . Check the counters.
   IF ( ( NATMIN > NATOMS ) .OR. ( NRESIN /= NRESID ) .OR. ( NSUBIN /= NSUBSYS ) ) THEN
      CALL PARSE_ERROR ( "ZMATRIX_READ", "Invalid atoms, residue or subsystem counter." )
   END IF

   ! . Reallocate the Z-matrix data structure.
   CALL ZMATRIX_ALLOCATE

   ! . Read in any symmetry definitions.
   CALL SYMMETRY_RECORD_READ

   ! . Initialization.
   NUMATM = 0

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Get the first line of the subsystem.
      IF ( ISUB > 1 ) THEN
         CALL GET_LINE
         CALL GET_WORD
      END IF
      IF ( WORD(1:WRDLEN) /= "SUBSYSTEM" ) THEN
         CALL PARSE_ERROR ( "ZMATRIX_READ", "SUBSYSTEM label invalid." )
      END IF
      CALL GET_INTEGER ( I )
      CALL GET_WORD
      IF ( ( WRDLEN <= 0 ) .OR. ( WRDLEN > SUBSYSTEM_NAME_LENGTH ) ) THEN
         CALL PARSE_ERROR ( "ZMATRIX_READ", "SUBSYSTEM name length invalid." )
      END IF

      ! . Check the subsystem name.
      IF ( SUBNAM(ISUB) /= WORD(1:WRDLEN) ) THEN
         CALL PARSE_ERROR ( "ZMATRIX_READ", "SUBSYSTEM name mismatch." )
      END IF

      ! . Get the number of residues in the subsystem.
      CALL GET_LINE
      CALL GET_INTEGER ( NRES )
      IF ( NRES /= ( SUBIND(ISUB+1) - SUBIND(ISUB) ) ) THEN
         CALL PARSE_ERROR ( "ZMATRIX_READ", "Invalid number of residues in a subsystem." )
      END IF

      ! . Loop over the residues in the subsystem.
      DO IRES = (SUBIND(ISUB)+1),SUBIND(ISUB+1)

         ! . Get the first line of the residue.
         CALL GET_LINE
         CALL GET_WORD
         IF ( WORD(1:WRDLEN) /= "RESIDUE" ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_READ", "RESIDUE label invalid." )
         END IF
         CALL GET_INTEGER ( I )
         IF ( I /= ( IRES - SUBIND(ISUB) ) ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_READ", "RESIDUE number invalid." )
         END IF
         CALL GET_WORD
         IF ( ( WRDLEN <= 0 ) .OR. ( WRDLEN > RESIDUE_NAME_LENGTH ) ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_READ", "RESIDUE name length invalid." )
         END IF

         ! . Check the residue name.
         IF ( RESNAM(IRES) /= WORD(1:WRDLEN) ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_READ", "RESIDUE name mismatch." )
         END IF

         ! . Get the number of atoms in the residue.
         CALL GET_LINE
         CALL GET_INTEGER ( NATM )
         IF ( NATM > ( RESIND(IRES+1) - RESIND(IRES) ) ) THEN
            CALL PARSE_ERROR ( "ZMATRIX_READ", "Too many atoms in a residue." )
         END IF

         ! . Read in the cards in the residue.
         DO IATOM = 1,NATM

            ! . Increment the atom counter.
            NUMATM = NUMATM + 1

            ! . Parse an atom line ( Atom's order, name and number ).
            CALL GET_LINE
            CALL GET_INTEGER ( I )
            CALL GET_WORD
            IF ( WRDLEN > ATOM_NAME_LENGTH ) THEN
               CALL PARSE_ERROR ( "ZMATRIX_READ", "ATOM name length invalid." )
            ENDIF

            ! . Loop over the atoms in the definition.
            DO JATOM = (RESIND(IRES)+1),RESIND(IRES+1)

               ! . Check the name.
               IF ( ATMNAM(JATOM) == WORD(1:WRDLEN) ) THEN

                  ! . Parse the remaining data.
                  CALL GET_INTEGER ( I )
                  IF ( ATMNUM(JATOM) /= I ) THEN
                     CALL PARSE_ERROR ( "ZMATRIX_READ", "ATOM number mismatch." )
                  END IF

                  ! . Read in the remaining data.
                  CALL GET_REAL    ( ZMCARDS(NUMATM)%BOND     )
                  CALL GET_REAL    ( ZMCARDS(NUMATM)%ANGLE    )
                  CALL GET_REAL    ( ZMCARDS(NUMATM)%DIHEDRAL )
                  CALL GET_INTEGER ( ZMCARDS(NUMATM)%J        )
                  CALL GET_INTEGER ( ZMCARDS(NUMATM)%K        )
                  CALL GET_INTEGER ( ZMCARDS(NUMATM)%L        )

                  ! . Save the order of the card.
                  ZMCARDS(NUMATM)%IORDER = JATOM

                  ! . Skip out of the loop.
                  GO TO 10

               END IF
            END DO

            ! . No match was found.
            CALL PARSE_ERROR ( "ZMATRIX_READ", "ATOM name mismatch." )

            ! . End of the loop.
            10 CONTINUE

         END DO
      END DO
   END DO

   ! . Close the input file if necessary.
   IF ( UNIT /= INPUT ) CLOSE ( UNIT )

   ! . Reset the parsing unit.
   CALL POP_UNIT

   ! . Check the Z-matrix syntax.
   CALL ZMATRIX_CHECK

   ! . Do some printing.
   IF ( PRESENT ( FILE ) ) THEN
      CALL PRINT_PARAGRAPH ( TEXT = "Z-matrix read from " // FILE )
   ELSE
      CALL PRINT_PARAGRAPH ( TEXT = "Z-matrix read from the input stream." )
   END IF

   END SUBROUTINE ZMATRIX_READ

   !--------------------------------
   SUBROUTINE ZMATRIX_WRITE ( FILE )
   !--------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Local scalars.
   INTEGER            :: IATOM, ICARD, IOSTAT, IRES, ISUB, JCARD, KCARD, LCARD, NATM, NRES, UNIT
   REAL ( KIND = DP ) :: PHI, R, THETA

   ! . Return if there are no Z-matrix cards or residues.
   IF ( ( NCARDS <= 0 ) .OR. ( NRESID <= 0 ) ) RETURN

   ! . Check that NATOMS and NCARDS are consistent.
   IF ( NATOMS /= NCARDS ) CALL PRINT_ERROR ( "ZMATRIX_WRITE", "The ATOMS and ZMATRIX modules are inconsistent." )

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "ZMATRIX_WRITE", "I/O Error.", IOSTAT )

      ! . Do some printing.
      CALL PRINT_PARAGRAPH ( TEXT = "Z-matrix written to " // FILE )

   ! . The unit for writing is the output stream.
   ELSE

      ! . Assign the unit number.
      UNIT = OUTPUT

      ! . Do some printing.
      CALL PRINT_PARAGRAPH ( TEXT = "Z-matrix written to the output stream." )
      CALL PRINT_BLANKLINE
      CALL PRINT_NOFORMAT_START

   END IF

   ! . Write out the separator.
   WRITE ( UNIT, "('!',79('='))" )

   ! . Write out the number of atoms, residues and subsystems.
   WRITE ( UNIT, "(3I6,A)" ) NATOMS, NRESID, NSUBSYS, " ! # of atoms, residues and subsystems."

   ! . Write out any symmetry information.
   CALL SYMMETRY_RECORD_WRITE ( UNIT )

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Get the number of residues.
      NRES = SUBIND(ISUB+1) - SUBIND(ISUB)

      ! . Write out a subsystem header.
      WRITE ( UNIT, "('!',79('='))" )
      WRITE ( UNIT, "(A,I6,2X,A)" ) "Subsystem", ISUB, SUBNAM(ISUB)
      WRITE ( UNIT, "(I6,A)" ) NRES, " ! # of residues."
      WRITE ( UNIT, "('!',79('='))" )

      ! . Loop over the residues.
      DO IRES = (SUBIND(ISUB)+1),SUBIND(ISUB+1)

         ! . Get the number of atoms.
         NATM = RESIND(IRES+1) - RESIND(IRES)

         ! . Write out a residue header.
         WRITE ( UNIT, "(A,I6,2X,A)" ) "Residue", ( IRES - SUBIND(ISUB) ), RESNAM(IRES)
         WRITE ( UNIT, "(I6,A)" ) NATM, " ! # of atoms."

         ! . Loop over the cards in the residue.
         DO ICARD = (RESIND(IRES)+1),RESIND(IRES+1)

            ! . Get the card indices.
            JCARD = ZMCARDS(ICARD)%J
            KCARD = ZMCARDS(ICARD)%K
            LCARD = ZMCARDS(ICARD)%L

            ! . Get the atom index.
            IATOM = ZMCARDS(ICARD)%IORDER

            ! . Get the internal coordinates.
            R     = ZMCARDS(ICARD)%BOND
            THETA = ZMCARDS(ICARD)%ANGLE
            PHI   = ZMCARDS(ICARD)%DIHEDRAL

            ! . Write out the card data.
            WRITE ( UNIT, "(I6,3X,A8,3X,I4,2X,3F12.4,3I6)" ) ICARD, ATMNAM(IATOM), ATMNUM(IATOM), &
                                                             R, THETA, PHI, JCARD, KCARD, LCARD

         END DO

         ! . Write out a comment separator.
         IF ( IRES /= SUBIND(ISUB+1) ) THEN
            WRITE ( UNIT, "('!',79('-'))" )
         END IF
      END DO
   END DO

   ! . Terminate the file.
   WRITE ( UNIT, "('!',79('='))" )

   ! . Finish up.
   IF ( UNIT == OUTPUT ) THEN
      CALL PRINT_NOFORMAT_STOP
   ELSE
      CLOSE ( UNIT )
   END IF

   END SUBROUTINE ZMATRIX_WRITE

END MODULE ZMATRIX_IO
