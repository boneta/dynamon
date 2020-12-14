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
!                               The XYZ I/O Module
!===============================================================================
!
! . Public subroutines:
!
!   XYZ_DEFINE             Define a structure from a XYZ file.
!   XYZ_READ               Read the coordinates from a XYZ file.
!   XYZ_WRITE              Write a XYZ file.
!
! . Private subroutines:
!
!   XYZ_ATOM_NAME_PARSE    Parse an atom name in XYZ format.
!
! . Notes:
!
!   XYZ_DEFINE is used to define a system from a XYZ file. All existing atom
!   data is lost. XYZ_READ reads in a new set of coordinates but the system
!   composition (the number, type and the order of the atoms) must be the same
!   as that which exists in the module ATOMS. XYZ_WRITE writes out a XYZ file.
!
!   For all three of the above subroutines, the argument FILE gives the file
!   name to read from or write to. If it is not present the standard input or
!   output streams are used (whichever is appropriate). For XYZ_READ and
!   XYZ_WRITE there is also a second argument which is an optional array that
!   can contain the data to read into or write out from. This is useful when
!   dealing with other atom data (such as derivatives or velocities).
!
!   XYZ_ATOM_NAME recognizes pure element names, atomic numbers or names of
!   the form "Element_Name"//"Integer".
!
!   XYZ_DEFINE will allocate subsystem and residue information. All atoms are
!   placed in the same subsystem with the name "SYSTEM" and the same residue
!   with the name "RESIDUE". This allows an XYZ to be read and the system
!   coordinates to be printed in either DYNAMO or PDB format. Atom names are
!   assigned of the form "Element Symbol" + "Atom Sequence Number". Unknown
!   elements are given the symbol "XX".
!
!===============================================================================
MODULE XYZ_IO

! . Module declarations.
USE DEFINITIONS,  ONLY : DP
USE ELEMENTS,     ONLY : MASS
USE FILES,        ONLY : NEXT_UNIT
USE IO_UNITS,     ONLY : INPUT, OUTPUT
USE PARSING
USE PRINTING,     ONLY : PRINT_BLANKLINE, PRINT_ERROR, PRINT_NOFORMAT_START, PRINT_NOFORMAT_STOP, &
                         PRINT_PARAGRAPH
USE STRING,       ONLY : DECODE_INTEGER, ENCODE_INTEGER, TO_UPPER_CASE

USE ATOMS
USE SEQUENCE
USE SYMMETRY

IMPLICIT NONE
PRIVATE
PUBLIC :: XYZ_DEFINE, XYZ_READ, XYZ_WRITE

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------
   SUBROUTINE XYZ_DEFINE ( FILE, TITLE )
   !------------------------------------

   ! . Optional scalar arugments.
   CHARACTER ( LEN = * ), INTENT(IN),  OPTIONAL :: FILE
   CHARACTER ( LEN = * ), INTENT(OUT), OPTIONAL :: TITLE

   ! . Local scalars.
   CHARACTER ( LEN = 20 ) :: NUMSTR
   INTEGER                :: IATOM, IOSTAT, NATMIN, UNIT

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "XYZ_DEFINE", "I/O Error.", IOSTAT )

   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Initialize the atom, sequence and symmetry data structures.
   CALL ATOMS_INITIALIZE
   CALL SEQUENCE_INITIALIZE
   CALL SYMMETRY_INITIALIZE

   ! . Read in the number of atoms in the file.
   CALL GET_LINE
   CALL GET_INTEGER ( NATMIN )

   ! . Check the counters.
   IF ( NATMIN <= 0 ) CALL PARSE_ERROR ( "XYZ_DEFINE", "Invalid number of atoms." )

   ! . Allocate the atom and sequence data arrays.
   CALL ATOMS_ALLOCATE ( NATMIN )
   CALL SEQUENCE_ALLOCATE ( 1, 1 )

   ! . Set all the sequence information.
   RESIND(1) = 0 ; RESIND(2) = NATOMS ; RESNAM(1) = "RESIDUE"
   SUBIND(1) = 0 ; SUBIND(2) = 1      ; SUBNAM(1) = "SYSTEM"

   ! . Read in any symmetry information or a title.
   CALL SYMMETRY_RECORD_READ ( COMPACT = .TRUE. )
   IF ( PRESENT ( TITLE ) ) CALL COPY_LINE ( TITLE )

   ! . Loop over the number of atoms.
   DO IATOM = 1,NATOMS

      ! . Read in the next line and the first word.
      CALL GET_LINE
      CALL GET_WORD

      ! . Get the atomic number of the element.
      CALL XYZ_ATOM_NAME_PARSE ( ATMNUM(IATOM) )

      ! . Read in the coordinate data.
      CALL GET_REAL ( ATMCRD(1,IATOM) )
      CALL GET_REAL ( ATMCRD(2,IATOM) )
      CALL GET_REAL ( ATMCRD(3,IATOM) ) 

   END DO

   ! . Fill the atom mass and atom name arrays.
   DO IATOM = 1,NATOMS

      ! . Get the atom number as a string.
      CALL ENCODE_INTEGER ( IATOM, NUMSTR )

      ! . Fill the arrays.
      IF ( ( ATMNUM(IATOM) > 0 ) .AND. ( ATMNUM(IATOM) <= NELEMENTS ) ) THEN
         ATMMAS(IATOM) = MASS(ATMNUM(IATOM))
	 ATMNAM(IATOM) = TRIM ( TO_UPPER_CASE ( SYMBOL(ATMNUM(IATOM)) ) ) // TRIM ( NUMSTR )
      ELSE
         ATMMAS(IATOM) = 0.0_DP
	 ATMNAM(IATOM) = "XX" // TRIM ( NUMSTR )
      END IF

   END DO

   ! . Close the input file if necessary.
   IF ( UNIT /= INPUT ) CLOSE ( UNIT )

   ! . Reset the parsing unit.
   CALL POP_UNIT

   ! . Do some printing.
   IF ( PRESENT ( FILE ) ) THEN
      CALL PRINT_PARAGRAPH ( TEXT = "System read from " // FILE )
   ELSE
      CALL PRINT_PARAGRAPH ( TEXT = "System read from the input stream." )
   END IF

   ! . Print out a summary of the atoms, sequence and symmetry data structures.
   CALL ATOMS_SUMMARY
   CALL SEQUENCE_SUMMARY
   CALL SYMMETRY_SUMMARY

   END SUBROUTINE XYZ_DEFINE

   !----------------------------------------
   SUBROUTINE XYZ_READ ( FILE, DATA, TITLE )
   !----------------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN),  OPTIONAL :: FILE
   CHARACTER ( LEN = * ), INTENT(OUT), OPTIONAL :: TITLE

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(OUT), OPTIONAL :: DATA

   ! . Local scalars.
   CHARACTER ( LEN = 16 ) :: TAG
   INTEGER                :: IATOM, IOSTAT, LTAG, NATMIN, NUMBER, UNIT

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: ATMDAT

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . Check for an error.
      IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "XYZ_READ", "I/O Error.", IOSTAT )

   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Read in and check the number of atoms in the file.
   CALL GET_LINE
   CALL GET_INTEGER ( NATMIN )
   IF ( NATMIN > NATOMS ) CALL PARSE_ERROR ( "XYZ_READ", "Invalid number of atoms." )

   ! . Read in any symmetry information or a title.
   CALL SYMMETRY_RECORD_READ ( COMPACT = .TRUE. )
   IF ( PRESENT ( TITLE ) ) CALL COPY_LINE ( TITLE )

   ! . Loop over the number of atoms.
   DO IATOM = 1,NATOMS

      ! . Read in the next line and the first word.
      CALL GET_LINE
      CALL GET_WORD

      ! . Get the atomic number of the element.
      CALL XYZ_ATOM_NAME_PARSE ( NUMBER )

      ! . Check the atom number.
      IF ( NUMBER /= ATMNUM(IATOM) ) CALL PARSE_ERROR ( "XYZ_READ", "Atom number mismatch." )

      ! . Read in the coordinate data.
      CALL GET_REAL ( ATMDAT(1,IATOM) )
      CALL GET_REAL ( ATMDAT(2,IATOM) )
      CALL GET_REAL ( ATMDAT(3,IATOM) ) 

   END DO

   ! . Close the input file if necessary.
   IF ( UNIT /= INPUT ) CLOSE ( UNIT )

   ! . Reset the parsing unit.
   CALL POP_UNIT

   ! . Put data in DATA if the DATA argument is present.
   IF ( PRESENT ( DATA ) ) THEN
      DATA = ATMDAT
      TAG  = "Data"
   ! . By default put the data in ATMCRD.
   ELSE
      ! . Rellocate and copy ATMDAT to the ATMCRD array.
      IF ( ALLOCATED ( ATMCRD ) ) DEALLOCATE ( ATMCRD )
      ALLOCATE ( ATMCRD(1:3,1:NATOMS) )
      ATMCRD = ATMDAT
      TAG    = "Coordinates"
   END IF

   ! . Get the length of the tag.
   LTAG = LEN_TRIM ( TAG )

   ! . Do some printing.
   IF ( PRESENT ( FILE ) ) THEN
      CALL PRINT_PARAGRAPH ( TEXT = TAG(1:LTAG)//" read from " // FILE )
   ELSE
      CALL PRINT_PARAGRAPH ( TEXT = TAG(1:LTAG)//" read from the input stream." )
   END IF

   END SUBROUTINE XYZ_READ

   !----------------------------------------------------
   SUBROUTINE XYZ_WRITE ( FILE, DATA, SELECTION, TITLE )
   !----------------------------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE, TITLE

   ! . Optional array arguments.
   LOGICAL,            DIMENSION(1:NATOMS),     INTENT(IN), OPTIONAL :: SELECTION
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(IN), OPTIONAL :: DATA

   ! . Local scalars.
   CHARACTER ( LEN = 16 ) :: TAG
   INTEGER                :: IATOM, IOSTAT, LTAG, NATM, UNIT

   ! . Local arrays.
   LOGICAL,            DIMENSION(1:NATOMS)     :: QATOM
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: ATMDAT

   ! . Write out DATA if the DATA argument is present.
   IF ( PRESENT ( DATA ) ) THEN
      ATMDAT = DATA
      TAG    = "Data"
   ! . By default write out the contents of ATMCRD.
   ELSE
      ATMDAT = ATMCRD
      TAG    = "Coordinates"
   END IF

   ! . Get the length of the tag.
   LTAG = LEN_TRIM ( TAG )

   ! . Check to see if the selection array is present.
   IF ( PRESENT ( SELECTION ) ) THEN

      ! . Save the SELECTION array in QATOM.
      QATOM = SELECTION

   ! . All data is to be printed.
   ELSE
      QATOM = .TRUE.
   END IF

   ! . Get the number of selected atoms.
   NATM = COUNT ( QATOM )

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "XYZ_WRITE", "I/O Error.", IOSTAT )

      ! . Do some printing.
      CALL PRINT_PARAGRAPH ( TEXT = TAG(1:LTAG)//" written to " // FILE )

      ! . Write out the number of atoms to the file.
      WRITE ( UNIT, "(I6)" ) NATM

   ! . The unit for writing is the output stream.
   ELSE

      ! . Assign the unit number.
      UNIT = OUTPUT

      ! . Do some printing.
      CALL PRINT_PARAGRAPH ( TEXT = TAG(1:LTAG)//" written to the output stream." )
      CALL PRINT_BLANKLINE
      CALL PRINT_NOFORMAT_START

      ! . Write out the separator.
      WRITE ( UNIT, "('!',69('='))" )

      ! . Write out the number of atoms.
      WRITE ( UNIT, "(I6,A)" ) NATM, " ! # of atoms."

   END IF

   ! . Write out a symmetry definition.
   IF ( PRESENT ( TITLE ) ) THEN
      WRITE ( UNIT, "(A)" ) TITLE
   ELSE
      CALL SYMMETRY_RECORD_WRITE ( UNIT, COMPACT = .TRUE. )
   END IF

   ! . Write out the selected atoms.
   DO IATOM = 1,NATOMS
      IF ( QATOM(IATOM) ) THEN
         IF ( ( ATMNUM(IATOM) > 0 ) .AND. ( ATMNUM(IATOM) <= NELEMENTS ) ) THEN
            WRITE ( UNIT, "(A,8X,3F20.12)" ) SYMBOL(ATMNUM(IATOM)), ATMDAT(1:3,IATOM)
         ELSE
            WRITE ( UNIT, "(I10,3F20.12)" ) ATMNUM(IATOM), ATMDAT(1:3,IATOM)
         END IF
      END IF
   END DO

   ! . Finish off.
   IF ( UNIT == OUTPUT ) THEN

      ! . Write out the terminator.
      WRITE ( UNIT, "('!',69('='))" )

      ! . Finish noformatting.
      CALL PRINT_NOFORMAT_STOP

   ! . Close the file.
   ELSE
      CLOSE ( UNIT )
   END IF

   END SUBROUTINE XYZ_WRITE

!------------------------------------------------------------------------------
! . Private subroutines.
!------------------------------------------------------------------------------

   !----------------------------------------
   SUBROUTINE XYZ_ATOM_NAME_PARSE ( NUMBER )
   !----------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(OUT) :: NUMBER

   ! . Local parameters.
   INTEGER, PARAMETER :: ZERO = ICHAR ( "0" ), NINE = ICHAR ( "9" )

   ! . Local scalars.
   INTEGER :: I, ICHR, NC

   ! . Initialization.
   NUMBER = -1

   ! . Check WRDLEN.
   IF ( WRDLEN <= 0 ) RETURN

   ! . Find out how many non-integers there are at the beginning of the word (maximum of two).
   NC = 0
   DO I = 1,2
      ICHR = ICHAR ( WORD(I:I) )
      IF ( ( ICHR >= ZERO ) .AND. ( ICHR <= NINE ) ) EXIT
      NC = NC + 1
   END DO

   ! . If there are characters check for a valid element name.
   IF ( NC > 0 ) THEN

      ! . Loop over the element names.
      DO I = 1,NELEMENTS
         IF ( WORD(1:NC) == TO_UPPER_CASE ( SYMBOL(I) ) ) THEN
            NUMBER = I ; EXIT
         END IF
      END DO

   ! . The word is an integer.
   ELSE
      NUMBER = DECODE_INTEGER ( WORD(1:WRDLEN) )
   END IF

   END SUBROUTINE XYZ_ATOM_NAME_PARSE

END MODULE XYZ_IO

