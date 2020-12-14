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
!                              The Symmetry Module
!===============================================================================
!
! . Parameter data:
!
!   BOX_NAME_LENGTH                 The length of the box type name.
!
! . Scalar data:
!
!   BOX_TYPE                        The type of periodic box.
!   QBOX                            The flag to indicate if a box exists.
!
! . Array data:
!
!   BOXL                            The sides of the rectangular box.
!
! . Subroutines:
!
!   SYMMETRY_CUBIC_BOX              Define a cubic periodic box.
!   SYMMETRY_INITIALIZE             Initialize the symmetry variables.
!   SYMMETRY_ORTHORHOMBIC_BOX       Define an orthorhombic periodic box.
!   SYMMETRY_RECORD_READ            Read a symmetry record.
!   SYMMETRY_RECORD_WRITE           Write a symmetry record.
!   SYMMETRY_SUMMARY                Write out a summary of the symmetry data.
!
! . Notes:
!
!   The symmetry module defines the symmetry of a system. In this version
!   there are only a limited number of options. It is only possible to use
!   translational symmetry and, at present, only two possibilities are
!   allowed - cubic and orthorhombic systems. These are particularly useful
!   for generating periodic systems for use with periodic boundary conditions.
!
!===============================================================================
MODULE SYMMETRY

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE PARSING
USE PRINTING,    ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, &
                        PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, PRINT_SUMMARY_STOP

IMPLICIT NONE
PUBLIC
#ifndef PGPC
SAVE
#endif

! . Parameter data.
INTEGER, PARAMETER :: BOX_NAME_LENGTH = 16

! . Scalar data.
CHARACTER ( LEN = BOX_NAME_LENGTH ) :: BOX_TYPE = " "
LOGICAL                             :: QBOX     = .FALSE.

! . Array data.
REAL ( KIND = DP ), DIMENSION(1:3) :: BOXL = 0.0_DP

!==============================================================================
CONTAINS
!==============================================================================

   !----------------------------------
   SUBROUTINE SYMMETRY_CUBIC_BOX ( A )
   !----------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: A

   ! . Initialize the symmetry variables.
   CALL SYMMETRY_INITIALIZE

   ! . Set the box parameters.
   BOXL = ABS ( A ) ; BOX_TYPE = "CUBIC" ; QBOX = .TRUE.

   ! . Write out the symmetry data.
   CALL SYMMETRY_SUMMARY

   END SUBROUTINE SYMMETRY_CUBIC_BOX

   !-----------------------------
   SUBROUTINE SYMMETRY_INITIALIZE
   !-----------------------------

   ! . Reset the box parameters.
   BOXL = 0.0_DP ; BOX_TYPE = " " ; QBOX = .FALSE.

   END SUBROUTINE SYMMETRY_INITIALIZE

   !-----------------------------------------------
   SUBROUTINE SYMMETRY_ORTHORHOMBIC_BOX ( A, B, C )
   !-----------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: A, B, C


   ! . Initialize the symmetry variables.
   CALL SYMMETRY_INITIALIZE

   ! . Set the box parameters.
   BOXL(1) = ABS ( A ) ; BOXL(2) = ABS ( B ) ; BOXL(3) = ABS ( C ) ; BOX_TYPE = "ORTHORHOMBIC" ; QBOX = .TRUE.

   ! . Write out the data.
   CALL SYMMETRY_SUMMARY

   END SUBROUTINE SYMMETRY_ORTHORHOMBIC_BOX

   !------------------------------------------
   SUBROUTINE SYMMETRY_RECORD_READ ( COMPACT )
   !------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: COMPACT

   ! . Local scalars.
   CHARACTER ( LEN = LEN ( BOX_TYPE ) ) :: TMPBOX
   INTEGER                              :: I, ICARD
   LOGICAL                              :: QCOMPACT

   ! . Check the format.
   IF ( PRESENT ( COMPACT ) ) THEN
      QCOMPACT = COMPACT
   ELSE
      QCOMPACT = .FALSE.
   END IF

   ! . Get the first word of the next line.
   IF ( QCOMPACT ) QADVANCE = .FALSE.
   CALL GET_LINE
   IF ( QCOMPACT ) QADVANCE = .TRUE.
   CALL GET_WORD

   ! . Check for symmetry information.
   IF ( WORD(1:WRDLEN) == "SYMMETRY" ) THEN

      ! . Full format.
      IF ( .NOT. QCOMPACT ) THEN

         ! . Get the rest of the information on the line.
         CALL GET_INTEGER ( ICARD )

         ! . Check the number of cards.
         IF ( ICARD /= 1 ) CALL PARSE_ERROR ( "SYMMETRY_RECORD_READ", "Invalid number of symmetry cards." )

         ! . Get the first word on the next line.
         CALL GET_LINE

      END IF

      ! . Interpret the next word on the line.
      CALL GET_WORD ; TMPBOX = WORD(1:WRDLEN)

      ! . Branch on the box type.
      SELECT CASE ( TMPBOX )
      CASE ( "CUBIC           " ) ; CALL GET_REAL ( BOXL(1) ) ; BOXL(2:3) = BOXL(1)
      CASE ( "ORTHORHOMBIC    " ) ; DO I = 1,3 ; CALL GET_REAL ( BOXL(I) ) ; END DO
      CASE DEFAULT ; CALL PRINT_ERROR ( "SYMMETRY_RECORD_READ", "Unknown periodic box type." )
      END SELECT

      ! . Set some symmetry variables.
      BOX_TYPE = TMPBOX
      QBOX     = .TRUE.

      ! . Reposition the file correctly by getting the first word on the next line.
      IF ( .NOT. QCOMPACT ) THEN
         CALL GET_LINE
         CALL GET_WORD
      END IF

      ! . Print a summary of the symmetry data.
      CALL SYMMETRY_SUMMARY

   END IF

   END SUBROUTINE SYMMETRY_RECORD_READ

   !-------------------------------------------------
   SUBROUTINE SYMMETRY_RECORD_WRITE ( UNIT, COMPACT )
   !-------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN)           :: UNIT
   LOGICAL, INTENT(IN), OPTIONAL :: COMPACT

   ! . Local scalars.
   INTEGER :: ICARD
   LOGICAL :: QCOMPACT

   ! . Check the format.
   IF ( PRESENT ( COMPACT ) ) THEN
      QCOMPACT = COMPACT
   ELSE
      QCOMPACT = .FALSE.
   END IF

   ! . Compact format.
   IF ( QCOMPACT ) THEN

      ! . There is symmetry data.
      IF ( QBOX ) THEN

         ! . Write out the periodic box information.
         SELECT CASE ( BOX_TYPE )
         CASE ( "CUBIC           " ) ; WRITE ( UNIT, "(A,F18.10)"  ) "SYMMETRY "//BOX_TYPE, BOXL(1)
         CASE ( "ORTHORHOMBIC    " ) ; WRITE ( UNIT, "(A,3F18.10)" ) "SYMMETRY "//BOX_TYPE, BOXL(1:3)
         CASE DEFAULT ; CALL PRINT_ERROR ( "SYMMETRY_RECORD_WRITE", "Unknown periodic box type." )
         END SELECT

      ! . There is no symmetry data so write a blank line.
      ELSE
         WRITE ( UNIT, "(1X)" )
      END IF

   ! . Full format.
   ELSE

      ! . There is symmetry data.
      IF ( QBOX ) THEN

         ! . Initialize the number of cards.
         ICARD = 1

         ! . Write out the symmetry header.
         WRITE ( UNIT, "('!',79('='))" )
         WRITE ( UNIT, "(A,I6)" ) "Symmetry", ICARD

         ! . Write out the periodic box information.
         SELECT CASE ( BOX_TYPE )
         CASE ( "CUBIC           " ) ; WRITE ( UNIT, "(A,F18.10)"  ) BOX_TYPE, BOXL(1)
         CASE ( "ORTHORHOMBIC    " ) ; WRITE ( UNIT, "(A,3F18.10)" ) BOX_TYPE, BOXL(1:3)
         CASE DEFAULT ; CALL PRINT_ERROR ( "SYMMETRY_RECORD_WRITE", "Unknown periodic box type." )
         END SELECT

      END IF
   END IF

   END SUBROUTINE SYMMETRY_RECORD_WRITE

   !--------------------------
   SUBROUTINE SYMMETRY_SUMMARY
   !--------------------------

   ! . Local scalars.
   REAL ( KIND = DP ) :: VOLUME

   ! . Check that there is symmetry.
   IF ( QBOX ) THEN

      ! . Calculate the volume.
      VOLUME = PRODUCT ( BOXL )

      ! . Branch on the box type.
      SELECT CASE ( BOX_TYPE )
      CASE ( "CUBIC           " )
         CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF5500" )
	 CALL PRINT_SUMMARY_START ( "Summary of Cubic Box Data" )
	 WRITE ( PRINT_LINE, "(F14.4)" ) BOXL(1) ; CALL PRINT_SUMMARY_ELEMENT ( "Box Length" )
	 WRITE ( PRINT_LINE, "(F14.4)" ) VOLUME  ; CALL PRINT_SUMMARY_ELEMENT ( "Box Volume" )
      CASE ( "ORTHORHOMBIC    " )
         CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF5500" )
	 CALL PRINT_SUMMARY_START ( "Summary of Orthorhombic Box Data" )
	 WRITE ( PRINT_LINE, "(F14.4)" ) BOXL(1) ; CALL PRINT_SUMMARY_ELEMENT ( "Box Length - X" )
	 WRITE ( PRINT_LINE, "(F14.4)" ) BOXL(2) ; CALL PRINT_SUMMARY_ELEMENT ( "Box Length - Y" )
	 WRITE ( PRINT_LINE, "(F14.4)" ) BOXL(3) ; CALL PRINT_SUMMARY_ELEMENT ( "Box Length - Z" )
	 WRITE ( PRINT_LINE, "(F14.4)" ) VOLUME  ; CALL PRINT_SUMMARY_ELEMENT ( "Box Volume" )
      END SELECT

      ! . Write out the terminator.
      CALL PRINT_SUMMARY_STOP

   ! . There is no symmetry.
   ELSE

      CALL PRINT_PARAGRAPH ( TEXT = "No symmetry is defined." )

   END IF

   END SUBROUTINE SYMMETRY_SUMMARY

END MODULE SYMMETRY
