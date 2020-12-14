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
!                         The Character Handling Module
!===============================================================================
!
! . Functions:
!
!   DECODE_INTEGER                  Decode a string into an integer.
!   DECODE_LOGICAL                  Decode a string into a logical.
!   DECODE_REAL                     Decode a string into a real.
!   TO_UPPER_CASE                   Convert a string to upper case.
!
! . Subroutines:
!
!   ENCODE_INTEGER                  Encode an integer in a string.
!
!===============================================================================
MODULE STRING

! . Module declarations.
USE DEFINITIONS, ONLY : DP, LINE_LENGTH
USE PRINTING,    ONLY : PRINT_ERROR

IMPLICIT NONE
PUBLIC
PRIVATE :: WORK, WRKLEN

! . Parameter data.
CHARACTER ( LEN = 5 ), PARAMETER :: DEFAULT_INTEGER_FORMAT = "(I20)"

! . Character data.
CHARACTER ( LEN = LINE_LENGTH ) :: WORK

! . Integer data.
INTEGER :: WRKLEN = 0

!===============================================================================
CONTAINS
!===============================================================================

   !---------------------------------
   FUNCTION DECODE_INTEGER ( STRING )
   !---------------------------------

   ! . Function declarations.
   INTEGER :: DECODE_INTEGER

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: STRING

   ! . Local scalars.
   CHARACTER ( LEN = 6 ) :: FMT
   INTEGER               :: IOSTAT, LENGTH, NUMBER

   ! . Initialize the integer.
   NUMBER = 0

   ! . Find the length of the string.
   LENGTH = LEN ( STRING )

   ! . The string has a non-zero length.
   IF ( LENGTH > 0 ) THEN

      ! . Write the format to the character string FMT.
      WRITE ( FMT, "(A,I2,A)", IOSTAT = IOSTAT ) "(I", LENGTH, ")"

      ! . Perform the conversion.
      IF ( IOSTAT == 0 ) READ ( STRING(1:LENGTH), FMT, IOSTAT = IOSTAT ) NUMBER

      ! . Check for an I/O error.
      IF ( IOSTAT /= 0 ) THEN
         CALL PRINT_ERROR ( "DECODE_INTEGER", "I/O Error.", IOSTAT )
         NUMBER = HUGE ( 0 )
      END IF

   END IF

   ! . Assign a value to the function.
   DECODE_INTEGER = NUMBER

   END FUNCTION DECODE_INTEGER

   !---------------------------------
   FUNCTION DECODE_LOGICAL ( STRING )
   !---------------------------------

   ! . Function declarations.
   LOGICAL :: DECODE_LOGICAL

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: STRING

   ! . Local scalars.
   LOGICAL :: QFLAG

   ! . Initialize the variable.
   QFLAG = .FALSE.

   ! . Branch on the string.
   SELECT CASE ( TRIM ( STRING ) )
   CASE ( "F"       ) ; QFLAG = .FALSE.
   CASE ( ".F."     ) ; QFLAG = .FALSE.
   CASE ( ".FALSE." ) ; QFLAG = .FALSE.
   CASE ( "T"       ) ; QFLAG = .TRUE.
   CASE ( ".T."     ) ; QFLAG = .TRUE.
   CASE ( ".TRUE."  ) ; QFLAG = .TRUE.
   CASE DEFAULT ; CALL PRINT_ERROR ( "DECODE_LOGICAL", "Unable to decode string to a logical." )
   END SELECT

   ! . Set the function value.
   DECODE_LOGICAL = QFLAG

   END FUNCTION DECODE_LOGICAL

   !------------------------------
   FUNCTION DECODE_REAL ( STRING )
   !------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: DECODE_REAL

   ! . Argument declarations.
   CHARACTER ( LEN = * ), INTENT(IN) :: STRING

   ! . Local scalars.
   CHARACTER ( LEN = 7 ) :: FMT
   INTEGER               :: EXPNT, IEXPNT, IOSTAT, LENGTH
   REAL ( KIND = DP )    :: EXPFAC, NUMBER, VALUE

   ! . Initialize the real.
   NUMBER = 0.0_DP

   ! . Get the length of the string.
   LENGTH = LEN ( STRING )

   ! . The string has a non-zero length.
   IF ( LENGTH > 0 ) THEN

      ! . Check for an exponent.
      IEXPNT = INDEX ( TO_UPPER_CASE ( STRING(1:LENGTH) ), "E" )

      ! . There is an exponent.
      IF ( IEXPNT > 0 ) THEN

         ! . Calculate the exponent factor.
         EXPNT  = DECODE_INTEGER ( STRING((IEXPNT+1):LENGTH) )
         EXPFAC = 10.0_DP ** EXPNT

         ! . Reset the length of the real number string.
         LENGTH = IEXPNT - 1

      ! . There is no exponent.
      ELSE
         EXPFAC = 1.0_DP
      END IF

      ! . Write the format to the character string FMT.
      WRITE ( FMT, "(A,I2,A)", IOSTAT = IOSTAT ) "(F", LENGTH, ".0)"

      ! . Perform the conversion.
      IF ( IOSTAT == 0 ) READ ( STRING(1:LENGTH), FMT, IOSTAT = IOSTAT ) VALUE
      NUMBER = EXPFAC * VALUE

      ! . Check for an I/O error.
      IF ( IOSTAT /= 0 ) THEN
         CALL PRINT_ERROR ( "DECODE_REAL", "I/O Error.", IOSTAT )
         NUMBER = HUGE ( 0.0_DP )
      END IF

   END IF

   ! . Assign a value to the function.
   DECODE_REAL = NUMBER

   END FUNCTION DECODE_REAL

   !---------------------------------------------------
   SUBROUTINE ENCODE_INTEGER ( NUMBER, STRING, FORMAT )
   !---------------------------------------------------

   ! . Argument declarations.
   CHARACTER ( LEN = * ), INTENT(OUT)          :: STRING
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FORMAT
   INTEGER,               INTENT(IN)           :: NUMBER

   ! . Write the integer to the work string.
   IF ( PRESENT ( FORMAT ) ) THEN
      WRITE ( WORK, FORMAT ) NUMBER
   ELSE
      WRITE ( WORK, DEFAULT_INTEGER_FORMAT ) NUMBER
   END IF

   ! . Remove any leading or trailing blanks.
   WRKLEN = LEN ( WORK )
   WORK(1:WRKLEN) = ADJUSTL  ( WORK(1:WRKLEN) )
   WRKLEN         = LEN_TRIM ( WORK(1:WRKLEN) ) 

   ! . Copy the scratch string to the output string.
   IF ( WRKLEN <= LEN ( STRING ) ) THEN
      STRING = WORK(1:WRKLEN)
   ! . There is an error.
   ELSE
      CALL PRINT_ERROR ( "ENCODE_INTEGER", "String is too short to hold the encoded integer." )
      STRING = " "
   END IF

   END SUBROUTINE ENCODE_INTEGER

   !--------------------------------
   FUNCTION TO_UPPER_CASE ( STRING )
   !--------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: STRING

   ! . Function declarations.
   CHARACTER ( LEN = LEN ( STRING ) ) :: TO_UPPER_CASE

   ! . Local parameters.
   INTEGER, PARAMETER :: BIG_A = ICHAR ( "A" ), LITTLE_A = ICHAR ( "a" ), LITTLE_Z = ICHAR ( "z" )

   ! . Local scalars.
   INTEGER :: I, ICHR

   ! . Loop over the characters in the string.
   DO I = 1,LEN ( STRING )

      ! . Get the ASCII order for the character to be converted.
      ICHR = ICHAR ( STRING(I:I) )

      ! . Use the order to change the case of the character.
      IF ( ( ICHR >= LITTLE_A ) .AND. ( ICHR <= LITTLE_Z ) ) THEN
         TO_UPPER_CASE(I:I) = CHAR ( ICHR + BIG_A - LITTLE_A )
      ELSE
         TO_UPPER_CASE(I:I) = STRING(I:I)
      END IF
   END DO

   END FUNCTION TO_UPPER_CASE

END MODULE STRING
