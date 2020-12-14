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
!                              The Printing Module
!===============================================================================
!
! . Scalars:
!
!   PRINT_LINE                   The line for printing.
!   PRINT_LINE_LENGTH            The length of the line.
!
!   QHTML                        The HTML flag.
!   QNOFORMAT                    The no-formatting flag.
!   QPARAGRAPH                   The paragraph flag.
!   QPRINT                       The flag to indicate that printing is on.
!   QSUMMARY                     The summary flag.
!   QTABLE                       The table flag.
!
! . Subroutines:
!
!   PRINT_ERROR                  Print an error.
!
!   PRINT_BLANKLINE              Print a blank line.
!   PRINT_HEADING                Print a heading.
!   PRINT_LINE_INITIALIZE        Initialize PRINT_LINE.
!   PRINT_LINEBREAK              Print a line break.   
!
!   PRINT_MATRIX                 Print a matrix stored in full form.
!   PRINT_SYMMETRIC_MATRIX       Print a symmetric matrix (upper triangular).
!
!   PRINT_PARAGRAPH              Print an entire paragraph.
!   PRINT_PARAGRAPH_START        Open a paragraph.
!   PRINT_PARAGRAPH_STOP         Close a paragraph.
!
!   PRINT_START                  Start the printing.
!   PRINT_STOP                   Stop the printing.
!   PRINT_STYLE                  Set the print style.
!
!   PRINT_SUMMARY_ELEMENT        Process an element in a summary.
!   PRINT_SUMMARY_ENDLINE        Finish a row with blank entries.
!   PRINT_SUMMARY_INITIALIZE     Initialize the summary data.
!   PRINT_SUMMARY_OPTIONS        Set the summary options.
!   PRINT_SUMMARY_START          Open a summary.
!   PRINT_SUMMARY_STOP           Close a summary.
!
!   PRINT_TABLE_ELEMENT          Process an element in a table.
!   PRINT_TABLE_ENDLINE          Finish a row with blank entries.
!   PRINT_TABLE_INITIALIZE       Initialize the table data.
!   PRINT_TABLE_OPTIONS          Set the table options.
!   PRINT_TABLE_START            Open a table.
!   PRINT_TABLE_STOP             Close a table.
!
!   PRINT_TEXT                   Print a text string.
!
! . Notes:
!
!   The allowed text styles are "HTML" or "TEXT". The latter is the default.
!
!===============================================================================
MODULE PRINTING

! . Module declarations.
USE DEFINITIONS, ONLY : DP, LINE_LENGTH
USE IO_UNITS,    ONLY : OUTPUT

IMPLICIT NONE
PRIVATE
PUBLIC :: PRINT_LINE, &
          PRINT_ERROR, PRINT_HEADING, PRINT_BLANKLINE, PRINT_LINEBREAK, PRINT_MATRIX,        &
	  PRINT_NOFORMAT_START, PRINT_NOFORMAT_STOP, PRINT_PARAGRAPH, PRINT_PARAGRAPH_START, &
	  PRINT_PARAGRAPH_STOP, PRINT_START, PRINT_STOP, PRINT_STYLE, PRINT_SUMMARY_ELEMENT, &
	  PRINT_SUMMARY_ENDLINE, PRINT_SUMMARY_INITIALIZE, PRINT_SUMMARY_OPTIONS,            &
	  PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_SYMMETRIC_MATRIX,                   &
          PRINT_TABLE_ELEMENT, PRINT_TABLE_ENDLINE, PRINT_TABLE_INITIALIZE,                  &
	  PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP, PRINT_TEXT
#ifndef PGPC
SAVE
#endif

! . HTML parameters.
CHARACTER ( LEN = 7 ), PARAMETER :: BODY_ALINK                 = "#00FF00", &
                                    BODY_BGCOLOR               = "#FFFFED", &
                                    BODY_LINK                  = "#0000FF", &
                                    BODY_TEXT                  = "#000000", &
                                    BODY_VLINK                 = "#800080", &
                                    SUMMARY_HEADER_COLOR       = "#0099FF", &
				    SUMMARY_TAG_COLOR          = "#AAAAAA", &
				    SUMMARY_TEXT_COLOR         = "#FFFFF7", &
				    TABLE_ENTRY_COLOR_DEFAULT  = "#FFFFF7", &
				    TABLE_HEADER_COLOR_DEFAULT = "#FF9900"

! . Text parameters.
INTEGER, PARAMETER :: SUMMARY_PAGEWIDTH_DEFAULT = 80, SUMMARY_VARIABLEWIDTH_DEFAULT = 14, &
                      TABLE_PAGEWIDTH_DEFAULT   = 80, TABLE_NCOLUMNS_DEFAULT        = 2

! . Scalars.
CHARACTER ( LEN = LINE_LENGTH ) :: PRINT_LINE = " "
INTEGER                         :: PRINT_LINE_LENGTH = 0
LOGICAL                         :: QHTML  = .FALSE., QNOFORMAT = .FALSE., QPARAGRAPH = .FALSE., &
                                   QPRINT = .FALSE., QSUMMARY  = .FALSE., QTABLE     = .FALSE.

! . Summary variables.
CHARACTER ( LEN =           7 ) :: SUMMARY_COLOR1 = SUMMARY_HEADER_COLOR, &
                                   SUMMARY_COLOR2 = SUMMARY_TAG_COLOR,    &
				   SUMMARY_COLOR3 = SUMMARY_TEXT_COLOR
CHARACTER ( LEN = LINE_LENGTH ) :: SUMMARY_LINE = " "
INTEGER                         :: NSUMMARY = 0, SUMMARY_LINE_LENGTH = 0, &
                                   SUMMARY_PAGEWIDTH     =   SUMMARY_PAGEWIDTH_DEFAULT,                                           &
                                   SUMMARY_TAGWIDTH      = ( SUMMARY_PAGEWIDTH_DEFAULT - 2 ) / 2 - SUMMARY_VARIABLEWIDTH_DEFAULT, &
				   SUMMARY_VARIABLEWIDTH =   SUMMARY_VARIABLEWIDTH_DEFAULT

! . Table variables.
CHARACTER ( LEN =           7 )    :: TABLE_ENTRY_COLOR = TABLE_ENTRY_COLOR_DEFAULT, TABLE_HEADER_COLOR = TABLE_HEADER_COLOR_DEFAULT
CHARACTER ( LEN = LINE_LENGTH )    :: TABLE_LINE = " "
INTEGER                            :: NCOLUMNS = TABLE_NCOLUMNS_DEFAULT, NTABLE = 0, TABLE_LINE_LENGTH = 0, &
                                      TABLE_PAGEWIDTH = TABLE_PAGEWIDTH_DEFAULT
INTEGER, ALLOCATABLE, DIMENSION(:) :: TABLE_VARIABLEWIDTHS

!===============================================================================
CONTAINS
!===============================================================================

!-------------------------------------------------------------------------------
! . Errors.
!-------------------------------------------------------------------------------

   !------------------------------------------------
   SUBROUTINE PRINT_ERROR ( ROUTINE, MESSAGE, CODE )
   !------------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: MESSAGE, ROUTINE

   ! . Optional scalar arguments.
   INTEGER, INTENT(IN), OPTIONAL :: CODE

   ! . Local scalars.
   INTEGER :: LENM, LENR

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check here for any unfinished environments, etc.
   IF ( QNOFORMAT  ) CALL PRINT_NOFORMAT_STOP
   IF ( QPARAGRAPH ) CALL PRINT_PARAGRAPH_STOP
   IF ( QSUMMARY   ) CALL PRINT_SUMMARY_STOP
   IF ( QTABLE     ) CALL PRINT_TABLE_STOP

   ! . Find the lengths of the character strings.
   LENM = LEN_TRIM ( MESSAGE )
   LENR = LEN_TRIM ( ROUTINE )

   WRITE ( OUTPUT, "(/A)" ) "Error in " // ROUTINE(1:LENR) // ":" // MESSAGE(1:LENM)
   IF ( PRESENT ( CODE ) ) WRITE ( OUTPUT, "(A,I6)" ) "Error code = ", CODE
   WRITE ( OUTPUT, "(/A)" ) "Program terminating."

   ! . Finish printing.
   CALL PRINT_STOP

   ! . Stop execution of the program.
   STOP

   END SUBROUTINE PRINT_ERROR

!-------------------------------------------------------------------------------
! . General printing.
!-------------------------------------------------------------------------------

   !-------------------------
   SUBROUTINE PRINT_BLANKLINE
   !-------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Print the blank line.
   WRITE ( OUTPUT, "(1X)" )

   END SUBROUTINE PRINT_BLANKLINE

   !------------------------------------------------------------------------------------
   SUBROUTINE PRINT_HEADING ( COLOR, TAG, RULE_SIZE, RULE_WIDTH, PAGEWIDTH, QBLANKLINE )
   !------------------------------------------------------------------------------------

   ! . HTML arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: COLOR, TAG
   INTEGER,               INTENT(IN) :: RULE_SIZE, RULE_WIDTH

   ! . Text arguments.
   INTEGER, INTENT(IN) :: PAGEWIDTH
   LOGICAL, INTENT(IN) :: QBLANKLINE

   ! . Local scalars.
   INTEGER :: N, NSPACES

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Get and check the print line length.
   PRINT_LINE_LENGTH = LEN_TRIM ( PRINT_LINE ) ; IF ( PRINT_LINE_LENGTH < 0 ) RETURN

      ! . Determine the spacing.
      N       = MIN ( PRINT_LINE_LENGTH, PAGEWIDTH )
      NSPACES = ( PAGEWIDTH - N + 1 ) / 2

      ! . Write out the header.
      IF ( QBLANKLINE ) WRITE ( OUTPUT, "(1X)" )
      WRITE ( OUTPUT, "(A)" ) REPEAT ( "-", PAGEWIDTH )
      WRITE ( OUTPUT, "(A)" ) REPEAT ( " ", NSPACES   ) // PRINT_LINE(1:N)
      WRITE ( OUTPUT, "(A)" ) REPEAT ( "-", PAGEWIDTH )

   ! . Initialize PRINT_LINE.
   CALL PRINT_LINE_INITIALIZE

   END SUBROUTINE PRINT_HEADING

   !-------------------------
   SUBROUTINE PRINT_LINEBREAK
   !-------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   END SUBROUTINE PRINT_LINEBREAK

   !-------------------------------
   SUBROUTINE PRINT_LINE_INITIALIZE
   !-------------------------------

   PRINT_LINE = REPEAT ( " ", LINE_LENGTH ) ; PRINT_LINE_LENGTH = 0

   END SUBROUTINE PRINT_LINE_INITIALIZE

   !---------------------------------------------------
   SUBROUTINE PRINT_TEXT ( TEXT, BOLD, CENTER, ITALIC )
   !---------------------------------------------------

   ! . Arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: TEXT
   LOGICAL,               INTENT(IN), OPTIONAL :: BOLD, CENTER, ITALIC

   ! . Local scalars.
   LOGICAL :: QBOLD, QCENTER, QITALIC

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . If TEXT is present copy it to PRINT_LINE
   IF ( PRESENT ( TEXT ) ) PRINT_LINE = TEXT

   ! . Get and check the print line length.
   PRINT_LINE_LENGTH = LEN_TRIM ( PRINT_LINE ) ; IF ( PRINT_LINE_LENGTH < 0 ) RETURN

   ! . Write out the text.
   WRITE ( OUTPUT, "(A)" ) PRINT_LINE(1:PRINT_LINE_LENGTH)

   ! . Initialize PRINT_LINE.
   CALL PRINT_LINE_INITIALIZE

   END SUBROUTINE PRINT_TEXT

!-------------------------------------------------------------------------------
! . Matrix printing.
!-------------------------------------------------------------------------------

   !----------------------------------------
   SUBROUTINE PRINT_MATRIX ( MATRIX, TITLE )
   !----------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: TITLE

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: MATRIX

   ! . Local scalars.
   INTEGER :: COLINC, I, ICOL, IROW, NCOL, NROW

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Get the dimensions of the matrix.
   NCOL = SIZE ( MATRIX, 2 )
   NROW = SIZE ( MATRIX, 1 )

   ! . Check that the matrix has a non-zero size.
   IF ( ( NCOL == 0 ) .OR. ( NROW == 0 ) ) RETURN

   ! . Print out the title.
   IF ( PRESENT ( TITLE ) ) THEN
      WRITE ( OUTPUT, "(/A)" ) TITLE
   ELSE
      WRITE ( OUTPUT, "(/)" )
   END IF

   ! . Initialization.
   ICOL = 0

   ! . Top of the loop over rows.
   10 COLINC = MIN ( 12, ( NCOL - ICOL ) )

   ! . Print out the column numbers.
   WRITE ( OUTPUT, '(12(1X,I6,4X))' ) ( ICOL+I, I = 1,COLINC )

   ! . Print out the matrix.
   DO IROW = 1,NROW
      WRITE ( OUTPUT, "(12(1X,F10.4))" ) MATRIX(IROW,ICOL+1:ICOL+COLINC)
   END DO

   ! . Increment the column number.
   ICOL = ICOL + COLINC

   ! . Check to see if more columns need to be printed.
   IF ( ICOL < NCOL ) GO TO 10

   END SUBROUTINE PRINT_MATRIX

   !--------------------------------------------------
   SUBROUTINE PRINT_SYMMETRIC_MATRIX ( MATRIX, TITLE )
   !--------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: TITLE

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: MATRIX

   ! . Local scalars.
   INTEGER :: COLINC, I, ICOL, II, IROW, N

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Get the dimension of the matrix.
   N = ( NINT ( SQRT ( REAL ( ( 8 * SIZE ( MATRIX ) + 1 ), DP ) ) ) - 1 ) / 2

   ! . Check that the matrix has the correct size.
   IF ( ( N <= 0 ) .OR. ( ( ( N * ( N + 1 ) ) / 2 ) /= SIZE ( MATRIX ) ) ) RETURN

   ! . Print out the title.
   IF ( PRESENT ( TITLE ) ) THEN
      WRITE ( OUTPUT, "(/A)" ) TITLE
   ELSE
      WRITE ( OUTPUT, "(/)" )
   END IF

   ! . Initialization.
   ICOL = 0

   ! . Top of the loop over rows.
   10 COLINC = MIN ( 12, ( N - ICOL ) )

   ! . Print out the column numbers.
   WRITE ( OUTPUT, "(12(1X,I6,4X))" ) ( ICOL+I, I = 1,COLINC )

   ! . Print out the diagonal part of the matrix.
   II = ((ICOL+1)*(ICOL+2))/2 - 1
   DO IROW = (ICOL+1),(ICOL+COLINC)
      WRITE ( OUTPUT, "(12(1X,F10.4))" ) MATRIX(II+1:II+(IROW-ICOL))
      II = II + IROW
   END DO

   ! . Print out the square part of the matrix.
   II = ((ICOL+COLINC)*(ICOL+COLINC+1))/2 + ICOL
   DO IROW = (ICOL+COLINC+1),N
      WRITE ( OUTPUT, "(12(1X,F10.4))" ) MATRIX(II+1:II+COLINC)
      II = II + IROW
   END DO

   ! . Increment the column number.
   ICOL = ICOL + COLINC

   ! . Check to see if more columns need to be printed.
   IF ( ICOL < N ) GO TO 10

   END SUBROUTINE PRINT_SYMMETRIC_MATRIX

!-------------------------------------------------------------------------------
! . Noformat printing.
!-------------------------------------------------------------------------------

   !------------------------------
   SUBROUTINE PRINT_NOFORMAT_START
   !------------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Set QNOFORMAT.
   QNOFORMAT = .TRUE.

   END SUBROUTINE PRINT_NOFORMAT_START

   !-----------------------------
   SUBROUTINE PRINT_NOFORMAT_STOP
   !-----------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Set QNOFORMAT.
   QNOFORMAT = .FALSE.

   END SUBROUTINE PRINT_NOFORMAT_STOP

!-------------------------------------------------------------------------------
! . Paragraph printing.
!-------------------------------------------------------------------------------

   !--------------------------------------------------------
   SUBROUTINE PRINT_PARAGRAPH ( TEXT, BOLD, CENTER, ITALIC )
   !--------------------------------------------------------

   ! . Arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: TEXT
   LOGICAL,               INTENT(IN), OPTIONAL :: BOLD, CENTER, ITALIC

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Print the paragraph.
   CALL PRINT_PARAGRAPH_START
   CALL PRINT_TEXT ( TEXT, BOLD, CENTER, ITALIC )
   CALL PRINT_PARAGRAPH_STOP

   END SUBROUTINE PRINT_PARAGRAPH

   !-------------------------------
   SUBROUTINE PRINT_PARAGRAPH_START
   !-------------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

      WRITE ( OUTPUT, "(1X)" )

   ! . Set QPARAGRAPH.
   QPARAGRAPH = .TRUE.

   END SUBROUTINE PRINT_PARAGRAPH_START

   !------------------------------
   SUBROUTINE PRINT_PARAGRAPH_STOP
   !------------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Set QPARAGRAPH.
   QPARAGRAPH = .FALSE.

   END SUBROUTINE PRINT_PARAGRAPH_STOP

!-------------------------------------------------------------------------------
! . Control and style print subroutines.
!-------------------------------------------------------------------------------

   !-------------------------------
   SUBROUTINE PRINT_START ( TITLE )
   !-------------------------------

   ! . Arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: TITLE

   ! . Check QPRINT.
   IF ( QPRINT ) CALL PRINT_ERROR ( "PRINT_START", "Printing has already started." )

   ! . Set the print flag.
   QPRINT = .TRUE.

   ! . Initialize PRINT_LINE.
   CALL PRINT_LINE_INITIALIZE

   END SUBROUTINE PRINT_START

   !--------------------
   SUBROUTINE PRINT_STOP
   !--------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check here for any unfinished paragraphs, summaries or tables.
   IF ( QNOFORMAT  ) CALL PRINT_NOFORMAT_STOP
   IF ( QPARAGRAPH ) CALL PRINT_PARAGRAPH_STOP
   IF ( QSUMMARY   ) CALL PRINT_SUMMARY_STOP
   IF ( QTABLE     ) CALL PRINT_TABLE_STOP

   END SUBROUTINE PRINT_STOP

   !------------------------------
   SUBROUTINE PRINT_STYLE ( HTML )
   !------------------------------

   ! . This subroutine can only be called once right at the beginning of the program.

   ! . Arguments.
   LOGICAL, INTENT(IN) :: HTML

   ! . Check QPRINT.
   IF ( QPRINT ) THEN
      CALL PRINT_ERROR ( "PRINT_STYLE", "The print style cannot be changed once printing has started." )
   END IF

   ! . Set the HTML flag.
   QHTML = HTML

   END SUBROUTINE PRINT_STYLE

!-------------------------------------------------------------------------------
! . Summary printing.
!-------------------------------------------------------------------------------

   !---------------------------------------------
   SUBROUTINE PRINT_SUMMARY_ELEMENT ( TAG, TEXT )
   !---------------------------------------------

   ! . Arguments.
   CHARACTER ( LEN = * ), INTENT(IN)           :: TAG
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: TEXT

   ! . Local scalars.
   INTEGER :: N

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QSUMMARY.
   IF ( .NOT. QSUMMARY ) CALL PRINT_ERROR ( "PRINT_SUMMARY_ELEMENT", "A SUMMARY block is not active." )

   ! . If TEXT is present copy it to PRINT_LINE
   IF ( PRESENT ( TEXT ) ) PRINT_LINE = TEXT

   ! . Get and check the print line length.
   PRINT_LINE_LENGTH = LEN_TRIM ( PRINT_LINE ) ; IF ( PRINT_LINE_LENGTH < 0 ) RETURN

      ! . Store the tag in the line (left justified).
      N = MIN ( LEN ( TAG ), SUMMARY_TAGWIDTH )
      SUMMARY_LINE(SUMMARY_LINE_LENGTH+1:SUMMARY_LINE_LENGTH+N) = TAG(1:N)

      ! . Increment SUMMARY_LINE_LENGTH.
      SUMMARY_LINE_LENGTH = SUMMARY_LINE_LENGTH + SUMMARY_TAGWIDTH

      ! . Add the equals sign.
      SUMMARY_LINE(SUMMARY_LINE_LENGTH-1:SUMMARY_LINE_LENGTH-1) = "="

      ! . Increment SUMMARY_LINE_LENGTH.
      SUMMARY_LINE_LENGTH = SUMMARY_LINE_LENGTH + SUMMARY_VARIABLEWIDTH

      ! . Store the variable in the line (right justified).
      N = MIN ( PRINT_LINE_LENGTH, SUMMARY_VARIABLEWIDTH )
      SUMMARY_LINE(SUMMARY_LINE_LENGTH-N+1:SUMMARY_LINE_LENGTH) = PRINT_LINE(1:N)

      ! . Write out the line if necessary.
      IF ( NSUMMARY == 1 ) WRITE ( OUTPUT, "(A)" ) SUMMARY_LINE(1:SUMMARY_LINE_LENGTH)

   ! . Reset NSUMMARY and some other variables.
   IF ( NSUMMARY == 0 ) THEN
      NSUMMARY            = 1
      SUMMARY_LINE_LENGTH = SUMMARY_TAGWIDTH + SUMMARY_VARIABLEWIDTH + 2
   ELSE
      NSUMMARY            = 0
      SUMMARY_LINE        = REPEAT ( " ", LINE_LENGTH )
      SUMMARY_LINE_LENGTH = 0
   END IF

   ! . Initialize PRINT_LINE.
   CALL PRINT_LINE_INITIALIZE

   END SUBROUTINE PRINT_SUMMARY_ELEMENT

   !-------------------------------
   SUBROUTINE PRINT_SUMMARY_ENDLINE
   !-------------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QSUMMARY.
   IF ( .NOT. QSUMMARY ) CALL PRINT_ERROR ( "PRINT_SUMMARY_ENDLINE", "A SUMMARY block is not active." )

   ! . Check NSUMMARY.
   IF ( NSUMMARY == 0 ) RETURN

      WRITE ( OUTPUT, "(A)" ) SUMMARY_LINE(1:SUMMARY_LINE_LENGTH)

   ! . Initialize some summary variables.
   NSUMMARY            = 0
   SUMMARY_LINE        = REPEAT ( " ", LINE_LENGTH )
   SUMMARY_LINE_LENGTH = 0

   END SUBROUTINE PRINT_SUMMARY_ENDLINE

   !----------------------------------
   SUBROUTINE PRINT_SUMMARY_INITIALIZE
   !----------------------------------

   ! . Initialize the data.
   NSUMMARY              = 0
   SUMMARY_LINE          = REPEAT ( " ", LINE_LENGTH )
   SUMMARY_LINE_LENGTH   = 0

   QSUMMARY              = .FALSE.

   SUMMARY_COLOR1        = SUMMARY_HEADER_COLOR
   SUMMARY_COLOR2        = SUMMARY_TAG_COLOR
   SUMMARY_COLOR3        = SUMMARY_TEXT_COLOR

   SUMMARY_PAGEWIDTH     = SUMMARY_PAGEWIDTH_DEFAULT
   SUMMARY_VARIABLEWIDTH = SUMMARY_VARIABLEWIDTH_DEFAULT

   SUMMARY_TAGWIDTH      = ( SUMMARY_PAGEWIDTH - 2 ) / 2 - SUMMARY_VARIABLEWIDTH

   END SUBROUTINE PRINT_SUMMARY_INITIALIZE

   !-------------------------------------------------------------------------------------------------
   SUBROUTINE PRINT_SUMMARY_OPTIONS ( HEADER_COLOR, TAG_COLOR, TEXT_COLOR, PAGEWIDTH, VARIABLEWIDTH )
   !-------------------------------------------------------------------------------------------------

   ! . HTML arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: HEADER_COLOR, TAG_COLOR, TEXT_COLOR

   ! . Text arguments.
   INTEGER, INTENT(IN), OPTIONAL :: PAGEWIDTH, VARIABLEWIDTH

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QSUMMARY.
   IF ( QSUMMARY ) CALL PRINT_ERROR ( "PRINT_SUMMARY_OPTIONS", "Cannot change options when a SUMMARY block is active." )

   ! . Initialize the SUMMARY data.
   CALL PRINT_SUMMARY_INITIALIZE

   ! . Set the summary options.
   IF ( PRESENT ( HEADER_COLOR  ) ) SUMMARY_COLOR1        = HEADER_COLOR
   IF ( PRESENT ( TAG_COLOR     ) ) SUMMARY_COLOR2        = TAG_COLOR
   IF ( PRESENT ( TEXT_COLOR    ) ) SUMMARY_COLOR3        = TEXT_COLOR
   IF ( PRESENT ( PAGEWIDTH     ) ) SUMMARY_PAGEWIDTH     = 2 * ( ( PAGEWIDTH + 1 ) / 2 )
   IF ( PRESENT ( VARIABLEWIDTH ) ) SUMMARY_VARIABLEWIDTH = VARIABLEWIDTH

   ! . Set TAGWIDTH.
   SUMMARY_TAGWIDTH = ( SUMMARY_PAGEWIDTH - 2 ) / 2 - SUMMARY_VARIABLEWIDTH

   END SUBROUTINE PRINT_SUMMARY_OPTIONS

   !---------------------------------------
   SUBROUTINE PRINT_SUMMARY_START ( TITLE )
   !---------------------------------------

   ! . Common arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: TITLE

   ! . Local scalars.
   INTEGER :: N, N1, N2

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QSUMMARY.
   IF ( QSUMMARY ) CALL PRINT_ERROR ( "PRINT_SUMMARY_START", "A SUMMARY block is already active." )

   ! . Initialize some summary variables.
   NSUMMARY            = 0
   SUMMARY_LINE        = REPEAT ( " ", LINE_LENGTH )
   SUMMARY_LINE_LENGTH = 0

      ! . Determine the spacing.
      N  = LEN ( TITLE ) + 2
      N1 = ( SUMMARY_PAGEWIDTH - N + 1 ) / 2
      N2 = SUMMARY_PAGEWIDTH - N - N1

      ! . Write out the header.
      WRITE ( OUTPUT, "(/A)" ) REPEAT ( "-", N1 ) // " " // TITLE // " " // REPEAT ( "-", N2 )

   ! . Set QSUMMARY.
   QSUMMARY = .TRUE.

   END SUBROUTINE PRINT_SUMMARY_START

   !----------------------------
   SUBROUTINE PRINT_SUMMARY_STOP
   !----------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QSUMMARY.
   IF ( .NOT. QSUMMARY ) CALL PRINT_ERROR ( "PRINT_SUMMARY_STOP", "A SUMMARY block is not active." )

   ! . Write out the remainder of the line if necessary.
   CALL PRINT_SUMMARY_ENDLINE

   ! . Write out the terminator.
      WRITE ( OUTPUT, "(A)" ) REPEAT ( "-", SUMMARY_PAGEWIDTH )

   ! . Initialize the SUMMARY variables.
   CALL PRINT_SUMMARY_INITIALIZE

   END SUBROUTINE PRINT_SUMMARY_STOP

!-------------------------------------------------------------------------------
! . Summary printing.
!-------------------------------------------------------------------------------

   !-----------------------------------------------------------------------------------
   SUBROUTINE PRINT_TABLE_ELEMENT ( ALIGN, COLOR, TEXT, COLSPAN, BOLD, HEADER, ITALIC )
   !-----------------------------------------------------------------------------------

   ! . Arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: ALIGN, COLOR, TEXT
   INTEGER,               INTENT(IN), OPTIONAL :: COLSPAN
   LOGICAL,               INTENT(IN), OPTIONAL :: BOLD, HEADER, ITALIC

   ! . Local parameters.
   CHARACTER ( LEN = 2 ) :: ENTRY_TAG = "TD", HEADER_TAG = "TH"

   ! . Local scalars.
   INTEGER :: N, NSTART, NSTOP, N1

   ! . Option scalars.
   CHARACTER ( LEN = 2 ) :: TAG
   CHARACTER ( LEN = 6 ) :: LOCAL_ALIGN
   CHARACTER ( LEN = 7 ) :: LOCAL_COLOR
   INTEGER               :: NCOLSPAN
   LOGICAL               :: QBOLD, QHEADER, QITALIC

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QTABLE.
   IF ( .NOT. QTABLE ) CALL PRINT_ERROR ( "PRINT_TABLE_ELEMENT", "A TABLE block is not active." )

   ! . If TEXT is present copy it to PRINT_LINE
   IF ( PRESENT ( TEXT ) ) PRINT_LINE = TEXT

   ! . Get and check the print line length.
   PRINT_LINE_LENGTH = LEN_TRIM ( PRINT_LINE ) ; IF ( PRINT_LINE_LENGTH < 0 ) RETURN

   ! . Get the number of spanning columns.
   NCOLSPAN = 1 ; IF ( PRESENT ( COLSPAN ) ) NCOLSPAN = COLSPAN

   ! . Find the starting and stopping columns.
   NSTART = NTABLE + 1
   NSTOP  = MIN ( NCOLUMNS, ( NTABLE + NCOLSPAN ) )

   ! . Get the HEADER option.
   QHEADER = .FALSE. ; IF ( PRESENT ( HEADER ) ) QHEADER = HEADER

   ! . Get the text alignment.
   IF ( PRESENT ( ALIGN ) ) THEN
      LOCAL_ALIGN = ALIGN
   ELSE
      IF ( QHEADER ) THEN
         LOCAL_ALIGN = "CENTER"
      ELSE
         LOCAL_ALIGN = "RIGHT"
      END IF
   END IF

   ! . HTML.
   IF ( QHTML ) THEN

      ! . Set the HTML options.
      QBOLD   = .FALSE. ; IF ( PRESENT ( BOLD	) ) QBOLD   = BOLD
      QITALIC = .FALSE. ; IF ( PRESENT ( ITALIC ) ) QITALIC = ITALIC

      ! . Get the color and appropriate tag.
      IF ( QHEADER ) THEN
         LOCAL_COLOR = TABLE_HEADER_COLOR
	 TAG         = HEADER_TAG
      ELSE
         LOCAL_COLOR = TABLE_ENTRY_COLOR
	 TAG         = ENTRY_TAG
      END IF
      IF ( PRESENT ( COLOR ) ) LOCAL_COLOR = COLOR

      ! . Write out the row beginning.
      IF ( NSTART == 1 ) WRITE ( OUTPUT, "(A)" ) "<TR>"

      ! . Write out the entry with the various options.
      WRITE ( OUTPUT, "(A,I2,A)" ) '<' // TAG // ' ALIGN = ' // LOCAL_ALIGN // ' BGCOLOR = "' // LOCAL_COLOR // &
                                   '" COLSPAN = ', NCOLSPAN , ' >'
      IF ( QBOLD   ) WRITE ( OUTPUT, "(A)" ) "<B>"
      IF ( QITALIC ) WRITE ( OUTPUT, "(A)" ) "<I>"
      WRITE ( OUTPUT, "(A)" ) PRINT_LINE(1:PRINT_LINE_LENGTH)
      IF ( QITALIC ) WRITE ( OUTPUT, "(A)" ) "</I>"
      IF ( QBOLD   ) WRITE ( OUTPUT, "(A)" ) "</B>"
      WRITE ( OUTPUT, "(A)" ) '</' // TAG // '>'

      ! . Write out the row ending.
      IF ( NSTOP == NCOLUMNS ) WRITE ( OUTPUT, "(A)" ) "</TR>"

   ! . Text.
   ELSE

      ! . Get the size of the column.
      N = SUM ( TABLE_VARIABLEWIDTHS(NSTART:NSTOP) )

      ! . Reduce PRINT_LINE_LENGTH if necessary.
      PRINT_LINE_LENGTH = MIN ( PRINT_LINE_LENGTH, N )

      ! . Store the entry (right justified is the default).
      SELECT CASE ( LOCAL_ALIGN )
      CASE ( "CENTER  " )
         N1 = ( N - PRINT_LINE_LENGTH + 1 ) / 2
         TABLE_LINE(TABLE_LINE_LENGTH+N1+1:TABLE_LINE_LENGTH+N1+PRINT_LINE_LENGTH) = PRINT_LINE(1:PRINT_LINE_LENGTH)
      CASE ( "LEFT    " )
         TABLE_LINE(TABLE_LINE_LENGTH+1:TABLE_LINE_LENGTH+PRINT_LINE_LENGTH)       = PRINT_LINE(1:PRINT_LINE_LENGTH)
      CASE DEFAULT
         TABLE_LINE(TABLE_LINE_LENGTH+N-PRINT_LINE_LENGTH+1:TABLE_LINE_LENGTH+N)   = PRINT_LINE(1:PRINT_LINE_LENGTH)
      END SELECT

      ! . Increment TABLE_LINE_LENGTH.
      TABLE_LINE_LENGTH = TABLE_LINE_LENGTH + N

      ! . Write out the line if necessary.
      IF ( NSTOP == NCOLUMNS ) THEN
         WRITE ( OUTPUT, "(A)" ) TABLE_LINE(1:TABLE_LINE_LENGTH)
	 IF ( QHEADER ) WRITE ( OUTPUT, "(A)" ) REPEAT ( "-", TABLE_PAGEWIDTH )
      END IF

   END IF

   ! . Reset NTABLE.
   NTABLE = NSTOP
   IF ( NTABLE == NCOLUMNS ) THEN
      NTABLE            = 0
      TABLE_LINE        = REPEAT ( " ", LINE_LENGTH )
      TABLE_LINE_LENGTH = 0
   END IF

   ! . Initialize PRINT_LINE.
   CALL PRINT_LINE_INITIALIZE

   END SUBROUTINE PRINT_TABLE_ELEMENT

   !-----------------------------
   SUBROUTINE PRINT_TABLE_ENDLINE
   !-----------------------------

   ! . Local scalars.
   INTEGER :: I

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QTABLE.
   IF ( .NOT. QTABLE ) CALL PRINT_ERROR ( "PRINT_TABLE_ENDLINE", "A TABLE block is not active." )

   ! . Check NTABLE.
   IF ( NTABLE == 0 ) RETURN

   ! . HTML.
   IF ( QHTML ) THEN

      DO I = (NTABLE+1),NCOLUMNS
         WRITE ( OUTPUT, "(A)" ) '<TD ALIGN = CENTER BGCOLOR = "' // TABLE_ENTRY_COLOR // '" >&#150;</TD>'
      END DO
      WRITE ( OUTPUT, "(A)" ) "</TR>"

   ! . Text.
   ELSE
      WRITE ( OUTPUT, "(A)" ) TABLE_LINE(1:TABLE_LINE_LENGTH)

   END IF

   ! . Initialize some table variables.
   NTABLE            = 0
   TABLE_LINE        = REPEAT ( " ", LINE_LENGTH )
   TABLE_LINE_LENGTH = 0

   END SUBROUTINE PRINT_TABLE_ENDLINE

   !--------------------------------
   SUBROUTINE PRINT_TABLE_INITIALIZE
   !--------------------------------

   ! . Initialize the data.
   NTABLE              = 0
   TABLE_LINE          = REPEAT ( " ", LINE_LENGTH )
   TABLE_LINE_LENGTH   = 0

   QTABLE              = .FALSE.

   NCOLUMNS            = TABLE_NCOLUMNS_DEFAULT

   TABLE_ENTRY_COLOR   = TABLE_ENTRY_COLOR_DEFAULT
   TABLE_HEADER_COLOR  = TABLE_HEADER_COLOR_DEFAULT

   TABLE_PAGEWIDTH     = TABLE_PAGEWIDTH_DEFAULT

   IF ( ALLOCATED ( TABLE_VARIABLEWIDTHS ) ) DEALLOCATE ( TABLE_VARIABLEWIDTHS )

   END SUBROUTINE PRINT_TABLE_INITIALIZE

   !-----------------------------------------------------------------------------------------------
   SUBROUTINE PRINT_TABLE_OPTIONS ( COLUMNS, ENTRY_COLOR, HEADER_COLOR, PAGEWIDTH, VARIABLEWIDTHS )
   !-----------------------------------------------------------------------------------------------

   ! . General arguments.
   INTEGER, INTENT(IN), OPTIONAL :: COLUMNS

   ! . HTML arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: ENTRY_COLOR, HEADER_COLOR

   ! . Text arguments.
   INTEGER,               INTENT(IN), OPTIONAL :: PAGEWIDTH
   INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: VARIABLEWIDTHS

   ! . Check QTABLE.
   IF ( QTABLE ) CALL PRINT_ERROR ( "PRINT_TABLE_OPTIONS", "Cannot change options when a TABLE block is active." )

   ! . Initialize the TABLE data.
   CALL PRINT_TABLE_INITIALIZE

   ! . Set the TABLE options.
   IF ( PRESENT ( COLUMNS      ) ) NCOLUMNS           = COLUMNS
   IF ( PRESENT ( ENTRY_COLOR  ) ) TABLE_ENTRY_COLOR  = ENTRY_COLOR
   IF ( PRESENT ( HEADER_COLOR ) ) TABLE_HEADER_COLOR = HEADER_COLOR
   IF ( PRESENT ( PAGEWIDTH    ) ) TABLE_PAGEWIDTH    = 2 * ( ( PAGEWIDTH + 1 ) / 2 )

   ! . Set TABLE_VARIABLEWIDTHS.
   IF ( PRESENT ( VARIABLEWIDTHS ) ) THEN
      IF ( ALLOCATED ( TABLE_VARIABLEWIDTHS ) ) DEALLOCATE ( TABLE_VARIABLEWIDTHS )
      ALLOCATE ( TABLE_VARIABLEWIDTHS(1:NCOLUMNS) ) ; TABLE_VARIABLEWIDTHS = VARIABLEWIDTHS
      IF ( SUM ( TABLE_VARIABLEWIDTHS ) /= TABLE_PAGEWIDTH ) DEALLOCATE ( TABLE_VARIABLEWIDTHS )
   END IF

   END SUBROUTINE PRINT_TABLE_OPTIONS

   !---------------------------
   SUBROUTINE PRINT_TABLE_START
   !---------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QTABLE.
   IF ( QTABLE ) CALL PRINT_ERROR ( "PRINT_TABLE_START", "A TABLE block is already active." )

   ! . Check NCOLUMNS.
   IF ( NCOLUMNS <= 0 ) CALL PRINT_ERROR ( "PRINT_TABLE_START", "This table has no columns." )

   ! . Initialize some table variables.
   NTABLE            = 0
   TABLE_LINE        = REPEAT ( " ", LINE_LENGTH )
   TABLE_LINE_LENGTH = 0

   ! . HTML.
   IF ( QHTML ) THEN
      WRITE ( OUTPUT, "(/A)" ) "<P>"
      WRITE ( OUTPUT, "(/A)" ) "<CENTER>"
      WRITE ( OUTPUT,  "(A)" ) "<TABLE BORDER = 1 CELLPADDING = 5>"

   ! . Text.
   ELSE

      ! . Make sure that TABLE_VARIABLEWIDTHS exists and is full.
      IF ( .NOT. ALLOCATED ( TABLE_VARIABLEWIDTHS ) ) THEN
         ALLOCATE ( TABLE_VARIABLEWIDTHS(1:NCOLUMNS) )
         TABLE_VARIABLEWIDTHS(1:NCOLUMNS-1) = ( TABLE_PAGEWIDTH + NCOLUMNS - 1 ) / NCOLUMNS
	 TABLE_VARIABLEWIDTHS(NCOLUMNS) = TABLE_PAGEWIDTH - SUM ( TABLE_VARIABLEWIDTHS(1:NCOLUMNS-1) )
      END IF

      ! . Write out the header.
      WRITE ( OUTPUT, "(/A)" ) REPEAT ( "-", TABLE_PAGEWIDTH )

   END IF

   ! . Set QTABLE.
   QTABLE = .TRUE.

   END SUBROUTINE PRINT_TABLE_START

   !--------------------------
   SUBROUTINE PRINT_TABLE_STOP
   !--------------------------

   ! . Check QPRINT.
   IF ( .NOT. QPRINT ) RETURN

   ! . Check QTABLE.
   IF ( .NOT. QTABLE ) CALL PRINT_ERROR ( "PRINT_TABLE_STOP", "A TABLE block is not active." )

   ! . Write out the remainder of the line if necessary.
   CALL PRINT_TABLE_ENDLINE

   ! . Write out the terminator.
   ! . HTML.
   IF ( QHTML ) THEN
      WRITE ( OUTPUT, "(/A)" ) "</TABLE>"
      WRITE ( OUTPUT, "(/A)" ) "</CENTER>"
      WRITE ( OUTPUT, "(/A)" ) "</P>"

   ! . Text.
   ELSE
      WRITE ( OUTPUT, "(A)" ) REPEAT ( "-", TABLE_PAGEWIDTH )

   END IF

   ! . Initialize the TABLE variables.
   CALL PRINT_TABLE_INITIALIZE

   END SUBROUTINE PRINT_TABLE_STOP

END MODULE PRINTING

