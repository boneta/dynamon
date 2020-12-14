MODULE URC

#ifdef	HACK_URC

! .Module declarations.
USE CONSTANTS,    ONLY : AMUA2PS2_TO_KJMOL
USE DEFINITIONS,  ONLY : DP
USE FILES,        ONLY : NEXT_UNIT
USE ATOMS,        ONLY : NATOMS, ATMMAS, ATMFIX
USE VELOCITY,     ONLY : ATMVEL
USE PRINTING,     ONLY : PRINT_ERROR, PRINT_LINE, PRINT_TEXT

IMPLICIT NONE
PRIVATE
PUBLIC :: URC_DEFINE, GRADIENT_URC, URC_WRITING_START, &
          URC_WRITING_STOP, URC_INITIALIZE
SAVE

REAL( KIND=DP ) :: EQ, FC

LOGICAL :: URC_WRITE = .FALSE., URC_CALC = .FALSE.
INTEGER :: URC_UNIT = -1

!===============================================================================
CONTAINS
!===============================================================================

   !----------------------------------------
   SUBROUTINE GRADIENT_URC( EPOT, GRADIENT )
   !----------------------------------------

   REAL( KIND = DP ), INTENT(INOUT) :: EPOT
   REAL( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT) :: GRADIENT

   INTEGER :: I
   REAL( KIND = DP ) :: CORR, ETOT

   IF( .NOT. URC_CALC ) RETURN

   CORR = .0_DP
   DO I = 1, NATOMS
       IF( .NOT. ATMFIX(I) ) CORR = CORR + ATMMAS(I) * DOT_PRODUCT( ATMVEL(1:3,I), ATMVEL(1:3,I) )
   END DO
   ETOT = EPOT + 0.5_DP * AMUA2PS2_TO_KJMOL * CORR

   CORR = 0.5_DP * FC * ( ETOT - EQ ) * ( ETOT - EQ )
   GRADIENT = GRADIENT + FC * ( ETOT - EQ ) * GRADIENT
   EPOT = EPOT + CORR

   IF( URC_WRITE ) THEN
       WRITE( URC_UNIT, "(F20.10)" ) ETOT + CORR
       CALL FLUSH( URC_UNIT )
   END IF


   END SUBROUTINE GRADIENT_URC


   !----------------------------------
   SUBROUTINE URC_WRITING_START
   !----------------------------------

   IF( URC_UNIT > -1 ) THEN
       URC_WRITE = .TRUE.
       WRITE( URC_UNIT, "(2F20.10)" ) FC, EQ
   END IF

   END SUBROUTINE URC_WRITING_START


   !----------------------------------
   SUBROUTINE URC_WRITING_STOP
   !----------------------------------

   IF( URC_WRITE ) THEN
       CALL FLUSH( URC_UNIT )
       CLOSE( URC_UNIT )
       URC_UNIT  = -1
       URC_WRITE = .FALSE.
   END IF

   END SUBROUTINE URC_WRITING_STOP

   !--------------------------------------------------------
   SUBROUTINE URC_DEFINE ( REFERENCE, FORCE_CONSTANT, FILE )
   !--------------------------------------------------------

   ! . Scalar arguments.
   REAL( KIND=DP ), INTENT(IN) :: REFERENCE, FORCE_CONSTANT
   CHARACTER( LEN=* ), INTENT(IN), OPTIONAL :: FILE

   EQ = REFERENCE
   FC = FORCE_CONSTANT

   IF( PRESENT( FILE ) ) THEN
       URC_UNIT = NEXT_UNIT()
       OPEN( UNIT = URC_UNIT, FILE = TRIM( FILE ), ACTION = "WRITE", FORM = "FORMATTED" )
   END IF

   URC_CALC = .TRUE.

   WRITE( PRINT_LINE, "(A)" ) "- Internal Energy RC Setup -------------------------------------------------"
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,F20.10)" ) "Reference   :", EQ
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,F20.10)" ) "Force Const.:", FC
   CALL PRINT_TEXT
   IF( URC_UNIT > -1 ) THEN
       WRITE( PRINT_LINE, "(A,A)" )  "File        :", TRIM( FILE )
       CALL PRINT_TEXT
   END IF
   WRITE( PRINT_LINE, "(A)" ) "----------------------------------------------------------------------------"
   CALL PRINT_TEXT

   END SUBROUTINE URC_DEFINE

   !------------------------
   SUBROUTINE URC_INITIALIZE
   !------------------------

   EQ        = .0_DP
   FC        = .0_DP
   URC_CALC  = .FALSE.

   END SUBROUTINE URC_INITIALIZE

#endif

END MODULE URC
