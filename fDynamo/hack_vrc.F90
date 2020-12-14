MODULE VRC

#ifdef   HACK_VRC

! .Module declarations.
USE CONSTANTS,          ONLY : ELECT_CONST
USE DEFINITIONS,        ONLY : DP
USE FILES,              ONLY : NEXT_UNIT
USE ATOMS,              ONLY : ATMCRD, NATOMS
USE MM_TERMS,           ONLY : ATMCHG
USE PRINTING,           ONLY : PRINT_ERROR, PRINT_LINE, PRINT_TEXT
USE ENERGY_NON_BONDING, ONLY : EPSILON, CUT_OFF, CUT_ON, NBLIST_TYPE, NBLISTQM_FIRST, QIMAGE
USE SYMMETRY,           ONLY : BOXL

IMPLICIT NONE
PRIVATE
PUBLIC :: VRC_DEFINE, CALC_VRC, VRC_WRITING_START, &
          VRC_WRITING_STOP, VRC_INITIALIZE
SAVE

REAL( KIND=DP ) :: EQ, FC
LOGICAL :: VRC_WRITE = .FALSE., VRC_CALC = .FALSE.
INTEGER :: VRC_UNIT = -1, VRC_ITMS
INTEGER, DIMENSION(:), ALLOCATABLE :: VRC_LIST
REAL( KIND=DP ), DIMENSION(:), ALLOCATABLE :: VRC_COEF, VRC_VALU

!===============================================================================
CONTAINS
!===============================================================================

   !------------------------------------
   SUBROUTINE CALC_VRC( ENER, GRADIENT )
   !------------------------------------

   REAL( KIND = DP ), INTENT(INOUT) :: ENER
   REAL( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT), OPTIONAL :: GRADIENT

   INTEGER :: I, J, QM, MM
   TYPE( NBLIST_TYPE ), POINTER :: NBLIST


   REAL( KIND=DP ) :: EPSFAC, R2OFF, R2ON, GAMMA, A, B, CC, D, SHIFT_EL1, SHIFT_EL2
   REAL( KIND=DP ) :: QIJ, RIJ2, RR, S, TMP, DF, R3, R5, SS, SQ
   REAL( KIND=DP ), DIMENSION(1:3) :: DR

   ENER = .0_DP
   IF( .NOT. VRC_CALC ) RETURN

   EPSFAC    = ELECT_CONST / EPSILON
   R2OFF     = CUT_OFF**2
   R2ON      = CUT_ON**2
   GAMMA     = ( R2OFF - R2ON )**3
   A         = R2OFF * R2OFF * ( R2OFF - 3.0_DP * R2ON ) / GAMMA
   B         = 6.0_DP * R2OFF * R2ON / GAMMA
   CC        = - ( R2OFF + R2ON ) / GAMMA
   D         = 0.4_DP / GAMMA
   SHIFT_EL1 = 8.0_DP * ( R2OFF * R2ON * ( CUT_OFF - CUT_ON ) - 0.2_DP * &
            ( CUT_OFF * R2OFF * R2OFF - CUT_ON * R2ON * R2ON ) ) / GAMMA
   SHIFT_EL2 = - ( A / CUT_OFF ) + B * CUT_OFF + CC * CUT_OFF * R2OFF + D * CUT_OFF * R2OFF * R2OFF
   VRC_VALU  = .0_DP

   NBLIST => NBLISTQM_FIRST
   NBLOOP_ENE: DO WHILE( ASSOCIATED( NBLIST ) )
      DO I = 1, VRC_ITMS
         ! CURRENT ATOM IS SUPPOSED TO BE QUANTUM ONE
         IF( NBLIST%ATOM == VRC_LIST(I) ) THEN
            QM = NBLIST%ATOM
            DO J = 1, SIZE( NBLIST%INTERACTIONS )
               MM = NBLIST%INTERACTIONS(J)
!write(*,"(2i6)") qm, mm
! -----------------------------------------------------------------------
! - INTERACTION - INTERACTION - INTERACTION - INTERACTION - INTERACTION -
               QIJ = EPSFAC * ATMCHG(MM)
               DR  = ATMCRD(1:3,QM) - ATMCRD(1:3,MM)
               IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
               RIJ2 = DOT_PRODUCT ( DR, DR )
               IF ( RIJ2 > R2OFF ) CYCLE
               RR = SQRT ( RIJ2 )
               S  = 1.0_DP / RR
               IF ( RIJ2 <= R2ON ) THEN
                  TMP         = QIJ * S
                  VRC_VALU(I) = VRC_VALU(I) + TMP + QIJ * SHIFT_EL1
                  DF          = - TMP / RIJ2
               ELSE
                  R3          = RR * RIJ2
                  R5          = R3 * RIJ2
                  VRC_VALU(I) = VRC_VALU(I) + QIJ * ( A * S - B * RR - CC * R3 - D * R5 + SHIFT_EL2 )
               END IF
! - INTERACTION - INTERACTION - INTERACTION - INTERACTION - INTERACTION -
! -----------------------------------------------------------------------
            END DO
         ELSE
         ! SEARCH IN THE INTERACTION VECTOR FOR THE PRESENCE OF THE QM ATOM
            MM = NBLIST%ATOM
            DO J = 1, SIZE( NBLIST%INTERACTIONS )
               QM = ABS( NBLIST%INTERACTIONS(J) )
               IF( QM == VRC_LIST(I) ) THEN
!write(*,"(2i6)") nblist%interactions(j), mm
! -----------------------------------------------------------------------
! - INTERACTION - INTERACTION - INTERACTION - INTERACTION - INTERACTION -
                  QIJ = EPSFAC * ATMCHG(MM)
                  DR  = ATMCRD(1:3,QM) - ATMCRD(1:3,MM)
                  IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
                  RIJ2 = DOT_PRODUCT ( DR, DR )
                  IF ( RIJ2 > R2OFF ) CYCLE
                  RR = SQRT ( RIJ2 )
                  S  = 1.0_DP / RR
                  IF ( RIJ2 <= R2ON ) THEN
                     TMP         = QIJ * S
                     VRC_VALU(I) = VRC_VALU(I) + TMP + QIJ * SHIFT_EL1
                  ELSE
                     R3          = RR * RIJ2
                     R5          = R3 * RIJ2
                     VRC_VALU(I) = VRC_VALU(I) + QIJ * ( A * S - B * RR - CC * R3 - D * R5 + SHIFT_EL2 )
                  END IF
! - INTERACTION - INTERACTION - INTERACTION - INTERACTION - INTERACTION -
! -----------------------------------------------------------------------
               END IF
            END DO
         END IF
      END DO
      NBLIST => NBLIST%NEXT_LIST
   END DO NBLOOP_ENE

   SQ   = DOT_PRODUCT( VRC_COEF, VRC_VALU )
   SS   = SQ - EQ
   ENER = 0.5_DP * FC * SS * SS

   IF( VRC_WRITE ) THEN
       WRITE( VRC_UNIT, "(F20.10)" ) SQ
       CALL FLUSH( VRC_UNIT )
   END IF

   IF( PRESENT( GRADIENT ) ) THEN

      NBLIST => NBLISTQM_FIRST
      NBLOOP_GRAD: DO WHILE( ASSOCIATED( NBLIST ) )
         DO I = 1, VRC_ITMS
            ! CURRENT ATOM IS SUPPOSED TO BE QUANTUM ONE
            IF( NBLIST%ATOM == VRC_LIST(I) ) THEN
               QM = NBLIST%ATOM
               DO J = 1, SIZE( NBLIST%INTERACTIONS )
                  MM = NBLIST%INTERACTIONS(J)
! -----------------------------------------------------------------------
! - INTERACTION - INTERACTION - INTERACTION - INTERACTION - INTERACTION -
                  QIJ = EPSFAC * ATMCHG(MM)
                  DR  = ATMCRD(1:3,QM) - ATMCRD(1:3,MM)
                  IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
                  RIJ2 = DOT_PRODUCT ( DR, DR )
                  IF ( RIJ2 > R2OFF ) CYCLE
                  RR = SQRT ( RIJ2 )
                  S  = 1.0_DP / RR
                  IF ( RIJ2 <= R2ON ) THEN
                     TMP = QIJ * S
                     DF  = - TMP / RIJ2
                  ELSE
                     R3  = RR * RIJ2
                     R5  = R3 * RIJ2
                     DF  = - QIJ * ( A / R3 + B * S + 3.0_DP * CC * RR + 5.0_DP * D * R3 ) 
                  END IF
                  GRADIENT(1:3,QM) = GRADIENT(1:3,QM) + DF * DR * FC * SS * VRC_COEF(I)
                  GRADIENT(1:3,MM) = GRADIENT(1:3,MM) - DF * DR * FC * SS * VRC_COEF(I)
! - INTERACTION - INTERACTION - INTERACTION - INTERACTION - INTERACTION -
! -----------------------------------------------------------------------
               END DO
            ELSE
            ! SEARCH IN THE INTERACTION VECTOR FOR THE PRESENCE OF THE QM ATOM
               MM = NBLIST%ATOM
               DO J = 1, SIZE( NBLIST%INTERACTIONS )
                  QM = ABS( NBLIST%INTERACTIONS(J) )
                  IF( QM == VRC_LIST(I) ) THEN
! -----------------------------------------------------------------------
! - INTERACTION - INTERACTION - INTERACTION - INTERACTION - INTERACTION -
                     QIJ = EPSFAC * ATMCHG(MM)
                     DR  = ATMCRD(1:3,QM) - ATMCRD(1:3,MM)
                     IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
                     RIJ2 = DOT_PRODUCT ( DR, DR )
                     IF ( RIJ2 > R2OFF ) CYCLE
                     RR = SQRT ( RIJ2 )
                     S  = 1.0_DP / RR
                     IF ( RIJ2 <= R2ON ) THEN
                        TMP = QIJ * S
                        DF  = - TMP / RIJ2
                     ELSE
                        R3  = RR * RIJ2
                        R5  = R3 * RIJ2
                        DF  = - QIJ * ( A / R3 + B * S + 3.0_DP * CC * RR + 5.0_DP * D * R3 ) 
                     END IF
                     GRADIENT(1:3,QM) = GRADIENT(1:3,QM) + DF * DR * FC * SS * VRC_COEF(I)
                     GRADIENT(1:3,MM) = GRADIENT(1:3,MM) - DF * DR * FC * SS * VRC_COEF(I)
! - INTERACTION - INTERACTION - INTERACTION - INTERACTION - INTERACTION -
! -----------------------------------------------------------------------
                  END IF
               END DO
            END IF
         END DO
         NBLIST => NBLIST%NEXT_LIST
      END DO NBLOOP_GRAD

   END IF

   END SUBROUTINE CALC_VRC


   !----------------------------------
   SUBROUTINE VRC_WRITING_START
   !----------------------------------

   IF( VRC_UNIT > -1 ) THEN
       VRC_WRITE = .TRUE.
       WRITE( VRC_UNIT, "(2F20.10)" ) FC, EQ
   END IF

   END SUBROUTINE VRC_WRITING_START


   !----------------------------------
   SUBROUTINE VRC_WRITING_STOP
   !----------------------------------

   IF( VRC_WRITE ) THEN
       CALL FLUSH( VRC_UNIT )
       CLOSE( VRC_UNIT )
       VRC_UNIT  = -1
       VRC_WRITE = .FALSE.
   END IF

   END SUBROUTINE VRC_WRITING_STOP

   !--------------------------------------------------------------------------
   SUBROUTINE VRC_DEFINE ( ITMS, LIST, COEF, REFERENCE, FORCE_CONSTANT, FILE )
   !--------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: ITMS
   INTEGER, DIMENSION(1:ITMS) :: LIST
   REAL( KIND=DP ), DIMENSION(1:ITMS) :: COEF
   REAL( KIND=DP ), INTENT(IN) :: REFERENCE, FORCE_CONSTANT
   CHARACTER( LEN=* ), INTENT(IN), OPTIONAL :: FILE

   EQ = REFERENCE
   FC = FORCE_CONSTANT

   IF( PRESENT( FILE ) ) THEN
       VRC_UNIT = NEXT_UNIT()
       OPEN( UNIT = VRC_UNIT, FILE = TRIM( FILE ), ACTION = "WRITE", FORM = "FORMATTED" )
   END IF

   VRC_CALC = .TRUE.

   VRC_ITMS = ITMS
   ALLOCATE( VRC_LIST(1:ITMS), VRC_COEF(1:ITMS), VRC_VALU(1:ITMS) )
   VRC_LIST = LIST
   VRC_COEF = COEF

   WRITE( PRINT_LINE, "(A)" ) "- Potential Energy RC Setup ------------------------------------------------"
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,I6)" ) "Num. Atoms  :", ITMS
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,F20.10)" ) "Reference   :", EQ
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,F20.10)" ) "Force Const.:", FC
   CALL PRINT_TEXT
   IF( VRC_UNIT > -1 ) THEN
       WRITE( PRINT_LINE, "(A,A)" )  "File        :", TRIM( FILE )
       CALL PRINT_TEXT
   END IF
   WRITE( PRINT_LINE, "(A)" ) "----------------------------------------------------------------------------"
   CALL PRINT_TEXT

   END SUBROUTINE VRC_DEFINE

   !------------------------
   SUBROUTINE VRC_INITIALIZE
   !------------------------

   EQ       = .0_DP
   FC       = .0_DP
   VRC_CALC = .FALSE.
   VRC_ITMS = 0 

   IF( ALLOCATED( VRC_LIST ) ) DEALLOCATE( VRC_LIST )
   IF( ALLOCATED( VRC_COEF ) ) DEALLOCATE( VRC_COEF )
   IF( ALLOCATED( VRC_VALU ) ) DEALLOCATE( VRC_VALU )

   END SUBROUTINE VRC_INITIALIZE

#endif

END MODULE VRC
