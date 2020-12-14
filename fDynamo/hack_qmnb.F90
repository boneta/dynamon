MODULE QMNB

#ifdef   HACK_QMNB

! .Module declarations.
USE DEFINITIONS,        ONLY : DP
USE ATOMS,              ONLY : ATMCRD, NATOMS, NATOMSQM, ATMQMI
USE ENERGY_NON_BONDING, ONLY : CUT_OFF, CUT_ON, NBLIST_TYPE, NBLISTQM_FIRST, QIMAGE, EPSILON
USE MM_TERMS,           ONLY : ATMEPS, ATMSIG, ATMCHG
USE SYMMETRY,           ONLY : BOXL
USE SEQUENCE,           ONLY : NRESID, RESIND

IMPLICIT NONE
PRIVATE
PUBLIC :: CALC_QMLJ, GRIM_INTERACTIONS, QMLJ_INTERACTIONS, SETUP_QMLJ_INTERACTIONS
SAVE

INTEGER, DIMENSION(:), ALLOCATABLE :: QM_IDX, QM_RES
REAL(KIND=DP) :: S_FACTOR = 1._DP, D_FACTOR = 1._dp
LOGICAL :: QMNB_CALC = .FALSE., GRIM_CALC = .FALSE.

!===============================================================================
CONTAINS
!===============================================================================

   !----------------------------------------------
   SUBROUTINE CALC_QMLJ ( QMLJ, GRADIENT )
   !----------------------------------------------
   REAL ( KIND = DP ), INTENT(INOUT) :: QMLJ
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT), OPTIONAL :: GRADIENT
   
   IF( PRESENT( GRADIENT ) ) THEN
       IF( QMNB_CALC  ) CALL QMLJ_INTERACTIONS( QMLJ, GRADIENT )
       IF( GRIM_CALC ) CALL GRIM_INTERACTIONS( QMLJ, GRADIENT )
   ELSE
       IF( QMNB_CALC  ) CALL QMLJ_INTERACTIONS( QMLJ )
       IF( GRIM_CALC ) CALL GRIM_INTERACTIONS( QMLJ )
   END IF

   END SUBROUTINE CALC_QMLJ

   !----------------------------------------------
   SUBROUTINE QMLJ_INTERACTIONS ( QMLJ, GRADIENT )
   !----------------------------------------------

   REAL ( KIND = DP ), INTENT(INOUT) :: QMLJ
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT), OPTIONAL :: GRADIENT

   INTEGER            :: II, I, JJ, J, K 
   REAL ( KIND = DP ) :: EI, EIJ, DF, R, RIJ2, S, S3, S6, SI, SIJ, S12
   REAL ( KIND = DP ) :: K6, K12, R2OFF, R2ON, SHIFT_LJ6, SHIFT_LJ12
   LOGICAL            :: QGRAD
   REAL ( KIND = DP ), DIMENSION(1:3) :: DR

   IF( .NOT. QMNB_CALC ) RETURN
   QGRAD = PRESENT( GRADIENT )
   QMLJ  = 0.0_DP
   R2OFF = CUT_OFF ** 2
   R2ON  = CUT_ON ** 2
   IF ( CUT_OFF > CUT_ON ) THEN
      K6    = ( CUT_OFF * R2OFF ) / ( CUT_OFF * R2OFF - CUT_ON * R2ON )
   ELSE
      K6    = 0.0_DP
      K12   = 0.0_DP
      R2OFF = 0.0_DP
      R2ON  = 0.0_DP
   END IF
   DO II = 1, NATOMSQM - 1
        I = QM_IDX(II)
       EI = ATMEPS(I)
       SI = ATMSIG(I)
       DO JJ = II + 1, NATOMSQM
           J = QM_IDX(JJ)
           IF( QM_RES(II) == QM_RES(JJ) ) CYCLE
           DR = ATMCRD(1:3,I) - ATMCRD(1:3,J)
           IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
           RIJ2 = DOT_PRODUCT ( DR, DR )
           IF ( RIJ2 > R2OFF ) CYCLE
           EIJ = EI * ATMEPS(J)
           SIJ = SI * ATMSIG(J)
           R   = SQRT ( RIJ2 )
           S   = 1.0_DP / R
           S3  = ( SIJ * S ) ** 3
           S6  = S3 * S3
           IF ( RIJ2 <= R2ON ) THEN
              S12 = S6 * S6
              SHIFT_LJ6  = ( ( SIJ / CUT_OFF ) * ( SIJ / CUT_ON ) ) ** 3
              SHIFT_LJ12 = SHIFT_LJ6 * SHIFT_LJ6
              QMLJ = QMLJ + EIJ * ( ( S12 - SHIFT_LJ12 ) * S_FACTOR - ( S6 - SHIFT_LJ6 ) * D_FACTOR )
              IF( .NOT. QGRAD ) CYCLE
              DF = 6.0_DP * EIJ * ( S6 * D_FACTOR - 2.0_DP * S12 * S_FACTOR  ) / RIJ2
           ELSE
              SHIFT_LJ6  = ( SIJ / CUT_OFF ) ** 3
              SHIFT_LJ12 = SHIFT_LJ6 * SHIFT_LJ6
              QMLJ = QMLJ + EIJ * ( K12 * ( S6 - SHIFT_LJ12 )**2 * S_FACTOR - K6 * ( S3 - SHIFT_LJ6 )**2  * D_FACTOR)
              IF( .NOT. QGRAD ) CYCLE
              DF = - 6.0_DP * EIJ * ( 2.0_DP * K12 * S6 * ( S6 - SHIFT_LJ12 ) * S_FACTOR - K6 * S3 * ( S3 - SHIFT_LJ6 ) * D_FACTOR ) / RIJ2
           END IF
           GRADIENT(1:3,I) = GRADIENT(1:3,I) + DF * DR
           GRADIENT(1:3,J) = GRADIENT(1:3,J) - DF * DR
       END DO
   END DO
   END SUBROUTINE QMLJ_INTERACTIONS

   !----------------------------------------------
   SUBROUTINE GRIM_INTERACTIONS ( QMLJ, GRADIENT )
   !----------------------------------------------

   REAL ( KIND = DP ), INTENT(INOUT) :: QMLJ
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT), OPTIONAL :: GRADIENT

   INTEGER            :: ii, I, jj, J, K 
   REAL ( KIND = DP ) :: EI, EIJ, DF, R, RIJ2, S, S3, S6, SI, SIJ, S12
   REAL ( KIND = DP ) :: K6, K12, R2OFF, R2ON, SHIFT_LJ6, SHIFT_LJ12
   REAL ( KIND = DP ) :: FR, DFR, SG, R0, EXD
   LOGICAL            :: QGRAD
   REAL ( KIND = DP ), DIMENSION(1:3) :: DR

   IF( .NOT. GRIM_CALC ) RETURN
   QGRAD = PRESENT( GRADIENT )
   QMLJ  = 0.0_DP
   R2OFF = CUT_OFF ** 2
   R2ON  = CUT_ON ** 2
   IF ( CUT_OFF > CUT_ON ) THEN
      K6    = ( CUT_OFF * R2OFF ) / ( CUT_OFF * R2OFF - CUT_ON * R2ON )
   ELSE
      K6    = 0.0_DP
      K12   = 0.0_DP
      R2OFF = 0.0_DP
      R2ON  = 0.0_DP
   END IF
   DO II = 1, NATOMSQM - 1
        I = QM_IDX(II)
       EI = ATMEPS(I)
       SI = ATMSIG(I)
       DO JJ = II + 1, NATOMSQM
           J = QM_IDX(JJ)
!           IF( QM_RES(II) == QM_RES(JJ) ) CYCLE
           DR = ATMCRD(1:3,I) - ATMCRD(1:3,J)
           IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
           RIJ2 = DOT_PRODUCT ( DR, DR )
           IF ( RIJ2 > R2OFF ) CYCLE
           EIJ = EI * ATMEPS(J)
           SIJ = SI * ATMSIG(J)
           IF (SIJ .EQ. 0.0_dp) CYCLE
           R   = SQRT ( RIJ2 )
           S   = 1.0_DP / R
           S3  = ( SIJ * S ) ** 3
           S6  = S3 * S3
   !GRIME FUNCTION
           SG  = S_FACTOR
           R0  = SIJ / ( 2._dp ** 6 )
           EXD = EXP( - D_FACTOR * ( R / R0 - 1._DP ) )
           FR  = SG / ( 1._dp + EXD )
           DFR = ( ( D_FACTOR / R0 ) * EXD ) / ( 1._DP + EXD ) ** 2
           IF ( RIJ2 <= R2ON ) THEN
              S12 = S6 * S6
              SHIFT_LJ6  = ( ( SIJ / CUT_OFF ) * ( SIJ / CUT_ON ) ) ** 3
              QMLJ = QMLJ + EIJ * (  - ( S6 - SHIFT_LJ6 ) ) * FR
              !write(100,*) QMLJ, QMLJ_TMP, FR
              IF( .NOT. QGRAD ) CYCLE
              DF = ( 6.0_DP * EIJ * ( S6 ) / RIJ2 ) * FR  +  EIJ * (  - ( S6 - SHIFT_LJ6 ) ) * DFR
           ELSE
              SHIFT_LJ6  = ( SIJ / CUT_OFF ) ** 3
              QMLJ = QMLJ + EIJ * ( - K6 * ( S3 - SHIFT_LJ6 )**2 ) * FR
              IF( .NOT. QGRAD ) CYCLE
              DF = ( 6.0_DP * EIJ * ( K6 * S3 * ( S3 - SHIFT_LJ6 ) ) / RIJ2 ) * FR + &
                   ( - K6 * ( S3 - SHIFT_LJ6 )**2 ) * DFR
           END IF
           GRADIENT(1:3,I) = GRADIENT(1:3,I) + DF * DR
           GRADIENT(1:3,J) = GRADIENT(1:3,J) - DF * DR
       END DO
   END DO

   END SUBROUTINE GRIM_INTERACTIONS     

   !----------------------------------------------
   SUBROUTINE SETUP_QMLJ_INTERACTIONS( GRIM, ALPHA, BETA )
   !----------------------------------------------
   IMPLICIT NONE
   REAL(KIND=DP),INTENT(IN),OPTIONAL :: ALPHA, BETA
   LOGICAL,INTENT(IN),OPTIONAL       :: GRIM 
   INTEGER       :: I, J, K
   LOGICAL       :: F

   IF( PRESENT  ( GRIM ) ) THEN
     GRIM_CALC = GRIM
     IF ( GRIM_CALC ) WRITE(*,*) &
        "---------------------------- GRIME FUCNTION ------------------------------------"
   END IF
     QMNB_CALC = .NOT. GRIM_CALC
     IF ( QMNB_CALC ) WRITE(*,*) &
        "---------------------------- QM NON BONDING ------------------------------------"
   

   IF( PRESENT  ( ALPHA  ) ) S_FACTOR = ALPHA
   IF( PRESENT  ( BETA   ) ) D_FACTOR = BETA 
   IF( ALLOCATED( QM_IDX ) ) DEALLOCATE( QM_IDX )
   IF( ALLOCATED( QM_RES ) ) DEALLOCATE( QM_RES )
   ALLOCATE( QM_IDX(1:NATOMSQM), QM_RES(1:NATOMSQM) )
   J = 0
   DO I = 1, NATOMS
      IF( ATMQMI(I) > 0 ) THEN
          J = J + 1
          QM_IDX(J) = I
          K = 1
          F = .FALSE.
          DO WHILE( .NOT. F .AND. K <= NRESID )
              F = I <= RESIND(K+1)
              IF( .NOT. F ) K = K + 1
          END DO
          QM_RES(J) = K
      END IF
   END DO
   END SUBROUTINE SETUP_QMLJ_INTERACTIONS

#endif

END MODULE QMNB
