MODULE HTSC

#ifdef	HACK_HTSC

! .Module declarations.
USE CONSTANTS,          ONLY : TO_DEGREES
USE DEFINITIONS,        ONLY : DP
USE FILES,              ONLY : NEXT_UNIT

USE ATOMS,              ONLY : ATMCRD, NATOMS, atmnum
use elements,           only : symbol

USE CONJUGATE_GRADIENT, ONLY : CONJUGATE_GRADIENT_MINIMIZE

USE QUANTUM_ENERGY,     ONLY : ENERGY_QUANTUM
USE ENERGY_NON_BONDING, ONLY : CUT_OFF, CUT_ON, NBLIST_TYPE, NBLISTQM_FIRST, QIMAGE
USE MM_TERMS,           ONLY : ATMEPS, ATMSIG
USE SYMMETRY,           ONLY : BOXL
USE PRINTING,           ONLY : PRINT_ERROR, PRINT_LINE, PRINT_TEXT
USE MOPAC_DATA,         ONLY : DENMAT, NBASTR

IMPLICIT NONE
PRIVATE
PUBLIC :: HTRANSFER_DEFINE, ENERGY_HTRANSFER, GRADIENT_HTRANSFER, HTRANSFER_WRITING_START, &
          HTRANSFER_WRITING_STOP, HTRANSFER_CLEAN, HTRANSFER_REFERENCE, HTSC_METHOD, HTSC_NSTEPS
SAVE

INTEGER :: HTSC_METHOD, HTSC_NSTEPS
INTEGER :: D_ATOM, A_ATOM, H_ATOM
INTEGER :: D_IDX, A_IDX, H_IDX
REAL( KIND=DP ) :: EQ, FC
INTEGER, DIMENSION(:), ALLOCATABLE :: HTSC_IDX
REAL( KIND=DP ), DIMENSION(:,:), ALLOCATABLE :: HTSC_CRD, HTSC_GRD1, HTSC_GRD2
LOGICAL :: HTSC_WRITE = .FALSE., HTSC_CALC = .FALSE.
INTEGER :: HTSC_UNIT = -1, HTSC_NAT
INTEGER :: WHOSYOURDADDY
REAL( KIND=DP ), ALLOCATABLE, DIMENSION(:) :: xDENMAT, d_DENMAT, a_DENMAT

!===============================================================================
CONTAINS
!===============================================================================

! -
! - current doesn't perform 1-4 interactions
! - the reason is that the LA should be far enough from the reactive atoms...
! - this could be changed "easily" through subroutine parameters: ATMEPS14, ...
! -
   !----------------------------------------------
   SUBROUTINE QMLJ_INTERACTIONS ( QMLJ, GRADIENT )
   !----------------------------------------------

   REAL ( KIND = DP ), INTENT(INOUT) :: QMLJ
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT), OPTIONAL :: GRADIENT

   INTEGER            :: I, IINT, J
   REAL ( KIND = DP ) :: EI, EIJ, DF, R, RIJ2, S, S3, S6, S12, SI, SIJ
   REAL ( KIND = DP ) :: K6, K12, R2OFF, R2ON, SHIFT_LJ6, SHIFT_LJ12
   LOGICAL            :: QGRAD

   REAL ( KIND = DP ), DIMENSION(1:3) :: DR

   TYPE(NBLIST_TYPE), POINTER :: NBLIST

   QGRAD = PRESENT( GRADIENT )

   QMLJ  = 0.0_DP
   R2OFF = CUT_OFF ** 2
   R2ON  = CUT_ON ** 2

   IF ( CUT_OFF > CUT_ON ) THEN

      K6    = ( CUT_OFF * R2OFF ) / ( CUT_OFF * R2OFF - CUT_ON * R2ON )
      K12   = ( R2OFF ** 3 ) / ( R2OFF ** 3 - R2ON ** 3 )

   ELSE

      K6        = 0.0_DP ; K12       = 0.0_DP
      R2OFF     = 0.0_DP ; R2ON      = 0.0_DP

   END IF

   NBLIST => NBLISTQM_FIRST
   DO WHILE( ASSOCIATED( NBLIST ) )

      I  = ABS( NBLIST%ATOM )
      EI = ATMEPS(I)
      SI = ATMSIG(I)

      DO IINT = 1, SIZE(NBLIST%INTERACTIONS)

         J = ABS( NBLIST%INTERACTIONS(IINT) )

         IF( I /= H_ATOM .AND. J /= H_ATOM ) CYCLE

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

            S12        = S6 * S6
            SHIFT_LJ6  = ( ( SIJ / CUT_OFF ) * ( SIJ / CUT_ON ) ) ** 3
            SHIFT_LJ12 = SHIFT_LJ6 * SHIFT_LJ6

            QMLJ = QMLJ + EIJ * ( ( S12 - SHIFT_LJ12 ) - ( S6 - SHIFT_LJ6 ) )

            IF( .NOT. QGRAD ) CYCLE
            DF = 6.0_DP * EIJ * ( S6 - 2.0_DP * S12 ) / RIJ2

         ELSE

            SHIFT_LJ6  = ( SIJ / CUT_OFF ) ** 3
            SHIFT_LJ12 = SHIFT_LJ6 * SHIFT_LJ6

            QMLJ = QMLJ + EIJ * ( K12 * ( S6 - SHIFT_LJ12 ) ** 2 - K6 * ( S3 - SHIFT_LJ6 ) ** 2 )

            IF( .NOT. QGRAD ) CYCLE
			DF = - 6.0_DP * EIJ * ( 2.0_DP * K12 * S6 * ( S6 - SHIFT_LJ12 ) - K6 * S3 * ( S3 - SHIFT_LJ6 ) ) / RIJ2

         END IF

         GRADIENT(1:3,I) = GRADIENT(1:3,I) + DF * DR
         GRADIENT(1:3,J) = GRADIENT(1:3,J) - DF * DR

      END DO

      NBLIST => NBLIST%NEXT_LIST

   END DO

   END SUBROUTINE QMLJ_INTERACTIONS

!--------------------------------------------------------------------------------------
SUBROUTINE LINEARIZE( X, F, G, THETA, KUMB )
   REAL( KIND=DP ), INTENT(IN) :: THETA, KUMB
   REAL( KIND=DP ), INTENT(INOUT) :: F
   REAL( KIND=DP ), DIMENSION(:), INTENT(IN) :: X
   REAL( KIND=DP ), DIMENSION(:), INTENT(INOUT) :: G

   INTEGER :: I
   REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP - 1.0E-6_DP
   REAL ( KIND = DP ) :: DF, DIFF, DOTFAC, RIJ, RKJ, VALUE
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ, DRKJ, DTI, DTJ, DTK
   
   DRIJ(1:3) = (/ X( 3 * ( D_IDX - 1 ) + 1 ) - X(  3 * ( H_IDX - 1 ) + 1 ), &
                  X( 3 * ( D_IDX - 1 ) + 2 ) - X(  3 * ( H_IDX - 1 ) + 2 ), &
                  X( 3 * ( D_IDX - 1 ) + 3 ) - X(  3 * ( H_IDX - 1 ) + 3 ) /)
   RIJ       = SQRT ( DOT_PRODUCT ( DRIJ, DRIJ ) )
   DRIJ      = DRIJ / RIJ

   DRKJ(1:3) = (/ X( 3 * ( A_IDX - 1 ) + 1 ) - X(  3 * ( H_IDX - 1 ) + 1 ), &
                  X( 3 * ( A_IDX - 1 ) + 2 ) - X(  3 * ( H_IDX - 1 ) + 2 ), &
                  X( 3 * ( A_IDX - 1 ) + 3 ) - X(  3 * ( H_IDX - 1 ) + 3 ) /)
   RKJ       = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )
   DRKJ      = DRKJ / RKJ
   
   DOTFAC = DOT_PRODUCT ( DRIJ, DRKJ )
   DOTFAC = SIGN ( MIN ( ABS ( DOTFAC ), DOT_LIMIT ), DOTFAC )
   VALUE  = ACOS ( DOTFAC ) * TO_DEGREES
   
   DIFF = VALUE - THETA
   DF   = KUMB * DIFF
   F    = F + 0.5_DP * DF * DIFF
   DF   = - TO_DEGREES * DF / SQRT ( 1.0_DP - DOTFAC * DOTFAC )
   DTI  = ( DRKJ - DOTFAC * DRIJ ) / RIJ
   DTK  = ( DRIJ - DOTFAC * DRKJ ) / RKJ
   DTJ  = - ( DTI + DTK )
   
   DO I = 1, 3
      G( 3 * ( D_IDX - 1 ) + I ) = G( 3 * ( D_IDX - 1 ) + I ) + DF * DTI(I)
      G( 3 * ( H_IDX - 1 ) + I ) = G( 3 * ( H_IDX - 1 ) + I ) + DF * DTJ(I)
      G( 3 * ( A_IDX - 1 ) + I ) = G( 3 * ( A_IDX - 1 ) + I ) + DF * DTK(I)
   END DO
END SUBROUTINE LINEARIZE

SUBROUTINE NORMALIZE( X, F, G, DIST, KUMB )
   REAL( KIND=DP ), INTENT(IN) :: DIST, KUMB
   REAL( KIND=DP ), INTENT(INOUT) :: F
   REAL( KIND=DP ), DIMENSION(:), INTENT(IN) :: X
   REAL( KIND=DP ), DIMENSION(:), INTENT(INOUT) :: G

   INTEGER :: I
   REAL ( KIND = DP ) :: DF, DIFF, VALUE
   REAL ( KIND = DP ), DIMENSION(1:3) :: DR
   
   DR(1:3) = (/ X( 3 * ( WHOSYOURDADDY - 1 ) + 1 ) - X(  3 * ( H_IDX - 1 ) + 1 ), &
                X( 3 * ( WHOSYOURDADDY - 1 ) + 2 ) - X(  3 * ( H_IDX - 1 ) + 2 ), &
                X( 3 * ( WHOSYOURDADDY - 1 ) + 3 ) - X(  3 * ( H_IDX - 1 ) + 3 ) /)
   VALUE   = SQRT ( DOT_PRODUCT ( DR, DR ) )

   DIFF = VALUE - DIST
   DF   = KUMB * DIFF
   F    = F + 0.5_DP * DF * DIFF
   
   DO I = 1, 3
      G( 3 * ( WHOSYOURDADDY - 1 ) + I ) = G( 3 * ( WHOSYOURDADDY - 1 ) + I ) + DF * DR(I)
      G( 3 * ( H_IDX - 1 ) + I )         = G( 3 * ( H_IDX - 1 ) + I )         - DF * DR(I)
   END DO
END SUBROUTINE NORMALIZE

SUBROUTINE EGCALC( X, F, G )
   REAL( KIND=DP ), INTENT(OUT) :: F
   REAL( KIND=DP ), DIMENSION(:), INTENT(IN) :: X
   REAL( KIND=DP ), DIMENSION(:), INTENT(OUT) :: G
   
   INTEGER :: I, J
   REAL( KIND=DP ) :: E_LJ, E_QM, T1, T2
   REAL( KIND=DP ), DIMENSION(1:3,1:NATOMS) :: LG
   
   DO I = 1, HTSC_NAT
       J = ( I - 1 ) * 3
       ATMCRD(1:3,HTSC_IDX(I)) = X(J+1:J+3)
   END DO
write(900,"(i6/)") htsc_nat
do i = 1, htsc_nat
write(900,"(a2,3f20.10)") symbol(atmnum(htsc_idx(i))), atmcrd(1:3,htsc_idx(i))
end do
   LG = .0_DP
   CALL QMLJ_INTERACTIONS( E_LJ, LG )
   CALL ENERGY_QUANTUM( E_QM, T1, T2, GRADIENT = LG, PRINT = .FALSE. )
   DO I = 1, HTSC_NAT
       J = ( I - 1 ) * 3
       G(J+1:J+3)  = LG(1:3,HTSC_IDX(I))
   END DO
   F = E_QM + E_LJ
   CALL LINEARIZE( X, F, G, 179.0_DP,     5.0_DP )
   CALL NORMALIZE( X, F, G,   1.0_DP, 10000.0_DP )
END SUBROUTINE EGCALC


SUBROUTINE LOCAL_CG( NSTEP, TOLG, MAXSTP, ENERGY, GRADIENT )
   INTEGER, INTENT(IN) :: NSTEP
   REAL( KIND=DP ), INTENT(IN) :: TOLG, MAXSTP
   REAL( KIND=DP ), INTENT(INOUT) :: ENERGY
   REAL( KIND=DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT) :: GRADIENT
   
   INTEGER :: I, J, IFLG
   REAL( KIND=DP ) :: RR, FF
   REAL( KIND=DP ), DIMENSION(:), ALLOCATABLE  :: X
   
   ALLOCATE(X(1:3*HTSC_NAT))
   DO I = 1, HTSC_NAT
       J = ( I - 1 ) * 3
       X(J+1:J+3) = ATMCRD(1:3,HTSC_IDX(I))
   END DO
   CALL CONJUGATE_GRADIENT_MINIMIZE( EGCALC, X, IFLG, 0, NSTEP, MAXSTP, TOLG )
   DO I = 1, HTSC_NAT
       J = ( I - 1 ) * 3
       ATMCRD(1:3,HTSC_IDX(I)) = X(J+1:J+3)
   END DO
   ENERGY   = .0_DP
   GRADIENT = .0_DP
   CALL QMLJ_INTERACTIONS( ENERGY, GRADIENT )
   CALL ENERGY_QUANTUM( ENERGY, RR, FF, GRADIENT = GRADIENT, PRINT = .FALSE. )
   DEALLOCATE(X)
END SUBROUTINE LOCAL_CG


SUBROUTINE LOCAL_LBFGSB( NSTEP, TOLG, ENERGY, GRADIENT )
   INTEGER, INTENT(IN) :: NSTEP
   REAL( KIND=DP ), INTENT(IN) :: TOLG
   REAL( KIND=DP ), INTENT(INOUT) :: ENERGY
   REAL( KIND=DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT) :: GRADIENT
   
   CHARACTER( LEN=60 ) :: TASK, CSAVE
   LOGICAL, DIMENSION(1:4) :: LSAVE
   INTEGER :: I, J, K, N, MX
   INTEGER, DIMENSION(1:44) :: ISAVE
   INTEGER, DIMENSION(:), ALLOCATABLE :: NBD, IWA
   REAL( KIND=DP ) :: F, FF, GG, RR
   REAL( KIND=DP ), DIMENSION(1:29) :: DSAVE
   REAL( KIND=DP ), DIMENSION(:), ALLOCATABLE  :: X, G, L, U, WA
   
   N = 3 * HTSC_NAT
   MX = 5
   RR = DSQRT( REAL( N, KIND=DP ) )
   ALLOCATE( X(1:N), G(1:N), NBD(1:N), L(1:N), U(1:N), IWA(1:3*N), WA(1:2*MX*N+4*N+12*MX*MX+12*MX) )
   DO I = 1, HTSC_NAT
       J = ( I - 1 ) * 3
       X(J+1:J+3) = ATMCRD(1:3,HTSC_IDX(I))
   END DO
   CALL EGCALC( X, F, G )
   GG = DSQRT( DOT_PRODUCT( G, G ) ) / RR
   FF = 1.E+3_DP
   NBD = 0
   L = .0_DP
   U = .0_DP
   K = 0
   TASK = "START"
   DO WHILE( K < NSTEP .AND. GG > TOLG .AND. &
   		TASK(1:4) /= "CONV" .AND. TASK(1:4) /= "STOP" .AND. TASK(1:4) /= "ERRO" .AND. TASK(1:4) /= "ABNO" )
   	CALL SETULB( N, MX, X, L, U, NBD, F, G, FF, TOLG, WA, IWA, TASK, -1, CSAVE, LSAVE, ISAVE, DSAVE )
   	IF( TASK(1:2) == "FG" ) THEN
   		CALL EGCALC( X, F, G )
   		GG = DSQRT( DOT_PRODUCT( G, G ) ) / RR
   	END IF
   	K = K + 1
   END DO
   DO I = 1, HTSC_NAT
       J = ( I - 1 ) * 3
       ATMCRD(1:3,HTSC_IDX(I)) = X(J+1:J+3)
   END DO
   ENERGY   = .0_DP
   GRADIENT = .0_DP
   CALL QMLJ_INTERACTIONS( ENERGY, GRADIENT )
   CALL ENERGY_QUANTUM( ENERGY, RR, FF, GRADIENT = GRADIENT, PRINT = .FALSE. )
   DEALLOCATE( X, G, NBD, L, U, IWA, WA )
END SUBROUTINE LOCAL_LBFGSB
!--------------------------------------------------------------------------------------

   !------------------------------------
   SUBROUTINE HTRANSFER_REFERENCE
   !------------------------------------
   INTEGER :: I
   REAL( KIND = DP ) :: E1, E2, DE
   REAL( KIND = DP ), DIMENSION(1:3) :: VEC_AD

   IF( .NOT. HTSC_CALC ) RETURN

   xDENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   VEC_AD = ATMCRD(1:3,D_ATOM) - ATMCRD(1:3,A_ATOM)
   VEC_AD = VEC_AD / DSQRT( SUM( VEC_AD ** 2 ) )

   DO I = 1, HTSC_NAT
       HTSC_CRD(1:3,I) = ATMCRD(1:3,HTSC_IDX(I))
   END DO

   ATMCRD(1:3,H_ATOM) = ATMCRD(1:3,D_ATOM) - VEC_AD(1:3) * 1.0_DP
   DENMAT(1:NBASTR,1) = d_DENMAT(1:NBASTR)
   WHOSYOURDADDY = D_IDX
   IF( HTSC_METHOD == 1 ) THEN
      CALL LOCAL_LBFGSB( HTSC_NSTEPS, 1._DP, E1, HTSC_GRD1 )
   ELSE
      CALL LOCAL_CG( HTSC_NSTEPS, 1._DP, 0.01_DP, E1, HTSC_GRD1 )
   END IF
   d_DENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   DO I = 1, HTSC_NAT
       ATMCRD(1:3,HTSC_IDX(I)) = HTSC_CRD(1:3,I)
   END DO

   ATMCRD(1:3,H_ATOM) = ATMCRD(1:3,A_ATOM) + VEC_AD(1:3) * 1.0_DP
   DENMAT(1:NBASTR,1) = a_DENMAT(1:NBASTR)
   WHOSYOURDADDY = A_IDX
   IF( HTSC_METHOD == 1 ) THEN
      CALL LOCAL_LBFGSB( HTSC_NSTEPS, 1._DP, E2, HTSC_GRD2 )
   ELSE
      CALL LOCAL_CG( HTSC_NSTEPS, 1._DP, 0.01_DP, E2, HTSC_GRD2 )
   END IF
   a_DENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   DENMAT(1:NBASTR,1) = xDENMAT(1:NBASTR)

   DO I = 1, HTSC_NAT
       ATMCRD(1:3,HTSC_IDX(I)) = HTSC_CRD(1:3,I)
   END DO

   DE = E1 - E2
   WRITE(*,"(A,F20.10)" ) "Energy_Transfer: ", DE

   END SUBROUTINE HTRANSFER_REFERENCE

   !------------------------------------
   SUBROUTINE ENERGY_HTRANSFER ( ECORR )
   !------------------------------------
   REAL ( KIND = DP ), INTENT(INOUT) :: ECORR

   INTEGER :: I
   REAL( KIND = DP ) :: E1, E2, DE
   REAL( KIND = DP ), DIMENSION(1:3) :: VEC_AD

   ECORR = 0.0_DP
   IF( .NOT. HTSC_CALC ) RETURN

   xDENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   VEC_AD = ATMCRD(1:3,D_ATOM) - ATMCRD(1:3,A_ATOM)
   VEC_AD = VEC_AD / DSQRT( SUM( VEC_AD ** 2 ) )

   DO I = 1, HTSC_NAT
       HTSC_CRD(1:3,I) = ATMCRD(1:3,HTSC_IDX(I))
   END DO

   ATMCRD(1:3,H_ATOM) = ATMCRD(1:3,D_ATOM) - VEC_AD(1:3) * 1.0_DP
   DENMAT(1:NBASTR,1) = d_DENMAT(1:NBASTR)
   WHOSYOURDADDY = D_IDX
   IF( HTSC_METHOD == 1 ) THEN
      CALL LOCAL_LBFGSB( HTSC_NSTEPS, 1._DP, E1, HTSC_GRD1 )
   ELSE
      CALL LOCAL_CG( HTSC_NSTEPS, 1._DP, 0.01_DP, E1, HTSC_GRD1 )
   END IF

   DO I = 1, HTSC_NAT
       ATMCRD(1:3,HTSC_IDX(I)) = HTSC_CRD(1:3,I)
   END DO
   d_DENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   ATMCRD(1:3,H_ATOM) = ATMCRD(1:3,A_ATOM) + VEC_AD(1:3) * 1.0_DP
   DENMAT(1:NBASTR,1) = a_DENMAT(1:NBASTR)
   WHOSYOURDADDY = A_IDX
   IF( HTSC_METHOD == 1 ) THEN
      CALL LOCAL_LBFGSB( HTSC_NSTEPS, 1._DP, E2, HTSC_GRD2 )
   ELSE
      CALL LOCAL_CG( HTSC_NSTEPS, 1._DP, 0.01_DP, E2, HTSC_GRD2 )
   END IF
   a_DENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   DENMAT(1:NBASTR,1) = xDENMAT(1:NBASTR)

   DO I = 1, HTSC_NAT
       ATMCRD(1:3,HTSC_IDX(I)) = HTSC_CRD(1:3,I)
   END DO

   DE       = E1 - E2
   ECORR    = 0.5_DP * FC * ( DE - EQ ) * ( DE - EQ )

   IF( HTSC_WRITE ) THEN
       WRITE( HTSC_UNIT, "(F20.10)" ) DE
       CALL FLUSH( HTSC_UNIT )
   END IF

   END SUBROUTINE ENERGY_HTRANSFER


   !----------------------------------------------
   SUBROUTINE GRADIENT_HTRANSFER ( ECORR, GRADIENT )
   !----------------------------------------------
   REAL ( KIND = DP ), INTENT(INOUT) :: ECORR
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT) :: GRADIENT

   INTEGER :: I
   REAL( KIND = DP ) :: E1, E2, DE
   REAL( KIND = DP ), DIMENSION(1:3) :: VEC_AD

   ECORR = 0.0_DP
   IF( .NOT. HTSC_CALC ) RETURN

   xDENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   VEC_AD = ATMCRD(1:3,D_ATOM) - ATMCRD(1:3,A_ATOM)
   VEC_AD = VEC_AD / DSQRT( SUM( VEC_AD ** 2 ) )

   DO I = 1, HTSC_NAT
       HTSC_CRD(1:3,I) = ATMCRD(1:3,HTSC_IDX(I))
   END DO

   ATMCRD(1:3,H_ATOM) = ATMCRD(1:3,D_ATOM) - VEC_AD(1:3) * 1.0_DP
   DENMAT(1:NBASTR,1) = d_DENMAT(1:NBASTR)
   WHOSYOURDADDY = D_IDX
   IF( HTSC_METHOD == 1 ) THEN
      CALL LOCAL_LBFGSB( HTSC_NSTEPS, 1._DP, E1, HTSC_GRD1 )
   ELSE
      CALL LOCAL_CG( HTSC_NSTEPS, 1._DP, 0.01_DP, E1, HTSC_GRD1 )
   END IF
   d_DENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

write(1024,"(i6/)") htsc_nat
do i = 1, htsc_nat
write(1024,"(a2,3f20.10)") symbol(atmnum(htsc_idx(i))), atmcrd(1:3,htsc_idx(i))
end do
call flush( 1024 )

   DO I = 1, HTSC_NAT
       ATMCRD(1:3,HTSC_IDX(I)) = HTSC_CRD(1:3,I)
   END DO

   ATMCRD(1:3,H_ATOM) = ATMCRD(1:3,A_ATOM) + VEC_AD(1:3) * 1.0_DP
   DENMAT(1:NBASTR,1) = a_DENMAT(1:NBASTR)
   WHOSYOURDADDY = A_IDX
   IF( HTSC_METHOD == 1 ) THEN
      CALL LOCAL_LBFGSB( HTSC_NSTEPS, 1._DP, E2, HTSC_GRD2 )
   ELSE
      CALL LOCAL_CG( HTSC_NSTEPS, 1._DP, 0.01_DP, E2, HTSC_GRD2 )
   END IF
   a_DENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   DENMAT(1:NBASTR,1) = xDENMAT(1:NBASTR)

write(1025,"(i6/)") htsc_nat
do i = 1, htsc_nat
write(1025,"(a2,3f20.10)") symbol(atmnum(htsc_idx(i))), atmcrd(1:3,htsc_idx(i))
end do
call flush( 1025 )

   DO I = 1, HTSC_NAT
       ATMCRD(1:3,HTSC_IDX(I)) = HTSC_CRD(1:3,I)
   END DO

   ! - cleanup involved atoms gradient...
   DO I = 1, HTSC_NAT
       HTSC_GRD1(1:3,HTSC_IDX(I)) = .0_DP
       HTSC_GRD2(1:3,HTSC_IDX(I)) = .0_DP
   END DO

   DE       = E1 - E2
   ECORR    = 0.5_DP * FC * ( DE - EQ ) * ( DE - EQ )
   GRADIENT = GRADIENT + FC * ( DE - EQ ) * ( HTSC_GRD1 - HTSC_GRD2 )

   IF( HTSC_WRITE ) THEN
       WRITE( HTSC_UNIT, "(F20.10)" ) DE
       CALL FLUSH( HTSC_UNIT )
   END IF

   END SUBROUTINE GRADIENT_HTRANSFER


   !----------------------------------
   SUBROUTINE HTRANSFER_WRITING_START
   !----------------------------------

   IF( HTSC_UNIT > -1 ) THEN
       HTSC_WRITE = .TRUE.
       WRITE( HTSC_UNIT, "(2F20.10)" ) FC, EQ
   END IF

   END SUBROUTINE HTRANSFER_WRITING_START


   !----------------------------------
   SUBROUTINE HTRANSFER_WRITING_STOP
   !----------------------------------

   IF( HTSC_WRITE ) THEN
       CALL FLUSH( HTSC_UNIT )
       CLOSE( HTSC_UNIT )
       HTSC_UNIT  = -1
       HTSC_WRITE = .FALSE.
   END IF

   END SUBROUTINE HTRANSFER_WRITING_STOP


   !-----------------------------------------------------------------------------------------------------
   SUBROUTINE HTRANSFER_DEFINE ( DONNOR, ACCEPTOR, HYDROGEN, REFERENCE, FORCE_CONSTANT, SELECTION, FILE )
   !-----------------------------------------------------------------------------------------------------
   INTEGER, INTENT(IN) :: DONNOR, ACCEPTOR, HYDROGEN
   REAL( KIND=DP ), INTENT(IN) :: REFERENCE, FORCE_CONSTANT
   CHARACTER( LEN=* ), INTENT(IN), OPTIONAL :: FILE
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION

   INTEGER :: I, J

   D_ATOM = DONNOR
   A_ATOM = ACCEPTOR
   H_ATOM = HYDROGEN
   EQ     = REFERENCE
   FC     = FORCE_CONSTANT

   ALLOCATE( xDENMAT(1:NBASTR) )
   ALLOCATE( d_DENMAT(1:NBASTR), a_DENMAT(1:NBASTR) )

   HTSC_NAT = COUNT( SELECTION )
   ALLOCATE( HTSC_IDX(1:HTSC_NAT), HTSC_CRD(1:3,1:HTSC_NAT) )
   J = 1
   DO I = 1, NATOMS
       IF( SELECTION(I) ) THEN
           HTSC_IDX(J) = I
           IF( I == D_ATOM ) D_IDX = J
           IF( I == A_ATOM ) A_IDX = J
           IF( I == H_ATOM ) H_IDX = J
           J = J + 1
       END IF
   END DO

   IF( PRESENT( FILE ) ) THEN
       HTSC_UNIT = NEXT_UNIT()
       OPEN( UNIT = HTSC_UNIT, FILE = TRIM( FILE ), ACTION = "WRITE", FORM = "FORMATTED" )
   END IF

   HTSC_CALC  = .TRUE.

   WRITE( PRINT_LINE, "(A)" ) "- Solvent Coordinate Setup -------------------------------------------------"
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,I20)" ) "Donnor: ", D_ATOM
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,I20)" ) "Acceptor: ", A_ATOM
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,I20)" ) "Hydrogen: ", H_ATOM
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,F20.10)" ) "Reference Energy: ", EQ
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,F20.10)" ) "Force Constant: ", FC
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,I20)" ) "Selected Atoms: ", HTSC_NAT
   CALL PRINT_TEXT
   IF( HTSC_UNIT > -1 ) THEN
       WRITE( PRINT_LINE, "(A,A)" ) "File: ", TRIM( FILE )
       CALL PRINT_TEXT
   END IF
   WRITE( PRINT_LINE, "(A)" ) "----------------------------------------------------------------------------"
   CALL PRINT_TEXT

   ALLOCATE( HTSC_GRD1(1:3,1:NATOMS), HTSC_GRD2(1:3,1:NATOMS) )

   ! - default to Conjugate Gradient
   HTSC_METHOD = 0
   HTSC_NSTEPS = 1000

open( unit=1024, file="dh.opt", action = "write", form = "formatted" )
open( unit=1025, file="ah.opt", action = "write", form = "formatted" )

   d_DENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)
   a_DENMAT(1:NBASTR) = DENMAT(1:NBASTR,1)

   END SUBROUTINE HTRANSFER_DEFINE


   SUBROUTINE HTRANSFER_CLEAN

close( 1024 )
close( 1025 )

   IF( ALLOCATED( HTSC_IDX ) )  DEALLOCATE( HTSC_IDX )
   IF( ALLOCATED( HTSC_CRD ) )  DEALLOCATE( HTSC_CRD )
   IF( ALLOCATED( HTSC_GRD1 ) ) DEALLOCATE( HTSC_GRD1 )
   IF( ALLOCATED( HTSC_GRD2 ) ) DEALLOCATE( HTSC_GRD2 )

   IF( ALLOCATED(  xDENMAT ) ) DEALLOCATE( xDENMAT )
   IF( ALLOCATED( d_DENMAT ) ) DEALLOCATE( d_DENMAT )
   IF( ALLOCATED( a_DENMAT ) ) DEALLOCATE( a_DENMAT )

   END SUBROUTINE HTRANSFER_CLEAN

#endif

END MODULE HTSC
