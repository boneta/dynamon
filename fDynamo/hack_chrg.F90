! ==============================================================================
! DYNAMO / WHAT-EVER Interface
! ==============================================================================
MODULE ABINITIO


USE ATOMS,              ONLY : NATOMS, ATMCRD, ATMQMI, ATMNUM, NATOMSMM, NATOMSQM, NFREE, ATMFIX, ATMMAS, atmnam
USE ENERGY_NON_BONDING, ONLY : NBLIST_TYPE, NBLISTQM_FIRST, QIMAGE, CUT_ON, CUT_OFF
USE DEFINITIONS,        ONLY : DP, SKIP_ABINITIO
USE MM_TERMS,           ONLY : ANGLES, ATMEPS, ATMEPS14, ATMCHG, ATMCHG14, &
                               ATME14I, ATME14J, BONDS, DIHEDRALS, IMPROPERS, &
                               NANGLES, NBONDS, NDIHEDRALS, NIMPROPERS
USE CONSTANTS,          ONLY : AU_TO_KJ, BOHRS_TO_ANGSTROMS
USE SYMMETRY,           ONLY : BOXL
USE MOPAC_DATA,         ONLY : QMATOM, MOPAC_DATA_ATOM_POSITION, QUSEBAQM
USE PRINTING,           ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
                               PRINT_SUMMARY_START, PRINT_SUMMARY_STOP

IMPLICIT NONE

PRIVATE
PUBLIC :: ABINITIO_INIT, ABINITIO_SETUP, ABINITIO_ENERGY, ABINITIO_GRADIENT, ABINITIO_HESSIAN, ABINITIO_EXIT

CONTAINS


SUBROUTINE ABINITIO_NON_BONDED( SELECTION )
    IMPLICIT NONE
    LOGICAL, DIMENSION(1:NATOMS), INTENT(INOUT) :: SELECTION

    INTEGER                       :: I
    TYPE( NBLIST_TYPE ), POINTER :: NBLIST

    SELECTION = .FALSE.

!	selection based on cut-off (no switching function or so...)
!    NBLIST => NBLISTQM_FIRST
!    NBLOOP: DO WHILE( ASSOCIATED( NBLIST ) )
!        DO I = 1, SIZE( NBLIST%INTERACTIONS )
!            ! . Negative values are designed to QM atoms
!            IF( NBLIST%INTERACTIONS(I) > 0 ) THEN
!                SELECTION(NBLIST%INTERACTIONS(I)) = .TRUE.
!            ELSE
!                SELECTION(NBLIST%ATOM) = .TRUE.
!            END IF
!        END DO
!        NBLIST => NBLIST%NEXT_LIST
!    END DO NBLOOP

	! - trunctaion scheme (take a big one in order to avoid problems... ideally all the mobile atoms)
	CALL MY_SELE( SELECTION )

	! - skip MM-Link counterpart from the interaction...
	WHERE( ATMQMI > 0 ) SELECTION = .FALSE.
END SUBROUTINE

! ==============================================================================
! ==============================================================================

SUBROUTINE ABINITIO_INIT
    IMPLICIT NONE
    SKIP_ABINITIO = .FALSE.
END SUBROUTINE

! ==============================================================================
! ==============================================================================

SUBROUTINE ABINITIO_REMOVE_MM_TERMS( SELECTION )
    IMPLICIT NONE
    LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION

    INTEGER                                  :: I, IINT, J, N, NANG, NBND, NDIH, NE14, NIMP, TMP
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: EX14I, EX14J

    N = 0
    DO I = 1,NANGLES
       TMP = 0
       IF( SELECTION( ANGLES(I)%I ) ) TMP = TMP + 1
       IF( SELECTION( ANGLES(I)%J ) ) TMP = TMP + 1
       IF( SELECTION( ANGLES(I)%K ) ) TMP = TMP + 1
       IF( TMP == 3 .OR. TMP == 2 ) CYCLE
!       IF( TMP == 3 ) CYCLE
       N = N + 1
       ANGLES(N) = ANGLES(I)
    END DO
    NANG    = NANGLES - N
    NANGLES = N

    N = 0
    DO I = 1,NBONDS
       TMP = 0
       IF( SELECTION( BONDS(I)%I ) ) TMP = TMP + 1
       IF( SELECTION( BONDS(I)%J ) ) TMP = TMP + 1
       IF( TMP == 2 ) CYCLE
       N = N + 1
       BONDS(N) = BONDS(I)
    END DO
    NBND   = NBONDS - N
    NBONDS = N

    N = 0
    DO I = 1,NDIHEDRALS
       TMP = 0
       IF( SELECTION( DIHEDRALS(I)%I ) ) TMP = TMP + 1
       IF( SELECTION( DIHEDRALS(I)%J ) ) TMP = TMP + 1
       IF( SELECTION( DIHEDRALS(I)%K ) ) TMP = TMP + 1
       IF( SELECTION( DIHEDRALS(I)%L ) ) TMP = TMP + 1
       IF( TMP == 4 .OR. TMP == 3 ) CYCLE
!       IF( TMP == 4 ) CYCLE
       N = N + 1
       DIHEDRALS(N) = DIHEDRALS(I)
    END DO
    NDIH       = NDIHEDRALS - N
    NDIHEDRALS = N

    N = 0
    DO I = 1,NIMPROPERS
       TMP = 0
       IF( SELECTION( IMPROPERS(I)%I ) ) TMP = TMP + 1
       IF( SELECTION( IMPROPERS(I)%J ) ) TMP = TMP + 1
       IF( SELECTION( IMPROPERS(I)%K ) ) TMP = TMP + 1
       IF( SELECTION( IMPROPERS(I)%L ) ) TMP = TMP + 1
       IF( TMP == 4 .OR. TMP == 3 ) CYCLE
!       IF( TMP == 4 ) CYCLE
       N = N + 1
       IMPROPERS(N) = IMPROPERS(I)
    END DO
    NIMP       = NIMPROPERS - N
    NIMPROPERS = N

    N = 0
    ALLOCATE( EX14I(1:NATOMS+1), EX14J(1:SIZE(ATME14J)) )
    DO I = 1,NATOMS
       EX14I(I) = N
       DO IINT = (ATME14I(I)+1),ATME14I(I+1)
          J = ATME14J(IINT)
          IF( SELECTION(I) .AND. SELECTION(J) ) CYCLE
          N        = N + 1
          EX14J(N) = J
       END DO
    END DO
    EX14I(NATOMS+1) = N
    NE14 = SIZE( ATME14J ) - N
    DEALLOCATE( ATME14I, ATME14J )
    ALLOCATE( ATME14I(1:NATOMS+1), ATME14J(1:N) )
    ATME14I(1:NATOMS+1) = EX14I(1:NATOMS+1)
    ATME14J(1:N)        = EX14J(1:N)
    DEALLOCATE( EX14I, EX14J )

    CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800" )
    CALL PRINT_SUMMARY_START ( "MM Terms Removed" )
    WRITE ( PRINT_LINE, "(I14)" ) NBND ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Bonds"     )
    WRITE ( PRINT_LINE, "(I14)" ) NANG ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Angles"    )
    WRITE ( PRINT_LINE, "(I14)" ) NDIH ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Dihedrals" )
    WRITE ( PRINT_LINE, "(I14)" ) NIMP ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Impropers" )
    WRITE ( PRINT_LINE, "(I14)" ) NE14 ; CALL PRINT_SUMMARY_ELEMENT ( "Number of 1-4 Excl." )
    WRITE ( PRINT_LINE, "(F14.4)" ) SUM( ATMCHG(1:NATOMS) ); CALL PRINT_SUMMARY_ELEMENT ( "Total MM Charge " )
    CALL PRINT_SUMMARY_STOP

END SUBROUTINE

! ==============================================================================
! ==============================================================================

SUBROUTINE ABINITIO_SETUP( SELECTION )
    IMPLICIT NONE
    LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION

    INTEGER                                  :: I, K1, K2
    INTEGER, DIMENSION(:), ALLOCATABLE       :: IDX

    ALLOCATE( IDX(1:NATOMS) )
    IDX = 0
    DO I = 1, NATOMS
        IF( SELECTION(I) ) IDX(I) = I
    END DO
    K1 = 0
    DO I = 1, NBONDS
        IF( SELECTION(BONDS(I)%I) .AND. .NOT. SELECTION(BONDS(I)%J) ) THEN
            K1 = K1 + 1
            IDX(BONDS(I)%J) = - BONDS(I)%I
        ELSE IF( .NOT. SELECTION(BONDS(I)%I) .AND. SELECTION(BONDS(I)%J) ) THEN
            K1 = K1 + 1
            IDX(BONDS(I)%I) = - BONDS(I)%J
        END IF
    END DO

    NATOMSQM = COUNT( SELECTION ) + K1
    NATOMSMM = NATOMS - NATOMSQM

    IF( ALLOCATED( QMATOM ) ) DEALLOCATE( QMATOM )
    ALLOCATE( QMATOM(1:NATOMSQM) )

    ! . the idea is placing the link-atoms at the end of the QMATOM list
    ATMQMI = 0
    K1     = 1
    K2     = NATOMSQM
    DO I = 1, NATOMS
        IF( IDX(I) > 0 ) THEN
            ATMQMI(I)            = K1
            ATMCHG(I)            = 0._DP
            ATMCHG14(I )         = 0._DP
            QMATOM(K1)%ATOM      = I
            QMATOM(K1)%BFIRST    = 0
            QMATOM(K1)%BLAST     = 0
            QMATOM(K1)%COPY      = 0
            QMATOM(K1)%NUMBER    = ATMNUM(I)
            QMATOM(K1)%LENGTH    = 0._DP
            QMATOM(K1)%PARTNER   = 0
            QMATOM(K1)%QBOUNDARY = .FALSE.
            K1                   = K1 + 1
        END IF
        IF( IDX(I) < 0 ) THEN
            ATMQMI(I)            = K2
            QMATOM(K2)%ATOM      = I
            QMATOM(K2)%BFIRST    = 0
            QMATOM(K2)%BLAST     = 0
            QMATOM(K2)%COPY      = 0
            QMATOM(K2)%NUMBER    = 1
            QMATOM(K2)%LENGTH    = 1._DP
            QMATOM(K2)%PARTNER   = - IDX(I)
            QMATOM(K2)%QBOUNDARY = .TRUE.
            K2                   = K2 - 1
        END IF
    END DO
    DEALLOCATE( IDX )
    CALL ABINITIO_REMOVE_MM_TERMS( SELECTION )
     QUSEBAQM     = .FALSE.
    SKIP_ABINITIO = .FALSE.
END SUBROUTINE

! ==============================================================================
! ==============================================================================

SUBROUTINE SHOW_CHARGES
	IMPLICIT NONE
	INTEGER :: I, J
	CHARACTER( LEN=100 ) :: MSG_LBL, MSG_CHG

	RETURN   ! SABOTAGE OF CHARGES PRINTING

	WRITE( *, "(/10('-'),A)" ) " AB INITIO CHARGES"
	J = 0
	MSG_LBL = ""
	MSG_CHG = ""
	DO I = 1, NATOMSQM
		WRITE( MSG_LBL(J*10+1:J*10+10),   "(A10)" ) ATMNAM(QMATOM(I)%ATOM)
		WRITE( MSG_CHG(J*10+1:J*10+10), "(F10.4)" ) ATMCHG(QMATOM(I)%ATOM)
		J = J + 1
		IF( J == 8 ) THEN
			WRITE( *, "(A80)" ) MSG_LBL(1:80)
			WRITE( *, "(A80)" ) MSG_CHG(1:80)
			J = 0
			MSG_LBL = ""
			MSG_CHG = ""
		END IF
	END DO
	IF( J /= 0 ) THEN
		WRITE( *, "(A80)" ) MSG_LBL(1:80)
		WRITE( *, "(A80)" ) MSG_CHG(1:80)
	END IF
	WRITE( *, "(80('-'))" )
END SUBROUTINE SHOW_CHARGES

! ==============================================================================
! ==============================================================================

SUBROUTINE ABINITIO_ENERGY( ENERGY )
    IMPLICIT NONE

    REAL ( KIND = DP ), INTENT(INOUT)          :: ENERGY

    REAL( KIND=DP )                            :: EXT_ENER
    INTEGER                                    :: I, J, K, NT, CODE
    LOGICAL, ALLOCATABLE, DIMENSION(:)         :: SELE
    REAL( KIND=DP ), DIMENSION(1:3)            :: TMP_CRD, QMGC
    REAL( KIND=DP ), DIMENSION(:), ALLOCATABLE :: EXT_CORD, EXT_ATMN, EXT_QFIT, EXT_GRAD, EXT_HESS

    ALLOCATE( SELE(1:NATOMS) )
    CALL ABINITIO_NON_BONDED( SELE )

    NT = NATOMSQM + COUNT( SELE )

    ALLOCATE( EXT_CORD(1:3*NT), EXT_ATMN(1:NT), EXT_QFIT(1:NATOMSQM) )
	ALLOCATE( EXT_GRAD(1), EXT_HESS(1) )

	QMGC = .0_DP
    DO I = 1, NATOMSQM
        EXT_ATMN(I) = QMATOM(I)%NUMBER * 1._DP
        TMP_CRD = MOPAC_DATA_ATOM_POSITION( I, ATMCRD )
        J = ( I - 1 ) * 3
        EXT_CORD(J + 1) = TMP_CRD(1)
        EXT_CORD(J + 2) = TMP_CRD(2)
        EXT_CORD(J + 3) = TMP_CRD(3)
		QMGC(1:3) = QMGC(1:3) + TMP_CRD(1:3)
    END DO
	QMGC(1:3) = QMGC(1:3) / REAL( NATOMSQM, KIND = DP )

    J = NATOMSQM
    DO I = 1, NATOMS
        IF( SELE(I) ) THEN
            J = J + 1
!            IF( QIMAGE ) THEN
!				TMP_CRD(1:3) = ATMCRD(1:3,I) - QMGC(1:3)
!                TMP_CRD(1:3) = TMP_CRD(1:3) - BOXL(1:3) * ANINT( TMP_CRD(1:3) / BOXL(1:3), DP )
!				TMP_CRD(1:3) = TMP_CRD(1:3) + QMGC(1:3)
!            ELSE
                TMP_CRD(1:3) = ATMCRD(1:3,I)
!            END IF
            K = ( J - 1 ) * 3
            EXT_CORD(K + 1) = TMP_CRD(1)
            EXT_CORD(K + 2) = TMP_CRD(2)
            EXT_CORD(K + 3) = TMP_CRD(3)
            EXT_ATMN(J) = ATMCHG(I)
        END IF
    END DO

    CODE = 0

    CALL CHRG_SERVER( CODE, NT, NATOMSQM, EXT_ATMN, EXT_CORD, EXT_ENER, EXT_QFIT, EXT_GRAD, EXT_HESS )

    ! . Map MK Charges into ATMCHG
    DO I = 1, NATOMSQM
        IF( .NOT. QMATOM(I)%QBOUNDARY ) ATMCHG(QMATOM(I)%ATOM) = EXT_QFIT(I)
    END DO

    ! . Accumulate Link-Atom charge on the QM counterpart... (not a good idea...)
    ! . Accumulate Link-Atom charge on the QM counterpart...
!    DO I = 1, NATOMSQM
!        IF( QMATOM(I)%QBOUNDARY ) ATMCHG(QMATOM(I)%PARTNER) = ATMCHG(QMATOM(I)%PARTNER) + EXT_QFIT(I)
!    END DO

	CALL SHOW_CHARGES

    ENERGY = ENERGY + EXT_ENER * AU_TO_KJ

    DEALLOCATE( SELE, EXT_ATMN, EXT_CORD, EXT_QFIT )
	DEALLOCATE( EXT_GRAD, EXT_HESS )
END SUBROUTINE

! ==============================================================================
! ==============================================================================

SUBROUTINE ABINITIO_GRADIENT( ENERGY, GRADIENT )
    IMPLICIT NONE

    REAL ( KIND = DP ), INTENT(INOUT)                          :: ENERGY
    REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT) :: GRADIENT

    REAL( KIND=DP )                            :: EXT_ENER, F1, TMP, GCI
    INTEGER                                    :: I, J, K, NT, NH, CODE, L
    LOGICAL, ALLOCATABLE, DIMENSION(:)         :: SELE
    REAL( KIND=DP ), DIMENSION(1:3)            :: TMP_CRD, QMGC
    REAL( KIND=DP ), DIMENSION(:), ALLOCATABLE :: EXT_CORD, EXT_ATMN, EXT_QFIT, EXT_GRAD, EXT_HESS

    ALLOCATE( SELE(1:NATOMS) )
    CALL ABINITIO_NON_BONDED( SELE )

    NT = NATOMSQM + COUNT( SELE )

    ALLOCATE( EXT_CORD(1:3*NT), EXT_ATMN(1:NT), EXT_QFIT(1:NATOMSQM), EXT_GRAD(1:3*NATOMSQM) )
	ALLOCATE( EXT_HESS(1) )

    DO I = 1, NATOMSQM
        EXT_ATMN(I) = QMATOM(I)%NUMBER * 1._DP
        TMP_CRD = MOPAC_DATA_ATOM_POSITION( I, ATMCRD )
        J = ( I - 1 ) * 3
        EXT_CORD(J + 1) = TMP_CRD(1)
        EXT_CORD(J + 2) = TMP_CRD(2)
        EXT_CORD(J + 3) = TMP_CRD(3)
		QMGC(1:3) = QMGC(1:3) + TMP_CRD(1:3)
    END DO
	QMGC(1:3) = QMGC(1:3) / REAL( NATOMSQM, KIND = DP )

    J = NATOMSQM
    DO I = 1, NATOMS
        IF( SELE(I) ) THEN
            J = J + 1
!            IF( QIMAGE ) THEN
!				TMP_CRD(1:3) = ATMCRD(1:3,I) - QMGC(1:3)
!                TMP_CRD(1:3) = TMP_CRD(1:3) - BOXL(1:3) * ANINT( TMP_CRD(1:3) / BOXL(1:3), DP )
!				TMP_CRD(1:3) = TMP_CRD(1:3) + QMGC(1:3)
!            ELSE
                TMP_CRD(1:3) = ATMCRD(1:3,I)
!            END IF
            K = ( J - 1 ) * 3
            EXT_CORD(K + 1) = TMP_CRD(1)
            EXT_CORD(K + 2) = TMP_CRD(2)
            EXT_CORD(K + 3) = TMP_CRD(3)
            EXT_ATMN(J) = ATMCHG(I)
        END IF
    END DO

    CODE = 1
    CALL CHRG_SERVER( CODE, NT, NATOMSQM, EXT_ATMN, EXT_CORD, EXT_ENER, EXT_QFIT, EXT_GRAD, EXT_HESS )

    ! . Map MK Charges into ATMCHG
    DO I = 1, NATOMSQM
        IF( .NOT. QMATOM(I)%QBOUNDARY ) ATMCHG(QMATOM(I)%ATOM) = EXT_QFIT(I)
    END DO
!    DO I = 1, NATOMSQM
!        IF( QMATOM(I)%QBOUNDARY ) ATMCHG(QMATOM(I)%PARTNER) = ATMCHG(QMATOM(I)%PARTNER) + EXT_QFIT(I)
!    END DO

	CALL SHOW_CHARGES

! ------------------------------------------------------------------------------------------------------
DO I = 1, NATOMSQM
	IF( QMATOM(I)%QBOUNDARY ) THEN
		J = ( ATMQMI( QMATOM(I)%PARTNER ) - 1 ) * 3
		K = ( I - 1 ) * 3
		DO L = 1, 3
			TMP_CRD(L) = EXT_CORD(J+L) - EXT_CORD(K+L)
		END DO
		TMP = SQRT( SUM( TMP_CRD ** 2 ) )
		GCI = .0_DP
		DO L = 1, 3
			GCI = GCI + TMP_CRD(L) * ( EXT_GRAD(J+L) - EXT_GRAD(K+L) )
		END DO
		GCI = GCI / TMP * 0.5_DP
		DO L = 1, 3
			EXT_GRAD(J+L) = EXT_GRAD(J+L) - GCI * TMP_CRD(L) / TMP
		END DO
	END IF
END DO
! ------------------------------------------------------------------------------------------------------

	! . Map gradient vector into ATMDER
    F1 = AU_TO_KJ / BOHRS_TO_ANGSTROMS
    DO I = 1, NATOMSQM
        J = ( I - 1 ) * 3
        IF( .NOT. QMATOM(I)%QBOUNDARY ) THEN
            DO K = 1, 3
                GRADIENT(K,QMATOM(I)%ATOM) = GRADIENT(K,QMATOM(I)%ATOM) + EXT_GRAD(J + K) * F1
            END DO
        END IF
    END DO

    ENERGY = ENERGY + EXT_ENER * AU_TO_KJ

    DEALLOCATE( SELE, EXT_ATMN, EXT_CORD, EXT_QFIT, EXT_GRAD )
	DEALLOCATE( EXT_HESS )
END SUBROUTINE

! ==============================================================================
! ==============================================================================

SUBROUTINE ABINITIO_HESSIAN( ENERGY, GRADIENT, HESSIAN )
    IMPLICIT NONE

    REAL ( KIND = DP ), INTENT(INOUT)                          :: ENERGY
    REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT) :: GRADIENT
    REAL ( KIND = DP ), DIMENSION(1:*), INTENT(INOUT)          :: HESSIAN
!    REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT) :: HESSIAN

    REAL( KIND=DP )                            :: EXT_ENER, F1, F2, TMP, GCI
    INTEGER                                    :: I, J, K, L, M, NT, NH, NF, CODE, KI, KJ
    LOGICAL, ALLOCATABLE, DIMENSION(:)         :: SELE
    REAL( KIND=DP ), DIMENSION(1:3)            :: TMP_CRD, QMGC
    REAL( KIND=DP ), DIMENSION(:), ALLOCATABLE :: EXT_CORD, EXT_ATMN, EXT_QFIT, EXT_GRAD, EXT_HESS

    ALLOCATE( SELE(1:NATOMS) )
    CALL ABINITIO_NON_BONDED( SELE )

    NT = NATOMSQM + COUNT( SELE )
    NH = 3*NATOMSQM*(3*NATOMSQM+1)/2
    NF = 3*NFREE*(3*NFREE+1)/2

    ALLOCATE( EXT_CORD(1:3*NT), EXT_ATMN(1:NT), EXT_QFIT(1:NATOMSQM), EXT_GRAD(1:3*NATOMSQM), EXT_HESS(1:NH) )

    DO I = 1, NATOMSQM
        EXT_ATMN(I) = QMATOM(I)%NUMBER * 1._DP
        TMP_CRD = MOPAC_DATA_ATOM_POSITION( I, ATMCRD )
        J = ( I - 1 ) * 3
        EXT_CORD(J + 1) = TMP_CRD(1)
        EXT_CORD(J + 2) = TMP_CRD(2)
        EXT_CORD(J + 3) = TMP_CRD(3)
		QMGC(1:3) = QMGC(1:3) + TMP_CRD(1:3)
    END DO
	QMGC(1:3) = QMGC(1:3) / REAL( NATOMSQM, KIND = DP )

    J = NATOMSQM
    DO I = 1, NATOMS
        IF( SELE(I) ) THEN
            J = J + 1
!            IF( QIMAGE ) THEN
!				TMP_CRD(1:3) = ATMCRD(1:3,I) - QMGC(1:3)
!                TMP_CRD(1:3) = TMP_CRD(1:3) - BOXL(1:3) * ANINT( TMP_CRD(1:3) / BOXL(1:3), DP )
!				TMP_CRD(1:3) = TMP_CRD(1:3) + QMGC(1:3)
!            ELSE
                TMP_CRD(1:3) = ATMCRD(1:3,I)
!            END IF
            K = ( J - 1 ) * 3
            EXT_CORD(K + 1) = TMP_CRD(1)
            EXT_CORD(K + 2) = TMP_CRD(2)
            EXT_CORD(K + 3) = TMP_CRD(3)
            EXT_ATMN(J) = ATMCHG(I)
        END IF
    END DO

    CODE = 2
    CALL CHRG_SERVER( CODE, NT, NATOMSQM, EXT_ATMN, EXT_CORD, EXT_ENER, EXT_QFIT, EXT_GRAD, EXT_HESS )

    ! . Map MK Charges into ATMCHG
    DO I = 1, NATOMSQM
        IF( .NOT. QMATOM(I)%QBOUNDARY ) ATMCHG(QMATOM(I)%ATOM) = EXT_QFIT(I)
    END DO
    ! . Accumulate Link-Atom charge on the QM counterpart...
!    DO I = 1, NATOMSQM
!        IF( QMATOM(I)%QBOUNDARY ) ATMCHG(QMATOM(I)%PARTNER) = ATMCHG(QMATOM(I)%PARTNER) + EXT_QFIT(I)
!    END DO

	CALL SHOW_CHARGES

! ------------------------------------------------------------------------------------------------------
DO I = 1, NATOMSQM
	IF( QMATOM(I)%QBOUNDARY ) THEN
		J = ( ATMQMI( QMATOM(I)%PARTNER ) - 1 ) * 3
		K = ( I - 1 ) * 3
		DO L = 1, 3
			TMP_CRD(L) = EXT_CORD(J+L) - EXT_CORD(K+L)
		END DO
		TMP = SQRT( SUM( TMP_CRD ** 2 ) )
		GCI = .0_DP
		DO L = 1, 3
			GCI = GCI + TMP_CRD(L) * ( EXT_GRAD(J+L) - EXT_GRAD(K+L) )
		END DO
		GCI = GCI / TMP * 0.5_DP
		DO L = 1, 3
			EXT_GRAD(J+L) = EXT_GRAD(J+L) - GCI * TMP_CRD(L) / TMP
		END DO
	END IF
END DO
! ------------------------------------------------------------------------------------------------------

	! . Map gradient vector into ATMDER
    F1 = AU_TO_KJ / BOHRS_TO_ANGSTROMS
    DO I = 1, NATOMSQM
        J = ( I - 1 ) * 3
        IF( .NOT. QMATOM(I)%QBOUNDARY ) THEN
            DO K = 1, 3
                GRADIENT(K,QMATOM(I)%ATOM) = GRADIENT(K,QMATOM(I)%ATOM) + EXT_GRAD(J + K) * F1
            END DO
        END IF
    END DO

    ! . Map Hessian matrix (truncated for link atoms)
    F2 = F1 / BOHRS_TO_ANGSTROMS
    DO I = 1, NF
	    HESSIAN(I) = HESSIAN(I) + EXT_HESS(I) * F2
    END DO

    ENERGY = ENERGY + EXT_ENER * AU_TO_KJ

    DEALLOCATE( SELE, EXT_ATMN, EXT_CORD, EXT_QFIT, EXT_GRAD, EXT_HESS )
END SUBROUTINE

! ==============================================================================
! ==============================================================================

SUBROUTINE ABINITIO_EXIT
    IMPLICIT NONE
END SUBROUTINE

END MODULE
