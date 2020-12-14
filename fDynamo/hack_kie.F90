MODULE KIE

#ifdef   HACK_KIE

! .Module declarations.
USE CONSTANTS
USE DEFINITIONS,        ONLY : DP
USE ATOMS,              ONLY : NATOMS, ATMCRD, ATMMAS, ATMFIX, NFREE
USE POTENTIAL_ENERGY,   ONLY : ATMHES
USE NORMAL_MODE,        ONLY : NORMAL_MODE_FREQUENCIES, NMODES, FREQUENCIES
USE TRANSFORMATION,     ONLY : MOMENTS_OF_INERTIA_AT_CENTER

IMPLICIT NONE
PUBLIC
SAVE

!===============================================================================
CONTAINS
!===============================================================================

   !------------------------------------------------------------------
   SUBROUTINE GIBBS_ENERGY( PRESS, TEMP, SYMM, NMSKP, GIBBS, WIGNER )
   !------------------------------------------------------------------
   REAL( KIND=DP ), INTENT( IN ) :: PRESS, TEMP, SYMM
   INTEGER, INTENT( IN ) :: NMSKP
   REAL( KIND=DP ), INTENT( INOUT ) :: GIBBS
   REAL( KIND=DP ), INTENT( INOUT ), OPTIONAL :: WIGNER

   INTEGER :: I, J
   REAL( KIND=DP ) :: MT, QT, QR, QV, GG, ZZ, KK
   REAL( KIND=DP ), DIMENSION(1:3) :: MI
   REAL( KIND=DP ), DIMENSION(:,:), ALLOCATABLE :: CV
   REAL( KIND=DP ), DIMENSION(:), ALLOCATABLE :: MV

   CALL NORMAL_MODE_FREQUENCIES( ATMHES, "PROJECT", PRINT = .TRUE. )

   ALLOCATE( CV(1:3,1:NFREE), MV(1:NFREE) )
   J = 1
   DO I = 1, NATOMS
       IF( .NOT. ATMFIX(I) ) THEN
           MV(J)     = ATMMAS(I)
           CV(1:3,J) = ATMCRD(1:3,I)
           J         = J + 1
       END IF
   END DO
   MT = SUM( MV ) * AMU_TO_KG
   WRITE(*,"(/A,G20.10)") "MASS (KG):                   ", MT

   ! -- TRANSLATIONAL COMPONENT (DIVIDED BY NAVOGADRO)
   QT = ( ( 2.0D0 * PI * MT * KBOLTZ * TEMP / ( PLANCK * PLANCK ) )**1.5 ) * KBOLTZ * TEMP / ( PRESS * ATM_TO_PASCALS )
   QT = DLOG( QT )
   WRITE(*,"(A,F20.4)") "LN( QTRA ):                  ", QT

   ! -- ROTATIONAL COMPONENT
   CALL MOMENTS_OF_INERTIA_AT_CENTER( CV, MI, WEIGHTS = MV )
   MI = MI * AMU_TO_KG * 1.0D-20
   KK = ( 8.0D0 * PI * PI * KBOLTZ * TEMP ) / ( PLANCK * PLANCK )
   QR = SQRT ( PI * KK**3 * PRODUCT( MI ) ) / SYMM
   QR = DLOG( QR )
   WRITE(*,"(A,F20.4)") "LN( QROT ):                  ", QR


   ! - VIBRATIONAL COMPONENT
   KK = 100.D0 * C * PLANCK / ( KBOLTZ * TEMP )
   IF( PRESENT( WIGNER ) ) WIGNER = 1.0D0 + 1.0D0 / 24.0D0 * ( FREQUENCIES(1) * KK ) ** 2
   QV = 1.0D0
   I = COUNT( FREQUENCIES < 10.0D0 )
   IF( NMSKP /= I ) &
   WRITE(*,"(A,I0,A,I0,A)") "Warning: Told to skip ", NMSKP, " frequencies, but found ", i, " (< 10 cm-1)."
       
   DO I = NMSKP + 1, NMODES
      QV = QV / ( 1.0D0 - DEXP( - FREQUENCIES(I) * KK ) )
!   WRITE(*,"(I6,F12.3,2G20.10)") I, FREQUENCIES(I), 1.0D0 / ( 1.0D0 - DEXP( - FREQUENCIES(I) * KK ) ), QV
   END DO
   QV = DLOG( QV )
   WRITE(*,"(A,F20.4)") "LN( QVIB ):                  ", QV

   ! - GIBBS CONTRIBUTION
   GG = - R * TEMP * ( QT + QR + QV )
   WRITE(*,"(A,F20.4)") "GIBBS CONTRIBUTION (KJ/MOL): ", GG

   ! - ZPE
   ZZ = SUM( FREQUENCIES(NMSKP+1:NMODES) ) * 0.5D0 * 100.D0 * C * PLANCK * NAVOGADRO * 1.0D-3
   WRITE(*,"(A,F20.4)") "ZPE (KJ/MOL):                ", ZZ

   GIBBS = GG + ZZ
   ! - G = F + PV = -RT Ln Q + nRT					(cancels out when calculating relative terms)
!   GIBBS = GIBBS + R * TEMP * DLOG( NAVOGADRO )
   WRITE(*,"(A,F20.4)") "TOTAL (KJ/MOL):              ", GIBBS

   DEALLOCATE( CV, MV )

   END SUBROUTINE GIBBS_ENERGY

#endif

END MODULE KIE
