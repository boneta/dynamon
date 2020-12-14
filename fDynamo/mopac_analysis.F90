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
!                        The MOPAC Energy Analysis Module
!===============================================================================
!
! . Subroutines:
!
!   MOPAC_ENERGY_ANALYSIS                Analyse a MOPAC energy.
!   MOPAC_ENERGY_ANALYSIS_CLEAR          Initialize the analysis arrays.
!   MOPAC_ENERGY_ANALYSIS_PRINT          Selectively print the analysis.
!
!   MOPAC_IONIZATION_POTENTIAL           Determine the ionization potential.
!
! . Notes:
!
!   This module needs to be updated to cope with PI atoms.
!
!===============================================================================
MODULE MOPAC_ANALYSIS

! . Module declarations.
USE CONSTANTS,          ONLY : ANGSTROMS_TO_BOHRS, AU_TO_EV, EV_TO_KJ, KCAL_TO_KJ
USE DEFINITIONS,        ONLY : DP
USE DIAGONALIZATION,    ONLY : SYMMETRIC_UPPER
USE PRINTING,           ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, &
                               PRINT_SUMMARY_STOP, PRINT_PARAGRAPH

USE ATOMS,              ONLY : ATMCRD, ATMQMI, NATOMS, NATOMSQM
USE ENERGY_NON_BONDING, ONLY : CUT_OFF, CUT_ON, ENERGY_NON_BONDING_INTERACTION, ENERGY_NON_BONDING_SELF, &
                               NBLIST_TYPE, NBLISTQM_FIRST, NNBLISTQM, QIMAGE
USE MM_TERMS,           ONLY : ATMCHG
USE MOPAC_DATA,         ONLY : ATHEAT, DENMAT, DENMATA, DENMATB, HAMILTONIAN, HCORE, MOPAC_DATA_ATOM_POSITION, &
                               NALPHA, NBASIS, NBASTR, NBETA, NPIATOMS, NPIBEADS, QMATOM, SETEI
USE MOPAC_FOCK_MATRIX,  ONLY : FOCK_MATRIX_RESTRICTED, FOCK_MATRIX_UNRESTRICTED
USE MOPAC_INTEGRALS,    ONLY : INTEGRALS_HCORE, SEINTC_INTEGRALAB, SEINTC_INTEGRALB, SEINTC_SETUP
USE MOPAC_PARAMETERS
USE SYMMETRY,           ONLY : BOXL

IMPLICIT NONE
PRIVATE
PUBLIC :: MOPAC_ENERGY_ANALYSIS, MOPAC_ENERGY_ANALYSIS_CLEAR, MOPAC_ENERGY_ANALYSIS_PRINT, MOPAC_IONIZATION_POTENTIAL
#ifndef PGPC
SAVE
#endif

! . Module scalars.
LOGICAL :: QANALYSIS = .FALSE.

! . Module arrays.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: EA1, EA2
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: EEEREPULSION, EENATTRACTION, EEXCHANGE, ENNREPULSION, ERESONANCE
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: EQMMMEL

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------------------
   SUBROUTINE MOPAC_ENERGY_ANALYSIS ( PRINT )
   !-----------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Local scalars.
   INTEGER            :: NATQMTR
   LOGICAL            :: QAM1PM3, QPRINT
   REAL ( KIND = DP ) :: EABE, EABEE, EABEN, EABNN, EABR, EABRX, EABX, EAE, EAU, EQMMM, EQMP, ET, TONE, TTWO

   ! . Initialize the data structure.
   CALL MOPAC_ENERGY_ANALYSIS_CLEAR

   ! . Check that there are QM atoms and the appropriate arrays are present.
   IF ( ( NATOMSQM <= 0 )              .OR. &
        ( .NOT. ALLOCATED ( DENMAT ) ) .OR. &
        ( .NOT. ALLOCATED ( HCORE  ) ) .OR. &
        ( .NOT. ALLOCATED ( SETEI  ) ) ) RETURN

   ! . Check that there are no PI atoms.
   IF ( NPIATOMS > 0 ) RETURN

   ! . Check for a UHF calculation.
   IF ( ALLOCATED ( DENMATA ) .OR. ALLOCATED ( DENMATB ) ) THEN
      CALL PRINT_ERROR ( "MOPAC_ENERGY_ANALYSIS", "Unable to analyse a UHF wavefunction at present." )
   END IF

   ! . Set the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Set the analysis flag.
   QANALYSIS = .TRUE.

   ! . Allocate the analysis arrays.
   NATQMTR = ( NATOMSQM * ( NATOMSQM + 1 ) ) / 2
   ALLOCATE ( EA1(1:NATOMSQM), EA2(1:NATOMSQM), EEEREPULSION(1:NATQMTR), EENATTRACTION(1:NATQMTR), EEXCHANGE(1:NATQMTR), &
                                            ENNREPULSION(1:NATQMTR), ERESONANCE(1:NATQMTR), EQMMMEL(1:NATOMSQM,1:NATOMS) )

   ! . Set the AM1PM3 flag.
   QAM1PM3 = ( HAMILTONIAN == "AM1 " ) .OR. ( HAMILTONIAN == "PDDG" ) .OR. ( HAMILTONIAN == "PM3 " ) .OR. ( HAMILTONIAN == "RM1 " )

   ! . Do the QM energy partitioning.
   CALL QMPARTITION

   ! . Do the QM/MM electrostatic energy partitioning.
   CALL QMPARTITION_ELECTROSTATIC

   ! . Print out the results if requested.
   IF ( QPRINT ) THEN

      ! . Determine the totals of the one-electron terms.
      EAU   = SUM ( EA1 )
      EAE   = SUM ( EA2 )
      TONE  = EAU + EAE

      ! . Determine the totals of the two-electron terms.
      EABR  = SUM ( ERESONANCE    )
      EABX  = SUM ( EEXCHANGE     )
      EABEE = SUM ( EEEREPULSION  )
      EABEN = SUM ( EENATTRACTION )
      EABNN = SUM ( ENNREPULSION  )
      EABRX = EABR  + EABX
      EABE  = EABEE + EABEN + EABNN
      TTWO  = EABRX + EABE

      ! . Calculate the total QM/MM electrostatic energy.
      EQMMM = SUM ( EQMMMEL )

      ! . Calculate the total energies.
      ET    = TONE  + TTWO + EQMMM
      EQMP  = EV_TO_KJ * ET + ATHEAT

      ! . Write out the energy terms.
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800", PAGEWIDTH = 112, VARIABLEWIDTH = 20 )
      CALL PRINT_SUMMARY_START ( "Summary of QM Energy Partitioning (eV)" )
      WRITE ( PRINT_LINE, "(F20.4)" ) EAU   ; CALL PRINT_SUMMARY_ELEMENT ( "Electron-Nuclear  (One-Electron)" )
      WRITE ( PRINT_LINE, "(F20.4)" ) EAE   ; CALL PRINT_SUMMARY_ELEMENT ( "Electron-Electron (Two-Electron)" )
      WRITE ( PRINT_LINE, "(F20.4)" ) TONE  ; CALL PRINT_SUMMARY_ELEMENT ( "Total One-Center Terms"           )
      WRITE ( PRINT_LINE, "(F20.4)" ) EABR  ; CALL PRINT_SUMMARY_ELEMENT ( "Resonance Energy"                 )
      WRITE ( PRINT_LINE, "(F20.4)" ) EABX  ; CALL PRINT_SUMMARY_ELEMENT ( "Exchange Energy"                  )
      WRITE ( PRINT_LINE, "(F20.4)" ) EABRX ; CALL PRINT_SUMMARY_ELEMENT ( "Exchange+Resonance Energy"        )
      WRITE ( PRINT_LINE, "(F20.4)" ) EABEE ; CALL PRINT_SUMMARY_ELEMENT ( "Electron-Electron Repulsion"      )
      WRITE ( PRINT_LINE, "(F20.4)" ) EABEN ; CALL PRINT_SUMMARY_ELEMENT ( "Electron-Nuclear Attraction"      )
      WRITE ( PRINT_LINE, "(F20.4)" ) EABNN ; CALL PRINT_SUMMARY_ELEMENT ( "Nuclear-Nuclear Repulsion"        )
      WRITE ( PRINT_LINE, "(F20.4)" ) EABE  ; CALL PRINT_SUMMARY_ELEMENT ( "Total Electrostatic Interaction"  )
      WRITE ( PRINT_LINE, "(F20.4)" ) TTWO  ; CALL PRINT_SUMMARY_ELEMENT ( "Total Two-Center Terms"           )
      WRITE ( PRINT_LINE, "(F20.4)" ) EQMMM ; CALL PRINT_SUMMARY_ELEMENT ( "QM/MM Energy"                     )
      WRITE ( PRINT_LINE, "(F20.4)" ) ET    ; CALL PRINT_SUMMARY_ELEMENT ( "Total Energy"                     )
      WRITE ( PRINT_LINE, "(F20.4)" ) EQMP  ; CALL PRINT_SUMMARY_ELEMENT ( "Total Energy (kJ mol^-1)"         )
      CALL PRINT_SUMMARY_STOP

   END IF

   !===========================================================================
   CONTAINS
   !===========================================================================

      !---------------------
      SUBROUTINE QMPARTITION
      !---------------------

      ! . Local scalars.
      INTEGER            :: I, IA, IAP1, IA1, IA2, IB, IG, II, IJ, IK, IL, IMINUS, &
                            ISS, ISX, ISY, ISZ, IXX, IXY, IXZ, IYY, IYZ, IZZ,      &
                            J, JA, JAP1, JB, JJ, JK, JL, JSS, K, KA, KB, KC, KINC, &
                            KK, KL, L, N, NI, NJ, NT
      REAL ( KIND = DP ) :: AA, BB, G, ONEII, ONEJJ, PIJ, R, SCALE, SS1, SS2, SS3, SS4, SS5, T, TT1, TT2, TT3, TT4, TT5

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3)      :: RI, RJ
      REAL ( KIND = DP ), DIMENSION(1:NBASTR) :: ALPHA, BETA

      !-------------------------------------------------------------------------
      ! . Calculate the one-center and two-center energies.
      !-------------------------------------------------------------------------
      ! . Define the ALPHA and BETA densities.
      ALPHA = 0.5_DP * DENMAT(:,1)
      BETA  = 0.5_DP * DENMAT(:,1)

      ! . Initialization.
      EA1           = 0.0_DP
      EA2           = 0.0_DP
      EEEREPULSION  = 0.0_DP
      EENATTRACTION = 0.0_DP
      EEXCHANGE     = 0.0_DP
      ENNREPULSION  = 0.0_DP
      ERESONANCE    = 0.0_DP

      ! . One-centre energies.
      K = 0
      DO I = 1,NATOMSQM

         ! . Get some information for the atom.
         IA = QMATOM(I)%BFIRST
         IB = QMATOM(I)%BLAST
         NI = QMATOM(I)%NUMBER

         ! . Loop over the orbitals.
         DO J = IA,IB
            K = K + J
            T = UPP(NI)
            IF ( J == IA ) T = USS(NI)
            EA1(I) = EA1(I) + DENMAT(K,1) * T
         END DO

         ISS = ( IA * ( IA + 1 ) ) / 2
         EA2(I) = 0.5_DP * GSS(NI) * DENMAT(ISS,1) * DENMAT(ISS,1) - &
                  0.5_DP * GSS(NI) * ( ALPHA(ISS) * ALPHA(ISS) + BETA(ISS) * BETA(ISS) )
         IF ( IA == IB ) CYCLE

         IA1 = IA+1
         IA2 = IA+2
         IXX = IA1*IA2/2
         IYY = IA2*IB/2
         IZZ = (IB*(IB+1))/2
         IXY = IA1+IA2*IA1/2
         IXZ = IA1+IB*IA2/2
         IYZ = IA2+IB*IA2/2
         ISX = IA+IA1*IA/2
         ISY = IA+IA2*IA1/2
         ISZ = IA+IB*IA2/2
         SS1 = DENMAT(IXX,1) *   DENMAT(IXX,1) + DENMAT(IYY,1) * DENMAT(IYY,1) + DENMAT(IZZ,1) * DENMAT(IZZ,1)
         SS2 = DENMAT(ISS,1) * ( DENMAT(IXX,1) + DENMAT(IYY,1) + DENMAT(IZZ,1) )
         SS3 = DENMAT(IXX,1) *   DENMAT(IYY,1) + DENMAT(IXX,1) * DENMAT(IZZ,1) + DENMAT(IYY,1) * DENMAT(IZZ,1)
         SS4 = DENMAT(ISX,1) *   DENMAT(ISX,1) + DENMAT(ISY,1) * DENMAT(ISY,1) + DENMAT(ISZ,1) * DENMAT(ISZ,1)
         SS5 = DENMAT(IXY,1) *   DENMAT(IXY,1) + DENMAT(IXZ,1) * DENMAT(IXZ,1) + DENMAT(IYZ,1) * DENMAT(IYZ,1)
         TT1 = ALPHA(IXX) * ALPHA(IXX) + ALPHA(IYY) * ALPHA(IYY) + ALPHA(IZZ) * ALPHA(IZZ) + &
               BETA(IXX)  * BETA(IXX)  + BETA(IYY)  * BETA(IYY)  + BETA(IZZ)  * BETA(IZZ)
         TT2 = ALPHA(ISS) * ( ALPHA(IXX) + ALPHA(IYY) + ALPHA(IZZ) ) + &
               BETA(ISS)  * ( BETA(IXX)  + BETA(IYY)  + BETA(IZZ)  )
         TT3 = ALPHA(IXX) * ALPHA(IYY) + ALPHA(IXX) * ALPHA(IZZ) + ALPHA(IYY) * ALPHA(IZZ) + &
               BETA(IXX)  * BETA(IYY)  + BETA(IXX)  * BETA(IZZ)  + BETA(IYY)  * BETA(IZZ)
         TT4 = ALPHA(ISX) * ALPHA(ISX) + ALPHA(ISY) * ALPHA(ISY) + ALPHA(ISZ) * ALPHA(ISZ) + &
               BETA(ISX)  * BETA(ISX)  + BETA(ISY)  * BETA(ISY)  + BETA(ISZ)  * BETA(ISZ)
         TT5 = ALPHA(IXY) * ALPHA(IXY) + ALPHA(IXZ) * ALPHA(IXZ) + ALPHA(IYZ) * ALPHA(IYZ) + &
               BETA(IXY)  * BETA(IXY)  + BETA(IXZ)  * BETA(IXZ)  + BETA(IYZ)  * BETA(IYZ)

         EA2(I) = EA2(I) + 0.5_DP*GPP(NI)*SS1 + GSP(NI)*SS2 + GP2(NI)*SS3           + &
                           HSP(NI)*SS4*2.0_DP + 0.5_DP*(GPP(NI)-GP2(NI))*SS5*2.0_DP - &
                           0.5_DP*GPP(NI)*TT1 - GSP(NI)*TT4 - GP2(NI)*TT5           - &
                           HSP(NI)*(TT2+TT4)  - 0.5_DP*(GPP(NI)-GP2(NI))*(TT3+TT5)

      END DO

      ! . Two-centre energies.
      ! . Resonance terms.
      N=1
      DO II = 2,NATOMSQM
         IA     = QMATOM(II)%BFIRST
         IB     = QMATOM(II)%BLAST
         IMINUS = II - 1
         ONEII  = 1.0_DP
         DO JJ = 1,IMINUS
            N     = N + 1
            JA    = QMATOM(JJ)%BFIRST
            JB    = QMATOM(JJ)%BLAST
            ONEJJ = 1.0_DP
            DO I = IA,IB
               KA = (I*(I-1))/2
               DO K = JA,JB
                  IK = KA + K
                  ERESONANCE(N) = ERESONANCE(N) + 2.0_DP * DENMAT(IK,1) * HCORE(IK) * ONEII * ONEJJ
               END DO
            END DO
         END DO
         N=N+1
      END DO

      ! . Core-core repulsion and core-electron attraction.
      N  = 1
      KK = 0
      DO II = 2,NATOMSQM

         ! . Get some information for the atom.
         IA     = QMATOM(II)%BFIRST
         IB     = QMATOM(II)%BLAST
         NI     = QMATOM(II)%NUMBER
         ISS    = (IA*(IA+1))/2
         IMINUS = II-1
         RI     = MOPAC_DATA_ATOM_POSITION ( II, ATMCRD )

         DO JJ = 1,IMINUS

            N = N + 1

            ! . Get some information for the atom.
            JA    = QMATOM(JJ)%BFIRST
            JB    = QMATOM(JJ)%BLAST
            NJ    = QMATOM(JJ)%NUMBER
            JSS   = (JA*(JA+1))/2
            KK    = KK+1
            RJ     = MOPAC_DATA_ATOM_POSITION ( JJ, ATMCRD )

            G     = SETEI(KK)
            R     = SQRT ( SUM ( ( RI - RJ )**2 ) )
            SCALE = 1.0_DP + EXP(-ALP(NI)*R) + EXP(-ALP(NJ)*R)
            NT    = NI + NJ
            IF ( ( NT == 8 ) .OR. ( NT == 9 ) ) THEN
               IF ( ( NI == 7 ) .OR. ( NI == 8 ) ) SCALE = SCALE + (R-1.0_DP)*EXP(-ALP(NI)*R)
               IF ( ( NJ == 7 ) .OR. ( NJ == 8 ) ) SCALE = SCALE + (R-1.0_DP)*EXP(-ALP(NJ)*R)
            END IF
            ENNREPULSION(N) = CORE(NI) * CORE(NJ) * G * SCALE

            IF ( QAM1PM3 ) THEN
               SCALE = 0.0_DP
               DO IG = 1,4
                  IF ( ABS ( FN1(NI,IG) ) > 0.0_DP ) SCALE = SCALE + CORE(NI)*CORE(NJ)/R * &
                                                                     FN1(NI,IG)*EXP(-FN2(NI,IG)*(R-FN3(NI,IG))**2)
                  IF ( ABS ( FN1(NJ,IG) ) > 0.0_DP ) SCALE = SCALE + CORE(NI)*CORE(NJ)/R * &
                                                                     FN1(NJ,IG)*EXP(-FN2(NJ,IG)*(R-FN3(NJ,IG))**2)
               END DO
               ENNREPULSION(N) = ENNREPULSION(N) + SCALE
            END IF

            EENATTRACTION(N) = - ( DENMAT(ISS,1)*CORE(NJ) + DENMAT(JSS,1)*CORE(NI) ) * G

            IF ( NJ >= 3 ) THEN
               KINC = 9
               JAP1 = JA+1
               DO K=JAP1,JB
                  KC = (K*(K-1))/2
                  DO L = JA,K
                     KL = KC+L
                     BB = 2.0_DP
                     IF ( K == L ) BB = 1.0_DP
                     KK = KK+1
                     EENATTRACTION(N) = EENATTRACTION(N) - DENMAT(KL,1)*CORE(NI)*BB*SETEI(KK)
                  END DO
               END DO
            ELSE
               KINC = 0
            END IF
            IF ( NI >= 3 ) THEN
               IAP1 = IA + 1
               DO I=IAP1,IB
                  KA = (I*(I-1))/2
                  DO J = IA,I
                     IJ = KA+J
                     AA = 2.0_DP
                     IF ( I == J ) AA = 1.0_DP
                     KK = KK+1
                     EENATTRACTION(N) = EENATTRACTION(N) - DENMAT(IJ,1)*CORE(NJ)*AA*SETEI(KK)
                     KK = KK + KINC
                  END DO
               END DO
            END IF

         END DO
         N = N + 1
      END DO

      ! . Coulomb and exchange terms.
      N  = 1
      KK = 0
      DO II = 2,NATOMSQM
         IA     = QMATOM(II)%BFIRST
         IB     = QMATOM(II)%BLAST
         IMINUS = II - 1
         DO JJ = 1,IMINUS
            JA = QMATOM(JJ)%BFIRST
            JB = QMATOM(JJ)%BLAST
            N  = N + 1
            DO I = IA,IB
               KA = (I*(I-1))/2
               DO J = IA,I
                  KB = (J*(J-1))/2
                  IJ = KA + J
                  AA = 2.0_DP
                  IF ( I == J ) AA=1.0_DP
                  PIJ = DENMAT(IJ,1)
                  DO K = JA,JB
                     KC = (K*(K-1))/2
                     IK = KA+K
                     JK = KB+K
                     DO L = JA,K
                        IL = KA+L
                        JL = KB+L
                        KL = KC+L
                        BB = 2.0_DP
                        IF ( K == L ) BB = 1.0_DP
                        KK = KK+1
                        G  = SETEI(KK)
                        EEEREPULSION(N) = EEEREPULSION(N) + AA * BB * G * PIJ * DENMAT(KL,1)
                        EEXCHANGE(N)    = EEXCHANGE(N) - 0.5_DP * AA * BB * G * &
                                                       ( ALPHA(IK) * ALPHA(JL) + ALPHA(IL) * ALPHA(JK) + &
                                                         BETA(IK)  * BETA(JL)  + BETA(IL)  * BETA(JK) )
                     END DO
                  END DO
               END DO
            END DO
         END DO
         N = N + 1
      END DO

      END SUBROUTINE QMPARTITION

      !-----------------------------------
      SUBROUTINE QMPARTITION_ELECTROSTATIC
      !-----------------------------------

      ! . Local variables.
      INTEGER            :: I, IATOM, II, IINT, IREC, J, JATOM, JJ, N1, N2, MATOM, NQM, NTERM, QATOM, QINDX
      REAL ( KIND = DP ) :: ALPHA, CGQM, D1, D2, D4, ENB, F1, F2, F3, GAMMA, RHO0, RHO1, RHO2, R2OFF, R2ON, SWITCH, XQ, YQ, ZQ

      ! . Local integral scalars.
      REAL ( KIND = DP ) :: SS, SPZ, PZPZ, PPPP, DSS, DSPZ, DPZPZ, DPPPP,         & 
                            SFCT1M,   SFCT1P,  SFCT20,  SFCT2M,  SFCT2P,  SFCT4P, &
                            DSFCT1M, DSFCT1P, DSFCT20, DSFCT2M, DSFCT2P, DSFCT4P

      ! . Local intermediate scalars.
      REAL ( KIND = DP ) :: CGL, ENUC, RQM, RQM2, RQMB, RQMI, TEMP1, TEMP2, TEMP3, TEMP4, &
                            XN, YN, ZN, XN2, YN2, ZN2, XQM, YQM, ZQM

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3)  :: RI, RQ
      REAL ( KIND = DP ), DIMENSION(1:10) :: DENFAC, E1B, PFAC

      ! . Local types.
      TYPE(NBLIST_TYPE), POINTER :: NBCURRENT, NBNEXT

      ! . QM/MM parameters.
      REAL ( KIND = DP ), PARAMETER :: ALPMM = 5.0_DP, RHO0MM = 0.0_DP

      ! . Calculate some factors.
      R2OFF = CUT_OFF**2
      R2ON  = CUT_ON**2
      GAMMA = ( R2OFF - R2ON )**3
      CALL SEINTC_SETUP

      ! . Fill the array of density factors (off-diagonal elements * 2).
      DENFAC     = 2.0_DP
      DENFAC(1)  = 1.0_DP
      DENFAC(3)  = 1.0_DP
      DENFAC(6)  = 1.0_DP
      DENFAC(10) = 1.0_DP

      ! . Initialization.
      EQMMMEL = 0.0_DP

      ! . Loop over the records in the non-bonding interaction lists.
      DO IREC = 1,NNBLISTQM

	 ! . Get the next list.
	 IF ( IREC == 1 ) THEN
            NBCURRENT => NBLISTQM_FIRST
            NBNEXT    => NBLISTQM_FIRST%NEXT_LIST
	 ELSE
            NBCURRENT => NBNEXT
            NBNEXT    => NBNEXT%NEXT_LIST
	 END IF

         ! . Get the first atom of the interaction.
	 IATOM = NBCURRENT%ATOM

         ! . Get the position vector for IATOM if it is a quantum or a boundary atom.
	 IF ( ATMQMI(IATOM) > 0 ) RI = MOPAC_DATA_ATOM_POSITION ( ATMQMI(IATOM), ATMCRD )

         ! . Loop over the interactions for the atom.
         DO IINT = 1,SIZE(NBCURRENT%INTERACTIONS)

            ! . Get the second atom of the interaction.
            JATOM = NBCURRENT%INTERACTIONS(IINT)

            ! . Determine which is the QM atom. A negative J indicates J is quantum.
            IF ( JATOM > 0 ) THEN
               QATOM =   IATOM
               MATOM =   JATOM
	       QINDX =   ATMQMI(QATOM)
	       RQ    =   RI
            ELSE
               QATOM = - JATOM
               MATOM =   IATOM
	       QINDX =   ATMQMI(QATOM)
               RQ    =   MOPAC_DATA_ATOM_POSITION ( QINDX, ATMCRD )
            END IF

            ! . Define some constants for the QM atom.
            NQM   = QMATOM(QINDX)%NUMBER
            ALPHA = ALP(NQM)
            CGQM  = CORE(NQM)
            N1    = QMATOM(QINDX)%BFIRST
            N2    = QMATOM(QINDX)%BLAST
            XQ    = RQ(1)
            YQ    = RQ(2)
            ZQ    = RQ(3)
            RHO0  = ( 0.5_DP / AM(NQM) + RHO0MM )**2

            IF ( QAM1PM3 ) NTERM = COUNT ( ABS ( FN1(NQM,1:4) ) > 0.0_DP )

            ! . Get the MM atom charge and scale it for a 1-4 interaction if necessary.
            CGL = ATMCHG(MATOM)

            ! . Define the interatomic distance vector.
            XQM = ( XQ - ATMCRD(1,MATOM) )
            YQM = ( YQ - ATMCRD(2,MATOM) )
            ZQM = ( ZQ - ATMCRD(3,MATOM) )
            IF ( QIMAGE ) THEN
               XQM = XQM - BOXL(1) * ANINT ( XQM / BOXL(1), DP )
               YQM = YQM - BOXL(2) * ANINT ( YQM / BOXL(2), DP )
               ZQM = ZQM - BOXL(3) * ANINT ( ZQM / BOXL(3), DP )
            END IF
            RQM2 = XQM * XQM + YQM * YQM + ZQM * ZQM

            ! . Check the interaction distance.
            IF ( RQM2 > R2OFF ) CYCLE

            ! . Calculate some more distance dependent terms.
            RQM  = SQRT ( RQM2 )
            RQMI = 1.0_DP / RQM
            XN   = RQMI * XQM
            YN   = RQMI * YQM
            ZN   = RQMI * ZQM

            ! . Calculate the switching function.
            IF ( RQM2 > R2ON ) THEN
               SWITCH = ( R2OFF - RQM2 )**2 * ( R2OFF + 2.0_DP * RQM2 - 3.0_DP * R2ON ) / GAMMA
            ELSE
               SWITCH = 1.0_DP
            END IF

            ! . QM atoms with one orbital.
            IF ( NATORB(NQM) == 1 ) THEN

               ! . Calculate some factors.
               RQMB = ANGSTROMS_TO_BOHRS * RQM

               ! . Calculate the SS integral and its derivative.
               CALL SEINTC_INTEGRALB ( RQMB, RHO0, SWITCH, SS, DSS )

               ! . Calculate the integrals.
               E1B(1)    = AU_TO_EV * SS
               E1B(2:10) = 0.0_DP

            ! . QM atoms with more than one orbital.
            ELSE

               ! . Calculate some intermediate quantities.
               RHO1 = ( 0.5_DP / AD(NQM) + RHO0MM )**2
               RHO2 = ( 0.5_DP / AQ(NQM) + RHO0MM )**2
               D1   = DD(NQM)
               D2   = 2.0_DP * QQ(NQM)
               D4   = D2 * D2
               XN2  = XN * XN
               YN2  = YN * YN
               ZN2  = ZN * ZN

               ! . Calculate some factors.
               RQMB = ANGSTROMS_TO_BOHRS * RQM

               ! . Calculate the SS integral and its derivative.
               CALL SEINTC_INTEGRALB ( RQMB, RHO0, SWITCH, SS, DSS )

               ! . Calculate the SPZ integral and its derivative.
               CALL SEINTC_INTEGRALAB ( RQMB,  D1, RHO1, SWITCH, SFCT1P, DSFCT1P )
               CALL SEINTC_INTEGRALAB ( RQMB, -D1, RHO1, SWITCH, SFCT1M, DSFCT1M )
               SPZ  = 0.5_DP * (  SFCT1P -  SFCT1M )
               DSPZ = 0.5_DP * ( DSFCT1P - DSFCT1M )

               ! . Calculate the PZPZ integral and its derivative.
               CALL SEINTC_INTEGRALB  ( RQMB,      RHO2, SWITCH, SFCT20, DSFCT20 )
               CALL SEINTC_INTEGRALAB ( RQMB,  D2, RHO2, SWITCH, SFCT2P, DSFCT2P )
               CALL SEINTC_INTEGRALAB ( RQMB, -D2, RHO2, SWITCH, SFCT2M, DSFCT2M )
               PZPZ  =  SS + 0.25_DP * (  SFCT2P +  SFCT2M ) - 0.5_DP *  SFCT20
               DPZPZ = DSS + 0.25_DP * ( DSFCT2P + DSFCT2M ) - 0.5_DP * DSFCT20

               ! . Calculate the PPPP integral and its derivative.
               CALL SEINTC_INTEGRALB ( RQMB, ( D4 + RHO2 ), SWITCH, SFCT4P, DSFCT4P )
               PPPP  =  SS + 0.5_DP * (  SFCT4P -  SFCT20 )
               DPPPP = DSS + 0.5_DP * ( DSFCT4P - DSFCT20 )

               ! . Fill arrays with the calculated integrals.
               E1B(1)  = SS
               E1B(2)  = XN * SPZ
               E1B(3)  = XN2 * PZPZ + ( YN2 + ZN2 ) * PPPP
               E1B(4)  = YN * SPZ
               E1B(5)  = XN * YN * ( PZPZ - PPPP )
               E1B(6)  = YN2 * PZPZ + ( XN2 + ZN2 ) * PPPP
               E1B(7)  = ZN * SPZ
               E1B(8)  = XN * ZN * ( PZPZ - PPPP )
               E1B(9)  = YN * ZN * ( PZPZ - PPPP )
               E1B(10) = ZN2 * PZPZ + ( XN2 + YN2 ) * PPPP

               ! . Scale the integrals.
               E1B(1)    = AU_TO_EV * E1B(1)
               E1B(2:10) = AU_TO_EV * CGL * E1B(2:10)

            END IF

            ! . Calculate the nuclear/nuclear terms.
            TEMP1 = E1B(1)
            TEMP2 = CGQM * CGL
            TEMP3 = ABS ( TEMP2 ) * EXP ( - ALPHA * RQM )
            TEMP4 = ABS ( TEMP2 ) * EXP ( - ALPMM * RQM )
            ENUC  = TEMP1 * ( TEMP2 + TEMP3 + TEMP4 )

            IF ( QAM1PM3 ) THEN
               DO I = 1,NTERM
                  F1    = FN1(NQM,I)
                  F2    = FN2(NQM,I)
                  F3    = FN3(NQM,I)
                  TEMP1 = EXP ( MAX ( -30.0_DP, - F2 * ( RQM - F3 )**2 ) )
                  ENUC  = ENUC + CGQM * CGL * RQMI * F1 * TEMP1
               END DO
            END IF

            ! . Scale the E1B(1) integral by CGL now ENUC is finished.
            E1B(1) = CGL * E1B(1)

            ! . Get the density elements for the atom.
            JJ = 0
            DO I = N1,N2
               II = (I * (I - 1)) / 2 + N1 - 1
               DO J = N1,I
                  II = II + 1
                  JJ = JJ + 1
                  PFAC(JJ) = DENFAC(JJ) * DENMAT(II,1)
               END DO
            END DO

            ! . Calculate the interaction energy.
            ENB = ENUC - E1B(1) * PFAC(1)
            IF ( NATORB(NQM) > 1 ) ENB = ENB - DOT_PRODUCT ( E1B(2:10), PFAC(2:10) )

            ! . Save the value.
            EQMMMEL(QINDX,MATOM) = ENB

         END DO
      END DO

      END SUBROUTINE QMPARTITION_ELECTROSTATIC

   END SUBROUTINE MOPAC_ENERGY_ANALYSIS

   !-------------------------------------
   SUBROUTINE MOPAC_ENERGY_ANALYSIS_CLEAR
   !-------------------------------------

   ! . Scalars.
   QANALYSIS = .FALSE.

   ! . Arrays.
   IF ( ALLOCATED ( EA1           ) ) DEALLOCATE ( EA1           )
   IF ( ALLOCATED ( EA2           ) ) DEALLOCATE ( EA2           )
   IF ( ALLOCATED ( EEEREPULSION  ) ) DEALLOCATE ( EEEREPULSION  )
   IF ( ALLOCATED ( EENATTRACTION ) ) DEALLOCATE ( EENATTRACTION )
   IF ( ALLOCATED ( EEXCHANGE     ) ) DEALLOCATE ( EEXCHANGE     )
   IF ( ALLOCATED ( ENNREPULSION  ) ) DEALLOCATE ( ENNREPULSION  )
   IF ( ALLOCATED ( ERESONANCE    ) ) DEALLOCATE ( ERESONANCE    )
   IF ( ALLOCATED ( EQMMMEL       ) ) DEALLOCATE ( EQMMMEL       )

   END SUBROUTINE MOPAC_ENERGY_ANALYSIS_CLEAR

   !-----------------------------------------------------------------------------------------------------
   SUBROUTINE MOPAC_ENERGY_ANALYSIS_PRINT ( SELECTION1, SELECTION2, PRINT, EINTERACTION, ESELF1, ESELF2 )
   !-----------------------------------------------------------------------------------------------------

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION1, SELECTION2

   ! . Optional scalar arguments.
   LOGICAL,            INTENT(IN),  OPTIONAL :: PRINT
   REAL ( KIND = DP ), INTENT(OUT), OPTIONAL :: EINTERACTION, ESELF1, ESELF2

   ! . Local scalars.
   LOGICAL            :: QPRINT
   REAL ( KIND = DP ) :: E1EL1, E1EL2, E1LJ1, E1LJ2, E1MM1, E1MM2, E1QM1, E1QM2, E2EL, E2LJ, E2MM, E2QM

   ! . Check to see if an analysis exists (only if there are QM atoms).
   IF ( NATOMSQM > 0 ) THEN
      IF ( .NOT. QANALYSIS ) CALL PRINT_ERROR ( "MOPAC_ENERGY_ANALYSIS_PRINT", "An analysis has not been done." )
   END IF

   ! . Check to see that some atoms have been selected.
   IF ( ( COUNT ( SELECTION1 ) <= 0 ) .OR. ( COUNT ( SELECTION2 ) <= 0 ) ) THEN
      CALL PRINT_ERROR ( "MOPAC_ENERGY_ANALYSIS_PRINT", "There is a null atom selection." )
   END IF

   ! . Check for overlaps between the selections.
   IF ( ANY ( SELECTION1 .AND. SELECTION2 ) ) THEN
      CALL PRINT_ERROR ( "MOPAC_ENERGY_ANALYSIS_PRINT", "An atom appears in both selections." )
   END IF

   ! . Get the contributions from the QM and the QM/MM electrostatic energies.
   CALL ENERGY_NON_BONDING_SELF ( E1EL1, E1LJ1, SELECTION1 )
   CALL ENERGY_NON_BONDING_SELF ( E1EL2, E1LJ2, SELECTION2 )
   CALL ENERGY_NON_BONDING_INTERACTION ( E2EL, E2LJ, SELECTION1, SELECTION2 )

   ! . Sum the total energies.
   E1MM1 = E1EL1 + E1LJ1
   E1MM2 = E1EL2 + E1LJ2
   E2MM  = E2EL  + E2LJ

   ! . Get the contributions from the QM and the QM/MM electrostatic energies.
   CALL GET_QM_VALUES_SELF ( E1QM1, SELECTION1 )
   CALL GET_QM_VALUES_SELF ( E1QM2, SELECTION2 )
   CALL GET_QM_VALUES_INTERACTION

   ! . Set the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Print out the results.
   IF ( QPRINT ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800", PAGEWIDTH = 90, VARIABLEWIDTH = 20 )
      CALL PRINT_SUMMARY_START ( "Self and Interaction Energies (kJ mol^-1)" )
      WRITE ( PRINT_LINE, "(I20)"   ) COUNT ( SELECTION1 ) ; CALL PRINT_SUMMARY_ELEMENT ( "Atoms in Selection 1" )
      WRITE ( PRINT_LINE, "(F20.4)" ) E1MM1 + E1QM1        ; CALL PRINT_SUMMARY_ELEMENT ( "Self Energy 1" )
      WRITE ( PRINT_LINE, "(I20)"   ) COUNT ( SELECTION2 ) ; CALL PRINT_SUMMARY_ELEMENT ( "Atoms in Selection 2" )
      WRITE ( PRINT_LINE, "(F20.4)" ) E1MM2 + E1QM2        ; CALL PRINT_SUMMARY_ELEMENT ( "Self Energy 2" )
      WRITE ( PRINT_LINE, "(F20.4)" ) E2MM  + E2QM         ; CALL PRINT_SUMMARY_ELEMENT ( "Interaction Energy" )
      CALL PRINT_SUMMARY_STOP
  END IF

   ! . Set the output arguments if they are present.
   IF ( PRESENT ( EINTERACTION ) ) EINTERACTION = E2MM  + E2QM
   IF ( PRESENT ( ESELF1       ) ) ESELF1       = E1MM1 + E1QM1
   IF ( PRESENT ( ESELF2       ) ) ESELF2       = E1MM2 + E1QM2

   !===========================================================================
   CONTAINS
   !===========================================================================

      !-----------------------------------
      SUBROUTINE GET_QM_VALUES_INTERACTION
      !-----------------------------------

      ! . Local scalars.
      INTEGER :: IATOM, IJ, JATOM, QINDXI, QINDXJ

      ! . Initialization.
      E2QM = 0.0_DP

      ! . Loop over the atoms of the first selection.
      DO IATOM = 1,NATOMS
         IF ( SELECTION1(IATOM) ) THEN

            ! . Get the QM index for the atom.
            QINDXI = ATMQMI(IATOM)

            ! . Loop over the atoms of the second selection.
            DO JATOM = 1,NATOMS
               IF ( SELECTION2(JATOM) ) THEN

                  ! . Get the QM index for the atom.
                  QINDXJ = ATMQMI(JATOM)

                  ! . A QM/QM interaction.
                  IF ( ( QINDXI > 0 ) .AND. ( QINDXJ > 0 ) ) THEN

                     ! . Get the interaction address.
                     IF ( QINDXI > QINDXJ ) THEN
                        IJ = ( QINDXI * ( QINDXI - 1 ) ) / 2 + QINDXJ
                     ELSE
                        IJ = ( QINDXJ * ( QINDXJ - 1 ) ) / 2 + QINDXI
                     END IF
                     E2QM = E2QM + EEEREPULSION(IJ) + EENATTRACTION(IJ) + EEXCHANGE(IJ) + ENNREPULSION(IJ) + ERESONANCE(IJ)

                  ! . A QM/MM interaction.
                  ELSE IF ( QINDXI > 0 ) THEN
                     E2QM = E2QM + EQMMMEL(QINDXI,JATOM)

                  ! . A MM/QM interaction.
                  ELSE IF ( QINDXJ > 0 ) THEN
                     E2QM = E2QM + EQMMMEL(QINDXJ,IATOM)

                  ! . There are no MM interactions here.
                  ENDIF
               END IF
            END DO
         END IF
      END DO

      ! . Convert E2QM to kJ mol^-1.
      E2QM = EV_TO_KJ * E2QM

      END SUBROUTINE GET_QM_VALUES_INTERACTION

      !-------------------------------------------------
      SUBROUTINE GET_QM_VALUES_SELF ( ESELF, SELECTION )
      !-------------------------------------------------

      ! . Scalar arguments.
      REAL ( KIND = DP ), INTENT(OUT) :: ESELF

      ! . Array arguments.
      LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION

      ! . Local scalars.
      INTEGER            :: IATOM, IJ, JATOM, NI, QINDXI, QINDXJ
      REAL ( KIND = DP ) :: EATHEAT

      ! . Initialization.
      EATHEAT = 0.0_DP
      ESELF   = 0.0_DP

      ! . Outer loop over atoms.
      DO IATOM = 1,NATOMS

         ! . The atom is in the selection.
         IF ( SELECTION(IATOM) ) THEN

            ! . The atom is a QM atom.
            IF ( ATMQMI(IATOM) > 0 ) THEN

               ! . Get the QM atom index and its atomic number.
               QINDXI = ATMQMI(IATOM)
               NI     = QMATOM(QINDXI)%NUMBER

               ! . Add in the contributions from the one-centre terms.
               EATHEAT = EATHEAT + KCAL_TO_KJ * EHEAT(NI) - EV_TO_KJ * EISOL(NI)
               ESELF   = ESELF   + EA1(QINDXI) + EA2(QINDXI)

               ! . Inner loop over the atoms.
               DO JATOM = 1,NATOMS

                  ! . The atom is in the selection.
                  IF ( SELECTION(JATOM) ) THEN

                      ! . The atom is QM.
                      IF ( ATMQMI(JATOM) > 0 ) THEN

                         ! . Get the QM atom index.
                         QINDXJ = ATMQMI(JATOM)

                         ! . Reject all QM interactions if QINDXJ >= QINDXI.
                         IF ( QINDXJ >= QINDXI ) CYCLE

                         ! . Add in the contributions from the two-centre terms.
                         IJ    = ( QINDXI * ( QINDXI - 1 ) ) / 2 + QINDXJ
                         ESELF = ESELF + EEEREPULSION(IJ) + EENATTRACTION(IJ) + EEXCHANGE(IJ) + ENNREPULSION(IJ) + ERESONANCE(IJ)

                      ! . The atom is MM.
                      ELSE

                         ! . Add in the contribution from the QM/MM interaction.
                         ESELF = ESELF + EQMMMEL(QINDXI,JATOM)

                      END IF
                  END IF
               END DO
            ENDIF
         END IF
      END DO

      ! . Convert ESELF to kJ mol^-1.
      ESELF = EV_TO_KJ * ESELF + EATHEAT

      END SUBROUTINE GET_QM_VALUES_SELF

   END SUBROUTINE MOPAC_ENERGY_ANALYSIS_PRINT

   !--------------------------------------------------
   SUBROUTINE MOPAC_IONIZATION_POTENTIAL ( IP, PRINT )
   !--------------------------------------------------

   ! . Scalar arguments.
   LOGICAL,            INTENT(IN),  OPTIONAL :: PRINT
   REAL ( KIND = DP ), INTENT(OUT), OPTIONAL :: IP

   ! . Other local scalars.
   INTEGER            :: IBEAD
   LOGICAL            :: QPRINT, QUHF
   REAL ( KIND = DP ) :: IPOT

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: FOCK1, FOCK2, ORBENE1, ORBENE2
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: MORBS1, MORBS2

   ! . Check to see that everything is there, etc.
   IF ( ( NBASIS <= 0 ) .OR. ( ( NALPHA <= 0 ) .AND. ( NBETA <= 0 ) ) .OR. ( .NOT. ALLOCATED ( DENMAT ) ) ) RETURN

   ! . Set the UHF flag.
   QUHF = ALLOCATED ( DENMATA ) .AND. ALLOCATED ( DENMATB )

   ! . Allocate some space.
   ALLOCATE ( FOCK1(1:NBASTR), MORBS1(1:NBASIS,1:NBASIS), ORBENE1(1:NBASIS) )
   IF ( QUHF ) ALLOCATE ( FOCK2(1:NBASTR), MORBS2(1:NBASIS,1:NBASIS), ORBENE2(1:NBASIS) )

   ! . Initialization.
   IPOT = 0.0_DP

   ! . Loop over the number of PI beads.
   DO IBEAD = 1,NPIBEADS

      ! . Set up the one-electron matrix.
      CALL INTEGRALS_HCORE  ( IBEAD )

      ! . Build the two-electron parts of the Fock matrices.
      IF ( QUHF ) THEN
         CALL FOCK_MATRIX_UNRESTRICTED ( IBEAD, DENMAT(:,IBEAD), DENMATA(:,IBEAD), DENMATB(:,IBEAD), FOCK1, FOCK2 )
      ELSE
         CALL FOCK_MATRIX_RESTRICTED ( IBEAD, DENMAT(:,IBEAD), FOCK1 )
      END IF

      ! . Produce the full FOCK matrices by adding in the one-electron matrix.
      FOCK1 = FOCK1 + HCORE
      IF ( QUHF ) FOCK2 = FOCK2 + HCORE

      ! . Diagonalize the Fock matrices.
      CALL SYMMETRIC_UPPER ( FOCK1, ORBENE1, MORBS1 )
      IF ( QUHF ) CALL SYMMETRIC_UPPER ( FOCK2, ORBENE2, MORBS2 )

      ! . Get IPOT.
      IPOT = IPOT + ORBENE1(NALPHA)

   END DO

   ! . Deallocate space.
   DEALLOCATE ( FOCK1, MORBS1, ORBENE1 )
   IF ( QUHF ) DEALLOCATE ( FOCK2, MORBS2, ORBENE2 )

   ! . Scale IPOT.
   IPOT = ABS ( IPOT ) / REAL ( NPIBEADS, DP )

   ! . Get the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Print the result if necessary.
   IF ( QPRINT ) THEN
      WRITE ( PRINT_LINE, "(A,F7.2,A)" ) "Ionization Potential = ", IPOT, " eV."
      CALL PRINT_PARAGRAPH
   END IF

   ! . Save IPOT if necessary.
   IF ( PRESENT ( IP ) ) IP = IPOT

   END SUBROUTINE MOPAC_IONIZATION_POTENTIAL

END MODULE MOPAC_ANALYSIS
