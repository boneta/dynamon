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
!                           The MOPAC Parameter Module
!===============================================================================
!
! . Integer arrays:
!
!     NATORB             the number of orbitals for each element.
!
! . Logical arrays:
!
!     SEPAR              flags for each element to show whether the
!                        parameterisation exists.
!
! . Real arrays:
!
!     AD, AM, AQ,        parameters for the multipole approximation to
!     DD, QQ             the two-centre, two-electron integrals.
!     ALP                the exponents used to calculate the core-core
!                        repulsion terms.
!     BETAD, BETAP,      parameters for the two-centre, one-electron core
!     BETAS              resonance integrals.
!     CORE               the CORE charge on each element.
!     EHEAT              the experimental heat of formation for each element.
!     EISOL              the calculated electronic energies for each element.
!     FN1, FN2, FN3      AM1/PM3 specific parameters for the core-core repulsion
!                        interaction.
!     GDD, GPD, GSD,     parameters for the Coulomb and exchange one-centre,
!     GP2, GPP, GSP,     two-electron integrals.
!     GSS, HSP
!     USS, UPP, UDD      electron kinetic energy integrals.
!     ZD, ZP, ZS         the Slater exponents of the basis functions.
!
! . Subroutines:
!
!     MOPAC_PARAMETERS_INITIALIZE     initialize the MOPAC parameter arrays.
!
! . Notes:
!
!   Enthalpies of formation of gaseous atoms from the CRC Handbook 1981-1982
!   and Annual Reports 71B, p119 (1974).
!
!===============================================================================
MODULE MOPAC_PARAMETERS

! . Module declarations.
USE CONSTANTS,   ONLY : AU_TO_EV
USE DEFINITIONS, ONLY : DP
USE ELEMENTS,    ONLY : NELEMENTS

IMPLICIT NONE
PUBLIC
#ifndef PGPC
SAVE
#endif

! . Options.
LOGICAL :: QNEWAM1PARAMETERS = .FALSE.

! . Integer arrays.
INTEGER, DIMENSION(1:NELEMENTS) :: NATORB = 0

! . Logical arrays.
LOGICAL, DIMENSION(1:NELEMENTS) :: SEPAR = .FALSE.

! . Real arrays.
REAL ( KIND = DP ), DIMENSION(1:NELEMENTS) :: AD    = 0.0_DP, ALP   = 0.0_DP, AM    = 0.0_DP, AQ   = 0.0_DP, &
                                              BETAD = 0.0_DP, BETAP = 0.0_DP, BETAS = 0.0_DP, CORE = 0.0_DP, &
                                              DD    = 0.0_DP, EHEAT = 0.0_DP, EISOL = 0.0_DP, GDD  = 0.0_DP, &
                                              GPD   = 0.0_DP, GSD   = 0.0_DP, GP2   = 0.0_DP, GPP  = 0.0_DP, &
                                              GSP   = 0.0_DP, GSS   = 0.0_DP, HSP   = 0.0_DP, QQ   = 0.0_DP, &
                                              USS   = 0.0_DP, UPP   = 0.0_DP, UDD   = 0.0_DP, ZD   = 0.0_DP, &
                                              ZP    = 0.0_DP, ZS    = 0.0_DP

! . AM1 and PM3 specific arrays.
REAL ( KIND = DP ), DIMENSION(1:NELEMENTS,1:4) :: FN1 = 0.0_DP, FN2 = 0.0_DP, FN3 = 0.0_DP

!   PDDG specific arrays
REAL ( KIND = DP ), DIMENSION(1:NELEMENTS,1:2) :: PDDGC = 0.0_DP, PDDGE = 0.0_DP

!===============================================================================
CONTAINS
!===============================================================================

   !------------------------------------------------
   SUBROUTINE MOPAC_PARAMETERS_INITIALIZE ( METHOD )
   !------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: METHOD

   ! . Initialize the common data elements.
   CALL PARAMETERS_COMMON

   ! . Branch on the method.
   IF ( METHOD == "AM1" ) THEN
      IF ( QNEWAM1PARAMETERS ) THEN
         CALL PARAMETERS_AM1TEST
      ELSE
         CALL PARAMETERS_AM1
      END IF
   ELSE IF ( METHOD == "MNDO" ) THEN
      CALL PARAMETERS_MNDO
   ELSE IF ( METHOD == "PDDG" ) THEN
      CALL PARAMETERS_PDDG
   ELSE IF ( METHOD == "PM3" ) THEN
      CALL PARAMETERS_PM3
   ELSE IF ( METHOD == "RM1" ) THEN
      CALL PARAMETERS_RM1
   END IF

   !============================================================================
   CONTAINS
   !============================================================================

      !---------------------------
      SUBROUTINE PARAMETERS_COMMON
      !---------------------------

      ! . Initialize the data arrays.
      SEPAR = .FALSE.

      AD    = 0.0_DP ; ALP   = 0.0_DP ; AM    = 0.0_DP ; AQ   = 0.0_DP
      BETAD = 0.0_DP ; BETAP = 0.0_DP ; BETAS = 0.0_DP ; CORE = 0.0_DP
      DD    = 0.0_DP ; EHEAT = 0.0_DP ; EISOL = 0.0_DP ; GDD  = 0.0_DP
      GPD   = 0.0_DP ; GSD   = 0.0_DP ; GP2   = 0.0_DP ; GPP  = 0.0_DP
      GSP   = 0.0_DP ; GSS   = 0.0_DP ; HSP   = 0.0_DP ; QQ   = 0.0_DP
      USS   = 0.0_DP ; UPP   = 0.0_DP ; UDD   = 0.0_DP ; ZD   = 0.0_DP
      ZP    = 0.0_DP ; ZS    = 0.0_DP

      FN1   = 0.0_DP ; FN2   = 0.0_DP ; FN3   = 0.0_DP

      ! . The number of orbitals per atom, the core charges and the enthalpies of formation.

      NATORB (  1) = 1 ; CORE (  1) =  1.0_DP ; EHEAT (  1) =  52.102_DP
      NATORB (  2) = 1 ; CORE (  2) =  0.0_DP ; EHEAT (  2) =   0.000_DP
      NATORB (  3) = 4 ; CORE (  3) =  1.0_DP ; EHEAT (  3) =  38.410_DP
      NATORB (  4) = 4 ; CORE (  4) =  2.0_DP ; EHEAT (  4) =  76.960_DP
      NATORB (  5) = 4 ; CORE (  5) =  3.0_DP ; EHEAT (  5) = 135.700_DP
      NATORB (  6) = 4 ; CORE (  6) =  4.0_DP ; EHEAT (  6) = 170.890_DP
      NATORB (  7) = 4 ; CORE (  7) =  5.0_DP ; EHEAT (  7) = 113.000_DP
      NATORB (  8) = 4 ; CORE (  8) =  6.0_DP ; EHEAT (  8) =  59.559_DP
      NATORB (  9) = 4 ; CORE (  9) =  7.0_DP ; EHEAT (  9) =  18.890_DP
      NATORB ( 10) = 4 ; CORE ( 10) =  0.0_DP ; EHEAT ( 10) =   0.000_DP

      NATORB ( 11) = 0 ; CORE ( 11) =  1.0_DP ; EHEAT ( 11) =  25.850_DP
      NATORB ( 12) = 4 ; CORE ( 12) =  2.0_DP ; EHEAT ( 12) =  35.000_DP
      NATORB ( 13) = 4 ; CORE ( 13) =  3.0_DP ; EHEAT ( 13) =  79.490_DP
      NATORB ( 14) = 4 ; CORE ( 14) =  4.0_DP ; EHEAT ( 14) = 108.390_DP
      NATORB ( 15) = 4 ; CORE ( 15) =  5.0_DP ; EHEAT ( 15) =  75.570_DP
      NATORB ( 16) = 4 ; CORE ( 16) =  6.0_DP ; EHEAT ( 16) =  66.400_DP
      NATORB ( 17) = 4 ; CORE ( 17) =  7.0_DP ; EHEAT ( 17) =  28.990_DP
      NATORB ( 18) = 4 ; CORE ( 18) =  0.0_DP ; EHEAT ( 18) =   0.000_DP
      NATORB ( 19) = 0 ; CORE ( 19) =  1.0_DP ; EHEAT ( 19) =  21.420_DP
      NATORB ( 20) = 4 ; CORE ( 20) =  2.0_DP ; EHEAT ( 20) =  42.600_DP

      NATORB ( 21) = 9 ; CORE ( 21) =  3.0_DP ; EHEAT ( 21) =  90.300_DP
      NATORB ( 22) = 9 ; CORE ( 22) =  4.0_DP ; EHEAT ( 22) = 112.300_DP
      NATORB ( 23) = 9 ; CORE ( 23) =  5.0_DP ; EHEAT ( 23) = 122.900_DP
      NATORB ( 24) = 9 ; CORE ( 24) =  6.0_DP ; EHEAT ( 24) =  95.000_DP
      NATORB ( 25) = 9 ; CORE ( 25) =  7.0_DP ; EHEAT ( 25) =  67.700_DP
      NATORB ( 26) = 9 ; CORE ( 26) =  8.0_DP ; EHEAT ( 26) =  99.300_DP
      NATORB ( 27) = 9 ; CORE ( 27) =  9.0_DP ; EHEAT ( 27) = 102.400_DP
      NATORB ( 28) = 9 ; CORE ( 28) = 10.0_DP ; EHEAT ( 28) = 102.800_DP
      NATORB ( 29) = 9 ; CORE ( 29) = 11.0_DP ; EHEAT ( 29) =  80.700_DP
      NATORB ( 30) = 4 ; CORE ( 30) =  2.0_DP ; EHEAT ( 30) =  31.170_DP

      NATORB ( 31) = 4 ; CORE ( 31) =  3.0_DP ; EHEAT ( 31) =  65.400_DP
      NATORB ( 32) = 4 ; CORE ( 32) =  4.0_DP ; EHEAT ( 32) =  89.500_DP
      NATORB ( 33) = 4 ; CORE ( 33) =  5.0_DP ; EHEAT ( 33) =  72.300_DP
      NATORB ( 34) = 4 ; CORE ( 34) =  6.0_DP ; EHEAT ( 34) =  54.300_DP
      NATORB ( 35) = 4 ; CORE ( 35) =  7.0_DP ; EHEAT ( 35) =  26.740_DP
      NATORB ( 36) = 4 ; CORE ( 36) =  0.0_DP ; EHEAT ( 36) =   0.000_DP
      NATORB ( 37) = 4 ; CORE ( 37) =  1.0_DP ; EHEAT ( 37) =  19.600_DP
      NATORB ( 38) = 4 ; CORE ( 38) =  2.0_DP ; EHEAT ( 38) =  39.100_DP
      NATORB ( 39) = 9 ; CORE ( 39) =  3.0_DP ; EHEAT ( 39) = 101.500_DP
      NATORB ( 40) = 9 ; CORE ( 40) =  4.0_DP ; EHEAT ( 40) = 145.500_DP

      NATORB ( 41) = 9 ; CORE ( 41) =  5.0_DP ; EHEAT ( 41) = 172.400_DP
      NATORB ( 42) = 9 ; CORE ( 42) =  6.0_DP ; EHEAT ( 42) = 157.300_DP 
      NATORB ( 43) = 9 ; CORE ( 43) =  7.0_DP ; EHEAT ( 43) =   0.000_DP
      NATORB ( 44) = 9 ; CORE ( 44) =  8.0_DP ; EHEAT ( 44) = 155.500_DP
      NATORB ( 45) = 9 ; CORE ( 45) =  9.0_DP ; EHEAT ( 45) = 133.000_DP
      NATORB ( 46) = 9 ; CORE ( 46) = 10.0_DP ; EHEAT ( 46) =  90.000_DP
      NATORB ( 47) = 9 ; CORE ( 47) = 11.0_DP ; EHEAT ( 47) =  68.100_DP
      NATORB ( 48) = 4 ; CORE ( 48) =  2.0_DP ; EHEAT ( 48) =  26.720_DP
      NATORB ( 49) = 4 ; CORE ( 49) =  3.0_DP ; EHEAT ( 49) =  58.000_DP
      NATORB ( 50) = 4 ; CORE ( 50) =  4.0_DP ; EHEAT ( 50) =  72.200_DP

      NATORB ( 51) = 4 ; CORE ( 51) =  5.0_DP ; EHEAT ( 51) =  63.200_DP
      NATORB ( 52) = 4 ; CORE ( 52) =  6.0_DP ; EHEAT ( 52) =  47.000_DP 
      NATORB ( 53) = 4 ; CORE ( 53) =  7.0_DP ; EHEAT ( 53) =  25.517_DP
      NATORB ( 54) = 4 ; CORE ( 54) =  0.0_DP ; EHEAT ( 54) =   0.000_DP
      NATORB ( 55) = 2 ; CORE ( 55) =  1.0_DP ; EHEAT ( 55) =  18.700_DP
      NATORB ( 56) = 2 ; CORE ( 56) =  2.0_DP ; EHEAT ( 56) =  42.500_DP
      NATORB ( 57) = 8 ; CORE ( 57) =  3.0_DP ; EHEAT ( 57) =   0.000_DP
      NATORB ( 58) = 8 ; CORE ( 58) =  4.0_DP ; EHEAT ( 58) = 101.300_DP
      NATORB ( 59) = 8 ; CORE ( 59) =  5.0_DP ; EHEAT ( 59) =   0.000_DP
      NATORB ( 60) = 8 ; CORE ( 60) =  6.0_DP ; EHEAT ( 60) =   0.000_DP

      NATORB ( 61) = 8 ; CORE ( 61) =  7.0_DP ; EHEAT ( 61) =   0.000_DP
      NATORB ( 62) = 8 ; CORE ( 62) =  8.0_DP ; EHEAT ( 62) =  49.400_DP
      NATORB ( 63) = 8 ; CORE ( 63) =  9.0_DP ; EHEAT ( 63) =   0.000_DP
      NATORB ( 64) = 8 ; CORE ( 64) = 10.0_DP ; EHEAT ( 64) =   0.000_DP
      NATORB ( 65) = 8 ; CORE ( 65) = 11.0_DP ; EHEAT ( 65) =   0.000_DP
      NATORB ( 66) = 8 ; CORE ( 66) = 12.0_DP ; EHEAT ( 66) =   0.000_DP
      NATORB ( 67) = 8 ; CORE ( 67) = 13.0_DP ; EHEAT ( 67) =   0.000_DP
      NATORB ( 68) = 8 ; CORE ( 68) = 14.0_DP ; EHEAT ( 68) =  75.800_DP
      NATORB ( 69) = 8 ; CORE ( 69) = 15.0_DP ; EHEAT ( 69) =   0.000_DP
      NATORB ( 70) = 8 ; CORE ( 70) = 16.0_DP ; EHEAT ( 70) =  36.350_DP

      NATORB ( 71) = 9 ; CORE ( 71) =  3.0_DP ; EHEAT ( 71) =   0.000_DP
      NATORB ( 72) = 9 ; CORE ( 72) =  4.0_DP ; EHEAT ( 72) = 148.000_DP
      NATORB ( 73) = 9 ; CORE ( 73) =  5.0_DP ; EHEAT ( 73) = 186.900_DP
      NATORB ( 74) = 9 ; CORE ( 74) =  6.0_DP ; EHEAT ( 74) = 203.100_DP
      NATORB ( 75) = 9 ; CORE ( 75) =  7.0_DP ; EHEAT ( 75) = 185.000_DP
      NATORB ( 76) = 9 ; CORE ( 76) =  8.0_DP ; EHEAT ( 76) = 188.000_DP
      NATORB ( 77) = 9 ; CORE ( 77) =  9.0_DP ; EHEAT ( 77) = 160.000_DP
      NATORB ( 78) = 9 ; CORE ( 78) = 10.0_DP ; EHEAT ( 78) = 135.200_DP
      NATORB ( 79) = 9 ; CORE ( 79) = 11.0_DP ; EHEAT ( 79) =  88.000_DP
      NATORB ( 80) = 4 ; CORE ( 80) =  2.0_DP ; EHEAT ( 80) =  14.690_DP

      NATORB ( 81) = 4 ; CORE ( 81) =  3.0_DP ; EHEAT ( 81) =  43.550_DP
      NATORB ( 82) = 4 ; CORE ( 82) =  4.0_DP ; EHEAT ( 82) =  46.620_DP
      NATORB ( 83) = 4 ; CORE ( 83) =  5.0_DP ; EHEAT ( 83) =  50.100_DP
      NATORB ( 84) = 4 ; CORE ( 84) =  6.0_DP ; EHEAT ( 84) =   0.000_DP
      NATORB ( 85) = 4 ; CORE ( 85) =  7.0_DP ; EHEAT ( 85) =   0.000_DP

      END SUBROUTINE PARAMETERS_COMMON

      !------------------------
      SUBROUTINE PARAMETERS_AM1
      !------------------------

      ! . Hydrogen.
      SEPAR(1) = .TRUE.
      ALP(1)   =   2.8823240_DP
      EISOL(1) = -11.3964270_DP
      BETAS(1) =  -6.1737870_DP
      ZS(1)    =   1.1880780_DP
      AM(1)    =   0.4721793_DP
      AD(1)    =   0.4721793_DP
      AQ(1)    =   0.4721793_DP
      USS(1)   = -11.3964270_DP
      GSS(1)   =  12.8480000_DP

      FN1(1,1) =  0.1227960_DP ; FN2(1,1) = 5.0000000_DP ; FN3(1,1) = 1.2000000_DP
      FN1(1,2) =  0.0050900_DP ; FN2(1,2) = 5.0000000_DP ; FN3(1,2) = 1.8000000_DP
      FN1(1,3) = -0.0183360_DP ; FN2(1,3) = 2.0000000_DP ; FN3(1,3) = 2.1000000_DP

      ! . Lithium.
      SEPAR(3) = .TRUE.
      ALP(3)   =   1.2501400_DP
      EISOL(3) =  -5.1280000_DP
      BETAS(3) =  -1.3500400_DP
      BETAP(3) =  -1.3500400_DP
      ZS(3)    =   0.7023800_DP
      ZP(3)    =   0.7023800_DP
      DD(3)    =   2.0549783_DP
      QQ(3)    =   1.7437069_DP
      AM(3)    =   0.2682837_DP
      AD(3)    =   0.2269793_DP
      AQ(3)    =   0.2614581_DP
      USS(3)   =  -5.1280000_DP
      UPP(3)   =  -2.7212000_DP
      GSS(3)   =   7.3000000_DP
      GSP(3)   =   5.4200000_DP
      GPP(3)   =   5.0000000_DP
      GP2(3)   =   4.5200000_DP
      HSP(3)   =   0.8300000_DP

      ! . Beryllium.
      SEPAR(4) = .TRUE.
      ALP(4)   =   1.6694340_DP
      EISOL(4) = -24.2047560_DP
      BETAS(4) =  -4.0170960_DP
      BETAP(4) =  -4.0170960_DP
      ZS(4)    =   1.0042100_DP
      ZP(4)    =   1.0042100_DP
      DD(4)    =   1.4373245_DP
      QQ(4)    =   1.2196103_DP
      AM(4)    =   0.3307607_DP
      AD(4)    =   0.3356142_DP
      AQ(4)    =   0.3846373_DP
      USS(4)   = -16.6023780_DP
      UPP(4)   = -10.7037710_DP
      GSS(4)   =   9.0000000_DP
      GSP(4)   =   7.4300000_DP
      GPP(4)   =   6.9700000_DP
      GP2(4)   =   6.2200000_DP
      HSP(4)   =   1.2800000_DP

      ! . Boron.
      SEPAR(5) = .TRUE.
      ALP(5)   =   2.4469090_DP
      EISOL(5) = -63.7172650_DP
      BETAS(5) =  -9.5991140_DP
      BETAP(5) =  -6.2737570_DP
      ZS(5)    =   1.6117090_DP
      ZP(5)    =   1.5553850_DP
      DD(5)    =   0.9107622_DP
      QQ(5)    =   0.7874223_DP
      AM(5)    =   0.3891951_DP
      AD(5)    =   0.5045152_DP
      AQ(5)    =   0.5678856_DP
      USS(5)   = -34.4928700_DP
      UPP(5)   = -22.6315250_DP
      GSS(5)   =  10.5900000_DP
      GSP(5)   =   9.5600000_DP
      GPP(5)   =   8.8600000_DP
      GP2(5)   =   7.8600000_DP
      HSP(5)   =   1.8100000_DP

      ! . Carbon.
      SEPAR(6) = .TRUE.
      ALP(6)   =    2.6482740_DP
      EISOL(6) = -120.8157940_DP
      BETAS(6) =  -15.7157830_DP
      BETAP(6) =   -7.7192830_DP
      ZS(6)    =    1.8086650_DP
      ZP(6)    =    1.6851160_DP
      DD(6)    =    0.8236736_DP
      QQ(6)    =    0.7268015_DP
      AM(6)    =    0.4494671_DP
      AD(6)    =    0.6082946_DP
      AQ(6)    =    0.6423492_DP
      USS(6)   =  -52.0286580_DP
      UPP(6)   =  -39.6142390_DP
      GSS(6)   =   12.2300000_DP
      GSP(6)   =   11.4700000_DP
      GPP(6)   =   11.0800000_DP
      GP2(6)   =    9.8400000_DP
      HSP(6)   =    2.4300000_DP

      FN1(6,1) =  0.0113550_DP ; FN2(6,1) =  5.0000000_DP ; FN3(6,1) =  1.6000000_DP
      FN1(6,2) =  0.0459240_DP ; FN2(6,2) =  5.0000000_DP ; FN3(6,2) =  1.8500000_DP
      FN1(6,3) = -0.0200610_DP ; FN2(6,3) =  5.0000000_DP ; FN3(6,3) =  2.0500000_DP
      FN1(6,4) = -0.0012600_DP ; FN2(6,4) =  5.0000000_DP ; FN3(6,4) =  2.6500000_DP

      ! . Nitrogen.
      SEPAR(7) = .TRUE.
      ALP(7)   =    2.9472860_DP
      EISOL(7) = -202.4077430_DP
      BETAS(7) =  -20.2991100_DP
      BETAP(7) =  -18.2386660_DP
      ZS(7)    =    2.3154100_DP
      ZP(7)    =    2.1579400_DP
      DD(7)    =    0.6433247_DP
      QQ(7)    =    0.5675528_DP
      AM(7)    =    0.4994487_DP
      AD(7)    =    0.7820840_DP
      AQ(7)    =    0.7883498_DP
      USS(7)   =  -71.8600000_DP
      UPP(7)   =  -57.1675810_DP
      GSS(7)   =   13.5900000_DP
      GSP(7)   =   12.6600000_DP
      GPP(7)   =   12.9800000_DP
      GP2(7)   =   11.5900000_DP
      HSP(7)   =    3.1400000_DP

      FN1(7,1) =  0.0252510_DP ; FN2(7,1) =  5.0000000_DP ; FN3(7,1) =  1.5000000_DP
      FN1(7,2) =  0.0289530_DP ; FN2(7,2) =  5.0000000_DP ; FN3(7,2) =  2.1000000_DP
      FN1(7,3) = -0.0058060_DP ; FN2(7,3) =  2.0000000_DP ; FN3(7,3) =  2.4000000_DP

      ! . Oxygen.
      SEPAR(8) = .TRUE.
      ALP(8)   =    4.4553710_DP
      EISOL(8) = -316.0995200_DP
      BETAS(8) =  -29.2727730_DP
      BETAP(8) =  -29.2727730_DP
      ZS(8)    =    3.1080320_DP
      ZP(8)    =    2.5240390_DP
      DD(8)    =    0.4988896_DP
      QQ(8)    =    0.4852322_DP
      AM(8)    =    0.5667034_DP
      AD(8)    =    0.9961066_DP
      AQ(8)    =    0.9065223_DP
      USS(8)   =  -97.8300000_DP
      UPP(8)   =  -78.2623800_DP
      GSS(8)   =   15.4200000_DP
      GSP(8)   =   14.4800000_DP
      GPP(8)   =   14.5200000_DP
      GP2(8)   =   12.9800000_DP
      HSP(8)   =    3.9400000_DP

      FN1(8,1) = 0.2809620_DP ; FN2(8,1) = 5.0000000_DP ; FN3(8,1) = 0.8479180_DP
      FN1(8,2) = 0.0814300_DP ; FN2(8,2) = 7.0000000_DP ; FN3(8,2) = 1.4450710_DP

      ! . Fluorine.
      SEPAR(9) = .TRUE.
      ALP(9)   =    5.5178000_DP
      EISOL(9) = -482.2905830_DP
      BETAS(9) =  -69.5902770_DP
      BETAP(9) =  -27.9223600_DP
      ZS(9)    =    3.7700820_DP
      ZP(9)    =    2.4946700_DP
      DD(9)    =    0.4145203_DP
      QQ(9)    =    0.4909446_DP
      AM(9)    =    0.6218302_DP
      AD(9)    =    1.2088792_DP
      AQ(9)    =    0.9449355_DP
      USS(9)   = -136.1055790_DP
      UPP(9)   = -104.8898850_DP
      GSS(9)   =   16.9200000_DP
      GSP(9)   =   17.2500000_DP
      GPP(9)   =   16.7100000_DP
      GP2(9)   =   14.9100000_DP
      HSP(9)   =    4.8300000_DP

      FN1(9,1) = 0.2420790_DP ; FN2(9,1) = 4.8000000_DP ; FN3(9,1) = 0.9300000_DP
      FN1(9,2) = 0.0036070_DP ; FN2(9,2) = 4.6000000_DP ; FN3(9,2) = 1.6600000_DP

      ! . Sodium.
      SEPAR(11) = .TRUE.
      ALP(11)   =  1.32_DP
      EISOL(11) =  0.00_DP
      AM(11)    =  0.50_DP

      ! . Magnesium (Hutter et al., JPCB 102, 8080-8090, 1998).
      SEPAR(12) = .TRUE.
      USS  ( 12) =  -14.96959313_DP
      UPP  ( 12) =  -11.56229248_DP
      BETAS( 12) =   -1.25974355_DP
      BETAP( 12) =   -0.77836604_DP
      ZS   ( 12) =    1.22339270_DP
      ZP   ( 12) =    1.02030798_DP
      ALP  ( 12) =    1.67049799_DP
      EISOL( 12) =  -22.43786349_DP ! From CALPAR formula. Otherwise PM3 value?
      GSS  ( 12) =    7.50132277_DP
      GSP  ( 12) =    6.34591536_DP
      GPP  ( 12) =    4.77534467_DP
      GP2  ( 12) =    4.34017279_DP
      HSP  ( 12) =    0.48930466_DP
      DD   ( 12) =    1.75012115_DP ! From formula (15) page 94 Thiel + Dewar.
      QQ   ( 12) =    1.64001467_DP ! From formula (16) page 95 Thiel + Dewar.
      AM   ( 12) =    0.27568257_DP ! From formula (21) page 95 Thiel + Dewar (rho0 = 1/2 Gss = 1/2 AM so AM = Gss/AU_TO_EV).
      AD   ( 12) =    0.19957236_DP ! By solving numerically equation (19) page 95 Thiel + Dewar.
      AQ   ( 12) =    0.26440152_DP ! By solving numerically equation (20) page 95 Thiel + Dewar.

      FN1( 12,1) =  2.55017735_DP ; FN2( 12,1) = 4.29397225_DP ; FN3( 12,1) = 0.79989601_DP
      FN1( 12,2) = -0.00565806_DP ; FN2( 12,2) = 2.96053910_DP ; FN3( 12,2) = 1.47499983_DP
      FN1( 12,3) = -0.00610286_DP ; FN2( 12,3) = 2.61416919_DP ; FN3( 12,3) = 2.42604040_DP

      ! . Aluminium.
      SEPAR(13) = .TRUE.
      ALP(13)   =   1.9765860_DP
      EISOL(13) = -46.4208150_DP
      BETAS(13) =  -3.8668220_DP
      BETAP(13) =  -2.3171460_DP
      ZS(13)    =   1.5165930_DP
      ZP(13)    =   1.3063470_DP
      ZD(13)    =   1.0000000_DP
      DD(13)    =   1.4040443_DP
      QQ(13)    =   1.2809154_DP
      AM(13)    =   0.2973172_DP
      AD(13)    =   0.2630229_DP
      AQ(13)    =   0.3427832_DP
      USS(13)   = -24.3535850_DP
      UPP(13)   = -18.3636450_DP
      GSS(13)   =   8.0900000_DP
      GSP(13)   =   6.6300000_DP
      GPP(13)   =   5.9800000_DP
      GP2(13)   =   5.4000000_DP
      HSP(13)   =   0.7000000_DP

      FN1(13,1) = 0.0900000_DP ; FN2(13,1) = 12.3924430_DP ; FN3(13,1) = 2.0503940_DP

      ! . Silicon.
      SEPAR(14) = .TRUE.
      ALP(14)   =   2.257816_DP
      EISOL(14) = -79.0017420_DP
      BETAS(14) =  -3.784852_DP
      BETAP(14) =  -1.968123_DP
      ZS(14)    =   1.830697_DP
      ZP(14)    =   1.2849530_DP
      ZD(14)    =   1.0000000_DP
      DD(14)    =   1.1631107_DP
      QQ(14)    =   1.3022422_DP
      AM(14)    =   0.3608967_DP
      AD(14)    =   0.3829813_DP
      AQ(14)    =   0.3999442_DP
      USS(14)   = -33.9536220_DP
      UPP(14)   = -28.9347490_DP
      GSS(14)   =   9.8200000_DP
      GSP(14)   =   8.3600000_DP
      GPP(14)   =   7.3100000_DP
      GP2(14)   =   6.5400000_DP
      HSP(14)   =   1.3200000_DP

      ! . Phosphorus.
      SEPAR(15) = .TRUE.
      ALP(15)   =    2.4553220_DP
      EISOL(15) = -124.4368355_DP
      BETAS(15) =   -6.3537640_DP
      BETAP(15) =   -6.5907090_DP
      ZS(15)    =    1.9812800_DP
      ZP(15)    =    1.8751500_DP
      ZD(15)    =    1.0000000_DP
      DD(15)    =    1.0452022_DP
      QQ(15)    =    0.8923660_DP
      AM(15)    =    0.4248440_DP
      AD(15)    =    0.3275319_DP
      AQ(15)    =    0.4386854_DP
      USS(15)   =  -42.0298630_DP
      UPP(15)   =  -34.0307090_DP
      GSS(15)   =   11.5600050_DP
      GSP(15)   =   5.2374490_DP
      GPP(15)   =    7.8775890_DP
      GP2(15)   =    7.3076480_DP
      HSP(15)   =    0.7792380_DP

      ! . Sulphur.
      SEPAR(16) = .TRUE.
      ALP(16)   =    2.4616480_DP
      EISOL(16) = -191.7321930_DP
      BETAS(16) =  -3.9205660_DP
      BETAP(16) =  -7.9052780_DP
      ZS(16)    =    2.3665150_DP
      ZP(16)    =    1.6672630_DP
      ZD(16)    =    1.0000000_DP
      DD(16)    =    0.9004265_DP
      QQ(16)    =    1.0036329_DP
      AM(16)    =    0.4331617_DP
      AD(16)    =    0.5907115_DP
      AQ(16)    =    0.6454943_DP
      USS(16)   =  -56.6940560_DP
      UPP(16)   =  -48.7170490_DP
      GSS(16)   =   11.7863290_DP
      GSP(16)   =    8.6631270_DP
      GPP(16)   =   10.0393080_DP
      GP2(16)   =    7.7816880_DP
      HSP(16)   =    2.5321370_DP

      ! . Chlorine.
      SEPAR(17) = .TRUE.
      ALP(17)   =    2.9193680_DP
      EISOL(17) = -372.1984310_DP
      BETAS(17) =  -24.5946700_DP
      BETAP(17) =  -14.6372160_DP
      ZS(17)    =    3.6313760_DP
      ZP(17)    =    2.0767990_DP
      ZD(17)    =    1.0000000_DP
      DD(17)    =    0.5406286_DP
      QQ(17)    =    0.8057208_DP
      AM(17)    =    0.5523705_DP
      AD(17)    =    0.7693200_DP
      AQ(17)    =    0.6133369_DP
      USS(17)   = -111.6139480_DP
      UPP(17)   =  -76.6401070_DP
      GSS(17)   =   15.0300000_DP
      GSP(17)   =   13.1600000_DP
      GPP(17)   =   11.3000000_DP
      GP2(17)   =    9.9700000_DP
      HSP(17)   =    2.4200000_DP

      FN1(17,1) = 0.0942430_DP ; FN2(17,1) = 4.0000000_DP ; FN3(17,1) = 1.3000000_DP
      FN1(17,2) = 0.0271680_DP ; FN2(17,2) = 4.0000000_DP ; FN3(17,2) = 2.1000000_DP

      ! . Potassium.
      SEPAR(19) = .TRUE.
      ALP(19)   =  1.16_DP
      EISOL(19) =  0.00_DP
      AM(19)    =  0.50_DP

      ! . Zinc (added by Aline).
      SEPAR(30)  = .TRUE.
      USS(30)    =    -21.0400080_DP
      UPP(30)    =    -17.6555740_DP
      BETAS(30)  =     -1.9974290_DP
      BETAP(30)  =     -4.7581190_DP
      ZS(30)     =      1.9542990_DP
      ZP(30)     =      1.3723650_DP
      ZD(30)     =      1.0000000_DP
      ALP(30)    =      1.4845630_DP
      EISOL(30)  =    -30.2800160_DP
      GSS(30)    =     11.8000000_DP
      GSP(30)    =     11.1820180_DP
      GPP(30)    =     13.3000000_DP
      GP2(30)    =     12.9305200_DP
      HSP(30)    =      0.4846060_DP
      DD(30)     =      1.3581113_DP
      QQ(30)     =      1.5457406_DP
      AM(30)     =      0.4336641_DP
      AD(30)     =      0.2317423_DP
      AQ(30)     =      0.2621165_DP

      ! . Bromine.
      SEPAR(35) = .TRUE.
      ALP(35)   =    2.5765460_DP
      EISOL(35) = -352.3142087_DP
      BETAS(35) =  -19.3998800_DP
      BETAP(35) =   -8.9571950_DP
      ZS(35)    =    3.0641330_DP
      ZP(35)    =    2.0383330_DP
      ZD(35)    =    1.0000000_DP
      DD(35)    =    0.8458104_DP
      QQ(35)    =    1.0407133_DP
      AM(35)    =    0.5526071_DP
      AD(35)    =    0.6024598_DP
      AQ(35)    =    0.5307555_DP
      USS(35)   = -104.6560630_DP
      UPP(35)   =  -74.9300520_DP
      GSS(35)   =   15.0364395_DP
      GSP(35)   =   13.0346824_DP
      GPP(35)   =   11.2763254_DP
      GP2(35)   =    9.8544255_DP
      HSP(35)   =    2.4558683_DP

      FN1(35,1) = 0.0666850_DP ; FN2(35,1) = 4.0000000_DP ; FN3(35,1) = 1.5000000_DP
      FN1(35,2) = 0.0255680_DP ; FN2(35,2) = 4.0000000_DP ; FN3(35,2) = 2.3000000_DP

      ! . Iodine.
      SEPAR(53) = .TRUE.
      USS(53)   = -103.5896630_DP
      UPP(53)   =  -74.4299970_DP
      BETAS(53) =   -8.4433270_DP
      BETAP(53) =   -6.3234050_DP
      ZS(53)    =    2.1028580_DP
      ZP(53)    =    2.1611530_DP
      ZD(53)    =    1.0000000_DP
      ALP(53)   =    2.2994240_DP
      EISOL(53) = -346.8642857_DP
      GSS(53)   =   15.0404486_DP
      GSP(53)   =   13.0565580_DP
      GPP(53)   =   11.1477837_DP
      GP2(53)   =    9.9140907_DP
      HSP(53)   =    2.4563820_DP
      DD(53)    =    1.4878778_DP
      QQ(53)    =    1.1887388_DP
      AM(53)    =    0.5527544_DP
      AD(53)    =    0.4497523_DP
      AQ(53)    =    0.4631775_DP

      FN1(53,1) = 0.0043610_DP ; FN2(53,1) = 2.3000000_DP ; FN3(53,1) = 1.8000000_DP
      FN1(53,2) = 0.0157060_DP ; FN2(53,2) = 3.0000000_DP ; FN3(53,2) = 2.2400000_DP

      END SUBROUTINE PARAMETERS_AM1

      !----------------------------
      SUBROUTINE PARAMETERS_AM1TEST
      !----------------------------

      ! . Hydrogen.
      SEPAR(1) = .TRUE.
      ALP(1)   =   2.8823240_DP
      EISOL(1) = -11.3964270_DP
      BETAS(1) =  -6.1737870_DP
      ZS(1)    =   1.1880780_DP
      AM(1)    =   0.4721793_DP
      AD(1)    =   0.4721793_DP
      AQ(1)    =   0.4721793_DP
      USS(1)   = -11.3964270_DP
      GSS(1)   =  12.8480000_DP

      FN1(1,1) =  0.1227960_DP ; FN2(1,1) = 5.0000000_DP ; FN3(1,1) = 1.2000000_DP
      FN1(1,2) =  0.0050900_DP ; FN2(1,2) = 5.0000000_DP ; FN3(1,2) = 1.8000000_DP
      FN1(1,3) = -0.0183360_DP ; FN2(1,3) = 2.0000000_DP ; FN3(1,3) = 2.1000000_DP

      ! . Lithium.
      SEPAR(3) = .TRUE.
      ALP(3)   =   1.2501400_DP
      EISOL(3) =  -5.1280000_DP
      BETAS(3) =  -1.3500400_DP
      BETAP(3) =  -1.3500400_DP
      ZS(3)    =   0.7023800_DP
      ZP(3)    =   0.7023800_DP
      DD(3)    =   2.0549783_DP
      QQ(3)    =   1.7437069_DP
      AM(3)    =   0.2682837_DP
      AD(3)    =   0.2269793_DP
      AQ(3)    =   0.2614581_DP
      USS(3)   =  -5.1280000_DP
      UPP(3)   =  -2.7212000_DP
      GSS(3)   =   7.3000000_DP
      GSP(3)   =   5.4200000_DP
      GPP(3)   =   5.0000000_DP
      GP2(3)   =   4.5200000_DP
      HSP(3)   =   0.8300000_DP

      ! . Beryllium.
      SEPAR(4) = .TRUE.
      ALP(4)   =   1.6694340_DP
      EISOL(4) = -24.2047560_DP
      BETAS(4) =  -4.0170960_DP
      BETAP(4) =  -4.0170960_DP
      ZS(4)    =   1.0042100_DP
      ZP(4)    =   1.0042100_DP
      DD(4)    =   1.4373245_DP
      QQ(4)    =   1.2196103_DP
      AM(4)    =   0.3307607_DP
      AD(4)    =   0.3356142_DP
      AQ(4)    =   0.3846373_DP
      USS(4)   = -16.6023780_DP
      UPP(4)   = -10.7037710_DP
      GSS(4)   =   9.0000000_DP
      GSP(4)   =   7.4300000_DP
      GPP(4)   =   6.9700000_DP
      GP2(4)   =   6.2200000_DP
      HSP(4)   =   1.2800000_DP

      ! . Boron.
      SEPAR(5) = .TRUE.
      ALP(5)   =   2.4469090_DP
      EISOL(5) = -63.7172650_DP
      BETAS(5) =  -9.5991140_DP
      BETAP(5) =  -6.2737570_DP
      ZS(5)    =   1.6117090_DP
      ZP(5)    =   1.5553850_DP
      DD(5)    =   0.9107622_DP
      QQ(5)    =   0.7874223_DP
      AM(5)    =   0.3891951_DP
      AD(5)    =   0.5045152_DP
      AQ(5)    =   0.5678856_DP
      USS(5)   = -34.4928700_DP
      UPP(5)   = -22.6315250_DP
      GSS(5)   =  10.5900000_DP
      GSP(5)   =   9.5600000_DP
      GPP(5)   =   8.8600000_DP
      GP2(5)   =   7.8600000_DP
      HSP(5)   =   1.8100000_DP

      ! . Carbon.
      SEPAR(6) = .TRUE.
      ALP(6)   =    2.6482740_DP
      EISOL(6) = -120.8157940_DP
      BETAS(6) =  -15.7157830_DP
      BETAP(6) =   -7.7192830_DP
      ZS(6)    =    1.8086650_DP
      ZP(6)    =    1.6851160_DP
      DD(6)    =    0.8236736_DP
      QQ(6)    =    0.7268015_DP
      AM(6)    =    0.4494671_DP
      AD(6)    =    0.6082946_DP
      AQ(6)    =    0.6423492_DP
      USS(6)   =  -52.0286580_DP
      UPP(6)   =  -39.6142390_DP
      GSS(6)   =   12.2300000_DP
      GSP(6)   =   11.4700000_DP
      GPP(6)   =   11.0800000_DP
      GP2(6)   =    9.8400000_DP
      HSP(6)   =    2.4300000_DP

      FN1(6,1) =  0.0113550_DP ; FN2(6,1) =  5.0000000_DP ; FN3(6,1) =  1.6000000_DP
      FN1(6,2) =  0.0459240_DP ; FN2(6,2) =  5.0000000_DP ; FN3(6,2) =  1.8500000_DP
      FN1(6,3) = -0.0200610_DP ; FN2(6,3) =  5.0000000_DP ; FN3(6,3) =  2.0500000_DP
      FN1(6,4) = -0.0012600_DP ; FN2(6,4) =  5.0000000_DP ; FN3(6,4) =  2.6500000_DP

      ! . Nitrogen.
      SEPAR(7) = .TRUE.
      ALP(7)   =    2.9472860_DP
      EISOL(7) = -202.4077430_DP
      BETAS(7) =  -20.2991100_DP
      BETAP(7) =  -18.2386660_DP
      ZS(7)    =    2.3154100_DP
      ZP(7)    =    2.1579400_DP
      DD(7)    =    0.6433247_DP
      QQ(7)    =    0.5675528_DP
      AM(7)    =    0.4994487_DP
      AD(7)    =    0.7820840_DP
      AQ(7)    =    0.7883498_DP
      USS(7)   =  -71.8600000_DP
      UPP(7)   =  -57.1675810_DP
      GSS(7)   =   13.5900000_DP
      GSP(7)   =   12.6600000_DP
      GPP(7)   =   12.9800000_DP
      GP2(7)   =   11.5900000_DP
      HSP(7)   =    3.1400000_DP

      FN1(7,1) =  0.0252510_DP ; FN2(7,1) =  5.0000000_DP ; FN3(7,1) =  1.5000000_DP
      FN1(7,2) =  0.0289530_DP ; FN2(7,2) =  5.0000000_DP ; FN3(7,2) =  2.1000000_DP
      FN1(7,3) = -0.0058060_DP ; FN2(7,3) =  2.0000000_DP ; FN3(7,3) =  2.4000000_DP

      ! . Oxygen.
      SEPAR(8) = .TRUE.
      ALP(8)   =    4.4553710_DP
      EISOL(8) = -316.0995200_DP
      BETAS(8) =  -29.2727730_DP
      BETAP(8) =  -29.2727730_DP
      ZS(8)    =    3.1080320_DP
      ZP(8)    =    2.5240390_DP
      DD(8)    =    0.4988896_DP
      QQ(8)    =    0.4852322_DP
      AM(8)    =    0.5667034_DP
      AD(8)    =    0.9961066_DP
      AQ(8)    =    0.9065223_DP
      USS(8)   =  -97.8300000_DP
      UPP(8)   =  -78.2623800_DP
      GSS(8)   =   15.4200000_DP
      GSP(8)   =   14.4800000_DP
      GPP(8)   =   14.5200000_DP
      GP2(8)   =   12.9800000_DP
      HSP(8)   =    3.9400000_DP

      FN1(8,1) = 0.2809620_DP ; FN2(8,1) = 5.0000000_DP ; FN3(8,1) = 0.8479180_DP
      FN1(8,2) = 0.0814300_DP ; FN2(8,2) = 7.0000000_DP ; FN3(8,2) = 1.4450710_DP

      ! . Fluorine.
      SEPAR(9) = .TRUE.
      ALP(9)   =    5.5178000_DP
      EISOL(9) = -482.2905830_DP
      BETAS(9) =  -69.5902770_DP
      BETAP(9) =  -27.9223600_DP
      ZS(9)    =    3.7700820_DP
      ZP(9)    =    2.4946700_DP
      DD(9)    =    0.4145203_DP
      QQ(9)    =    0.4909446_DP
      AM(9)    =    0.6218302_DP
      AD(9)    =    1.2088792_DP
      AQ(9)    =    0.9449355_DP
      USS(9)   = -136.1055790_DP
      UPP(9)   = -104.8898850_DP
      GSS(9)   =   16.9200000_DP
      GSP(9)   =   17.2500000_DP
      GPP(9)   =   16.7100000_DP
      GP2(9)   =   14.9100000_DP
      HSP(9)   =    4.8300000_DP

      FN1(9,1) = 0.2420790_DP ; FN2(9,1) = 4.8000000_DP ; FN3(9,1) = 0.9300000_DP
      FN1(9,2) = 0.0036070_DP ; FN2(9,2) = 4.6000000_DP ; FN3(9,2) = 1.6600000_DP

      ! . Sodium.
      SEPAR(11) = .TRUE.
      ALP(11)   =  1.32_DP
      EISOL(11) =  0.00_DP
      AM(11)    =  0.50_DP

      ! . Magnesium (Hutter et al., JPCB 102, 8080-8090, 1998).
      SEPAR(12) = .TRUE.
      USS  ( 12) =  -14.96959313_DP
      UPP  ( 12) =  -11.56229248_DP
      BETAS( 12) =   -1.25974355_DP
      BETAP( 12) =   -0.77836604_DP
      ZS   ( 12) =    1.22339270_DP
      ZP   ( 12) =    1.02030798_DP
      ALP  ( 12) =    1.67049799_DP
      EISOL( 12) =  -22.43786349_DP ! From CALPAR formula. Otherwise PM3 value?
      GSS  ( 12) =    7.50132277_DP
      GSP  ( 12) =    6.34591536_DP
      GPP  ( 12) =    4.77534467_DP
      GP2  ( 12) =    4.34017279_DP
      HSP  ( 12) =    0.48930466_DP
      DD   ( 12) =    1.75012115_DP ! From formula (15) page 94 Thiel + Dewar.
      QQ   ( 12) =    1.64001467_DP ! From formula (16) page 95 Thiel + Dewar.
      AM   ( 12) =    0.27568257_DP ! From formula (21) page 95 Thiel + Dewar (rho0 = 1/2 Gss = 1/2 AM so AM = Gss/AU_TO_EV).
      AD   ( 12) =    0.19957236_DP ! By solving numerically equation (19) page 95 Thiel + Dewar.
      AQ   ( 12) =    0.26440152_DP ! By solving numerically equation (20) page 95 Thiel + Dewar.

      FN1( 12,1) =  2.55017735_DP ; FN2( 12,1) = 4.29397225_DP ; FN3( 12,1) = 0.79989601_DP
      FN1( 12,2) = -0.00565806_DP ; FN2( 12,2) = 2.96053910_DP ; FN3( 12,2) = 1.47499983_DP
      FN1( 12,3) = -0.00610286_DP ; FN2( 12,3) = 2.61416919_DP ; FN3( 12,3) = 2.42604040_DP

      ! . Aluminium.
      SEPAR(13) = .TRUE.
      ALP(13)   =   1.9765860_DP
      EISOL(13) = -46.4208150_DP
      BETAS(13) =  -3.8668220_DP
      BETAP(13) =  -2.3171460_DP
      ZS(13)    =   1.5165930_DP
      ZP(13)    =   1.3063470_DP
      ZD(13)    =   1.0000000_DP
      DD(13)    =   1.4040443_DP
      QQ(13)    =   1.2809154_DP
      AM(13)    =   0.2973172_DP
      AD(13)    =   0.2630229_DP
      AQ(13)    =   0.3427832_DP
      USS(13)   = -24.3535850_DP
      UPP(13)   = -18.3636450_DP
      GSS(13)   =   8.0900000_DP
      GSP(13)   =   6.6300000_DP
      GPP(13)   =   5.9800000_DP
      GP2(13)   =   5.4000000_DP
      HSP(13)   =   0.7000000_DP

      FN1(13,1) = 0.0900000_DP ; FN2(13,1) = 12.3924430_DP ; FN3(13,1) = 2.0503940_DP

      ! . Silicon.
      SEPAR(14) = .TRUE.
      ALP   (14) =      2.2578160_DP
      EISOL (14) =    -79.0017420_DP
      BETAS (14) =     -3.7848520_DP
      BETAP (14) =     -1.9681230_DP
      ZS    (14) =      1.8306970_DP
      ZP    (14) =      1.2849530_DP
      DD    (14) =      1.1631107_DP
      QQ    (14) =      1.3022422_DP
      AM    (14) =      0.3608967_DP
      AD    (14) =      0.3829813_DP
      AQ    (14) =      0.3712106_DP
      USS   (14) =    -33.9536220_DP
      UPP   (14) =    -28.9347490_DP
      GSS   (14) =      9.8200000_DP
      GSP   (14) =      8.3600000_DP
      GPP   (14) =      7.3100000_DP
      GP2   (14) =      6.5400000_DP
      HSP   (14) =      1.3200000_DP

      FN1(14, 1) =      0.2500000_DP ; FN2(14, 1) =      9.0000000_DP ; FN3(14, 1) =      0.9114530_DP
      FN1(14, 2) =      0.0615130_DP ; FN2(14, 2) =      5.0000000_DP ; FN3(14, 2) =      1.9955690_DP
      FN1(14, 3) =      0.0207890_DP ; FN2(14, 3) =      5.0000000_DP ; FN3(14, 3) =      2.9906100_DP

      ! . Phosphorus.
      SEPAR(15) = .TRUE.
      ALP   (15) =      2.4553220_DP
      EISOL (15) =   -124.4368355_DP
      BETAS (15) =     -6.3537640_DP
      BETAP (15) =     -6.5907090_DP
      ZS    (15) =      1.9812800_DP
      ZP    (15) =      1.8751500_DP
      DD    (15) =      1.0452022_DP
      QQ    (15) =      0.8923660_DP
      AM    (15) =      0.4248440_DP
      AD    (15) =      0.3275319_DP
      AQ    (15) =      0.4386854_DP
      USS   (15) =    -42.0298630_DP
      UPP   (15) =    -34.0307090_DP
      GSS   (15) =     11.5600050_DP
      GSP   (15) =      5.2374490_DP
      GPP   (15) =      7.8775890_DP
      GP2   (15) =      7.3076480_DP
      HSP   (15) =      0.7792380_DP

      FN1(15, 1) =     -0.0318270_DP ; FN2(15, 1) =      6.0000000_DP ; FN3(15, 1) =      1.4743230_DP
      FN1(15, 2) =      0.0184700_DP ; FN2(15, 2) =      7.0000000_DP ; FN3(15, 2) =      1.7793540_DP
      FN1(15, 3) =      0.0332900_DP ; FN2(15, 3) =      9.0000000_DP ; FN3(15, 3) =      3.0065760_DP

      ! . Sulphur.
      SEPAR(16) = .TRUE.
      ALP   (16) =      2.4616480_DP
      EISOL (16) =   -191.7321930_DP
      BETAS (16) =     -3.9205660_DP
      BETAP (16) =     -7.9052780_DP
      ZS    (16) =      2.3665150_DP
      ZP    (16) =      1.6672630_DP
      DD    (16) =      0.9004265_DP
      QQ    (16) =      1.0036329_DP
      AM    (16) =      0.4331617_DP
      AD    (16) =      0.5907115_DP
      AQ    (16) =      0.6454943_DP
      USS   (16) =    -56.6940560_DP
      UPP   (16) =    -48.7170490_DP
      GSS   (16) =     11.7863290_DP
      GSP   (16) =      8.6631270_DP
      GPP   (16) =     10.0393080_DP
      GP2   (16) =      7.7816880_DP
      HSP   (16) =      2.5321370_DP

      FN1(16, 1) =     -0.5091950_DP ; FN2(16, 1) =      4.5936910_DP ; FN3(16, 1) =      0.7706650_DP
      FN1(16, 2) =     -0.0118630_DP ; FN2(16, 2) =      5.8657310_DP ; FN3(16, 2) =      1.5033130_DP
      FN1(16, 3) =      0.0123340_DP ; FN2(16, 3) =     13.5573360_DP ; FN3(16, 3) =      2.0091730_DP

      ! . Chlorine.
      SEPAR(17) = .TRUE.
      ALP(17)   =    2.9193680_DP
      EISOL(17) = -372.1984310_DP
      BETAS(17) =  -24.5946700_DP
      BETAP(17) =  -14.6372160_DP
      ZS(17)    =    3.6313760_DP
      ZP(17)    =    2.0767990_DP
      ZD(17)    =    1.0000000_DP
      DD(17)    =    0.5406286_DP
      QQ(17)    =    0.8057208_DP
      AM(17)    =    0.5523705_DP
      AD(17)    =    0.7693200_DP
      AQ(17)    =    0.6133369_DP
      USS(17)   = -111.6139480_DP
      UPP(17)   =  -76.6401070_DP
      GSS(17)   =   15.0300000_DP
      GSP(17)   =   13.1600000_DP
      GPP(17)   =   11.3000000_DP
      GP2(17)   =    9.9700000_DP
      HSP(17)   =    2.4200000_DP

      FN1(17,1) = 0.0942430_DP ; FN2(17,1) = 4.0000000_DP ; FN3(17,1) = 1.3000000_DP
      FN1(17,2) = 0.0271680_DP ; FN2(17,2) = 4.0000000_DP ; FN3(17,2) = 2.1000000_DP

      ! . Potassium.
      SEPAR(19) = .TRUE.
      ALP(19)   =  1.16_DP
      EISOL(19) =  0.00_DP
      AM(19)    =  0.50_DP

      ! . Zinc (added by Aline).
      SEPAR(30)  = .TRUE.
      USS(30)    =    -21.0400080_DP
      UPP(30)    =    -17.6555740_DP
      BETAS(30)  =     -1.9974290_DP
      BETAP(30)  =     -4.7581190_DP
      ZS(30)     =      1.9542990_DP
      ZP(30)     =      1.3723650_DP
      ZD(30)     =      1.0000000_DP
      ALP(30)    =      1.4845630_DP
      EISOL(30)  =    -30.2800160_DP
      GSS(30)    =     11.8000000_DP
      GSP(30)    =     11.1820180_DP
      GPP(30)    =     13.3000000_DP
      GP2(30)    =     12.9305200_DP
      HSP(30)    =      0.4846060_DP
      DD(30)     =      1.3581113_DP
      QQ(30)     =      1.5457406_DP
      AM(30)     =      0.4336641_DP
      AD(30)     =      0.2317423_DP
      AQ(30)     =      0.2621165_DP

      ! . Bromine.
      SEPAR(35) = .TRUE.
      ALP(35)   =    2.5765460_DP
      EISOL(35) = -352.3142087_DP
      BETAS(35) =  -19.3998800_DP
      BETAP(35) =   -8.9571950_DP
      ZS(35)    =    3.0641330_DP
      ZP(35)    =    2.0383330_DP
      ZD(35)    =    1.0000000_DP
      DD(35)    =    0.8458104_DP
      QQ(35)    =    1.0407133_DP
      AM(35)    =    0.5526071_DP
      AD(35)    =    0.6024598_DP
      AQ(35)    =    0.5307555_DP
      USS(35)   = -104.6560630_DP
      UPP(35)   =  -74.9300520_DP
      GSS(35)   =   15.0364395_DP
      GSP(35)   =   13.0346824_DP
      GPP(35)   =   11.2763254_DP
      GP2(35)   =    9.8544255_DP
      HSP(35)   =    2.4558683_DP

      FN1(35,1) = 0.0666850_DP ; FN2(35,1) = 4.0000000_DP ; FN3(35,1) = 1.5000000_DP
      FN1(35,2) = 0.0255680_DP ; FN2(35,2) = 4.0000000_DP ; FN3(35,2) = 2.3000000_DP

      ! . Iodine.
      SEPAR(53) = .TRUE.
      USS(53)   = -103.5896630_DP
      UPP(53)   =  -74.4299970_DP
      BETAS(53) =   -8.4433270_DP
      BETAP(53) =   -6.3234050_DP
      ZS(53)    =    2.1028580_DP
      ZP(53)    =    2.1611530_DP
      ZD(53)    =    1.0000000_DP
      ALP(53)   =    2.2994240_DP
      EISOL(53) = -346.8642857_DP
      GSS(53)   =   15.0404486_DP
      GSP(53)   =   13.0565580_DP
      GPP(53)   =   11.1477837_DP
      GP2(53)   =    9.9140907_DP
      HSP(53)   =    2.4563820_DP
      DD(53)    =    1.4878778_DP
      QQ(53)    =    1.1887388_DP
      AM(53)    =    0.5527544_DP
      AD(53)    =    0.4497523_DP
      AQ(53)    =    0.4631775_DP

      FN1(53,1) = 0.0043610_DP ; FN2(53,1) = 2.3000000_DP ; FN3(53,1) = 1.8000000_DP
      FN1(53,2) = 0.0157060_DP ; FN2(53,2) = 3.0000000_DP ; FN3(53,2) = 2.2400000_DP

      END SUBROUTINE PARAMETERS_AM1TEST

      !-------------------------
      SUBROUTINE PARAMETERS_MNDO
      !-------------------------

      ! . Hydrogen.
      SEPAR(1) = .TRUE.
      ALP(1)   =   2.5441341_DP
      EISOL(1) = -11.9062760_DP
      BETAS(1) =  -6.9890640_DP
      ZS(1)    =   1.3319670_DP
      AM(1)    =   0.4721793_DP
      AD(1)    =   0.4721793_DP
      AQ(1)    =   0.4721793_DP
      USS(1)   = -11.9062760_DP
      GSS(1)   =  12.8480000_DP

      ! . Lithium.
      SEPAR(3) = .TRUE.
      ALP(3)   =   1.2501400_DP
      EISOL(3) =  -5.1280000_DP
      BETAS(3) =  -1.3500400_DP
      BETAP(3) =  -1.3500400_DP
      ZS(3)    =   0.7023800_DP
      ZP(3)    =   0.7023800_DP
      DD(3)    =   2.0549783_DP
      QQ(3)    =   1.7437069_DP
      AM(3)    =   0.2682837_DP
      AD(3)    =   0.2269793_DP
      AQ(3)    =   0.2614581_DP
      USS(3)   =  -5.1280000_DP
      UPP(3)   =  -2.7212000_DP
      GSS(3)   =   7.3000000_DP
      GSP(3)   =   5.4200000_DP
      GPP(3)   =   5.0000000_DP
      GP2(3)   =   4.5200000_DP
      HSP(3)   =   0.8300000_DP

      ! . Beryllium.
      SEPAR(4) = .TRUE.
      ALP(4)   =   1.6694340_DP
      EISOL(4) = -24.2047560_DP
      BETAS(4) =  -4.0170960_DP
      BETAP(4) =  -4.0170960_DP
      ZS(4)    =   1.0042100_DP
      ZP(4)    =   1.0042100_DP
      DD(4)    =   1.4373245_DP
      QQ(4)    =   1.2196103_DP
      AM(4)    =   0.3307607_DP
      AD(4)    =   0.3356142_DP
      AQ(4)    =   0.3846373_DP
      USS(4)   = -16.6023780_DP
      UPP(4)   = -10.7037710_DP
      GSS(4)   =   9.0000000_DP
      GSP(4)   =   7.4300000_DP
      GPP(4)   =   6.9700000_DP
      GP2(4)   =   6.2200000_DP
      HSP(4)   =   1.2800000_DP

      ! . Boron.
      SEPAR(5) = .TRUE.
      ALP(5)   =   2.1349930_DP
      EISOL(5) = -64.3159500_DP
      BETAS(5) =  -8.2520540_DP
      BETAP(5) =  -8.2520540_DP
      ZS(5)    =   1.5068010_DP
      ZP(5)    =   1.5068010_DP
      DD(5)    =   0.9579073_DP
      QQ(5)    =   0.8128113_DP
      AM(5)    =   0.3891951_DP
      AD(5)    =   0.4904730_DP
      AQ(5)    =   0.5556979_DP
      USS(5)   = -34.5471300_DP
      UPP(5)   = -23.1216900_DP
      GSS(5)   =  10.5900000_DP
      GSP(5)   =   9.5600000_DP
      GPP(5)   =   8.8600000_DP
      GP2(5)   =   7.8600000_DP
      HSP(5)   =   1.8100000_DP

      ! . Carbon.
      SEPAR(6) = .TRUE.
      ALP(6)   =    2.5463800_DP
      EISOL(6) = -120.5006060_DP
      BETAS(6) =  -18.9850440_DP
      BETAP(6) =   -7.9341220_DP
      ZS(6)    =    1.7875370_DP
      ZP(6)    =    1.7875370_DP
      DD(6)    =    0.8074662_DP
      QQ(6)    =    0.6851578_DP
      AM(6)    =    0.4494671_DP
      AD(6)    =    0.6149474_DP
      AQ(6)    =    0.6685897_DP
      USS(6)   =  -52.2797450_DP
      UPP(6)   =  -39.2055580_DP
      GSS(6)   =   12.2300000_DP
      GSP(6)   =   11.4700000_DP
      GPP(6)   =   11.0800000_DP
      GP2(6)   =    9.8400000_DP
      HSP(6)   =    2.4300000_DP

      ! . Nitrogen.
      SEPAR(7) = .TRUE.
      ALP(7)   =    2.8613420_DP
      EISOL(7) = -202.5812010_DP
      BETAS(7) =  -20.4957580_DP
      BETAP(7) =  -20.4957580_DP
      ZS(7)    =    2.2556140_DP
      ZP(7)    =    2.2556140_DP
      DD(7)    =    0.6399037_DP
      QQ(7)    =    0.5429763_DP
      AM(7)    =    0.4994487_DP
      AD(7)    =    0.7843643_DP
      AQ(7)    =    0.8144720_DP
      USS(7)   =  -71.9321220_DP
      UPP(7)   =  -57.1723190_DP
      GSS(7)   =   13.5900000_DP
      GSP(7)   =   12.6600000_DP
      GPP(7)   =   12.9800000_DP
      GP2(7)   =   11.5900000_DP
      HSP(7)   =    3.1400000_DP

      ! . Oxygen.
      SEPAR(8) = .TRUE.
      ALP(8)   =    3.1606040_DP
      EISOL(8) = -317.8685060_DP
      BETAS(8) =  -32.6880820_DP
      BETAP(8) =  -32.6880820_DP
      ZS(8)    =    2.6999050_DP
      ZP(8)    =    2.6999050_DP
      DD(8)    =    0.5346024_DP
      QQ(8)    =    0.4536252_DP
      AM(8)    =    0.5667034_DP
      AD(8)    =    0.9592562_DP
      AQ(8)    =    0.9495934_DP
      USS(8)   =  -99.6443090_DP
      UPP(8)   =  -77.7974720_DP
      GSS(8)   =   15.4200000_DP
      GSP(8)   =   14.4800000_DP
      GPP(8)   =   14.5200000_DP
      GP2(8)   =   12.9800000_DP
      HSP(8)   =    3.9400000_DP

      ! . Fluorine.
      SEPAR(9) = .TRUE.
      ALP(9)   =    3.4196606_DP
      EISOL(9) = -476.6837810_DP
      BETAS(9) =  -48.2904660_DP
      BETAP(9) =  -36.5085400_DP
      ZS(9)    =    2.8484870_DP
      ZP(9)    =    2.8484870_DP
      DD(9)    =    0.5067166_DP
      QQ(9)    =    0.4299633_DP
      AM(9)    =    0.6218302_DP
      AD(9)    =    1.0850301_DP
      AQ(9)    =    1.0343643_DP
      USS(9)   = -131.0715480_DP
      UPP(9)   = -105.7821370_DP
      GSS(9)   =   16.9200000_DP
      GSP(9)   =   17.2500000_DP
      GPP(9)   =   16.7100000_DP
      GP2(9)   =   14.9100000_DP
      HSP(9)   =    4.8300000_DP

      ! . Sodium.
      SEPAR(11) = .TRUE.
      ALP(11)   =  1.66_DP
      EISOL(11) =  0.00_DP
      AM(11)    =  0.50_DP

      ! . Aluminium.
      SEPAR(13) = .TRUE.
      ALP(13)   =   1.8688394_DP
      EISOL(13) = -44.4840720_DP
      BETAS(13) =  -2.6702840_DP
      BETAP(13) =  -2.6702840_DP
      ZS(13)    =   1.4441610_DP
      ZP(13)    =   1.4441610_DP
      ZD(13)    =   1.0000000_DP
      DD(13)    =   1.3992387_DP
      QQ(13)    =   1.1586797_DP
      AM(13)    =   0.2973172_DP
      AD(13)    =   0.2635574_DP
      AQ(13)    =   0.3673560_DP
      USS(13)   = -23.8070970_DP
      UPP(13)   = -17.5198780_DP
      GSS(13)   =   8.0900000_DP
      GSP(13)   =   6.6300000_DP
      GPP(13)   =   5.9800000_DP
      GP2(13)   =   5.4000000_DP
      HSP(13)   =   0.7000000_DP

      ! . Silicon.
      SEPAR(14) = .TRUE.
      ALP(14)   =   2.2053160_DP
      EISOL(14) = -82.8394220_DP
      BETAS(14) =  -9.0868040_DP
      BETAP(14) =  -1.0758270_DP
      ZS(14)    =   1.3159860_DP
      ZP(14)    =   1.7099430_DP
      ZD(14)    =   1.0000000_DP
      DD(14)    =   1.2580349_DP
      QQ(14)    =   0.9785824_DP
      AM(14)    =   0.3608967_DP
      AD(14)    =   0.3664244_DP
      AQ(14)    =   0.4506740_DP
      USS(14)   = -37.0375330_DP
      UPP(14)   = -27.7696780_DP
      GSS(14)   =   9.8200000_DP
      GSP(14)   =   8.3600000_DP
      GPP(14)   =   7.3100000_DP
      GP2(14)   =   6.5400000_DP
      HSP(14)   =   1.3200000_DP

      ! . Phosphorus.
      SEPAR(15) = .TRUE.
      ALP(15)   =    2.4152800_DP
      EISOL(15) = -152.9599600_DP
      BETAS(15) =   -6.7916000_DP
      BETAP(15) =   -6.7916000_DP
      ZS(15)    =    2.1087200_DP
      ZP(15)    =    1.7858100_DP
      ZD(15)    =    1.0000000_DP
      DD(15)    =    1.0129699_DP
      QQ(15)    =    0.9370090_DP
      AM(15)    =    0.4248438_DP
      AD(15)    =    0.4882420_DP
      AQ(15)    =    0.4979406_DP
      USS(15)   =  -56.1433600_DP
      UPP(15)   =  -42.8510800_DP
      GSS(15)   =   11.5600000_DP
      GSP(15)   =   10.0800000_DP
      GPP(15)   =    8.6400000_DP
      GP2(15)   =    7.6800000_DP
      HSP(15)   =    1.9200000_DP

      ! . Sulphur.
      SEPAR(16) = .TRUE.
      ALP(16)   =    2.4780260_DP
      EISOL(16) = -226.0123900_DP
      BETAS(16) =  -10.7616700_DP
      BETAP(16) =  -10.1084330_DP
      ZS(16)    =    2.3129620_DP
      ZP(16)    =    2.0091460_DP
      ZD(16)    =    1.0000000_DP
      DD(16)    =    0.9189935_DP
      QQ(16)    =    0.8328514_DP
      AM(16)    =    0.4733554_DP
      AD(16)    =    0.5544502_DP
      AQ(16)    =    0.5585244_DP
      USS(16)   =  -72.2422810_DP
      UPP(16)   =  -56.9732070_DP
      GSS(16)   =   12.8800000_DP
      GSP(16)   =   11.2600000_DP
      GPP(16)   =    9.9000000_DP
      GP2(16)   =    8.8300000_DP
      HSP(16)   =    2.2600000_DP

      ! . Chlorine.
      SEPAR(17) = .TRUE.
      ALP(17)   =    2.5422010_DP
      EISOL(17) = -353.1176670_DP
      BETAS(17) =  -14.2623200_DP
      BETAP(17) =  -14.2623200_DP
      ZS(17)    =    3.7846450_DP
      ZP(17)    =    2.0362630_DP
      ZD(17)    =    1.0000000_DP
      DD(17)    =    0.4986870_DP
      QQ(17)    =    0.8217603_DP
      AM(17)    =    0.5523705_DP
      AD(17)    =    0.8061220_DP
      AQ(17)    =    0.6053435_DP
      USS(17)   = -100.2271660_DP
      UPP(17)   =  -77.3786670_DP
      GSS(17)   =   15.0300000_DP
      GSP(17)   =   13.1600000_DP
      GPP(17)   =   11.3000000_DP
      GP2(17)   =    9.9700000_DP
      HSP(17)   =    2.4200000_DP

      ! . Potassium.
      SEPAR(19) = .TRUE.
      ALP(19)   =  1.16_DP
      EISOL(19) =  0.00_DP
      AM(19)    =  0.50_DP

      ! . Chromium.
      SEPAR(24) = .TRUE.
      ALP(24)   =    3.0683070_DP
      EISOL(24) = -134.8187920_DP
      BETAS(24) =   -0.1000000_DP
      BETAP(24) =   -0.1000000_DP
      BETAD(24) =   -8.7766360_DP
      ZS(24)    =    1.5000000_DP
      ZP(24)    =    1.5000000_DP
      ZD(24)    =    2.8845490_DP
      DD(24)    =    1.7320508_DP
      QQ(24)    =    1.4142136_DP
      AM(24)    =    0.2205072_DP
      AD(24)    =    0.2711332_DP
      AQ(24)    =    0.4464656_DP
      USS(24)   =  -17.5170270_DP
      UPP(24)   =  -12.5337290_DP
      UDD(24)   =  -44.1249280_DP
      GSS(24)   =    6.0000000_DP
      GSP(24)   =    4.1500000_DP
      GPP(24)   =    5.0000000_DP
      GP2(24)   =    3.5000000_DP
      HSP(24)   =    1.0000000_DP
      GSD(24)   =    2.8746410_DP
      GPD(24)   =    3.0000000_DP
      GDD(24)   =    8.8949670_DP

      ! . Zinc (added by Aline).
      SEPAR(30) = .TRUE.
      ALP(30)   =    1.5064570_DP
      EISOL(30) =  -29.8794320_DP
      BETAS(30) =   -1.0000000_DP
      BETAP(30) =   -2.0000000_DP
      ZS(30)    =    2.0473590_DP
      ZP(30)    =    1.4609460_DP
      ZD(30)    =    1.0000000_DP
      DD(30)    =    1.3037826_DP
      QQ(30)    =    1.4520183_DP
      AM(30)    =    0.4336641_DP
      AD(30)    =    0.2375912_DP
      AQ(30)    =    0.2738858_DP
      USS(30)   =  -20.8397160_DP
      UPP(30)   =  -19.6252240_DP
      GSS(30)   =   11.8000000_DP
      GSP(30)   =   11.1820180_DP
      GPP(30)   =   13.3000000_DP
      GP2(30)   =   12.9305200_DP
      HSP(30)   =    0.4846060_DP

      ! . Germanium.
      SEPAR(32) = .TRUE.
      ALP(32)   =   1.9784980_DP
      EISOL(32) = -76.2489440_DP
      BETAS(32) =  -4.5164790_DP
      BETAP(32) =  -1.7555170_DP
      ZS(32)    =   1.2931800_DP
      ZP(32)    =   2.0205640_DP
      DD(32)    =   1.2556091_DP
      QQ(32)    =   1.0498655_DP
      AM(32)    =   0.3601617_DP
      AD(32)    =   0.3643722_DP
      AQ(32)    =   0.4347337_DP
      USS(32)   = -33.9493670_DP
      UPP(32)   = -27.4251050_DP
      GSS(32)   =   9.8000000_DP
      GSP(32)   =   8.3000000_DP
      GPP(32)   =   7.3000000_DP
      GP2(32)   =   6.5000000_DP
      HSP(32)   =   1.3000000_DP

      ! . Bromine.
      SEPAR(35) = .TRUE.
      ALP(35)   =    2.44570510_DP
      EISOL(35) = -346.68125000_DP
      BETAS(35) =   -8.91710700_DP
      BETAP(35) =   -9.94374000_DP
      ZS(35)    =    3.85430190_DP
      ZP(35)    =    2.19920910_DP
      ZD(35)    =    1.00000000_DP
      DD(35)    =    0.60510740_DP
      QQ(35)    =    0.96458730_DP
      AM(35)    =    0.55260680_DP
      AD(35)    =    0.72583300_DP
      AQ(35)    =    0.55745890_DP
      USS(35)   =  -99.98644050_DP
      UPP(35)   =  -75.67130750_DP
      GSS(35)   =   15.03643948_DP
      GSP(35)   =   13.03468242_DP
      GPP(35)   =   11.27632539_DP
      GP2(35)   =    9.85442552_DP
      HSP(35)   =    2.45586832_DP

      ! . Tin.
      SEPAR(50) = .TRUE.
      ALP(50)   =   1.8008140_DP
      EISOL(50) = -92.3241020_DP
      BETAS(50) =  -3.2351470_DP
      BETAP(50) =  -4.2904160_DP
      ZS(50)    =   2.0803800_DP
      ZP(50)    =   1.9371060_DP
      DD(50)    =   1.5697766_DP
      QQ(50)    =   1.3262292_DP
      AM(50)    =   0.3601617_DP
      AD(50)    =   0.3219998_DP
      AQ(50)    =   0.3713827_DP
      USS(50)   = -40.8518020_DP
      UPP(50)   = -28.5602490_DP
      GSS(50)   =   9.8000000_DP
      GSP(50)   =   8.3000000_DP
      GPP(50)   =   7.3000000_DP
      GP2(50)   =   6.5000000_DP
      HSP(50)   =   1.3000000_DP

      ! . Iodine.
      SEPAR(53) = .TRUE.
      USS(53)   = -100.00305380_DP
      UPP(53)   =  -74.61146920_DP
      BETAS(53) =   -7.41445100_DP
      BETAP(53) =   -6.19678100_DP
      ZS(53)    =    2.27296100_DP
      ZP(53)    =    2.16949800_DP
      ZD(53)    =    1.00000000_DP
      ALP(53)   =    2.20732000_DP
      EISOL(53) = -340.59836000_DP
      GSS(53)   =   15.04044855_DP
      GSP(53)   =   13.05655798_DP
      GPP(53)   =   11.14778369_DP
      GP2(53)   =    9.91409071_DP
      HSP(53)   =    2.45638202_DP
      DD(53)    =    1.42532330_DP
      QQ(53)    =    1.18417070_DP
      AM(53)    =    0.55275410_DP
      AD(53)    =    0.45934510_DP
      AQ(53)    =    0.45853760_DP

      ! . Mercury.
      SEPAR(80) = .TRUE.
      ALP(80)   =   1.3356410_DP
      EISOL(80) = -28.8191480_DP
      BETAS(80) =  -0.4045250_DP
      BETAP(80) =  -6.2066830_DP
      ZS(80)    =   2.2181840_DP
      ZP(80)    =   2.0650380_DP
      DD(80)    =   1.7378048_DP
      QQ(80)    =   1.4608064_DP
      AM(80)    =   0.3969129_DP
      AD(80)    =   0.3047694_DP
      AQ(80)    =   0.3483102_DP
      USS(80)   = -19.8095740_DP
      UPP(80)   = -13.1025300_DP
      GSS(80)   =  10.8000000_DP
      GSP(80)   =   9.3000000_DP
      GPP(80)   =  14.3000000_DP
      GP2(80)   =  13.5000000_DP
      HSP(80)   =   1.3000000_DP

      ! . Lead.
      SEPAR(82) = .TRUE.
      ALP(82)   =    1.7283330_DP
      EISOL(82) = -105.8345040_DP
      BETAS(82) =   -8.0423870_DP
      BETAP(82) =   -3.0000000_DP
      ZS(82)    =    2.4982860_DP
      ZP(82)    =    2.0820710_DP
      DD(82)    =    1.5526624_DP
      QQ(82)    =    1.4488558_DP
      AM(82)    =    0.3601617_DP
      AD(82)    =    0.3239309_DP
      AQ(82)    =    0.3502057_DP
      USS(82)   =  -47.3196920_DP
      UPP(82)   =  -28.8475600_DP
      GSS(82)   =    9.8000000_DP
      GSP(82)   =    8.3000000_DP
      GPP(82)   =    7.3000000_DP
      GP2(82)   =    6.5000000_DP
      HSP(82)   =    1.3000000_DP

      END SUBROUTINE PARAMETERS_MNDO
     
      !------------------------
      SUBROUTINE PARAMETERS_PDDG
      !------------------------

      ! . Hydrogen.
      SEPAR(1)   = .TRUE.
      USS  (  1) =  -12.8932720_DP 
      BETAS(  1) =   -6.1526540_DP
      ZS   (  1) =    0.9727860_DP
      ALP  (  1) =    3.3816860_DP
      EISOL(  1) =  -13.1205660_DP
      GSS  (  1) =   14.7942080_DP
      AM   (  1) =    0.5437048_DP
      AD   (  1) =    0.5437048_DP
      AQ   (  1) =    0.5437048_DP

      FN1(  1,1) =  1.1222440_DP ; FN2(  1,1) = 4.7077900_DP ; FN3(  1,1) = 1.5470990_DP
      FN1(  1,2) = -1.0697370_DP ; FN2(  1,2) = 5.8579950_DP ; FN3(  1,2) = 1.5678930_DP
      PDDGC ( 1,1) =  0.0571930_DP ; PDDGE(1,1) =  0.6633950_DP
      PDDGC ( 1,2) = -0.0348230_DP ; PDDGE(1,2) =  1.0819010_DP

      ! . Carbon.
      SEPAR(6)   = .TRUE.
      USS  (  6) =  -48.2412410_DP
      UPP  (  6) =  -36.4612560_DP
      BETAS(  6) =  -11.9528180_DP
      BETAP(  6) =   -9.9224110_DP
      ZS   (  6) =    1.5678640_DP
      ZP   (  6) =    1.8466590_DP
      ALP  (  6) =    2.7257720_DP
      EISOL(  6) = -113.4282420_DP
      GSS  (  6) =   11.2007080_DP
      GSP  (  6) =   10.2650270_DP
      GPP  (  6) =   10.7962920_DP
      GP2  (  6) =    9.0425660_DP
      HSP  (  6) =    2.2909800_DP
      DD   (  6) =    0.8314130_DP
      QQ   (  6) =    0.6632220_DP
      AM   (  6) =    0.4116388_DP
      AD   (  6) =    0.5892981_DP
      AQ   (  6) =    0.7659489_DP

      FN1(  6,1) = 0.0489060_DP ; FN2(  6,1) = 5.7653400_DP ; FN3(  6,1) = 1.6822320_DP
      FN1(  6,2) = 0.0476970_DP ; FN2(  6,2) = 5.9737210_DP ; FN3(  6,2) = 0.8944060_DP
      PDDGC ( 6,1) = -0.0007430_DP ; PDDGE(6,1) =  0.8369150_DP
      PDDGC ( 6,2) =  0.0009850_DP ; PDDGE(6,2) =  1.5852360_DP

      ! . Nitrogen.                      
      SEPAR(7)  = .TRUE.
      USS  (  7) =  -49.4545460_DP
      UPP  (  7) =  -47.7574060_DP
      BETAS(  7) =  -14.1172300_DP
      BETAP(  7) =  -19.9385090_DP
      ZS   (  7) =    2.0358070_DP
      ZP   (  7) =    2.3243270_DP
      ALP  (  7) =    2.8491240_DP
      EISOL(  7) = -158.4162050_DP
      GSS  (  7) =   11.9047870_DP
      GSP  (  7) =    7.3485650_DP
      GPP  (  7) =   11.7546720_DP
      GP2  (  7) =   10.8072770_DP
      HSP  (  7) =    1.1367130_DP
      DD   (  7) =    0.6548550_DP
      QQ   (  7) =    0.5269240_DP
      AM   (  7) =    0.4375150_DP
      AD   (  7) =    0.5044213_DP
      AQ   (  7) =    0.7388754_DP

      FN1(  7,1) =  1.5133200_DP ; FN2(  7,1) = 5.9043940_DP ; FN3(  7,1) = 1.7283760_DP
      FN1(  7,2) = -1.5118920_DP ; FN2(  7,2) = 6.0300140_DP ; FN3(  7,2) = 1.7341080_DP
      PDDGC( 7,1) = -0.0031600_DP ; PDDGE(7,1) =  1.0041720_DP
      PDDGC( 7,2) =  0.0125010_DP ; PDDGE(7,2) =  1.5163360_DP
    
    
      ! . Oxygen.
      SEPAR(8)   = .TRUE.
      USS  (  8) =  -87.4125050_DP
      UPP  (  8) =  -72.1830698_DP
      BETAS(  8) =  -44.8745530_DP
      BETAP(  8) =  -24.6019390_DP
      ZS   (  8) =    3.8145650_DP
      ZP   (  8) =    2.3180110_DP
      ALP  (  8) =    3.2253090_DP
      EISOL(  8) = -292.1887660_DP
      GSS  (  8) =   15.7557600_DP
      GSP  (  8) =   10.6211600_DP
      GPP  (  8) =   13.6540160_DP
      GP2  (  8) =   12.4060950_DP
      HSP  (  8) =    0.5938830_DP
      DD   (  8) =    0.4037410_DP
      QQ   (  8) =    0.5283600_DP
      AM   (  8) =    0.5790428_DP
      AD   (  8) =    0.5340363_DP
      AQ   (  8) =    0.8009086_DP

      FN1(  8,1) = -1.1384550_DP ; FN2(  8,1) = 6.0000430_DP ; FN3(  8,1) = 1.6223620_DP
      FN1(  8,2) =  1.1460070_DP ; FN2(  8,2) = 5.9634940_DP ; FN3(  8,2) = 1.6147880_DP
      PDDGC( 8,1) = -0.0010000_DP ; PDDGE(8,1) =  1.3606850_DP
      PDDGC( 8,2) = -0.0015220_DP ; PDDGE(8,2) =  1.3664070_DP

      ! . Fluorine.
      SEPAR(9)   = .TRUE.
      USS  (  9) = -111.4004320_DP
      UPP  (  9) = -106.3952640_DP
      BETAS(  9) =  -50.9373010_DP
      BETAP(  9) =  -31.6369760_DP
      ZS   (  9) =    5.5380330_DP
      ZP   (  9) =    2.5380660_DP
      ALP  (  9) =    3.2005710_DP
      EISOL(  9) = -442.4571330_DP
      GSS  (  9) =   10.4966670_DP
      GSP  (  9) =   16.0736890_DP
      GPP  (  9) =   14.8172560_DP
      GP2  (  9) =   14.4183930_DP
      HSP  (  9) =    0.7277630_DP
      DD   (  9) =    0.2466010_DP
      QQ   (  9) =    0.4825510_DP
      AM   (  9) =    0.3857650_DP
      AD   (  9) =    0.7878570_DP
      AQ   (  9) =    0.6205000_DP

      FN1(  9,1) = -0.0080790_DP ; FN2(  9,1) = 5.9389690_DP ; FN3(  9,1) = 1.8639490_DP
      FN1(  9,2) = -0.0026590_DP ; FN2(  9,2) = 5.9251050_DP ; FN3(  9,2) = 2.3888640_DP
      
      PDDGC(  9,1) =    -0.0128660000_DP
      PDDGC(  9,2) =     0.0073150000_DP
      PDDGE(  9,1) =     1.3056810000_DP
      PDDGE(  9,2) =     1.8425720000_DP

      ! . Silicon.
      SEPAR( 14)   = .TRUE.
      USS  ( 14)   =   -26.3325220000_DP
      UPP  ( 14)   =   -22.6025400000_DP
      BETAS( 14)   =    -3.3764450000_DP
      BETAP( 14)   =    -3.6209690000_DP
      ZS   ( 14)   =     1.5863890000_DP
      ZP   ( 14)   =     1.4859580000_DP
      ALP  ( 14)   =     2.2151570000_DP
      EISOL( 14)   =   -66.8390000000_DP
      DD   ( 14)   =     1.3105150000_DP
      QQ   ( 14)   =     1.1260890000_DP
      AM   ( 14)   =     0.1854904888_DP
      AD   ( 14)   =     0.3066060731_DP
      AQ   ( 14)   =     0.5267593763_DP
      FN1  ( 14,1) =    -0.0713140000_DP
      FN2  ( 14,1) =     6.0000000000_DP
      FN3  ( 14,1) =     0.2379950000_DP
      FN1  ( 14,2) =     0.0894510000_DP
      FN2  ( 14,2) =     6.0000000000_DP
      FN3  ( 14,2) =     1.8977280000_DP
      PDDGC( 14,1) =    -0.0919280000_DP
      PDDGC( 14,2) =    -0.0407530000_DP
      PDDGE( 14,1) =     1.1631900000_DP
      PDDGE( 14,2) =     2.1905260000_DP

      GSS  ( 14) =    5.0471960_DP
      GSP  ( 14) =    5.9490570_DP
      GPP  ( 14) =    6.7593670_DP
      GP2  ( 14) =    5.1612970_DP
      HSP  ( 14) =    0.9198320_DP

      ! . Phosphorus.
      SEPAR( 15)   = .TRUE.
      USS  ( 15)   =   -37.8821130000_DP
      UPP  ( 15)   =   -30.3129790000_DP
      BETAS( 15)   =   -12.6762970000_DP
      BETAP( 15)   =    -7.0933180000_DP
      ZS   ( 15)   =     2.3958820000_DP
      ZP   ( 15)   =     1.7422130000_DP
      ALP  ( 15)   =     2.0052940000_DP
      EISOL( 15)   =  -117.2128540000_DP
      DD   ( 15)   =     0.8939780000_DP
      QQ   ( 15)   =     0.9604570000_DP
      AM   ( 15)   =     0.2867186201_DP
      AD   ( 15)   =     0.4758048477_DP
      AQ   ( 15)   =     0.4135967448_DP
      FN1  ( 15,1) =    -0.3980550000_DP
      FN2  ( 15,1) =     1.9972720000_DP
      FN3  ( 15,1) =     0.9500730000_DP
      FN1  ( 15,2) =    -0.0796530000_DP
      FN2  ( 15,2) =     1.9983600000_DP
      FN3  ( 15,2) =     2.3369590000_DP
      PDDGC( 15,1) =     0.4627410000_DP
      PDDGC( 15,2) =    -0.0204440000_DP
      PDDGE( 15,1) =     0.7142960000_DP
      PDDGE( 15,2) =     2.0412090000_DP

      GSS  ( 15) =    7.8016150_DP
      GSP  ( 15) =    5.1869490_DP
      GPP  ( 15) =    6.6184780_DP
      GP2  ( 15) =    6.0620020_DP
      HSP  ( 15) =    1.5428090_DP

      ! . Sulfur.
      SEPAR( 16)   = .TRUE.
      USS  ( 16)   =   -43.9063660000_DP
      UPP  ( 16)   =   -43.4613480000_DP
      BETAS( 16)   =    -2.9539120000_DP
      BETAP( 16)   =    -8.5077790000_DP
      ZS   ( 16)   =     1.0120020000_DP
      ZP   ( 16)   =     1.8769990000_DP
      ALP  ( 16)   =     2.5397510000_DP
      EISOL( 16)   =  -166.3365540000_DP
      DD   ( 16)   =     1.0069890000_DP
      QQ   ( 16)   =     0.8914870000_DP
      AM   ( 16)   =     0.3294621530_DP
      AD   ( 16)   =     0.7025708472_DP
      AQ   ( 16)   =     0.6628345989_DP
      FN1  ( 16,1) =    -0.3306920000_DP
      FN2  ( 16,1) =     6.0000000000_DP
      FN3  ( 16,1) =     0.8238370000_DP
      FN1  ( 16,2) =     0.0241710000_DP
      FN2  ( 16,2) =     6.0000000000_DP
      FN3  ( 16,2) =     2.0177560000_DP
      PDDGC( 16,1) =     0.1204340000_DP
      PDDGC( 16,2) =    -0.0026630000_DP
      PDDGE( 16,1) =     0.6728700000_DP
      PDDGE( 16,2) =     2.0323400000_DP

      GSS  ( 16) =    8.9646670_DP
      GSP  ( 16) =    6.7859360_DP
      GPP  ( 16) =    9.9681640_DP
      GP2  ( 16) =    7.9702470_DP
      HSP  ( 16) =    4.0418360_DP

      ! . Chlorine.
      SEPAR(17)  = .TRUE.
      USS  ( 17) =  -95.0944340_DP
      UPP  ( 17) =  -53.9216510_DP
      BETAS( 17) =  -26.9131290_DP
      BETAP( 17) =  -14.9911780_DP
      ZS   ( 17) =    2.5482680_DP
      ZP   ( 17) =    2.2846240_DP
      ZD   ( 17) =    1.0000000_DP
      ALP  ( 17) =    2.4976170_DP
      EISOL( 17) = -305.7152010_DP
      GSS  ( 17) =   16.0136010_DP
      GSP  ( 17) =    8.0481150_DP
      GPP  ( 17) =    7.5222150_DP
      GP2  ( 17) =    7.5041540_DP
      HSP  ( 17) =    3.4811530_DP
      DD   ( 17) =    0.8275610_DP
      QQ   ( 17) =    0.7324270_DP
      AM   ( 17) =    0.5885190_DP
      AD   ( 17) =    0.7182216_DP
      AQ   ( 17) =    0.2174760_DP

      FN1( 17,1) = -0.1122220_DP ; FN2( 17,1) = 5.9637190_DP ; FN3( 17,1) = 1.0277190_DP
      FN1( 17,2) = -0.0130610_DP ; FN2( 17,2) = 1.9995560_DP ; FN3( 17,2) = 2.2863770_DP
      
      PDDGC(  17,1) =    -0.0165520000_DP
      PDDGC(  17,2) =    -0.0166460000_DP
      PDDGE(  17,1) =     1.7276900000_DP
      PDDGE(  17,2) =     1.7846550000_DP

      END SUBROUTINE PARAMETERS_PDDG

      !------------------------
      SUBROUTINE PARAMETERS_PM3
      !------------------------

      ! . Hydrogen.
      SEPAR(1)   = .TRUE.
      USS  (  1) =  -13.0733210_DP
      BETAS(  1) =   -5.6265120_DP
      ZS   (  1) =    0.9678070_DP
      ALP  (  1) =    3.3563860_DP
      EISOL(  1) =  -13.0733210_DP
      GSS  (  1) =   14.7942080_DP
      AM   (  1) =    0.5437048_DP
      AD   (  1) =    0.5437048_DP
      AQ   (  1) =    0.5437048_DP

      FN1(  1,1) =  1.1287500_DP ; FN2(  1,1) = 5.0962820_DP ; FN3(  1,1) = 1.5374650_DP
      FN1(  1,2) = -1.0603290_DP ; FN2(  1,2) = 6.0037880_DP ; FN3(  1,2) = 1.5701890_DP

      ! . Beryllium.
      SEPAR(4)   = .TRUE.
      USS  (  4) =  -17.2647520_DP
      UPP  (  4) =  -11.3042430_DP
      BETAS(  4) =   -3.9620530_DP
      BETAP(  4) =   -2.7806840_DP
      ZS   (  4) =    0.8774390_DP
      ZP   (  4) =    1.5087550_DP
      ALP  (  4) =    1.5935360_DP
      EISOL(  4) =  -25.5166530_DP
      GSS  (  4) =    9.0128510_DP
      GSP  (  4) =    6.5761990_DP
      GPP  (  4) =    6.0571820_DP
      GP2  (  4) =    9.0052190_DP
      HSP  (  4) =    0.5446790_DP
      DD   (  4) =    1.0090531_DP
      QQ   (  4) =    0.8117586_DP
      AM   (  4) =    0.3312330_DP
      AD   (  4) =    0.2908996_DP
      AQ   (  4) =    0.3530008_DP

      FN1(  4,1) =  1.6315720_DP ; FN2(  4,1) = 2.6729620_DP ; FN3(  4,1) = 1.7916860_DP
      FN1(  4,2) = -2.1109590_DP ; FN2(  4,2) = 1.9685940_DP ; FN3(  4,2) = 1.7558710_DP

      ! . Carbon.
      SEPAR(6)   = .TRUE.
      USS  (  6) =  -47.2703200_DP
      UPP  (  6) =  -36.2669180_DP
      BETAS(  6) =  -11.9100150_DP
      BETAP(  6) =   -9.8027550_DP
      ZS   (  6) =    1.5650850_DP
      ZP   (  6) =    1.8423450_DP
      ALP  (  6) =    2.7078070_DP
      EISOL(  6) = -111.2299170_DP
      GSS  (  6) =   11.2007080_DP
      GSP  (  6) =   10.2650270_DP
      GPP  (  6) =   10.7962920_DP
      GP2  (  6) =    9.0425660_DP
      HSP  (  6) =    2.2909800_DP
      DD   (  6) =    0.8332396_DP
      QQ   (  6) =    0.6647750_DP
      AM   (  6) =    0.4116394_DP
      AD   (  6) =    0.5885862_DP
      AQ   (  6) =    0.7647667_DP

      FN1(  6,1) = 0.0501070_DP ; FN2(  6,1) = 6.0031650_DP ; FN3(  6,1) = 1.6422140_DP
      FN1(  6,2) = 0.0507330_DP ; FN2(  6,2) = 6.0029790_DP ; FN3(  6,2) = 0.8924880_DP

      ! . Nitrogen.
      SEPAR(7)  = .TRUE.
      USS  (  7) =  -49.3356720_DP
      UPP  (  7) =  -47.5097360_DP
      BETAS(  7) =  -14.0625210_DP
      BETAP(  7) =  -20.0438480_DP
      ZS   (  7) =    2.0280940_DP
      ZP   (  7) =    2.3137280_DP
      ALP  (  7) =    2.8305450_DP
      EISOL(  7) = -157.6137755_DP
      GSS  (  7) =   11.9047870_DP
      GSP  (  7) =    7.3485650_DP
      GPP  (  7) =   11.7546720_DP
      GP2  (  7) =   10.8072770_DP
      HSP  (  7) =    1.1367130_DP
      DD   (  7) =    0.6577006_DP
      QQ   (  7) =    0.5293383_DP
      AM   (  7) =    0.4375151_DP
      AD   (  7) =    0.5030995_DP
      AQ   (  7) =    0.7364933_DP

      FN1(  7,1) =  1.5016740_DP ; FN2(  7,1) = 5.9011480_DP ; FN3(  7,1) = 1.7107400_DP
      FN1(  7,2) = -1.5057720_DP ; FN2(  7,2) = 6.0046580_DP ; FN3(  7,2) = 1.7161490_DP

      ! . Oxygen.
      SEPAR(8)   = .TRUE.
      USS  (  8) =  -86.9930020_DP
      UPP  (  8) =  -71.8795800_DP
      BETAS(  8) =  -45.2026510_DP
      BETAP(  8) =  -24.7525150_DP
      ZS   (  8) =    3.7965440_DP
      ZP   (  8) =    2.3894020_DP
      ALP  (  8) =    3.2171020_DP
      EISOL(  8) = -289.3422065_DP
      GSS  (  8) =   15.7557600_DP
      GSP  (  8) =   10.6211600_DP
      GPP  (  8) =   13.6540160_DP
      GP2  (  8) =   12.4060950_DP
      HSP  (  8) =    0.5938830_DP
      DD   (  8) =    0.4086173_DP
      QQ   (  8) =    0.5125738_DP
      AM   (  8) =    0.5790430_DP
      AD   (  8) =    0.5299630_DP
      AQ   (  8) =    0.8179630_DP

      FN1(  8,1) = -1.1311280_DP ; FN2(  8,1) = 6.0024770_DP ; FN3(  8,1) = 1.6073110_DP
      FN1(  8,2) =  1.1378910_DP ; FN2(  8,2) = 5.9505120_DP ; FN3(  8,2) = 1.5983950_DP

      ! . Fluorine.
      SEPAR(9)   = .TRUE.
      USS  (  9) = -110.4353030_DP
      UPP  (  9) = -105.6850470_DP
      BETAS(  9) =  -48.4059390_DP
      BETAP(  9) =  -27.7446600_DP
      ZS   (  9) =    4.7085550_DP
      ZP   (  9) =    2.4911780_DP
      ALP  (  9) =    3.3589210_DP
      EISOL(  9) = -437.5171690_DP
      GSS  (  9) =   10.4966670_DP
      GSP  (  9) =   16.0736890_DP
      GPP  (  9) =   14.8172560_DP
      GP2  (  9) =   14.4183930_DP
      HSP  (  9) =    0.7277630_DP
      DD   (  9) =    0.3125302_DP
      QQ   (  9) =    0.4916328_DP
      AM   (  9) =    0.3857650_DP
      AD   (  9) =    0.6768503_DP
      AQ   (  9) =    0.6120047_DP

      FN1(  9,1) = -0.0121660_DP ; FN2(  9,1) = 6.0235740_DP ; FN3(  9,1) = 1.8568590_DP
      FN1(  9,2) = -0.0028520_DP ; FN2(  9,2) = 6.0037170_DP ; FN3(  9,2) = 2.6361580_DP

      ! . Sodium.
      SEPAR(11) = .TRUE.
      ALP  (11) = 1.681_DP
      EISOL(11) = 0.000_DP
      AM   (11) = 0.500_DP

      ! . Magnesium.
      SEPAR(12) = .TRUE.
      USS  ( 12) =  -14.6236880_DP
      UPP  ( 12) =  -14.1734600_DP
      BETAS( 12) =   -2.0716910_DP
      BETAP( 12) =   -0.5695810_DP
      ZS   ( 12) =    0.6985520_DP
      ZP   ( 12) =    1.4834530_DP
      ALP  ( 12) =    1.3291470_DP
      EISOL( 12) =  -22.5530760_DP
      GSS  ( 12) =    6.6943000_DP
      GSP  ( 12) =    6.7939950_DP
      GPP  ( 12) =    6.9104460_DP
      GP2  ( 12) =    7.0908230_DP
      HSP  ( 12) =    0.5433000_DP
      DD   ( 12) =    1.1403950_DP
      QQ   ( 12) =    1.1279899_DP
      AM   ( 12) =    0.2460235_DP
      AD   ( 12) =    0.2695751_DP
      AQ   ( 12) =    0.2767522_DP

      FN1( 12,1) =  2.1170500_DP ; FN2( 12,1) = 6.0094770_DP ; FN3( 12,1) = 2.0844060_DP
      FN1( 12,2) = -2.5477670_DP ; FN2( 12,2) = 4.3953700_DP ; FN3( 12,2) = 2.0636740_DP

      ! . Aluminium.
      SEPAR(13)  = .TRUE.
      USS  ( 13) =  -24.8454040_DP
      UPP  ( 13) =  -22.2641590_DP
      BETAS( 13) =   -0.5943010_DP
      BETAP( 13) =   -0.9565500_DP
      ZS   ( 13) =    1.7028880_DP
      ZP   ( 13) =    1.0736290_DP
      ZD   ( 13) =    1.0000000_DP
      ALP  ( 13) =    1.5217030_DP
      EISOL( 13) =  -46.8647630_DP
      GSS  ( 13) =    5.7767370_DP
      GSP  ( 13) =   11.6598560_DP
      GPP  ( 13) =    6.3477900_DP
      GP2  ( 13) =    6.1210770_DP
      HSP  ( 13) =    4.0062450_DP
      DD   ( 13) =    1.2102799_DP
      QQ   ( 13) =    1.5585645_DP
      AM   ( 13) =    0.2123020_DP
      AD   ( 13) =    0.6418584_DP
      AQ   ( 13) =    0.2262838_DP

      FN1( 13,1) = -0.4730900_DP ; FN2( 13,1) = 1.9158250_DP ; FN3( 13,1) = 1.4517280_DP
      FN1( 13,2) = -0.1540510_DP ; FN2( 13,2) = 6.0050860_DP ; FN3( 13,2) = 2.5199970_DP

      ! . Silicon.
      SEPAR(14)  = .TRUE.
      USS  ( 14) =  -26.7634830_DP
      UPP  ( 14) =  -22.8136350_DP
      BETAS( 14) =   -2.8621450_DP
      BETAP( 14) =   -3.9331480_DP
      ZS   ( 14) =    1.6350750_DP
      ZP   ( 14) =    1.3130880_DP
      ZD   ( 14) =    1.0000000_DP
      ALP  ( 14) =    2.1358090_DP
      EISOL( 14) =  -67.7882140_DP
      GSS  ( 14) =    5.0471960_DP
      GSP  ( 14) =    5.9490570_DP
      GPP  ( 14) =    6.7593670_DP
      GP2  ( 14) =    5.1612970_DP
      HSP  ( 14) =    0.9198320_DP
      DD   ( 14) =    1.3144550_DP
      QQ   ( 14) =    1.2743396_DP
      AM   ( 14) =    0.1854905_DP
      AD   ( 14) =    0.3060715_DP
      AQ   ( 14) =    0.4877432_DP

      FN1( 14,1) = -0.3906000_DP ; FN2( 14,1) = 6.0000540_DP ; FN3( 14,1) = 0.6322620_DP
      FN1( 14,2) =  0.0572590_DP ; FN2( 14,2) = 6.0071830_DP ; FN3( 14,2) = 2.0199870_DP

      ! . Phosphorus.
      SEPAR(15)  = .TRUE.
      USS  ( 15) =  -40.4130960_DP
      UPP  ( 15) =  -29.5930520_DP
      BETAS( 15) =  -12.6158790_DP
      BETAP( 15) =   -4.1600400_DP
      ZS   ( 15) =    2.0175630_DP
      ZP   ( 15) =    1.5047320_DP
      ZD   ( 15) =    1.0000000_DP
      ALP  ( 15) =    1.9405340_DP
      EISOL( 15) = -117.9591740_DP
      GSS  ( 15) =    7.8016150_DP
      GSP  ( 15) =    5.1869490_DP
      GPP  ( 15) =    6.6184780_DP
      GP2  ( 15) =    6.0620020_DP
      HSP  ( 15) =    1.5428090_DP
      DD   ( 15) =    1.0644947_DP
      QQ   ( 15) =    1.1120386_DP
      AM   ( 15) =    0.2867187_DP
      AD   ( 15) =    0.4309446_DP
      AQ   ( 15) =    0.3732517_DP

      FN1( 15,1) = -0.6114210_DP ; FN2( 15,1) = 1.9972720_DP ; FN3( 15,1) = 0.7946240_DP
      FN1( 15,2) = -0.0939350_DP ; FN2( 15,2) = 1.9983600_DP ; FN3( 15,2) = 1.9106770_DP

      ! . Sulphur.
      SEPAR(16)  = .TRUE.
      USS  ( 16) =  -49.8953710_DP
      UPP  ( 16) =  -44.3925830_DP
      BETAS( 16) =   -8.8274650_DP
      BETAP( 16) =   -8.0914150_DP
      ZS   ( 16) =    1.8911850_DP
      ZP   ( 16) =    1.6589720_DP
      ZD   ( 16) =    1.0000000_DP
      ALP  ( 16) =    2.2697060_DP
      EISOL( 16) = -183.4537395_DP
      GSS  ( 16) =    8.9646670_DP
      GSP  ( 16) =    6.7859360_DP
      GPP  ( 16) =    9.9681640_DP
      GP2  ( 16) =    7.9702470_DP
      HSP  ( 16) =    4.0418360_DP
      DD   ( 16) =    1.1214313_DP
      QQ   ( 16) =    1.0086488_DP
      AM   ( 16) =    0.3294622_DP
      AD   ( 16) =    0.6679118_DP
      AQ   ( 16) =    0.6137472_DP

      FN1( 16,1) = -0.3991910_DP ; FN2( 16,1) = 6.0006690_DP ; FN3( 16,1) = 0.9621230_DP
      FN1( 16,2) = -0.0548990_DP ; FN2( 16,2) = 6.0018450_DP ; FN3( 16,2) = 1.5799440_DP

      ! . Chlorine.
      SEPAR(17)  = .TRUE.
      USS  ( 17) = -100.6267470_DP
      UPP  ( 17) =  -53.6143960_DP
      BETAS( 17) =  -27.5285600_DP
      BETAP( 17) =  -11.5939220_DP
      ZS   ( 17) =    2.2462100_DP
      ZP   ( 17) =    2.1510100_DP
      ZD   ( 17) =    1.0000000_DP
      ALP  ( 17) =    2.5172960_DP
      EISOL( 17) = -315.1949480_DP
      GSS  ( 17) =   16.0136010_DP
      GSP  ( 17) =    8.0481150_DP
      GPP  ( 17) =    7.5222150_DP
      GP2  ( 17) =    7.5041540_DP
      HSP  ( 17) =    3.4811530_DP
      DD   ( 17) =    0.9175856_DP
      QQ   ( 17) =    0.7779230_DP
      AM   ( 17) =    0.5885190_DP
      AD   ( 17) =    0.6814522_DP
      AQ   ( 17) =    0.3643694_DP

      FN1( 17,1) = -0.1715910_DP ; FN2( 17,1) = 6.0008020_DP ; FN3( 17,1) = 1.0875020_DP
      FN1( 17,2) = -0.0134580_DP ; FN2( 17,2) = 1.9666180_DP ; FN3( 17,2) = 2.2928910_DP

      ! . Potassium.
      SEPAR(19) = .TRUE.
      ALP  (19) = 1.400_DP
      EISOL(19) = 0.000_DP
      AM   (19) = 0.500_DP

!      ! . Zinc.
!      SEPAR(30)  = .TRUE.
!      USS  ( 30) =  -18.5321980_DP
!      UPP  ( 30) =  -11.0474090_DP
!      BETAS( 30) =   -0.7155780_DP
!      BETAP( 30) =   -6.3518640_DP
!      ZS   ( 30) =    1.8199890_DP
!      ZP   ( 30) =    1.5069220_DP
!      ZD   ( 30) =    1.0000000_DP
!      ALP  ( 30) =    1.3501260_DP
!      EISOL( 30) =  -27.3872000_DP
!      GSS  ( 30) =    9.6771960_DP
!      GSP  ( 30) =    7.7362040_DP
!      GPP  ( 30) =    4.9801740_DP
!      GP2  ( 30) =    4.6696560_DP
!      HSP  ( 30) =    0.6004130_DP
!      DD   ( 30) =    1.5005758_DP
!      QQ   ( 30) =    1.4077174_DP
!      AM   ( 30) =    0.3556485_DP
!      AD   ( 30) =    0.2375689_DP
!      AQ   ( 30) =    0.2661069_DP
!
!      FN1( 30,1) = -0.1112340_DP ; FN2( 30,1) = 6.0014780_DP ; FN3( 30,1) = 1.5160320_DP
!      FN1( 30,2) = -0.1323700_DP ; FN2( 30,2) = 1.9958390_DP ; FN3( 30,2) = 2.5196420_DP

      ! . Zinc (J. Comput. Chem., 25: 1677-1692, 2004, Merz)
      SEPAR(30)  = .TRUE.
      USS  ( 30) =  -16.9746360_DP
      UPP  ( 30) =   -9.7941560_DP
      BETAS( 30) =   -0.5705140_DP
      BETAP( 30) =   -4.1245870_DP
      ZS   ( 30) =    1.4341890_DP
      ZP   ( 30) =    1.4545820_DP
      ZD   ( 30) =    1.0000000_DP
      ALP  ( 30) =    1.3602520_DP
      EISOL( 30) =  -25.4527010_DP
      GSS  ( 30) =    8.4965710_DP
      GSP  ( 30) =    8.3019450_DP
      GPP  ( 30) =    6.4853290_DP
      GP2  ( 30) =    6.2358930_DP
      HSP  ( 30) =    0.8023580_DP
      DD   ( 30) =    1.7983380_DP
      QQ   ( 30) =    1.4583710_DP
      AM   ( 30) =    0.3122590_DP
      AD   ( 30) =    0.2411480_DP
      AQ   ( 30) =    0.2438070_DP

      FN1( 30,1) = -0.2621700_DP ; FN2( 30,1) = 4.7309390_DP ; FN3( 30,1) = 1.8022450_DP
      FN1( 30,2) = -0.1329170_DP ; FN2( 30,2) = 0.9599290_DP ; FN3( 30,2) = 2.3824630_DP

      ! . Gallium.
      SEPAR(31)  = .TRUE.
      USS  ( 31) =  -29.8555930_DP
      UPP  ( 31) =  -21.8753710_DP
      BETAS( 31) =   -4.9456180_DP
      BETAP( 31) =   -0.4070530_DP
      ZS   ( 31) =    1.8470400_DP
      ZP   ( 31) =    0.8394110_DP
      ALP  ( 31) =    1.6051150_DP
      EISOL( 31) =  -57.3280250_DP
      GSS  ( 31) =    8.4585540_DP
      GSP  ( 31) =    8.9256190_DP
      GPP  ( 31) =    5.0868550_DP
      GP2  ( 31) =    4.9830450_DP
      HSP  ( 31) =    2.0512600_DP
      DD   ( 31) =    0.9776692_DP
      QQ   ( 31) =    2.5271534_DP
      AM   ( 31) =    0.3108620_DP
      AD   ( 31) =    0.5129360_DP
      AQ   ( 31) =    0.1546208_DP

      FN1( 31,1) = -0.5601790_DP ; FN2( 31,1) = 5.6232730_DP ; FN3( 31,1) = 1.5317800_DP
      FN1( 31,2) = -0.2727310_DP ; FN2( 31,2) = 1.9918430_DP ; FN3( 31,2) = 2.1838640_DP

      ! . Germanium.
      SEPAR(32) = .TRUE.
      USS  ( 32) =  -35.4671955_DP
      UPP  ( 32) =  -31.5863583_DP
      BETAS( 32) =   -5.3250024_DP
      BETAP( 32) =   -2.2501567_DP
      ZS   ( 32) =    2.2373526_DP
      ZP   ( 32) =    1.5924319_DP
      ALP  ( 32) =    1.9723370_DP
      EISOL( 32) =  -84.0156006_DP
      GSS  ( 32) =    5.3769635_DP
      GSP  ( 32) =   10.2095293_DP
      GPP  ( 32) =    7.6718647_DP
      GP2  ( 32) =    6.9242663_DP
      HSP  ( 32) =    1.3370204_DP
      DD   ( 32) =    1.1920304_DP
      QQ   ( 32) =    1.3321263_DP
      AM   ( 32) =    0.1976098_DP
      AD   ( 32) =    0.3798182_DP
      AQ   ( 32) =    0.3620669_DP

      FN1( 32,1) =  0.9631726_DP ; FN2( 32,1) = 6.0120134_DP ; FN3( 32,1) = 2.1633655_DP
      FN1( 32,2) = -0.9593891_DP ; FN2( 32,2) = 5.7491802_DP ; FN3( 32,2) = 2.1693724_DP

      ! . Arsenic.
      SEPAR(33)  = .TRUE.
      USS  ( 33) =  -38.5074240_DP
      UPP  ( 33) =  -35.1524150_DP
      BETAS( 33) =   -8.2321650_DP
      BETAP( 33) =   -5.0173860_DP
      ZS   ( 33) =    2.6361770_DP
      ZP   ( 33) =    1.7038890_DP
      ALP  ( 33) =    1.7944770_DP
      EISOL( 33) = -122.6326140_DP
      GSS  ( 33) =    8.7890010_DP
      GSP  ( 33) =    5.3979830_DP
      GPP  ( 33) =    8.2872500_DP
      GP2  ( 33) =    8.2103460_DP
      HSP  ( 33) =    1.9510340_DP
      DD   ( 33) =    0.9679655_DP
      QQ   ( 33) =    1.2449874_DP
      AM   ( 33) =    0.3230063_DP
      AD   ( 33) =    0.5042239_DP
      AQ   ( 33) =    0.2574219_DP

      FN1( 33,1) = -0.4600950_DP ; FN2( 33,1) = 1.9831150_DP ; FN3( 33,1) = 1.0867930_DP
      FN1( 33,2) = -0.0889960_DP ; FN2( 33,2) = 1.9929440_DP ; FN3( 33,2) = 2.1400580_DP

      ! . Selenium.
      SEPAR(34)  = .TRUE.
      USS  ( 34) =  -55.3781350_DP
      UPP  ( 34) =  -49.8230760_DP
      BETAS( 34) =   -6.1578220_DP
      BETAP( 34) =   -5.4930390_DP
      ZS   ( 34) =    2.8280510_DP
      ZP   ( 34) =    1.7325360_DP
      ALP  ( 34) =    3.0439570_DP
      EISOL( 34) = -192.7748115_DP
      GSS  ( 34) =    7.4325910_DP
      GSP  ( 34) =   10.0604610_DP
      GPP  ( 34) =    9.5683260_DP
      GP2  ( 34) =    7.7242890_DP
      HSP  ( 34) =    4.0165580_DP
      DD   ( 34) =    0.8719813_DP
      QQ   ( 34) =    1.2244019_DP
      AM   ( 34) =    0.2731566_DP
      AD   ( 34) =    0.7509697_DP
      AQ   ( 34) =    0.5283737_DP

      FN1( 34,1) = 0.0478730_DP ; FN2( 34,1) = 6.0074000_DP ; FN3( 34,1) = 2.0817170_DP
      FN1( 34,2) = 0.1147200_DP ; FN2( 34,2) = 6.0086720_DP ; FN3( 34,2) = 1.5164230_DP

      ! . Bromine.
      SEPAR(35)  = .TRUE.
      USS  ( 35) = -116.6193110_DP
      UPP  ( 35) =  -74.2271290_DP
      BETAS( 35) =  -31.1713420_DP
      BETAP( 35) =   -6.8140130_DP
      ZS   ( 35) =    5.3484570_DP
      ZP   ( 35) =    2.1275900_DP
      ZD   ( 35) =    1.0000000_DP
      ALP  ( 35) =    2.5118420_DP
      EISOL( 35) = -352.5398970_DP
      GSS  ( 35) =   15.9434250_DP
      GSP  ( 35) =   16.0616800_DP
      GPP  ( 35) =    8.2827630_DP
      GP2  ( 35) =    7.8168490_DP
      HSP  ( 35) =    0.5788690_DP
      DD   ( 35) =    0.2759025_DP
      QQ   ( 35) =    0.9970532_DP
      AM   ( 35) =    0.5859399_DP
      AD   ( 35) =    0.6755383_DP
      AQ   ( 35) =    0.3823719_DP

      FN1( 35,1) =  0.9604580_DP ; FN2( 35,1) = 5.9765080_DP ; FN3( 35,1) = 2.3216540_DP
      FN1( 35,2) = -0.9549160_DP ; FN2( 35,2) = 5.9447030_DP ; FN3( 35,2) = 2.3281420_DP

      ! . Cadmium.
      SEPAR(48)  = .TRUE.
      USS  ( 48) =  -15.8285840_DP
      UPP  ( 48) =    8.7497950_DP
      BETAS( 48) =   -8.5819440_DP
      BETAP( 48) =   -0.6010340_DP
      ZS   ( 48) =    1.6793510_DP
      ZP   ( 48) =    2.0664120_DP
      ALP  ( 48) =    1.5253820_DP
      EISOL( 48) =  -22.4502080_DP
      GSS  ( 48) =    9.2069600_DP
      GSP  ( 48) =    8.2315390_DP
      GPP  ( 48) =    4.9481040_DP
      GP2  ( 48) =    4.6696560_DP
      HSP  ( 48) =    1.6562340_DP
      DD   ( 48) =    1.5982681_DP
      QQ   ( 48) =    1.2432402_DP
      AM   ( 48) =    0.3383668_DP
      AD   ( 48) =    0.3570290_DP
      AQ   ( 48) =    0.2820582_DP

      ! . Indium.
      SEPAR(49)  = .TRUE.
      USS  ( 49) =  -26.1762050_DP
      UPP  ( 49) =  -20.0058220_DP
      BETAS( 49) =   -2.9933190_DP
      BETAP( 49) =   -1.8289080_DP
      ZS   ( 49) =    2.0161160_DP
      ZP   ( 49) =    1.4453500_DP
      ALP  ( 49) =    1.4183850_DP
      EISOL( 49) =  -51.9750470_DP
      GSS  ( 49) =    6.5549000_DP
      GSP  ( 49) =    8.2298730_DP
      GPP  ( 49) =    6.2992690_DP
      GP2  ( 49) =    4.9842110_DP
      HSP  ( 49) =    2.6314610_DP
      DD   ( 49) =    1.5766241_DP
      QQ   ( 49) =    1.7774563_DP
      AM   ( 49) =    0.2409004_DP
      AD   ( 49) =    0.4532655_DP
      AQ   ( 49) =    0.3689812_DP

      FN1( 49,1) = -0.3431380_DP ; FN2( 49,1) = 1.9940340_DP ; FN3( 49,1) = 1.6255160_DP
      FN1( 49,2) = -0.1095320_DP ; FN2( 49,2) = 5.6832170_DP ; FN3( 49,2) = 2.8670090_DP

      ! . Tin.
      SEPAR(50) = .TRUE.
      USS  ( 50) =  -34.5501920_DP
      UPP  ( 50) =  -25.8944190_DP
      BETAS( 50) =   -2.7858020_DP
      BETAP( 50) =   -2.0059990_DP
      ZS   ( 50) =    2.3733280_DP
      ZP   ( 50) =    1.6382330_DP
      ALP  ( 50) =    1.6996500_DP
      EISOL( 50) =  -78.8877790_DP
      GSS  ( 50) =   10.1900330_DP
      GSP  ( 50) =    7.2353270_DP
      GPP  ( 50) =    5.6738100_DP
      GP2  ( 50) =    5.1822140_DP
      HSP  ( 50) =    1.0331570_DP
      DD   ( 50) =    1.3120038_DP
      QQ   ( 50) =    1.5681814_DP
      AM   ( 50) =    0.3744959_DP
      AD   ( 50) =    0.3218163_DP
      AQ   ( 50) =    0.2832529_DP

      FN1( 50,1) = -0.1503530_DP ; FN2( 50,1) = 6.0056940_DP ; FN3( 50,1) = 1.7046420_DP
      FN1( 50,2) = -0.0444170_DP ; FN2( 50,2) = 2.2573810_DP ; FN3( 50,2) = 2.4698690_DP

      ! . Antimony.
      SEPAR(51) = .TRUE.
      USS  ( 51) =  -56.4321960_DP
      UPP  ( 51) =  -29.4349540_DP
      BETAS( 51) =  -14.7942170_DP
      BETAP( 51) =   -2.8179480_DP
      ZS   ( 51) =    2.3430390_DP
      ZP   ( 51) =    1.8999920_DP
      ALP  ( 51) =    2.0343010_DP
      EISOL( 51) = -148.9382890_DP
      GSS  ( 51) =    9.2382770_DP
      GSP  ( 51) =    5.2776800_DP
      GPP  ( 51) =    6.3500000_DP
      GP2  ( 51) =    6.2500000_DP
      HSP  ( 51) =    2.4244640_DP
      DD   ( 51) =    1.4091903_DP
      QQ   ( 51) =    1.3521354_DP
      AM   ( 51) =    0.3395177_DP
      AD   ( 51) =    0.4589010_DP
      AQ   ( 51) =    0.2423472_DP

      FN1( 51,1) =  3.0020280_DP ; FN2( 51,1) = 6.0053420_DP ; FN3( 51,1) = 0.8530600_DP
      FN1( 51,2) = -0.0188920_DP ; FN2( 51,2) = 6.0114780_DP ; FN3( 51,2) = 2.7933110_DP

      ! . Tellurium.
      SEPAR(52) = .TRUE.
      USS  ( 52) =  -44.9380360_DP
      UPP  ( 52) =  -46.3140990_DP
      BETAS( 52) =   -2.6651460_DP
      BETAP( 52) =   -3.8954300_DP
      ZS   ( 52) =    4.1654920_DP
      ZP   ( 52) =    1.6475550_DP
      ALP  ( 52) =    2.4850190_DP
      EISOL( 52) = -168.0945925_DP
      GSS  ( 52) =   10.2550730_DP
      GSP  ( 52) =    8.1691450_DP
      GPP  ( 52) =    7.7775920_DP
      GP2  ( 52) =    7.7551210_DP
      HSP  ( 52) =    3.7724620_DP
      DD   ( 52) =    0.3484177_DP
      QQ   ( 52) =    1.5593085_DP
      AM   ( 52) =    0.3768862_DP
      AD   ( 52) =    1.1960743_DP
      AQ   ( 52) =    0.2184786_DP

      FN1( 52,1) =  0.0333910_DP ; FN2( 52,1) = 5.9563790_DP ; FN3( 52,1) = 2.2775750_DP 
      FN1( 52,2) = -1.9218670_DP ; FN2( 52,2) = 4.9732190_DP ; FN3( 52,2) = 0.5242430_DP

      ! . Iodine.
      SEPAR(53)  = .TRUE.
      USS  ( 53) =  -96.4540370_DP
      UPP  ( 53) =  -61.0915820_DP
      BETAS( 53) =  -14.4942340_DP
      BETAP( 53) =   -5.8947030_DP
      ZS   ( 53) =    7.0010130_DP
      ZP   ( 53) =    2.4543540_DP
      ZD   ( 53) =    1.0000000_DP
      ALP  ( 53) =    1.9901850_DP
      EISOL( 53) = -288.3160860_DP
      GSS  ( 53) =   13.6319430_DP
      GSP  ( 53) =   14.9904060_DP
      GPP  ( 53) =    7.2883300_DP
      GP2  ( 53) =    5.9664070_DP
      HSP  ( 53) =    2.6300350_DP
      DD   ( 53) =    0.1581469_DP
      QQ   ( 53) =    1.0467302_DP
      AM   ( 53) =    0.5009902_DP
      AD   ( 53) =    1.6699104_DP
      AQ   ( 53) =    0.5153082_DP

      FN1( 53,1) = -0.1314810_DP ; FN2( 53,1) = 5.2064170_DP ; FN3( 53,1) = 1.7488240_DP
      FN1( 53,2) = -0.0368970_DP ; FN2( 53,2) = 6.0101170_DP ; FN3( 53,2) = 2.7103730_DP

      ! . Mercury.
      SEPAR(80)  = .TRUE.
      USS  ( 80) =  -17.7622290_DP
      UPP  ( 80) =  -18.3307510_DP
      BETAS( 80) =   -3.1013650_DP
      BETAP( 80) =   -3.4640310_DP
      ZS   ( 80) =    1.4768850_DP
      ZP   ( 80) =    2.4799510_DP
      ALP  ( 80) =    1.5293770_DP
      EISOL( 80) =  -28.8997380_DP
      GSS  ( 80) =    6.6247200_DP
      GSP  ( 80) =   10.6392970_DP
      GPP  ( 80) =   14.7092830_DP
      GP2  ( 80) =   16.0007400_DP
      HSP  ( 80) =    2.0363110_DP
      DD   ( 80) =    1.2317811_DP
      QQ   ( 80) =    1.2164033_DP
      AM   ( 80) =    0.2434664_DP
      AD   ( 80) =    0.4515472_DP
      AQ   ( 80) =    0.2618394_DP

      FN1( 80,1) =  1.0827200_DP ; FN2( 80,1) = 6.4965980_DP ; FN3( 80,1) = 1.1951460_DP
      FN1( 80,2) = -0.0965530_DP ; FN2( 80,2) = 3.9262810_DP ; FN3( 80,2) = 2.6271600_DP

      ! . Thallium.
      SEPAR(81)  = .TRUE.
      USS  ( 81) =  -30.0531700_DP
      UPP  ( 81) =  -26.9206370_DP
      BETAS( 81) =   -1.0844950_DP
      BETAP( 81) =   -7.9467990_DP
      ZS   ( 81) =    6.8679210_DP
      ZP   ( 81) =    1.9694450_DP
      ALP  ( 81) =    1.3409510_DP
      EISOL( 81) =  -56.6492050_DP
      GSS  ( 81) =   10.4604120_DP
      GSP  ( 81) =   11.2238830_DP
      GPP  ( 81) =    4.9927850_DP
      GP2  ( 81) =    8.9627270_DP
      HSP  ( 81) =    2.5304060_DP
      DD   ( 81) =    0.0781362_DP
      QQ   ( 81) =    1.5317110_DP
      AM   ( 81) =    0.3844326_DP
      AD   ( 81) =    2.5741815_DP
      AQ   ( 81) =    0.2213264_DP

      FN1( 81,1) = -1.3613990_DP ; FN2( 81,1) = 3.5572260_DP ; FN3( 81,1) = 1.0928020_DP
      FN1( 81,2) = -0.0454010_DP ; FN2( 81,2) = 2.3069950_DP ; FN3( 81,2) = 2.9650290_DP

      ! . Lead.
      SEPAR(82)  = .TRUE.
      USS  ( 82) =  -30.3227560_DP
      UPP  ( 82) =  -24.4258340_DP
      BETAS( 82) =   -6.1260240_DP
      BETAP( 82) =   -1.3954300_DP
      ZS   ( 82) =    3.1412890_DP
      ZP   ( 82) =    1.8924180_DP
      ALP  ( 82) =    1.6200450_DP
      EISOL( 82) =  -73.4660775_DP
      GSS  ( 82) =    7.0119920_DP
      GSP  ( 82) =    6.7937820_DP
      GPP  ( 82) =    5.1837800_DP
      GP2  ( 82) =    5.0456510_DP
      HSP  ( 82) =    1.5663020_DP
      DD   ( 82) =    0.9866290_DP
      QQ   ( 82) =    1.5940562_DP
      AM   ( 82) =    0.2576991_DP
      AD   ( 82) =    0.4527678_DP
      AQ   ( 82) =    0.2150175_DP

      FN1( 82,1) = -0.1225760_DP ; FN2( 82,1) = 6.0030620_DP ; FN3( 82,1) = 1.9015970_DP
      FN1( 82,2) = -0.0566480_DP ; FN2( 82,2) = 4.7437050_DP ; FN3( 82,2) = 2.8618790_DP

      ! . Bismuth.
      SEPAR( 83) = .TRUE.
      USS  ( 83) =  -33.4959380_DP
      UPP  ( 83) =  -35.5210260_DP
      BETAS( 83) =   -5.6072830_DP
      BETAP( 83) =   -5.8001520_DP
      ZS   ( 83) =    4.9164510_DP
      ZP   ( 83) =    1.9349350_DP
      ALP  ( 83) =    1.8574310_DP
      EISOL( 83) = -109.2774910_DP
      GSS  ( 83) =    4.9894800_DP
      GSP  ( 83) =    6.1033080_DP
      GPP  ( 83) =    8.6960070_DP
      GP2  ( 83) =    8.3354470_DP
      HSP  ( 83) =    0.5991220_DP
      DD   ( 83) =    0.2798609_DP
      QQ   ( 83) =    1.5590294_DP
      AM   ( 83) =    0.1833693_DP
      AD   ( 83) =    0.6776013_DP
      AQ   ( 83) =    0.2586520_DP

      FN1( 83,1) = 2.5816930_DP ; FN2( 83,1) = 5.0940220_DP ; FN3( 83,1) = 0.4997870_DP
      FN1( 83,2) = 0.0603200_DP ; FN2( 83,2) = 6.0015380_DP ; FN3( 83,2) = 2.4279700_DP

      END SUBROUTINE PARAMETERS_PM3

   END SUBROUTINE MOPAC_PARAMETERS_INITIALIZE

   !------------------------
   SUBROUTINE PARAMETERS_RM1
   !------------------------

      ! . Hydrogen.
      SEPAR ( 1) = .TRUE.
      ALP   ( 1) =      3.0683595_DP
      EISOL ( 1) =    -11.9606770_DP
      BETAS ( 1) =     -5.7654447_DP
      ZS    ( 1) =      1.0826737_DP
      AM    ( 1) =      0.5138998_DP
      AD    ( 1) =      0.5138998_DP
      AQ    ( 1) =      0.5138998_DP
      USS   ( 1) =    -11.9606770_DP
      GSS   ( 1) =     13.9832130_DP

      FN1( 1, 1) =      0.1028888_DP ; FN2( 1, 1) =      5.9017227_DP ; FN3( 1, 1) =      1.1750118_DP
      FN1( 1, 2) =      0.0645745_DP ; FN2( 1, 2) =      6.4178567_DP ; FN3( 1, 2) =      1.9384448_DP
      FN1( 1, 3) =     -0.0356739_DP ; FN2( 1, 3) =      2.8047313_DP ; FN3( 1, 3) =      1.6365524_DP

      ! . Carbon.
      SEPAR ( 6) = .TRUE.
      ALP   ( 6) =      2.7928208_DP
      EISOL ( 6) =   -117.8673444_DP
      BETAS ( 6) =    -15.4593243_DP
      BETAP ( 6) =     -8.2360864_DP
      ZS    ( 6) =      1.8501880_DP
      ZP    ( 6) =      1.7683009_DP
      DD    ( 6) =      0.7967571_DP
      QQ    ( 6) =      0.6926111_DP
      AM    ( 6) =      0.4797179_DP
      AD    ( 6) =      0.5097516_DP
      AQ    ( 6) =      0.6614863_DP
      USS   ( 6) =    -51.7255603_DP
      UPP   ( 6) =    -39.4072894_DP
      GSS   ( 6) =     13.0531244_DP
      GSP   ( 6) =     11.3347939_DP
      GPP   ( 6) =     10.9511374_DP
      GP2   ( 6) =      9.7239510_DP
      HSP   ( 6) =      1.5521513_DP

      FN1( 6, 1) =      0.0746227_DP ; FN2( 6, 1) =      5.7392160_DP ; FN3( 6, 1) =      1.0439698_DP
      FN1( 6, 2) =      0.0117705_DP ; FN2( 6, 2) =      6.9240173_DP ; FN3( 6, 2) =      1.6615957_DP
      FN1( 6, 3) =      0.0372066_DP ; FN2( 6, 3) =      6.2615894_DP ; FN3( 6, 3) =      1.6315872_DP
      FN1( 6, 4) =     -0.0027066_DP ; FN2( 6, 4) =      9.0000373_DP ; FN3( 6, 4) =      2.7955790_DP

      ! . Nitrogen.
      SEPAR ( 7) = .TRUE.
      ALP   ( 7) =      2.9642254_DP
      EISOL ( 7) =   -205.0876419_DP
      BETAS ( 7) =    -20.8712455_DP
      BETAP ( 7) =    -16.6717185_DP
      ZS    ( 7) =      2.3744716_DP
      ZP    ( 7) =      1.9781257_DP
      DD    ( 7) =      0.6495620_DP
      QQ    ( 7) =      0.6191441_DP
      AM    ( 7) =      0.4809762_DP
      AD    ( 7) =      0.9706271_DP
      AQ    ( 7) =      0.8023683_DP
      USS   ( 7) =    -70.8512372_DP
      UPP   ( 7) =    -57.9773092_DP
      GSS   ( 7) =     13.0873623_DP
      GSP   ( 7) =     13.2122683_DP
      GPP   ( 7) =     13.6992432_DP
      GP2   ( 7) =     11.9410395_DP
      HSP   ( 7) =      5.0000085_DP

      FN1( 7, 1) =      0.0607338_DP ; FN2( 7, 1) =      4.5889295_DP ; FN3( 7, 1) =      1.3787388_DP
      FN1( 7, 2) =      0.0243856_DP ; FN2( 7, 2) =      4.6273052_DP ; FN3( 7, 2) =      2.0837070_DP
      FN1( 7, 3) =     -0.0228343_DP ; FN2( 7, 3) =      2.0527466_DP ; FN3( 7, 3) =      1.8676382_DP

      ! . Oxygen.
      SEPAR ( 8) = .TRUE.
      ALP   ( 8) =      4.1719672_DP
      EISOL ( 8) =   -312.0403540_DP
      BETAS ( 8) =    -29.8510121_DP
      BETAP ( 8) =    -29.1510131_DP
      ZS    ( 8) =      3.1793691_DP
      ZP    ( 8) =      2.5536191_DP
      DD    ( 8) =      0.4886701_DP
      QQ    ( 8) =      0.4796114_DP
      AM    ( 8) =      0.5146059_DP
      AD    ( 8) =      1.0065837_DP
      AQ    ( 8) =      0.8953404_DP
      USS   ( 8) =    -96.9494807_DP
      UPP   ( 8) =    -77.8909298_DP
      GSS   ( 8) =     14.0024279_DP
      GSP   ( 8) =     14.9562504_DP
      GPP   ( 8) =     14.1451514_DP
      GP2   ( 8) =     12.7032550_DP
      HSP   ( 8) =      3.9321716_DP

      FN1( 8, 1) =      0.2309355_DP ; FN2( 8, 1) =      5.2182874_DP ; FN3( 8, 1) =      0.9036355_DP
      FN1( 8, 2) =      0.0585987_DP ; FN2( 8, 2) =      7.4293293_DP ; FN3( 8, 2) =      1.5175461_DP

      ! . Fluorine.
      SEPAR ( 9) = .TRUE.
      ALP   ( 9) =      6.0000006_DP
      EISOL ( 9) =   -484.5957024_DP
      BETAS ( 9) =    -70.0000051_DP
      BETAP ( 9) =    -32.6798271_DP
      ZS    ( 9) =      4.4033791_DP
      ZP    ( 9) =      2.6484156_DP
      DD    ( 9) =      0.3488927_DP
      QQ    ( 9) =      0.4624444_DP
      AM    ( 9) =      0.6145135_DP
      AD    ( 9) =      0.9225494_DP
      AQ    ( 9) =      0.6236567_DP
      USS   ( 9) =   -134.1836959_DP
      UPP   ( 9) =   -107.8466092_DP
      GSS   ( 9) =     16.7209132_DP
      GSP   ( 9) =     16.7614263_DP
      GPP   ( 9) =     15.2258103_DP
      GP2   ( 9) =     14.8657868_DP
      HSP   ( 9) =      1.9976617_DP

      FN1( 9, 1) =      0.4030203_DP ; FN2( 9, 1) =      7.2044196_DP ; FN3( 9, 1) =      0.8165301_DP
      FN1( 9, 2) =      0.0708583_DP ; FN2( 9, 2) =      9.0000156_DP ; FN3( 9, 2) =      1.4380238_DP

      ! . Phosphorus.
      SEPAR (15) = .TRUE.
      ALP   (15) =      1.9099329_DP
      EISOL (15) =   -123.1797789_DP
      BETAS (15) =     -6.1351497_DP
      BETAP (15) =     -5.9444213_DP
      ZS    (15) =      2.1224012_DP
      ZP    (15) =      1.7432795_DP
      DD    (15) =      1.0106955_DP
      QQ    (15) =      0.9598690_DP
      AM    (15) =      0.4072250_DP
      AD    (15) =      0.3932223_DP
      AQ    (15) =      0.3123478_DP
      USS   (15) =    -41.8153318_DP
      UPP   (15) =    -34.3834253_DP
      GSS   (15) =     11.0805926_DP
      GSP   (15) =      5.6833920_DP
      GPP   (15) =      7.6041756_DP
      GP2   (15) =      7.4026518_DP
      HSP   (15) =      1.1618179_DP

      FN1(15, 1) =     -0.4106347_DP ; FN2(15, 1) =      6.0875283_DP ; FN3(15, 1) =      1.3165026_DP
      FN1(15, 2) =     -0.1629929_DP ; FN2(15, 2) =      7.0947260_DP ; FN3(15, 2) =      1.9072132_DP
      FN1(15, 3) =     -0.0488712_DP ; FN2(15, 3) =      8.9997931_DP ; FN3(15, 3) =      2.6585778_DP

      ! . Sulphur.
      SEPAR (16) = .TRUE.
      ALP   (16) =      2.4401564_DP
      EISOL (16) =   -185.3861382_DP
      BETAS (16) =     -1.9591072_DP
      BETAP (16) =     -8.7743065_DP
      ZS    (16) =      2.1334431_DP
      ZP    (16) =      1.8746065_DP
      DD    (16) =      0.9936921_DP
      QQ    (16) =      0.8926247_DP
      AM    (16) =      0.4589594_DP
      AD    (16) =      0.6930934_DP
      AQ    (16) =      0.4959415_DP
      USS   (16) =    -55.1677512_DP
      UPP   (16) =    -46.5293042_DP
      GSS   (16) =     12.4882841_DP
      GSP   (16) =      8.5691057_DP
      GPP   (16) =      8.5230117_DP
      GP2   (16) =      7.6686330_DP
      HSP   (16) =      3.8897893_DP

      FN1(16, 1) =     -0.7460106_DP ; FN2(16, 1) =      4.8103800_DP ; FN3(16, 1) =      0.5938013_DP
      FN1(16, 2) =     -0.0651929_DP ; FN2(16, 2) =      7.2076086_DP ; FN3(16, 2) =      1.2949201_DP
      FN1(16, 3) =     -0.0065598_DP ; FN2(16, 3) =      9.0000018_DP ; FN3(16, 3) =      1.8006015_DP

      ! . Chlorine.
      SEPAR (17) = .TRUE.
      ALP   (17) =      3.6935883_DP
      EISOL (17) =   -382.4700938_DP
      BETAS (17) =    -19.9243043_DP
      BETAP (17) =    -11.5293520_DP
      ZS    (17) =      3.8649107_DP
      ZP    (17) =      1.8959314_DP
      DD    (17) =      0.4541788_DP
      QQ    (17) =      0.8825847_DP
      AM    (17) =      0.5645068_DP
      AD    (17) =      0.7489737_DP
      AQ    (17) =      0.7707320_DP
      USS   (17) =   -118.4730692_DP
      UPP   (17) =    -76.3533034_DP
      GSS   (17) =     15.3602310_DP
      GSP   (17) =     13.3067117_DP
      GPP   (17) =     12.5650264_DP
      GP2   (17) =      9.6639708_DP
      HSP   (17) =      1.7648990_DP

      FN1(17, 1) =      0.1294711_DP ; FN2(17, 1) =      2.9772442_DP ; FN3(17, 1) =      1.4674978_DP
      FN1(17, 2) =      0.0028890_DP ; FN2(17, 2) =      7.0982759_DP ; FN3(17, 2) =      2.5000272_DP

      ! . Bromine.
      SEPAR (35) = .TRUE.
      ALP   (35) =      2.8671053_DP
      EISOL (35) =   -357.1164272_DP
      BETAS (35) =     -1.3413984_DP
      BETAP (35) =     -8.2022599_DP
      ZS    (35) =      5.7315721_DP
      ZP    (35) =      2.0314758_DP
      DD    (35) =      0.2099005_DP
      QQ    (35) =      1.0442262_DP
      AM    (35) =      0.6290199_DP
      AD    (35) =      1.3165755_DP
      AQ    (35) =      0.5863476_DP
      USS   (35) =   -113.4839818_DP
      UPP   (35) =    -76.1872002_DP
      GSS   (35) =     17.1156307_DP
      GSP   (35) =     15.6241925_DP
      GPP   (35) =     10.7354629_DP
      GP2   (35) =      8.8605620_DP
      HSP   (35) =      2.2351276_DP

      FN1(35, 1) =      0.9868994_DP ; FN2(35, 1) =      4.2848419_DP ; FN3(35, 1) =      2.0001970_DP
      FN1(35, 2) =     -0.9273125_DP ; FN2(35, 2) =      4.5400591_DP ; FN3(35, 2) =      2.0161770_DP

      ! . Iodine.
      SEPAR (53) = .TRUE.
      ALP   (53) =      2.1415709_DP
      EISOL (53) =   -248.4933241_DP
      BETAS (53) =     -4.1931615_DP
      BETAP (53) =     -4.4003841_DP
      ZS    (53) =      2.5300375_DP
      ZP    (53) =      2.3173868_DP
      DD    (53) =      1.2963425_DP
      QQ    (53) =      1.1085963_DP
      AM    (53) =      0.7350144_DP
      AD    (53) =      0.3717412_DP
      AQ    (53) =      0.3512697_DP
      USS   (53) =    -74.8999784_DP
      UPP   (53) =    -51.4102380_DP
      GSS   (53) =     19.9997413_DP
      GSP   (53) =      7.6895767_DP
      GPP   (53) =      7.3048834_DP
      GP2   (53) =      6.8542461_DP
      HSP   (53) =      1.4160294_DP

      FN1(53, 1) =     -0.0814772_DP ; FN2(53, 1) =      1.5606507_DP ; FN3(53, 1) =      2.0000206_DP
      FN1(53, 2) =      0.0591499_DP ; FN2(53, 2) =      5.7611127_DP ; FN3(53, 2) =      2.2048880_DP

   END SUBROUTINE PARAMETERS_RM1

END MODULE MOPAC_PARAMETERS
