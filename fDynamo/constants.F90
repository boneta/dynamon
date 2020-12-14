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
!                              The Constants Module
!===============================================================================
!
! . Mathematical constants:
!
!   PI                             Pi.
!
! . Physical constants:
!
!   AMU_TO_KG                      a.m.u. to kg.
!   ANGSTROMS_TO_BOHRS             Angstroms to bohrs.
!   ATM_TO_PASCALS                 Atmospheres to Pascals.
!   AU_TO_DB                       Atomic units to Debyes.
!   AU_TO_EV                       Atomic units to eVs.
!   AU_TO_KJ                       Atomic units to kJ mol^-1.
!   BOHRS_TO_ANGSTROMS             Bohrs to Angstroms.
!   C                              The speed of light.
!   EV_TO_KJ                       eV to kJ mol^-1.
!   E0                             Permittivity of vacuum.
!   KBOLTZ                         Boltzmann's constant.
!   KCAL_TO_KJ                     kcal/mole to kJ/mole.
!   NAVOGADRO                      Avogadro's number.
!   PI                             Pi.
!   PLANCK                         Planck's constant.
!   TO_COULOMBS                    Elementary charge in coulombs.
!   UNDEFINED                      An undefined coordinate.
!
! . Constants derived from those defined above:
!
!   ELECT_CONST                    e^2/Angstroms to kJ mol^-1.
!   PV_TO_KJ_MOLE                  Atm. Ang.^3 to kJ mol^-1.
!   R                              The gas constant.
!   TO_DEGREES                     Radians to degrees.
!   TO_HZ                          cm^-1 to ps^-1.
!   TO_KM_PER_MOLE                 To km/mole from internal units.
!   TO_RADIANS                     Degrees to radians.
!   TO_WAVENUMBERS                 To cm^-1 from internal units.
!
! . Constants for the dynamics:
!
!   AMUA2PS2_TO_K                  From amu (A ps^-1)^2 to Kelvin.
!   AMUA2PS2_TO_KJMOL              From amu (A ps^-1)^2 to kJ mol^-1.
!   MS_TO_APS                      From m s^-1 to A ps^-1.
!   KJMOL_TO_AMUA2PS2              From kJ mol^-1 to amu (A ps^-1)^2.
!
! . Notes:
!
!   S.I. units are used by default, although more chemically common units
!   often override these. For example, lengths are in Angstroms, energies
!   in kJ mol^-1 and frequencies in cm^-1.
!
!   The values of all physical constants were obtained from the "CRC Handbook
!   of Chemistry and Physics", D.R. Lide (Editor in Chief), 73rd Edition,
!   1992-1993, CRC Press.
!
!===============================================================================
MODULE CONSTANTS

! . Module declarations.
USE DEFINITIONS, ONLY : DP

IMPLICIT NONE
PUBLIC

!===============================================================================
! . Mathematical and physical constants.
!===============================================================================
! . Atomic mass units to kg.
REAL ( KIND = DP ), PARAMETER :: AMU_TO_KG = 1.6605402E-27_DP

! . Atmospheres to Pascals.
REAL ( KIND = DP ), PARAMETER :: ATM_TO_PASCALS = 1.013250E+5_DP

! . The conversion factor from atomic units to Debyes.
REAL ( KIND = DP ), PARAMETER :: AU_TO_DB = 2.54176568_DP

! . The conversion factor from atomic units to electron-volts.
REAL ( KIND = DP ), PARAMETER :: AU_TO_EV = 27.21_DP

! . The conversion factor from atomic units to kJ/mole.
! -- Este valor (2625.41736504_DP) discrepa de: 2625.49962955_DP
REAL ( KIND = DP ), PARAMETER :: AU_TO_KJ = 27.21_DP * 4.184_DP * 23.061_DP

! . The conversion factor from Bohrs to Angstroms.
REAL ( KIND = DP ), PARAMETER :: BOHRS_TO_ANGSTROMS = 0.529177249_DP

! . The speed of light.
REAL ( KIND = DP ), PARAMETER :: C = 299792458.0_DP

! . The conversion factor from electron-volts to kJ/mole.
REAL ( KIND = DP ), PARAMETER :: EV_TO_KJ = 4.184_DP * 23.061_DP

! . The permittivity of the vacuum (F m^-1).
REAL ( KIND = DP ), PARAMETER :: E0 = 8.854187817E-12_DP

! . Boltzmann's constant.
REAL ( KIND = DP ), PARAMETER :: KBOLTZ = 1.380658E-23_DP

! . kcal mol^-1 to kJ mol^-1.
REAL ( KIND = DP ), PARAMETER :: KCAL_TO_KJ = 4.184_DP

! . Avogadro's number.
REAL ( KIND = DP ), PARAMETER :: NAVOGADRO = 6.0221367E+23_DP

! . Pi.
REAL ( KIND = DP ), PARAMETER :: PI = 3.141592653589793_DP

! . Planck's constant.
REAL ( KIND = DP ), PARAMETER :: PLANCK = 6.6260755E-34_DP

! . Elementary charge in coulombs.
REAL ( KIND = DP ), PARAMETER :: TO_COULOMBS = 1.60217733E-19_DP 

! . The constant to use for an undefined coordinate.
REAL ( KIND = DP ), PARAMETER :: UNDEFINED = 999999.0_DP

!========================================================================
! . Derived constants.
!========================================================================
! . The conversion factor from Angstroms to Bohrs.
REAL ( KIND = DP ), PARAMETER :: ANGSTROMS_TO_BOHRS = 1.0_DP / BOHRS_TO_ANGSTROMS

! . The conversion factor from e^2/Angstroms to kJ mol^-1.
REAL ( KIND = DP ), PARAMETER :: ELECT_CONST = ( 1.0E+7_DP * NAVOGADRO * TO_COULOMBS**2 ) / ( 4.0_DP * PI * E0 )

! . The conversion factor from atm. Angstroms^3 to kJ mol^-1.
REAL ( KIND = DP ), PARAMETER :: PV_TO_KJ_MOLE = ATM_TO_PASCALS * NAVOGADRO * 1.0E-33_DP

! . The gas constant (kJ mol^-1 K^-1).
REAL ( KIND = DP ), PARAMETER :: R = KBOLTZ * NAVOGADRO * 1.0E-3_DP

! . The conversion factor from cm^-1 to ps^-1.
REAL ( KIND = DP ), PARAMETER :: TO_HZ = C / 1.0E+10_DP

! . Radians to degrees.
REAL ( KIND = DP ), PARAMETER :: TO_DEGREES = 180.0_DP / PI

! . The conversion factor from internal units to km/mole.
REAL ( KIND = DP ), PARAMETER :: TO_KM_PER_MOLE = ( NAVOGADRO * TO_COULOMBS**2 ) / ( 12.0_DP * 1.0E+3_DP * E0 * C**2 * AMU_TO_KG )

! . Degrees to radians.
REAL ( KIND = DP ), PARAMETER :: TO_RADIANS = PI / 180.0_DP

! . The conversion factor from internal units to cm^-1.
REAL ( KIND = DP ), PARAMETER :: TO_WAVENUMBERS = 1.0E+11_DP / ( 2.0_DP * PI * C )

!========================================================================
! . Dynamics constants.
!========================================================================
! . The conversion factor from amu A^2 ps^-2 to kJ mol^-1 (equivalent to
! . AMU_TO_KG * NAVOGADRO * 10^-3 / MS_TO_APS^2 ).
REAL ( KIND = DP ), PARAMETER :: AMUA2PS2_TO_KJMOL = 1.0E-2_DP

! . The conversion factor from m s^-1 to A ps^-1.
REAL ( KIND = DP ), PARAMETER :: MS_TO_APS = 1.0E-2_DP

! . The conversion factor from amu A^2 ps^-2 to Kelvin.
REAL ( KIND = DP ), PARAMETER :: AMUA2PS2_TO_K = AMU_TO_KG / ( KBOLTZ * MS_TO_APS**2 )

! . The conversion factor from kJ mol^-1 to amu A^2 ps^-2.
REAL ( KIND = DP ), PARAMETER :: KJMOL_TO_AMUA2PS2 = 1.0_DP / AMUA2PS2_TO_KJMOL

END MODULE CONSTANTS
