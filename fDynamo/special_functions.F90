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
!========================================================================
!                       The Special Functions Module
!========================================================================
!
! . Functions:
!
!   ERFC                           Calculate erfc.
!
! . Notes:
!
!   This subroutine was modified from the subroutine ERRINT from the
!   package SPECFN and uses the method from Irene A. Stegun and Ruth
!   Zucker in J. RES. NBS, VOL 74B, NO. 3, 211-224, (1970).
!
!   It is probably not the most efficient ERFC subroutine available but
!   ERFC is obtained to the maximum machine accuracy possible.
!
!   The precision of the subroutine is determined by the parameters
!   NBC and NBM. NBC is the number of binary digits in the characteristic
!   of a floating point number while NBM is the accuracy desired or the
!   maximum number of binary digits in the mantissa of a floating point
!   number.
!
!
! Method.    Power series for abs(x) .le. 1,ulps,upper limit for series
!            continued fraction for abs(x) .gt. 1 and .le. ulcf (upper
!                              limit of continued fraction) 
!
! Range.     abs(x) .le. ulcf           ulcf=.83*(2.**((nbc-1)/2))
!            abs(ulcf) approx.  9.3, nbc=8
!                              26.5, nbc=11
!            Beyond this range the limiting values of the functions are
!            supplied - erf(inf)=1.0,   erfc(inf)=0.
!
! Accuracy.  Routine will return (nbm-i-3) significant binary digits
!            where i is the number of binary digits representing the
!            integer part of x**2. (This is essentially the accuracy of
!            the exponential routine.)
!
! Precision. Variable - by setting desired NBM.
!
!========================================================================
MODULE SPECIAL_FUNCTIONS

! . Module declarations.
USE DEFINITIONS, ONLY : DP

IMPLICIT NONE
PRIVATE
PUBLIC :: ERFC

!========================================================================
CONTAINS
!========================================================================

   !------------------
   FUNCTION ERFC ( X )
   !------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: ERFC

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: X

   ! . Local parameter.
   INTEGER,            PARAMETER :: NBC = 11, NBM = 60
   REAL ( KIND = DP ), PARAMETER :: ZERO = 0.0_DP, ONE = 1.0_DP, TWO = 2.0_DP, FOUR = 4.0_DP, ULPS = 1.0_DP, &
                                    CONS = 0.83_DP, TRRTPI = 1.128379167095512574_DP

   ! . Local scalars.
   REAL ( KIND = DP ) :: AN, BN, C1, DN, F, FN, FNM1, FNM2, GN, GNM1, GNM2, PREV, PWR, RNBC, SCF, SUM, TN, &
                         TOLER, ULCF, WN, Y, YSQ

   RNBC=NBC
   TOLER = TWO ** (-NBM)

   ! . Test on zero.
   IF ( X == ZERO ) THEN

      ERFC = ONE
      RETURN

   END IF

   Y=ABS(X)
   YSQ=Y ** 2

   IF ( ( Y - ULPS ) <= ZERO ) GO TO 3

      ! . Maximum argument.
      C1=TWO ** ((RNBC-ONE)/TWO)
      ULCF=CONS * C1

      ! . Scale factor.
      SCF=TWO ** (C1 ** 2 - RNBC)

   ! . Limiting value.
   IF ( ( Y - ULCF ) <= ZERO ) GO TO 10

      ERFC=ZERO
      GO TO 7

      ! . Power series method.
  3   SUM=ZERO
      DN=ONE
      TN=ONE
      PWR=TWO*YSQ
  6   DN=DN + TWO
      TN=PWR*TN/DN
      SUM=TN+SUM

      ! . Tolerance check.
      IF ( ( TN - TOLER ) >= ZERO ) GO TO 6

      ERFC = ONE - ( SUM + ONE ) * TRRTPI * Y * EXP ( - YSQ ) 

      ! . Negative argument.
  7   CONTINUE
      IF ( X < ZERO ) ERFC = TWO - ERFC
      RETURN

      ! . Continued fraction method.
  10  FNM2=ZERO
      GNM2=ONE
      FNM1=TWO * Y
      GNM1=TWO * YSQ + ONE

      PREV=FNM1/GNM1

      WN=ONE
      BN=GNM1 + FOUR

  14  AN=-WN * (WN + ONE)

      FN=BN * FNM1 + AN * FNM2
      GN=BN * GNM1 + AN * GNM2

      F=FN/GN

      ! . Tolerance check.
      IF ( ABS ( ONE - ( F / PREV ) ) - TOLER <= ZERO ) GO TO 12
      IF ( ( PREV - F ) > ZERO ) GO TO 18

      ! . Both FN and GN must be tested if ABS ( X ) < .61.
      IF ( GN < SCF ) GO TO 16

      ! . Scaling.
      FN=FN/SCF
      GN=GN/SCF
      FNM1=FNM1/SCF 
      GNM1=GNM1/SCF 

  16  FNM2=FNM1
      GNM2=GNM1
      FNM1=FN
      GNM1=GN

      WN=WN + TWO
      BN=BN + FOUR
      PREV=F
      GO TO 14

  18  F=PREV
  12  ERFC=F*EXP(-YSQ)*TRRTPI/TWO

   GO TO 7

   END FUNCTION ERFC

END MODULE SPECIAL_FUNCTIONS
