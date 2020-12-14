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
!                             The Statistics Module
!===============================================================================
!
! . Subroutines:
!
!   STATISTICS_ACF_DIRECT          Calculate an autocorrelation function with a
!                                  direct method.
!   STATISTICS_CCF_DIRECT          Calculate a cross correlation function with a
!                                  direct method.
!
! . Notes:
!
!   Both subroutines calculate an autocorrelation using a direct scheme (see, for
!   example, the book by Allen & Tildesley). This is not the most efficient
!   scheme for large systems but it is simple to implement.
!
!===============================================================================
MODULE STATISTICS

! . Module declarations.
USE DEFINITIONS, ONLY : DP

IMPLICIT NONE
PRIVATE
PUBLIC :: STATISTICS_ACF_DIRECT, STATISTICS_CCF_DIRECT

!===============================================================================
CONTAINS
!===============================================================================

   !------------------------------------------------------
   SUBROUTINE STATISTICS_ACF_DIRECT ( X, ACF, NORMALIZED )
   !------------------------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),  INTENT(IN)  :: X
   REAL ( KIND = DP ), DIMENSION(0:), INTENT(OUT) :: ACF

   ! . Scalar optional arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: NORMALIZED

   ! . Local scalars.
   INTEGER            :: T, TCOR, TMAX, TRUN
   REAL ( KIND = DP ) :: AVERAGE

   ! . Get the sizes of the arrays.
   TCOR = UBOUND ( ACF, 1 )
   TRUN = SIZE ( X )

   ! . Check the array sizes.
   IF ( TRUN < ( TCOR + 1 ) ) THEN
      ACF  = HUGE ( 0.0_DP )
      TCOR = TRUN - 1
   END IF

   ! . Loop over the elements in the ACF.
   DO T = 0,TCOR

      ! . Calculate the number of points in the sum.
      TMAX = TRUN - T

      ! . Calculate the element of the ACF.
      ACF(T) = DOT_PRODUCT ( X(1:TMAX), X(T+1:T+TMAX) ) / REAL ( TMAX, DP )

   END DO

   ! . Check to see whether the function is to be normalized.
   IF ( PRESENT ( NORMALIZED ) ) THEN

      ! . Normalization is required.
      IF ( NORMALIZED ) THEN

         ! . Calculate the average.
         AVERAGE = SUM ( X(1:TRUN) ) / REAL ( TRUN, DP )

         ! . Normalize the function.
         ACF(0:TCOR) = ( ACF(0:TCOR) - AVERAGE * AVERAGE ) / ( ACF(0) - AVERAGE * AVERAGE )

      END IF

   END IF

   END SUBROUTINE STATISTICS_ACF_DIRECT

   !---------------------------------------------------------
   SUBROUTINE STATISTICS_CCF_DIRECT ( X, Y, CCF, NORMALIZED )
   !---------------------------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),  INTENT(IN)  :: X, Y
   REAL ( KIND = DP ), DIMENSION(0:), INTENT(OUT) :: CCF

   ! . Scalar optional arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: NORMALIZED

   ! . Local scalars.
   INTEGER            :: T, TCOR, TMAX, TRUN
   REAL ( KIND = DP ) :: XAVERAGE, XSD, YAVERAGE, YSD

   ! . Get the sizes of the arrays.
   TCOR = UBOUND ( CCF, 1 )
   TRUN = MIN ( SIZE ( X ), SIZE ( Y ) )

   ! . Check the array sizes.
   IF ( TRUN < ( TCOR + 1 ) ) THEN
      CCF  = HUGE ( 0.0_DP )
      TCOR = TRUN - 1
   END IF

   ! . Loop over the elements in the ACF.
   DO T = 0,TCOR

      ! . Calculate the number of points in the sum.
      TMAX = TRUN - T

      ! . Calculate the element of the ACF.
      CCF(T) = DOT_PRODUCT ( X(1:TMAX), Y(T+1:T+TMAX) ) / REAL ( TMAX, DP )

   END DO

   ! . Check to see whether the function is to be normalized.
   IF ( PRESENT ( NORMALIZED ) ) THEN

      ! . Normalization is required.
      IF ( NORMALIZED ) THEN

         ! . Calculate the averages of X and Y.
         XAVERAGE = SUM ( X(1:TRUN) ) / REAL ( TRUN, DP )
         YAVERAGE = SUM ( Y(1:TRUN) ) / REAL ( TRUN, DP )

         ! . Calculate the standard deviations of X and Y.
	 XSD = SQRT ( SUM ( X(1:TRUN)**2 ) / REAL ( TRUN, DP ) - XAVERAGE * XAVERAGE )
	 YSD = SQRT ( SUM ( Y(1:TRUN)**2 ) / REAL ( TRUN, DP ) - YAVERAGE * YAVERAGE )

         ! . Normalize the function.
         CCF(0:TCOR) = ( CCF(0:TCOR) - XAVERAGE * YAVERAGE ) / ( XSD * YSD )

      END IF

   END IF

   END SUBROUTINE STATISTICS_CCF_DIRECT

END MODULE STATISTICS
