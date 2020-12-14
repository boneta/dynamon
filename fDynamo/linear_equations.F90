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
!                          The Linear Equations Module
!===============================================================================
!
! . Subroutines:
!
!   LINEAR_EQUATIONS_SOLVE               Solve a set of linear equations.
!
! . Notes:
!
!   The subroutine in this module acts a front end to the LAPACK subroutine
!   DGESV.
!
!===============================================================================
MODULE LINEAR_EQUATIONS

! . Module declarations.
USE DEFINITIONS, ONLY : DP

! . Remaining statements.
IMPLICIT NONE
PRIVATE
PUBLIC :: LINEAR_EQUATIONS_SOLVE

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------------------------
   SUBROUTINE LINEAR_EQUATIONS_SOLVE ( N, A, IA, B, IERR )
   !------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN)  :: IA, N
   INTEGER, INTENT(OUT) :: IERR

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:IA,1:N), INTENT(INOUT) :: A
   REAL ( KIND = DP ), DIMENSION(1:N),      INTENT(INOUT) :: B

   ! . Local scalars.
   INTEGER :: INFO

   ! . Local arrays.
   INTEGER, DIMENSION(1:N) :: IPIV

   ! . Initialize INFO.
   INFO = 0

   ! . Solve the linear equations.
   CALL DGESV ( N, 1, A(1:IA,1:N), IA, IPIV(1:N), B(1:N), N, INFO )

   ! . Set IERR.
   IERR = INFO

   END SUBROUTINE LINEAR_EQUATIONS_SOLVE

END MODULE LINEAR_EQUATIONS
