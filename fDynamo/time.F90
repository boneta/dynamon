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
!                                The Time Module
!===============================================================================
!
! . Subroutines:
!
!   TIME_CPU_INITIALIZE            Initialize the CPU time.
!   TIME_CPU_PRINT                 Print the CPU time.
!   TIME_PRINT                     Print out the current date and time.
!
! . Notes:
!
!   The standard Fortran 90 subroutine DATE_AND_TIME is used to get the current
!   date and time strings.
!
!   The standard Fortran 95 subroutine is used to get the CPU time. As such, it
!   can only be used when Fortran 95 compilers are available.
!
!===============================================================================
MODULE TIME

! . Module declarations.
USE PRINTING, ONLY : PRINT_LINE, PRINT_LINEBREAK, PRINT_PARAGRAPH, PRINT_PARAGRAPH_START, &
                     PRINT_PARAGRAPH_STOP, PRINT_TEXT

IMPLICIT NONE
PRIVATE
#ifdef F95
PUBLIC :: TIME_CPU_INITIALIZE, TIME_CPU_PRINT, TIME_PRINT
#else
PUBLIC :: TIME_PRINT
#endif
#ifndef PGPC
SAVE
#endif

#ifdef F95
! . Module scalars.
REAL :: CPUSTART = 0.0
#endif

!==============================================================================
CONTAINS
!==============================================================================

#ifdef F95

   !-----------------------------
   SUBROUTINE TIME_CPU_INITIALIZE
   !-----------------------------

   ! . Initialize the CPU time counter.
   CALL CPU_TIME ( CPUSTART )

   END SUBROUTINE TIME_CPU_INITIALIZE

   !------------------------
   SUBROUTINE TIME_CPU_PRINT
   !------------------------

   ! . Local scalars.
   REAL :: CPUACTUAL

   ! . Get the CPU time counter.
   CALL CPU_TIME ( CPUACTUAL )

   ! . Print the CPU time.
   WRITE ( PRINT_LINE, "(A,F12.3,A)" ) "CPU time used = ", CPUACTUAL - CPUSTART, " seconds."
   CALL PRINT_PARAGRAPH

   END SUBROUTINE TIME_CPU_PRINT

#endif

   !--------------------
   SUBROUTINE TIME_PRINT
   !--------------------

   ! . Local scalars.
   CHARACTER ( LEN =  8 ) :: DATSTR
   CHARACTER ( LEN = 10 ) :: TIMSTR

   ! . Get the current date and time.
   CALL DATE_AND_TIME ( DATSTR, TIMSTR )

   ! . Write out the date and time.
   CALL PRINT_PARAGRAPH_START
   CALL PRINT_TEXT ( TEXT = "Date = " // DATSTR(7:8) // "/" // DATSTR(5:6) // "/" // DATSTR(1:4)  )
   CALL PRINT_LINEBREAK
   CALL PRINT_TEXT ( TEXT = "Time = " // TIMSTR(1:2) // ":" // TIMSTR(3:4) // ":" // TIMSTR(5:10) )
   CALL PRINT_PARAGRAPH_STOP

   END SUBROUTINE TIME_PRINT

END MODULE TIME
