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
!                                The Files Module
!===============================================================================
!
! . Functions:
!
!   NEXT_UNIT                      Get the next available FORTRAN stream number.
!
! . Notes:
!
!   The Fortran 90 INQUIRE statement is used to find the next unit which is not
!   open. The input and output streams are always assumed to be assigned and
!   are not checked.
!
!===============================================================================
MODULE FILES

! . Module declarations.
USE DEFINITIONS, ONLY : MAX_UNITS
USE IO_UNITS,    ONLY : INPUT, OUTPUT

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !---------------------
   FUNCTION NEXT_UNIT ( )
   !---------------------

   ! . Function declarations.
   INTEGER :: NEXT_UNIT

   ! . Local scalars.
   INTEGER :: I, UNIT
   LOGICAL :: QOPEN

   ! . Initialize the unit number.
   UNIT = -1

   ! . Loop over the number of units available.
   DO I = 10, MAX_UNITS

      ! . Check all units except for the standard input and output streams.
      IF ( ( I /= INPUT ) .AND. ( I /= OUTPUT ) ) THEN

         ! . Inquire if the unit is connected.
         INQUIRE ( I, OPENED = QOPEN )

         ! . If the unit is unconnected assign the unit and exit the loop.
         IF ( .NOT. QOPEN ) THEN
            UNIT = I ; EXIT
         END IF

      END IF

   END DO

   ! . Assign a value to the function.
   NEXT_UNIT = UNIT

   END FUNCTION NEXT_UNIT

END MODULE FILES
