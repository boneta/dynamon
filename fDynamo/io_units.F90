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
!                              The I/O Units Module
!===============================================================================
!
! . I/O Channels:
!
!   INPUT                          The input  stream unit number.
!   OUTPUT                         The output stream unit number.
!
! . Notes:
!
!   The I/O unit values are parameters and so cannot be changed.
!
!===============================================================================
MODULE IO_UNITS

! . Module declarations.
USE DEFINITIONS, ONLY : DEFAULT_INPUT, DEFAULT_OUTPUT

IMPLICIT NONE
PUBLIC

! . Integer data.
INTEGER, PARAMETER :: INPUT = DEFAULT_INPUT, OUTPUT = DEFAULT_OUTPUT

END MODULE IO_UNITS
