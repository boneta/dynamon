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
!                             The Definitions Module
!===============================================================================
!
! . Parameter data:
!
!   DEFAULT_INPUT                   The default input stream.
!   DEFAULT_OUTPUT                  The default output stream.
!   DP                              The Default Precision for reals.
!   FORCE_FIELD                     The module library force field.
!   LINE_LENGTH                     The maximum line length.
!   MAX_RECORD_LENGTH               The maximum binary file record length.
!   MAX_UNITS                       The maximum number of files.
!   SP                              Single Precision.
!   VERSION                         The version number of the program.
!
! . Notes:
!
!   This version has been tested with the NAG compiler on PC/Linux and HP/HP-UX
!   machines and the DEC compiler on DEC Alpha systems.
!
!   A 32 bit precision is sufficient for integers. For real numbers a 64 bit
!   precision is used throughout.
!
!===============================================================================
MODULE DEFINITIONS

IMPLICIT NONE
PUBLIC

! . The default program input stream.
INTEGER, PARAMETER :: DEFAULT_INPUT = 5

! . The default program output stream.
INTEGER, PARAMETER :: DEFAULT_OUTPUT = 6

! . The default precision parameter for real numbers.
INTEGER, PARAMETER :: DP = 8

! . The module library force field.
CHARACTER ( LEN = 16 ), PARAMETER :: FORCE_FIELD = "OPLS_AA"

! . The input line length.
INTEGER, PARAMETER :: LINE_LENGTH = 132

! . The maximum record length for a binary file.
INTEGER, PARAMETER :: MAX_RECORD_LENGTH = 2**20 ! Linux.

! . The maximum FORTRAN stream number.
INTEGER, PARAMETER :: MAX_UNITS = 100

! . The single precision parameter for real numbers.
INTEGER, PARAMETER :: SP = 4

! . The version number of the module library.
REAL ( KIND = DP ), PARAMETER :: VERSION = 2.2_DP

LOGICAL :: SKIP_ABINITIO = .TRUE.

END MODULE DEFINITIONS
