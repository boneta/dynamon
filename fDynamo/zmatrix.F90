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
!                              The Z-Matrix Module
!===============================================================================
!
! . Scalar data:
!
!   NCARDS                          The number of Z-matrix cards.
!
! . Array data:
!
!   ZMCARDS                         The Z-matrix cards.
!
! . Subroutines:
!
!   ZMATRIX_ALLOCATE                Allocate the atoms variables.
!   ZMATRIX_BUILD                   Build the coordinates from the Z-matrix.
!   ZMATRIX_CHECK                   Check the syntax of a Z-matrix.
!   ZMATRIX_FILL                    Fill the Z-matrix from the coordinates.
!   ZMATRIX_INITIALIZE              Initialize the atoms variables.
!
! . Notes:
!
!   There is no dummy atom facility as is common in many Z-matrices. All
!   the atoms must be explicitly defined.
!
!   The J, K and L counters for each Z-matrix card refers to the card
!   order. The IORDER counter gives the index of the atom referred to
!   on the card in the ATOMS data structure.
!
!===============================================================================
MODULE ZMATRIX

! . Module declarations.
USE CONSTANTS,      ONLY : TO_RADIANS, UNDEFINED
USE DEFINITIONS,    ONLY : DP
USE PRINTING,       ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH

USE ATOMS
USE GEOMETRY,       ONLY : GEOMETRY_ANGLE, GEOMETRY_COORDINATE_IJ, GEOMETRY_COORDINATE_IJK, &
                           GEOMETRY_COORDINATE_IJKL, GEOMETRY_DIHEDRAL, GEOMETRY_DISTANCE
USE TRANSFORMATION, ONLY : TO_PRINCIPAL_AXES

IMPLICIT NONE
PUBLIC
#ifndef PGPC
SAVE
#endif

! . Type definitions.
TYPE ZMATRIX_CARD
   INTEGER            :: IORDER, J, K, L
   REAL ( KIND = DP ) :: BOND, ANGLE, DIHEDRAL
END TYPE ZMATRIX_CARD

! . Scalar data.
INTEGER :: NCARDS = 0

! . Array data.
TYPE(ZMATRIX_CARD), ALLOCATABLE, DIMENSION(:) :: ZMCARDS

!==============================================================================
CONTAINS
!==============================================================================

   !--------------------------
   SUBROUTINE ZMATRIX_ALLOCATE
   !--------------------------

   ! . Set the NCARDS counter.
   NCARDS = NATOMS

   ! . Allocate ZMCARDS.
   ALLOCATE ( ZMCARDS(1:NCARDS) )

   END SUBROUTINE ZMATRIX_ALLOCATE

   !-----------------------
   SUBROUTINE ZMATRIX_BUILD
   !-----------------------

   ! . Local scalars.
   INTEGER            :: I, ICARD, J, JCARD, K, KCARD, L, LCARD
   REAL ( KIND = DP ) :: PHI, R, THETA

   ! . Check the number of atoms and cards.
   IF ( ( NATOMS <= 0 ) .OR. ( NCARDS <= 0 ) ) RETURN

   ! . Reallocate and initialize ATMCRD.
   IF ( ALLOCATED ( ATMCRD ) ) DEALLOCATE ( ATMCRD )
   ALLOCATE ( ATMCRD(1:3,1:NATOMS) ) ; ATMCRD = UNDEFINED

   ! . Place the first atom at the origin.
   I = ZMCARDS(1)%IORDER
   ATMCRD(1:3,I) = 0.0_DP

   ! . There is more than one card.
   IF ( NCARDS > 1 ) THEN

      ! . Place the remaining atoms.
      DO ICARD = 2,NCARDS

         ! . Get the card indices.
         JCARD = ZMCARDS(ICARD)%J
         KCARD = ZMCARDS(ICARD)%K
         LCARD = ZMCARDS(ICARD)%L

         ! . Get the atom indices.
         I = ZMCARDS(ICARD)%IORDER
         IF ( JCARD > 0 ) THEN ; J = ZMCARDS(JCARD)%IORDER ; ELSE ; J = 0 ; END IF
         IF ( KCARD > 0 ) THEN ; K = ZMCARDS(KCARD)%IORDER ; ELSE ; K = 0 ; END IF
         IF ( LCARD > 0 ) THEN ; L = ZMCARDS(LCARD)%IORDER ; ELSE ; L = 0 ; END IF

         ! . Get the internal coordinates.
         R     = ZMCARDS(ICARD)%BOND
         THETA = ZMCARDS(ICARD)%ANGLE
         PHI   = ZMCARDS(ICARD)%DIHEDRAL

         ! . Calculate the coordinate for atom I.
         SELECT CASE ( ICARD )
         CASE ( 2 )   ; CALL GEOMETRY_COORDINATE_IJ   ( ATMCRD, I, J,       R,        "X" )
         CASE ( 3 )   ; CALL GEOMETRY_COORDINATE_IJK  ( ATMCRD, I, J, K,    R, THETA, "Y" )
         CASE DEFAULT ; CALL GEOMETRY_COORDINATE_IJKL ( ATMCRD, I, J, K, L, R, THETA, PHI )
         END SELECT

      END DO
   END IF

   ! . Transform the coordinates to their principal axes if all coordinates are defined.
   IF ( .NOT. ANY ( ATMCRD == UNDEFINED ) ) CALL TO_PRINCIPAL_AXES ( ATMCRD )

   END SUBROUTINE ZMATRIX_BUILD

   !-----------------------
   SUBROUTINE ZMATRIX_CHECK
   !-----------------------

   ! . Local scalars.
   INTEGER            :: I, ICARD, JCARD, KCARD, LCARD
   LOGICAL            :: QERROR, QWRONG
   REAL ( KIND = DP ) :: ANGLE, BOND, DIHEDRAL

   ! . Check the number of cards.
   IF ( NCARDS <= 0 ) RETURN

   ! . Initialize the error flag.
   QERROR = .FALSE.

   ! . Loop over the cards.
   DO ICARD = 1,NCARDS

      ! . Get the atom indices.
      I = ZMCARDS(ICARD)%IORDER

      ! . Get the card indices.
      JCARD = ZMCARDS(ICARD)%J
      KCARD = ZMCARDS(ICARD)%K
      LCARD = ZMCARDS(ICARD)%L

      ! . Get the internal coordintes.
      BOND     = ZMCARDS(ICARD)%BOND
      ANGLE    = ZMCARDS(ICARD)%ANGLE
      DIHEDRAL = ZMCARDS(ICARD)%DIHEDRAL

      ! . Do some general error checking.
      QWRONG = ( I <= 0 ) .OR. ( I > NCARDS ) .OR. ( BOND < 0.0_DP ) .OR. ( ANGLE < 0.0_DP ) .OR. &
                                    ( ( ANGLE    > 180.0_DP ) .AND. ( ANGLE /= UNDEFINED ) ) .OR. &
                                                                    ( DIHEDRAL < -180.0_DP ) .OR. &
                                    ( ( DIHEDRAL > 180.0_DP ) .AND. ( DIHEDRAL /= UNDEFINED ) )

      ! . Do some card specific error checking.
      SELECT CASE ( ICARD )
      CASE ( 1 )
         QWRONG = QWRONG .OR. ( JCARD /= 0 ) .OR. ( KCARD /= 0 ) .OR. ( LCARD /= 0 ) .OR. &
                    ( BOND /= 0.0_DP ) .OR. ( ANGLE /= 0.0_DP ) .OR. ( DIHEDRAL /= 0.0_DP )
      CASE ( 2 )
         QWRONG = QWRONG .OR. ( JCARD /= 1 ) .OR. ( KCARD /= 0 ) .OR. ( LCARD /= 0 ) .OR. &
                                            ( ANGLE /= 0.0_DP ) .OR. ( DIHEDRAL /= 0.0_DP )
      CASE ( 3 )
         QWRONG = QWRONG .OR. .NOT. ( ( ( JCARD /= 1 ) .AND. ( KCARD /= 2 ) )   .OR. &
                                      ( ( JCARD /= 2 ) .AND. ( KCARD /= 1 ) ) ) .OR. &
                                            ( LCARD /= 0 ) .OR. ( DIHEDRAL /= 0.0_DP )
      CASE DEFAULT
         QWRONG = QWRONG .OR. ( JCARD >= ICARD ) .OR. ( KCARD >= ICARD ) .OR. ( LCARD >= ICARD ) .OR. &
                              ( JCARD == KCARD ) .OR. ( JCARD == LCARD ) .OR. ( KCARD == LCARD )
      END SELECT

      ! . Write out a message if there is an error.
      IF ( QWRONG ) THEN
         WRITE ( PRINT_LINE, "(A,I6,A)" ) "Syntax error in card = ", ICARD, "."
	 CALL PRINT_PARAGRAPH
      END IF

      ! . Set the error flag.
      QERROR = QERROR .OR. QWRONG

   END DO

   ! . Check for an error.
   IF ( QERROR ) CALL PRINT_ERROR ( "ZMATRIX_CHECK", "Invalid Z-matrix." )

   END SUBROUTINE ZMATRIX_CHECK

   !----------------------
   SUBROUTINE ZMATRIX_FILL
   !----------------------

   ! . Local scalars.
   INTEGER :: I, ICARD, J, JCARD, K, KCARD, L, LCARD

   ! . Check the number of atoms and cards.
   IF ( ( NATOMS <= 0 ) .OR. ( NCARDS <= 0 ) ) RETURN

   ! . Loop over the cards.
   DO ICARD = 1,NCARDS

      ! . Get the card indices.
      JCARD = ZMCARDS(ICARD)%J
      KCARD = ZMCARDS(ICARD)%K
      LCARD = ZMCARDS(ICARD)%L

      ! . Get the atom indices.
      I = ZMCARDS(ICARD)%IORDER
      IF ( JCARD > 0 ) THEN ; J = ZMCARDS(JCARD)%IORDER ; ELSE ; J = 0 ; END IF
      IF ( KCARD > 0 ) THEN ; K = ZMCARDS(KCARD)%IORDER ; ELSE ; K = 0 ; END IF
      IF ( LCARD > 0 ) THEN ; L = ZMCARDS(LCARD)%IORDER ; ELSE ; L = 0 ; END IF

      ! . Get the internal coordinates.
      ZMCARDS(ICARD)%BOND     = GEOMETRY_DISTANCE ( ATMCRD, I, J       )
      ZMCARDS(ICARD)%ANGLE    = GEOMETRY_ANGLE    ( ATMCRD, I, J, K    )
      ZMCARDS(ICARD)%DIHEDRAL = GEOMETRY_DIHEDRAL ( ATMCRD, I, J, K, L )

   END DO

   ! . Zero the non-existent internal coordinates on the first three cards.
   ZMCARDS(1)%BOND     = 0.0_DP
   ZMCARDS(1)%ANGLE    = 0.0_DP
   ZMCARDS(1)%DIHEDRAL = 0.0_DP
   IF ( NCARDS > 1 ) THEN
      ZMCARDS(2)%ANGLE    = 0.0_DP
      ZMCARDS(2)%DIHEDRAL = 0.0_DP
      IF ( NCARDS > 2 ) ZMCARDS(3)%DIHEDRAL = 0.0_DP
   END IF

   END SUBROUTINE ZMATRIX_FILL

   !----------------------------
   SUBROUTINE ZMATRIX_INITIALIZE
   !----------------------------

   ! . Initialize NCARDS.
   NCARDS = 0

   ! . Initialize the ZMCARDS array.
   IF ( ALLOCATED ( ZMCARDS ) ) DEALLOCATE ( ZMCARDS )

   END SUBROUTINE ZMATRIX_INITIALIZE

END MODULE ZMATRIX
