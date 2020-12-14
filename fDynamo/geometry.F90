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
!                              The Geometry Module
!===============================================================================
!
! . Functions:
!
!   GEOMETRY_ANGLE                 Calculate an angle.
!   GEOMETRY_ANGLES                Calculate some angles.
!   GEOMETRY_DIHEDRAL              Calculate a dihedral.
!   GEOMETRY_DIHEDRALS             Calculate some dihedrals.
!   GEOMETRY_DISTANCE              Calculate a distance.
!   GEOMETRY_DISTANCES             Calculate some distances.
!
! . Subroutines:
!
!   GEOMETRY_COORDINATE_IJ         Calculate the coordinates of an atom given
!                                  the coordinates of another atom and a
!                                  distance. The atom is put in a random
!                                  direction or along one of the axes.
!   GEOMETRY_COORDINATE_IJK        Calculate the coordinates of an atom given
!                                  the coordinates of two other atoms, a
!                                  distance and an angle. The atom is placed
!                                  so that its dihedral is 0 with respect to
!                                  a random direction or one of the axes.
!   GEOMETRY_COORDINATE_IJKL       Calculate the coordinates of an atom given
!                                  the coordinates of three other atoms and a
!                                  distance, angle and dihedral.
!
! . Notes:
!
!   GEOMETRY_ANGLE, GEOMETRY_DIHEDRAL and GEOMETRY_DISTANCE calculate the
!   angle, dihedral and distance for a list of atoms (3, 4 and 2 respectively). 
!
!   GEOMETRY_ANGLES, GEOMETRY_DIHEDRALS and GEOMETRY_DISTANCES carry out the
!   same function except that they are vector functions are return many such
!   values. Each takes an array, LIST(1:n,1:x),  which is a list of atom
!   indices where n is the number of terms to be calculated and x is 2, 3
!   and 4 for distances, angles and dihedrals respectively.
!
!   Any errors result in a nonsensical value, UNDEFINED, being returned for
!   the internal coordinate in question.
!
!===============================================================================
MODULE GEOMETRY

! . Module declarations.
USE CONSTANTS,      ONLY : TO_DEGREES, TO_RADIANS, UNDEFINED
USE DEFINITIONS,    ONLY : DP
USE LINEAR_ALGEBRA, ONLY : CROSS_PRODUCT, NORMALIZE
USE RANDOM_NUMBERS, ONLY : RANDOM_VECTOR

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------------------------
   FUNCTION GEOMETRY_ANGLE ( COORDINATES, I, J, K )
   !-----------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: GEOMETRY_ANGLE

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: I, J, K

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Local scalars.
   INTEGER            :: N
   REAL ( KIND = DP ) :: DOTFAC, RIJ, RKJ

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ, DRKJ

   ! . Initialize the angle.
   GEOMETRY_ANGLE = UNDEFINED

   ! . Check the first dimension of COORDINATES.
   IF ( SIZE ( COORDINATES, 1 ) /= 3 ) RETURN

   ! . Get the second dimension.
   N = SIZE ( COORDINATES, 2 )

   ! . Check I, J and K.
   IF ( ( I <= 0 ) .OR. ( I > N ) .OR. ( J <= 0 ) .OR. ( J > N ) .OR. ( K <= 0 ) .OR. ( K > N ) ) RETURN

   ! . Check the coordinates.
   IF ( ( ANY ( COORDINATES(1:3,I) == UNDEFINED ) ) .OR. &
        ( ANY ( COORDINATES(1:3,J) == UNDEFINED ) ) .OR. &
        ( ANY ( COORDINATES(1:3,K) == UNDEFINED ) ) ) RETURN

   ! . Calculate some displacement vectors.
   DRIJ = COORDINATES(1:3,I) - COORDINATES(1:3,J)
   DRKJ = COORDINATES(1:3,K) - COORDINATES(1:3,J)

   ! . Calculate the size of the vectors.
   RIJ = SQRT ( DOT_PRODUCT ( DRIJ, DRIJ ) )
   RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

   ! . The distances are non-zero.
   IF ( ( RIJ /= 0.0_DP ) .AND. ( RKJ /= 0.0_DP ) ) THEN

      ! . Normalize the vectors.
      DRIJ = DRIJ / RIJ
      DRKJ = DRKJ / RKJ

      ! . Calculate the dot product of the vectors.
      DOTFAC = DOT_PRODUCT ( DRIJ, DRKJ )

      ! . Ensure DOTFAC is between -1 and 1.
      DOTFAC = SIGN ( MIN ( ABS ( DOTFAC ), 1.0_DP ), DOTFAC )

      ! . Calculate the angle.
      GEOMETRY_ANGLE = TO_DEGREES * ACOS ( DOTFAC )

   END IF

   END FUNCTION GEOMETRY_ANGLE

   !---------------------------------------------
   FUNCTION GEOMETRY_ANGLES ( COORDINATES, LIST )
   !---------------------------------------------

   ! . Array arguments.
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: LIST
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:SIZE(LIST,2)) :: GEOMETRY_ANGLES

   ! . Local scalars.
   INTEGER :: I, J, K, N, P

   ! . Initialize GEOMETRY_ANGLES.
   GEOMETRY_ANGLES = UNDEFINED

   ! . Check the dimensions of LIST.
   N = SIZE ( LIST, 2 )
   IF ( SIZE ( LIST, 1 ) /= 3 ) RETURN

   ! . Loop over the elements in the list.
   DO P = 1,N

      ! . Get the atom indices.
      I = LIST(1,P)
      J = LIST(2,P)
      K = LIST(3,P)

      ! . Calculate the angle.
      GEOMETRY_ANGLES(P) = GEOMETRY_ANGLE ( COORDINATES, I, J, K )

   END DO

   END FUNCTION GEOMETRY_ANGLES

   !-----------------------------------------------------
   FUNCTION GEOMETRY_DIHEDRAL ( COORDINATES, I, J, K, L )
   !-----------------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: GEOMETRY_DIHEDRAL

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: I, J, K, L

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Local scalars.
   INTEGER            :: N
   REAL ( KIND = DP ) :: DOTFAC, R, RKJ, S, SGNFAC

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DN, DR, DRIJ, DRKJ, DRLK, DS

   ! . Initialize the dihedral.
   GEOMETRY_DIHEDRAL = UNDEFINED

   ! . Check the first dimension of COORDINATES.
   IF ( SIZE ( COORDINATES, 1 ) /= 3 ) RETURN

   ! . Get the second dimension.
   N = SIZE ( COORDINATES, 2 )

   ! . Check I, J, K and L.
   IF ( ( I <= 0 ) .OR. ( I > N ) .OR. ( J <= 0 ) .OR. ( J > N ) .OR. &
        ( K <= 0 ) .OR. ( K > N ) .OR. ( L <= 0 ) .OR. ( L > N ) ) RETURN

   ! . Check the coordinates.
   IF ( ( ANY ( COORDINATES(1:3,I) == UNDEFINED ) ) .OR. &
        ( ANY ( COORDINATES(1:3,J) == UNDEFINED ) ) .OR. &
        ( ANY ( COORDINATES(1:3,K) == UNDEFINED ) ) .OR. &
        ( ANY ( COORDINATES(1:3,L) == UNDEFINED ) ) ) RETURN

   ! . Calculate some displacement vectors.
   DRIJ = COORDINATES(1:3,I) - COORDINATES(1:3,J)
   DRKJ = COORDINATES(1:3,K) - COORDINATES(1:3,J)
   DRLK = COORDINATES(1:3,L) - COORDINATES(1:3,K)

   ! . Calculate the size of the RKJ vector.
   RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

   ! . The distance is non-zero.
   IF ( RKJ /= 0.0_DP ) THEN

      ! . Normalize the vector.
      DRKJ = DRKJ / RKJ

      ! . Calculate the DR and DS vectors.
      DR = DRIJ - DOT_PRODUCT ( DRIJ, DRKJ ) * DRKJ
      DS = DRLK - DOT_PRODUCT ( DRLK, DRKJ ) * DRKJ

      ! . Calculate the magnitudes of DR and DS.
      R = SQRT ( DOT_PRODUCT ( DR, DR ) )
      S = SQRT ( DOT_PRODUCT ( DS, DS ) )

      ! . The distances are non-zero.
      IF ( ( R /= 0.0_DP ) .AND. ( S /= 0.0_DP ) ) THEN

         ! . Normalize DR and DS.
         DR = DR / R
         DS = DS / S

         ! . Calculate the cross-product of DRKJ and DRLK.
         DN = CROSS_PRODUCT ( DRKJ, DRLK )

         ! . Calculate the factor to determine the sign of the dihedral.
         SGNFAC = 1.0_DP
         IF ( - DOT_PRODUCT ( DRIJ, DN ) < 0.0_DP ) SGNFAC = -1.0_DP

         ! . Calculate the dot product of the vectors.
         DOTFAC = DOT_PRODUCT ( DR, DS )

         ! . Ensure DOTFAC is between -1 and 1.
         DOTFAC = SIGN ( MIN ( ABS ( DOTFAC ), 1.0_DP ), DOTFAC )

         ! . Calculate the dihedral.
         GEOMETRY_DIHEDRAL = TO_DEGREES * SGNFAC * ACOS ( DOTFAC )

      END IF

   END IF

   END FUNCTION GEOMETRY_DIHEDRAL

   !------------------------------------------------
   FUNCTION GEOMETRY_DIHEDRALS ( COORDINATES, LIST )
   !------------------------------------------------

   ! . Array arguments.
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: LIST
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:SIZE(LIST,2)) :: GEOMETRY_DIHEDRALS

   ! . Local scalars.
   INTEGER :: I, J, K, L, N, P

   ! . Initialize GEOMETRY_DIHEDRALS.
   GEOMETRY_DIHEDRALS = UNDEFINED

   ! . Check the dimensions of LIST.
   N = SIZE ( LIST, 2 )
   IF ( SIZE ( LIST, 1 ) /= 4 ) RETURN

   ! . Loop over the elements in the list.
   DO P = 1,N

      ! . Get the atom indices.
      I = LIST(1,P)
      J = LIST(2,P)
      K = LIST(3,P)
      L = LIST(4,P)

      ! . Calculate the angle.
      GEOMETRY_DIHEDRALS(P) = GEOMETRY_DIHEDRAL ( COORDINATES, I, J, K, L )

   END DO

   END FUNCTION GEOMETRY_DIHEDRALS

   !-----------------------------------------------
   FUNCTION GEOMETRY_DISTANCE ( COORDINATES, I, J )
   !-----------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: GEOMETRY_DISTANCE

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: I, J

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Local scalars.
   INTEGER :: N

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ

   ! . Initialize the angle.
   GEOMETRY_DISTANCE = UNDEFINED

   ! . Check the first dimension of COORDINATES.
   IF ( SIZE ( COORDINATES, 1 ) /= 3 ) RETURN

   ! . Get the second dimension.
   N = SIZE ( COORDINATES, 2 )

   ! . Check I and J.
   IF ( ( I <= 0 ) .OR. ( I > N ) .OR. ( J <= 0 ) .OR. ( J > N ) ) RETURN

   ! . Check the coordinates.
   IF ( ( ANY ( COORDINATES(1:3,I) == UNDEFINED ) ) .OR. &
        ( ANY ( COORDINATES(1:3,J) == UNDEFINED ) ) ) RETURN

   ! . Calculate the displacement vector.
   DRIJ = COORDINATES(1:3,I) - COORDINATES(1:3,J)

   ! . Calculate the size of the vector.
   GEOMETRY_DISTANCE = SQRT ( DOT_PRODUCT ( DRIJ, DRIJ ) )

   END FUNCTION GEOMETRY_DISTANCE

   !------------------------------------------------
   FUNCTION GEOMETRY_DISTANCES ( COORDINATES, LIST )
   !------------------------------------------------

   ! . Array arguments.
   INTEGER,            DIMENSION(:,:), INTENT(IN) :: LIST
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:SIZE(LIST,2)) :: GEOMETRY_DISTANCES

   ! . Local scalars.
   INTEGER :: I, J, N, P

   ! . Initialize GEOMETRY_DISTANCES.
   GEOMETRY_DISTANCES = UNDEFINED

   ! . Check the dimensions of LIST.
   N = SIZE ( LIST, 2 )
   IF ( SIZE ( LIST, 1 ) /= 2 ) RETURN

   ! . Loop over the elements in the list.
   DO P = 1,N

      ! . Get the atom indices.
      I = LIST(1,P)
      J = LIST(2,P)

      ! . Calculate the angle.
      GEOMETRY_DISTANCES(P) = GEOMETRY_DISTANCE ( COORDINATES, I, J )

   END DO

   END FUNCTION GEOMETRY_DISTANCES

   !-----------------------------------------------------------------
   SUBROUTINE GEOMETRY_COORDINATE_IJ ( COORDINATES, I, J, R, OPTION )
   !-----------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = 1 ), INTENT(IN) :: OPTION
   INTEGER,               INTENT(IN) :: I, J
   REAL ( KIND = DP ),    INTENT(IN) :: R

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: COORDINATES

   ! . Local scalars.
   INTEGER :: N

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DIRECTION

   ! . Get the length of the coordinate array.
   N = SIZE ( COORDINATES, 2 )

   ! . If COORDINATES or I is invalid destroy the whole coordinate array and return.
   IF ( ( SIZE ( COORDINATES, 1 ) /= 3 ) .OR. ( I <= 0 ) .OR. ( I > N ) ) THEN
      COORDINATES = UNDEFINED
      RETURN
   END IF

   ! . Initialize the coordinates for I.
   COORDINATES(1:3,I) = UNDEFINED

   ! . Check J.
   IF ( ( J <= 0 ) .OR. ( J > N ) .OR. ( ANY ( COORDINATES(1:3,J) == UNDEFINED ) ) ) RETURN

   ! . Initialize the direction.
   DIRECTION = 0.0_DP

   ! . Branch on the option.
   SELECT CASE ( OPTION )
   CASE ( "R" ) ; DIRECTION = NORMALIZE ( RANDOM_VECTOR ( 3 ) )
   CASE ( "X" ) ; DIRECTION(1) =  1.0_DP
   CASE ( "Y" ) ; DIRECTION(2) =  1.0_DP
   CASE ( "Z" ) ; DIRECTION(3) =  1.0_DP
   CASE DEFAULT ; RETURN
   END SELECT

   ! . Calculate the coordinates of I.
   COORDINATES(1:3,I) = COORDINATES(1:3,J) + R * DIRECTION

   END SUBROUTINE GEOMETRY_COORDINATE_IJ

   !----------------------------------------------------------------------------
   SUBROUTINE GEOMETRY_COORDINATE_IJK ( COORDINATES, I, J, K, R, THETA, OPTION )
   !----------------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = 1 ), INTENT(IN) :: OPTION
   INTEGER,               INTENT(IN) :: I, J, K
   REAL ( KIND = DP ),    INTENT(IN) :: R, THETA

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: COORDINATES

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: SMALL = 1.0E-9_DP

   ! . Local scalars.
   INTEGER            :: N
   REAL ( KIND = DP ) :: COST, RA, RKJ, SINT, WA, WB

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRA, DRKJ, DRLK

   ! . Get the length of the coordinate array.
   N = SIZE ( COORDINATES, 2 )

   ! . If COORDINATES or I is invalid destroy the whole coordinate array and return.
   IF ( ( SIZE ( COORDINATES, 1 ) /= 3 ) .OR. ( I <= 0 ) .OR. ( I > N ) ) THEN
      COORDINATES = UNDEFINED
      RETURN
   END IF

   ! . Initialize the coordinates for I.
   COORDINATES(1:3,I) = UNDEFINED

   ! . Check J and K.
   IF ( ( J <= 0 ) .OR. ( J > N ) .OR. ( ANY ( COORDINATES(1:3,J) == UNDEFINED ) ) .OR. &
        ( K <= 0 ) .OR. ( K > N ) .OR. ( ANY ( COORDINATES(1:3,K) == UNDEFINED ) ) ) RETURN

   ! . Initialize the direction DRLK.
   DRLK = 0.0_DP

   ! . Branch on the option.
   SELECT CASE ( OPTION )
   CASE ( "R" ) ; DRLK    = NORMALIZE ( RANDOM_VECTOR ( 3 ) )
   CASE ( "X" ) ; DRLK(1) =  1.0_DP
   CASE ( "Y" ) ; DRLK(2) =  1.0_DP
   CASE ( "Z" ) ; DRLK(3) =  1.0_DP
   CASE DEFAULT ; RETURN
   END SELECT

   ! . Calculate some trigonometric functions.
   COST = COS ( TO_RADIANS * THETA )
   SINT = SIN ( TO_RADIANS * THETA )

   ! . Calculate the JK displacement vector.
   DRKJ = COORDINATES(1:3,K) - COORDINATES(1:3,J)

   ! . Calculate the length of DRKJ.
   RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

   ! . Check that J and K are not too close.
   IF ( RKJ > SMALL ) THEN

      ! . Normalize DRKJ.
      DRKJ = DRKJ / RKJ

      ! . Calculate the cross product of DRLK and DRKJ.
      DRA = CROSS_PRODUCT ( DRLK, DRKJ )

      ! . Calculate the length of DRA.
      RA = SQRT ( DOT_PRODUCT ( DRA, DRA ) )

      ! . Check the length of DRA.
      IF ( RA > SMALL ) THEN

         ! . Normalize DRA.
         DRA = DRA / RA

         ! . Calculate the cross product of DRKJ and DRA (put in DRLK).
         DRLK = CROSS_PRODUCT ( DRKJ, DRA )

         ! . Calculate the coordinate displacements (spherical polars).
         WA = R * COST
         WB = R * SINT

         ! . Calculate the coordinates of I.
         COORDINATES(1:3,I) = COORDINATES(1:3,J) + WA * DRKJ + WB * DRLK

      END IF
   END IF

   END SUBROUTINE GEOMETRY_COORDINATE_IJK

   !-----------------------------------------------------------------------------
   SUBROUTINE GEOMETRY_COORDINATE_IJKL ( COORDINATES, I, J, K, L, R, THETA, PHI )
   !-----------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: I, J, K, L
   REAL ( KIND = DP ), INTENT(IN) :: R, THETA, PHI 

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: COORDINATES

   ! . Parameter declarations.
   REAL ( KIND = DP ), PARAMETER :: SMALL = 1.0E-9_DP

   ! . Local scalars.
   INTEGER            :: N
   REAL ( KIND = DP ) :: COST, SINT, COSP, SINP, RA, RKJ, WA, WB, WC

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRA, DRKJ, DRLK

   ! . Get the length of the coordinate array.
   N = SIZE ( COORDINATES, 2 )

   ! . If COORDINATES or I is invalid destroy the whole coordinate array and return.
   IF ( ( SIZE ( COORDINATES, 1 ) /= 3 ) .OR. ( I <= 0 ) .OR. ( I > N ) ) THEN
      COORDINATES = UNDEFINED
      RETURN
   END IF

   ! . Initialize the coordinates for I.
   COORDINATES(1:3,I) = UNDEFINED

   ! . Check J, K and L.
   IF ( ( J <= 0 ) .OR. ( J > N ) .OR. ( ANY ( COORDINATES(1:3,J) == UNDEFINED ) ) .OR. &
        ( K <= 0 ) .OR. ( K > N ) .OR. ( ANY ( COORDINATES(1:3,K) == UNDEFINED ) ) .OR. &
        ( L <= 0 ) .OR. ( L > N ) .OR. ( ANY ( COORDINATES(1:3,L) == UNDEFINED ) ) ) RETURN

   ! . Calculate some trigonometric functions.
   COSP = COS ( TO_RADIANS * PHI   )
   SINP = SIN ( TO_RADIANS * PHI   )
   COST = COS ( TO_RADIANS * THETA )
   SINT = SIN ( TO_RADIANS * THETA )

   ! . Calculate some displacement vectors.
   DRKJ = COORDINATES(1:3,K) - COORDINATES(1:3,J)
   DRLK = COORDINATES(1:3,L) - COORDINATES(1:3,K)

   ! . Calculate the length of DRKJ.
   RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

   ! . Check that J and K are not too close.
   IF ( RKJ > SMALL ) THEN

      ! . Normalize DRKJ.
      DRKJ = DRKJ / RKJ

      ! . Calculate the cross product of DRLK and DRKJ.
      DRA = CROSS_PRODUCT ( DRLK, DRKJ )

      ! . Calculate the length of DRA.
      RA = SQRT ( DOT_PRODUCT ( DRA, DRA ) )

      ! . Check the length of DRA.
      IF ( RA > SMALL ) THEN

         ! . Normalize DRA.
         DRA = DRA / RA

         ! . Calculate the cross product of DRKJ and RA (put in DRLK).
         DRLK = CROSS_PRODUCT ( DRKJ, DRA )

         ! . Calculate the coordinate displacements (spherical polars).
         WA = R * COST
         WB = R * SINT * COSP
         WC = R * SINT * SINP

         ! . Calculate the coordinates of I.
         COORDINATES(1:3,I) = COORDINATES(1:3,J) + WA * DRKJ + WB * DRLK + WC * DRA

      END IF

   END IF

   END SUBROUTINE GEOMETRY_COORDINATE_IJKL

END MODULE GEOMETRY
