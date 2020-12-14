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
!                          The Build Coordinates Module
!===============================================================================
!
! . Subroutines:
!
!   BUILD_HYDROGENS                 Build hydrogen coordinates.
!
! . Notes:
!
!   Only a reasonable guess at building the coordinates of the desired atoms is
!   made. Refinement using energy minimization or dynamics will be required
!   afterwards.
!
!   An MM TERMS data structure is required.
!
!   Hydrogens are built using only bond and angle information - no dihedral or
!   non-bonding information. Bonds to other hydrogens are ignored and hydrogens
!   linked to more than one heavy atom will not be built.
!
!===============================================================================
MODULE BUILD_COORDINATES

! . Module declarations.
USE CONSTANTS,      ONLY : PI, UNDEFINED
USE DEFINITIONS,    ONLY : DP
USE LINEAR_ALGEBRA, ONLY : CROSS_PRODUCT, NORMALIZE
USE PRINTING,       ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, &
                           PRINT_SUMMARY_STOP
USE RANDOM_NUMBERS, ONLY : RANDOM_VECTOR

USE ATOMS,          ONLY : ATMCRD, ATMNUM, NATOMS
USE MM_TERMS,       ONLY : ANGLES, BONDS, NANGLES, NBONDS

IMPLICIT NONE
PRIVATE
PUBLIC :: BUILD_HYDROGENS

! . Module parameters.
REAL ( KIND = DP ), PARAMETER :: LINEAR          = PI,                       &
                                 TETRAHEDRAL     = 1.9106332362490185563_DP, & ! 109.47 degrees.
                                 TRIGONAL_PLANAR = 2.0_DP * PI / 3.0_DP

!==============================================================================
CONTAINS
!==============================================================================

   !-------------------------
   SUBROUTINE BUILD_HYDROGENS
   !-------------------------

   ! . Local parameters.
   INTEGER, PARAMETER :: MAX_VALENCY = 4

   ! . Local type definitions.
   TYPE CENTER_TYPE
      CHARACTER ( LEN = 16 )                       :: TYPE
      INTEGER                                      :: CENTER, NDEFINED, NUNDEFINED
      REAL ( KIND = DP )                           :: ANGLE
      INTEGER,            DIMENSION(1:MAX_VALENCY) :: ATOMS
      REAL ( KIND = DP ), DIMENSION(1:MAX_VALENCY) :: BONDS
   END TYPE CENTER_TYPE

   ! . Local scalars.
   INTEGER :: I, NCENTERS, NUNDEFA, NUNDEFH, NUNDEFH2

   ! . Local arrays.
   LOGICAL, DIMENSION(1:NATOMS) :: QUNDEFINED

   ! . Local types.
   TYPE(CENTER_TYPE), ALLOCATABLE, DIMENSION(:) :: CENTERS

   !---------------------------------------------------------------------------
   ! . Initialization.
   !---------------------------------------------------------------------------
   ! . Check the number of atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Find the number of atoms (hydrogen and non-hydrogen) with undefined coordinates.
   DO I = 1,NATOMS
      QUNDEFINED(I) = ANY ( ATMCRD(1:3,I) == UNDEFINED )
   END DO
   NUNDEFH = COUNT ( QUNDEFINED .AND. ( ATMNUM == 1 ) )
   NUNDEFA = COUNT ( QUNDEFINED ) - NUNDEFH

   ! . There are hydrogens that need building.
   IF ( NUNDEFH > 0 ) THEN

      ! . Set up the connectivity data structures.
      CALL GET_CENTERS

      ! . Use the connectivity data to build the coordinates.
      CALL GET_COORDINATES

      ! . Deallocate temporary space.
      DEALLOCATE ( CENTERS )

   END IF

   ! . Write out data about the building.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF0000", VARIABLEWIDTH = 12 )
   CALL PRINT_SUMMARY_START ( "Summary of Hydrogen-Building Data" )
   WRITE ( PRINT_LINE, "(I12)" ) NUNDEFH ; CALL PRINT_SUMMARY_ELEMENT ( "No. Undefined Hydrogens" )
   WRITE ( PRINT_LINE, "(I12)" ) NUNDEFA ; CALL PRINT_SUMMARY_ELEMENT ( "No. Undefined Non-Hs"    )

   ! . If there were some hydrogens were built write out more data.
   IF ( NUNDEFH > 0 ) THEN
      NUNDEFH2 = COUNT ( QUNDEFINED .AND. ( ATMNUM == 1 ) )
      WRITE ( PRINT_LINE, "(I12)" ) NUNDEFH - NUNDEFH2 ; CALL PRINT_SUMMARY_ELEMENT ( "No. Hydrogens Built"   )
      WRITE ( PRINT_LINE, "(I12)" ) NCENTERS           ; CALL PRINT_SUMMARY_ELEMENT ( "No. Centers"           )
      WRITE ( PRINT_LINE, "(I12)" ) NUNDEFH2           ; CALL PRINT_SUMMARY_ELEMENT ( "No. Hydrogens Unbuilt" )
   END IF

   ! . Write out the terminator.
   CALL PRINT_SUMMARY_STOP

   !===========================================================================
   CONTAINS
   !===========================================================================

      !---------------------
      SUBROUTINE GET_CENTERS
      !---------------------

      ! . Local type definitions.
      TYPE CBOND_TYPE
         INTEGER            :: CENTER, J
         LOGICAL            :: QUNDEFINED
         REAL ( KIND = DP ) :: EQ
      END TYPE CBOND_TYPE

      ! . Local scalars.
      INTEGER :: I, IANGLE, IBOND, ICENTER, J, JATOM, K, N, NCBONDS
      LOGICAL :: QI, QJ

      ! . Local arrays.
      INTEGER,          DIMENSION(1:NATOMS)   :: CINDEX, CNDEFINED, CNUNDEFINED, HVALENCY
      LOGICAL,          DIMENSION(1:NATOMS)   :: QOK
      TYPE(CBOND_TYPE), DIMENSION(1:2*NBONDS) :: CBONDS

      !------------------------------------------------------------------------
      ! . Initialization.
      !------------------------------------------------------------------------
      ! . Check that a BONDS array exists.
      IF ( ( NBONDS <= 0 ) .AND. ( .NOT. ALLOCATED ( BONDS ) ) ) THEN
         CALL PRINT_ERROR ( "BUILD_HYDROGENS", "Cannot build hydrogen coordinates without appropriate MM definitions." )
      END IF

      !------------------------------------------------------------------------
      ! . Find some center and valency data for the undefined hydrogens.
      !------------------------------------------------------------------------
      ! . Initialization.
      CNUNDEFINED = 0
      HVALENCY    = 0

      ! . Loop over the bonds to find the heavy atom centers to which the undefined hydrogens are bound.
      DO IBOND = 1,NBONDS

         ! . Initialization.
         ICENTER = 0

         ! . Get the atom indices.
         I = BONDS(IBOND)%I
         J = BONDS(IBOND)%J

         ! . Are the atoms in the bond undefined hydrogens?
         QI = QUNDEFINED(I) .AND. ( ATMNUM(I) == 1 )
         QJ = QUNDEFINED(J) .AND. ( ATMNUM(J) == 1 )

         ! . Two undefined hydrogens.
         IF ( QI .AND. QJ ) THEN
            CONTINUE
         ! . One undefined hydrogen.
         ELSE IF ( QI ) THEN
            IF ( ATMNUM(J) /= 1 ) THEN
               JATOM   = I
               ICENTER = J
            END IF
         ELSE IF ( QJ ) THEN
            IF ( ATMNUM(I) /= 1 ) THEN
               JATOM   = J
               ICENTER = I
            END IF
         END IF

         ! . A center has been found so save some data about the atoms.
         IF ( ICENTER > 0 ) THEN
            CNUNDEFINED(ICENTER) = CNUNDEFINED(ICENTER) + 1
            HVALENCY(JATOM)      = HVALENCY(JATOM)      + 1
         END IF
      END DO

      !------------------------------------------------------------------------
      ! . Create a reduced bond list for the undefined-hydrogens' centers.
      !------------------------------------------------------------------------
      ! . Initialization.
      CNDEFINED = 0
      NCBONDS   = 0

      ! . Initialize QOK to be all centers with undefined hydrogens.
      QOK = ( CNUNDEFINED > 0 )

      ! . Loop over the bond list again to create a reduced bond list.
      DO IBOND = 1,NBONDS

         ! . Get the atom indices.
         I = BONDS(IBOND)%I
         J = BONDS(IBOND)%J

         ! . The first atom in the bond is a valid center.
         IF ( QOK(I) ) THEN

            ! . Save some bond data.
            NCBONDS = NCBONDS + 1
            CBONDS(NCBONDS)%CENTER     = I
            CBONDS(NCBONDS)%J          = J
            CBONDS(NCBONDS)%QUNDEFINED = .FALSE.
            CBONDS(NCBONDS)%EQ         = BONDS(IBOND)%EQ

         END IF

         ! . The second atom in the bond is a valid center.
         IF ( QOK(J) ) THEN

            ! . Save some bond data.
            NCBONDS = NCBONDS + 1
            CBONDS(NCBONDS)%CENTER     = J
            CBONDS(NCBONDS)%J          = I
            CBONDS(NCBONDS)%QUNDEFINED = .FALSE.
            CBONDS(NCBONDS)%EQ         = BONDS(IBOND)%EQ

         END IF

      END DO

      ! . Loop over the reduced bond list to find centers which are unsuitable.
      DO IBOND = 1,NCBONDS

         ! . Get the center and second atom indices.
         ICENTER = CBONDS(IBOND)%CENTER
         JATOM   = CBONDS(IBOND)%J

         ! . The center has undefined coordinates.
         IF ( QUNDEFINED(ICENTER) ) THEN
            QOK(ICENTER) = .FALSE.

         ! . The non-hydrogen atom to which the center is bound has undefined coordinates.
         ELSE IF ( ( ATMNUM(JATOM) /= 1 ) .AND. QUNDEFINED(JATOM) ) THEN
            QOK(ICENTER) = .FALSE.

         ! . The second atom is an undefined hydrogen atom.
         ELSE IF ( ( ATMNUM(JATOM) == 1 ) .AND. QUNDEFINED(JATOM) ) THEN

            ! . The hydrogen atom is bridging.
            IF ( HVALENCY(JATOM) /= 1 ) THEN
               QOK(ICENTER) = .FALSE.

            ! . Set the QUNDEFINED flag.
            ELSE
               CBONDS(IBOND)%QUNDEFINED = .TRUE.
            END IF

         ! . The second atom is defined.
         ELSE
            CNDEFINED(ICENTER) = CNDEFINED(ICENTER) + 1
         END IF

      END DO

      ! . Take out from the center list all those with a valency that is too high.
      QOK = QOK .AND. ( ( CNDEFINED + CNUNDEFINED ) <= MAX_VALENCY )

      !------------------------------------------------------------------------
      ! . Reorganize the center data into a more usable form.
      !------------------------------------------------------------------------
      ! . Determine the number of centers.
      NCENTERS = COUNT ( QOK )

      ! . Allocate the centers data structure.
      ALLOCATE ( CENTERS(1:NCENTERS) )

      ! . Initialization.
      CINDEX  = -1
      ICENTER =  0

      ! . Loop over the atoms.
      DO I = 1,NATOMS

         ! . Check for a valid center.
         IF ( QOK(I) ) THEN

            ! . Increment the center counter.
            ICENTER = ICENTER + 1

            ! . Set the index of the center and the number of defined and undefined atoms.
            CENTERS(ICENTER)%CENTER     = I
            CENTERS(ICENTER)%NDEFINED   = CNDEFINED(I)
            CENTERS(ICENTER)%NUNDEFINED = CNUNDEFINED(I)

            ! . Initialize the center type and angle depending upon the valency of the center.
            SELECT CASE ( CNDEFINED(I) + CNUNDEFINED(I) )
            CASE ( 2 )   ; CENTERS(ICENTER)%ANGLE = LINEAR          ; CENTERS(ICENTER)%TYPE = "LINEAR          "
            CASE ( 3 )   ; CENTERS(ICENTER)%ANGLE = TRIGONAL_PLANAR ; CENTERS(ICENTER)%TYPE = "PLANAR          "
            CASE ( 4 )   ; CENTERS(ICENTER)%ANGLE = TETRAHEDRAL     ; CENTERS(ICENTER)%TYPE = "TETRAHEDRAL     "
            CASE DEFAULT ; CENTERS(ICENTER)%ANGLE = 0.0_DP          ; CENTERS(ICENTER)%TYPE = "                "
            END SELECT

            ! . Set the index into the atom array.
            CINDEX(I) = ICENTER

         END IF
      END DO

      ! . Reinitialize CNDEFINED and CNUNDEFINED.
      CNDEFINED = 0 ; CNUNDEFINED = 0

      ! . Loop over the reduced bond list.
      DO IBOND = 1,NCBONDS

         ! . Get the center index.
         I = CBONDS(IBOND)%CENTER

         ! . The center is valid.
         IF ( QOK(I) ) THEN

            ! . Get the second atom index and the center index.
            J       = CBONDS(IBOND)%J
            ICENTER = CINDEX(I)

            ! . Get the position of the new entry (defined atoms first, undefined afterwards).
            IF ( CBONDS(IBOND)%QUNDEFINED ) THEN
               CNUNDEFINED(I) = CNUNDEFINED(I) + 1
               N              = CENTERS(ICENTER)%NDEFINED + CNUNDEFINED(I)
            ELSE
               CNDEFINED(I)   = CNDEFINED(I)   + 1
               N              = CNDEFINED(I)
            END IF

            ! . Put the bond in the appropriate place.
            CENTERS(ICENTER)%ATOMS(N) = J
            CENTERS(ICENTER)%BONDS(N) = CBONDS(IBOND)%EQ

         END IF

      END DO

      !------------------------------------------------------------------------
      ! . For centers with a valency of two check ANGLE.
      !------------------------------------------------------------------------
      ! . Check that an ANGLES array exists.
      IF ( ( NANGLES > 0 ) .AND. ( ALLOCATED ( ANGLES ) ) ) THEN

         ! . Set QOK to only those centers with a valency of two.
         QOK = QOK .AND. ( ( CNDEFINED + CNUNDEFINED ) == 2 )

         ! . Loop over the angles.
         DO IANGLE = 1,NANGLES

            ! . Get the central atom of the angle.
            J = ANGLES(IANGLE)%J

            ! . The central atom of the angle is a valid center.
            IF ( QOK ( J ) ) THEN

               ! . Get some information for the center.
               ICENTER = CINDEX(J)
               I       = CENTERS(ICENTER)%ATOMS(1)
               K       = CENTERS(ICENTER)%ATOMS(2)

               ! . The other two atoms of the angle match those of the center.
               IF ( ( ( I == ANGLES(IANGLE)%I ) .AND. ( K == ANGLES(IANGLE)%K ) ) .OR. &
                    ( ( I == ANGLES(IANGLE)%K ) .AND. ( K == ANGLES(IANGLE)%I ) ) ) THEN

                  ! . Reset the angle definition terms for the center.
                  CENTERS(ICENTER)%ANGLE = ANGLES(IANGLE)%EQ
                  CENTERS(ICENTER)%TYPE  = "                "

               END IF
            END IF
         END DO

      END IF

      END SUBROUTINE GET_CENTERS

      !-------------------------
      SUBROUTINE GET_COORDINATES
      !-------------------------

      ! . Local scalars.
      INTEGER :: C, I, ICENTER, J, K, L, NDEF, NUND

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: RA, RB, RC, RIC, RJC, RKC

      ! . Loop over the centers for which hydrogens are to be built.
      DO ICENTER = 1,NCENTERS

         ! . Get some information for the center.
         C    = CENTERS(ICENTER)%CENTER
         NDEF = CENTERS(ICENTER)%NDEFINED
         NUND = CENTERS(ICENTER)%NUNDEFINED

         ! . Get the atom indices.
         I = CENTERS(ICENTER)%ATOMS(1)
         J = CENTERS(ICENTER)%ATOMS(2)
         K = CENTERS(ICENTER)%ATOMS(3)
         L = CENTERS(ICENTER)%ATOMS(4)

         !---------------------------------------------------------------------
         ! . No atoms on the center are defined.
         !---------------------------------------------------------------------
         IF ( ( NDEF == 0 ) .AND. ( NUND > 0 ) ) THEN

            ! . Put the atom in a random direction from the center.
            ATMCRD(1:3,I) = ATMCRD(1:3,C) + CENTERS(ICENTER)%BONDS(1) * RANDOM_VECTOR ( 3, NORMALIZE = .TRUE. )

            ! . Reset some variables.
            NDEF          = 1
            NUND          = NUND - 1
            QUNDEFINED(I) = .FALSE.

         END IF

         !---------------------------------------------------------------------
         ! . One atom on the center is defined.
         !---------------------------------------------------------------------
         IF ( ( NDEF == 1 ) .AND. ( NUND > 0 ) ) THEN

            ! . Get the normalized CI vector.
            RIC = NORMALIZE ( ATMCRD(1:3,I) - ATMCRD(1:3,C) )

            ! . Get a random vector perpendicular to RIC.
            RA  = RANDOM_VECTOR ( 3, NORMALIZE = .TRUE., ORTHOGONAL_TO = RESHAPE ( RIC, (/ 3, 1 /) ) )

            ! . Put the atom in a random direction from the center.
            ATMCRD(1:3,J) = ATMCRD(1:3,C) + CENTERS(ICENTER)%BONDS(2) * ( COS ( CENTERS(ICENTER)%ANGLE ) * RIC + &
                                                                          SIN ( CENTERS(ICENTER)%ANGLE ) * RA    )

            ! . Reset some variables.
            NDEF          = 2
            NUND          = NUND - 1
            QUNDEFINED(J) = .FALSE.

         END IF

         !---------------------------------------------------------------------
         ! . Two atoms on the center are defined.
         !---------------------------------------------------------------------
         IF ( ( NDEF == 2 ) .AND. ( NUND > 0 ) ) THEN

            ! . Get the RIC and RJC vectors.
            RIC = NORMALIZE ( ATMCRD(1:3,I) - ATMCRD(1:3,C) )
            RJC = NORMALIZE ( ATMCRD(1:3,J) - ATMCRD(1:3,C) )

            ! . Get the bisector of RIC and RJC pointing in the opposite direction.
            RA = NORMALIZE ( RIC + RJC )

            ! . Branch on the center type.
            SELECT CASE ( CENTERS(ICENTER)%TYPE )

            ! . Add a planar atom.
            CASE ( "PLANAR          " ) ; ATMCRD(1:3,K) = ATMCRD(1:3,C) - CENTERS(ICENTER)%BONDS(3) * RA

            ! . Add a tetrahedral atom.
            CASE ( "TETRAHEDRAL     " )

            ! . Get the in-plane vector perpendicular to RC.
            RB = NORMALIZE ( RJC - RIC )

            ! . Get the out-of-plane vector perpendicular to RA and RB.
            RC = NORMALIZE ( CROSS_PRODUCT ( RA, RB ) )

            ! . Construct the coordinates of the tetrahedral atom.
            ATMCRD(1:3,K) = ATMCRD(1:3,C) + CENTERS(ICENTER)%BONDS(3) * ( COS ( TETRAHEDRAL ) * RA + &
                                                                          SIN ( TETRAHEDRAL ) * RC   )

            END SELECT

            ! . Reset some variables.
            NDEF          = 3
            NUND          = NUND - 1
            QUNDEFINED(K) = .FALSE.

         END IF

         !---------------------------------------------------------------------
         ! . Three atoms on the center are defined.
         !---------------------------------------------------------------------
         IF ( ( NDEF == 3 ) .AND. ( NUND > 0 ) ) THEN

            ! . Get the normalized RIC, RJC and RKC vectors.
            RIC = NORMALIZE ( ATMCRD(1:3,I) - ATMCRD(1:3,C) )
            RJC = NORMALIZE ( ATMCRD(1:3,J) - ATMCRD(1:3,C) )
            RKC = NORMALIZE ( ATMCRD(1:3,K) - ATMCRD(1:3,C) )

            ! . Get the vector from the mid-point of the (normalized) IJK plane to C.
            RA = - NORMALIZE ( RIC + RJC + RKC )

            ! . Construct the coordinates of the tetrahedral atom.
            ATMCRD(1:3,L) = ATMCRD(1:3,C) + CENTERS(ICENTER)%BONDS(4) * RA

            ! . Reset some variables.
            NDEF          = 4
            NUND          = NUND - 1
            QUNDEFINED(L) = .FALSE.

         END IF

      END DO

      END SUBROUTINE GET_COORDINATES

   END SUBROUTINE BUILD_HYDROGENS

END MODULE BUILD_COORDINATES
