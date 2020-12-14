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
!                         The Atoms Manipulation Module
!===============================================================================
!
! . Functions:
!
!   ATOM_NUMBER                 Get an atom's sequence number given its identity.
!   ATOM_SELECTION              Get an atom selection from sequence information.
!   ATOM_SELECTION_BY_DISTANCE  Get an atom selection from distance information.
!
!===============================================================================
MODULE ATOM_MANIPULATION

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_ERROR
USE STRING,      ONLY : TO_UPPER_CASE

USE ATOMS,       ONLY : ATMCRD, ATMNAM, NATOMS
USE SEQUENCE,    ONLY : NRESID, NSUBSYS, RESIND, RESNAM, SUBIND, SUBNAM

IMPLICIT NONE
PRIVATE
PUBLIC :: ATOM_NUMBER, ATOM_SELECTION, ATOM_SELECTION_BY_DISTANCE

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------------------------------
   FUNCTION ATOM_NUMBER ( ATOM_NAME, RESIDUE_NUMBER, SUBSYSTEM )
   !------------------------------------------------------------

   ! . Function declarations.
   INTEGER :: ATOM_NUMBER

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: ATOM_NAME, SUBSYSTEM
   INTEGER,               INTENT(IN) :: RESIDUE_NUMBER

   ! . Local scalars.
   INTEGER :: I, IATOM, IRESID, ISUBSYS

   ! . Initialize the function value.
   ATOM_NUMBER = -1

   ! . Check the number of atoms, residues or subsystems.
   IF ( ( NATOMS <= 0 ) .OR. ( NRESID <= 0 ) .OR. ( NSUBSYS <= 0 ) ) RETURN

   ! . Initialize some counters.
   IATOM   = -1
   ISUBSYS = -1

   ! . Get the subsystem number.
   DO I = 1,NSUBSYS
      IF ( SUBSYSTEM == SUBNAM(I) ) THEN
         ISUBSYS = I ; EXIT
      END IF
   END DO
   IF ( ISUBSYS <= 0 ) RETURN

   ! . Get the residue number in the subsystem.
   IRESID = RESIDUE_NUMBER + SUBIND(ISUBSYS)

   ! . Check the residue number.
   IF ( ( IRESID <= SUBIND(ISUBSYS) ) .OR. ( IRESID > SUBIND(ISUBSYS+1) ) ) RETURN

   ! . Check the atom names in the residue.
   DO I = (RESIND(IRESID)+1),RESIND(IRESID+1)
      IF ( ATOM_NAME == ATMNAM(I) ) THEN
         IATOM = I ; EXIT
      END IF
   END DO

   ! . Set the function value.
   ATOM_NUMBER = IATOM

   END FUNCTION ATOM_NUMBER

   !------------------------------------------------------------------------------------------
   FUNCTION ATOM_SELECTION ( ATOM_NAME, ATOM_NUMBER, RESIDUE_NAME, RESIDUE_NUMBER, SUBSYSTEM )
   !------------------------------------------------------------------------------------------

   ! . Function declarations.
   LOGICAL, DIMENSION(1:NATOMS) :: ATOM_SELECTION

   ! . Array arguments.
   CHARACTER ( LEN = * ), DIMENSION(:), INTENT(IN), OPTIONAL :: ATOM_NAME, RESIDUE_NAME, SUBSYSTEM
   INTEGER,               DIMENSION(:), INTENT(IN), OPTIONAL :: ATOM_NUMBER, RESIDUE_NUMBER

   ! . Local scalars.
   INTEGER :: I, J, N

   ! . Local arrays.
   LOGICAL, DIMENSION(1:NATOMS) :: QATOM_NAME, QATOM_NUMBER, QRESIDUE_NAME, QRESIDUE_NUMBER, QSUBSYSTEM

   ! . Initialize the function value.
   ATOM_SELECTION = .FALSE.

   ! . Check the number of atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Initialize the individual selection arrays.
   QATOM_NAME      = .TRUE.
   QATOM_NUMBER    = .TRUE.
   QRESIDUE_NAME   = .TRUE.
   QRESIDUE_NUMBER = .TRUE.
   QSUBSYSTEM      = .TRUE.

   ! . The ATOM_NAME argument is present.
   IF ( PRESENT ( ATOM_NAME ) ) THEN

      ! . Initialize the selection array.
      QATOM_NAME = .FALSE.

      ! . Loop over the atom names.
      DO I = 1,SIZE ( ATOM_NAME )
         ! . Loop over the atom names.
         DO J = 1,NATOMS
            IF ( ATMNAM(J) == ATOM_NAME(I) ) QATOM_NAME(J) = .TRUE.
         END DO
      END DO

   END IF

   ! . The ATOM_NUMBER argument is present.
   IF ( PRESENT ( ATOM_NUMBER ) ) THEN

      ! . Initialize the selection array.
      QATOM_NUMBER = .FALSE.

      ! . Loop over the atom numbers.
      DO I = 1,SIZE ( ATOM_NUMBER )
         IF ( ( ATOM_NUMBER(I) > 0 ) .AND. ( ATOM_NUMBER(I) <= NATOMS ) ) QATOM_NUMBER(ATOM_NUMBER(I)) = .TRUE.
      END DO

   END IF

   ! . The RESIDUE_NAME argument is present.
   IF ( PRESENT ( RESIDUE_NAME ) .AND. ( NRESID > 0 ) ) THEN

      ! . Initialize the selection array.
      QRESIDUE_NAME = .FALSE.

      ! . Loop over the residue names in the argument.
      DO I = 1,SIZE ( RESIDUE_NAME )
         ! . Loop over the residue names in the sequence.
         DO J = 1,NRESID
            IF ( RESNAM(J) == RESIDUE_NAME(I) ) QRESIDUE_NAME(RESIND(J)+1:RESIND(J+1)) = .TRUE.
         END DO
      END DO

   END IF

   ! . The RESIDUE_NUMBER argument is present.
   IF ( PRESENT ( RESIDUE_NUMBER ) .AND. ( NRESID > 0 ) .AND. ( NSUBSYS > 0 ) ) THEN

      ! . Initialize the selection array.
      QRESIDUE_NUMBER = .FALSE.

      ! . Loop over the residue numbers in the argument.
      DO I = 1,SIZE ( RESIDUE_NUMBER )

         ! . Loop over the subsystems.
         DO J = 1,NSUBSYS

            ! . Get the residue number in the sequence.
            N = RESIDUE_NUMBER(I) + SUBIND(J)

            ! . Check to see if a residue of this number exists.
            IF ( ( N > SUBIND(J) ) .AND. ( N <= SUBIND(J+1) ) ) QRESIDUE_NUMBER(RESIND(N)+1:RESIND(N+1)) = .TRUE.

         END DO

      END DO

   END IF

   ! . The SUBSYSTEM argument is present.
   IF ( PRESENT ( SUBSYSTEM ) .AND. ( NSUBSYS > 0 )  ) THEN

      ! . Initialize the selection array.
      QSUBSYSTEM = .FALSE.

      ! . Loop over the subsystem names in the argument.
      DO I = 1,SIZE ( SUBSYSTEM )
         ! . Loop over the subsystem names in the sequence.
         DO J = 1,NSUBSYS
            IF ( SUBNAM(J) == SUBSYSTEM(I) ) QSUBSYSTEM(RESIND(SUBIND(J)+1)+1:RESIND(SUBIND(J+1)+1)) = .TRUE.
         END DO
      END DO

   END IF

   ! . Set the function value.
   ATOM_SELECTION = QATOM_NAME .AND. QATOM_NUMBER .AND. QRESIDUE_NAME .AND. QRESIDUE_NUMBER .AND. QSUBSYSTEM

   END FUNCTION ATOM_SELECTION

   !----------------------------------------------------------------------
   FUNCTION ATOM_SELECTION_BY_DISTANCE ( TEST, CUTOFF, POINTS, INCLUSION )
   !----------------------------------------------------------------------

   ! . Function declarations.
   LOGICAL, DIMENSION(1:NATOMS) :: ATOM_SELECTION_BY_DISTANCE

   ! . Scalar arguments.
   CHARACTER ( LEN = 1 ), INTENT(IN) :: TEST
   REAL ( KIND = DP ),    INTENT(IN) :: CUTOFF

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: POINTS

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: INCLUSION

   ! . Local scalars.
   INTEGER            :: I, J
   LOGICAL            :: QLESSTHAN
   REAL ( KIND = DP ) :: RIJ2

   ! . Local arrays.
   LOGICAL, DIMENSION(1:NATOMS) :: QATOM

   ! . Initialize the function value.
   ATOM_SELECTION_BY_DISTANCE = .FALSE.

   ! . Check the number of atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Check the dimensions of POINTS.
   IF ( SIZE ( POINTS, 1 ) /= 3 ) CALL PRINT_ERROR ( "ATOM_SELECTION_BY_DISTANCE", "Invalid POINTS array." )

   ! . Check the TEST argument.
   IF ( ( TEST == ">" ) .OR. ( TEST == "+" ) ) THEN
      QLESSTHAN = .FALSE.
   ELSE IF ( ( TEST == "<" ) .OR. ( TEST == "-" ) ) THEN
      QLESSTHAN = .TRUE.
   ELSE
      CALL PRINT_ERROR ( "ATOM_SELECTION_BY_DISTANCE", "Invalid TEST argument." )
   END IF

   ! . Initialize QATOM.
   QATOM = .FALSE.

   ! . Loop over the points.
   DO I = 1,SIZE ( POINTS, 2 )

      ! . Loop over the atoms.
      DO J = 1,NATOMS

         ! . Calculate the distance.
         RIJ2 = SUM ( ( POINTS(1:3,I) - ATMCRD(1:3,J) )**2 )

         ! . Set QATOM.
         QATOM(J) = QATOM(J) .OR. ( RIJ2 < CUTOFF**2 )

      END DO
   END DO

   ! . Negate QATOM if the test is less than.
   IF ( .NOT. QLESSTHAN ) QATOM = .NOT. QATOM

   ! . The inclusion argument is present.
   IF ( PRESENT ( INCLUSION ) ) THEN

      ! . Branch on the option.
#ifdef PGPC
      ! . This is a PGPC compiler error.
      SELECT CASE ( INCLUSION )
#else
      SELECT CASE ( TO_UPPER_CASE ( INCLUSION ) )
#endif

      ! . Inclusion by residues.
      CASE ( "RESIDUE" )

         ! . Loop over the residues.
         DO I = 1,NRESID
            IF ( ANY ( QATOM(RESIND(I)+1:RESIND(I+1)) ) ) QATOM(RESIND(I)+1:RESIND(I+1)) = .TRUE.
         END DO

      ! . Inclusion by subsystems.
      CASE ( "SUBSYSTEM" )

         ! . Loop over the subsystems.
         DO I = 1,NSUBSYS
            IF ( ANY ( QATOM(RESIND(SUBIND(I)+1)+1:RESIND(SUBIND(I+1)+1)) ) ) THEN
               QATOM(RESIND(SUBIND(I)+1)+1:RESIND(SUBIND(I+1)+1)) = .TRUE.
            END IF
         END DO

      CASE DEFAULT ; CALL PRINT_ERROR ( "ATOM_SELECTION_BY_DISTANCE", "Unknown INCLUSION option." )
      END SELECT

   END IF

   ! . Set the function value.
   ATOM_SELECTION_BY_DISTANCE = QATOM

   END FUNCTION ATOM_SELECTION_BY_DISTANCE

END MODULE ATOM_MANIPULATION
