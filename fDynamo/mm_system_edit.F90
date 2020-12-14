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
!                   The Molecular Mechanics System Edit Module
!===============================================================================
!
! . Subroutines (public):
!
!   MM_SYSTEM_APPEND               Add a systems to an existing one.
!   MM_SYSTEM_DELETE               Delete some parts of a system.
!
! . Subroutines (private):
!
!   MM_SYSTEM_WRITE_SCRATCH        Write the MM system to a temporary file.
!
! . Notes:
!
!   Both subroutines leave unaltered the current symmetry information.
!
!===============================================================================
MODULE MM_SYSTEM_EDIT

! . Module declarations.
USE CONSTANTS,   ONLY : UNDEFINED
USE DEFINITIONS, ONLY : DP, FORCE_FIELD, MAX_RECORD_LENGTH, VERSION
USE FILES,       ONLY : NEXT_UNIT
USE PRINTING,    ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH

USE ATOMS
USE MM_TERMS
USE SEQUENCE
USE SYMMETRY,    ONLY : SYMMETRY_SUMMARY

IMPLICIT NONE
PRIVATE
PUBLIC :: MM_SYSTEM_APPEND, MM_SYSTEM_DELETE

!===============================================================================
CONTAINS
!===============================================================================

   !------------------------------------------
   SUBROUTINE MM_SYSTEM_APPEND ( FILE, PRINT )
   !------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Local scalars.
   CHARACTER ( LEN = LEN ( FORCE_FIELD ) ) :: HEADER
   INTEGER                                 :: I, IOSTAT, SCRATCH, UNIT
   LOGICAL                                 :: QPRINT
   REAL ( KIND = DP )                      :: VERNUM

   ! . New counters.
   INTEGER            :: NANG, NATM, NBND, NDIH, NEXCL, NEX14, NIMP, NRES, NSUB
   REAL ( KIND = DP ) :: SCELNEW, SCLJNEW

   ! . Old counters.
   INTEGER            :: NANGOLD, NATMOLD, NBNDOLD, NDIHOLD, NEXCOLD, NE14OLD, NIMPOLD, NRESOLD, NSUBOLD
   REAL ( KIND = DP ) :: SCELOLD, SCLJOLD

   ! . Get the next unit number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, ACTION = "READ", FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "OLD", IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "MM_SYSTEM_APPEND", "I/O Error in System File.", IOSTAT )

   ! . Header.
   READ ( UNIT ) HEADER, VERNUM
   IF ( ( HEADER /= FORCE_FIELD ) .OR. ( VERNUM > VERSION ) ) CALL PRINT_ERROR ( "MM_SYSTEM_APPEND", "Invalid system file header." )

   ! . Open a scratch file.
   SCRATCH = NEXT_UNIT ( )
   OPEN ( SCRATCH, FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "SCRATCH", IOSTAT = IOSTAT )
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "MM_SYSTEM_APPEND", "I/O Error in Scratch file.", IOSTAT )

   ! . Save the old numbers of atoms, residues and subsystems.
   NATMOLD = NATOMS ; NRESOLD = NRESID ; NSUBOLD = NSUBSYS

   ! . Save the old MM terms counters.
   NANGOLD = NANGLES ; NBNDOLD = NBONDS ; NDIHOLD = NDIHEDRALS ; NIMPOLD = NIMPROPERS
   NEXCOLD = ATMEXCI(NATOMS+1) ; NE14OLD = ATME14I(NATOMS+1)
   SCELOLD = SCALE_EL14 ; SCLJOLD = SCALE_LJ14

   ! . Write all the old data to the scratch file.
   CALL MM_SYSTEM_WRITE_SCRATCH ( SCRATCH, NRESOLD, NSUBOLD, NATMOLD, NANGOLD, NBNDOLD, NDIHOLD, NIMPOLD, &
                                           NEXCOLD, NE14OLD, SCELOLD, SCLJOLD )

   ! . Rewind the scratch file.
   REWIND ( SCRATCH )

   ! . Read in the number of residues and subsystems to append.
   READ ( UNIT ) NRES, NSUB

   ! . Initialize and reallocate the sequence data structure.
   CALL SEQUENCE_INITIALIZE
   CALL SEQUENCE_ALLOCATE ( ( NRES + NRESOLD ), ( NSUB + NSUBOLD ) )

   ! . Read in the old data from the scratch file.
   READ ( SCRATCH )
   DO I = 1,NRESOLD+1,MM_IO_BATCH ; READ ( SCRATCH ) RESIND(I:MIN(I+MM_IO_BATCH,NRESOLD+1))   ; ENDDO
   DO I = 1,NRESOLD,  MM_IO_BATCH ; READ ( SCRATCH ) RESNAM(I:MIN(I+MM_IO_BATCH,NRESOLD)  )   ; ENDDO
   DO I = 1,NSUBOLD+1,MM_IO_BATCH ; READ ( SCRATCH ) SUBIND(I:MIN(I+MM_IO_BATCH,NSUBOLD+1))   ; ENDDO
   DO I = 1,NSUBOLD,  MM_IO_BATCH ; READ ( SCRATCH ) SUBNAM(I:MIN(I+MM_IO_BATCH,NSUBOLD)  )   ; ENDDO

   ! . Read in the new data.
   READ ( UNIT ) RESIND(NRESOLD+1:NRESID+1)
   READ ( UNIT ) RESNAM(NRESOLD+1:NRESID)
   READ ( UNIT ) SUBIND(NSUBOLD+1:NSUBSYS+1)
   READ ( UNIT ) SUBNAM(NSUBOLD+1:NSUBSYS)

   ! . Increment the new elements of the residue and subsystem index arrays.
   RESIND(NRESOLD+1:NRESID+1)  = RESIND(NRESOLD+1:NRESID+1)  + NATMOLD
   SUBIND(NSUBOLD+1:NSUBSYS+1) = SUBIND(NSUBOLD+1:NSUBSYS+1) + NRESOLD

   ! . Read in the number of atoms to append.
   READ ( UNIT ) NATM

   ! . Initialize and reallocate the atom data structure.
   CALL ATOMS_INITIALIZE
   CALL ATOMS_ALLOCATE ( NATM + NATMOLD )

   ! . Read in the old data from the scratch file.
   READ ( SCRATCH )
   DO I = 1,NATMOLD,MM_IO_BATCH ; READ ( SCRATCH ) ATMCRD(1:3,I:MIN(I+MM_IO_BATCH,NATMOLD)) ; ENDDO
   DO I = 1,NATMOLD,MM_IO_BATCH ; READ ( SCRATCH ) ATMMAS(    I:MIN(I+MM_IO_BATCH,NATMOLD)) ; ENDDO
   DO I = 1,NATMOLD,MM_IO_BATCH ; READ ( SCRATCH ) ATMNAM(    I:MIN(I+MM_IO_BATCH,NATMOLD)) ; ENDDO
   DO I = 1,NATMOLD,MM_IO_BATCH ; READ ( SCRATCH ) ATMNUM(    I:MIN(I+MM_IO_BATCH,NATMOLD)) ; ENDDO

   ! . Read in the new data.
   ATMCRD(1:3,NATMOLD+1:NATOMS) = UNDEFINED
   READ ( UNIT ) ATMMAS(NATMOLD+1:NATOMS)
   READ ( UNIT ) ATMNAM(NATMOLD+1:NATOMS)
   READ ( UNIT ) ATMNUM(NATMOLD+1:NATOMS)

   ! . Read in the new MM terms counters.
   READ ( UNIT ) NANG, NBND, NDIH, NIMP, NEXCL, NEX14, SCELNEW, SCLJNEW

   ! . Print a warning if the scale factors are different.
   IF ( ( SCELOLD /= SCELNEW ) .OR. ( SCLJOLD /= SCLJNEW ) ) THEN
      CALL PRINT_PARAGRAPH ( TEXT = "Warning: Electrostatic and/or LJ scaling factors are different." )
   END IF

   ! . Initialize and reallocate the MM terms data structure.
   CALL MM_TERMS_INITIALIZE
   CALL MM_TERMS_ALLOCATE ( ( NBND + NBNDOLD ), ( NANG + NANGOLD ), ( NDIH + NDIHOLD ), ( NIMP + NIMPOLD ) )
   ALLOCATE ( ATMEXCJ(1:(NEXCL+NEXCOLD)), ATME14J(1:(NEX14+NE14OLD)) )

   ! . Reset the scaling factors.
   SCALE_EL14 = SCELOLD
   SCALE_LJ14 = SCLJOLD

   ! . Read in the old data from the scratch file.
   READ ( SCRATCH )
   DO I = 1,NATMOLD,  MM_IO_BATCH ; READ ( SCRATCH ) ATMCHG  (I:MIN(I+MM_IO_BATCH,NATMOLD))   ; ENDDO
   DO I = 1,NATMOLD,  MM_IO_BATCH ; READ ( SCRATCH ) ATMCHG14(I:MIN(I+MM_IO_BATCH,NATMOLD))   ; ENDDO
   DO I = 1,NATMOLD,  MM_IO_BATCH ; READ ( SCRATCH ) ATMEPS  (I:MIN(I+MM_IO_BATCH,NATMOLD))   ; ENDDO
   DO I = 1,NATMOLD,  MM_IO_BATCH ; READ ( SCRATCH ) ATMEPS14(I:MIN(I+MM_IO_BATCH,NATMOLD))   ; ENDDO
   DO I = 1,NATMOLD,  MM_IO_BATCH ; READ ( SCRATCH ) ATMSIG  (I:MIN(I+MM_IO_BATCH,NATMOLD))   ; ENDDO
   DO I = 1,NATMOLD,  MM_IO_BATCH ; READ ( SCRATCH ) ATMTYP  (I:MIN(I+MM_IO_BATCH,NATMOLD))   ; ENDDO

   DO I = 1,NATMOLD+1,MM_IO_BATCH ; READ ( SCRATCH ) ATMEXCI(I:MIN(I+MM_IO_BATCH,NATMOLD+1))  ; ENDDO
   DO I = 1,NEXCOLD,  MM_IO_BATCH ; READ ( SCRATCH ) ATMEXCJ(I:MIN(I+MM_IO_BATCH,NEXCOLD  ))  ; ENDDO
   DO I = 1,NATMOLD+1,MM_IO_BATCH ; READ ( SCRATCH ) ATME14I(I:MIN(I+MM_IO_BATCH,NATMOLD+1))  ; ENDDO
   DO I = 1,NE14OLD,  MM_IO_BATCH ; READ ( SCRATCH ) ATME14J(I:MIN(I+MM_IO_BATCH,NE14OLD  ))  ; ENDDO

   DO I = 1,NANGOLD,  MM_IO_BATCH ; READ ( SCRATCH )    ANGLES(I:MIN(I+MM_IO_BATCH,NANGOLD))  ; ENDDO
   DO I = 1,NBNDOLD,  MM_IO_BATCH ; READ ( SCRATCH )     BONDS(I:MIN(I+MM_IO_BATCH,NBNDOLD))  ; ENDDO
   DO I = 1,NDIHOLD,  MM_IO_BATCH ; READ ( SCRATCH ) DIHEDRALS(I:MIN(I+MM_IO_BATCH,NDIHOLD))  ; ENDDO
   DO I = 1,NIMPOLD,  MM_IO_BATCH ; READ ( SCRATCH ) IMPROPERS(I:MIN(I+MM_IO_BATCH,NIMPOLD))  ; ENDDO

   ! . Read in the new data.
   READ ( UNIT ) ATMCHG(NATMOLD+1:NATOMS)
   READ ( UNIT ) ATMCHG14(NATMOLD+1:NATOMS)
   READ ( UNIT ) ATMEPS(NATMOLD+1:NATOMS)
   READ ( UNIT ) ATMEPS14(NATMOLD+1:NATOMS)
   READ ( UNIT ) ATMSIG(NATMOLD+1:NATOMS)
   READ ( UNIT ) ATMTYP(NATMOLD+1:NATOMS)

   READ ( UNIT ) ATMEXCI(NATMOLD+1:NATOMS+1)
   READ ( UNIT ) ATMEXCJ(NEXCOLD+1:NEXCOLD+NEXCL)
   READ ( UNIT ) ATME14I(NATMOLD+1:NATOMS+1)
   READ ( UNIT ) ATME14J(NE14OLD+1:NE14OLD+NEX14)

   DO I = NANGOLD+1,NANGLES,    MM_IO_BATCH ; READ ( UNIT )    ANGLES(I:MIN(I+MM_IO_BATCH,NANGLES))    ; ENDDO
   DO I = NBNDOLD+1,NBONDS,     MM_IO_BATCH ; READ ( UNIT )     BONDS(I:MIN(I+MM_IO_BATCH,NBONDS))     ; ENDDO
   DO I = NDIHOLD+1,NDIHEDRALS, MM_IO_BATCH ; READ ( UNIT ) DIHEDRALS(I:MIN(I+MM_IO_BATCH,NDIHEDRALS)) ; ENDDO
   DO I = NIMPOLD+1,NIMPROPERS, MM_IO_BATCH ; READ ( UNIT ) IMPROPERS(I:MIN(I+MM_IO_BATCH,NIMPROPERS)) ; ENDDO

   ! . Increment the new elements of the exclusion arrays.
   ATMEXCI(NATMOLD+1:NATOMS+1) = ATMEXCI(NATMOLD+1:NATOMS+1) + NEXCOLD
   ATME14I(NATMOLD+1:NATOMS+1) = ATME14I(NATMOLD+1:NATOMS+1) + NE14OLD
   ATMEXCJ(NEXCOLD+1:NEXCOLD+NEXCL) = ATMEXCJ(NEXCOLD+1:NEXCOLD+NEXCL) + NATMOLD
   ATME14J(NE14OLD+1:NE14OLD+NEX14) = ATME14J(NE14OLD+1:NE14OLD+NEX14) + NATMOLD

   ! . Increment the atom indices in the angle, bond, dihedral and improper arrays.
   ANGLES(NANGOLD+1:NANGLES)%I = ANGLES(NANGOLD+1:NANGLES)%I + NATMOLD
   ANGLES(NANGOLD+1:NANGLES)%J = ANGLES(NANGOLD+1:NANGLES)%J + NATMOLD
   ANGLES(NANGOLD+1:NANGLES)%K = ANGLES(NANGOLD+1:NANGLES)%K + NATMOLD

   BONDS(NBNDOLD+1:NBONDS)%I = BONDS(NBNDOLD+1:NBONDS)%I + NATMOLD
   BONDS(NBNDOLD+1:NBONDS)%J = BONDS(NBNDOLD+1:NBONDS)%J + NATMOLD

   DIHEDRALS(NDIHOLD+1:NDIHEDRALS)%I = DIHEDRALS(NDIHOLD+1:NDIHEDRALS)%I + NATMOLD
   DIHEDRALS(NDIHOLD+1:NDIHEDRALS)%J = DIHEDRALS(NDIHOLD+1:NDIHEDRALS)%J + NATMOLD
   DIHEDRALS(NDIHOLD+1:NDIHEDRALS)%K = DIHEDRALS(NDIHOLD+1:NDIHEDRALS)%K + NATMOLD
   DIHEDRALS(NDIHOLD+1:NDIHEDRALS)%L = DIHEDRALS(NDIHOLD+1:NDIHEDRALS)%L + NATMOLD

   IMPROPERS(NIMPOLD+1:NIMPROPERS)%I = IMPROPERS(NIMPOLD+1:NIMPROPERS)%I + NATMOLD
   IMPROPERS(NIMPOLD+1:NIMPROPERS)%J = IMPROPERS(NIMPOLD+1:NIMPROPERS)%J + NATMOLD
   IMPROPERS(NIMPOLD+1:NIMPROPERS)%K = IMPROPERS(NIMPOLD+1:NIMPROPERS)%K + NATMOLD
   IMPROPERS(NIMPOLD+1:NIMPROPERS)%L = IMPROPERS(NIMPOLD+1:NIMPROPERS)%L + NATMOLD

   ! . Close the files.
   CLOSE ( SCRATCH )
   CLOSE ( UNIT    )

   ! . Write out some information.
   CALL PRINT_PARAGRAPH ( TEXT = "System data to append read from " // FILE )

   ! . Set the print flag.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Print out a summary of the system data.
   IF ( QPRINT ) THEN
      CALL    ATOMS_SUMMARY
      CALL MM_TERMS_SUMMARY
      CALL SEQUENCE_SUMMARY
      CALL SYMMETRY_SUMMARY
   END IF

   END SUBROUTINE MM_SYSTEM_APPEND

   !---------------------------------------------
   SUBROUTINE MM_SYSTEM_DELETE ( QSELECT, PRINT )
   !---------------------------------------------

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: QSELECT

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Local scalars.
   INTEGER :: IOSTAT, NDELETE, SCRATCH
   LOGICAL :: QPRINT

   ! . Counters.
   INTEGER :: NANG, NATM, NBND, NDIH, NEXC, NE14, NIMP, NRES, NSUB

   ! . Local arrays.
   INTEGER, DIMENSION(1:NATOMS) :: OLDTONEW

   ! . Check the number of atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Get the number of atoms to delete.
   NDELETE = COUNT ( QSELECT )

   ! . Return if there are no atoms to delete.
   IF ( NDELETE <= 0 ) RETURN

   ! . Contract the atom data.
!PGF . The PGF compiler may require OLDTONEW to be initialized here rather than in the subroutine.
   CALL CONTRACT_ATOM_DATA

   ! . Contract the sequence data structures.
   CALL CONTRACT_SEQUENCE_DATA

   ! . Contract the exclusion lists.
   CALL CONTRACT_EXCLUSION_LIST ( ATMEXCI, ATMEXCJ, NEXC )
   CALL CONTRACT_EXCLUSION_LIST ( ATME14I, ATME14J, NE14 )

   ! . Contract the angle, bond, dihedral and improper terms.
   CALL CONTRACT_MM_TERMS

   ! . Open a scratch file.
   SCRATCH = NEXT_UNIT ( )
   OPEN ( SCRATCH, FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "SCRATCH", IOSTAT = IOSTAT )
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "MM_SYSTEM_DELETE", "I/O Error in Scratch file.", IOSTAT )

   ! . Save the contracted data on the scratch file.
   CALL MM_SYSTEM_WRITE_SCRATCH ( SCRATCH, NRES, NSUB, NATM, NANG, NBND, NDIH, NIMP, NEXC, NE14, SCALE_EL14, SCALE_LJ14 )

   ! . Initialize the data structures.
   CALL    ATOMS_INITIALIZE
   CALL MM_TERMS_INITIALIZE
   CALL SEQUENCE_INITIALIZE

   ! . Rewind the scratch file.
   REWIND ( SCRATCH )

   ! . Read in the contracted data.
   CALL MM_SYSTEM_READ_SCRATCH ( SCRATCH )

   ! . Close the scratch file.
   CLOSE ( SCRATCH )

   ! . Write out some data.
   WRITE ( PRINT_LINE, "(A,I6,A)" ) "Number of atoms deleted from system file = ", NDELETE, "."
   CALL PRINT_PARAGRAPH

   ! . Set the print flag.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Print out a summary of the system data.
   IF ( QPRINT ) THEN
      CALL    ATOMS_SUMMARY
      CALL MM_TERMS_SUMMARY
      CALL SEQUENCE_SUMMARY
      CALL SYMMETRY_SUMMARY
   END IF

   !===============================================================================
   CONTAINS
   !===============================================================================

      !----------------------------
      SUBROUTINE CONTRACT_ATOM_DATA
      !----------------------------

      ! . Local scalars.
      INTEGER :: I, N

      ! . Initialize the atom index array.
      OLDTONEW = 0

      ! . Loop over the atoms.
      N = 0
      DO I = 1,NATOMS

         ! . The atom is to be saved.
         IF ( .NOT. QSELECT(I) ) THEN

            ! . Increment the atom number.
            N = N + 1

            ! . Save the atom index.
            OLDTONEW(I) = N

            ! . Keep the atom data.
            ATMCRD(1:3,N) = ATMCRD(1:3,I)
            ATMMAS(N)     = ATMMAS(I)
            ATMNAM(N)     = ATMNAM(I)
            ATMNUM(N)     = ATMNUM(I)

            ! . Keep the MM terms data.
            ATMCHG(N)   = ATMCHG(I)
            ATMCHG14(N) = ATMCHG14(I)
            ATMEPS(N)   = ATMEPS(I)
            ATMEPS14(N) = ATMEPS14(I)
            ATMSIG(N)   = ATMSIG(I)
            ATMTYP(N)   = ATMTYP(I)

         END IF
      END DO
      NATM = N

      END SUBROUTINE CONTRACT_ATOM_DATA

      !-------------------------------------------------------
      SUBROUTINE CONTRACT_EXCLUSION_LIST ( LISTI, LISTJ, NEX )
      !-------------------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(OUT) :: NEX

      ! . Array arguments.
      INTEGER, DIMENSION(:), INTENT(INOUT) :: LISTI, LISTJ

      ! . Local scalars.
      INTEGER :: I, IINT, J, N, START, STOP

      ! . Loop over the exclusions.
      N = 0
      DO I = 1,NATOMS

         ! . The atom is to be saved.
         IF ( .NOT. QSELECT(I) ) THEN

            ! . Get the starting and stopping points in the list.
            START = LISTI(I)+1
            STOP  = LISTI(I+1)

            ! . Save the new index.
            LISTI(OLDTONEW(I)) = N

            ! . Loop over the exclusions in the old list.
            DO IINT = START,STOP

               ! . Get the second atom in the list.
               J = LISTJ(IINT)

               ! . The exclusion is to be saved.
               IF ( .NOT. QSELECT(J) ) THEN

                  ! . Increment the exclusion count.
                  N = N + 1

                  ! . Save the exclusion.
                  LISTJ(N) = OLDTONEW(J)
  
               END IF
            END DO
         END IF
      END DO
      LISTI(NATM+1) = N

      ! . Save the exclusion counter.
      NEX = N

      END SUBROUTINE CONTRACT_EXCLUSION_LIST

      !---------------------------
      SUBROUTINE CONTRACT_MM_TERMS
      !---------------------------

      ! . Local scalars.
      INTEGER :: I, N

      ! . Loop over the angles.
      N = 0
      DO I = 1,NANGLES

         ! . Skip any angle that has a deleted atom.
         IF ( QSELECT(ANGLES(I)%I) .OR. QSELECT(ANGLES(I)%J) .OR. QSELECT(ANGLES(I)%K) ) CYCLE

         ! . Increment the angle number.
         N = N + 1

         ! . Keep the angle.
         ANGLES(N) = ANGLES(I)

         ! . Alter the atom indices for the angle.
         ANGLES(N)%I = OLDTONEW(ANGLES(I)%I)
         ANGLES(N)%J = OLDTONEW(ANGLES(I)%J)
         ANGLES(N)%K = OLDTONEW(ANGLES(I)%K)

      END DO
      NANG = N

      ! . Loop over the bonds.
      N = 0
      DO I = 1,NBONDS

         ! . Skip any bond that has a deleted atom.
         IF ( QSELECT(BONDS(I)%I) .OR. QSELECT(BONDS(I)%J) ) CYCLE

         ! . Increment the bond number.
         N = N + 1

         ! . Keep the bond.
         BONDS(N) = BONDS(I)

         ! . Alter the atom indices for the bond.
         BONDS(N)%I = OLDTONEW(BONDS(I)%I)
         BONDS(N)%J = OLDTONEW(BONDS(I)%J)

      END DO
      NBND = N

      ! . Loop over the dihedrals.
      N = 0
      DO I = 1,NDIHEDRALS

         ! . Skip any dihedral that has a deleted atom.
         IF ( QSELECT(DIHEDRALS(I)%I) .OR. QSELECT(DIHEDRALS(I)%J) .OR. QSELECT(DIHEDRALS(I)%K) .OR. QSELECT(DIHEDRALS(I)%L) ) CYCLE

         ! . Increment the dihedral number.
         N = N + 1

         ! . Keep the dihedral.
         DIHEDRALS(N) = DIHEDRALS(I)

         ! . Alter the atom indices for the dihedral.
         DIHEDRALS(N)%I = OLDTONEW(DIHEDRALS(I)%I)
         DIHEDRALS(N)%J = OLDTONEW(DIHEDRALS(I)%J)
         DIHEDRALS(N)%K = OLDTONEW(DIHEDRALS(I)%K)
         DIHEDRALS(N)%L = OLDTONEW(DIHEDRALS(I)%L)

      END DO
      NDIH = N

      ! . Loop over the impropers.
      N = 0
      DO I = 1,NIMPROPERS

         ! . Skip any improper that has a deleted atom.
         IF ( QSELECT(IMPROPERS(I)%I) .OR. QSELECT(IMPROPERS(I)%J) .OR. QSELECT(IMPROPERS(I)%K) .OR. QSELECT(IMPROPERS(I)%L) ) CYCLE

         ! . Increment the improper number.
         N = N + 1

         ! . Keep the improper.
         IMPROPERS(N) = IMPROPERS(I)

         ! . Alter the atom indices for the improper.
         IMPROPERS(N)%I = OLDTONEW(IMPROPERS(I)%I)
         IMPROPERS(N)%J = OLDTONEW(IMPROPERS(I)%J)
         IMPROPERS(N)%K = OLDTONEW(IMPROPERS(I)%K)
         IMPROPERS(N)%L = OLDTONEW(IMPROPERS(I)%L)

      END DO
      NIMP = N

      END SUBROUTINE CONTRACT_MM_TERMS

      !--------------------------------
      SUBROUTINE CONTRACT_SEQUENCE_DATA
      !--------------------------------

      ! . Local scalars.
      INTEGER :: ATMLFT, I, IRES, ISUB, N, NRESOLD, START, STOP

      ! . Local arrays.
      INTEGER, DIMENSION(1:NRESID+1)  :: RESINDNEW
      INTEGER, DIMENSION(1:NSUBSYS+1) :: SUBINDNEW

      ! . Initialize the new index arrays.
      RESINDNEW = 0
      SUBINDNEW = 0

      ! . Initialize some counters.
      N    = 0
      NRES = 0
      NSUB = 0

      ! . Loop over the subsystems.
      DO ISUB = 1,NSUBSYS

         ! . Save NRES.
         NRESOLD = NRES

         ! . Get the residue indices for the subsystem.
         START = SUBIND(ISUB)+1
         STOP  = SUBIND(ISUB+1)

         ! . Loop over the residues.
         DO IRES = START,STOP

            ! . Get the number of undeleted atoms in the residue.
	    ATMLFT = 0
            DO I = (RESIND(IRES)+1),RESIND(IRES+1)
               IF ( .NOT. QSELECT(I) ) ATMLFT = ATMLFT + 1
            END DO

            ! . Include the residue if there are atoms left.
            IF ( ATMLFT > 0 ) THEN

               ! . Increment the residue number.
               NRES = NRES + 1

               ! . Save the residue index and name.
               RESINDNEW(NRES) = N
               RESNAM(NRES)    = RESNAM(IRES)

               ! . Increment the atom number.
               N = N + ATMLFT

            END IF
         END DO

         ! . Include the subsystem if there are residues left.
         IF ( NRES - NRESOLD >  0 ) THEN

            ! . Increment the subsystem number.
            NSUB = NSUB + 1

            ! . Save the subsystem index and name.
            SUBINDNEW(NSUB) = NRESOLD
            SUBNAM(NSUB)    = SUBNAM(ISUB)

         END IF
      END DO

      ! . Fill the last elements of RESIND and SUBIND.
      RESINDNEW(NRES+1) = N
      SUBINDNEW(NSUB+1) = NRES

      ! . Copy the new to the old index arrays.
      RESIND = RESINDNEW
      SUBIND = SUBINDNEW

      END SUBROUTINE CONTRACT_SEQUENCE_DATA

   END SUBROUTINE MM_SYSTEM_DELETE

   !----------------------------------------------------------------------------
   ! . Private I/O Subroutines.
   !----------------------------------------------------------------------------

   !-----------------------------------------
   SUBROUTINE MM_SYSTEM_READ_SCRATCH ( UNIT )
   !-----------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: UNIT

   ! . Local scalars.
   INTEGER :: I, NRES, NSUB, NATM, NANG, NBND, NDIH, NIMP, NEXC, NE14

   ! . Residue and subsystem data.
   READ ( UNIT ) NRES, NSUB
   CALL SEQUENCE_ALLOCATE ( NRES, NSUB )
   DO I = 1,NRES+1,MM_IO_BATCH ; READ ( UNIT ) RESIND(I:MIN(I+MM_IO_BATCH,NRES+1))   ; ENDDO
   DO I = 1,NRES,  MM_IO_BATCH ; READ ( UNIT ) RESNAM(I:MIN(I+MM_IO_BATCH,NRES)  )   ; ENDDO
   DO I = 1,NSUB+1,MM_IO_BATCH ; READ ( UNIT ) SUBIND(I:MIN(I+MM_IO_BATCH,NSUB+1))   ; ENDDO
   DO I = 1,NSUB,  MM_IO_BATCH ; READ ( UNIT ) SUBNAM(I:MIN(I+MM_IO_BATCH,NSUB)  )   ; ENDDO

   ! . Atom data.
   READ ( UNIT ) NATM
   CALL ATOMS_ALLOCATE ( NATM )
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMCRD(1:3,I:MIN(I+MM_IO_BATCH,NATM)) ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMMAS(    I:MIN(I+MM_IO_BATCH,NATM)) ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMNAM(    I:MIN(I+MM_IO_BATCH,NATM)) ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMNUM(    I:MIN(I+MM_IO_BATCH,NATM)) ; ENDDO

   ! . MM terms data.
   READ ( UNIT ) NANG, NBND, NDIH, NIMP, NEXC, NE14, SCALE_EL14, SCALE_LJ14
   CALL MM_TERMS_ALLOCATE ( NBND, NANG, NDIH, NIMP )
   ALLOCATE ( ATMEXCJ(1:NEXC), ATME14J(1:NE14) )
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMCHG  (I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMCHG14(I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMEPS  (I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMEPS14(I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMSIG  (I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; READ ( UNIT ) ATMTYP  (I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO

   DO I = 1,NATM+1,MM_IO_BATCH ; READ ( UNIT ) ATMEXCI(I:MIN(I+MM_IO_BATCH,NATM+1))  ; ENDDO
   DO I = 1,NEXC,  MM_IO_BATCH ; READ ( UNIT ) ATMEXCJ(I:MIN(I+MM_IO_BATCH,NEXC  ))  ; ENDDO
   DO I = 1,NATM+1,MM_IO_BATCH ; READ ( UNIT ) ATME14I(I:MIN(I+MM_IO_BATCH,NATM+1))  ; ENDDO
   DO I = 1,NE14,  MM_IO_BATCH ; READ ( UNIT ) ATME14J(I:MIN(I+MM_IO_BATCH,NE14  ))  ; ENDDO

   DO I = 1,NANG,  MM_IO_BATCH ; READ ( UNIT )    ANGLES(I:MIN(I+MM_IO_BATCH,NANG))  ; ENDDO
   DO I = 1,NBND,  MM_IO_BATCH ; READ ( UNIT )     BONDS(I:MIN(I+MM_IO_BATCH,NBND))  ; ENDDO
   DO I = 1,NDIH,  MM_IO_BATCH ; READ ( UNIT ) DIHEDRALS(I:MIN(I+MM_IO_BATCH,NDIH))  ; ENDDO
   DO I = 1,NIMP,  MM_IO_BATCH ; READ ( UNIT ) IMPROPERS(I:MIN(I+MM_IO_BATCH,NIMP))  ; ENDDO

   END SUBROUTINE MM_SYSTEM_READ_SCRATCH

   !------------------------------------------------------------------------------------------------------------
   SUBROUTINE MM_SYSTEM_WRITE_SCRATCH ( UNIT, NRES, NSUB, NATM, NANG, NBND, NDIH, NIMP, NEXC, NE14, SCEL, SCLJ )
   !------------------------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: UNIT, NRES, NSUB, NATM, NANG, NBND, NDIH, NIMP, NEXC, NE14
   REAL ( KIND = DP ), INTENT(IN) :: SCEL, SCLJ

   ! . Local scalars.
   INTEGER :: I

   ! . Residue and subsystem data.
   WRITE ( UNIT ) NRES, NSUB
   DO I = 1,NRES+1,MM_IO_BATCH ; WRITE ( UNIT ) RESIND(I:MIN(I+MM_IO_BATCH,NRES+1))   ; ENDDO
   DO I = 1,NRES,  MM_IO_BATCH ; WRITE ( UNIT ) RESNAM(I:MIN(I+MM_IO_BATCH,NRES)  )   ; ENDDO
   DO I = 1,NSUB+1,MM_IO_BATCH ; WRITE ( UNIT ) SUBIND(I:MIN(I+MM_IO_BATCH,NSUB+1))   ; ENDDO
   DO I = 1,NSUB,  MM_IO_BATCH ; WRITE ( UNIT ) SUBNAM(I:MIN(I+MM_IO_BATCH,NSUB)  )   ; ENDDO

   ! . Atom data.
   WRITE ( UNIT ) NATM
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMCRD(1:3,I:MIN(I+MM_IO_BATCH,NATM)) ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMMAS(    I:MIN(I+MM_IO_BATCH,NATM)) ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMNAM(    I:MIN(I+MM_IO_BATCH,NATM)) ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMNUM(    I:MIN(I+MM_IO_BATCH,NATM)) ; ENDDO

   ! . MM terms data.
   WRITE ( UNIT ) NANG, NBND, NDIH, NIMP, NEXC, NE14, SCEL, SCLJ
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMCHG  (I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMCHG14(I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMEPS  (I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMEPS14(I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMSIG  (I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO
   DO I = 1,NATM,  MM_IO_BATCH ; WRITE ( UNIT ) ATMTYP  (I:MIN(I+MM_IO_BATCH,NATM))   ; ENDDO

   DO I = 1,NATM+1,MM_IO_BATCH ; WRITE ( UNIT ) ATMEXCI(I:MIN(I+MM_IO_BATCH,NATM+1))  ; ENDDO
   DO I = 1,NEXC,  MM_IO_BATCH ; WRITE ( UNIT ) ATMEXCJ(I:MIN(I+MM_IO_BATCH,NEXC  ))  ; ENDDO
   DO I = 1,NATM+1,MM_IO_BATCH ; WRITE ( UNIT ) ATME14I(I:MIN(I+MM_IO_BATCH,NATM+1))  ; ENDDO
   DO I = 1,NE14,  MM_IO_BATCH ; WRITE ( UNIT ) ATME14J(I:MIN(I+MM_IO_BATCH,NE14  ))  ; ENDDO

   DO I = 1,NANG,  MM_IO_BATCH ; WRITE ( UNIT )    ANGLES(I:MIN(I+MM_IO_BATCH,NANG))  ; ENDDO
   DO I = 1,NBND,  MM_IO_BATCH ; WRITE ( UNIT )     BONDS(I:MIN(I+MM_IO_BATCH,NBND))  ; ENDDO
   DO I = 1,NDIH,  MM_IO_BATCH ; WRITE ( UNIT ) DIHEDRALS(I:MIN(I+MM_IO_BATCH,NDIH))  ; ENDDO
   DO I = 1,NIMP,  MM_IO_BATCH ; WRITE ( UNIT ) IMPROPERS(I:MIN(I+MM_IO_BATCH,NIMP))  ; ENDDO

   END SUBROUTINE MM_SYSTEM_WRITE_SCRATCH

END MODULE MM_SYSTEM_EDIT










