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
!                            The MOPAC Density Module
!===============================================================================
!
! . Subroutines.
!
!   DENSITY_CALCULATE        calculate the density matrix from the MOs.
!   DENSITY_GUESS            guess an initial density matrix.
!   DENSITY_READ             read the density matrix from a file.
!   DENSITY_WRITE            write the density matrix to a file.
!
! . Notes:
!
!   The algorithm in DENSITY_GUESS needs improvement.
!
!===============================================================================
MODULE MOPAC_DENSITY

! . Module declarations.
USE DEFINITIONS,      ONLY : DP
USE FILES,            ONLY : NEXT_UNIT
USE PARSING
USE PRINTING,         ONLY : PRINT_ERROR, PRINT_PARAGRAPH
USE RANDOM_NUMBERS,   ONLY : RANDOM_VECTOR

USE ATOMS,            ONLY : NATOMSQM
USE MOPAC_DATA,       ONLY : DENMAT, DENMATA, DENMATB, NALPHA, NBETA, NBASIS, NBASTR, NPIBEADS, QMATOM, TOTCHG
USE MOPAC_PARAMETERS, ONLY : CORE, NATORB

IMPLICIT NONE
PUBLIC

!===============================================================================
CONTAINS
!===============================================================================

   !------------------------------------------------------------
   SUBROUTINE DENSITY_CALCULATE ( DENSITY, MORBS, NOCC, OCCNUM )
   !------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: NOCC
   REAL ( KIND = DP ), INTENT(IN) :: OCCNUM

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NBASTR),        INTENT(OUT) :: DENSITY
   REAL ( KIND = DP ), DIMENSION(1:NBASIS,1:NOCC), INTENT(IN)  :: MORBS

   ! . Local scalars.
   INTEGER :: I, IBASIS, JBASIS

   ! . Loop over each element in the density matrix.
   I = 0
   DO IBASIS = 1,NBASIS
      DO JBASIS = 1,IBASIS
         I = I + 1
         DENSITY(I) = SUM ( MORBS(IBASIS,1:NOCC) * MORBS(JBASIS,1:NOCC) )
      END DO
   END DO

   ! . Scale the density matrix by the occupation number.
   IF ( OCCNUM /= 1.0_DP ) DENSITY = OCCNUM * DENSITY

   END SUBROUTINE DENSITY_CALCULATE

   !-----------------------
   SUBROUTINE DENSITY_GUESS
   !-----------------------

   ! . Local scalars.
   INTEGER :: NELEBAS

   ! . Deallocate space for the density matrices.
   IF ( ALLOCATED ( DENMAT  ) ) DEALLOCATE ( DENMAT  )
   IF ( ALLOCATED ( DENMATA ) ) DEALLOCATE ( DENMATA )
   IF ( ALLOCATED ( DENMATB ) ) DEALLOCATE ( DENMATB )

   ! . Allocate the density matrices.
   ALLOCATE ( DENMAT(1:NBASTR,1:NPIBEADS) )
   IF ( NALPHA /= NBETA ) ALLOCATE ( DENMATA(1:NBASTR,1:NPIBEADS), DENMATB(1:NBASTR,1:NPIBEADS) )

   ! . Calculate the number of electrons on the atoms.
   NELEBAS = NALPHA + NBETA + TOTCHG

   ! . Total density only.
   IF ( NALPHA == NBETA ) THEN

      CALL GET_DENSITY ( DENMAT, NALPHA + NBETA, NELEBAS, 1.0_DP )

   ! . Separate densities.
   ELSE

      CALL GET_DENSITY ( DENMATA, NALPHA, NELEBAS, 0.5_DP )
      CALL GET_DENSITY ( DENMATB, NBETA,  NELEBAS, 0.5_DP )

      ! . Get the total density.
      DENMAT = DENMATA + DENMATB

   END IF

   !============================================================================
   CONTAINS
   !============================================================================

      !-----------------------------------------------------------
      SUBROUTINE GET_DENSITY ( DMATRIX, NELETOT, NELEBAS, OCCFAC )
      !-----------------------------------------------------------

      ! . Scalar arguments.
      INTEGER,            INTENT(IN) :: NELEBAS, NELETOT
      REAL ( KIND = DP ), INTENT(IN) :: OCCFAC

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: DMATRIX

      ! . Local parameters.
      REAL ( KIND = DP ), PARAMETER :: SMALL = 0.01_DP

      ! . Local scalars.
      INTEGER            :: I, IBASIS, ICOPY, II, IQM, NI, NORBS
      REAL ( KIND = DP ) :: FRACTQ, ORBCHG

      ! . Initialization.
      DMATRIX(:,1) = SMALL * RANDOM_VECTOR ( NBASTR )

      ! . Get the excess fractional charge per orbital.
      FRACTQ = ( REAL ( NELETOT, DP ) - OCCFAC * NELEBAS ) / REAL ( NBASIS, DP )

      ! . Loop over the QM atoms.
      DO IQM = 1,NATOMSQM

	 ! . Get some information for the atom.
	 ICOPY = QMATOM(IQM)%COPY
	 NI    = QMATOM(IQM)%NUMBER
	 NORBS = NATORB(NI)

	 ! . Skip atoms of the wrong index.
	 IF ( ICOPY > 1 ) CYCLE

	 ! . The atom has orbitals.
	 IF ( NORBS > 0 ) THEN
            ORBCHG = ( OCCFAC * CORE(NI) ) / REAL ( NORBS, DP ) - FRACTQ
            DO IBASIS = QMATOM(IQM)%BFIRST,QMATOM(IQM)%BLAST
               II = ( IBASIS * ( IBASIS + 1 ) ) / 2
               DMATRIX(II,1) = ORBCHG
            END DO
	 END IF
      END DO

      ! . Copy the first density matrix to the remaining arrays.
      IF ( NPIBEADS > 1 ) THEN
	 DO I = 2,NPIBEADS
            DMATRIX(:,I) = DMATRIX(:,1)
	 END DO
      END IF

      END SUBROUTINE GET_DENSITY

   END SUBROUTINE DENSITY_GUESS

   !-------------------------------
   SUBROUTINE DENSITY_READ ( FILE, PRINT )
   !-------------------------------
   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Local scalars.
   INTEGER :: I, IMAT, IOSTAT, NBASIN, NMATIN, NPIMIN, UNIT
   LOGICAL :: QPRINT
   
   ! . Check the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Get the next unit number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "DENSITY_READ", "I/O Error.", IOSTAT )

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Read in and check the number of basis functions and matrices in the file.
   CALL GET_LINE
   CALL GET_INTEGER ( NBASIN )
   CALL GET_INTEGER ( NPIMIN )
   CALL GET_INTEGER ( NMATIN )
   IF ( NBASIN /= NBASIS   ) CALL PARSE_ERROR ( "DENSITY_READ", "Invalid number of basis functions." )
   IF ( NPIMIN /= NPIBEADS ) CALL PARSE_ERROR ( "DENSITY_READ", "Invalid number of PI matrices." )
   IF ( ( NMATIN /= 1 ) .AND. &
        ( NMATIN /= 2 ) )    CALL PARSE_ERROR ( "DENSITY_READ", "Invalid number of input matrices." )

   ! . Check to see if the density matrix exists.
   IF ( .NOT. ALLOCATED ( DENMAT ) ) CALL PRINT_ERROR ( "DENSITY_READ", "Total density matrix not allocated." )
   IF ( ( NMATIN == 2 ) ) THEN
      IF ( .NOT. ALLOCATED ( DENMATA ) .OR. .NOT. ALLOCATED ( DENMATB ) ) &
                                     CALL PRINT_ERROR ( "DENSITY_READ", "Alpha or beta density matrices not allocated." )
   END IF

   ! . Read in the density matrix elements in free format (with no comment lines!).
   DO IMAT = 1,NPIBEADS
      IF ( NMATIN == 2 ) THEN
         READ ( UNIT,   *    ) ( DENMATA(I,IMAT), I = 1,NBASTR )
         READ ( UNIT,   *    ) ( DENMATB(I,IMAT), I = 1,NBASTR )
      ELSE
         READ ( UNIT,   *    ) ( DENMAT(I,IMAT), I = 1,NBASTR )
      END IF
      READ ( UNIT, "(1X)" )
   END DO

   ! . Fill the total density matrix in the case where there are two matrices.
   IF ( NMATIN == 2 ) DENMAT = DENMATA + DENMATB

   ! . Close the input file.
   CLOSE ( UNIT )

   ! . Reset the parsing unit.
   CALL POP_UNIT

   ! . Do some printing.
   IF (QPRINT) CALL PRINT_PARAGRAPH ( TEXT = "Density matrix read from " // TRIM ( FILE ) // "." )

   END SUBROUTINE DENSITY_READ

   !--------------------------------
   SUBROUTINE DENSITY_WRITE ( FILE, PRINT )
   !--------------------------------
   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Local scalars.
   INTEGER :: I, IMAT, IOSTAT, NMATRIX, UNIT
   LOGICAL :: QPRINT
   
   ! . Check the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF


   ! . A density matrix does not exist.
   IF ( ( NBASIS <= 0 ) .OR. ( .NOT. ALLOCATED ( DENMAT ) ) ) RETURN

   ! . Check for the presence of alpha and beta matrices.
   IF ( ALLOCATED ( DENMATA ) .AND. ALLOCATED ( DENMATB ) ) THEN
      NMATRIX = 2
   ELSE
      NMATRIX = 1
   END IF

   ! . Get the next unit number.
   UNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "UNKNOWN", IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "DENSITY_WRITE", "I/O Error.", IOSTAT )

   ! . Write out the separator.
   WRITE ( UNIT, "('!',119('='))" )

   ! . Write out the number of basis functions, PI beads and matrices.
   WRITE ( UNIT, "(3I6,A/)" ) NBASIS, NPIBEADS, NMATRIX, "  ! # of basis functions, PI beads and matrices."

   ! . Loop over the matrices.
   DO IMAT = 1,NPIBEADS

      ! . Two matrices.
      IF ( NMATRIX == 2 ) THEN
         WRITE ( UNIT, "(6G20.10)" ) ( DENMATA(I,IMAT), I = 1,NBASTR )
         WRITE ( UNIT, "(6G20.10)" ) ( DENMATB(I,IMAT), I = 1,NBASTR )
      ! . One matrix.
      ELSE
         WRITE ( UNIT, "(6G20.10)" ) ( DENMAT(I,IMAT), I = 1,NBASTR )
      END IF

      ! . Write out the terminator.
      WRITE ( UNIT, "('!',119('='))" )

   END DO

   ! . Close the input file.
   CLOSE ( UNIT )

   ! . Do some printing.
   IF (QPRINT) CALL PRINT_PARAGRAPH ( TEXT = "Density matrix written to " // TRIM ( FILE ) // "." )

   END SUBROUTINE DENSITY_WRITE

END MODULE MOPAC_DENSITY
