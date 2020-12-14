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
!                      The MOPAC Quantum Properties Module
!===============================================================================
!
! . Subroutines.
!
!   QUANTUM_CHARGES                  Calculate the atom charges.
!   QUANTUM_DIPOLE                   Calculate the dipole moment.
!   QUANTUM_DIPOLE_DERIVATIVES       Calculate the dipole moment derivatives.
!   QUANTUM_SPIN_DENSITIES           Calculate the atom spin densities.
!
!===============================================================================
MODULE QUANTUM_PROPERTIES

! . Module declarations.
USE DEFINITIONS,      ONLY : DP
USE ELEMENTS,         ONLY : SYMBOL
USE PRINTING,         ONLY : PRINT_ERROR, PRINT_LINE, PRINT_TABLE_ELEMENT, PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP
USE TRANSFORMATION,   ONLY : CENTER_OF_MASS => CENTER

USE ATOMS,            ONLY : ATMCRD, ATMMAS, ATMNUM, NATOMS, NATOMSQM
USE MOPAC_DATA,       ONLY : DENMAT, DENMATA, DENMATB, MOPAC_DATA_ATOM_POSITION, NBASIS, NPIBEADS, QMATOM
USE MOPAC_PARAMETERS, ONLY : CORE, DD

IMPLICIT NONE
PRIVATE
PUBLIC :: QUANTUM_CHARGES, QUANTUM_DIPOLE, QUANTUM_DIPOLE_DERIVATIVES, QUANTUM_SPIN_DENSITIES

!===============================================================================
CONTAINS
!===============================================================================

   !-------------------------------------
   SUBROUTINE QUANTUM_CHARGES ( CHARGES )
   !-------------------------------------

   ! . Optional array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:NATOMS), INTENT(OUT) :: CHARGES

   ! . Local scalars.
   INTEGER :: IATOM, I, ICOPY, IFIRST, ILAST, IQM, NUMBER, P
 
   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NBASIS) :: BQTEMP

   ! . Check that there are QM atoms and a density matrix.
   IF ( ( NATOMSQM <= 0 ) .OR. ( NBASIS <= 0 ) .OR. ( .NOT. ALLOCATED ( DENMAT ) ) ) RETURN

   ! . Initialize the charge arrays.
   CHARGES = 0.0_DP

   ! . Loop over the density matrices.
   DO P = 1,NPIBEADS

      ! . Get the basis function charges.
      DO I = 1,NBASIS
         BQTEMP(I) = DENMAT((I*(I+1))/2,P)
      END DO

      ! . Loop over the QM atoms.
      DO IQM = 1,NATOMSQM

         ! . Get some information for the atom.
         IATOM  = QMATOM(IQM)%ATOM
         ICOPY  = QMATOM(IQM)%COPY
         IFIRST = QMATOM(IQM)%BFIRST
         ILAST  = QMATOM(IQM)%BLAST
         NUMBER = QMATOM(IQM)%NUMBER

         ! . Skip inappropriate atoms.
         IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

         ! . Set the charge.
         CHARGES(IATOM) = CHARGES(IATOM) + ( CORE(NUMBER) - SUM ( BQTEMP(IFIRST:ILAST) ) )

      END DO
   END DO

   ! . Scale the charge arrays.
   IF ( NPIBEADS > 1 ) THEN
      DO IQM = 1,NATOMSQM
         IATOM = QMATOM(IQM)%ATOM
         IF ( QMATOM(IQM)%COPY == 0 ) THEN
            CHARGES(IATOM) = CHARGES(IATOM) / REAL ( NPIBEADS, DP )
         END IF
      END DO
   END IF

   END SUBROUTINE QUANTUM_CHARGES

   !-------------------------------------------
   SUBROUTINE QUANTUM_DIPOLE ( DIPOLE, CENTER )
   !-------------------------------------------

   ! . Optional array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(OUT)          :: DIPOLE
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN), OPTIONAL :: CENTER

   ! . Local parameter declarations.
   ! . Note: DIPFC1 = AU_to_Debyes * ANGSTROMS_TO_BOHRS and DIPFC2 = 2 * AU_to_Debyes.
   REAL ( KIND = DP ), PARAMETER :: DIPFC1 = 4.803_DP, DIPFC2 = 5.0832_DP

   ! . Local scalars.
   INTEGER            :: I, IATOM, ICOPY, IELMNT, IFIRST, II, ILAST, IQM, NUMBER, P
   REAL ( KIND = DP ) :: QTOT, SUM

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3)   :: CENMAS, DIPQME, DIPQMN
   REAL ( KIND = DP ), DIMENSION(1:100) :: HYF

   ! . Check that there are QM atoms and a density matrix.
   IF ( ( NATOMSQM <= 0 ) .OR. ( NBASIS <= 0 ) .OR. ( .NOT. ALLOCATED ( DENMAT ) ) ) RETURN

   ! . Initialization.
   DIPOLE = 0.0_DP
   DIPQME = 0.0_DP
   DIPQMN = 0.0_DP

   ! . A center is present.
   IF ( PRESENT ( CENTER ) ) THEN
      CENMAS = CENTER
   ! . Calculate the centre of mass of the system.
   ELSE
      CENMAS = CENTER_OF_MASS ( ATMCRD, ATMMAS )
   END IF

   ! . Calculate the hybridisation terms.
   HYF(1)     = 0.0_DP
   HYF(2:100) = DIPFC2 * DD(2:100)

   ! . Loop over the density matrices.
   DO P = 1,NPIBEADS

      ! . Loop over the QM atoms.
      DO IQM = 1,NATOMSQM

         ! . Get some information for the atom.
         IATOM  = QMATOM(IQM)%ATOM
         ICOPY  = QMATOM(IQM)%COPY
         IFIRST = QMATOM(IQM)%BFIRST
         ILAST  = QMATOM(IQM)%BLAST
         NUMBER = QMATOM(IQM)%NUMBER

         ! . Skip an inappropriate atom.
         IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

         ! . Calculate the total charge on the atom.
         SUM    = 0.0_DP
         IELMNT = (IFIRST*(IFIRST-1))/2
         DO I = IFIRST,ILAST
            IELMNT = IELMNT + I
            SUM    = SUM + DENMAT(IELMNT,P)
         END DO
         QTOT = CORE(NUMBER) - SUM

         ! . Calculate the electronic dipole moment.
         ILAST  = ILAST - IFIRST
         DO I = 1,ILAST
            II = ((IFIRST+I) * (IFIRST+I-1))/2 + IFIRST
            DIPQME(I) = DIPQME(I) - ( HYF(NUMBER) * DENMAT(II,P) )
         END DO

         ! . Calculate the nuclear dipole moment.
         DIPQMN = DIPQMN + ( DIPFC1 * QTOT * ( MOPAC_DATA_ATOM_POSITION ( IQM, ATMCRD ) - CENMAS ) )

      END DO
   END DO

   ! . Scale the dipole moment vectors.
   DIPQME = DIPQME / REAL ( NPIBEADS, DP )
   DIPQMN = DIPQMN / REAL ( NPIBEADS, DP )

   ! . Sum up the total dipole moment.
   DIPOLE = DIPQME + DIPQMN

   END SUBROUTINE QUANTUM_DIPOLE

   !----------------------------------------------
   SUBROUTINE QUANTUM_DIPOLE_DERIVATIVES ( DMUDR )
   !----------------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3*NATOMS,1:3), INTENT(OUT) :: DMUDR

   ! . Check that there are QM atoms and a density matrix.
   IF ( ( NATOMSQM <= 0 ) .OR. ( NBASIS <= 0 ) .OR. ( .NOT. ALLOCATED ( DENMAT ) ) ) RETURN

   ! . Initialization.
   DMUDR = 0.0_DP

   ! . Call an error.
   CALL PRINT_ERROR ( "QUANTUM_DIPOLE_DERIVATIVES", "MOPAC dipole derivatives are unavailable." )

   END SUBROUTINE QUANTUM_DIPOLE_DERIVATIVES

   !---------------------------------------------------
   SUBROUTINE QUANTUM_SPIN_DENSITIES ( SPIN_DENSITIES )
   !---------------------------------------------------

   ! . Optional array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:NATOMSQM), INTENT(OUT), OPTIONAL :: SPIN_DENSITIES

   ! . Local scalars.
   INTEGER :: I, ICOPY, IFIRST, ILAST, IQM, P
 
   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NATOMSQM) :: SPIND
   REAL ( KIND = DP ), DIMENSION(1:NBASIS)   :: BQTEMP

   ! . Initialize the spin density arrays.
   SPIND = 0.0_DP

   ! . Check that there are QM atoms and a density matrix.
   IF ( ( NATOMSQM > 0 ) .AND. ( NBASIS > 0 ) .AND. ALLOCATED ( DENMATA ) .AND. ALLOCATED ( DENMATB ) ) THEN

      ! . Loop over the density matrices.
      DO P = 1,NPIBEADS

	 ! . Get the basis function spin densities.
	 DO I = 1,NBASIS
            BQTEMP(I) = DENMATA((I*(I+1))/2,P) - DENMATB((I*(I+1))/2,P)
	 END DO

	 ! . Loop over the QM atoms.
	 DO IQM = 1,NATOMSQM

            ! . Get some information for the atom.
            ICOPY  = QMATOM(IQM)%COPY
            IFIRST = QMATOM(IQM)%BFIRST
            ILAST  = QMATOM(IQM)%BLAST

            ! . Skip inappropriate atoms.
            IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

            ! . Set the spin density.
            SPIND(IQM) = SPIND(IQM) + SUM ( BQTEMP(IFIRST:ILAST) )

	 END DO
      END DO

      ! . Scale the spin density arrays.
      IF ( NPIBEADS > 1 ) THEN
	 DO IQM = 1,NATOMSQM
            IF ( QMATOM(IQM)%COPY == 0 ) SPIND(IQM) = SPIND(IQM) / REAL ( NPIBEADS, DP )
	 END DO
      END IF

      ! . Write out the header.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 4, HEADER_COLOR = "#FF8800", VARIABLEWIDTHS = (/ 20, 20, 20, 20 /) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Quantum Atom Spin Densities", COLSPAN = 4, HEADER = .TRUE. )

      ! . Print out the spin densities.
      DO I = 1,NATOMSQM
         WRITE ( PRINT_LINE, "(I4,2X,A,I3,2X,F7.3)" ) I, SYMBOL(ATMNUM(QMATOM(I)%ATOM)), ATMNUM(QMATOM(I)%ATOM), SPIND(I)
         CALL PRINT_TABLE_ELEMENT
      END DO

      ! . Write out the terminator.
      CALL PRINT_TABLE_STOP

   END IF

   ! . Copy SPIND to SPIN_DENSITIES if necessary.
   IF ( PRESENT ( SPIN_DENSITIES ) ) SPIN_DENSITIES = SPIND

   END SUBROUTINE QUANTUM_SPIN_DENSITIES

END MODULE QUANTUM_PROPERTIES
