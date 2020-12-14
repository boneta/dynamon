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
!                        The MOPAC Quantum Energy Module
!===============================================================================
!
!
! . Hessian calculation scalars:
!
!   DELTA                      The hessian finite difference step size.
!
! . Subroutines:
!
!   ENERGY_QUANTUM             Calculate the QM energy.
!   ENERGY_PATH_INTEGRALS      Calculate the path integral harmonic term.
!   MOPAC_ENERGY               Calculate the MOPAC energy and gradients.
!   MOPAC_HESSIAN              Calculate the MOPAC Hessian.
!   MOPAC_HESSIAN_OPTIONS      Set the options for calculation of the Hessian.
!
!===============================================================================
MODULE QUANTUM_ENERGY

! . Module declarations.
USE CONSTANTS,       ONLY : EV_TO_KJ
USE DEFINITIONS,     ONLY : DP
USE PRINTING,        ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH

USE ATOMS,           ONLY : ATMCRD, ATMFIX, ATMIND, ATMMAS, NATOMS, NATOMSQM, NFIXED, NFREE
USE MOPAC_DATA,      ONLY : ATHEAT, NPIATOMS, NPIBEADS, PIFC, PILIST
USE MOPAC_GRADIENTS, ONLY : GRADIENTS
USE MOPAC_INTEGRALS, ONLY : INTEGRALS
USE MOPAC_SCF,       ONLY : MOPAC_SCF_CALCULATE

IMPLICIT NONE
PRIVATE
PUBLIC :: ENERGY_QUANTUM
#ifndef PGPC
SAVE
#endif

! . Module scalars.
REAL ( KIND = DP ) :: DELTA  = 1.0E-4_DP

!===============================================================================
CONTAINS
!===============================================================================

   !-----------------------------------------------------------------------
   SUBROUTINE ENERGY_QUANTUM ( EQM, EPI, VIRIAL, GRADIENT, HESSIAN, PRINT )
   !-----------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: EPI, EQM
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Local scalars.
   LOGICAL :: QGRADIENT, QHESSIAN, QPRINT

   ! . Check the number of QM atoms.
   IF ( NATOMSQM <= 0 ) RETURN

   ! . Set the gradient flag.
   QGRADIENT = PRESENT ( GRADIENT )
   QHESSIAN  = PRESENT ( HESSIAN  )

   ! . Check the consistency of the derivative options.
   IF ( QHESSIAN .AND. .NOT. QGRADIENT ) CALL PRINT_ERROR ( "ENERGY_QUANTUM", "First derivative argument missing." )

   ! . Set the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Calculate an energy and, if necessary, the gradients.
   CALL MOPAC_ENERGY ( EQM, QPRINT, VIRIAL, GRADIENT )

   ! . Calculate the Hessian if necessary.
   IF ( QHESSIAN ) CALL MOPAC_HESSIAN ( HESSIAN )

   ! . Calculate the PI energy if necessary.
   CALL ENERGY_PATH_INTEGRALS ( EPI, VIRIAL, GRADIENT, HESSIAN )

   END SUBROUTINE ENERGY_QUANTUM

   !------------------------------------------------------------------
   SUBROUTINE ENERGY_PATH_INTEGRALS ( EPI, VIRIAL, GRADIENT, HESSIAN )
   !------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: EPI
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN
  
   ! . Local scalars.
   INTEGER            :: I, IATOM, IBEAD, IFAC, INDEX, J, JFAC, SWAP
   LOGICAL            :: QGRADIENT, QHESSIAN
   REAL ( KIND = DP ) :: D2F, E, FC

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ

   ! . Initialization.
   EPI = 0.0_DP

   ! . Check the number of path integral particles.
   IF ( NPIATOMS <= 0 ) RETURN

   ! . Set the gradient flag.
   QGRADIENT = PRESENT ( GRADIENT )
   QHESSIAN  = PRESENT ( HESSIAN  )

   ! . Check the consistency of the derivative options.
   IF ( QHESSIAN .AND. .NOT. QGRADIENT ) CALL PRINT_ERROR ( "ENERGY_PATH_INTEGRALS", "First derivative argument missing." )

   ! . Loop over the path integral atoms.
   DO IATOM = 1,NPIATOMS

      ! . Get the mass prefactor (the bead masses must be the same).
      FC = ATMMAS(PILIST(IATOM,1)) * PIFC

      ! . Loop over the beads.
      DO IBEAD = 1,NPIBEADS

         ! . Get the first and second beads of the interaction.
         I = PILIST(IATOM,IBEAD)
         IF ( IBEAD == NPIBEADS ) THEN
            J = PILIST(IATOM,1)
         ELSE
            J = PILIST(IATOM,IBEAD+1)
         END IF

         ! . Calculate the displacement vector.
         DRIJ = ATMCRD(1:3,I) - ATMCRD(1:3,J)

         ! . Calculate the energy.
	 E = FC * DOT_PRODUCT ( DRIJ, DRIJ )

         ! . Add in the contribution to the total energy.
         EPI = EPI + E

         ! . Calculate the gradients if necessary.
         IF ( .NOT. QGRADIENT ) CYCLE

	 ! . Calculate the contribution to the virial.
	 VIRIAL = VIRIAL + ( 2.0_DP * E )

         ! . Calculate the gradient vector.
         DRIJ = 2.0_DP * FC * DRIJ

         ! . Add in the contribution to the total gradient.
         GRADIENT(1:3,I) = GRADIENT(1:3,I) + DRIJ
         GRADIENT(1:3,J) = GRADIENT(1:3,J) - DRIJ

         ! . Check for a Hessian calculation.
         IF ( .NOT. QHESSIAN ) CYCLE

         ! . Determine the SD element.
	 D2F = 2.0_DP * FC

	 ! . Calculate some index factors.
	 IFAC = 3 * ( ATMIND(I) - 1 )
	 JFAC = 3 * ( ATMIND(J) - 1 )

	 ! . Calculate the II block of the hessian.
	 IF ( ATMIND(I) > 0 ) THEN
            INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + D2F
            INDEX = INDEX + IFAC + 1
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + D2F
            INDEX = INDEX + IFAC + 2
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + D2F
	 END IF

	 ! . Calculate the JJ block of the hessian.
	 IF ( ATMIND(J) > 0 ) THEN
	    INDEX = ( JFAC * ( JFAC + 1 ) ) / 2 + JFAC + 1
	    HESSIAN(INDEX)   = HESSIAN(INDEX)   + D2F
	    INDEX = INDEX + JFAC + 1
	    HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + D2F
	    INDEX = INDEX + JFAC + 2
	    HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + D2F
	 END IF

	 ! . Swap the I and J indices if I is less than J.
	 IF ( I < J ) THEN
            SWAP = IFAC ; IFAC = JFAC ; JFAC = SWAP
	 END IF

	 ! . Calculate the IJ block of the hessian.
	 IF ( ( ATMIND(I) > 0 ) .AND. ( ATMIND(J) > 0 ) ) THEN
	    INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + JFAC + 1
	    HESSIAN(INDEX)   = HESSIAN(INDEX)   - D2F
	    INDEX = INDEX + IFAC + 1
	    HESSIAN(INDEX+1) = HESSIAN(INDEX+1) - D2F
	    INDEX = INDEX + IFAC + 2
	    HESSIAN(INDEX+2) = HESSIAN(INDEX+2) - D2F
      END IF

      END DO
   END DO

   END SUBROUTINE ENERGY_PATH_INTEGRALS

   !--------------------------------------------------------
   SUBROUTINE MOPAC_ENERGY ( EQM, QPRINT, VIRIAL, GRADIENT )
   !--------------------------------------------------------

   ! . Scalar arguments.
   LOGICAL,            INTENT(IN)    :: QPRINT
   REAL ( KIND = DP ), INTENT(OUT)   :: EQM
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT), OPTIONAL :: GRADIENT

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NPIBEADS) :: EHF, EELEC, ENUCL

   ! . Calculate the integrals.
   IF ( PRESENT ( GRADIENT ) ) THEN
      CALL INTEGRALS ( ENUCL, VIRIAL, GRADIENT )
   ELSE
      CALL INTEGRALS ( ENUCL, VIRIAL )
   END IF

   ! . Perform the SCF calculations.
   CALL MOPAC_SCF_CALCULATE ( EELEC, QPRINT )
!write ( 6, "(a,3f25.15)" ) "Energies = ", EELEC, ENUCL, ATHEAT
   ! . Calculate EHF.
   EHF = EV_TO_KJ * ( EELEC + ENUCL ) + ATHEAT

   ! . Calculate the average of the energies.
   EQM = SUM ( EHF ) / REAL ( NPIBEADS, DP )

   ! . Calculate the derivatives if necessary.
   IF ( PRESENT ( GRADIENT ) ) CALL GRADIENTS ( VIRIAL, GRADIENT )

   END SUBROUTINE MOPAC_ENERGY

   !-----------------------------------
   SUBROUTINE MOPAC_HESSIAN ( HESSIAN )
   !-----------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT) :: HESSIAN

   ! . Local scalars.
   INTEGER            :: I, IATOM, IFREE, INDEX, J
   REAL ( KIND = DP ) :: EDUMMY, SCALE, TEMP, VDUMMY

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: GMINUS, GTEMP, SD_MATRIX

   ! . Return if there are no free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Allocate a square second derivative matrix.
   ALLOCATE ( GMINUS(1:3,1:NFREE), GTEMP(1:3,1:NATOMS), SD_MATRIX(1:3*NFREE,1:3*NFREE) ) ; SD_MATRIX = 0.0_DP

   ! . Initialization.
   IFREE = 0

   ! . Loop over the atoms.
   DO IATOM = 1,NATOMS

      ! . Skip fixed atoms.
      IF ( ATMFIX(IATOM) ) CYCLE

      ! . Increment the free atom counter.
      IFREE = IFREE + 1

      ! . Loop over the Cartesian components.
      DO I = 1,3

         ! . Save the coordinate.
         TEMP = ATMCRD(I,IATOM)

         ! . Do the backward step.
         ! . Decrement the coordinate.
         ATMCRD(I,IATOM) = TEMP - DELTA

         ! . Calculate and save the gradient.
         GTEMP = 0.0_DP ; CALL MOPAC_ENERGY ( EDUMMY, .FALSE., VDUMMY, GTEMP ) ; GMINUS = CONTRACT ( )

         ! . Do the forward step.
         ! . Increment the coordinate.
         ATMCRD(I,IATOM) = TEMP + DELTA

         ! . Calculate and save the gradient.
         GTEMP = 0.0_DP ; CALL MOPAC_ENERGY ( EDUMMY, .FALSE., VDUMMY, GTEMP )

         ! . Restore the coordinate.
         ATMCRD(I,IATOM) = TEMP

         ! . Calculate the components of the second derivative matrix.
         INDEX = 3 * ( IFREE - 1 ) + I
         SD_MATRIX(1:3*NFREE,INDEX) = PACK ( ( CONTRACT ( ) - GMINUS ), .TRUE. )

      END DO
   END DO

   ! . Calculate the scale factor.
   SCALE = 0.5_DP / ( 2.0_DP * DELTA )

   ! . Add in the elements of SD_MATRIX to HESSIAN.
   INDEX = 0
   DO I = 1,3*NFREE
      DO J = 1,I
         INDEX = INDEX + 1
         HESSIAN(INDEX) = HESSIAN(INDEX) + SCALE * ( SD_MATRIX(I,J) + SD_MATRIX(J,I) )
      END DO
   END DO

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( GMINUS, GTEMP, SD_MATRIX )

   !============================================================================
   CONTAINS
   !============================================================================

      !--------------------
      FUNCTION CONTRACT ( )
      !--------------------

      ! . Function declarations.
      REAL ( KIND = DP ), DIMENSION(1:3,1:NFREE) :: CONTRACT

      ! . Local scalars.
      INTEGER :: IATOM, IFREE

      ! . There are no fixed atoms.
      IF ( NFIXED == 0 ) THEN
         CONTRACT = GTEMP
      ! . There are fixed atoms.
      ELSE

         ! . Loop over the atoms.
         IFREE = 0
         DO IATOM = 1,NATOMS
            IF ( ATMFIX(IATOM) ) CYCLE
            IFREE = IFREE + 1
            CONTRACT(1:3,IFREE) = GTEMP(1:3,IATOM)
         END DO

      END IF

      END FUNCTION CONTRACT

   END SUBROUTINE MOPAC_HESSIAN

   !----------------------------------------
   SUBROUTINE MOPAC_HESSIAN_OPTIONS ( STEP )
   !----------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: STEP

   ! . Set the options.
   IF ( PRESENT ( STEP ) ) DELTA  = STEP

   ! . Print out information about the options.
   WRITE ( PRINT_LINE, "(A,G16.8,A)" ) "MOPAC Hessian to be calculated with a step of ", DELTA, "."
   CALL PRINT_PARAGRAPH

   END SUBROUTINE MOPAC_HESSIAN_OPTIONS

END MODULE QUANTUM_ENERGY

