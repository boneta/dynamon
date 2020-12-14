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
!                         The Monte Carlo Energy Module
!===============================================================================
!
! . Scalar data:
!
!   BUFFER                          The smoothing buffer distance.
!   CUT                             The cutoff distance.
!   EPSILON                         The dielectric constant.
!
! . Subroutines:
!
!   MONTE_CARLO_ENERGY_OPTIONS      Define the energy options.
!
! . Functions:
!
!   MONTE_CARLO_ENERGY_FULL         Calculate the full system energy.
!   MONTE_CARLO_ENERGY_ONE          Calculate the interaction energy between
!                                   a molecule and the rest of the system.
!
! . Notes:
!
!   Neither fixed atoms or quantum atoms can be handled by this module.
!
!===============================================================================
MODULE MONTE_CARLO_ENERGY

! . Module declarations.
USE CONSTANTS,      ONLY : ELECT_CONST
USE DEFINITIONS,    ONLY : DP
USE PRINTING,       ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
                           PRINT_SUMMARY_START, PRINT_SUMMARY_STOP

USE ATOMS,          ONLY : ATMCRD, ATMMAS
USE MM_TERMS,       ONLY : ATMCHG, ATMEPS, ATMSIG
USE SEQUENCE,       ONLY : MOLIND => RESIND, NMOLECULES => NRESID
USE SYMMETRY,       ONLY : BOXL
USE TRANSFORMATION, ONLY : CENTER

IMPLICIT NONE
PRIVATE
PUBLIC :: MONTE_CARLO_ENERGY_FULL, MONTE_CARLO_ENERGY_ONE, MONTE_CARLO_ENERGY_OPTIONS
#ifndef PGPC
SAVE
#endif

! . Module scalars.
REAL ( KIND = DP ) :: BUFFER = 0.5_DP, CUT = 8.5_DP, EPSILON = 1.0_DP

!==============================================================================
CONTAINS
!==============================================================================

   !--------------------------------------------------------------------------
   SUBROUTINE MONTE_CARLO_ENERGY_OPTIONS ( CUTOFF, DIELECTRIC, SMOOTH, PRINT )
   !--------------------------------------------------------------------------

   ! . Optional scalar arguments.
   LOGICAL,            INTENT(IN), OPTIONAL :: PRINT
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: CUTOFF, DIELECTRIC, SMOOTH

   ! . Local scalars.
   LOGICAL :: QPRINT

   ! . The options.
   IF ( PRESENT ( CUTOFF     ) ) CUT     = CUTOFF
   IF ( PRESENT ( DIELECTRIC ) ) EPSILON = DIELECTRIC
   IF ( PRESENT ( SMOOTH     ) ) BUFFER  = SMOOTH

   ! . Check BUFFER.
   IF ( ( BUFFER < 0.0_DP ) .OR. ( BUFFER >= CUT ) ) BUFFER = 0.0_DP

   ! . Check for printing.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Write out the options.
   IF ( QPRINT ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF5500", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Monte Carlo Energy Options" )
      WRITE ( PRINT_LINE, "(F16.2)" ) CUT     ; CALL PRINT_SUMMARY_ELEMENT ( "Non-Bond Cutoff"     )
      WRITE ( PRINT_LINE, "(F16.2)" ) EPSILON ; CALL PRINT_SUMMARY_ELEMENT ( "Dielectric Constant" )
      WRITE ( PRINT_LINE, "(F16.2)" ) BUFFER  ; CALL PRINT_SUMMARY_ELEMENT ( "Dielectric Constant" )
      CALL PRINT_SUMMARY_STOP
   END IF

   END SUBROUTINE MONTE_CARLO_ENERGY_OPTIONS

   !------------------------------------------------------
   REAL ( KIND = DP ) FUNCTION MONTE_CARLO_ENERGY_FULL ( )
   !------------------------------------------------------

   ! . Local scalars.
   INTEGER            :: I, IMOL, ISTART, ISTOP, J, JMOL
   REAL ( KIND = DP ) :: CUTSQ, EE,  EELECT, EI, EL, ELJ, LOWERSQ, QI, RIJ2, SCALE, SI, SIJ2, SIJ6

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3)              :: DCORR, DR
   REAL ( KIND = DP ), DIMENSION(1:3,1:NMOLECULES) :: MOLCEN

   ! . Check BOXL.
   IF ( 2.0_DP * CUT >= MINVAL ( BOXL ) ) THEN
      CALL PRINT_ERROR ( "MONTE_CARLO_ENERGY_FULL", "The cutoff is too large for the box size." )
   END IF

   ! . Calculate some cutoff constants.
   CUTSQ   = CUT * CUT
   LOWERSQ = ( CUT - BUFFER )**2

   ! . Initialize the energies.
   EELECT = 0.0_DP
   ELJ    = 0.0_DP

   ! . Outer loop over the molecules.
   DO IMOL = 1,NMOLECULES

      ! . Get the starting and stopping atoms of the molecule.
      ISTART = MOLIND(IMOL)+1
      ISTOP  = MOLIND(IMOL+1)

      ! . Calculate the center of the molecule.
      MOLCEN(1:3,IMOL) = CENTER ( ATMCRD(1:3,ISTART:ISTOP), ATMMAS(ISTART:ISTOP) )

      ! . Inner loop over the molecules.
      DO JMOL = 1,(IMOL-1)

         ! . Calculate the distance between the molecule centers (applying the minimum image convention).
         DR    = MOLCEN(1:3,IMOL) - MOLCEN(1:3,JMOL)
         DCORR = - BOXL * ANINT ( DR / BOXL, DP )
         RIJ2  = SUM ( ( DR + DCORR )**2 )

         ! . Check to see whether the molecules are within the cutoff distance.
         IF ( RIJ2 < CUTSQ ) THEN

            ! . Initialize the local energy accumulators.
            EE = 0.0_DP
            EL = 0.0_DP

            ! . Loop over the atoms in the IMOLth molecule.
            DO I = ISTART,ISTOP

               ! . Get some data for the Ith atom.
               EI = ATMEPS(I)
               SI = ATMSIG(I)
               QI = ATMCHG(I)

               ! . Loop over the atoms in the JMOLth molecule.
               DO J = (MOLIND(JMOL)+1),MOLIND(JMOL+1)

                  ! . Calculate the inverse distance squared between atoms.
                  SIJ2 = 1.0_DP / SUM ( ( ATMCRD(1:3,I) - ATMCRD(1:3,J) + DCORR )**2 )

                  ! . Calculate the charge-charge interaction.
                  EE = EE + QI * ATMCHG(J) * SQRT ( SIJ2 )

                  ! . Calculate the Lennard-Jones interaction.
                  SIJ2 = ( SI * ATMSIG(J) )**2 * SIJ2
                  SIJ6 = SIJ2 * SIJ2 * SIJ2
                  EL   = EL + EI * ATMEPS(J) * SIJ6 * ( SIJ6 - 1.0_DP )

               END DO
            END DO

            ! . Calculate the scale factor.
            IF ( RIJ2 > LOWERSQ ) THEN
               SCALE = ( CUTSQ - RIJ2 ) / ( CUTSQ - LOWERSQ )
            ELSE
               SCALE = 1.0_DP
            END IF

            ! . Accumulate the total energies.
            EELECT = EELECT + EE * SCALE
            ELJ    = ELJ    + EL * SCALE

         END IF
      END DO
   END DO

   ! . Calculate the total energy.
   MONTE_CARLO_ENERGY_FULL = ( ELECT_CONST * EELECT / EPSILON ) + ELJ

   END FUNCTION MONTE_CARLO_ENERGY_FULL

   !------------------------------------------------------------
   REAL ( KIND = DP ) FUNCTION MONTE_CARLO_ENERGY_ONE ( CHOSEN )
   !------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: CHOSEN

   ! . Local scalars.
   INTEGER            :: I, ISTART, ISTOP, J, JMOL, JSTART, JSTOP
   REAL ( KIND = DP ) :: CUTSQ, EE, EELECT, EI, EL, ELJ, LOWERSQ, QI, RIJ2, SCALE, SI, SIJ2, SIJ6

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DCORR, DR, MOLCEN

   ! . Calculate some cutoff constants.
   CUTSQ   = CUT * CUT
   LOWERSQ = ( CUT - BUFFER )**2

   ! . Initialize the energies.
   EELECT = 0.0_DP
   ELJ    = 0.0_DP

   ! . Get the starting and stopping atoms in the CHOSEN molecule.
   ISTART = MOLIND(CHOSEN)+1
   ISTOP  = MOLIND(CHOSEN+1)

   ! . Calculate the center of the molecule.
   MOLCEN = CENTER ( ATMCRD(1:3,ISTART:ISTOP), ATMMAS(ISTART:ISTOP) )

   ! . Loop over the molecules.
   DO JMOL = 1,NMOLECULES

      ! . Skip over the self-interaction.
      IF ( JMOL == CHOSEN ) CYCLE

      ! . Get the atom indices for the molecule.
      JSTART = MOLIND(JMOL)+1
      JSTOP  = MOLIND(JMOL+1)

      ! . Calculate the distance between the molecule centers (applying the minimum image convention).
      DR    = MOLCEN - CENTER ( ATMCRD(1:3,JSTART:JSTOP), ATMMAS(JSTART:JSTOP) )
      DCORR = - BOXL * ANINT ( DR / BOXL, DP )
      RIJ2  = SUM ( ( DR + DCORR )**2 )

      ! . Check to see whether the molecules are within the cutoff distance.
      IF ( RIJ2 <= CUTSQ ) THEN

         ! . Initialize the local energy accumulators.
         EE = 0.0_DP
         EL = 0.0_DP

         ! . Loop over the atoms in the CHOSEN molecule.
         DO I = ISTART,ISTOP

            ! . Get some data for the Ith atom.
            EI = ATMEPS(I)
            SI = ATMSIG(I)
            QI = ATMCHG(I)

            ! . Loop over the atoms in the IMOLth molecule.
            DO J = JSTART,JSTOP

               ! . Calculate the inverse distance squared between atoms.
               SIJ2 = 1.0_DP / SUM ( ( ATMCRD(1:3,I) - ATMCRD(1:3,J) + DCORR )**2 )

               ! . Calculate the charge-charge interaction.
               EE = EE + QI * ATMCHG(J) * SQRT ( SIJ2 )

               ! . Calculate the Lennard-Jones interaction.
               SIJ2 = ( SI * ATMSIG(J) )**2 * SIJ2
               SIJ6 = SIJ2 * SIJ2 * SIJ2
               EL   = EL + EI * ATMEPS(J) * SIJ6 * ( SIJ6 - 1.0_DP )

            END DO
         END DO

         ! . Calculate the scale factor.
         IF ( RIJ2 > LOWERSQ ) THEN
            SCALE = ( CUTSQ - RIJ2 ) / ( CUTSQ - LOWERSQ )
         ELSE
            SCALE = 1.0_DP
         END IF

         ! . Accumulate the total energies.
         EELECT = EELECT + EE * SCALE
         ELJ    = ELJ    + EL * SCALE

      END IF
   END DO

   ! . Calculate the total energy.
   MONTE_CARLO_ENERGY_ONE = ( ELECT_CONST * EELECT / EPSILON ) + ELJ

   END FUNCTION MONTE_CARLO_ENERGY_ONE

END MODULE MONTE_CARLO_ENERGY
