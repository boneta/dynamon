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
!========================================================================
!            The Non-Bonding Energy Module (All Interactions)
!========================================================================
!
! . Subroutines:
!
!   ENERGY_NON_BONDING               Calculate the full energy.
!   ENERGY_NON_BONDING_INTERACTION   Calculate an interaction energy
!                                    between two groups of atoms.
!   ENERGY_NON_BONDING_OPTIONS       Process the non-bonding energy
!                                    options.
!   ENERGY_NON_BONDING_SELF          Calculate the self energy for a
!                                    group of atoms.
!
!   INTERACTIONS                     This is the subroutine that does
!                                    the work.
!
! . Notes:
!
!   This module calculates the non-bonding energy for a system without
!   the use of cutoffs or truncation schemes. The interactions are
!   calculated using a simple (not very efficient) double loop O(N^2)
!   procedure. Note that interactions between excluded atoms or between
!   atoms that are fixed are ignored.
!
!   This module cannot handle periodic boundary conditions or quantum
!   atoms.
!
!========================================================================
MODULE ENERGY_NON_BONDING !_FULL

! . Module declarations.
USE CONSTANTS,   ONLY : ELECT_CONST
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH

USE ATOMS,       ONLY : ATMCRD, ATMFIX, ATMIND, NATOMS, NATOMSQM
USE MM_TERMS
USE SYMMETRY,    ONLY : QBOX

IMPLICIT NONE
PRIVATE
PUBLIC :: ENERGY_NON_BONDING_CALCULATE, ENERGY_NON_BONDING_INTERACTION, ENERGY_NON_BONDING_OPTIONS, &
          ENERGY_NON_BONDING_SELF, &
          CUT_OFF, CUT_ON, NQMINT, QIMAGE, NNBLISTQM, NBLISTQM_FIRST, NBLISTQM_LAST
#ifndef PGPC
SAVE
#endif

! . Module scalars.
REAL ( KIND = DP ) :: EPSILON = 1.0_DP

! . The non-bonding interaction list type definition.
TYPE NBLIST_TYPE
   INTEGER :: ATOM
   INTEGER,           DIMENSION(:), POINTER :: INTERACTIONS
   TYPE(NBLIST_TYPE),               POINTER :: NEXT_LIST
END TYPE NBLIST_TYPE

! . Module scalars which are not used in the module but which are needed by the QM modules.
INTEGER            :: NQMINT = 0
LOGICAL            :: QIMAGE = .FALSE.
REAL ( KIND = DP ) :: CUT_ON = 99998.0_DP

! . The QM list variables.
INTEGER                    :: NNBLISTQM
TYPE(NBLIST_TYPE), POINTER :: NBLISTQM_FIRST, NBLISTQM_LAST

!========================================================================
CONTAINS
!========================================================================

   !---------------------------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_CALCULATE ( EELECT, ELJ, VIRIAL, GRADIENT, HESSIAN )
   !---------------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: EELECT,  ELJ
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),                INTENT(OUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NATOMS*(3*NATOMS+1))/2), INTENT(OUT), OPTIONAL :: HESSIAN

   ! . Local arrays.
   LOGICAL, DIMENSION(1:NATOMS) :: SELECTION1, SELECTION2

   ! . Initialize the atom selection arrays.
   SELECTION1 = .TRUE.
   SELECTION2 = .TRUE.

   ! . Calculate the interactions.
   CALL INTERACTIONS ( SELECTION1, SELECTION2, EELECT, ELJ, VIRIAL, GRADIENT, HESSIAN )

   END SUBROUTINE ENERGY_NON_BONDING_CALCULATE

   !--------------------------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_INTERACTION ( EELECT, ELJ, SELECTION1, SELECTION2 )
   !--------------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EELECT, ELJ

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION1, SELECTION2

   ! . Local scalars.
   REAL ( KIND = DP ) :: EELECT1, EELECT2, ELJ1, ELJ2, VIRIAL

   ! . Check for overlaps between the atom selections.
   IF ( ANY ( SELECTION1 .AND. SELECTION2 ) ) THEN
      CALL PRINT_ERROR ( "ENERGY_NON_BONDING_INTERACTION", "There is overlap between the atom selections." )
   END IF

   ! . Initialization.
   VIRIAL = 0.0_DP

   ! . Calculate the two sets of interactions.
   CALL INTERACTIONS ( SELECTION1, SELECTION2, EELECT1, ELJ1, VIRIAL )
   CALL INTERACTIONS ( SELECTION2, SELECTION1, EELECT2, ELJ2, VIRIAL )

   ! . Return the total interaction.
   EELECT = EELECT1 + EELECT2
   ELJ    = ELJ1    + ELJ2

   END SUBROUTINE ENERGY_NON_BONDING_INTERACTION

   !---------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_OPTIONS ( DIELECTRIC )
   !---------------------------------------------------

   ! . Optional scalar arguments.
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: DIELECTRIC

   ! . Process the input options.
   IF ( PRESENT ( DIELECTRIC ) ) EPSILON = DIELECTRIC

   ! . Write out the header.
   WRITE ( PRINT_LINE, "(A,F16.2,A)" ) "Full Non-Bonding Energy Calculation with Dielectric = ", EPSILON, "."
   CALL PRINT_PARAGRAPH

   END SUBROUTINE ENERGY_NON_BONDING_OPTIONS

   !------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_SELF ( EELECT, ELJ, SELECTION )
   !------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EELECT, ELJ

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION

   ! . Local scalars.
   REAL ( KIND = DP ) :: VIRIAL = 0.0_DP

   ! . Local arrays.
   LOGICAL, DIMENSION(1:NATOMS) :: SELECTION2

   ! . Copy the first selection to the second one.
   SELECTION2 = SELECTION

   ! . Calculate the interactions.
   CALL INTERACTIONS ( SELECTION, SELECTION2, EELECT, ELJ, VIRIAL )

   END SUBROUTINE ENERGY_NON_BONDING_SELF

!------------------------------------------------------------------------
! . Private subroutines.
!------------------------------------------------------------------------

   !-----------------------------------------------------------------------------------------
   SUBROUTINE INTERACTIONS ( SELECTION1, SELECTION2, EELECT, ELJ, VIRIAL, GRADIENT, HESSIAN )
   !-----------------------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: EELECT,  ELJ
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION1, SELECTION2

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),                INTENT(OUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NATOMS*(3*NATOMS+1))/2), INTENT(OUT), OPTIONAL :: HESSIAN

   ! . Local scalars.
   INTEGER            :: I, IFAC, INDEX, J, JFAC
   LOGICAL            :: QGRADIENT, QHESSIAN
   REAL ( KIND = DP ) :: DF, D2F, EPSFAC, EI, EIJ, EI14, ETEMP, HXX, HXY, HXZ, HYY, HYZ, HZZ, &
                         QI, QIJ, QI14, RIJ2, S, S6, SI, SIJ

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ

   ! . Initialization.
   EELECT = 0.0_DP
   ELJ    = 0.0_DP

   ! . Check for quantum atoms and periodic boundary conditions.
   IF ( ( NATOMSQM > 0 ) .OR. QBOX ) THEN
      CALL PRINT_ERROR ( "INTERACTIONS", "This non-bonding energy module does not support quantum atoms or SYMMETRY." )
   END IF

   ! . Check for the presence of derivatives.
   QGRADIENT = PRESENT ( GRADIENT )
   QHESSIAN  = PRESENT ( HESSIAN  )

   ! . Check the consistency of the derivative options.
   IF ( QHESSIAN .AND. .NOT. QGRADIENT ) CALL PRINT_ERROR ( "INTERACTIONS", "First derivative argument missing." )

   ! . Calculate the conversion factor for the electrostatic interactions.
   EPSFAC = ELECT_CONST / EPSILON

   ! . Loop over the first atom of the interaction.
   DO I = 1,(NATOMS-1)

      ! . Check for a selected atom in the first selection.
      IF ( .NOT. SELECTION1(I) ) CYCLE

      ! . Get some data for the Ith atom.
      EI   = ATMEPS(I)
      EI14 = ATMEPS14(I)
      SI   = ATMSIG(I)
      QI   = EPSFAC * ATMCHG(I)
      QI14 = EPSFAC * ATMCHG14(I)

      ! . Loop over the second atom of the interaction.
      DO J = (I+1),NATOMS

         ! . Check for a selected atom in the second selection.
         IF ( .NOT. SELECTION2(J) ) CYCLE

         ! . Skip an interaction between fixed atoms.
	 IF ( ATMFIX(I) .AND. ATMFIX(J) ) CYCLE

         ! . Check for an exclusion.
         IF ( ANY ( ATMEXCJ(ATMEXCI(I)+1:ATMEXCI(I+1)) == J ) ) THEN

            ! . Get some data for a 1-4 I/J interaction.
            IF ( ANY ( ATME14J(ATME14I(I)+1:ATME14I(I+1)) == J ) ) THEN
               EIJ = EI14 * ATMEPS14(J)
               QIJ = QI14 * ATMCHG14(J)
            ELSE
               CYCLE
            END IF

         ELSE

            ! . Get some data for the I/J interaction.
            EIJ = EI * ATMEPS(J)
            QIJ = QI * ATMCHG(J)

         END IF

         ! . Get the product of the LJ radii.
         SIJ = SI * ATMSIG(J)

         ! . Get the interatomic vector.
         DRIJ = ATMCRD(1:3,I) - ATMCRD(1:3,J)

         ! . Calculate the distance.
         RIJ2 = DOT_PRODUCT ( DRIJ, DRIJ )

         ! . Calculate some distance factors.
         S  = 1.0_DP / SQRT ( RIJ2 )
         S6 = ( SIJ * S )**6

         ! . Calculate the electrostatic interaction.
         ETEMP  = QIJ * S
         EELECT = EELECT + ETEMP

         ! . Calculate the Lennard-Jones interaction.
         ELJ = ELJ + EIJ * S6 * ( S6 - 1.0_DP )

         ! . Check for a gradient calculation.
         IF ( .NOT. QGRADIENT ) CYCLE

         ! . Calculate an intermediate for the gradient calculation.
         DF = ( - ETEMP + 6.0_DP * EIJ * S6 * ( 1.0_DP - 2.0_DP * S6 ) ) / RIJ2

         ! . Calculate the contribution to the virial.
         VIRIAL = VIRIAL + DF * RIJ2

         ! . Calculate the gradients.
         GRADIENT(1:3,I) = GRADIENT(1:3,I) + DF * DRIJ
         GRADIENT(1:3,J) = GRADIENT(1:3,J) - DF * DRIJ

         ! . Check for a Hessian calculation.
         IF ( .NOT. QHESSIAN ) CYCLE

         ! . Calculate an intermediate factor for the Hessian calculation.
         D2F = 3.0_DP * ETEMP + 24.0_DP * EIJ * S6 * ( 7.0_DP * S6 - 2.0_DP )

         ! . Scale DRIJ.
         DRIJ = DRIJ / RIJ2

         ! . Calculate the elements of the hessian block.
         HXX = DRIJ(1) * DRIJ(1) * D2F + DF
         HXY = DRIJ(1) * DRIJ(2) * D2F
         HXZ = DRIJ(1) * DRIJ(3) * D2F
         HYY = DRIJ(2) * DRIJ(2) * D2F + DF
         HYZ = DRIJ(2) * DRIJ(3) * D2F
         HZZ = DRIJ(3) * DRIJ(3) * D2F + DF

         ! . Calculate some index factors.
         IFAC = 3 * ( ATMIND(I) - 1 )
         JFAC = 3 * ( ATMIND(J) - 1 )

         ! . Calculate the II block of the hessian.
	 IF ( ATMIND(I) > 0 ) THEN
            INDEX = ( IFAC * ( IFAC + 1 ) ) / 2 + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXX
            INDEX = INDEX + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXY
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + HYY
            INDEX = INDEX + IFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXZ
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + HYZ
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + HZZ
	 END IF

         ! . Calculate the JJ block of the hessian.
	 IF ( ATMIND(J) > 0 ) THEN
            INDEX = ( JFAC * ( JFAC + 1 ) ) / 2 + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXX
            INDEX = INDEX + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXY
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + HYY
            INDEX = INDEX + JFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   + HXZ
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) + HYZ
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) + HZZ
	 END IF

         ! . Calculate the IJ block of the hessian.
         IF ( ( ATMIND(I) > 0 ) .AND. ( ATMIND(J) > 0 ) ) THEN
            INDEX = ( JFAC * ( JFAC + 1 ) ) / 2 + IFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   - HXX
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) - HXY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) - HXZ
            INDEX = INDEX + JFAC + 1
            HESSIAN(INDEX)   = HESSIAN(INDEX)   - HXY
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) - HYY
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) - HYZ
            INDEX = INDEX + JFAC + 2
            HESSIAN(INDEX)   = HESSIAN(INDEX)   - HXZ
            HESSIAN(INDEX+1) = HESSIAN(INDEX+1) - HYZ
            HESSIAN(INDEX+2) = HESSIAN(INDEX+2) - HZZ
	 END IF

      END DO
   END DO

   END SUBROUTINE INTERACTIONS

END MODULE ENERGY_NON_BONDING !_FULL
