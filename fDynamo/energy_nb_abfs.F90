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
!       The Non-Bonding Energy Module (OPLS - Atom-Based Force-Switching)
!===============================================================================
!
! . Module scalars:
!
!   CUT_LIST                        The cutoff for inclusion on the list.
!   CUT_OFF                         The outer cutoff for the interactions.
!   CUT_ON                          The inner cutoff for the interactions.
!   EPSILON                         The dielectric constant.
!   NCALLS                          The number of calls.
!   NQMINT                          The number of QM/MM interactions.
!   NUPDATES                        The number of updates.
!   QIMAGE                          The minimum image flag.
!   QNBPRINT                        The update printing flag.
!   SUM_INTM, SUM_INTQ              The accumulators for the total number
!                                   of interactions.
!
! . Module arrays:
!
!   REFCRD                          The reference atom coordinates.
!   REFFIX                          The reference fixed atom array.
!
! . Module lists:
!
!   NBLIST14_FIRST, NBLIST14_LAST   The 1,4 MM/MM and QM/MM interaction lists.
!   NBLISTMM_FIRST, NBLISTMM_LAST   The MM/MM interaction lists (with QM/MM LJ).
!   NBLISTQM_FIRST, NBLISTQM_LAST   The QM/MM electrostatic interaction lists.
!
! . Subroutines (public):
!
!   ENERGY_NON_BONDING_CALCULATE    Calculate the non-bonding energy.
!   ENERGY_NON_BONDING_OPTIONS      Set the non-bonding options.
!   ENERGY_NON_BONDING_STATISTICS   Print out some non-bonding statistics.
!
! . Interaction subroutines (public):
!
!   ENERGY_NON_BONDING_INTERACTION  Calculate the interaction energy.
!   ENERGY_NON_BONDING_SELF         Calculate the self energy.
!
! . Subroutines (private):
!
!   GENERATE_REDUCED_LISTS          Reduce the non-bonding lists given two atom
!                                   selections.
!   INTERACTIONS                    Calculate the interactions given a list.
!   NBLIST_APPEND                   Append interactions to a list.
!   NBLIST_INITIALIZE               Initialize a non-bonding interaction list.
!   UPDATE                          Update the non-bonding interaction lists.
!
! . Notes:
!
!   Calculate the non-bonding energy using an atom-based force-switching
!   truncation scheme applied to both the electrostatic and LJ energies. The
!   non-bonding interaction lists are generated and updated automatically using
!   a residue-based search scheme.
!
!   Interactions between fixed atoms and between QM atoms are removed from the
!   interaction lists. It is assumed that the charges on the QM atoms are zero!
!
!   This module calculates the MM/MM interactions and updates both the QM/MM and
!   MM/MM lists.
!
!   The QM boundary atoms have both QM and MM electrostatic interactions and so
!   enter into both lists.
!
!===============================================================================
MODULE ENERGY_NON_BONDING_ABFS

! . Module declarations.
USE CONSTANTS,   ONLY : ELECT_CONST
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_ERROR, PRINT_LINE, PRINT_LINEBREAK, PRINT_PARAGRAPH_START,    &
                        PRINT_PARAGRAPH_STOP, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
                        PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TEXT

USE ATOMS,       ONLY : ATMCRD, ATMFIX, ATMIND, ATMQMI, NATOMS, NATOMSMM, NATOMSQM, NFREE
USE MM_TERMS
USE MOPAC_DATA,  ONLY : MOPAC_DATA_ATOM_POSITION, QMATOM, QUSEBAQM
USE SEQUENCE,    ONLY : NRESID, RESIND
USE SYMMETRY,    ONLY : BOXL, QBOX

IMPLICIT NONE
PRIVATE
PUBLIC :: CUT_OFF, CUT_ON, ENERGY_NON_BONDING_CALCULATE, ENERGY_NON_BONDING_INTERACTION,      &
          ENERGY_NON_BONDING_OPTIONS, ENERGY_NON_BONDING_SELF, ENERGY_NON_BONDING_STATISTICS, &
          NBLIST_TYPE, NBLISTQM_FIRST, NNBLISTQM, NQMINT, QIMAGE,                             &
          NBLIST14_FIRST, NBLIST14_LAST, NBLISTMM_FIRST, NBLISTMM_LAST, NNBLIST14, NNBLISTMM, EPSILON
#ifndef PGPC
SAVE
#endif

! . The non-bonding interaction list type definition.
TYPE NBLIST_TYPE
   INTEGER :: ATOM
   INTEGER,           DIMENSION(:), POINTER :: INTERACTIONS
   TYPE(NBLIST_TYPE),               POINTER :: NEXT_LIST
END TYPE NBLIST_TYPE

! . General module scalars.
INTEGER            :: NCALLS = 0, NQMINT = 0, NUPDATES = 0
LOGICAL            :: QIMAGE   = .FALSE., QNBPRINT = .FALSE., QNONBONDING = .TRUE.
REAL ( KIND = DP ) :: CUT_LIST = 999999.0_DP, CUT_OFF  = 990000.0_DP, CUT_ON   = 980000.0_DP, &
                      EPSILON  =      1.0_DP, SUM_INTM =      0.0_DP, SUM_INTQ =      0.0_DP

! . Module arrays.
LOGICAL,            ALLOCATABLE, DIMENSION(:)   :: REFFIX
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: REFCRD

! . The number of records in the 1,4 MM/MM, MM/MM (plus QM/MM LJ) and QM/MM electrostatics lists.
INTEGER :: NNBLIST14 = -1, NNBLISTMM = -1, NNBLISTQM = -1

! . The list type definitions.
TYPE(NBLIST_TYPE), POINTER :: NBLIST14_FIRST, NBLIST14_LAST, &
                              NBLISTMM_FIRST, NBLISTMM_LAST, &
                              NBLISTQM_FIRST, NBLISTQM_LAST

#define DEBUG

!===============================================================================
CONTAINS
!===============================================================================

   !---------------------------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_CALCULATE ( EELECT, ELJ, VIRIAL, GRADIENT, HESSIAN )
   !---------------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: EELECT, ELJ
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Local scalars.
   REAL ( KIND = DP ) :: EELECT14, ELJ14

   !---------------------------------------------------------------------
   ! . Initialization.
   !---------------------------------------------------------------------
   ! . Initialize the energy values.
   EELECT   = 0.0_DP
   EELECT14 = 0.0_DP
   ELJ      = 0.0_DP
   ELJ14    = 0.0_DP

   ! . Return if there are no MM atoms or no residues.
   IF ( ( NATOMSMM <= 0 ) .OR. ( NRESID <= 0 ) ) RETURN

   ! . Check that non-bonding interactions are to be calculated.
   IF ( .NOT. QNONBONDING ) RETURN

   ! . Update the non-bonding interaction lists if necessary.
   CALL UPDATE

   ! . Calculate the 1-4 interactions.
   IF ( NNBLIST14 > 0 ) THEN

      ! . Calculate the interactions.
      CALL INTERACTIONS ( EELECT14, ELJ14, VIRIAL, NNBLIST14, NBLIST14_FIRST, ATMCHG14, ATMEPS14, ATMSIG, &
                                                                                        GRADIENT, HESSIAN )

   END IF

   ! . Calculate the normal interactions.
   IF ( NNBLISTMM > 0 ) THEN

      ! . Calculate the interactions.
      CALL INTERACTIONS ( EELECT, ELJ, VIRIAL, NNBLISTMM, NBLISTMM_FIRST, ATMCHG, ATMEPS, ATMSIG, GRADIENT, HESSIAN )

   END IF

   ! . Sum the total energies.
   EELECT = EELECT + EELECT14
   ELJ    = ELJ    + ELJ14

   ! . Increment the NCALLS counter.
   NCALLS = NCALLS + 1

   END SUBROUTINE ENERGY_NON_BONDING_CALCULATE

   !------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_OPTIONS ( LIST_CUTOFF, OUTER_CUTOFF, INNER_CUTOFF, DIELECTRIC, MINIMUM_IMAGE, PRINT )
   !------------------------------------------------------------------------------------------------------------------

   ! . Optional scalar arguments.
   LOGICAL,            INTENT(IN), OPTIONAL :: MINIMUM_IMAGE, PRINT
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: LIST_CUTOFF, OUTER_CUTOFF, INNER_CUTOFF, DIELECTRIC

   ! . Reinitialize the statistics counters.
   NCALLS   = 0
   NQMINT   = 0
   NUPDATES = 0
   SUM_INTM = 0.0_DP
   SUM_INTQ = 0.0_DP

   ! . Deallocate the module arrays.
   IF ( ALLOCATED ( REFCRD ) ) DEALLOCATE ( REFCRD )
   IF ( ALLOCATED ( REFFIX ) ) DEALLOCATE ( REFFIX )

   ! . Deallocate the module non-bonding lists.
   CALL NBLIST_INITIALIZE ( NNBLIST14, NBLIST14_FIRST, NBLIST14_LAST )
   CALL NBLIST_INITIALIZE ( NNBLISTMM, NBLISTMM_FIRST, NBLISTMM_LAST )
   CALL NBLIST_INITIALIZE ( NNBLISTQM, NBLISTQM_FIRST, NBLISTQM_LAST )

   ! . The covalent energy terms.
   IF ( PRESENT ( LIST_CUTOFF   ) ) CUT_LIST   = LIST_CUTOFF
   IF ( PRESENT ( OUTER_CUTOFF  ) ) CUT_OFF    = OUTER_CUTOFF
   IF ( PRESENT ( INNER_CUTOFF  ) ) CUT_ON     = INNER_CUTOFF
   IF ( PRESENT ( DIELECTRIC    ) ) EPSILON    = DIELECTRIC
   IF ( PRESENT ( MINIMUM_IMAGE ) ) QIMAGE     = MINIMUM_IMAGE
   IF ( PRESENT ( PRINT         ) ) QNBPRINT   = PRINT

   ! . Check the values of the cutoffs.
   IF ( CUT_LIST < CUT_OFF ) CALL PRINT_ERROR ( "ENERGY_NON_BONDING_OPTIONS", "Invalid CUT_LIST/CUT_OFF values." )
   IF ( CUT_OFF  < CUT_ON  ) CUT_ON = CUT_OFF

   ! . Set the non-bonding interaction flag.
   QNONBONDING = ( CUT_LIST > 0.0_DP )

   ! . Write out some information.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#550055", VARIABLEWIDTH = 16 )
   CALL PRINT_SUMMARY_START ( "Atom-Based Force-Switching Non-Bonding Energy Options" )
   WRITE ( PRINT_LINE, "(F16.2)" ) CUT_LIST ; CALL PRINT_SUMMARY_ELEMENT ( "List Cutoff"          )
   WRITE ( PRINT_LINE, "(F16.2)" ) EPSILON  ; CALL PRINT_SUMMARY_ELEMENT ( "Dielectric Constant"  )
   WRITE ( PRINT_LINE, "(F16.2)" ) CUT_OFF  ; CALL PRINT_SUMMARY_ELEMENT ( "Outer Switch Cutoff"  )
   WRITE ( PRINT_LINE, "(F16.2)" ) CUT_ON   ; CALL PRINT_SUMMARY_ELEMENT ( "Inner Switch Cutoff"  )
   WRITE ( PRINT_LINE, "(L16)"   ) QIMAGE   ; CALL PRINT_SUMMARY_ELEMENT ( "Minimum Image Option" )
   WRITE ( PRINT_LINE, "(L16)"   ) QNBPRINT ; CALL PRINT_SUMMARY_ELEMENT ( "Update Printing"      )
   CALL PRINT_SUMMARY_STOP

   END SUBROUTINE ENERGY_NON_BONDING_OPTIONS

   !---------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_STATISTICS
   !---------------------------------------

   ! . Write out some information.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#550055", VARIABLEWIDTH = 16 )
   CALL PRINT_SUMMARY_START ( "Atom-Based Force-Switching Non-Bonding Energy Statistics" )
   WRITE ( PRINT_LINE, "(I16)" ) NCALLS   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Calls"   )
   WRITE ( PRINT_LINE, "(I16)" ) NUPDATES ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Updates" )
   IF ( NUPDATES > 0 ) THEN
      WRITE ( PRINT_LINE, "(F16.4)" ) REAL ( NCALLS, DP ) / REAL ( NUPDATES, DP )
      CALL PRINT_SUMMARY_ELEMENT ( "Calls per Update"   )
      WRITE ( PRINT_LINE, "(F16.4)" ) SUM_INTM / REAL ( NUPDATES, DP )
      CALL PRINT_SUMMARY_ELEMENT ( "<# of MM/MM Ints.>" )
      WRITE ( PRINT_LINE, "(F16.4)" ) SUM_INTQ / REAL ( NUPDATES, DP )
      CALL PRINT_SUMMARY_ELEMENT ( "<# of QM/MM Ints.>" )
   END IF
   CALL PRINT_SUMMARY_STOP

   ! . Reinitialize the integer counters.
   NCALLS   = 0
! . This is wrong as it is required in MOPAC_INTEGRALS.
! . The alternative is to remove all NB lists.
!   NQMINT   = 0
   NUPDATES = 0
   SUM_INTM = 0.0_DP
   SUM_INTQ = 0.0_DP

   END SUBROUTINE ENERGY_NON_BONDING_STATISTICS

!-------------------------------------------------------------------------------
! . Interaction and Self-Energy Subroutines.
!-------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_INTERACTION ( EELECT, ELJ, SELECTION1, SELECTION2 )
   !--------------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EELECT, ELJ

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION1, SELECTION2

   ! . Local scalars.
   INTEGER            :: NNBTEMP
   REAL ( KIND = DP ) :: EELECT14, ELJ14, VIRIAL

   ! . Local types.
   TYPE(NBLIST_TYPE), POINTER :: NBTEMP_FIRST, NBTEMP_LAST

   ! . Initialize the energies.
   EELECT   = 0.0_DP
   EELECT14 = 0.0_DP
   ELJ      = 0.0_DP
   ELJ14    = 0.0_DP

   ! . Return if there are no MM atoms or no residues.
   IF ( ( NATOMSMM <= 0 ) .OR. ( NRESID <= 0 ) ) RETURN

   ! . Check that non-bonding interactions are to be calculated.
   IF ( .NOT. QNONBONDING ) RETURN

   ! . Check for overlaps between the atom selections.
   IF ( ANY ( SELECTION1 .AND. SELECTION2 ) ) THEN
      CALL PRINT_ERROR ( "ENERGY_NON_BONDING_INTERACTION", "There is overlap between the atom selections." )
   END IF

   ! . Update the non-bond lists if necessary.
   CALL UPDATE

   ! . Calculate the 1-4 interactions.
   IF ( NNBLIST14 > 0 ) THEN

      ! . Initialize the temporary list record counter.
      NNBTEMP = 0 ; NULLIFY ( NBTEMP_FIRST, NBTEMP_LAST )

      ! . Get the reduced 1-4 non-bond lists.
      CALL GENERATE_REDUCED_LISTS ( SELECTION1, SELECTION2, NNBLIST14, NBLIST14_FIRST, NNBTEMP, NBTEMP_FIRST, NBTEMP_LAST )

      ! . Calculate the interactions.
      CALL INTERACTIONS ( EELECT14, ELJ14, VIRIAL, NNBTEMP, NBTEMP_FIRST, ATMCHG14, ATMEPS14, ATMSIG )

      ! . Initialize the temporary lists.
      CALL NBLIST_INITIALIZE ( NNBTEMP, NBTEMP_FIRST, NBTEMP_LAST )

   END IF

   ! . Calculate the normal interactions.
   IF ( NNBLISTMM > 0 ) THEN

      ! . Initialize the temporary list record counter.
      NNBTEMP = 0 ; NULLIFY ( NBTEMP_FIRST, NBTEMP_LAST )

      ! . Get the reduced 1-4 non-bond lists.
      CALL GENERATE_REDUCED_LISTS ( SELECTION1, SELECTION2, NNBLISTMM, NBLISTMM_FIRST, NNBTEMP, NBTEMP_FIRST, NBTEMP_LAST )

      ! . Calculate the interactions.
      CALL INTERACTIONS ( EELECT, ELJ, VIRIAL, NNBTEMP, NBTEMP_FIRST, ATMCHG, ATMEPS, ATMSIG )

      ! . Initialize the temporary lists.
      CALL NBLIST_INITIALIZE ( NNBTEMP, NBTEMP_FIRST, NBTEMP_LAST )

   END IF

   ! . Sum the total energies.
   EELECT = EELECT + EELECT14
   ELJ    = ELJ    + ELJ14

   END SUBROUTINE ENERGY_NON_BONDING_INTERACTION

   !------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_SELF ( EELECT, ELJ, SELECTION )
   !------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EELECT, ELJ

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION

   ! . Local scalars.
   INTEGER            :: NNBTEMP
   REAL ( KIND = DP ) :: EELECT14, ELJ14, VIRIAL

   ! . Local types.
   TYPE(NBLIST_TYPE), POINTER :: NBTEMP_FIRST, NBTEMP_LAST

   ! . Initialize the energies.
   EELECT   = 0.0_DP
   EELECT14 = 0.0_DP
   ELJ      = 0.0_DP
   ELJ14    = 0.0_DP

   ! . Return if there are no MM atoms or no residues.
   IF ( ( NATOMSMM <= 0 ) .OR. ( NRESID <= 0 ) ) RETURN

   ! . Check that non-bonding interactions are to be calculated.
   IF ( .NOT. QNONBONDING ) RETURN

   ! . Update the non-bond lists if necessary.
   CALL UPDATE

   ! . Calculate the 1-4 interactions.
   IF ( NNBLIST14 > 0 ) THEN

      ! . Initialize the temporary list record counter.
      NNBTEMP = 0 ; NULLIFY ( NBTEMP_FIRST, NBTEMP_LAST )

      ! . Get the reduced 1-4 non-bond lists.
      CALL GENERATE_REDUCED_LISTS ( SELECTION, SELECTION, NNBLIST14, NBLIST14_FIRST, NNBTEMP, NBTEMP_FIRST, NBTEMP_LAST )

      ! . Calculate the interactions.
      CALL INTERACTIONS ( EELECT14, ELJ14, VIRIAL, NNBTEMP, NBTEMP_FIRST, ATMCHG14, ATMEPS14, ATMSIG )

      ! . Initialize the temporary lists.
      CALL NBLIST_INITIALIZE ( NNBTEMP, NBTEMP_FIRST, NBTEMP_LAST )

   END IF

   ! . Calculate the normal interactions.
   IF ( NNBLISTMM > 0 ) THEN

      ! . Initialize the temporary list record counter.
      NNBTEMP = 0 ; NULLIFY ( NBTEMP_FIRST, NBTEMP_LAST )

      ! . Get the reduced 1-4 non-bond lists.
      CALL GENERATE_REDUCED_LISTS ( SELECTION, SELECTION, NNBLISTMM, NBLISTMM_FIRST, NNBTEMP, NBTEMP_FIRST, NBTEMP_LAST )

      ! . Calculate the interactions.
      CALL INTERACTIONS ( EELECT, ELJ, VIRIAL, NNBTEMP, NBTEMP_FIRST, ATMCHG, ATMEPS, ATMSIG )

      ! . Initialize the temporary lists.
      CALL NBLIST_INITIALIZE ( NNBTEMP, NBTEMP_FIRST, NBTEMP_LAST )

   END IF

   ! . Sum the total energies.
   EELECT = EELECT + EELECT14
   ELJ    = ELJ    + ELJ14

   END SUBROUTINE ENERGY_NON_BONDING_SELF

!-------------------------------------------------------------------------------
! . Private subroutines.
!-------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   SUBROUTINE GENERATE_REDUCED_LISTS ( SELECTION1, SELECTION2, NRECORDSI, NBLISTI_FIRST, &
                                                               NRECORDSO, NBLISTO_FIRST, &
							                  NBLISTO_LAST   )
   !--------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN)  :: NRECORDSI
   INTEGER, INTENT(OUT) :: NRECORDSO

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION1, SELECTION2

   ! . Type arguments.
   TYPE(NBLIST_TYPE), POINTER :: NBLISTI_FIRST, NBLISTO_FIRST, NBLISTO_LAST

   ! . Local scalars.
   INTEGER :: I, IINT, IREC, J, N, NSAME, NSIZE1, NSIZE2

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:) :: TEMPI

   ! . Local types.
   TYPE(NBLIST_TYPE), POINTER :: NBCURRENT, NBNEXT

   ! . Find the number of .TRUE. elements in each selection array and the number of common elements.
   NSAME  = COUNT ( SELECTION1 .AND. SELECTION2 )
   NSIZE1 = COUNT ( SELECTION1 )
   NSIZE2 = COUNT ( SELECTION2 )

   ! . Allocate temporary space for the interactions.
   ALLOCATE ( TEMPI(1:MAX(NSIZE1,NSIZE2)) )

   ! . Loop over the records in the input non-bonding interaction list.
   DO IREC = 1,NRECORDSI

      ! . Get the next list.
      IF ( IREC == 1 ) THEN
         NBCURRENT => NBLISTI_FIRST
         NBNEXT    => NBLISTI_FIRST%NEXT_LIST
      ELSE
         NBCURRENT => NBNEXT
         NBNEXT    => NBNEXT%NEXT_LIST
      END IF

      ! . Get the atom of the interaction.
      I = NBCURRENT%ATOM

      ! . Make sure that the atom is in one of the selections.
      IF ( .NOT. ( SELECTION1(I) .OR. SELECTION2(I) ) ) CYCLE

      ! . Initialize the number of interactions for the atom.
      N = 0

      ! . Loop over the interactions for the atom.
      DO IINT = 1,SIZE(NBCURRENT%INTERACTIONS)

         ! . Get the second atom of the interaction.
	 J = NBCURRENT%INTERACTIONS(IINT)

         ! . Check to see whether the interaction is to be included.
	 IF ( ( SELECTION1(I) .AND. SELECTION2(J) ) .OR. ( SELECTION1(J) .AND. SELECTION2(I) ) ) THEN

            ! . Save the interaction.
	    N = N + 1
	    TEMPI(N) = J

         END IF
      END DO

      ! . Append the interactions to the list.
      CALL NBLIST_APPEND ( N, I, TEMPI(1:N), NRECORDSO, NBLISTO_FIRST, NBLISTO_LAST )

   END DO

   ! . Deallocate TEMPI.
   DEALLOCATE ( TEMPI )

   END SUBROUTINE GENERATE_REDUCED_LISTS

   !----------------------------------------------------------------------------------------------------------------
   SUBROUTINE INTERACTIONS ( EELECT, ELJ, VIRIAL, NRECORDS, NBLIST_FIRST, CHARGE, EDEPTH, SIGMA, GRADIENT, HESSIAN )
   !----------------------------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN)    :: NRECORDS
   REAL ( KIND = DP ), INTENT(OUT)   :: EELECT, ELJ
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NATOMS),                  INTENT(IN)              :: CHARGE, EDEPTH, SIGMA
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Pointer arguments.
   TYPE(NBLIST_TYPE), POINTER :: NBLIST_FIRST

   ! . Local scalars.
   INTEGER            :: I, IFAC, IINT, INDEX, IREC, J, JFAC
   LOGICAL            :: QGRADIENT, QHESSIAN
   REAL ( KIND = DP ) :: DF, D2F, EI, EIJ, EPSFAC, ETEMP, HXX, HXY, HXZ, HYY, HYZ, HZZ, &
                         QI, QIJ, R, RIJ2, R3, R5, S, S3, S6, S12, SI, SIJ

   ! . Local cutoff scalars.
   REAL ( KIND = DP ) :: A, B, C, D, GAMMA, K6, K12, R2OFF, R2ON, SHIFT_EL1, SHIFT_EL2, SHIFT_LJ6, SHIFT_LJ12

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DR

   ! . Local types.
   TYPE(NBLIST_TYPE), POINTER :: NBCURRENT, NBNEXT

   !------------------------------------------------------------------------
   ! . Initialization.
   !------------------------------------------------------------------------
   ! . Initialize the energies.
   EELECT = 0.0_DP
   ELJ    = 0.0_DP

   ! . Check the derivative options.
   QGRADIENT = PRESENT ( GRADIENT )
   QHESSIAN  = PRESENT ( HESSIAN  )

   ! . Check the consistency of the derivative options.
   IF ( QHESSIAN .AND. .NOT. QGRADIENT ) CALL PRINT_ERROR ( "INTERACTIONS", "First derivative argument missing." )

   ! . Calculate the conversion factor for the electrostatic interactions.
   EPSFAC = ELECT_CONST / EPSILON

   ! . Calculate some constants needed for the ABFS option.
   ! . Calculate the cutoff distance squared.
   R2OFF = CUT_OFF**2
   R2ON  = CUT_ON**2

   ! . Calculate some intermediate cutoff factors.
   IF ( CUT_OFF > CUT_ON ) THEN

      ! . Calculate GAMMA, A, B, C and D.
      GAMMA = ( R2OFF - R2ON )**3
      A     = R2OFF * R2OFF * ( R2OFF - 3.0_DP * R2ON ) / GAMMA
      B     = 6.0_DP * R2OFF * R2ON / GAMMA
      C     = - ( R2OFF + R2ON ) / GAMMA
      D     = 0.4_DP / GAMMA

      ! . Calculate the electrostatic shift factors.
      SHIFT_EL1 = 8.0_DP * ( R2OFF * R2ON * ( CUT_OFF - CUT_ON ) - 0.2_DP * &
                  ( CUT_OFF * R2OFF * R2OFF - CUT_ON * R2ON * R2ON ) ) / GAMMA
      SHIFT_EL2 = - ( A / CUT_OFF ) + B * CUT_OFF + C * CUT_OFF * R2OFF + D * CUT_OFF * R2OFF * R2OFF

      ! . Calculate the Lennard-Jones shift factors.      
      K6    = ( CUT_OFF * R2OFF ) / ( CUT_OFF * R2OFF - CUT_ON * R2ON )
      K12   = ( R2OFF**3 ) / ( R2OFF**3 - R2ON**3 )

   ELSE

      ! . Initialize all values.
      GAMMA     = 1.0_DP
      A         = 0.0_DP ; B         = 0.0_DP
      C         = 0.0_DP ; D         = 0.0_DP
      K6        = 0.0_DP ; K12       = 0.0_DP
      R2OFF     = 0.0_DP ; R2ON      = 0.0_DP
      SHIFT_EL1 = 0.0_DP ; SHIFT_EL2 = 0.0_DP

   END IF

   !------------------------------------------------------------------------
   ! . Calculate the non-bonding interactions.
   !------------------------------------------------------------------------
   ! . Loop over the records in the non-bonding interaction lists.
   DO IREC = 1,NRECORDS

      ! . Get the next list.
      IF ( IREC == 1 ) THEN
         NBCURRENT => NBLIST_FIRST
         NBNEXT    => NBLIST_FIRST%NEXT_LIST
      ELSE
         NBCURRENT => NBNEXT
         NBNEXT    => NBNEXT%NEXT_LIST
      END IF

      ! . Get the atom of the interaction.
      I = NBCURRENT%ATOM

      ! . Get some information for the atom.
      EI = EDEPTH(I)
      SI = SIGMA(I)
      QI = EPSFAC * CHARGE(I)

      ! . Loop over the interactions for the atom.
      DO IINT = 1,SIZE(NBCURRENT%INTERACTIONS)

         ! . Get the second atom of the interaction.
         J = NBCURRENT%INTERACTIONS(IINT)

         ! . Get the distance between the atoms (applying the minimum image convention if necessary).
         DR = ATMCRD(1:3,I) - ATMCRD(1:3,J)
         IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
         RIJ2 = DOT_PRODUCT ( DR, DR )

         ! . Check the distance.
         IF ( RIJ2 > R2OFF ) CYCLE

         ! . Get some data for the I/J interaction.
         EIJ = EI * EDEPTH(J)
         SIJ = SI * SIGMA(J)
         QIJ = QI * CHARGE(J)

         ! . Calculate some distance factors.
         R   = SQRT ( RIJ2 )
         S   = 1.0_DP / R
         S3  = ( SIJ * S )**3
         S6  = S3 * S3

         ! . Calculate the full interaction.
         IF ( RIJ2 <= R2ON ) THEN

            ! . Calculate some intermediate quantities.
            S12        = S6 * S6
            SHIFT_LJ6  = ( ( SIJ / CUT_OFF ) * ( SIJ / CUT_ON ) )**3
            SHIFT_LJ12 = SHIFT_LJ6 * SHIFT_LJ6

            ! . Calculate the electrostatic interaction.
            ETEMP  = QIJ * S
            EELECT = EELECT + ETEMP + QIJ * SHIFT_EL1

            ! . Calculate the Lennard-Jones interaction.
            ELJ = ELJ + EIJ * ( ( S12 - SHIFT_LJ12 ) - ( S6 - SHIFT_LJ6 ) )

            ! . Check for a gradient calculation.
            IF ( .NOT. QGRADIENT ) CYCLE

            ! . Calculate an intermediate factor for the gradient calculation.
            DF  = ( - ETEMP + 6.0_DP * EIJ * ( S6 - 2.0_DP * S12 ) ) / RIJ2

         ! . Calculate the switched interaction.
         ELSE

            ! . Calculate some intermediate quantities.
            R3         = R  * RIJ2
            R5         = R3 * RIJ2
            SHIFT_LJ6  = ( SIJ / CUT_OFF )**3
            SHIFT_LJ12 = SHIFT_LJ6 * SHIFT_LJ6

            ! . Calculate the electrostatic interaction.
            EELECT = EELECT + QIJ * ( A * S - B * R - C * R3 - D * R5 + SHIFT_EL2 )

            ! . Calculate the Lennard-Jones interaction.
            ELJ = ELJ + EIJ * ( K12 * ( S6 - SHIFT_LJ12 )**2 - K6 * ( S3 - SHIFT_LJ6 )**2 )

            ! . Check for a gradient calculation.
            IF ( .NOT. QGRADIENT ) CYCLE

            ! . Calculate an intermediate factor for the gradient calculation.
            DF  = - QIJ * ( A / R3 + B * S + 3.0_DP * C * R + 5.0_DP * D * R3 ) - 6.0_DP * &
                    EIJ * ( 2.0_DP * K12 * S6 * ( S6 - SHIFT_LJ12 ) - K6 * S3 * ( S3 - SHIFT_LJ6 ) ) / RIJ2

         END IF

         ! . Calculate the contribution to the virial.
         VIRIAL = VIRIAL + DF * RIJ2

         ! . Add in the contribution to the gradients.
         GRADIENT(1:3,I) = GRADIENT(1:3,I) + DF * DR
         GRADIENT(1:3,J) = GRADIENT(1:3,J) - DF * DR

         ! . Check for a Hessian calculation.
         IF ( .NOT. QHESSIAN ) CYCLE

         ! . Calculate an intermediate factor for the Hessian calculation.
	 IF ( RIJ2 <= R2ON ) THEN
            D2F = ( 3.0_DP * ETEMP + 24.0_DP * EIJ * ( 7.0_DP * S12 - 2.0_DP * S6 ) ) / ( RIJ2 * RIJ2 )
	 ELSE
            D2F = ( QIJ * ( 2.0_DP * A / R3 - 6.0_DP * C * R - 20.0_DP * D * R3 ) - 6.0_DP * &
                    EIJ * ( 2.0_DP * K12 * S6 * ( 7.0_DP * SHIFT_LJ12 - 13.0_DP * S6 ) - &
                                     K6  * S3 * ( 4.0_DP * SHIFT_LJ6  -  7.0_DP * S3 ) ) / RIJ2 - DF ) / RIJ2
	 END IF

         ! . Calculate the elements of the hessian block.
         HXX = DR(1) * DR(1) * D2F + DF
         HXY = DR(1) * DR(2) * D2F
         HXZ = DR(1) * DR(3) * D2F
         HYY = DR(2) * DR(2) * D2F + DF
         HYZ = DR(2) * DR(3) * D2F
         HZZ = DR(3) * DR(3) * D2F + DF

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

   !------------------------------------------------------------------------------------------
   SUBROUTINE NBLIST_APPEND ( NINTS, ATOM, INTERACTIONS, NRECORDS, NBLIST_FIRST, NBLIST_LAST )
   !------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN)    :: ATOM, NINTS
   INTEGER, INTENT(INOUT) :: NRECORDS

   ! . Array arguments.
   INTEGER, DIMENSION(:), INTENT(IN) :: INTERACTIONS

   ! . Type arguments.
   TYPE(NBLIST_TYPE), POINTER :: NBLIST_FIRST, NBLIST_LAST

   ! . Local scalars.
   INTEGER :: FLAG

   ! . Local types.
   TYPE(NBLIST_TYPE), POINTER :: NLOCAL

   ! . There are no interactions.
   IF ( NINTS <= 0 ) RETURN

   ! . Allocate the new list.
   ALLOCATE ( NLOCAL )

   ! . Set the first atom of the interaction.
   NLOCAL%ATOM = ATOM

   ! . Allocate space for the second atoms of the interaction.
   ALLOCATE ( NLOCAL%INTERACTIONS(1:NINTS), STAT = FLAG )

   ! . Save the interactions for the atom.
   IF ( FLAG == 0 ) THEN
      NLOCAL%INTERACTIONS(1:NINTS) = INTERACTIONS(1:NINTS)
   ! . Space cannot be allocated for the interactions.
   ELSE
      CALL PRINT_ERROR ( "NBLIST_APPEND", "Unable to allocate non-bonding interaction lists.", ATOM )
   END IF

   ! . Nullify the pointer in NLOCAL.
   NULLIFY ( NLOCAL%NEXT_LIST )

   ! . Increment the number of lists.
   NRECORDS = NRECORDS + 1

   ! . This is the first list.
   IF ( NRECORDS == 1 ) THEN

      ! . Define the first list pointer.
      NBLIST_FIRST => NLOCAL

   ! . There are other lists.
   ELSE

      ! . Set NEXT_LIST of the old last list.
      NBLIST_LAST%NEXT_LIST => NLOCAL

   END IF

   ! . Update the last list pointer.
   NBLIST_LAST => NLOCAL

   END SUBROUTINE NBLIST_APPEND

   !-------------------------------------------------------------------
   SUBROUTINE NBLIST_INITIALIZE ( NRECORDS, NBLIST_FIRST, NBLIST_LAST )
   !-------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(INOUT) :: NRECORDS

   ! . Type arguments.
   TYPE(NBLIST_TYPE), POINTER :: NBLIST_FIRST, NBLIST_LAST

   ! . Local scalars.
   INTEGER :: I

   ! . Local types.
   TYPE(NBLIST_TYPE), POINTER :: NBCURRENT, NBNEXT

   ! . Check the number of records.
   IF ( NRECORDS < 0 ) RETURN

   ! . Loop over the records.
   DO I = 1,NRECORDS

      ! . Get the next list.
      IF ( I == 1 ) THEN
         NBCURRENT => NBLIST_FIRST
         NBNEXT    => NBLIST_FIRST%NEXT_LIST
      ELSE
         NBCURRENT => NBNEXT
         NBNEXT    => NBNEXT%NEXT_LIST
      END IF

      ! . Deallocate the interaction array.
      DEALLOCATE ( NBCURRENT%INTERACTIONS )

      ! . Nullify the pointer in NBCURRENT.
      NULLIFY ( NBCURRENT%NEXT_LIST )

      ! . Deallocate the list.
      DEALLOCATE ( NBCURRENT )

  END DO

   ! . Initialize the number of records.
   NRECORDS = -1

   ! . Nullify the first and last list pointers.
   NULLIFY ( NBLIST_FIRST, NBLIST_LAST )

   END SUBROUTINE NBLIST_INITIALIZE

   !----------------
   SUBROUTINE UPDATE
   !----------------

   ! . Local scalars.
   INTEGER            :: I, IINT, IRES, J, JRES, NINTM, NINTQ, NINTS, NTOTINTM, NTOTINTQ, QMI, QMJ
   LOGICAL            :: QOKCRD, QOKFIX
   REAL ( KIND = DP ) :: CUTSQ, RIJ2

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DR, RI, RJ, RMAX, RMIN

   ! . Local allocatable arrays.
   INTEGER,            ALLOCATABLE, DIMENSION(:)   :: TEMPI, TEMPM, TEMPQ
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: CENTER, EXTENT

   !---------------------------------------------------------------------------
   ! . Do some general checks.
   !---------------------------------------------------------------------------
   ! . There are no MM atoms.
   IF ( NATOMSMM < 0 ) RETURN

   ! . Check the consistency of the minimum image and symmetry options.
   IF ( QIMAGE ) THEN
      IF ( .NOT. QBOX ) CALL PRINT_ERROR ( "UPDATE", "A periodic box is not defined for minimum image option." )
      IF ( 2.0_DP * CUT_LIST > MINVAL ( BOXL ) ) THEN
         CALL PRINT_ERROR ( "UPDATE", "The box size is too small for the specified minimum image cutoff." )
      END IF
   ELSE
      IF ( QBOX ) CALL PRINT_ERROR ( "UPDATE", "A periodic box is defined without the minimum image option." )
   END IF

   !---------------------------------------------------------------------------
   ! . Check for changes in fixed atoms.
   !---------------------------------------------------------------------------
   ! . Check for changes in REFFIX.
   QOKFIX = ALLOCATED ( REFFIX )
   IF ( QOKFIX ) THEN
       IF ( SIZE ( REFFIX ) /= NATOMS ) THEN
          QOKFIX = .FALSE.
       ELSE IF ( .NOT. ALL ( ATMFIX .EQV. REFFIX ) ) THEN
          QOKFIX = .FALSE.
       END IF
   END IF

   ! . Reinitialize the relevant lists and arrays.
   IF ( .NOT. QOKFIX ) THEN
      CALL NBLIST_INITIALIZE ( NNBLIST14, NBLIST14_FIRST, NBLIST14_LAST )
      IF ( ALLOCATED ( REFFIX ) ) DEALLOCATE ( REFFIX )
      ALLOCATE ( REFFIX(1:NATOMS) )
      REFFIX = ATMFIX
   END IF

   !---------------------------------------------------------------------------
   ! . Update the 1,4 lists (MM/MM and QM/MM LJ) if necessary.
   !---------------------------------------------------------------------------
   ! . Check to see if an update is necessary.
   IF ( NNBLIST14 < 0 ) THEN

      ! . Initialize NNBLIST14.
      NNBLIST14 = 0

      ! . There are 1-4 interactions.
      IF ( SIZE ( ATME14J ) > 0 ) THEN

	 ! . Allocate the temporary interaction array.
	 ALLOCATE ( TEMPM(1:MIN(NATOMS-1,SIZE(ATME14J))) )

	 ! . Loop over the atoms.
	 DO I = 1,NATOMS

            ! . Initialize the number of 1-4 interactions for the atom.
	    NINTS = 0

            ! . Loop over the exclusions for the atom.
            DO IINT = (ATME14I(I)+1),ATME14I(I+1)

               ! . Get the second atom of the exclusion.
               J = ATME14J(IINT)

               ! . Skip interactions between two fixed atoms. Interactions between QM atoms (neither of which
	       ! . is boundary atoms) are also skipped but these are removed in MOPAC_HAMILTONIAN.
	       IF ( ATMFIX(I) .AND. ATMFIX(J) ) CYCLE
   
               ! . Save the interaction.
               NINTS = NINTS + 1
               TEMPM(NINTS) = J

            END DO

            ! . Append the interactions to the list.
	    CALL NBLIST_APPEND ( NINTS, I, TEMPM(1:NINTS), NNBLIST14, NBLIST14_FIRST, NBLIST14_LAST )

	 END DO

	 ! . Deallocate TEMPM.
	 DEALLOCATE ( TEMPM )

      END IF
   END IF

   !---------------------------------------------------------------------------
   ! . Check for changes in the coordinates.
   !---------------------------------------------------------------------------
   ! . The lists exist.
   IF ( ( NNBLISTMM >= 0 ) .AND. ( NNBLISTQM >= 0 ) ) THEN

      ! . Initialization.
      QOKCRD = .FALSE.

      ! . Check for valid REFCRD.
      IF ( ALLOCATED ( REFCRD ) ) THEN
         IF ( SIZE ( REFCRD, 2 ) == NATOMS ) THEN

            ! . Initialize the buffer size constant.
            CUTSQ = ( ( CUT_LIST - CUT_OFF ) / 2.0_DP )**2

            ! . Loop to see if some atoms have moved too far.
            QOKCRD = .TRUE.
            DO I = 1,NATOMS
               IF ( .NOT. ATMFIX(I) ) THEN
                  IF ( SUM ( ( ATMCRD(1:3,I) - REFCRD(1:3,I) )**2 ) > CUTSQ ) THEN
                     QOKCRD = .FALSE.
                     EXIT
                  END IF
               END IF
            END DO
         END IF
      END IF

      ! . Return if no update needs to be performed.
      IF ( QOKCRD ) RETURN

   END IF

   ! . Initialize the lists and reference arrays.
   CALL NBLIST_INITIALIZE ( NNBLISTMM, NBLISTMM_FIRST, NBLISTMM_LAST )
   CALL NBLIST_INITIALIZE ( NNBLISTQM, NBLISTQM_FIRST, NBLISTQM_LAST )
   IF ( ALLOCATED ( REFCRD ) ) DEALLOCATE ( REFCRD )

   ! . Reinitialize REFCRD.
   ALLOCATE ( REFCRD(1:3,1:NATOMS) )
   REFCRD = ATMCRD

   ! . Reset NNBLISTMM and NNBLISTQM.
   NNBLISTMM = 0
   NNBLISTQM = 0

   !---------------------------------------------------------------------------
   ! . Calculate the centers of the residues and their extents.
   !---------------------------------------------------------------------------
   ! . Allocate space for the residue information.
   ALLOCATE ( CENTER(1:3,1:NRESID), EXTENT(1:3,1:NRESID) )

   ! . Loop over the residues.
   DO IRES = 1,NRESID

      ! . Initialize the maximum and minimum arrays.
      RMAX = - HUGE ( 0.0_DP )
      RMIN =   HUGE ( 0.0_DP )

      ! . Find the maximum and minimum coordinates.
      DO I = (RESIND(IRES)+1),RESIND(IRES+1)
         RMAX = MAX ( ATMCRD(1:3,I), RMAX )
         RMIN = MIN ( ATMCRD(1:3,I), RMIN )
      END DO

      ! . Calculate the center and the extent.
      CENTER(1:3,IRES) = 0.5_DP * ( RMAX + RMIN )
      EXTENT(1:3,IRES) = 0.5_DP * ( RMAX - RMIN )

   END DO

   !---------------------------------------------------------------------------
   ! . Perform the update.
   !---------------------------------------------------------------------------
   ! . Allocate the temporary interaction arrays.
   ALLOCATE ( TEMPI(1:NATOMS), TEMPM(1:NATOMS) )
   IF ( NATOMSQM > 0 ) THEN
      ALLOCATE ( TEMPQ(1:2*NATOMS) ) ! . Depending upon the link-atom option this may need to be 2*NATOMS.
   ELSE
      ALLOCATE ( TEMPQ(1:0) )
   END IF

   ! . Initialize the cutoff constant.
   CUTSQ = CUT_LIST**2

   ! . Initialize the global number of interactions.
   NTOTINTM = 0
   NTOTINTQ = 0

   ! . Outer loop over residues.
   DO IRES = 1,NRESID

      ! . Initialize the local number of interactions.
      NINTS = 0

      ! . Inner loop over residues to determine which are in range.
      DO JRES = IRES,NRESID

         ! . Calculate the distance differences in each dimension (applying the minimum image convention).
         DR = CENTER(1:3,IRES) - CENTER(1:3,JRES)
         IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
         DR = MAX ( ( ABS ( DR ) - EXTENT(1:3,IRES) - EXTENT(1:3,JRES) ), 0.0_DP )

         ! . Calculate the distance squared.
         RIJ2 = DOT_PRODUCT ( DR, DR )

         ! . Check to see if the residue is in range.
         IF ( RIJ2 <= CUTSQ ) THEN

            ! . Loop over the atoms in JRES and fill the interaction array.
            DO J = (RESIND(JRES)+1),RESIND(JRES+1)
               NINTS = NINTS + 1
               TEMPI(NINTS) = J
            END DO

         END IF
      END DO

      ! . Loop over the atoms in the Ith residue.
      DO I = (RESIND(IRES)+1),RESIND(IRES+1)

         ! . Initialize the number of interactions for the atom.
	 NINTM = 0
	 NINTQ = 0

         ! . Set QMI to indicate whether I is MM, BA or QM.
	 IF ( ATMQMI(I) > 0 ) THEN
            IF ( QMATOM(ATMQMI(I))%QBOUNDARY ) THEN
	       QMI = 0
	       RI  = MOPAC_DATA_ATOM_POSITION ( ATMQMI(I), ATMCRD )
	    ELSE
	       QMI = 1
	    END IF
	 ELSE
	    QMI = -1
	 END IF

         ! . Loop over the atoms in the interaction array.
         DO IINT = 1,NINTS

            ! . Get the index of the second atom.
            J = TEMPI(IINT)

            ! . Check the index of J.
            IF ( I >= J ) CYCLE

            ! . Set QMJ to indicate whether J is MM, BA or QM.
	    IF ( ATMQMI(J) > 0 ) THEN
               IF ( QMATOM(ATMQMI(J))%QBOUNDARY ) THEN
		  QMJ = 0
	          RJ  = MOPAC_DATA_ATOM_POSITION ( ATMQMI(J), ATMCRD )
	       ELSE
		  QMJ = 1
	       END IF
	    ELSE
	       QMJ = -1
	    END IF

            ! . Skip interactions between two QM atoms.
            IF ( ( QMI > 0 ) .AND. ( QMJ > 0 ) ) CYCLE

            ! . Calculate the distance squared between the atoms (applying the minimum image convention).
            DR = ATMCRD(1:3,I) - ATMCRD(1:3,J)
            IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
            RIJ2 = DOT_PRODUCT ( DR, DR )

            ! . Adjust the distance squared if QM(BA)/MM interactions are to be used.
	    IF ( QUSEBAQM ) THEN
	       IF ( QMI == 0 ) THEN
        	  DR = RI - ATMCRD(1:3,J)
        	  IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
        	  RIJ2 = MIN ( DOT_PRODUCT ( DR, DR ), RIJ2 )
	       END IF
	       IF ( QMJ == 0 ) THEN
        	  DR = ATMCRD(1:3,I) - RJ
        	  IF ( QIMAGE ) DR = DR - BOXL * ANINT ( DR / BOXL, DP )
        	  RIJ2 = MIN ( DOT_PRODUCT ( DR, DR ), RIJ2 )
	       END IF
	    END IF

            ! . The interaction is in range.
            IF ( RIJ2 < CUTSQ ) THEN

               ! . There is a QM/MM electrostatic interaction (include all interactions and interactions
	       ! . between fixed atoms). A negative J indicates that J is the quantum atom.
               IF ( QMI > 0 ) THEN
                  NINTQ = NINTQ + 1
                  TEMPQ(NINTQ) = J
               ELSE IF ( QMJ > 0 ) THEN
                  NINTQ = NINTQ + 1
                  TEMPQ(NINTQ) = -J
               END IF

               ! . There is a BA(QM)/MM electrostatic interaction.
	       IF ( QUSEBAQM ) THEN
  	          IF ( ( QMI == 0 ) .AND. ( QMJ < 1 ) ) THEN
                     NINTQ = NINTQ + 1
                     TEMPQ(NINTQ) = J
		  END IF
		  IF ( ( QMJ == 0 ) .AND. ( QMI < 1 ) ) THEN
                     NINTQ = NINTQ + 1
                     TEMPQ(NINTQ) = -J
		  END IF
               END IF

               ! . Skip all excluded interactions.
	       IF ( ANY ( ATMEXCJ(ATMEXCI(I)+1:ATMEXCI(I+1)) == J ) ) CYCLE

               ! . Skip MM/MM and QM(BA and LJ)/MM interactions between fixed atoms.
               IF ( ATMFIX(I) .AND. ATMFIX(J) ) CYCLE

               ! . There is a MM/MM or QM/MM (LJ) interaction (note that there is never a BA(QM)/MM LJ interaction).
               NINTM = NINTM + 1
               TEMPM(NINTM) = J

            END IF
         END DO

         ! . Include the BA(QM)/BA(MM) self-interaction if necessary.
	 IF ( QUSEBAQM .AND. ( QMI == 0 ) ) THEN
	    NINTQ = NINTQ + 1
	    TEMPQ(NINTQ) = I
	 END IF

         ! . Append the interactions to the lists.
	 CALL NBLIST_APPEND ( NINTM, I, TEMPM(1:NINTM), NNBLISTMM, NBLISTMM_FIRST, NBLISTMM_LAST )
	 CALL NBLIST_APPEND ( NINTQ, I, TEMPQ(1:NINTQ), NNBLISTQM, NBLISTQM_FIRST, NBLISTQM_LAST )

         ! . Increment the total number of interactions.
	 NTOTINTM = NTOTINTM + NINTM
	 NTOTINTQ = NTOTINTQ + NINTQ

      END DO
   END DO

   !---------------------------------------------------------------------------
   ! . Finish up.
   !---------------------------------------------------------------------------
   ! . Deallocate the temporary arrays.
   DEALLOCATE ( CENTER, EXTENT, TEMPI, TEMPM, TEMPQ )

   ! . Print out some information if required.
   IF ( QNBPRINT ) THEN
      CALL PRINT_PARAGRAPH_START
      WRITE ( PRINT_LINE, "(A,I10,A)" ) "MM/MM non-bonding interaction lists created with ", NTOTINTM, " interactions."
      CALL PRINT_TEXT
      IF ( NATOMSQM > 0 ) THEN
         CALL PRINT_LINEBREAK
	 WRITE ( PRINT_LINE, "(A,I10,A)" ) "QM/MM non-bonding interaction lists created with ", NTOTINTQ, " interactions."
         CALL PRINT_TEXT
      END IF
      CALL PRINT_PARAGRAPH_STOP
   END IF

   ! . Increment the NUPDATES and the SUM_INTM and SUM_INTQ counters.
   NUPDATES = NUPDATES + 1
   SUM_INTM = SUM_INTM + REAL ( NTOTINTM, DP )
   SUM_INTQ = SUM_INTQ + REAL ( NTOTINTQ, DP )

   ! . Save the total number of QM/MM interactions.
   NQMINT = NTOTINTQ

   END SUBROUTINE UPDATE

END MODULE ENERGY_NON_BONDING_ABFS
