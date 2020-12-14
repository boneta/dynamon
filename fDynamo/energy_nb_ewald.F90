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
!              The Non-Bonding Energy Module (Ewald - OPLS)
!========================================================================
!
! . Module scalars:
!
!   BETA                           The beta/kappa value for the Ewald procedure.
!   CUT_LIST                       The cutoff for inclusion on the list.
!   CUT_OFF                        The cutoff for the interactions.
!   EPSILON                        The dielectric constant.
!   NCALLS                         The number of calls.
!   NKVEC                          The number of k-vectors.
!   NUPDATESK                      The number of k-vector updates.
!   NUPDATESL                      The number of list updates.
!   QPRINT                         The update printing flag.
!   QTINFOIL                       The tinfoil boundary condition flag.
!   SUM_INTS                       The accumulator for the total number
!                                  of interactions.
!   SUM_VEC                        The accumulator for the total number
!                                  of k-vectors.
!
! . Module arrays:
!
!   ECOMPEL, ECOMPLJ               The electrostatic and LJ energy components
!                                  (direct, reciprocal, self and tinfoil).
!   KMAX                           The maximum value of the wave-vector
!                                  in each direction (x,y,z).
!   LISTI, LISTJ                   The atom non-bond lists.
!   REFCRD                         The reference atom coordinates.
!   RKEXPEL, RKEXPLJ               The prefactors for the reciprocal
!                                  space calculation.
!   RKVEC                          The wave-vectors ( 2 Pi k / L ).
!
! . Subroutines (public):
!
!   ENERGY_NON_BONDING_CALCULATE   Calculate the non-bonding energy and
!                                  its derivatives.
!   ENERGY_NON_BONDING_OPTIONS     Set the non-bonding options.
!   ENERGY_NON_BONDING_STATISTICS  Print out some non-bonding statistics.
!
! . Subroutines (public but not implemented):
!
!   ENERGY_NON_BONDING_INTERACTION
!   ENERGY_NON_BONDING_SELF
!
! . Subroutines (internal to ENERGY_NON_BONDING_CALCULATE):
!
!   SELF_ENERGIES                  Calculate the self-energies.
!   UPDATE                         Update the non-bond lists.
!
! . Notes:
!
!   Calculate the non-bonding energy using a standard Ewald summation
!   method for both the electrostatic and Lennard-Jones energies. The
!   non-bonding interaction lists are generated and updated automatically
!   using a residue-based search scheme.
!
!   Other points:
!
!   1. "Tin-foil" boundary conditions are assumed.
!
!   2. Only the LJ dispersion energies are calculated using the Ewald
!      scheme. The repulsive terms are calculated directly using the
!      list (i.e. it is assumed they are sufficiently short-range).
!
!   3. The excluded interactions are treated in the way suggested by
!      Essmann et al. The 1,2, 1,3 and 1,4 interactions are skipped in
!      the direct sum and their contributions are subtracted separately.
!      The 1,4 terms are then calculated as normal.
!
!   4. The following features are not implemented at present: (i) the
!      calculation of the virial; (ii) second derivatives; (iii) use of
!      fixed atom information; (iv) quantum atoms; (v) interaction
!      energies.
!
!========================================================================
MODULE ENERGY_NON_BONDING !_EWALD

! . Module declarations.
USE CONSTANTS,         ONLY : ELECT_CONST, PI
USE DEFINITIONS,       ONLY : DP
USE PRINTING,          ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, &
                              PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, PRINT_SUMMARY_STOP
USE SPECIAL_FUNCTIONS, ONLY : ERFC

USE ATOMS,             ONLY : ATMCRD, NATOMS, NATOMSQM, NFREE
USE MM_TERMS
USE SEQUENCE,          ONLY : NRESID, RESIND
USE SYMMETRY,          ONLY : BOXL, QBOX

IMPLICIT NONE
PRIVATE
PUBLIC :: ENERGY_NON_BONDING_CALCULATE, ENERGY_NON_BONDING_INTERACTION, ENERGY_NON_BONDING_OPTIONS, &
          ENERGY_NON_BONDING_SELF, ENERGY_NON_BONDING_STATISTICS, ECOMPEL, ECOMPLJ,                 &
          CUT_OFF, CUT_ON, NQMINT, QIMAGE, NNBLISTQM, NBLISTQM_FIRST, NBLISTQM_LAST
#ifndef PGPC
SAVE
#endif

! . Module scalars.
INTEGER            :: NCALLS    =         0, NKVEC     =         0, &
                      NUPDATESK =         0, NUPDATESL =         0
LOGICAL            :: QPRINT    =   .FALSE., QTINFOIL  =    .TRUE.
REAL ( KIND = DP ) :: BETA      =    0.2_DP, CUT_LIST  = 9999.0_DP, &
                      CUT_OFF   = 9999.0_DP, EPSILON   =    1.0_DP, &
                      SUM_INTS  =    0.0_DP, SUM_VEC   =    0.0_DP

! . Module arrays.
INTEGER,            DIMENSION(1:3) :: KMAX    = 0
REAL ( KIND = DP ), DIMENSION(1:3) :: BOXOLD  = 0.0_DP
REAL ( KIND = DP ), DIMENSION(1:4) :: ECOMPEL = 0.0_DP, ECOMPLJ = 0.0_DP

! . Allocatable module arrays.
INTEGER,            ALLOCATABLE, DIMENSION(:)   :: LISTI,   LISTJ
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: RKEXPEL, RKEXPLJ
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: REFCRD,  RKVEC

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
   REAL ( KIND = DP ), INTENT(OUT) :: EELECT, ELJ, VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),                INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NATOMS*(3*NATOMS+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Local scalars.
   LOGICAL :: QGRADIENT, QHESSIAN

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NATOMS) :: BII

   !---------------------------------------------------------------------
   ! . Initialization.
   !---------------------------------------------------------------------
   ! . Initialize the energy values.
   EELECT = 0.0_DP
   ELJ    = 0.0_DP

   ! . Set the virial to a very large number!
   VIRIAL = HUGE ( 0.0_DP )

   ! . Check for no calculation.
   IF ( NFREE <= 0 ) RETURN

   ! . Check for the presence of a periodic box.
   IF ( .NOT. QBOX ) CALL PRINT_ERROR ( "ENERGY_NON_BONDING", "A periodic box has not been defined." )

   ! . Check the derivative arrays.
   QGRADIENT = PRESENT ( GRADIENT )
   QHESSIAN  = PRESENT ( HESSIAN  )

   ! . Ewald second derivatives cannot be calculated.
   IF ( QHESSIAN ) CALL PRINT_ERROR ( "ENERGY_NON_BONDING", "Cannot calculate Ewald second derivatives." )

   ! . Check for quantum atoms.
   IF ( NATOMSQM > 0 ) CALL PRINT_ERROR ( "ENERGY_NON_BONDING", "This module cannot handle quantum atoms." )

   ! . Initialize the energy components.
   ECOMPEL = 0.0_DP
   ECOMPLJ = 0.0_DP

   !---------------------------------------------------------------------
   ! . Calculate the real space interactions.
   !---------------------------------------------------------------------
   ! . Update the non-bond lists if necessary.
   CALL UPDATE

   ! . Calculate the direct interactions (minus the exclusions).
   IF ( SIZE ( LISTJ ) > 0 ) CALL DIRECT

   ! . Remove the contributions to the real space sum from the exclusions.
   IF ( SIZE ( ATMEXCJ ) > 0 ) CALL EXCLUSIONS

   ! . Calculate the full 1,4 interactions.
   IF ( SIZE ( ATME14J ) > 0 ) CALL INTERACTIONS14

   !---------------------------------------------------------------------
   ! . Fill the Bii array for the dispersion interactions.
   !---------------------------------------------------------------------
   ! . Calculate the array.
   BII = ATMEPS * ATMSIG**6

   !---------------------------------------------------------------------
   ! . Calculate the reciprocal space terms.
   !---------------------------------------------------------------------
   ! . Calculate the k-vectors if necessary.
   CALL KVECTORS

   ! . Calculate the reciprocal space terms.
   CALL RECIPROCAL

   !---------------------------------------------------------------------
   ! . Calculate the self-energies.
   !---------------------------------------------------------------------
   CALL SELF_ENERGIES

   !---------------------------------------------------------------------
   ! . Calculate the tin-foil boundary condition corrections.
   !---------------------------------------------------------------------
   CALL TINFOIL_ENERGIES

   ! . Increment the NCALLS counter.
   NCALLS = NCALLS + 1

   ! . Sum the energies.
   EELECT = SUM ( ECOMPEL )
   ELJ    = SUM ( ECOMPLJ )

   !=====================================================================
   CONTAINS
   !=====================================================================

      !----------------
      SUBROUTINE DIRECT
      !----------------

      ! . Local scalars.
      INTEGER            :: I, IINT, J
      REAL ( KIND = DP ) :: ERFFAC, EI, EIJ, EPSFAC, G6, QI, QIJ, R, RIJ2, R2OFF, S, S6, S12, SI, SIJ, X2

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      ! . Calculate the conversion factor for the electrostatic interactions.
      EPSFAC = ELECT_CONST / EPSILON

      ! . Calculate the cutoff distance squared.
      R2OFF = CUT_OFF**2

      ! . Outer loop over the atoms.
      DO I = 1,NATOMS

         ! . Get some information for the atom.
         EI = ATMEPS(I)
         SI = ATMSIG(I)
         QI = EPSFAC * ATMCHG(I)

         ! . Loop over the interactions for the atom.
         DO IINT = (LISTI(I)+1),LISTI(I+1)

            ! . Get the second atom of the interaction.
            J = LISTJ(IINT)

            ! . Get the distance between the atoms (applying the minimum image convention).
            DR = ATMCRD(1:3,I) - ATMCRD(1:3,J)
            DR = DR - BOXL * ANINT ( DR / BOXL, DP )
            RIJ2 = DOT_PRODUCT ( DR, DR )

            ! . Check the distance.
            IF ( RIJ2 > R2OFF ) CYCLE

            ! . Get some data for the I/J interaction.
            EIJ = EI * ATMEPS(J)
            SIJ = SI * ATMSIG(J)
            QIJ = QI * ATMCHG(J)

            ! . Calculate some distance factors.
            R = SQRT ( RIJ2 )
            S = 1.0_DP / R

            ! . Get ERFC.
            ERFFAC = ERFC ( BETA * R )

            ! . Calculate the electrostatic interaction.
            ECOMPEL(1) = ECOMPEL(1) + QIJ * ERFFAC * S

            ! . Calculate some factors for the LJ calculation.
            S6  = ( SIJ * S )**6
            S12 = S6 * S6

            ! . Calculate G6.
            X2 = BETA * BETA * RIJ2
            G6 = ( 1.0_DP + X2 + 0.5_DP * X2 * X2 ) * EXP ( - X2 )

            ! . Calculate the Lennard-Jones interaction.
            ECOMPLJ(1) = ECOMPLJ(1) + EIJ * ( S12 - S6 * G6 )

            ! . Calculate the derivatives.
            IF ( QGRADIENT ) THEN

               ! . Calculate some intermediate factors for the gradient calculation.
               DR = DR * ( - QIJ * ( ERFFAC * S + 2.0_DP * BETA * EXP ( - X2 ) / SQRT ( PI ) ) + &
                             EIJ * S6 * ( 6.0_DP * G6 + X2**3 * EXP ( -X2 ) - 12.0_DP * S6 ) ) / RIJ2

               ! . Calculate the gradients.
               GRADIENT(1:3,I) = GRADIENT(1:3,I) + DR
               GRADIENT(1:3,J) = GRADIENT(1:3,J) - DR

            END IF
         END DO
      END DO

      END SUBROUTINE DIRECT

      !--------------------
      SUBROUTINE EXCLUSIONS
      !--------------------

      ! . Local scalars.
      INTEGER            :: I, IINT, J
      REAL ( KIND = DP ) :: ERFFAC, EI, EIJ, EPSFAC, G6, QI, QIJ, R, RIJ2, S, S6, SI, SIJ, X2

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      ! . Calculate the conversion factor for the electrostatic interactions.
      EPSFAC = ELECT_CONST / EPSILON

      ! . Outer loop over the atoms.
      DO I = 1,NATOMS

         ! . Get some information for the atom.
         EI = ATMEPS(I)
         SI = ATMSIG(I)
         QI = EPSFAC * ATMCHG(I)

         ! . Loop over the interactions for the atom.
         DO IINT = (ATMEXCI(I)+1),ATMEXCI(I+1)

            ! . Get the second atom of the interaction.
            J = ATMEXCJ(IINT)

            ! . Get the distance between the atoms (applying the minimum image convention).
            DR = ATMCRD(1:3,I) - ATMCRD(1:3,J)
            DR = DR - BOXL * ANINT ( DR / BOXL, DP )
            RIJ2 = DOT_PRODUCT ( DR, DR )

            ! . Get some data for the I/J interaction.
            EIJ = EI * ATMEPS(J)
            SIJ = SI * ATMSIG(J)
            QIJ = QI * ATMCHG(J)

            ! . Calculate some distance factors.
            R = SQRT ( RIJ2 )
            S = 1.0_DP / R

            ! . Get ERF.
            ERFFAC = ERFC ( BETA * R ) - 1.0_DP

            ! . Calculate the electrostatic interaction.
            ECOMPEL(1) = ECOMPEL(1) + QIJ * ERFFAC * S

            ! . Calculate some factors for the LJ calculation (but leave out S12!).
            S6  = ( SIJ * S )**6

            ! . Calculate G6.
            X2 = BETA * BETA * RIJ2
            G6 = ( 1.0_DP + X2 + 0.5_DP * X2 * X2 ) * EXP ( - X2 ) - 1.0_DP

            ! . Calculate the Lennard-Jones interaction.
            ECOMPLJ(1) = ECOMPLJ(1) - EIJ * S6 * G6

            ! . Calculate the derivatives.
            IF ( QGRADIENT ) THEN

               ! . Calculate some intermediate factors for the gradient calculation.
               DR = DR * ( - QIJ * ( ERFFAC * S + 2.0_DP * BETA * EXP ( - X2 ) / SQRT ( PI ) ) + &
                           EIJ * S6 * ( 6.0_DP * G6 + X2**3 * EXP ( -X2 ) ) ) / RIJ2

               ! . Calculate the gradients.
               GRADIENT(1:3,I) = GRADIENT(1:3,I) + DR
               GRADIENT(1:3,J) = GRADIENT(1:3,J) - DR

            END IF
         END DO
      END DO

      END SUBROUTINE EXCLUSIONS

      !------------------------
      SUBROUTINE INTERACTIONS14
      !------------------------

      ! . Local scalars.
      INTEGER            :: I, IINT, J
      REAL ( KIND = DP ) :: EI, EIJ, EPSFAC, QI, QIJ, R, RIJ2, S, S6, S12, SI, SIJ

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      ! . Calculate the conversion factor for the electrostatic interactions.
      EPSFAC = ELECT_CONST / EPSILON

      ! . Outer loop over the atoms.
      DO I = 1,NATOMS

         ! . Get some information for the atom.
         EI = ATMEPS14(I)
         SI = ATMSIG(I)
         QI = EPSFAC * ATMCHG14(I)

         ! . Loop over the interactions for the atom.
         DO IINT = (ATME14I(I)+1),ATME14I(I+1)

            ! . Get the second atom of the interaction.
            J = ATME14J(IINT)

            ! . Get the distance between the atoms (applying the minimum image convention).
            DR = ATMCRD(1:3,I) - ATMCRD(1:3,J)
            DR = DR - BOXL * ANINT ( DR / BOXL, DP )
            RIJ2 = DOT_PRODUCT ( DR, DR )

            ! . Get some data for the I/J interaction.
            EIJ = EI * ATMEPS14(J)
            SIJ = SI * ATMSIG(J)
            QIJ = QI * ATMCHG14(J)

            ! . Calculate some distance factors.
            R = SQRT ( RIJ2 )
            S = 1.0_DP / R

            ! . Calculate the electrostatic interaction.
            ECOMPEL(1) = ECOMPEL(1) + QIJ * S

            ! . Calculate some factors for the LJ calculation.
            S6  = ( SIJ * S )**6
            S12 = S6 * S6

            ! . Calculate the Lennard-Jones interaction.
            ECOMPLJ(1) = ECOMPLJ(1) + EIJ * ( S12 - S6 )

            ! . Calculate the derivatives.
            IF ( QGRADIENT ) THEN

               ! . Calculate some intermediate factors for the gradient calculation.
               DR = DR * ( - QIJ * S + 6.0_DP * EIJ * ( S6 - 2.0_DP * S12 ) ) / RIJ2

               ! . Calculate the gradients.
               GRADIENT(1:3,I) = GRADIENT(1:3,I) + DR
               GRADIENT(1:3,J) = GRADIENT(1:3,J) - DR

            END IF
         END DO
      END DO

      END SUBROUTINE INTERACTIONS14

      !------------------
      SUBROUTINE KVECTORS
      !------------------

      ! . Local scalars.
      INTEGER            :: I, IX, IY, IZ, N
      REAL ( KIND = DP ) :: EXPFAC, K2, X

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: KTMP

      ! . Return if the wave-vectors do not need to be calculated.
      IF ( ALL ( BOXL == BOXOLD ) .AND. ( NKVEC > 0 ) ) RETURN

      ! . Calculate the number of k-vectors (without the central box and taking into account inverses).
      NKVEC = ( PRODUCT ( 2 * KMAX + 1 ) - 1 ) / 2

      ! . Deallocate the wave-vector arrays.
      IF ( ALLOCATED ( RKEXPEL ) ) DEALLOCATE ( RKEXPEL )
      IF ( ALLOCATED ( RKEXPLJ ) ) DEALLOCATE ( RKEXPLJ )
      IF ( ALLOCATED ( RKVEC   ) ) DEALLOCATE ( RKVEC   )

      ! . Allocate the wave-vector arrays.
      ALLOCATE ( RKEXPEL(1:NKVEC), RKEXPLJ(1:NKVEC), RKVEC(1:3,1:NKVEC) )

      ! . Loop over the k-vectors.
      N = 0
      DO IZ = -KMAX(3),KMAX(3)
         DO IY = -KMAX(2),KMAX(2)
            DO IX = -KMAX(1),KMAX(1)

               ! . Skip the central box.
               IF ( ( IX == 0 ) .AND. ( IY == 0 ) .AND. ( IZ == 0 ) ) CYCLE

               ! . Save the vector.
               KTMP(1) = REAL ( IX, DP )
               KTMP(2) = REAL ( IY, DP )
               KTMP(3) = REAL ( IZ, DP )

               ! . Check the previous k-vectors for the inverse.
               DO I = 1,N
                  IF ( ALL ( KTMP == - RKVEC(:,I) ) ) GOTO 10
               END DO

               ! . Save the vector.
               N = N + 1
               RKVEC(:,N) = KTMP

               ! . End of the loop.
               10 CONTINUE

            END DO
         END DO
      END DO

      ! . Loop over the k-vectors.
      DO I = 1,NKVEC

         ! . Scale RKVEC.
         RKVEC(:,I) = 2.0_DP * PI * RKVEC(:,I) / BOXL

         ! . Calculate K2 and X.
         K2 = DOT_PRODUCT ( RKVEC(:,I), RKVEC(:,I) )
         X  = SQRT ( K2 ) / ( 2.0_DP * BETA )

         ! . Calculate some additional factors.
         EXPFAC = EXP ( - X * X )

         ! . Calculate the electrostatic and LJ prefactors.
         RKEXPEL(I) = EXPFAC / K2
         RKEXPLJ(I) = ( 0.5_DP - X * X ) * EXPFAC + SQRT ( PI ) * X**3 * ERFC ( X )

      END DO

      ! . Scale RKEXPEL and RKEXPLJ (with a factor of 2 for inverses).
      RKEXPEL =   2.0_DP * ( ELECT_CONST / EPSILON ) * ( 2.0_DP * PI / PRODUCT ( BOXL ) ) * RKEXPEL
      RKEXPLJ = - 2.0_DP * ( SQRT ( PI**3 ) * BETA**3 ) / ( 3.0_DP * PRODUCT ( BOXL ) ) * RKEXPLJ  

      ! . Save the box-size.
      BOXOLD = BOXL

      ! . Print out some information if required.
      IF ( QPRINT ) THEN
         WRITE ( PRINT_LINE, "(A,I10)" ) "Number of reciprocal space k-vectors generated = ", NKVEC, "."
	 CALL PRINT_PARAGRAPH
      END IF

      ! . Increment the NUPDATESK and the SUM_VEC counters.
      NUPDATESK = NUPDATESK + 1
      SUM_VEC   = SUM_VEC + REAL ( NKVEC, DP )

      END SUBROUTINE KVECTORS

      !--------------------
      SUBROUTINE RECIPROCAL
      !--------------------

      ! . Local scalars.
      INTEGER            :: I, IKVEC
      REAL ( KIND = DP ) :: COSKR, KR, SINKR

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3)     :: SUMEL, SUMLJ
      REAL ( KIND = DP ), DIMENSION(1:NKVEC) :: COSEL, COSLJ, SINEL, SINLJ

      ! . Initialization.
      COSEL = 0.0_DP ; SINEL = 0.0_DP
      COSLJ = 0.0_DP ; SINLJ = 0.0_DP

      ! . Calculate the electrostatic and Lennard-Jones cosine and sine factors.
      DO IKVEC = 1,NKVEC

         ! . Loop over the atoms.
         DO I = 1,NATOMS

            ! . Calculate the dot product of the atom position and the k-vector.
            KR = DOT_PRODUCT ( RKVEC(:,IKVEC), ATMCRD(:,I) )

            ! . Calculate the cos and sin.
            COSKR = COS ( KR )
            SINKR = SIN ( KR )

            ! . Calculate the contributions to the cosine and sine factors.
            COSEL(IKVEC) = COSEL(IKVEC) + ATMCHG(I) * COSKR
            SINEL(IKVEC) = SINEL(IKVEC) + ATMCHG(I) * SINKR
            COSLJ(IKVEC) = COSLJ(IKVEC) + BII(I)    * COSKR
            SINLJ(IKVEC) = SINLJ(IKVEC) + BII(I)    * SINKR

         END DO
      END DO

      ! . Calculate the energy of the reciprocal space term.
      ECOMPEL(2) = DOT_PRODUCT ( RKEXPEL, ( COSEL * COSEL + SINEL * SINEL ) )
      ECOMPLJ(2) = DOT_PRODUCT ( RKEXPLJ, ( COSLJ * COSLJ + SINLJ * SINLJ ) )

      ! . Check for calculation of the gradient.
      IF ( QGRADIENT ) THEN

         ! . Calculate the derivatives for the atoms.
         DO I = 1,NATOMS

            ! . Initialization.
            SUMEL = 0.0_DP
            SUMLJ = 0.0_DP

            ! . Loop over the k-vectors.
            DO IKVEC = 1,NKVEC

               ! . Calculate the dot product of the atom position and the k-vector.
               KR = DOT_PRODUCT ( RKVEC(:,IKVEC), ATMCRD(:,I) )

               ! . Calculate the cos and sin.
               COSKR = COS ( KR )
               SINKR = SIN ( KR )

               ! . Calculate the contribution to the electrostatic gradient.
               SUMEL = SUMEL + RKEXPEL(IKVEC) * ( - COSEL(IKVEC) * SINKR + SINEL(IKVEC) * COSKR ) * RKVEC(:,IKVEC)
               SUMLJ = SUMLJ + RKEXPLJ(IKVEC) * ( - COSLJ(IKVEC) * SINKR + SINLJ(IKVEC) * COSKR ) * RKVEC(:,IKVEC)

            END DO

            ! . Add in the contribution to the gradient.
            GRADIENT(:,I) = GRADIENT(:,I) + 2.0_DP * ( ATMCHG(I) * SUMEL + BII(I) * SUMLJ )

         END DO
      END IF

      END SUBROUTINE RECIPROCAL

      !-----------------------
      SUBROUTINE SELF_ENERGIES
      !-----------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: FAC1, FAC2

      ! . Calculate the electrostatic self-energy.
      ECOMPEL(3) = - ( ELECT_CONST / EPSILON ) * ( BETA / SQRT ( PI ) ) * DOT_PRODUCT ( ATMCHG, ATMCHG )

      ! . Calculate some factors the LJ self-energy calculation.
      FAC1 = ( SQRT ( PI ) * BETA )**3 / ( 6.0_DP * PRODUCT ( BOXL ) )
      FAC2 = BETA**6 / 12.0_DP

      ! . Calculate the LJ self-energy.
      ECOMPLJ(3) = - ( FAC1 * ( SUM ( BII ) )**2 - FAC2 * DOT_PRODUCT ( BII, BII ) )

      END SUBROUTINE SELF_ENERGIES

      !--------------------------
      SUBROUTINE TINFOIL_ENERGIES
      !--------------------------

      ! . Local scalars.
      INTEGER            :: I
      REAL ( KIND = DP ) :: PREFAC

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DIPOLE

      ! . There are no corrections for tinfoil boundary conditions.
      IF ( QTINFOIL ) RETURN

      ! . Calculate the total dipole for the cell.
      DO I = 1,3
         DIPOLE(I) = SUM ( ATMCHG * ATMCRD(I,:) )
      END DO

      ! . Calculate some intermediate constants.
      PREFAC = ( ELECT_CONST / EPSILON ) * ( 2.0_DP * PI / ( 3.0_DP * PRODUCT ( BOXL ) ) )

      ! . Calculate the electrostatic correction (assume a sphere).
      ECOMPEL(4) = PREFAC * DOT_PRODUCT ( DIPOLE, DIPOLE )

      ! . Calculate the correction to the gradients.
      IF ( QGRADIENT ) THEN

         ! . Loop over the atoms.
         DO I = 1,NATOMS
            GRADIENT(1:3,I) = GRADIENT(1:3,I) + 2.0_DP * PREFAC * ATMCHG(I) * DIPOLE
         END DO

      END IF

      END SUBROUTINE TINFOIL_ENERGIES

      !----------------
      SUBROUTINE UPDATE
      !----------------

      ! . Local scalars.
      INTEGER            :: I, IINT, IRES, J, JRES, MAXINT, NINTL, NINTR
      LOGICAL            :: QOK
      REAL ( KIND = DP ) :: CUTSQ, RIJ2

      ! . Local arrays.
      INTEGER,            DIMENSION(1:NATOMS)     :: INDEX
      REAL ( KIND = DP ), DIMENSION(1:3)          :: DR, RMAX, RMIN
      REAL ( KIND = DP ), DIMENSION(1:3,1:NRESID) :: CENTER, EXTENT

      ! . Local allocatable arrays.
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TEMPJ

      !---------------------------------------------------------------------------
      ! . Check to see if an update is necessary.
      !---------------------------------------------------------------------------
      ! . Check the value of the list cutoff against the box-size (as the
      ! . minimum image convention must be satisfied).
      IF ( 2.0_DP * CUT_LIST > MINVAL ( BOXL ) ) THEN
         CALL PRINT_ERROR ( "UPDATE", "The box size is too small for the specified list cutoff." )
      END IF

      ! . Check to see if all the list arrays are allocated.
      IF ( ALLOCATED ( REFCRD ) .AND. ALLOCATED ( LISTI ) .AND. ALLOCATED ( LISTJ ) ) THEN

         ! . Initialize the buffer size constant.
         CUTSQ = ( ( CUT_LIST - CUT_OFF ) / 2.0_DP )**2

         ! . Determine if any atoms have moved by a given amount.
         QOK = .TRUE.
         DO I = 1,NATOMS
            IF ( SUM ( ( ATMCRD(1:3,I) - REFCRD(1:3,I) )**2 ) > CUTSQ ) THEN
               QOK = .FALSE.
               EXIT
            END IF
         END DO

      ! . An update needs to be done.
      ELSE

         ! . Set the update flag.
         QOK = .FALSE.

      END IF

      ! . Return if no update needs to be performed.
      IF ( QOK ) RETURN

      ! . Deallocate the module arrays.
      IF ( ALLOCATED ( REFCRD ) ) DEALLOCATE ( REFCRD )
      IF ( ALLOCATED ( LISTI  ) ) DEALLOCATE ( LISTI  )
      IF ( ALLOCATED ( LISTJ  ) ) DEALLOCATE ( LISTJ  )

      ! . Reallocate some of the module arrays.
      ALLOCATE ( REFCRD(1:3,1:NATOMS), LISTI(1:NATOMS+1) )

      ! . Save the current atom coordinates in the reference coordinate set.
      REFCRD = ATMCRD

      !---------------------------------------------------------------------------
      ! . Calculate the centers of the residues and their extents.
      !---------------------------------------------------------------------------
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
      ! . Initialize the maximum number of interactions.
      MAXINT = 40 * NATOMS

      ! . Allocate space for the temporary LISTJ array.
      ALLOCATE ( TEMPJ(1:MAXINT) )

      ! . Initialize the cutoff constant.
      CUTSQ = CUT_LIST**2

      ! . Initialize the global number of interactions.
      NINTR = 0

      ! . Outer loop over residues.
      DO IRES = 1,NRESID

         ! . Initialize the local number of interactions.
         NINTL = 0

         ! . Inner loop over residues to determine which are in range.
         DO JRES = IRES,NRESID

            ! . Calculate the distance differences in each dimension (applying the minimum image convention).
            DR = CENTER(1:3,IRES) - CENTER(1:3,JRES)
            DR = DR - BOXL * ANINT ( DR / BOXL, DP )
            DR = MAX ( ( ABS ( DR ) - EXTENT(1:3,IRES) - EXTENT(1:3,JRES) ), 0.0_DP )

            ! . Calculate the distance squared.
            RIJ2 = DOT_PRODUCT ( DR, DR )

            ! . Check to see if the residue is in range.
            IF ( RIJ2 <= CUTSQ ) THEN

               ! . Loop over the atoms in JRES and fill the interaction array.
               DO J = (RESIND(JRES)+1),RESIND(JRES+1)
                  NINTL = NINTL + 1
                  INDEX(NINTL) = J
               END DO

            END IF
         END DO

         ! . Check the size of the temporary interaction array.
         IF ( ( NINTR + NINTL * ( RESIND(IRES+1) - RESIND(IRES) ) ) > MAXINT ) THEN

            ! . Save the existing lists.
            ALLOCATE ( LISTJ(1:NINTR) ) ; LISTJ(1:NINTR) = TEMPJ(1:NINTR)

            ! . Increment MAXINT.
            MAXINT = MAXINT + 20 * NATOMS

            ! . Reallocate the temporary lists and resave the data.
            DEALLOCATE ( TEMPJ ) ; ALLOCATE ( TEMPJ(1:MAXINT) ) ; TEMPJ(1:NINTR) = LISTJ(1:NINTR)

            ! . Deallocate LISTJ.
            DEALLOCATE ( LISTJ )

         END IF

         ! . Loop over the atoms in the Ith residue.
         DO I = (RESIND(IRES)+1),RESIND(IRES+1)

            ! . Fill the element of the LISTI array for the atom.
            LISTI(I) = NINTR

            ! . Loop over the atoms in the interaction array.
            DO IINT = 1,NINTL

               ! . Get the index of the second atom.
               J = INDEX(IINT)

               ! . Check the index of J.
               IF ( I >= J ) CYCLE

               ! . Check to see if J is excluded.
               IF ( ANY ( ATMEXCJ(ATMEXCI(I)+1:ATMEXCI(I+1)) == J ) ) CYCLE

               ! . Calculate the distance squared between the atoms (applying the minimum image convention).
               DR   = ATMCRD(1:3,I) - ATMCRD(1:3,J)
               DR   = DR - BOXL * ANINT ( DR / BOXL, DP )
               RIJ2 = DOT_PRODUCT ( DR, DR )

               ! . Add the interaction to the list.
               IF ( RIJ2 < CUTSQ ) THEN
                  NINTR = NINTR + 1
                  TEMPJ(NINTR) = J
               END IF
            END DO
         END DO
      END DO

      ! . Fill the last LISTI element.
      LISTI(NATOMS+1) = NINTR

      !---------------------------------------------------------------------------
      ! . Finish up.
      !---------------------------------------------------------------------------
      ! . Allocate the LISTJ array.
      ALLOCATE ( LISTJ(1:NINTR) )

      ! . Save the connectivity data.
      LISTJ(1:NINTR) = TEMPJ(1:NINTR)

      ! . Deallocate the remaining temporary arrays.
      DEALLOCATE ( TEMPJ )

      ! . Print out some information if required.
      IF ( QPRINT ) THEN
         WRITE ( PRINT_LINE, "(A,I10,A)" ) "Non-bonding interaction lists created with ", NINTR, " interactions."
         CALL PRINT_PARAGRAPH
      END IF

      ! . Increment the NUPDATESL and the SUM_INTS counters.
      NUPDATESL = NUPDATESL + 1
      SUM_INTS  = SUM_INTS + REAL ( NINTR, DP )

      END SUBROUTINE UPDATE

   END SUBROUTINE ENERGY_NON_BONDING_CALCULATE

   !---------------------------------------------------------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_OPTIONS ( LIST_CUTOFF, OUTER_CUTOFF, KAPPA, DIELECTRIC, KMAXIMUM, TINFOIL, PRINT )
   !---------------------------------------------------------------------------------------------------------------

   ! . Optional scalar arguments.
   LOGICAL,            INTENT(IN), OPTIONAL :: PRINT, TINFOIL
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: LIST_CUTOFF, OUTER_CUTOFF, KAPPA, DIELECTRIC

   ! . Optional array arguments.
   INTEGER, DIMENSION(1:3), INTENT(IN), OPTIONAL :: KMAXIMUM

   ! . Reinitialize the statistics counters.
   NCALLS    = 0
   NUPDATESK = 0
   NUPDATESL = 0
   SUM_INTS  = 0.0_DP
   SUM_VEC   = 0.0_DP

   ! . Initialize the old box size and the number of wave-vectors.
   BOXOLD = 0.0_DP
   NKVEC  = 0

   ! . Deallocate the module arrays.
   IF ( ALLOCATED ( LISTI   ) ) DEALLOCATE ( LISTI   )
   IF ( ALLOCATED ( LISTJ   ) ) DEALLOCATE ( LISTJ   )
   IF ( ALLOCATED ( REFCRD  ) ) DEALLOCATE ( REFCRD  )
   IF ( ALLOCATED ( RKEXPEL ) ) DEALLOCATE ( RKEXPEL )
   IF ( ALLOCATED ( RKEXPLJ ) ) DEALLOCATE ( RKEXPLJ )
   IF ( ALLOCATED ( RKVEC   ) ) DEALLOCATE ( RKVEC   )

   ! . The covalent energy terms.
   IF ( PRESENT ( DIELECTRIC   ) ) EPSILON    = DIELECTRIC
   IF ( PRESENT ( KAPPA        ) ) BETA       = KAPPA
   IF ( PRESENT ( KMAXIMUM     ) ) KMAX       = KMAXIMUM
   IF ( PRESENT ( LIST_CUTOFF  ) ) CUT_LIST   = LIST_CUTOFF
   IF ( PRESENT ( OUTER_CUTOFF ) ) CUT_OFF    = OUTER_CUTOFF
   IF ( PRESENT ( PRINT        ) ) QPRINT     = PRINT
   IF ( PRESENT ( TINFOIL      ) ) QTINFOIL   = TINFOIL

   ! . Check the values of the cutoffs.
   IF ( CUT_LIST < CUT_OFF ) CALL PRINT_ERROR ( "ENERGY_NON_BONDING_OPTIONS", "Invalid CUT_LIST/CUT_OFF values." )

   ! . Write out some information.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#550055", VARIABLEWIDTH = 16 )
   CALL PRINT_SUMMARY_START ( "Ewald Non-Bonding Energy Options" )
   WRITE ( PRINT_LINE, "(F16.2)" ) CUT_LIST ; CALL PRINT_SUMMARY_ELEMENT ( "List Cutoff"          )
   WRITE ( PRINT_LINE, "(F16.2)" ) EPSILON  ; CALL PRINT_SUMMARY_ELEMENT ( "Dielectric Constant"  )
   WRITE ( PRINT_LINE, "(F16.2)" ) CUT_OFF  ; CALL PRINT_SUMMARY_ELEMENT ( "Outer Switch Cutoff"  )
   WRITE ( PRINT_LINE, "(F16.4)" ) BETA     ; CALL PRINT_SUMMARY_ELEMENT ( "Kappa Value"          )
   WRITE ( PRINT_LINE, "(L16)"   ) QPRINT   ; CALL PRINT_SUMMARY_ELEMENT ( "Update Printing"      )
   WRITE ( PRINT_LINE, "(L16)"   ) QTINFOIL ; CALL PRINT_SUMMARY_ELEMENT ( "Tinfoil Boundary"     )
   WRITE ( PRINT_LINE, "(I16)"   ) KMAX(1)  ; CALL PRINT_SUMMARY_ELEMENT ( "Maximum k-vector (x)" )
   WRITE ( PRINT_LINE, "(I16)"   ) KMAX(2)  ; CALL PRINT_SUMMARY_ELEMENT ( "Maximum k-vector (y)" )
   WRITE ( PRINT_LINE, "(I16)"   ) KMAX(3)  ; CALL PRINT_SUMMARY_ELEMENT ( "Maximum k-vector (z)" )
   WRITE ( PRINT_LINE, "(G16.4)" ) ERFC ( BETA * CUT_OFF )
   CALL PRINT_SUMMARY_ELEMENT ( "Erfc(kappa*cutoff)" )
   CALL PRINT_SUMMARY_STOP

   END SUBROUTINE ENERGY_NON_BONDING_OPTIONS

   !---------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_STATISTICS
   !---------------------------------------

   ! . Write out some information.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#550055", VARIABLEWIDTH = 16 )
   CALL PRINT_SUMMARY_START ( "Ewald Non-Bonding Energy Statistics" )
   WRITE ( PRINT_LINE, "(I16)" ) NCALLS    ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Calls"   )
   WRITE ( PRINT_LINE, "(I16)" ) NUPDATESL ; CALL PRINT_SUMMARY_ELEMENT ( "List Updates" )
   WRITE ( PRINT_LINE, "(I16)" ) NUPDATESK ; CALL PRINT_SUMMARY_ELEMENT ( "Wave-Vector Updates" )
   IF ( NUPDATESL > 0 ) THEN
      WRITE ( PRINT_LINE, "(F16.4)" ) REAL ( NCALLS, DP ) / REAL ( NUPDATESL, DP )
      CALL PRINT_SUMMARY_ELEMENT ( "List Update Freq."   )
      WRITE ( PRINT_LINE, "(F16.4)" ) SUM_INTS / REAL ( NUPDATESL, DP )
      CALL PRINT_SUMMARY_ELEMENT ( "<Number of Ints.>" )
   END IF
   IF ( NUPDATESK > 0 ) THEN
      WRITE ( PRINT_LINE, "(F16.4)" ) REAL ( NCALLS, DP ) / REAL ( NUPDATESK, DP )
      CALL PRINT_SUMMARY_ELEMENT ( "Vector Update Freq."   )
      WRITE ( PRINT_LINE, "(F16.4)" ) SUM_INTS / REAL ( NUPDATESK, DP )
      CALL PRINT_SUMMARY_ELEMENT ( "<Number of Vectors>" )
   END IF
   CALL PRINT_SUMMARY_STOP

   ! . Reinitialize the integer counters.
   NCALLS    = 0
   NUPDATESK = 0
   NUPDATESL = 0
   SUM_INTS  = 0.0_DP
   SUM_VEC   = 0.0_DP

   ! . Reinitialize the old box size and the number of k-vectors.
   BOXOLD = 0.0_DP
   NKVEC  = 0

   END SUBROUTINE ENERGY_NON_BONDING_STATISTICS

!========================================================================
! . Unimplemened subroutines.
!========================================================================

   !--------------------------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_INTERACTION ( EELECT, ELJ, SELECTION1, SELECTION2 )
   !--------------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EELECT, ELJ

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION1, SELECTION2

   ! . Initialize the energies.
   EELECT = 0.0_DP
   ELJ    = 0.0_DP

   ! . Check for overlaps between the atom selections.
   IF ( ANY ( SELECTION1 .AND. SELECTION2 ) ) THEN
      CALL PRINT_ERROR ( "ENERGY_NON_BONDING_INTERACTION", "There is overlap between the atom selections." )
   END IF

   ! . Write out a warning.
   CALL PRINT_ERROR ( "ENERGY_NON_BONDING_INTERACTION", "This subroutine is not implemented at present." )

   END SUBROUTINE ENERGY_NON_BONDING_INTERACTION

   !------------------------------------------------------------
   SUBROUTINE ENERGY_NON_BONDING_SELF ( EELECT, ELJ, SELECTION )
   !------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EELECT, ELJ

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION

   ! . Initialize the energies.
   EELECT = 0.0_DP
   ELJ    = 0.0_DP

   ! . Check for a null atom selection.
   IF ( .NOT. ANY ( SELECTION ) ) RETURN

   ! . Write out a warning.
   CALL PRINT_ERROR ( "ENERGY_NON_BONDING_SELF", "This subroutine is not implemented at present." )

   END SUBROUTINE ENERGY_NON_BONDING_SELF

END MODULE ENERGY_NON_BONDING !_EWALD
