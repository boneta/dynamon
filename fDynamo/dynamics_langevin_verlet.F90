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
!                          The Langevin Dynamics Module
!===============================================================================
!
! . Subroutines:
!
!   LANGEVIN_VERLET_DYNAMICS       Do a molecular dynamics simulation.
!
!===============================================================================
MODULE DYNAMICS_LANGEVIN_VERLET

! . Module declarations.
USE CONSTANTS,        ONLY : AMU_TO_KG, KBOLTZ, MS_TO_APS
USE DEFINITIONS,      ONLY : DP
USE PRINTING,         ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, &
                             PRINT_SUMMARY_STOP
USE RANDOM_NUMBERS,   ONLY : RANDOM_GAUSS

USE ATOMS,            ONLY : ATMCRD, ATMFIX, ATMMAS, NATOMS, NFREE
USE DYNAMICS_UTILITIES
USE VELOCITY,         ONLY : ATMVEL, EKE, TEMPERATURE, VELOCITY_TEMPERATURE

#ifdef	HACK_LAMB
use mopac_integrals,  only : mopac_lambda_value, mopac_lambda_delta
use quantum_energy,   only : energy_quantum
#endif

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !---------------------------------------------------
   SUBROUTINE LANGEVIN_VERLET_DYNAMICS ( TBATH, GAMMA )
   !---------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: GAMMA, TBATH

   ! . Local scalars.
   INTEGER            :: I, ISTEP
   REAL ( KIND = DP ) :: C0, C1, C2, CRV1, CRV2, EPOT, FACR1, FACR2, FACV1, FACV2, FACV3, FACT, SDR, SDV, TIME

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NATOMS)     :: MASFAC
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: ATMACC, ATMTMP

#ifdef	HACK_LAMB
   real( kind=dp ) :: f_lambda, prv_ene, cur_ene, nxt_ene, f_epi, f_virial
#endif

   ! . Check the number of free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Check the options.
   IF ( TBATH <= 0.0_DP ) CALL PRINT_ERROR ( "LANGEVIN_VERLET_DYNAMICS", "Invalid TBATH for Langevin dynamics." )
   IF ( GAMMA <= 0.0_DP ) CALL PRINT_ERROR ( "LANGEVIN_VERLET_DYNAMICS", "Invalid GAMMA for Langevin dynamics." )

   ! . Check that a velocity array exists.
   IF ( .NOT. ALLOCATED ( ATMVEL ) ) CALL PRINT_ERROR ( "LANGEVIN_VERLET_DYNAMICS", "No velocities exist." )

   ! . Print out the options.
   IF ( DYNAMICS_PRFRQ > 0 ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#005500", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Langevin Verlet Molecular Dynamics" )
      WRITE ( PRINT_LINE, "(G16.6)" ) TBATH ; CALL PRINT_SUMMARY_ELEMENT ( "Bath Temperature"    )
      WRITE ( PRINT_LINE, "(G16.6)" ) GAMMA ; CALL PRINT_SUMMARY_ELEMENT ( "Collision Frequency" )
      CALL PRINT_SUMMARY_STOP
   END IF

   ! . Calculate some quantities necessary for the integration.
   ! . Integration constants.
   FACT = DYNAMICS_DELTA * GAMMA
   C0   = EXP ( - FACT )
   C1   = ( 1.0_DP - C0 ) / FACT
   C2   = ( 1.0_DP - C1 ) / FACT

   ! . Random force constants.
   SDR  = SQRT ( DYNAMICS_DELTA**2 * ( 2.0_DP - ( 3.0_DP - 4.0_DP * C0 + C0 * C0 ) / FACT ) / FACT ) ! ps.
   SDV  = SQRT ( ( 1.0_DP - C0 * C0 ) )                                                              ! Dimensionless.
   CRV1 = DYNAMICS_DELTA * ( 1.0_DP - C0 )**2 / ( FACT * SDR * SDV )                                 ! Dimensionless.
   CRV2 = SQRT ( 1.0_DP - CRV1 * CRV1 )                                                              ! Dimensionless.

   ! . Temperature and mass factors.
   FACT = ( KBOLTZ * TBATH / AMU_TO_KG )
   WHERE ( .NOT. ATMFIX ) MASFAC = MS_TO_APS * SQRT ( FACT / ATMMAS ) ! A ps^-1.

   ! . Calculate some integration constants.
   FACR1 = C1 * DYNAMICS_DELTA          ! ps.
   FACR2 = C2 * DYNAMICS_DELTA**2       ! ps^2.
   FACV1 = C0                           ! Dimensionless.
   FACV2 = ( C1 - C2 ) * DYNAMICS_DELTA ! ps.
   FACV3 = C2 * DYNAMICS_DELTA          ! ps.

   ! . Activate the trajectories if necessary.
   CALL DYNAMICS_SAVE_ACTIVATE

   ! . Initialize the time.
   TIME = 0.0_DP

   ! . Calculate the initial energy and accelerations.
   CALL ENERGY_AND_ACCELERATIONS ( EPOT, ATMACC )

   ! . Calculate the temperature.
   CALL VELOCITY_TEMPERATURE

   ! . Initialize the statistics calculation.
   CALL DYNAMICS_STATISTICS_INITIALIZE ( TIME, EKE, EPOT, TEMPERATURE, 0.0_DP, 0.0_DP )

#ifdef	HACK_LAMB
if( mopac_lambda_delta /= 0.0_dp ) then

	write( 999, "(a,f19.10,2f20.10)" ) "#", mopac_lambda_value - mopac_lambda_delta, &
			mopac_lambda_value, mopac_lambda_value + mopac_lambda_delta
   do istep = 1, dynamics_steps
      time = time + dynamics_delta
      ! -----------------------------------------------------------------
      f_lambda = mopac_lambda_value
      call energy_quantum( cur_ene, f_epi, f_virial, print = .false. )
      if( f_lambda - mopac_lambda_delta < .0_dp ) then
         prv_ene = cur_ene
      else
         mopac_lambda_value = f_lambda - mopac_lambda_delta
         call energy_quantum( prv_ene, f_epi, f_virial, print = .false. )
      end if
      if( f_lambda + mopac_lambda_delta > 1._dp ) then
         nxt_ene = cur_ene
      else
         mopac_lambda_value = f_lambda + mopac_lambda_delta
         call energy_quantum( nxt_ene, f_epi, f_virial, print = .false. )
      end if
      mopac_lambda_value = f_lambda
      write( 999, "(3f20.10)" ) prv_ene, cur_ene, nxt_ene
      ! -----------------------------------------------------------------
      do i = 1, natoms
         if( .not. atmfix(i) ) atmcrd(1:3,i) = atmcrd(1:3,i) + facr1 * atmvel(1:3,i) + facr2 * atmacc(1:3,i)
      end do
      do i = 1, natoms
         if( .not. atmfix(i) ) atmtmp(1:3,i) = facv1 * atmvel(1:3,i) + facv2 * atmacc(1:3,i)
      end do
      call random_forces
      call energy_and_accelerations( epot, atmacc )
      do i = 1, natoms
         if( .not. atmfix(i) ) atmvel(1:3,i) = atmtmp(1:3,i) + facv3 * atmacc(1:3,i)
      end do
      call velocity_temperature
      call dynamics_statistics_accumulate( istep, time, eke, epot, temperature, 0.0_dp, 0.0_dp )
      call dynamics_save_write( istep )
   end do

else
#endif

   ! . Start the loop over steps.
   DO ISTEP = 1,DYNAMICS_STEPS

      ! . Increment the time.
      TIME = TIME + DYNAMICS_DELTA

      ! . Calculate the new position.
      DO I = 1,NATOMS
         IF ( .NOT. ATMFIX(I) ) ATMCRD(1:3,I) = ATMCRD(1:3,I) + FACR1 * ATMVEL(1:3,I) + FACR2 * ATMACC(1:3,I)
      END DO

      ! . Calculate the half-velocities.
      DO I = 1,NATOMS
         IF ( .NOT. ATMFIX(I) ) ATMTMP(1:3,I) = FACV1 * ATMVEL(1:3,I) + FACV2 * ATMACC(1:3,I)
      END DO

      ! . Add in the random forces.
      CALL RANDOM_FORCES

      ! . Calculate the new energy and accelerations.
      CALL ENERGY_AND_ACCELERATIONS ( EPOT, ATMACC )

      ! . Calculate the new velocities.
      DO I = 1,NATOMS
         IF ( .NOT. ATMFIX(I) ) ATMVEL(1:3,I) = ATMTMP(1:3,I) + FACV3 * ATMACC(1:3,I)
      END DO

      ! . Calculate the temperature.
      CALL VELOCITY_TEMPERATURE

      ! . Accumulate the statistics.
      CALL DYNAMICS_STATISTICS_ACCUMULATE ( ISTEP, TIME, EKE, EPOT, TEMPERATURE, 0.0_DP, 0.0_DP )

      ! . Save some data if necessary.
      CALL DYNAMICS_SAVE_WRITE ( ISTEP )

   END DO

#ifdef	HACK_LAMB
end if
#endif

   ! . Print the statistics.
   CALL DYNAMICS_STATISTICS_PRINT

   ! . Close the trajectories.
   CALL DYNAMICS_SAVE_DEACTIVATE

   ! . Reinitialize the dynamics data structure.
   CALL DYNAMICS_INITIALIZE

   !===========================================================================
   CONTAINS
   !===========================================================================

      !-----------------------
      SUBROUTINE RANDOM_FORCES
      !-----------------------

      ! . Local scalars.
      INTEGER            :: I, IATOM
      REAL ( KIND = DP ) :: RANR, RANV, R1, R2

      ! . Loop over the atoms.
      DO IATOM = 1,NATOMS

         ! . The atom is free.
         IF ( .NOT. ATMFIX(IATOM) ) THEN

            ! . Loop over the Cartesian components for the atom.
            DO I = 1,3

               ! . Generate two Gaussian random numbers.
               R1 = RANDOM_GAUSS ( 0.0_DP, 1.0_DP )
               R2 = RANDOM_GAUSS ( 0.0_DP, 1.0_DP )

               ! . Calculate the position and velocity changes.
               RANR = MASFAC(IATOM) * SDR * R1
               RANV = MASFAC(IATOM) * SDV * ( CRV1 * R1 + CRV2 * R2 )

               ! . Add in the components.
               ATMCRD(I,IATOM) = ATMCRD(I,IATOM) + RANR
               ATMTMP(I,IATOM) = ATMTMP(I,IATOM) + RANV

            END DO
         END IF
      END DO

      END SUBROUTINE RANDOM_FORCES

   END SUBROUTINE LANGEVIN_VERLET_DYNAMICS

END MODULE DYNAMICS_LANGEVIN_VERLET
