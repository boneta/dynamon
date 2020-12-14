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
!                      The Velocity Verlet Dynamics Module
!===============================================================================
!
! . Subroutines:
!
!   VELOCITY_VERLET_DYNAMICS       Do a molecular dynamics simulation.
!
!===============================================================================
MODULE DYNAMICS_VELOCITY_VERLET

! . Module declarations.
USE DEFINITIONS,      ONLY : DP
USE PRINTING,         ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, &
                             PRINT_SUMMARY_STOP

USE ATOMS,            ONLY : ATMCRD, ATMFIX, NATOMS, NFREE
USE DYNAMICS_UTILITIES
USE VELOCITY,         ONLY : ATMVEL, EKE, TEMPERATURE, VELOCITY_SCALE, VELOCITY_TEMPERATURE

#ifdef	HACK_GH
USE GH
#endif


IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !----------------------------------------------------------------------------------------
   SUBROUTINE VELOCITY_VERLET_DYNAMICS ( TARGET_TEMPERATURE, SCALE_FREQUENCY, SCALE_OPTION )
   !----------------------------------------------------------------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: SCALE_OPTION
   INTEGER,               INTENT(IN), OPTIONAL :: SCALE_FREQUENCY
   REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: TARGET_TEMPERATURE

   ! . Local scalars.
   CHARACTER ( LEN = 16 ) :: TOPTION
   INTEGER                :: I, ISTEP, MODFRQ
   REAL ( KIND = DP )     :: EPOT, FACR1, FACR2, FACV, TIME, TSCALE, TSTART, TSTOP

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: ATMACC, ATMTMP

   ! . Check the number of free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Initialize the input options.
   TOPTION = " "
   MODFRQ  = 0
   TSTART  = - HUGE ( 0.0_DP )
   TSTOP   = - HUGE ( 0.0_DP )

   ! . Assign the input options.
   IF ( PRESENT ( TARGET_TEMPERATURE  ) ) TSTOP   = TARGET_TEMPERATURE
   IF ( PRESENT ( SCALE_FREQUENCY     ) ) MODFRQ  = SCALE_FREQUENCY
   IF ( PRESENT ( SCALE_OPTION        ) ) TOPTION = SCALE_OPTION

   ! . Check the input options.
   IF ( ( MODFRQ < 0 ) .OR. ( MODFRQ > DYNAMICS_STEPS ) ) MODFRQ = 0

   ! . Check the velocity options.
   CALL CHECK_VELOCITY_OPTIONS

   ! . Print out the options.
   IF ( DYNAMICS_PRFRQ > 0 ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#005500", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Velocity Verlet Molecular Dynamics" )
      WRITE ( PRINT_LINE, "(G16.6)" ) TSTART ; CALL PRINT_SUMMARY_ELEMENT ( "Initial Temperature" )
      WRITE ( PRINT_LINE, "(I16)"   ) MODFRQ ; CALL PRINT_SUMMARY_ELEMENT ( "Scale Frequency"     )
      IF ( MODFRQ > 0 ) THEN
         CALL PRINT_SUMMARY_ELEMENT ( "Scale Option", TEXT = TOPTION )
         WRITE ( PRINT_LINE, "(G16.6)" ) TSTOP ; CALL PRINT_SUMMARY_ELEMENT ( "Target Temperature" )
      END IF
      CALL PRINT_SUMMARY_STOP
   END IF

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

   ! . Calculate some integration constants.
   FACR1 =          DYNAMICS_DELTA                  ! ps.
   FACR2 = 0.5_DP * DYNAMICS_DELTA * DYNAMICS_DELTA ! ps^2.
   FACV  = 0.5_DP * DYNAMICS_DELTA                  ! ps.

#ifdef	HACK_GH
   CALL GH_VV_INIT( FACV, FACR2 )
#endif

   ! . Start the loop over steps.
   DO ISTEP = 1,DYNAMICS_STEPS

      ! . Increment the time.
      TIME = TIME + DYNAMICS_DELTA

#ifdef	HACK_GH
   CALL GH_VV_A
#endif

      ! . Calculate the new position.
      DO I = 1,NATOMS
         IF ( .NOT. ATMFIX(I) ) ATMCRD(1:3,I) = ATMCRD(1:3,I) + FACR1 * ATMVEL(1:3,I) + FACR2 * ATMACC(1:3,I)
      END DO

#ifdef	HACK_GH
   CALL GH_VV_B( ISTEP )
#endif

      ! . Calculate the half-velocities.
      DO I = 1,NATOMS
         IF ( .NOT. ATMFIX(I) ) ATMTMP(1:3,I) = ATMVEL(1:3,I) + FACV * ATMACC(1:3,I)
      END DO

#ifdef	HACK_GH
   CALL GH_VV_C( ATMTMP )
#endif

      ! . Calculate the new energy and accelerations.
      CALL ENERGY_AND_ACCELERATIONS ( EPOT, ATMACC )

      ! . Calculate the new velocities.
      DO I = 1,NATOMS
         IF ( .NOT. ATMFIX(I) ) ATMVEL(1:3,I) = ATMTMP(1:3,I) + FACV * ATMACC(1:3,I)
      END DO

      ! . Calculate the temperature.
      CALL VELOCITY_TEMPERATURE

      ! . Check for velocity modification.
      IF ( MODFRQ > 0 ) THEN
         IF ( MOD ( ISTEP, MODFRQ ) == 0 ) CALL SCALE_VELOCITIES
      END IF

      ! . Accumulate the statistics.
      CALL DYNAMICS_STATISTICS_ACCUMULATE ( ISTEP, TIME, EKE, EPOT, TEMPERATURE, 0.0_DP, 0.0_DP )

      ! . Save some data if necessary.
      CALL DYNAMICS_SAVE_WRITE ( ISTEP )

   END DO

   ! . Print the statistics.
   CALL DYNAMICS_STATISTICS_PRINT

   ! . Close the trajectories.
   CALL DYNAMICS_SAVE_DEACTIVATE

   ! . Reinitialize the dynamics data structure.
   CALL DYNAMICS_INITIALIZE

   !===========================================================================
   CONTAINS
   !===========================================================================

      !--------------------------------
      SUBROUTINE CHECK_VELOCITY_OPTIONS
      !--------------------------------

      ! . Check that a velocity array exists.
      IF ( .NOT. ALLOCATED ( ATMVEL ) ) CALL PRINT_ERROR ( "VELOCITY_VERLET_DYNAMICS", "No velocities exist." )

      ! . Set the starting temperature.
      TSTART = TEMPERATURE

      ! . Initialize all options if the scaling frequency is zero.
      IF ( MODFRQ == 0 ) THEN
         TOPTION = " "
         TSTOP   = 0.0_DP

      ! . Check the scaling options.
      ELSE

         ! . Check TSTOP.
         IF ( TSTOP <= 0.0_DP ) THEN
            CALL PRINT_ERROR ( "VELOCITY_VERLET_DYNAMICS", "Target temperature not specified for scaling options." )
         END IF

         ! . Branch on the modify option and calculate the scaling factor.
         SELECT CASE ( TOPTION )
         CASE ( "CONSTANT        " ) ; TSCALE = 0.0_DP ; CALL VELOCITY_SCALE ( TSTOP ) ; TSTART = TSTOP
         CASE ( "EXPONENTIAL     " ) ; TSCALE = LOG ( TSTOP / TSTART ) / DYNAMICS_TIME
         CASE ( "LINEAR          " ) ; TSCALE = ( TSTOP - TSTART ) / DYNAMICS_TIME
         CASE DEFAULT ; CALL PRINT_ERROR ( "VELOCITY_VERLET_DYNAMICS", "Invalid SCALE OPTION." )
         END SELECT

      END IF

      END SUBROUTINE CHECK_VELOCITY_OPTIONS

      !--------------------------
      SUBROUTINE SCALE_VELOCITIES
      !--------------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: TNEEDED

      ! . Branch on the modify option.
      SELECT CASE ( TOPTION )
      CASE ( "CONSTANT   " ) ; TNEEDED = TSTART
      CASE ( "EXPONENTIAL" ) ; TNEEDED = TSTART * EXP ( TSCALE * TIME )
      CASE ( "LINEAR     " ) ; TNEEDED = TSCALE * TIME + TSTART 
      CASE DEFAULT ; RETURN
      END SELECT

      ! . Scale the velocities.
      CALL VELOCITY_SCALE ( TNEEDED )

      END SUBROUTINE SCALE_VELOCITIES

   END SUBROUTINE VELOCITY_VERLET_DYNAMICS

END MODULE DYNAMICS_VELOCITY_VERLET

