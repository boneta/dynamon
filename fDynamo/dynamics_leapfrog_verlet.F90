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
!                      The Leapfrog Verlet Dynamics Module
!===============================================================================
!
! . Subroutines:
!
!   LEAPFROG_VERLET_DYNAMICS       Do a molecular dynamics simulation.
!
!===============================================================================
MODULE DYNAMICS_LEAPFROG_VERLET

! . Module declarations.
USE CONSTANTS,        ONLY : PV_TO_KJ_MOLE
USE DEFINITIONS,      ONLY : DP
USE PRINTING,         ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, &
                             PRINT_SUMMARY_STOP

USE ATOMS,            ONLY : ATMCRD, ATMFIX, NATOMS, NATOMSQM, NFIXED, NFREE
USE CONSTRAINT,       ONLY : QTETHER
USE DYNAMICS_UTILITIES
USE POTENTIAL_ENERGY, ONLY : VIRIAL
USE SYMMETRY,         ONLY : BOXL, QBOX
USE VELOCITY,         ONLY : ATMVEL, EKE, TEMPERATURE, VELOCITY_TEMPERATURE

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !-------------------------------------------------------------------------------------------------------------------
   SUBROUTINE LEAPFROG_VERLET_DYNAMICS ( TARGET_TEMPERATURE, TEMPERATURE_COUPLING, TARGET_PRESSURE, PRESSURE_COUPLING )
   !-------------------------------------------------------------------------------------------------------------------

   ! . Optional scalar arguments.
   REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: PRESSURE_COUPLING, TARGET_PRESSURE, TARGET_TEMPERATURE, &
                                                  TEMPERATURE_COUPLING

   ! . Local scalars.
   INTEGER                :: I, ISTEP
   LOGICAL                :: QPCALC, QPCOUPLE, QTCOUPLE
   REAL ( KIND = DP )     :: EPOT, PCOUPLE, PRES, PTARGET, TCOUPLE, TIME, TTARGET, VOLU

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: ATMACC

   ! . Check the number of free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Initialize the input options.
   QPCOUPLE = .FALSE.
   QTCOUPLE = .FALSE.
   PCOUPLE  =          0.0_DP
   PTARGET  = - HUGE ( 0.0_DP )
   TCOUPLE  =          0.0_DP
   TTARGET  = - HUGE ( 0.0_DP )

   ! . Assign the input options.
   IF ( PRESENT ( PRESSURE_COUPLING    ) ) PCOUPLE = PRESSURE_COUPLING
   IF ( PRESENT ( TARGET_PRESSURE      ) ) PTARGET = TARGET_PRESSURE
   IF ( PRESENT ( TARGET_TEMPERATURE   ) ) TTARGET = TARGET_TEMPERATURE
   IF ( PRESENT ( TEMPERATURE_COUPLING ) ) TCOUPLE = TEMPERATURE_COUPLING

   ! . Check that a velocity array exists.
   IF ( .NOT. ALLOCATED ( ATMVEL ) ) CALL PRINT_ERROR ( "LEAPFROG_VERLET_DYNAMICS", "No velocities exist." )

   ! . Check the pressure options.
   IF ( PCOUPLE /= 0.0_DP ) THEN
      IF ( .NOT. QBOX ) THEN
         CALL PRINT_ERROR ( "LEAPFROG_VERLET_DYNAMICS", "Pressure control must be done with periodic boundary conditions." )
      END IF
      IF ( ( PCOUPLE <= 0.0_DP ) .OR. ( PTARGET <= 0.0_DP ) ) THEN
         CALL PRINT_ERROR ( "LEAPFROG_VERLET_DYNAMICS", "Invalid pressure coupling or target variables." )
      END IF
      IF ( ( NATOMSQM > 0 ) .OR. ( NFIXED > 0 ) .OR. QTETHER ) THEN
         CALL PRINT_ERROR ( "LEAPFROG_VERLET_DYNAMICS", "Pressure control not allowed with fixed, quantum or tethered atoms." )
      END IF
      QPCOUPLE = .TRUE.
   END IF

   ! . Check the temperature options.
   IF ( TCOUPLE /= 0.0_DP ) THEN
      IF ( ( TCOUPLE <= 0.0_DP ) .OR. ( TTARGET <= 0.0_DP ) ) THEN
         CALL PRINT_ERROR ( "LEAPFROG_VERLET_DYNAMICS", "Invalid temperature coupling or target variables." )
      END IF
      QTCOUPLE = .TRUE.
   END IF

   ! . Determine the flag for the pressure calculation.
   QPCALC = QPCOUPLE .OR. QBOX

   ! . Print out the options.
   IF ( DYNAMICS_PRFRQ > 0 ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#005500", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Leapfrog Verlet Molecular Dynamics" )
      WRITE ( PRINT_LINE, "(L16)" ) QPCOUPLE ; CALL PRINT_SUMMARY_ELEMENT ( "Pressure Control"    )
      WRITE ( PRINT_LINE, "(L16)" ) QTCOUPLE ; CALL PRINT_SUMMARY_ELEMENT ( "Temperature Control" )
      IF ( QPCOUPLE ) THEN
         WRITE ( PRINT_LINE, "(G16.6)" ) PTARGET ; CALL PRINT_SUMMARY_ELEMENT ( "Target Pressure"   )
         WRITE ( PRINT_LINE, "(G16.6)" ) PCOUPLE ; CALL PRINT_SUMMARY_ELEMENT ( "Pressure Coupling" )
      END IF
      IF ( QTCOUPLE ) THEN
         WRITE ( PRINT_LINE, "(G16.6)" ) TTARGET ; CALL PRINT_SUMMARY_ELEMENT ( "Target Temperature"   )
         WRITE ( PRINT_LINE, "(G16.6)" ) TCOUPLE ; CALL PRINT_SUMMARY_ELEMENT ( "Temperature Coupling" )
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

   ! . Calculate the pressure.
   CALL PRESSURE_CALCULATE

   ! . Initialize the statistics calculation.
   CALL DYNAMICS_STATISTICS_INITIALIZE ( TIME, EKE, EPOT, TEMPERATURE, PRES, VOLU )

   ! . Start the loop over steps.
   DO ISTEP = 1,DYNAMICS_STEPS

      ! . Increment the time.
      TIME = TIME + DYNAMICS_DELTA

      ! . Calculate the new velocities.
      DO I = 1,NATOMS
         IF ( .NOT. ATMFIX(I) ) ATMVEL(1:3,I) = ATMVEL(1:3,I) + DYNAMICS_DELTA * ATMACC(1:3,I)
      END DO

      ! . Scale the velocities if necessary.
      CALL TEMPERATURE_SCALING

      ! . Calculate the new position.
      DO I = 1,NATOMS
         IF ( .NOT. ATMFIX(I) ) ATMCRD(1:3,I) = ATMCRD(1:3,I) + DYNAMICS_DELTA * ATMVEL(1:3,I)
      END DO

      ! . Scale the coordinates and box size.
      CALL PRESSURE_SCALING

      ! . Calculate the initial energy and accelerations.
      CALL ENERGY_AND_ACCELERATIONS ( EPOT, ATMACC )

      ! . Calculate the temperature.
      CALL VELOCITY_TEMPERATURE

      ! . Calculate the pressure.
      CALL PRESSURE_CALCULATE

      ! . Accumulate the statistics.
      CALL DYNAMICS_STATISTICS_ACCUMULATE ( ISTEP, TIME, EKE, EPOT, TEMPERATURE, PRES, VOLU )

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

      !----------------------------
      SUBROUTINE PRESSURE_CALCULATE
      !----------------------------

      ! . The pressure is to be calculated.
      IF ( QPCALC ) THEN

         ! . Calculate the volume.
         VOLU = PRODUCT ( BOXL )

         ! . Calculate the pressure (in atmospheres).
         PRES = ( 2.0_DP * EKE - VIRIAL ) / ( 3.0_DP * VOLU * PV_TO_KJ_MOLE )

      ! . No pressure is to be calculated.
      ELSE
         PRES = 0.0_DP
         VOLU = 0.0_DP
      END IF

      END SUBROUTINE PRESSURE_CALCULATE

      !--------------------------
      SUBROUTINE PRESSURE_SCALING
      !--------------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: ZETAP

      ! . Pressure scaling is to be done.
      IF ( QPCOUPLE ) THEN

         ! . Calculate the scale factor.
         ZETAP = EXP ( LOG ( 1.0_DP - ( DYNAMICS_DELTA / PCOUPLE ) * ( PTARGET - PRES ) ) / 3.0_DP )

         ! . Scale the box size.
         BOXL = ZETAP * BOXL

         ! . Scale the coordinates.
         ATMCRD = ZETAP * ATMCRD

      END IF

      END SUBROUTINE PRESSURE_SCALING

      !-----------------------------
      SUBROUTINE TEMPERATURE_SCALING
      !-----------------------------

      ! . Local scalars.
      INTEGER            :: I
      REAL ( KIND = DP ) :: ZETAT

      ! . Pressure scaling is to be done.
      IF ( QTCOUPLE ) THEN

         ! . Calculate the scale factor.
         ZETAT = SQRT ( 1.0_DP + ( DYNAMICS_DELTA / TCOUPLE ) * ( ( TTARGET / TEMPERATURE ) - 1.0_DP ) )

         ! . Scale the velocities.
	 DO I = 1,NATOMS
	    IF ( .NOT. ATMFIX(I) ) ATMVEL(1:3,I) = ZETAT * ATMVEL(1:3,I)
	 END DO

         ! . Recalculate the temperature.
	 CALL VELOCITY_TEMPERATURE

      END IF

      END SUBROUTINE TEMPERATURE_SCALING

   END SUBROUTINE LEAPFROG_VERLET_DYNAMICS

END MODULE DYNAMICS_LEAPFROG_VERLET

