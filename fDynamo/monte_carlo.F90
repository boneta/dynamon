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
!                       The Monte Carlo Simulation Module
!===============================================================================
!
! . Subroutines:
!
!   MONTE_CARLO                     Do a Monte Carlo simulation.
!
! . Notes:
!
!   Neither fixed atoms or quantum atoms can be handled by this module.
!
!===============================================================================
MODULE MONTE_CARLO_SIMULATION

! . Module declarations.
USE CONSTANTS,      ONLY : PV_TO_KJ_MOLE, R, TO_DEGREES, TO_RADIANS
USE DEFINITIONS,    ONLY : DP
USE PRINTING,       ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
                           PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT,          &
			   PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP
USE RANDOM_NUMBERS, ONLY : RANDOM

USE ATOMS,          ONLY : ATMCRD, ATMMAS, NATOMS, NATOMSQM, NFIXED
USE DCD_IO
USE SEQUENCE,       ONLY : MOLIND => RESIND,  NMOLECULES => NRESID
USE SYMMETRY,       ONLY : BOXL, BOX_TYPE, QBOX
USE TRANSFORMATION, ONLY : CENTER

USE MONTE_CARLO_ENERGY, ONLY : MONTE_CARLO_ENERGY_FULL, MONTE_CARLO_ENERGY_ONE

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !--------------------------------------------------------------------------------------------
   SUBROUTINE MONTE_CARLO ( NBLOCK, NCONF, ADJUST_FREQUENCY, SAVE_FREQUENCY, VOLUME_FREQUENCY, &
                                           ACCEPTANCE, ROTATION, TRANSLATION, VOLUME_MOVE,     &
                                           PRESSURE, TEMPERATURE, FILE )
   !--------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: NBLOCK, NCONF

   ! . Optional arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE
   INTEGER,               INTENT(IN), OPTIONAL :: ADJUST_FREQUENCY, SAVE_FREQUENCY, VOLUME_FREQUENCY
   REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: ACCEPTANCE, ROTATION, TRANSLATION, VOLUME_MOVE, &
                                                  PRESSURE, TEMPERATURE

   ! . Option scalars.
   INTEGER            :: ADJFRQ, SAVFRQ, VOLFRQ
   REAL ( KIND = DP ) :: ACCRAT, RMAX, TMAX, VMAX, P, T

   ! . Miscellaneous scalars.
   INTEGER            :: IBLOCK, ICONF, IMOVE, NREJECT, NREJECTM, NREJECTV, NTRYM, NTRYV
   LOGICAL            :: QVMOVE
   REAL ( KIND = DP ) :: BETA, ECURRENT, PFACT, TFACT, VOLUME

   ! . Statistics scalars.
   INTEGER            :: NTOT
   REAL ( KIND = DP ) :: EAV, EAV2, ETOT, ETOT2, ETOTB2, &
                         HAV, HAV2, HTOT, HTOT2, HTOTB2, &
                         VAV, VAV2, VTOT, VTOT2, VTOTB2

   ! . The trajectory.
   TYPE(DCD_TYPE) :: MCTRAJECTORY

   !---------------------------------------------------------------------------
   ! . Check the input options.
   !---------------------------------------------------------------------------
   ! . Check that there are atoms and molecules.
   IF ( ( NATOMS <= 0 ) .OR. ( NMOLECULES <= 0 ) ) RETURN

   ! . Check the essential input arguments.
   IF ( ( NBLOCK * NCONF ) <= 0 ) RETURN

   ! . Check for fixed or quantum atoms.
   IF ( ( NFIXED > 0 ) .OR. ( NATOMSQM > 0 ) ) THEN
      CALL PRINT_ERROR ( "MONTE_CARLO", "The Monte Carlo module cannot handle fixed or quantum atoms." )
   END IF

   ! . Check that cubic periodic boundary conditions have been defined.
   IF ( .NOT. QBOX .OR. ( BOX_TYPE /= "CUBIC" ) ) THEN
      CALL PRINT_ERROR ( "MONTE_CARLO", "A cubic periodic box has not been defined." )
   END IF

   ! . Initialize the optional arguments (frequencies).
   ADJFRQ = 1000
   SAVFRQ = 0
   VOLFRQ = 500

   ! . Initialize the optional arguments (move parameters).
   ACCRAT = 0.4_DP
   RMAX   = 15.0_DP * TO_RADIANS
   TMAX   = 0.15_DP
   VMAX   = 400.0_DP

   ! . Initialize the optional arguments (pressure and temperature).
   P = 1.0_DP
   T = 300.0_DP

   ! . Check the optional arguments.
   IF ( PRESENT ( ADJUST_FREQUENCY ) ) ADJFRQ = ADJUST_FREQUENCY
   IF ( PRESENT ( SAVE_FREQUENCY   ) ) SAVFRQ = SAVE_FREQUENCY
   IF ( PRESENT ( VOLUME_FREQUENCY ) ) VOLFRQ = VOLUME_FREQUENCY

   IF ( PRESENT ( ACCEPTANCE  ) ) ACCRAT = ACCEPTANCE
   IF ( PRESENT ( ROTATION    ) ) RMAX   = ROTATION * TO_RADIANS
   IF ( PRESENT ( TRANSLATION ) ) TMAX   = TRANSLATION
   IF ( PRESENT ( VOLUME_MOVE ) ) VMAX   = VOLUME_MOVE

   IF ( PRESENT ( PRESSURE    ) ) P = PRESSURE
   IF ( PRESENT ( TEMPERATURE ) ) T = TEMPERATURE

   ! . Write out the Monte Carlo options.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF5500", VARIABLEWIDTH = 16 )
   CALL PRINT_SUMMARY_START ( "Monte Carlo Simulation" )
   WRITE ( PRINT_LINE, "(I16)"   ) NBLOCK            ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Blocks"     )
   WRITE ( PRINT_LINE, "(I16)"   ) NCONF             ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Configs."   )
   WRITE ( PRINT_LINE, "(G16.6)" ) T                 ; CALL PRINT_SUMMARY_ELEMENT ( "Temperature (K)"      )
   WRITE ( PRINT_LINE, "(G16.6)" ) P                 ; CALL PRINT_SUMMARY_ELEMENT ( "Pressure (atm.)"      )
   WRITE ( PRINT_LINE, "(I16)"   ) SAVFRQ            ; CALL PRINT_SUMMARY_ELEMENT ( "Save Frequency"       )
   WRITE ( PRINT_LINE, "(I16)"   ) VOLFRQ            ; CALL PRINT_SUMMARY_ELEMENT ( "Volume Frequency"     )
   WRITE ( PRINT_LINE, "(I16)"   ) ADJFRQ            ; CALL PRINT_SUMMARY_ELEMENT ( "Adjust Frequency"     )
   WRITE ( PRINT_LINE, "(G16.6)" ) ACCRAT            ; CALL PRINT_SUMMARY_ELEMENT ( "Acceptance Ratio"     )
   WRITE ( PRINT_LINE, "(G16.6)" ) RMAX * TO_DEGREES ; CALL PRINT_SUMMARY_ELEMENT ( "Max. Rotation (Deg.)" )
   WRITE ( PRINT_LINE, "(G16.6)" ) TMAX              ; CALL PRINT_SUMMARY_ELEMENT ( "Max. Translation (A)" )
   WRITE ( PRINT_LINE, "(G16.6)" ) VMAX              ; CALL PRINT_SUMMARY_ELEMENT ( "Max. Vol. Move (A^3)" )
   CALL PRINT_SUMMARY_STOP

   !----------------------------------------------------------------------------
   ! . Do the simulation.
   !----------------------------------------------------------------------------
   ! . Initialize some constants.
   BETA  = 1.0_DP / ( R * T )
   PFACT = PV_TO_KJ_MOLE * P
   TFACT = REAL ( NMOLECULES, DP ) / BETA

   ! . Calculate the volume of the box.
   VOLUME = PRODUCT ( BOXL )

   ! . Calculate the energy of the initial configuration.
   ECURRENT = MONTE_CARLO_ENERGY_FULL ( )

   ! . Initialize the trajectory.
   IF ( SAVFRQ > 0 ) THEN

      ! . Check the FILE argument.
      IF ( .NOT. PRESENT ( FILE ) ) CALL PRINT_ERROR ( "MONTE_CARLO", "The FILE argument is missing." )

      ! . Open the trajectory.
      CALL DCD_INITIALIZE ( MCTRAJECTORY )
      CALL DCD_ACTIVATE_WRITE ( FILE, MCTRAJECTORY, "CORD", NATOMS, 0, (NBLOCK*NCONF)/SAVFRQ+1, QCRYSTAL = .TRUE. )
      CALL DCD_WRITE ( MCTRAJECTORY, ATMCRD, BOXL )

   END IF

   ! . Initialize the statistics accumulators.
   CALL STATISTICS_START

   ! . Initialize the move counter.
   IMOVE = 0

   ! . Initialize the move adjustment counters.
   NREJECTM = 0 ; NTRYM = 0
   NREJECTV = 0 ; NTRYV = 0

   ! . Loop over the number of blocks.
   DO IBLOCK = 1,NBLOCK

      ! . Do some initialisation for the block statistics.
      CALL STATISTICS_BLOCK_START

      ! . Loop over the configurations.
      DO ICONF = 1,NCONF

         ! . Increment the IMOVE counter.
         IMOVE = IMOVE + 1

         ! . Determine whether this is a volume move.
         IF ( VOLFRQ > 0 ) THEN
            QVMOVE = MOD ( IMOVE, VOLFRQ ) == 0
         ELSE
            QVMOVE = .FALSE.
         END IF

         ! . Perform a volume move.
         IF ( QVMOVE ) THEN
            CALL MOVE_VOLUME
         ! . Perform a molecule move.
         ELSE
            CALL MOVE_MOLECULE
         END IF

         ! . Accumulate statistics for the configuration.
         CALL STATISTICS_BLOCK_ACCUMULATE

         ! . Check to see if the move sizes need to be adjusted.
         IF ( ADJFRQ > 0 ) THEN
            IF ( MOD ( IMOVE, ADJFRQ ) == 0 ) CALL ADJUST_MOVE_SIZES
         END IF

         ! . Save the configuration.
         IF ( SAVFRQ > 0 ) THEN
            IF ( MOD ( IMOVE, SAVFRQ ) == 0 ) CALL DCD_WRITE ( MCTRAJECTORY, ATMCRD, BOXL )
         END IF

      END DO

      ! . Do the statistics for the block.
      CALL STATISTICS_BLOCK_STOP

   END DO

   !----------------------------------------------------------------------------
   ! . Finish up.
   !----------------------------------------------------------------------------
   ! . Do the statistics for the run.
   CALL STATISTICS_STOP

   ! . Write out the final move sizes.
   IF ( ADJFRQ > 0 ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF5500", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Adjusted Move Sizes" )
      WRITE ( PRINT_LINE, "(G16.6)" ) RMAX * TO_DEGREES ; CALL PRINT_SUMMARY_ELEMENT ( "Max. Rotation (Deg.)" )
      WRITE ( PRINT_LINE, "(G16.6)" ) TMAX              ; CALL PRINT_SUMMARY_ELEMENT ( "Max. Translation (A)" )
      WRITE ( PRINT_LINE, "(G16.6)" ) VMAX              ; CALL PRINT_SUMMARY_ELEMENT ( "Max. Vol. Move (A^3)" )
      CALL PRINT_SUMMARY_STOP
   END IF

   ! . Deactivate the trajectory.
   IF ( SAVFRQ > 0 ) CALL DCD_DEACTIVATE ( MCTRAJECTORY )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !---------------------------
      SUBROUTINE ADJUST_MOVE_SIZES
      !---------------------------

      ! . Local parameters.
      REAL ( KIND = DP ), PARAMETER :: DOWN = 0.95_DP, UP = 1.05_DP

      ! . Adjust the rotation and translation move sizes.
      IF ( NTRYM > 0 ) THEN
         IF ( ( REAL ( NTRYM - NREJECTM, DP ) / REAL ( NTRYM, DP ) ) > ACCRAT ) THEN
            RMAX =   UP * RMAX
            TMAX =   UP * TMAX
         ELSE
            RMAX = DOWN * RMAX
            TMAX = DOWN * TMAX
         END IF

         ! . Reset the accumulators.
         NREJECTM = 0 ; NTRYM = 0

      END IF

      ! . Adjust the box move size.
      IF ( VOLFRQ > 0 ) THEN
         IF ( NTRYV > 0 ) THEN
            IF ( ( REAL ( NTRYV - NREJECTV, DP ) / REAL ( NTRYV, DP ) ) > ACCRAT ) THEN
               VMAX =   UP * VMAX
            ELSE
               VMAX = DOWN * VMAX
            END IF
         END IF

         ! . Reset the accumulators.
         NREJECTV = 0 ; NTRYV = 0

      END IF

      END SUBROUTINE ADJUST_MOVE_SIZES

      !-----------------------
      SUBROUTINE MOVE_MOLECULE
      !-----------------------

      ! . Local scalars.
      INTEGER            :: AXIS, CHOSEN, I, START, STOP
      REAL ( KIND = DP ) :: ANGLE, COSA, DX, DY, DZ, EAFTER, EBEFORE, OLDENE, SINA

      ! . Local arrays.
      REAL ( KIND = DP ),              DIMENSION(1:3) :: DR, MOLCEN, T
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: OLDCRD

      ! . Increment the number of tries.
      NTRYM = NTRYM + 1

      ! . Choose a molecule to move.
      CHOSEN = INT ( REAL ( NMOLECULES, DP ) * RANDOM ( ) ) + 1

      ! . Get the starting and stopping atoms of the molecule.
      START = MOLIND(CHOSEN)+1
      STOP  = MOLIND(CHOSEN+1)

      ! . Calculate the energy of the molecule at the old configuration.
      EBEFORE = MONTE_CARLO_ENERGY_ONE ( CHOSEN )

      ! . Allocate space for the old coordinates.
      ALLOCATE ( OLDCRD(1:3,1:STOP-START+1) )

      ! . Save the old configuration.
      OLDENE                     = ECURRENT
      OLDCRD(1:3,1:STOP-START+1) = ATMCRD(1:3,START:STOP)

      ! . Calculate the center of the molecule.
      MOLCEN = CENTER ( ATMCRD(1:3,START:STOP), ATMMAS(START:STOP) )

      ! . Calculate the translation for the molecule.
      T(1) = 2.0_DP * TMAX * ( RANDOM ( ) - 0.5_DP )
      T(2) = 2.0_DP * TMAX * ( RANDOM ( ) - 0.5_DP )
      T(3) = 2.0_DP * TMAX * ( RANDOM ( ) - 0.5_DP )

      ! . Find the overall displacement (applying the minimum image convention).
      DR = T - BOXL * ANINT ( ( MOLCEN + T ) / BOXL, DP )

      ! . Displace the molecule.
      MOLCEN = MOLCEN + DR
      DO I = START,STOP
         ATMCRD(1:3,I) = ATMCRD(1:3,I) + DR
      END DO

      ! . Do a rotation only if the molecule has more than one particle.
      IF ( ( STOP - START + 1 ) > 1 ) THEN

         ! . Determine the rotation angle and the rotation axis.
         ANGLE = 2.0_DP * RMAX * ( RANDOM ( ) - 0.5_DP )
         AXIS  = INT ( 3.0_DP * RANDOM ( ) ) + 1
         COSA  = COS ( ANGLE )
         SINA  = SIN ( ANGLE )

         ! . Rotate the molecule about the X, Y or Z axis.
         SELECT CASE ( AXIS )
         CASE ( 1 )
            DO I = START,STOP
               DY = ATMCRD(2,I) - MOLCEN(2)
               DZ = ATMCRD(3,I) - MOLCEN(3)
               ATMCRD(2,I) = MOLCEN(2) + ( COSA * DY - SINA * DZ )
	       ATMCRD(3,I) = MOLCEN(3) + ( SINA * DY + COSA * DZ )
            END DO
         CASE ( 2 )
            DO I = START,STOP
               DZ = ATMCRD(3,I) - MOLCEN(3)
               DX = ATMCRD(1,I) - MOLCEN(1)
               ATMCRD(3,I) = MOLCEN(3) + ( COSA * DZ - SINA * DX )
               ATMCRD(1,I) = MOLCEN(1) + ( SINA * DZ + COSA * DX )
            END DO
         CASE ( 3 )
            DO I = START,STOP
               DX = ATMCRD(1,I) - MOLCEN(1)
               DY = ATMCRD(2,I) - MOLCEN(2)
               ATMCRD(1,I) = MOLCEN(1) + ( COSA * DX - SINA * DY )
               ATMCRD(2,I) = MOLCEN(2) + ( SINA * DX + COSA * DY )
            END DO
         END SELECT
      END IF

      ! . Calculate the energy of the molecule at the new configuration.
      EAFTER = MONTE_CARLO_ENERGY_ONE ( CHOSEN )

      ! . Calculate the total energy of the new configuration.
      ECURRENT = OLDENE + EAFTER - EBEFORE

      ! . Check to see if the move is rejected.
      IF ( QREJECT ( BETA * ( ECURRENT - OLDENE ) ) ) THEN

         ! . Increment NREJECT and NREJECTM.
         NREJECT  = NREJECT  + 1
         NREJECTM = NREJECTM + 1

         ! . Reactivate the old configuration.
         ECURRENT               = OLDENE
         ATMCRD(1:3,START:STOP) = OLDCRD(1:3,1:STOP-START+1)

      END IF

      ! . Deallocate temporary space.
      DEALLOCATE ( OLDCRD )

      END SUBROUTINE MOVE_MOLECULE

      !---------------------
      SUBROUTINE MOVE_VOLUME
      !---------------------

      ! . Local scalars.
      INTEGER            :: I, IMOL, START, STOP
      REAL ( KIND = DP ) :: DELTAEB, OLDENE, OLDVOL

      ! . Local arrays.
      REAL ( KIND = DP ),              DIMENSION(1:3) :: DR, OLDBOX, SCALE
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: OLDCRD

      ! . Increment the number of tries.
      NTRYV = NTRYV + 1

      ! . Allocate space for the old coordinates.
      ALLOCATE ( OLDCRD(1:3,1:NATOMS) )

      ! . Save the old configuration.
      OLDENE = ECURRENT
      OLDBOX = BOXL
      OLDVOL = VOLUME
      OLDCRD = ATMCRD

      ! . Calculate the new volume and the new box length.
      VOLUME = OLDVOL + 2.0_DP * VMAX * ( RANDOM ( ) - 0.5_DP )
      BOXL   = EXP ( LOG ( VOLUME ) / 3.0_DP )
      IF ( ABS ( PRODUCT ( BOXL ) - VOLUME ) > 1.0E-6_DP ) CALL PRINT_ERROR ( "MONTE_CARLO", "Box length/volume error." )

      ! . Calculate the scale factor for changing the coordinates of the molecular centers.
      SCALE = BOXL / OLDBOX

      ! . Translate the coordinates of the atoms in each molecule.
      DO IMOL = 1,NMOLECULES

         ! . Find the atom indices for the molecule.
         START = MOLIND(IMOL)+1
         STOP  = MOLIND(IMOL+1)

         ! . Find the corrections for the other coordinates.
         DR = ( SCALE - 1.0_DP ) * CENTER ( ATMCRD(1:3,START:STOP), ATMMAS(START:STOP) )

         ! . Alter the coordinates of the atoms (maintaining the internal geometry).
         DO I = START,STOP
            ATMCRD(1:3,I) = ATMCRD(1:3,I) + DR
         END DO

      END DO

      !. Calculate the total energy for the configuration.
      ECURRENT = MONTE_CARLO_ENERGY_FULL ( )

      ! . Calculate the quantity required for the Metropolis MC acceptance test.
      DELTAEB = BETA * ( ECURRENT - OLDENE + PFACT * ( VOLUME - OLDVOL ) - TFACT * LOG ( VOLUME / OLDVOL ) )

      ! . Check to see if the move is rejected.
      IF ( QREJECT ( DELTAEB ) ) THEN

         ! . Increment NREJECT and NREJECTV.
         NREJECT  = NREJECT  + 1
         NREJECTV = NREJECTV + 1

         ! . Reactivate the old configuration.
         ECURRENT = OLDENE
         BOXL     = OLDBOX
         VOLUME   = OLDVOL
         ATMCRD   = OLDCRD

      END IF

      ! . Deallocate temporary space.
      DEALLOCATE ( OLDCRD )

      END SUBROUTINE MOVE_VOLUME

      !-----------------------------------
      LOGICAL FUNCTION QREJECT ( DELTAEB )
      !-----------------------------------

      ! . Scalar arguments.
      REAL ( KIND = DP ), INTENT(IN) :: DELTAEB

      ! . Local parameters.
      REAL ( KIND = DP ), PARAMETER :: UNDERFLOW = 75.0_DP

      ! . Local scalars.
      LOGICAL :: QACCEPT

      ! . Initialize QACCEPT.
      QACCEPT = .FALSE.

      ! . Check to see if the move has been accepted.
      IF ( DELTAEB < UNDERFLOW ) THEN
         IF ( DELTAEB <= 0.0_DP ) THEN
            QACCEPT = .TRUE.
         ELSE IF ( EXP ( - DELTAEB ) > RANDOM ( ) ) THEN
            QACCEPT = .TRUE.
         END IF
      END IF

      ! . Set QREJECT.
      QREJECT = .NOT. QACCEPT

      END FUNCTION QREJECT

      !-------------------------------------------------------------------------
      ! . Statistics Subroutines.
      !-------------------------------------------------------------------------

      !-------------------------------------
      SUBROUTINE STATISTICS_BLOCK_ACCUMULATE
      !-------------------------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: H

      ! . Accumulate the statistics.
      H    = ECURRENT + PFACT * VOLUME
      EAV  = EAV      + ECURRENT
      EAV2 = EAV2     + ECURRENT * ECURRENT
      HAV  = HAV      + H
      HAV2 = HAV2     + H * H
      VAV  = VAV      + VOLUME
      VAV2 = VAV2     + VOLUME * VOLUME

      END SUBROUTINE STATISTICS_BLOCK_ACCUMULATE

      !--------------------------------
      SUBROUTINE STATISTICS_BLOCK_START
      !--------------------------------

      ! . Initialization.
      NREJECT = 0
      EAV     = 0.0_DP ; EAV2    = 0.0_DP
      HAV     = 0.0_DP ; HAV2    = 0.0_DP
      VAV     = 0.0_DP ; VAV2    = 0.0_DP

      END SUBROUTINE STATISTICS_BLOCK_START

      !-------------------------------
      SUBROUTINE STATISTICS_BLOCK_STOP
      !-------------------------------

      ! . Local scalars.
      INTEGER            :: NACCEPT
      REAL ( KIND = DP ) :: N

      ! . Determine the number of accepted moves.
      NACCEPT = NCONF - NREJECT

      ! . Accumulate the run statistics.
      NTOT  = NTOT  + NACCEPT
      ETOT  = ETOT  + EAV ; ETOT2 = ETOT2 + EAV2
      HTOT  = HTOT  + HAV ; HTOT2 = HTOT2 + HAV2
      VTOT  = VTOT  + VAV ; VTOT2 = VTOT2 + VAV2

      ! . Calculate the block statistics.
      N    = REAL ( NCONF, DP )
      EAV  = EAV  / N ; EAV2 = EAV2 / N - EAV * EAV
      HAV  = HAV  / N ; HAV2 = HAV2 / N - HAV * HAV
      VAV  = VAV  / N ; VAV2 = VAV2 / N - VAV * VAV

      ! . Accumulate the run block averages.
      ETOTB2 = ETOTB2 + EAV * EAV
      HTOTB2 = HTOTB2 + HAV * HAV
      VTOTB2 = VTOTB2 + VAV * VAV

      ! . Write out the statistics for the block.
      WRITE ( PRINT_LINE, "(I8)"    ) IBLOCK                   ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F16.4)" ) REAL ( NACCEPT, DP ) / N ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F16.4)" ) EAV                      ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F16.4)" ) EAV2                     ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F16.4)" ) HAV                      ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F16.4)" ) HAV2                     ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F16.4)" ) VAV                      ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F16.4)" ) VAV2                     ; CALL PRINT_TABLE_ELEMENT

      END SUBROUTINE STATISTICS_BLOCK_STOP

      !--------------------------
      SUBROUTINE STATISTICS_START
      !--------------------------

      ! . Initialize the statistics counters.
      NTOT = 0
      ETOT = 0.0_DP ; ETOT2 = 0.0_DP ; ETOTB2 = 0.0_DP
      HTOT = 0.0_DP ; HTOT2 = 0.0_DP ; HTOTB2 = 0.0_DP
      VTOT = 0.0_DP ; VTOT2 = 0.0_DP ; VTOTB2 = 0.0_DP

      ! . Write out the header.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 8, HEADER_COLOR = "#FF5500", PAGEWIDTH = 120, &
                                              VARIABLEWIDTHS = (/ 8, 16, 16, 16, 16, 16, 16, 16 /) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Run Statistics", COLSPAN = 8, HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Block",   HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Accept.", HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "<E>",     HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "<dE^2>",  HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "<H>",     HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "<dH^2>",  HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "<V>",     HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "<dV^2>",  HEADER = .TRUE. )

      END SUBROUTINE STATISTICS_START

      !-------------------------
      SUBROUTINE STATISTICS_STOP
      !-------------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: N

      ! . Do the total run statistics.
      IF ( NBLOCK > 1 ) THEN

         ! . Calculate the run statistics.
         N    = REAL ( ( NBLOCK * NCONF ), DP )
         ETOT = ETOT / N ; ETOT2 = ETOT2 / N - ETOT * ETOT
         HTOT = HTOT / N ; HTOT2 = HTOT2 / N - HTOT * HTOT
         VTOT = VTOT / N ; VTOT2 = VTOT2 / N - VTOT * VTOT

         ! . Write out the run statistics.
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "     Run" )
	 WRITE ( PRINT_LINE, "(F16.4)" ) REAL ( NTOT, DP ) / N ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F16.4)" ) ETOT                  ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F16.4)" ) ETOT2                 ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F16.4)" ) HTOT                  ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F16.4)" ) HTOT2                 ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F16.4)" ) VTOT                  ; CALL PRINT_TABLE_ELEMENT
	 WRITE ( PRINT_LINE, "(F16.4)" ) VTOT2                 ; CALL PRINT_TABLE_ELEMENT

      END IF

      ! . Write out the terminator.
      CALL PRINT_TABLE_STOP

      ! . Calculate the standard deviations in the block averages.
      IF ( NBLOCK > 1 ) THEN

         ! . Calculate the run block statistics.
         N      = REAL ( NBLOCK, DP )
         ETOTB2 = ETOTB2 / N - ETOT * ETOT
         HTOTB2 = HTOTB2 / N - HTOT * HTOT
         VTOTB2 = VTOTB2 / N - VTOT * VTOT
         N      = REAL ( NBLOCK - 1, DP )
         ETOTB2 = SQRT ( MAX ( ETOTB2 / N, 0.0_DP ) )
         HTOTB2 = SQRT ( MAX ( HTOTB2 / N, 0.0_DP ) )
         VTOTB2 = SQRT ( MAX ( VTOTB2 / N, 0.0_DP ) )

         ! . Write out the run block statistics.
         CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF5500", VARIABLEWIDTH = 16 )
	 CALL PRINT_SUMMARY_START ( "Run Block Statistics" )
	 WRITE ( PRINT_LINE, "(F16.4)" ) ETOT   ; CALL PRINT_SUMMARY_ELEMENT ( "Energy   Average"     )
	 WRITE ( PRINT_LINE, "(F16.4)" ) ETOTB2 ; CALL PRINT_SUMMARY_ELEMENT ( "Energy   Stand. Dev." )
	 WRITE ( PRINT_LINE, "(F16.4)" ) HTOT   ; CALL PRINT_SUMMARY_ELEMENT ( "Enthalpy Average"     )
	 WRITE ( PRINT_LINE, "(F16.4)" ) HTOTB2 ; CALL PRINT_SUMMARY_ELEMENT ( "Enthalpy Stand. Dev." )
	 WRITE ( PRINT_LINE, "(F16.4)" ) VTOT   ; CALL PRINT_SUMMARY_ELEMENT ( "Volume   Average"     )
	 WRITE ( PRINT_LINE, "(F16.4)" ) VTOTB2 ; CALL PRINT_SUMMARY_ELEMENT ( "Volume   Stand. Dev." )
	 CALL PRINT_SUMMARY_STOP

      END IF

      END SUBROUTINE STATISTICS_STOP

   END SUBROUTINE MONTE_CARLO

END MODULE MONTE_CARLO_SIMULATION
