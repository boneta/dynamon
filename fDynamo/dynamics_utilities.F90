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
!                         The Dynamics Utilities Module
!===============================================================================
!
! . Subroutines:
!
!   DYNAMICS_INITIALIZE                 Initialize the dynamics options.
!   DYNAMICS_OPTIONS                    Set some basic dynamics options.
!   DYNAMICS_PRINT                      Print out dynamics information.
!   DYNAMICS_SAVE_ACTIVATE              Activate dynamics trajectories.
!   DYNAMICS_SAVE_DEACTIVATE            Deactivate dynamics trajectories.
!   DYNAMICS_SAVE_WRITE                 Write dynamics trajectories.
!   DYNAMICS_STATISTICS_ACCUMULATE      Accumulate the statistics.
!   DYNAMICS_STATISTICS_INITIALIZE      Initialize statistics accumulation.
!   DYNAMICS_STATISTICS_PRINT           Print the statistics.
!   ENERGY_AND_ACCELERATIONS            Calculate the energy and accelerations.
!
!===============================================================================
MODULE DYNAMICS_UTILITIES

! . Module declarations.
USE CONSTANTS,        ONLY : KJMOL_TO_AMUA2PS2
USE DEFINITIONS,      ONLY : DP, LINE_LENGTH
USE PRINTING,         ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
                             PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT,          &
			     PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP

USE ATOMS,            ONLY : ATMCRD, ATMFIX, ATMMAS, NATOMS, NFIXED
USE DCD_IO,           ONLY : DCD_ACTIVATE_WRITE, DCD_DEACTIVATE, DCD_INITIALIZE, DCD_TYPE, DCD_WRITE
USE POTENTIAL_ENERGY, ONLY : ATMDER, ETOTAL, GRADIENT
USE SYMMETRY,         ONLY : BOXL, QBOX
USE VELOCITY,         ONLY : ATMVEL, VELOCITY_PROJECT

USE FILES

#ifdef	HACK_URC
USE URC
#endif

#ifdef	HACK_GH
USE GH, ONLY : GH_FORCES
#endif

IMPLICIT NONE
PRIVATE
PUBLIC :: DYNAMICS_DELTA, DYNAMICS_PRFRQ, DYNAMICS_STEPS, DYNAMICS_TIME, &
          DYNAMICS_INITIALIZE, DYNAMICS_OPTIONS, DYNAMICS_SAVE_ACTIVATE, DYNAMICS_SAVE_DEACTIVATE, DYNAMICS_SAVE_WRITE, &
          DYNAMICS_STATISTICS_ACCUMULATE, DYNAMICS_STATISTICS_INITIALIZE, DYNAMICS_STATISTICS_PRINT, &
	  ENERGY_AND_ACCELERATIONS
#ifndef PGPC
SAVE
#endif

! . Public scalars.
INTEGER                         :: DYNAMICS_PRFRQ = 0,      DYNAMICS_STEPS = 0
REAL ( KIND = DP )              :: DYNAMICS_DELTA = 0.0_DP, DYNAMICS_TIME  = 0.0_DP

! . Private scalars.
CHARACTER ( LEN = LINE_LENGTH ) :: CFILE, VFILE
INTEGER                         :: SAVFRQ = 0
LOGICAL                         :: QCFILE = .FALSE., QVFILE = .FALSE., QRFILE = .FALSE.
TYPE(DCD_TYPE)                  :: CTRAJECTORY, VTRAJECTORY

! . Accumulators.
REAL ( KIND = DP ) :: EKE1, EKE2, EPOT1, EPOT2, ETOT1, ETOT2, PRES1, PRES2, TEMP1, TEMP2, VOLU1, VOLU2

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------
   SUBROUTINE DYNAMICS_INITIALIZE
   !-----------------------------

   ! . Initialize the dynamics options.
   CFILE = " "
   VFILE = " "

   DYNAMICS_PRFRQ = 0
   DYNAMICS_STEPS = 0
   SAVFRQ         = 0

   QCFILE         = .FALSE.
   QVFILE         = .FALSE.

   DYNAMICS_DELTA = 0.0_DP
   DYNAMICS_TIME  = 0.0_DP

   END SUBROUTINE DYNAMICS_INITIALIZE

   !----------------------------------------------------------------------------------------------------------------
   SUBROUTINE DYNAMICS_OPTIONS ( TIME_STEP, STEPS, PRINT_FREQUENCY, SAVE_FREQUENCY, COORDINATE_FILE, VELOCITY_FILE )
   !----------------------------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: PRINT_FREQUENCY, STEPS
   REAL ( KIND = DP ), INTENT(IN) :: TIME_STEP

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: COORDINATE_FILE, VELOCITY_FILE
   INTEGER,               INTENT(IN), OPTIONAL :: SAVE_FREQUENCY

   ! . Initialize the dynamics data structure.
   CALL DYNAMICS_INITIALIZE

   ! . Initialize the basic options.
   DYNAMICS_DELTA = TIME_STEP
   DYNAMICS_PRFRQ = PRINT_FREQUENCY
   DYNAMICS_STEPS = STEPS

   ! . Calculate the total simulation time.
   DYNAMICS_TIME = DYNAMICS_DELTA * REAL ( DYNAMICS_STEPS, DP )

   IF( DYNAMICS_PRFRQ < 0 ) THEN
       DYNAMICS_PRFRQ = - DYNAMICS_PRFRQ
       QRFILE = .TRUE.
   END IF

   ! . Check DYNAMICS_PRFRQ.
   IF ( ( DYNAMICS_PRFRQ < 0 ) .OR. ( DYNAMICS_PRFRQ > DYNAMICS_STEPS ) ) DYNAMICS_PRFRQ = 0

   ! . Get the save frequency.
   SAVFRQ = 0
   IF ( PRESENT ( SAVE_FREQUENCY ) ) SAVFRQ = SAVE_FREQUENCY

   ! . Initialize the input options.
   QCFILE = PRESENT ( COORDINATE_FILE )
   QVFILE = PRESENT (   VELOCITY_FILE )

   ! . Save the file names.
   IF ( QCFILE ) CFILE = COORDINATE_FILE
   IF ( QVFILE ) VFILE =   VELOCITY_FILE

   ! . Check the input save options.
   IF ( ( SAVFRQ < 0 ) .OR. ( SAVFRQ > DYNAMICS_STEPS ) ) SAVFRQ = 0
   IF ( ( SAVFRQ > 0 ) .AND. .NOT. ( QCFILE .OR. QVFILE ) ) THEN
      CALL PRINT_ERROR ( "DYNAMICS_OPTIONS", "Missing coordinate or velocity trajectory file name." )
   END IF

   ! . Print out the options.
   IF ( DYNAMICS_PRFRQ > 0 ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#005500", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Dynamics Options" )
      WRITE ( PRINT_LINE, "(G16.6)" ) DYNAMICS_DELTA ; CALL PRINT_SUMMARY_ELEMENT ( "Time Step (ps)"  )
      WRITE ( PRINT_LINE, "(G16.6)" ) DYNAMICS_TIME  ; CALL PRINT_SUMMARY_ELEMENT ( "Total Time"      )
      WRITE ( PRINT_LINE, "(I16)"   ) DYNAMICS_STEPS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Steps" )
      WRITE ( PRINT_LINE, "(I16)"   ) DYNAMICS_PRFRQ ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency" )
      WRITE ( PRINT_LINE, "(I16)"   ) SAVFRQ         ; CALL PRINT_SUMMARY_ELEMENT ( "Save Frequency"  )
      IF ( SAVFRQ > 0 ) THEN
         WRITE ( PRINT_LINE, "(L16)" ) QCFILE ; CALL PRINT_SUMMARY_ELEMENT ( "Save Coordinates" )
         WRITE ( PRINT_LINE, "(L16)" ) QVFILE ; CALL PRINT_SUMMARY_ELEMENT ( "Save Velocities"  )
      END IF
      CALL PRINT_SUMMARY_STOP
   END IF

   END SUBROUTINE DYNAMICS_OPTIONS

   !----------------------------------------------------------------------------
   SUBROUTINE DYNAMICS_PRINT ( OPTION, TIME, ETOT, EKE, EPOT, TEMP, PRES, VOLU )
   !----------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: OPTION
   REAL ( KIND = DP ), INTENT(IN) :: EKE, EPOT, ETOT, PRES, TEMP, TIME, VOLU

   ! . Local scalars.
   INTEGER :: NCOLUMNS, I

   ! . Start the table.
   IF ( OPTION == 0 ) THEN
      IF ( ( PRES /= 0.0_DP ) .AND. ( VOLU /= 0.0_DP ) ) THEN
         NCOLUMNS = 7
      ELSE
         NCOLUMNS = 5
      END IF
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = NCOLUMNS, HEADER_COLOR = "#005500", PAGEWIDTH = ( 16 * NCOLUMNS ) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Dynamics Results (ps,kJ mol^-1,K)", COLSPAN = NCOLUMNS, HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Time",           HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Total Energy",   HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Kinetic Energy", HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Potent. Energy", HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Temperature",    HEADER = .TRUE. )
      IF ( ( PRES /= 0.0_DP ) .AND. ( VOLU /= 0.0_DP ) ) THEN
         CALL PRINT_TABLE_ELEMENT ( TEXT = "Press. (atm.)", HEADER = .TRUE. )
         CALL PRINT_TABLE_ELEMENT ( TEXT = "Volume (A^3)",  HEADER = .TRUE. )
      END IF
   END IF

   ! . Print out the first entry in the table (depends on the option).
   SELECT CASE ( OPTION )
   CASE ( -1 ) ; CALL PRINT_TABLE_ELEMENT ( TEXT = "Averages",       ALIGN = "LEFT" )
   CASE ( -2 ) ; CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS Deviations", ALIGN = "LEFT" )
   CASE DEFAULT ; WRITE ( PRINT_LINE, "(G16.6)" ) TIME ; CALL PRINT_TABLE_ELEMENT
   END SELECT

   ! . Write out the remaining results.
   WRITE ( PRINT_LINE, "(G16.6)" ) ETOT ; CALL PRINT_TABLE_ELEMENT
   WRITE ( PRINT_LINE, "(G16.6)" ) EKE  ; CALL PRINT_TABLE_ELEMENT
   WRITE ( PRINT_LINE, "(G16.6)" ) EPOT ; CALL PRINT_TABLE_ELEMENT
   WRITE ( PRINT_LINE, "(G16.6)" ) TEMP ; CALL PRINT_TABLE_ELEMENT
   IF ( ( PRES /= 0.0_DP ) .AND. ( VOLU /= 0.0_DP ) ) THEN
      WRITE ( PRINT_LINE, "(G16.6)" ) PRES ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(G16.6)" ) VOLU ; CALL PRINT_TABLE_ELEMENT
   END IF

   ! . Print the terminator.
   IF ( OPTION == -2 ) CALL PRINT_TABLE_STOP
   
   IF( QRFILE ) THEN
       I = NEXT_UNIT()
       OPEN( I, FILE = "restart.bin", ACTION = "WRITE", STATUS = "REPLACE", FORM = "UNFORMATTED" )
       WRITE( I ) ATMCRD
       WRITE( I ) ATMVEL
       CLOSE( I )
   END IF

   END SUBROUTINE DYNAMICS_PRINT

   !--------------------------------
   SUBROUTINE DYNAMICS_SAVE_ACTIVATE
   !--------------------------------

   ! . Activate the trajectories if necessary.
   IF ( SAVFRQ > 0 ) THEN
      IF ( QCFILE ) THEN
         CALL DCD_INITIALIZE ( CTRAJECTORY )
         CALL DCD_ACTIVATE_WRITE ( CFILE, CTRAJECTORY, "CORD", NATOMS, NFIXED, (DYNAMICS_STEPS/SAVFRQ+1), &
	                                                                   QCRYSTAL = QBOX, QFIX = ATMFIX )
         CALL DCD_WRITE ( CTRAJECTORY, ATMCRD, BOXL )
      END IF
      IF ( QVFILE ) THEN
         CALL DCD_INITIALIZE ( VTRAJECTORY )
         CALL DCD_ACTIVATE_WRITE ( VFILE, VTRAJECTORY, "VELO", NATOMS, NFIXED, (DYNAMICS_STEPS/SAVFRQ+1), &
	                                                                                    QFIX = ATMFIX )
         CALL DCD_WRITE ( VTRAJECTORY, ATMVEL )
      END IF
   END IF

   END SUBROUTINE DYNAMICS_SAVE_ACTIVATE

   !----------------------------------
   SUBROUTINE DYNAMICS_SAVE_DEACTIVATE
   !----------------------------------

   ! . Close the trajectories.
   IF ( SAVFRQ > 0 ) THEN
      IF ( QCFILE ) CALL DCD_DEACTIVATE ( CTRAJECTORY )
      IF ( QVFILE ) CALL DCD_DEACTIVATE ( VTRAJECTORY )
   END IF

   END SUBROUTINE DYNAMICS_SAVE_DEACTIVATE

   !---------------------------------------
   SUBROUTINE DYNAMICS_SAVE_WRITE ( ISTEP )
   !---------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: ISTEP

   ! . Save some data if necessary.
   IF ( SAVFRQ > 0 ) THEN
      IF ( MOD ( ISTEP, SAVFRQ ) == 0 ) THEN
         IF ( QCFILE ) CALL DCD_WRITE ( CTRAJECTORY, ATMCRD, BOXL )
         IF ( QVFILE ) CALL DCD_WRITE ( VTRAJECTORY, ATMVEL )
      END IF
   END IF

   END SUBROUTINE DYNAMICS_SAVE_WRITE

   !-------------------------------------------------------------------------------------
   SUBROUTINE DYNAMICS_STATISTICS_ACCUMULATE ( ISTEP, TIME, EKE, EPOT, TEMP, PRES, VOLU )
   !-------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: ISTEP
   REAL ( KIND = DP ), INTENT(IN) :: EKE, EPOT, PRES, TEMP, TIME, VOLU

   ! . Local scalars.
   REAL ( KIND = DP ) :: ETOT

   ! . Do nothing if DYNAMICS_PRFRQ is zero.
   IF ( DYNAMICS_PRFRQ <= 0 ) RETURN

   ! . Calculate the total energy.
   ETOT = EKE + EPOT

   ! . Accumulate the average and fluctuation variables.
   EKE1  = EKE1  + EKE  ; EKE2  = EKE2  + EKE  * EKE
   EPOT1 = EPOT1 + EPOT ; EPOT2 = EPOT2 + EPOT * EPOT
   ETOT1 = ETOT1 + ETOT ; ETOT2 = ETOT2 + ETOT * ETOT
   PRES1 = PRES1 + PRES ; PRES2 = PRES2 + PRES * PRES
   TEMP1 = TEMP1 + TEMP ; TEMP2 = TEMP2 + TEMP * TEMP
   VOLU1 = VOLU1 + VOLU ; VOLU2 = VOLU2 + VOLU * VOLU

   ! . Output the results for the step.
   IF ( MOD ( ISTEP, DYNAMICS_PRFRQ ) == 0 ) CALL DYNAMICS_PRINT ( ISTEP, TIME, ETOT, EKE, EPOT, TEMP, PRES, VOLU )

   END SUBROUTINE DYNAMICS_STATISTICS_ACCUMULATE

   !------------------------------------------------------------------------------
   SUBROUTINE DYNAMICS_STATISTICS_INITIALIZE ( TIME, EKE, EPOT, TEMP, PRES, VOLU )
   !------------------------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: EKE, EPOT, PRES, TEMP, TIME, VOLU

   ! . Local scalars.
   REAL ( KIND = DP ) :: ETOT

   ! . Do nothing if DYNAMICS_PRFRQ is zero.
   IF ( DYNAMICS_PRFRQ <= 0 ) RETURN

   ! . Calculate the total energy.
   ETOT = EKE + EPOT

   ! . Initialize the average and fluctuation accumulators.
   EKE1  = EKE  ; EKE2  = EKE  * EKE
   EPOT1 = EPOT ; EPOT2 = EPOT * EPOT
   ETOT1 = ETOT ; ETOT2 = ETOT * ETOT
   PRES1 = PRES ; PRES2 = PRES * PRES
   TEMP1 = TEMP ; TEMP2 = TEMP * TEMP
   VOLU1 = VOLU ; VOLU2 = VOLU * VOLU

   ! . Start the printing.
   CALL DYNAMICS_PRINT ( 0, TIME, ETOT, EKE, EPOT, TEMP, PRES, VOLU )

   END SUBROUTINE DYNAMICS_STATISTICS_INITIALIZE

   !-----------------------------------
   SUBROUTINE DYNAMICS_STATISTICS_PRINT
   !-----------------------------------

   ! . Local scalars.
   REAL ( KIND = DP ) :: FACT

   ! . Do nothing if DYNAMICS_PRFRQ is zero.
   IF ( DYNAMICS_PRFRQ <= 0 ) RETURN

   ! . Calculate the averages and fluctuations.
   FACT  = REAL ( DYNAMICS_STEPS + 1, DP )
   EKE1  = EKE1  / FACT ; EKE2  = SQRT ( MAX ( ( EKE2  / FACT - EKE1  * EKE1  ), 0.0_DP ) )
   EPOT1 = EPOT1 / FACT ; EPOT2 = SQRT ( MAX ( ( EPOT2 / FACT - EPOT1 * EPOT1 ), 0.0_DP ) )
   ETOT1 = ETOT1 / FACT ; ETOT2 = SQRT ( MAX ( ( ETOT2 / FACT - ETOT1 * ETOT1 ), 0.0_DP ) )
   PRES1 = PRES1 / FACT ; PRES2 = SQRT ( MAX ( ( PRES2 / FACT - PRES1 * PRES1 ), 0.0_DP ) )
   TEMP1 = TEMP1 / FACT ; TEMP2 = SQRT ( MAX ( ( TEMP2 / FACT - TEMP1 * TEMP1 ), 0.0_DP ) )
   VOLU1 = VOLU1 / FACT ; VOLU2 = SQRT ( MAX ( ( VOLU2 / FACT - VOLU1 * VOLU1 ), 0.0_DP ) )

   ! . Print the averages and fluctuations.
   CALL DYNAMICS_PRINT ( -1, DYNAMICS_TIME, ETOT1, EKE1, EPOT1, TEMP1, PRES1, VOLU1 )
   CALL DYNAMICS_PRINT ( -2, DYNAMICS_TIME, ETOT2, EKE2, EPOT2, TEMP2, PRES2, VOLU2 )

   END SUBROUTINE DYNAMICS_STATISTICS_PRINT

   !---------------------------------------------------
   SUBROUTINE ENERGY_AND_ACCELERATIONS ( EPOT, ATMACC )
   !---------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: EPOT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(OUT) :: ATMACC

   ! . Local scalars.
   INTEGER :: I

   ! . Calculate the energy and first derivatives.
   CALL GRADIENT ( PRINT = .FALSE. )

   ! . Set the total energy.
   EPOT = ETOTAL

#ifdef	HACK_URC
   CALL GRADIENT_URC( EPOT, ATMDER )
#endif

#ifdef	HACK_GH
   CALL GH_FORCES
#endif

   ! . Convert the first derivatives into accelerations.
   DO I = 1,NATOMS
      IF ( .NOT. ATMFIX(I) ) ATMACC(1:3,I) = - KJMOL_TO_AMUA2PS2 * ATMDER(1:3,I) / ATMMAS(I)
   END DO

   ! . Project out the translational degrees of freedom.
   CALL VELOCITY_PROJECT ( ATMACC )

   END SUBROUTINE ENERGY_AND_ACCELERATIONS

END MODULE DYNAMICS_UTILITIES

