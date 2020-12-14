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
!                         The DCD Trajectory I/O Module
!===============================================================================
!
! . Basic subroutines:
!
!   DCD_ACTIVATE_READ        Activate the read trajectory.
!   DCD_ACTIVATE_WRITE       Activate the write trajectory.
!   DCD_DEACTIVATE           Deactivate a trajectory.
!   DCD_INITIALIZE           Initialize a trajectory.
!   DCD_READ                 Read a trajectory.
!   DCD_SUMMARY              Write out some information about the trajectory.
!   DCD_WRITE                Write a trajectory.
!
! . Subroutines handling more sophisticated tasks:
!
!   DCD_MERGE                Merge several trajectory files.
!
! . Notes:
!
!   The DCD trajectory format is the format used by the modeling programs CHARMM
!   and XPLOR and by the visualization program VMD.
!
!   In Fortran 90 it is not possible to initialize data inside types. This
!   problem is cured in Fortran 95 and will be made use of in a future
!   version. At present any declarations of a trajectory should be accompanied
!   by a call to DCD_INITIALIZE.
!
!===============================================================================
MODULE DCD_IO

! . Module declarations.
USE DEFINITIONS, ONLY : DP, MAX_RECORD_LENGTH, SP, VERSION
USE FILES,       ONLY : NEXT_UNIT
USE PRINTING,    ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, &
                        PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, PRINT_SUMMARY_STOP

IMPLICIT NONE
PUBLIC
PRIVATE :: DCD_SUMMARY, UNDEFINED

! . Trajectory type definitions.
TYPE DCD_TYPE
   CHARACTER ( LEN = 4 ) :: HEADER
   INTEGER               :: NATOMS, NFIXED, NFRAMES, NFREE, POSITION, UNIT
   LOGICAL               :: QCRYSTAL, QWRITE
   INTEGER, DIMENSION(:), POINTER :: FREE
END TYPE DCD_TYPE

! . The header for undefined trajectories.
CHARACTER ( LEN = 4 ), PARAMETER :: UNDEFINED = "????"

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------------------
   SUBROUTINE DCD_ACTIVATE_READ ( FILE, DCD, PRINT )
   !-----------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)    :: FILE
   TYPE(DCD_TYPE),        INTENT(INOUT) :: DCD
   
   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Local scalars.
   INTEGER :: IOSTAT
   LOGICAL :: QPRINT

   ! . Local arrays.
   INTEGER, DIMENSION(1:20) :: ICNTRL

   ! . Check the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Check to see that the trajectory has been initialized.
   IF ( DCD%HEADER /= UNDEFINED ) THEN
      CALL PRINT_ERROR ( "DCD_ACTIVATE_READ", "The trajectory has not been initialized." )
   ! . Get the next unit for the file.
   ELSE
      DCD%UNIT = NEXT_UNIT ( )
   END IF

   ! . Open the file.
   OPEN ( UNIT = DCD%UNIT, FILE = FILE, ACTION = "READ", FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, &
                                                                         STATUS = "OLD", IOSTAT = IOSTAT )
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "DCD_ACTIVATE_READ", "I/O Error.", IOSTAT )

   ! . Initialize the trajectory position and read/write flag.
   DCD%POSITION = 0
   DCD%QWRITE   = .FALSE.

   ! . Read in the header.
   READ ( DCD%UNIT ) DCD%HEADER, ICNTRL(1:20)
   READ ( DCD%UNIT )
   READ ( DCD%UNIT ) DCD%NATOMS

   ! . Fill the remaining trajectory fields.
   DCD%NFRAMES  = ICNTRL(1)
   DCD%NFIXED   = ICNTRL(9)
   DCD%NFREE    = DCD%NATOMS - DCD%NFIXED
   DCD%QCRYSTAL = ( ICNTRL(11) /= 0 )

   ! . Read in the fixed atom array if necessary.
   IF ( DCD%NFIXED > 0 ) THEN
      ALLOCATE ( DCD%FREE(1:DCD%NFREE) )
      READ ( DCD%UNIT ) DCD%FREE(1:DCD%NFREE)
   END IF

   ! . Write out a summary of the trajectory information.
   IF (QPRINT) CALL DCD_SUMMARY ( DCD )

   END SUBROUTINE DCD_ACTIVATE_READ

   !-------------------------------------------------------------------------------------------------
   SUBROUTINE DCD_ACTIVATE_WRITE ( FILE, DCD, TYPE, NATOMS, NFIXED, NFRAMES, QCRYSTAL, QFIX, PRINT )
   !-------------------------------------------------------------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)    :: FILE, TYPE
   INTEGER,               INTENT(IN)    :: NATOMS, NFIXED, NFRAMES
   TYPE(DCD_TYPE),        INTENT(INOUT) :: DCD
   

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: QCRYSTAL, PRINT

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN), OPTIONAL :: QFIX

   ! . Local scalars.
   INTEGER :: IATOM, IFREE, IOSTAT
   LOGICAL :: QPRINT

   ! . Local arrays.
   INTEGER, DIMENSION(1:20) :: ICNTRL = 0
   
   ! . Check the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Check to see that the trajectory has been initialized.
   IF ( DCD%HEADER /= UNDEFINED ) THEN
      CALL PRINT_ERROR ( "DCD_ACTIVATE_WRITE", "The trajectory has not been initialized." )
   ! . Get the next unit for the file.
   ELSE
      DCD%UNIT = NEXT_UNIT ( )
   END IF

   ! . Open the file.
   OPEN ( UNIT = DCD%UNIT, FILE = FILE, ACTION = "WRITE", FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, &
                                                                      STATUS = "UNKNOWN", IOSTAT = IOSTAT )
   IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "DCD_ACTIVATE_WRITE", "I/O Error.", IOSTAT )

   ! . Initialize the trajectory position and read/write flag.
   DCD%POSITION = 0
   DCD%QWRITE   = .TRUE.

   ! . Assign the values from the input options.
   DCD%HEADER  = TYPE
   DCD%NATOMS  = NATOMS
   DCD%NFIXED  = NFIXED
   DCD%NFRAMES = NFRAMES
   DCD%NFREE   = NATOMS - NFIXED

   ! . Set the QCRYSTAL option.
   IF ( PRESENT ( QCRYSTAL ) ) DCD%QCRYSTAL = QCRYSTAL

   ! . Create the free atom index array.
   IF ( NFIXED > 0 ) THEN

      ! . Check that QFIX exists.
      IF ( .NOT. PRESENT ( QFIX ) ) CALL PRINT_ERROR ( "DCD_ACTIVATE_WRITE", "Fixed atom array missing." )

      ! . Allocate the array.
      ALLOCATE ( DCD%FREE(1:NATOMS-NFIXED) )

      ! . Fill the array.
      IFREE = 0
      DO IATOM = 1,NATOMS
         IF ( .NOT. QFIX(IATOM) ) THEN
            IFREE = IFREE + 1
            DCD%FREE(IFREE) = IATOM
         END IF
      END DO

   END IF

   ! . Fill ICNTRL.
   ICNTRL(1)  = NFRAMES         ! Number of frames.
   ICNTRL(2)  = 1               ! Current step.
   ICNTRL(3)  = 1               ! Save frequency.
   ICNTRL(4)  = NFRAMES         ! Total number of steps.
   ICNTRL(8)  = 3 * DCD%NFREE   ! Number of degrees of freedom (not necessarily exact).
   ICNTRL(9)  = NFIXED          ! Number of fixed atoms.
   ICNTRL(10) = 0               ! Time step.
   IF ( DCD%QCRYSTAL ) THEN
      ICNTRL(11) = 1            ! Crystal present.
   ELSE
      ICNTRL(11) = 0            ! Crystal absent.
   END IF
   ICNTRL(20) = INT ( VERSION ) ! Version number.

   ! . Write out the header.
   WRITE ( DCD%UNIT ) DCD%HEADER, ICNTRL(1:20)
   WRITE ( DCD%UNIT ) 1, REPEAT ( " ", 80 )
   WRITE ( DCD%UNIT ) NATOMS
   IF ( NFIXED > 0 ) WRITE ( DCD%UNIT ) DCD%FREE(1:NATOMS-NFIXED)

   ! . Write out a summary of the trajectory information.
      IF (QPRINT) CALL DCD_SUMMARY ( DCD )

   END SUBROUTINE DCD_ACTIVATE_WRITE

   !--------------------------------
   SUBROUTINE DCD_DEACTIVATE ( DCD )
   !--------------------------------

   ! . Scalar arguments.
   TYPE(DCD_TYPE), INTENT(INOUT) :: DCD

   ! . Check that the read trajectory is active.
   IF ( DCD%UNIT <= 0 ) CALL PRINT_ERROR ( "DCD_DEACTIVATE", "The trajectory is not active." )

   ! . If the trajectory was open for writing check its position.
   IF ( DCD%QWRITE .AND. ( DCD%POSITION /= DCD%NFRAMES ) ) THEN
      CALL PRINT_ERROR ( "DCD_DEACTIVATE", "A trajectory open for writing is not full." )
   END IF

   ! . Close the trajectory file.
   CLOSE ( DCD%UNIT )

   ! . Deallocate the free atom array if present.
   IF ( DCD%NFIXED > 0 ) DEALLOCATE ( DCD%FREE )

   ! . Reinitialize the trajectory.
   CALL DCD_INITIALIZE ( DCD )

   END SUBROUTINE DCD_DEACTIVATE

   !--------------------------------
   SUBROUTINE DCD_INITIALIZE ( DCD )
   !--------------------------------

   ! . Scalar arguments.
   TYPE(DCD_TYPE), INTENT(OUT) :: DCD

   ! . Initialize the fields in the trajectory type.
   DCD%HEADER   = UNDEFINED
   DCD%NATOMS   = -1
   DCD%NFIXED   = -1
   DCD%NFRAMES  = -1
   DCD%NFREE    = -1
   DCD%POSITION = -1
   DCD%UNIT     = -1
   DCD%QCRYSTAL = .FALSE.
   DCD%QWRITE   = .FALSE.

   ! . Deallocate the pointer array.
   NULLIFY ( DCD%FREE )

   END SUBROUTINE DCD_INITIALIZE

   !------------------------------------------------------
   SUBROUTINE DCD_MERGE ( FILES_IN, FILE_OUT, SKIP_FIRST )
   !------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)           :: FILE_OUT
   LOGICAL,               INTENT(IN), OPTIONAL :: SKIP_FIRST

   ! . Array arguments.
   CHARACTER ( LEN = * ), DIMENSION(:), INTENT(IN) :: FILES_IN

   ! . Local scalars.
   INTEGER        :: I, IFILE, LOWER, NFILES, NFRAMES
   LOGICAL        :: QSKIP
   TYPE(DCD_TYPE) :: TRAJECTORY_OUT

   ! . Local arrays.
   LOGICAL,            ALLOCATABLE, DIMENSION(:)   :: QFIX
   REAL ( KIND = DP ),              DIMENSION(1:3) :: BOXDAT
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: ATMDAT
   TYPE(DCD_TYPE),     ALLOCATABLE, DIMENSION(:)   :: TRAJECTORY_IN

   ! . Check the SKIP option.
   IF ( PRESENT ( SKIP_FIRST ) ) THEN
      QSKIP = SKIP_FIRST
   ELSE
      QSKIP = .TRUE.
   END IF

   ! . Get the number of files.
   NFILES = SIZE ( FILES_IN )

   ! . Allocate the input trajectory array.
   ALLOCATE ( TRAJECTORY_IN(1:NFILES) )

   ! . Initialize the number of frames.
   NFRAMES = 0

   ! . Loop over the input trajectories.
   DO IFILE = 1,NFILES

      ! . Activate the trajectory.
      CALL DCD_INITIALIZE    (                  TRAJECTORY_IN(IFILE) )
      CALL DCD_ACTIVATE_READ ( FILES_IN(IFILE), TRAJECTORY_IN(IFILE) )

      ! . Check the trajectory parameters to ensure consistency.
      IF ( IFILE > 1 ) THEN

         ! . Check the scalar parameters.
         IF ( ( TRAJECTORY_IN(IFILE)%HEADER     /=   TRAJECTORY_IN(1)%HEADER   ) .OR. &
              ( TRAJECTORY_IN(IFILE)%NATOMS     /=   TRAJECTORY_IN(1)%NATOMS   ) .OR. &
              ( TRAJECTORY_IN(IFILE)%NFIXED     /=   TRAJECTORY_IN(1)%NFIXED   ) .OR. &
              ( TRAJECTORY_IN(IFILE)%QCRYSTAL .NEQV. TRAJECTORY_IN(1)%QCRYSTAL ) ) &
            CALL PRINT_ERROR ( "DCD_MERGE", "Inconsistent scalar data for input trajectories." )

         ! . Check the array parameters.
         IF ( TRAJECTORY_IN(1)%NFIXED > 0 ) THEN
            IF ( ANY ( TRAJECTORY_IN(IFILE)%FREE /= TRAJECTORY_IN(1)%FREE ) ) &
               CALL PRINT_ERROR ( "DCD_MERGE", "Inconsistent array data for input trajectories." )
         END IF

      END IF

      ! . Calculate the total number of frames.
      NFRAMES = NFRAMES + TRAJECTORY_IN(IFILE)%NFRAMES

   END DO

   ! . Adjust NFRAMES to account for QSKIP.
   IF ( QSKIP ) NFRAMES = NFRAMES - ( NFILES - 1 )

   ! . Allocate and initialize the atom arrays.
   ALLOCATE ( ATMDAT(1:3,1:TRAJECTORY_IN(1)%NATOMS), QFIX(1:TRAJECTORY_IN(1)%NATOMS) ) ; QFIX = .FALSE.

   ! . Set the required elements of QFIX.
   IF ( TRAJECTORY_IN(1)%NFIXED > 0 ) THEN
      QFIX = .TRUE.
      DO I = 1,TRAJECTORY_IN(1)%NFREE
         QFIX(TRAJECTORY_IN(1)%FREE(I)) = .FALSE.
      END DO
   END IF

   ! . Activate the trajectory for writing.
   CALL DCD_INITIALIZE     (           TRAJECTORY_OUT )
   CALL DCD_ACTIVATE_WRITE ( FILE_OUT, TRAJECTORY_OUT, TRAJECTORY_IN(1)%HEADER, TRAJECTORY_IN(1)%NATOMS, &
                                       TRAJECTORY_IN(1)%NFIXED, NFRAMES, TRAJECTORY_IN(1)%QCRYSTAL, QFIX )

   ! . Loop over the input trajectories again.
   DO IFILE = 1,NFILES

      ! . Determine the lower bound of the loop.
      IF ( ( IFILE == 1 ) .OR. ( .NOT. QSKIP ) ) THEN
         LOWER = 1
      ELSE
         LOWER = 2
         CALL DCD_READ ( TRAJECTORY_IN(IFILE), ATMDAT, BOXDAT )
      END IF

      ! . Copy the trajectory data to the output trajectory.
      DO I = LOWER,TRAJECTORY_IN(IFILE)%NFRAMES
         CALL DCD_READ  ( TRAJECTORY_IN(IFILE), ATMDAT, BOXDAT )
         CALL DCD_WRITE ( TRAJECTORY_OUT,       ATMDAT, BOXDAT )
      END DO

      ! . Deactivate the trajectory.
      CALL DCD_DEACTIVATE ( TRAJECTORY_IN(IFILE) )

   END DO

   ! . Deactivate the trajectory for writing.
   CALL DCD_DEACTIVATE ( TRAJECTORY_OUT )

   ! . Deallocate any temporary arrays.
   DEALLOCATE ( ATMDAT, QFIX, TRAJECTORY_IN )

   ! . Write out some information.
   CALL PRINT_PARAGRAPH ( TEXT = "Trajectory files merged into "//FILE_OUT(1:LEN_TRIM(FILE_OUT))//"." )

   END SUBROUTINE DCD_MERGE

   !------------------------------------------
   SUBROUTINE DCD_READ ( DCD, ATMDAT, BOXDAT )
   !------------------------------------------

   ! . Scalar arguments.
   TYPE(DCD_TYPE), INTENT(INOUT) :: DCD

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT)             :: ATMDAT
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(INOUT), OPTIONAL :: BOXDAT

   ! . Local scalars.
   INTEGER :: I

   ! . Local arrays.
   REAL ( KIND = DP ),              DIMENSION(1:6) :: TMPBOX
   REAL ( KIND = SP ), ALLOCATABLE, DIMENSION(:)   :: TEMP

   ! . Check the input trajectory unit number.
   IF ( DCD%UNIT <= 0 ) CALL PRINT_ERROR ( "DCD_READ", "The input trajectory is inactive." )

   ! . Check the array dimensions.
   IF ( ( SIZE ( ATMDAT, 1 ) /= 3 ) .OR. ( SIZE ( ATMDAT, 2 ) /= DCD%NATOMS ) ) THEN
      CALL PRINT_ERROR ( "DCD_READ", "Incompatible atom data array dimensions." )
   END IF

   ! . There are no more frames on the file.
   IF ( DCD%POSITION == DCD%NFRAMES ) THEN
      CALL PRINT_ERROR ( "DCD_READ", "There are no more frames on the file." )
   ! . Increment the file position.
   ELSE
      DCD%POSITION = DCD%POSITION + 1
   END IF

   ! . There is crystal information.
   IF ( DCD%QCRYSTAL ) THEN

      ! . BOXDAT must be present.
      IF ( .NOT. PRESENT ( BOXDAT ) ) CALL PRINT_ERROR ( "DCD_READ", "There is crystal information but no BOXDAT argument." )

      ! . Read in the data.
      READ ( DCD%UNIT ) TMPBOX

      ! . Fill BOXDAT correctly.
      BOXDAT(1) = TMPBOX(1)
      BOXDAT(2) = TMPBOX(3)
      BOXDAT(3) = TMPBOX(6)

   END IF

   ! . There are fixed atoms (except for the first frame).
   IF ( ( DCD%NFIXED > 0 ) .AND. ( DCD%POSITION > 1 ) ) THEN

      ! . Allocate TEMP.
      ALLOCATE ( TEMP(1:DCD%NFREE) )

      ! . Read in the data for the free atoms.
      READ ( DCD%UNIT ) TEMP ; DO I = 1,DCD%NFREE ; ATMDAT(1,DCD%FREE(I)) = REAL ( TEMP(I), DP ) ; ENDDO
      READ ( DCD%UNIT ) TEMP ; DO I = 1,DCD%NFREE ; ATMDAT(2,DCD%FREE(I)) = REAL ( TEMP(I), DP ) ; ENDDO
      READ ( DCD%UNIT ) TEMP ; DO I = 1,DCD%NFREE ; ATMDAT(3,DCD%FREE(I)) = REAL ( TEMP(I), DP ) ; ENDDO

   ! . There are no fixed atoms.
   ELSE

      ! . Allocate TEMP.
      ALLOCATE ( TEMP(1:DCD%NATOMS) )

      ! . Read in the data for all the atoms.
      READ ( DCD%UNIT ) TEMP ; ATMDAT(1,1:DCD%NATOMS) = REAL ( TEMP, DP )
      READ ( DCD%UNIT ) TEMP ; ATMDAT(2,1:DCD%NATOMS) = REAL ( TEMP, DP )
      READ ( DCD%UNIT ) TEMP ; ATMDAT(3,1:DCD%NATOMS) = REAL ( TEMP, DP )

   END IF

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( TEMP )

   END SUBROUTINE DCD_READ

   !-----------------------------
   SUBROUTINE DCD_SUMMARY ( DCD )
   !-----------------------------

   ! . Scalar arguments.
   TYPE(DCD_TYPE), INTENT(IN) :: DCD

   ! . Write out the header.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF00FF", VARIABLEWIDTH = 16 )
   CALL PRINT_SUMMARY_START ( "DCD Trajectory Summary" )

   ! . Write out the trajectory details.
   CALL PRINT_SUMMARY_ELEMENT ( "Trajectory Type", TEXT = "        "//DCD%HEADER )
   WRITE ( PRINT_LINE, "(L16)" ) DCD%QWRITE   ; CALL PRINT_SUMMARY_ELEMENT ( "Write Activation"    )
   WRITE ( PRINT_LINE, "(I16)" ) DCD%NATOMS   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Atoms"     )
   WRITE ( PRINT_LINE, "(I16)" ) DCD%NFRAMES  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Frames"    )
   WRITE ( PRINT_LINE, "(I16)" ) DCD%NFIXED   ; CALL PRINT_SUMMARY_ELEMENT ( "Fixed Atoms"         )
   WRITE ( PRINT_LINE, "(L16)" ) DCD%QCRYSTAL ; CALL PRINT_SUMMARY_ELEMENT ( "Crystal Information" )

   ! . Write out the terminator.
   CALL PRINT_SUMMARY_STOP

   END SUBROUTINE DCD_SUMMARY

   !-------------------------------------------
   SUBROUTINE DCD_WRITE ( DCD, ATMDAT, BOXDAT )
   !-------------------------------------------

   ! . Scalar arguments.
   TYPE(DCD_TYPE), INTENT(INOUT) :: DCD

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)           :: ATMDAT
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN), OPTIONAL :: BOXDAT

   ! . Local scalars.
   INTEGER :: I

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:6) :: TMPBOX

   ! . Check the output trajectory unit number.
   IF ( DCD%UNIT <= 0 ) CALL PRINT_ERROR ( "DCD_WRITE", "The write trajectory is inactive." )

   ! . Check the array dimensions.
   IF ( ( SIZE ( ATMDAT, 1 ) /= 3 ) .OR. ( SIZE ( ATMDAT, 2 ) /= DCD%NATOMS ) ) THEN
      CALL PRINT_ERROR ( "DCD_WRITE", "Incompatible atom data array dimensions." )
   END IF

   ! . The file is full.
   IF ( DCD%POSITION == DCD%NFRAMES ) THEN
      CALL PRINT_ERROR ( "DCD_WRITE", "The file is full." )
   ! . Increment the file position.
   ELSE
      DCD%POSITION = DCD%POSITION + 1
   END IF

   ! . There is crystal information.
   IF ( DCD%QCRYSTAL ) THEN

      ! . BOXDAT must be present.
      IF ( .NOT. PRESENT ( BOXDAT ) ) CALL PRINT_ERROR ( "DCD_WRITE", "There is crystal information but no BOXDAT argument." )

      ! . Fill TMPBOX correctly.
      TMPBOX    = 0.0_DP
      TMPBOX(1) = BOXDAT(1)
      TMPBOX(3) = BOXDAT(2)
      TMPBOX(6) = BOXDAT(3)

      ! . Write out the data.
      WRITE ( DCD%UNIT ) TMPBOX

   END IF

   ! . There are fixed atoms (except for the first frame).
   IF ( ( DCD%NFIXED > 0 ) .AND. ( DCD%POSITION > 1 ) ) THEN

      ! . Write out the data for the free atoms.
      WRITE ( DCD%UNIT ) ( REAL ( ATMDAT(1,DCD%FREE(I)), SP ), I = 1,DCD%NFREE )
      WRITE ( DCD%UNIT ) ( REAL ( ATMDAT(2,DCD%FREE(I)), SP ), I = 1,DCD%NFREE )
      WRITE ( DCD%UNIT ) ( REAL ( ATMDAT(3,DCD%FREE(I)), SP ), I = 1,DCD%NFREE )

   ! . There are no fixed atoms.
   ELSE

      ! . Write out the data for all the atoms.
      WRITE ( DCD%UNIT ) ( REAL ( ATMDAT(1,I), SP ), I = 1,DCD%NATOMS )
      WRITE ( DCD%UNIT ) ( REAL ( ATMDAT(2,I), SP ), I = 1,DCD%NATOMS )
      WRITE ( DCD%UNIT ) ( REAL ( ATMDAT(3,I), SP ), I = 1,DCD%NATOMS )

   END IF

   END SUBROUTINE DCD_WRITE

END MODULE DCD_IO
