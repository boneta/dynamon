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
!                         The Trajectory Analysis Module
!===============================================================================
!
! . Subroutines:
!
!   DCD_DSELF                Calculate the self-diffusion coefficient.
!   DCD_RDF_CROSS            Calculate the radial distribution function between
!                            two groups of atoms.
!   DCD_RDF_SELF             Calculate the radial distribution function for a
!                            single group of atoms.
!
!===============================================================================
MODULE DCD_ANALYSIS

! . Module declarations.
USE CONSTANTS,    ONLY : PI
USE DEFINITIONS,  ONLY : DP
USE IO_UNITS,     ONLY : OUTPUT
USE PRINTING,     ONLY : PRINT_ERROR, PRINT_NOFORMAT_START, PRINT_NOFORMAT_STOP

USE DCD_IO

IMPLICIT NONE
PRIVATE
PUBLIC :: DCD_DSELF, DCD_RDF_CROSS, DCD_RDF_SELF

!==============================================================================
CONTAINS
!==============================================================================

   !---------------------------------------------------------
   SUBROUTINE DCD_DSELF ( FILE, TSTOP, TIME_STEP, SELECTION )
   !---------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)           :: FILE
   INTEGER,               INTENT(IN)           :: TSTOP
   REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: TIME_STEP

   ! . Optional array arguments.
   LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: SELECTION

   ! . Local scalars.
   INTEGER            :: I, IATOM, IP, IFRAME, NATOMS, NFRAMES, NP
   LOGICAL            :: QSELECT
   REAL ( KIND = DP ) :: TSTEP

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: BOXDAT

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)     :: DSELF
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:)   :: ATMDAT
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:,:) :: R

   ! . Local types.
   TYPE(DCD_TYPE) :: TRAJECTORY

   ! . Activate the trajectory.
   CALL DCD_INITIALIZE ( TRAJECTORY )
   CALL DCD_ACTIVATE_READ ( FILE, TRAJECTORY )

   ! . Set some variables for the trajectory.
   NATOMS  = TRAJECTORY%NATOMS
   NFRAMES = TRAJECTORY%NFRAMES

   ! . Check NFRAMES and TSTOP.
   IF ( NFRAMES <= TSTOP ) THEN
      CALL PRINT_ERROR ( "DCD_DSELF", "The number of frames in the trajectory is too small for the input TSTOP value." )
   END IF

   ! . Get the time step.
   IF ( PRESENT ( TIME_STEP ) ) THEN
      TSTEP = TIME_STEP
   ELSE
      TSTEP = 1.0_DP
   END IF

   ! . Set the atom selection flag.
   QSELECT = PRESENT ( SELECTION )

   ! . An atom selection is present.
   IF ( QSELECT ) THEN

      ! . Check NATOMS and the dimension of SELECTION.
      IF ( SIZE ( SELECTION ) /= NATOMS ) THEN
         CALL PRINT_ERROR ( "DCD_DSELF", "The SELECTION array is incompatible with the trajectory." )
      END IF

      ! . Determine the number of particles.
      NP = COUNT ( SELECTION )

   ! . Assign the total number of particles.
   ELSE
      NP = NATOMS
   END IF

   ! . Allocate space for the temporary arrays.
   ALLOCATE ( ATMDAT(1:3,1:NATOMS), DSELF(0:TSTOP), R(1:3,1:NP,1:NFRAMES) )

   ! . Loop over the number of frames in the trajectory.
   DO IFRAME = 1,NFRAMES

      ! . Read in the data.
      CALL DCD_READ ( TRAJECTORY, ATMDAT, BOXDAT )

      ! . Contract the data if necessary.
      IF ( QSELECT ) THEN
         IP = 0
         DO IATOM = 1,NATOMS
            IF ( SELECTION(IATOM) ) THEN
               IP = IP + 1
               R(1:3,IP,IFRAME) = ATMDAT(1:3,IATOM)
            END IF
         END DO
      ! . Copy the data.
      ELSE
         R(1:3,1:NP,IFRAME) = ATMDAT
      END IF

   END DO

   ! . Deactivate the trajectory.
   CALL DCD_DEACTIVATE ( TRAJECTORY )

   ! . Calculate the self-diffusion function.
   CALL CALCULATE

   ! . Write out the values.
   CALL PRINT_NOFORMAT_START
   WRITE ( OUTPUT, "(/8('-'),A,7('-'))" ) " Self-Diffusion Function "
   WRITE ( OUTPUT, "(16X,A,15X,A)" ) "Time", "Dself"
   DO I = 0,TSTOP
      WRITE ( OUTPUT, "(2F20.4)" ) ( TSTEP * REAL ( I, DP ) ), DSELF(I)
   END DO
   WRITE ( OUTPUT, "(40('-'))" )
   CALL PRINT_NOFORMAT_STOP

   ! . Deallocate the temporary arrays.
   DEALLOCATE  ( ATMDAT, DSELF, R )

   !============================================================================
   CONTAINS
   !============================================================================

      !-------------------
      SUBROUTINE CALCULATE
      !-------------------

      ! . Local scalars.
      INTEGER            :: P, T, TMAX, T0
      REAL ( KIND = DP ) :: TOTAL

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      ! . Loop over the elements in the block.
      DO T = 0,TSTOP

         TMAX  = NFRAMES - T
         TOTAL = 0.0_DP

         ! . Loop over the coordinates in the array.
         DO T0 = 1,TMAX
            DO P = 1,NP
               DR = R(1:3,P,T0+T) - R(1:3,P,T0)
               TOTAL = TOTAL + DOT_PRODUCT ( DR, DR )
            END DO
         END DO

         ! . Add in the sums to DSELF.
         DSELF(T) = TOTAL / REAL ( 3 * NP * TMAX, DP )

      END DO

      END SUBROUTINE CALCULATE

   END SUBROUTINE DCD_DSELF

   !--------------------------------------------------------------------------------
   SUBROUTINE DCD_RDF_CROSS ( FILE, NBINS, UPPER, SELECTION1, SELECTION2, BOX_SIZE )
   !--------------------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE
   INTEGER,               INTENT(IN) :: NBINS
   REAL ( KIND = DP ),    INTENT(IN) :: UPPER

   ! . Array arguments.
   LOGICAL,            DIMENSION(:),   INTENT(IN)           :: SELECTION1, SELECTION2
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN), OPTIONAL :: BOX_SIZE

   ! . Local scalars.
   INTEGER            :: IBIN, FSTART, FSTOP, IFRAME, NATOMS, NFRAMES, NP1, NP2
   REAL ( KIND = DP ) :: FACT, LOWER, RLOWER, RUPPER, VAVE, WIDTH

   ! . Local static arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: BOXDAT

   ! . Local dynamic arrays.
   INTEGER,            ALLOCATABLE, DIMENSION(:)   :: HIST
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: GR
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: R, R1, R2

   ! . Local types.
   TYPE(DCD_TYPE) :: TRAJECTORY

   ! . Initialization.
   FSTART = 1
   LOWER  = 0.0_DP

   ! . Calculate the width of each bin.
   WIDTH = ( UPPER - LOWER ) / REAL ( NBINS, DP ) 

   ! . Activate the trajectory.
   CALL DCD_INITIALIZE ( TRAJECTORY )
   CALL DCD_ACTIVATE_READ ( FILE, TRAJECTORY )

   ! . Set some variables for the trajectory.
   NATOMS  = TRAJECTORY%NATOMS
   NFRAMES = TRAJECTORY%NFRAMES

   ! . There is crystal information on the trajectory.
   IF ( TRAJECTORY%QCRYSTAL ) THEN

      ! . BOX_SIZE must be absent.
      IF ( PRESENT ( BOX_SIZE ) ) THEN
         CALL PRINT_ERROR ( "DCD_RDF_CROSS", "There is trajectory crystal information and BOX_SIZE was specified." )
      END IF

   ! . There is no crystal information on the trajectory.
   ELSE

      ! . BOX_SIZE must be present.
      IF ( PRESENT ( BOX_SIZE ) ) THEN
         BOXDAT = BOX_SIZE
      ELSE
         CALL PRINT_ERROR ( "DCD_RDF_CROSS", "There is no trajectory crystal information and BOX_SIZE was missing." )
      END IF

   END IF

   ! . Initialize FSTOP.
   FSTOP = NFRAMES

   ! . Check the sizes of SELECTION1 and SELECTION2.
   IF ( ( SIZE ( SELECTION1 ) /= NATOMS ) .OR. ( SIZE ( SELECTION2 ) /= NATOMS ) ) THEN
      CALL PRINT_ERROR ( "DCD_RDF_CROSS", "The SELECTION arrays are incompatible with the trajectory." )
   END IF

   ! . Get the number of atoms in each selection.
   NP1 = COUNT ( SELECTION1 )
   NP2 = COUNT ( SELECTION2 )

   ! . Return if a selection is null.
   IF ( ( NP1 <= 0 ) .OR. ( NP2 <= 0 ) ) RETURN

   ! . Allocate space for the temporary arrays.
   ALLOCATE ( GR(1:NBINS), HIST(1:NBINS), R(1:3,1:NATOMS), R1(1:3,1:NP1), R2(1:3,1:NP2) )

   ! . Initialize the tabulation array.
   HIST = 0

   ! . Initialize the average volume variable.
   VAVE = 0.0_DP

   ! . Loop over the frames in the file.
   DO IFRAME = 1,FSTOP

      ! . Read in the data.
      CALL DCD_READ ( TRAJECTORY, R, BOXDAT )

      ! . Check the starting frame.
      IF ( IFRAME < FSTART ) CYCLE

      ! . Check UPPER.
      IF ( UPPER > ( MINVAL ( BOXDAT ) / 2.0_DP ) ) THEN
         CALL PRINT_ERROR ( "DCD_RDF_CROSS", "The box size is too small for upper bound of g(r)." )
      END IF

      ! . Accumulate the volume.
      VAVE = VAVE + PRODUCT ( BOXDAT )

      ! . Select the appropriate coordinates using the selection arrays.
      CALL SELECT ( SELECTION1, R, R1 )
      CALL SELECT ( SELECTION2, R, R2 )

      ! . Tabulate the distances.
      CALL TABULATE

   END DO

   ! . Calculate the average volume.
   VAVE = VAVE / REAL ( FSTOP - FSTART + 1, DP )

   ! . Calculate the density factor.
   FACT = 4.0_DP * PI * REAL ( NP1 * NP2 * ( FSTOP - FSTART + 1 ), DP ) / ( 3.0_DP * VAVE )

   ! . Calculate g(r).
   DO IBIN = 1,NBINS
      RLOWER = LOWER  + REAL ( IBIN - 1, DP ) * WIDTH
      RUPPER = RLOWER + WIDTH
      GR(IBIN) = REAL ( HIST(IBIN), DP ) / ( FACT * ( RUPPER**3 - RLOWER**3 ) )
   END DO

   ! . Write out the values.
   CALL PRINT_NOFORMAT_START
   WRITE ( OUTPUT, "(/5('-'),A,5('-'))" ) " Radial Distribution Function "
   WRITE ( OUTPUT, "(12X,A,16X,A)" ) "Distance", "g(r)"
   DO IBIN = 1,NBINS
      WRITE ( OUTPUT, "(2F20.4)" ) ( LOWER + ( 0.5_DP + REAL ( IBIN - 1, DP ) ) * WIDTH ), GR(IBIN)
   END DO
   WRITE ( OUTPUT, "(40('-'))" )
   CALL PRINT_NOFORMAT_STOP

   ! . Deactivate the trajectory.
   CALL DCD_DEACTIVATE ( TRAJECTORY )

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( GR, HIST, R, R1, R2 )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !-----------------------------------------
      SUBROUTINE SELECT ( SELECTION, RIN, ROUT )
      !-----------------------------------------

      ! . Array arguments.
      LOGICAL,            DIMENSION(:),   INTENT(IN)  :: SELECTION
      REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)  :: RIN
      REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: ROUT

      ! . Local scalars.
      INTEGER :: IATOM, IP

      ! . Select the appropriate coordinates.
      IP = 0
      DO IATOM = 1,NATOMS
         IF ( SELECTION(IATOM) ) THEN
            IP = IP + 1
            ROUT(1:3,IP) = RIN(1:3,IATOM)
         END IF
      END DO

      END SUBROUTINE SELECT

      !------------------
      SUBROUTINE TABULATE
      !------------------

      ! . Local scalars.
      INTEGER            :: BIN, I, J
      REAL ( KIND = DP ) :: RIJ

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      ! . Loop over the particles in the first selection.
      DO I = 1,NP1

         ! . Loop over the particles in the second selection.
         DO J = 1,NP2

            ! . Get the interparticle vector.
            DR(1:3) = R1(1:3,I) - R2(1:3,J)

            ! . Reset the vector by applying the minimum image convention.
            DR = DR - BOXDAT * ANINT ( DR / BOXDAT, DP )

            ! . Calculate the distance.
            RIJ = SQRT ( SUM ( DR * DR ) )

            ! . The RIJ is neither too small or too large.
            IF ( ( RIJ >= LOWER ) .AND. ( RIJ < UPPER ) ) THEN

               ! . Calculate the bin required.
               BIN = INT ( ( RIJ - LOWER ) / WIDTH ) + 1

               ! . Store the result.
               HIST(BIN) = HIST(BIN) + 1

            END IF 
         END DO
      END DO

      END SUBROUTINE TABULATE

   END SUBROUTINE DCD_RDF_CROSS

   !------------------------------------------------------------------
   SUBROUTINE DCD_RDF_SELF ( FILE, NBINS, UPPER, SELECTION, BOX_SIZE )
   !------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE
   INTEGER,               INTENT(IN) :: NBINS
   REAL ( KIND = DP ),    INTENT(IN) :: UPPER

   ! . Optional array arguments.
   LOGICAL,            DIMENSION(:),   INTENT(IN), OPTIONAL :: SELECTION
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN), OPTIONAL :: BOX_SIZE

   ! . Local scalars.
   INTEGER            :: IATOM, IBIN, IP, FSTART, FSTOP, IFRAME, NATOMS, NFRAMES, NP
   LOGICAL            :: QSELECT
   REAL ( KIND = DP ) :: FACT1, FACT2, LOWER, NIDEAL, RLOWER, RUPPER, VAVE, WIDTH

   ! . Local static arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: BOXDAT

   ! . Local dynamic arrays.
   INTEGER,            ALLOCATABLE, DIMENSION(:)   :: HIST
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: GR
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: R

   ! . Local types.
   TYPE(DCD_TYPE) :: TRAJECTORY

   ! . Initialization.
   FSTART = 1
   LOWER  = 0.0_DP

   ! . Calculate the width of each bin.
   WIDTH = ( UPPER - LOWER ) / REAL ( NBINS, DP ) 

   ! . Activate the trajectory.
   CALL DCD_INITIALIZE ( TRAJECTORY )
   CALL DCD_ACTIVATE_READ ( FILE, TRAJECTORY )

   ! . Set some variables for the trajectory.
   NATOMS  = TRAJECTORY%NATOMS
   NFRAMES = TRAJECTORY%NFRAMES

   ! . There is crystal information on the trajectory.
   IF ( TRAJECTORY%QCRYSTAL ) THEN

      ! . BOX_SIZE must be absent.
      IF ( PRESENT ( BOX_SIZE ) ) THEN
         CALL PRINT_ERROR ( "DCD_RDF_SELF", "There is trajectory crystal information and BOX_SIZE was specified." )
      END IF

   ! . There is no crystal information on the trajectory.
   ELSE

      ! . BOX_SIZE must be present.
      IF ( PRESENT ( BOX_SIZE ) ) THEN
         BOXDAT = BOX_SIZE
      ELSE
         CALL PRINT_ERROR ( "DCD_RDF_SELF", "There is no trajectory crystal information and BOX_SIZE was missing." )
      END IF

   END IF
   ! . Initialize FSTOP.
   FSTOP = NFRAMES

   ! . Set the atom selection flag.
   QSELECT = PRESENT ( SELECTION )

   ! . An atom selection is present.
   IF ( QSELECT ) THEN

      ! . Check NATOMS and the dimension of SELECTION.
      IF ( SIZE ( SELECTION ) /= NATOMS ) THEN
         CALL PRINT_ERROR ( "DCD_RDF_SELF", "The SELECTION array is incompatible with the trajectory." )
      END IF

      ! . Determine the number of particles.
      NP = COUNT ( SELECTION )

   ! . Assign the total number of particles.
   ELSE
      NP = NATOMS
   END IF

   ! . Allocate space for the temporary arrays.
   ALLOCATE ( GR(1:NBINS), HIST(1:NBINS), R(1:3,1:NATOMS) )

   ! . Initialize the tabulation array.
   HIST = 0

   ! . Initialize the average volume variable.
   VAVE = 0.0_DP

   ! . Loop over the frames in the file.
   DO IFRAME = 1,FSTOP

      ! . Read in the data.
      CALL DCD_READ ( TRAJECTORY, R, BOXDAT )

      ! . Check the starting frame.
      IF ( IFRAME < FSTART ) CYCLE

      ! . Check UPPER.
      IF ( UPPER > ( MINVAL ( BOXDAT ) / 2.0_DP ) ) THEN
         CALL PRINT_ERROR ( "DCD_RDF_SELF", "The box size is too small for upper bound of g(r)." )
      END IF

      ! . Accumulate the volume.
      VAVE = VAVE + PRODUCT ( BOXDAT )

      ! . Contract the data if necessary.
      IF ( QSELECT ) THEN
         IP = 0
         DO IATOM = 1,NATOMS
            IF ( SELECTION(IATOM) ) THEN
               IP = IP + 1
               R(1:3,IP) = R(1:3,IATOM)
            END IF
         END DO
      END IF

      ! . Tabulate the distances.
      CALL TABULATE

   END DO

   ! . Calculate the average volume.
   VAVE = VAVE / REAL ( FSTOP - FSTART + 1, DP )

   ! . Calculate some density factors.
   FACT1 = 4.0_DP * PI * REAL ( NP, DP ) / ( 3.0_DP * VAVE )
   FACT2 = REAL ( NP * ( FSTOP - FSTART + 1 ), DP )

   ! . Calculate g(r).
   DO IBIN = 1,NBINS
      RLOWER = LOWER  + REAL ( IBIN - 1, DP ) * WIDTH
      RUPPER = RLOWER + WIDTH
      NIDEAL = FACT1 * ( RUPPER**3 - RLOWER**3 )
      GR(IBIN) = REAL ( HIST(IBIN), DP ) / ( NIDEAL * FACT2 )
   END DO

   ! . Write out the values.
   CALL PRINT_NOFORMAT_START
   WRITE ( OUTPUT, "(/5('-'),A,5('-'))" ) " Radial Distribution Function "
   WRITE ( OUTPUT, "(12X,A,16X,A)" ) "Distance", "g(r)"
   DO IBIN = 1,NBINS
      WRITE ( OUTPUT, "(2F20.4)" ) ( LOWER + ( 0.5_DP + REAL ( IBIN - 1, DP ) ) * WIDTH ), GR(IBIN)
   END DO
   WRITE ( OUTPUT, "(40('-'))" )
   CALL PRINT_NOFORMAT_STOP

   ! . Deactivate the trajectory.
   CALL DCD_DEACTIVATE ( TRAJECTORY )

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( GR, HIST, R )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !------------------
      SUBROUTINE TABULATE
      !------------------

      ! . Local scalars.
      INTEGER            :: BIN, I, J
      REAL ( KIND = DP ) :: RIJ

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      ! . Outer loop over particles.
      DO I = 1,(NP-1)

         ! . Inner loop over particles.
         DO J = (I+1),NP

            ! . Get the interparticle vector.
            DR(1:3) = R(1:3,I) - R(1:3,J)

            ! . Reset the vector by applying the minimum image convention.
            DR = DR - BOXDAT * ANINT ( DR / BOXDAT, DP )

            ! . Calculate the distance.
            RIJ = SQRT ( SUM ( DR * DR ) )

            ! . The RIJ is neither too small or too large.
            IF ( ( RIJ >= LOWER ) .AND. ( RIJ < UPPER ) ) THEN

               ! . Calculate the bin required.
               BIN = INT ( ( RIJ - LOWER ) / WIDTH ) + 1

               ! . Store the result.
               HIST(BIN) = HIST(BIN) + 2

            END IF 
         END DO
      END DO

      END SUBROUTINE TABULATE

   END SUBROUTINE DCD_RDF_SELF

END MODULE DCD_ANALYSIS
