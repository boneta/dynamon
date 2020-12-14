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
!                                The WHAM Module
!===============================================================================
!
! . Subroutines:
!
!   WHAM_ANALYSE                   Analyse constraint data using the WHAM
!                                  technique.
!
! . Notes:
!
!   Only one-dimensional analyses can be performed at present.
!
!===============================================================================
MODULE WHAM

! . Module declarations.
USE CONSTANTS,   ONLY : R
USE DEFINITIONS, ONLY : DP
USE FILES,       ONLY : NEXT_UNIT
USE PRINTING,    ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH, PRINT_TABLE_ELEMENT, &
			PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP

IMPLICIT NONE
PRIVATE
PUBLIC :: WHAM_ANALYSE

!===============================================================================
CONTAINS
!===============================================================================

   !-------------------------------------------------------------------------
   SUBROUTINE WHAM_ANALYSE ( DATA, NBINS, T, CONVERGENCE, ITERATIONS, PRINT )
   !-------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: NBINS
   REAL ( KIND = DP ), INTENT(IN) :: T

   ! . Array arguments.
   CHARACTER ( LEN = * ), DIMENSION(:), INTENT(IN) :: DATA

   ! . Optional scalar arguments.
   INTEGER,            INTENT(IN), OPTIONAL :: ITERATIONS
   LOGICAL,            INTENT(IN), OPTIONAL :: PRINT
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: CONVERGENCE

   ! . Local scalars.
   INTEGER            :: I, NWINDOWS
   LOGICAL            :: QCONVERGE
   REAL ( KIND = DP ) :: RT

   ! . Local arrays.
   INTEGER,            ALLOCATABLE, DIMENSION(:)   :: FREQUENCY
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: COORDINATE, F, FORCE, PMF, RAWDATA, REFERENCE, RHO
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: UMBPOT, UMBPOTE

   ! . Determine the number of windows.
   NWINDOWS = SIZE ( DATA )

   ! . Allocate some temporary arrays.
   ALLOCATE ( COORDINATE(1:NBINS), F(1:NWINDOWS), FORCE(1:NWINDOWS), FREQUENCY(1:NWINDOWS), &
                       PMF(1:NBINS), RAWDATA(1:NBINS), REFERENCE(1:NWINDOWS), RHO(1:NBINS), &
                                    UMBPOT(1:NWINDOWS,1:NBINS), UMBPOTE(1:NWINDOWS,1:NBINS) )

   ! . Read in and bin the data for the windows.
   CALL PROCESS_DATA

   ! . Calculate RT.
   RT = R * T

   ! . Solve the WHAM equations.
   CALL WHAM_SOLVE

   ! . The calculation converged.
   IF ( QCONVERGE ) THEN

      ! . Write out a header.
      WRITE ( PRINT_LINE, "(A,F8.2,A)" ) "PMF generated at a temperature of ", T, "K."
      CALL PRINT_PARAGRAPH

      ! . Write out the window information.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#AAAA00", VARIABLEWIDTHS = (/ 10, 10, 20, 20, 20 /) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Umbrella Potential Windows", COLSPAN = 5, HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Window",         HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Points",         HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Reference",      HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Force Constant", HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Free Energy",    HEADER = .TRUE. )
      DO I = 1,NWINDOWS
         WRITE ( PRINT_LINE, "(I10)"    ) I            ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(I10)"    ) FREQUENCY(I) ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F20.10)" ) REFERENCE(I) ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F20.10)" ) FORCE(I)     ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F20.10)" ) F(I)         ; CALL PRINT_TABLE_ELEMENT
      END DO
      CALL PRINT_TABLE_STOP

      ! . Write out the PMF.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 3, HEADER_COLOR = "#AAAA00", PAGEWIDTH = 60, VARIABLEWIDTHS = (/ 20, 20, 20 /) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Umbrella Potential PMF", COLSPAN = 3, HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Coordinate", HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Density",    HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "PMF",        HEADER = .TRUE. )
      DO I = 1,NBINS
         WRITE ( PRINT_LINE, "(F20.10)" ) COORDINATE(I) ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G20.10)" ) RHO(I)        ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F20.10)" ) PMF(I)        ; CALL PRINT_TABLE_ELEMENT
      END DO
      CALL PRINT_TABLE_STOP

   ! . Unable to obtain convergence.
   ELSE
      CALL PRINT_PARAGRAPH ( TEXT = "Unable to obtain convergence of the WHAM equations." )
   END IF

   ! . Deallocate all temporary arrays.
   DEALLOCATE ( COORDINATE, F, FORCE, FREQUENCY, PMF, RAWDATA, REFERENCE, RHO, UMBPOT, UMBPOTE )

   !=============================================================================
   CONTAINS
   !=============================================================================

      !----------------------
      SUBROUTINE PROCESS_DATA
      !----------------------

      ! . Local parameters.
      REAL ( KIND = DP ), PARAMETER :: SMALL = 1.0E-6_DP

      ! . Local scalars.
      INTEGER            :: IBIN, IOSTAT, IWINDOW, N, NPOINTS, UNIT
      REAL ( KIND = DP ) :: DUM1, DUM2, VALUE, VMAX, VMIN, WIDTH

      ! . Initialize the total number of data points.
      NPOINTS = 0

      ! . Initialize the VMAX and VMIN variables.
      VMAX = - HUGE ( 0.0_DP )
      VMIN =   HUGE ( 0.0_DP )

      ! . Loop over the data files.
      DO IWINDOW = 1,NWINDOWS

         ! . Open the file.
         UNIT = NEXT_UNIT ( )
         OPEN ( UNIT, FILE = TRIM ( DATA(IWINDOW) ), FORM = "FORMATTED", STATUS = "OLD", IOSTAT = IOSTAT )
         IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "WHAM_ANALYSE", "I/O Error on data file.", IOSTAT )

         ! . Read in the force constant and reference value for the window.
         READ ( UNIT, * ) FORCE(IWINDOW), REFERENCE(IWINDOW)

         ! . Read in remaining data until the end of the file.
         N = 0
         DO

            ! . Read a line.
            READ ( UNIT, *, END = 10 ) VALUE
            N = N + 1

            ! . Check the value.
            VMAX = MAX ( VMAX, VALUE )
            VMIN = MIN ( VMIN, VALUE )

         END DO
         10 CONTINUE

         ! . Set the frequency.
         FREQUENCY(IWINDOW) = N

         ! . Accumulate the total number of points.
         NPOINTS = NPOINTS + N

         ! . Close the file.
         CLOSE ( UNIT )

      END DO

      ! . Adjust VMAX and VMIN by a small amount.
      VMAX = VMAX + SMALL
      VMIN = VMIN - SMALL

      ! . Determine the width of each bin.
      WIDTH = ( VMAX - VMIN ) / REAL ( NBINS, DP )

      ! . Determine the COORDINATE array.
      DO IBIN = 1,NBINS
         COORDINATE(IBIN) = VMIN + WIDTH * ( REAL ( IBIN, DP ) - 0.5_DP )
      END DO

      ! . Initialize RAWDATA.
      RAWDATA = 0.0_DP

      ! . Loop over the data files.
      DO IWINDOW = 1,NWINDOWS

         ! . Open the file.
         UNIT = NEXT_UNIT ( )
         OPEN ( UNIT, FILE = TRIM ( DATA(IWINDOW) ), FORM = "FORMATTED", STATUS = "OLD", IOSTAT = IOSTAT )
         IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "WHAM_ANALYSE", "I/O Error on data file.", IOSTAT )

         ! . Read in the force constant and reference value for the window.
         READ ( UNIT, * ) DUM1, DUM2

         ! . Read in remaining data until the end of the file.
         DO

            ! . Read a line.
            READ ( UNIT, *, END = 20 ) VALUE

            ! . Find the appropriate bin.
            IBIN = INT ( ( VALUE - VMIN ) / WIDTH ) + 1

            ! . Increment RAWDATA.
            RAWDATA(IBIN) = RAWDATA(IBIN) + 1.0_DP

         END DO
         20 CONTINUE

         ! . Close the file.
         CLOSE ( UNIT )

      END DO

      END SUBROUTINE PROCESS_DATA

      !--------------------
      SUBROUTINE WHAM_SOLVE
      !--------------------

      ! . Local scalars.
      INTEGER            :: I, ITER, MAXIT
      LOGICAL            :: QPRINT
      REAL ( KIND = DP ) :: DIFF, TOLERANCE

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:NWINDOWS) :: FOLD

      ! . Initialize some algorithmic defaults.
      MAXIT     = 500
      QPRINT    = .FALSE.
      TOLERANCE = 1.0E-3_DP

      ! . Check the convergence parameters.
      IF ( PRESENT ( ITERATIONS      ) ) MAXIT     = ITERATIONS
      IF ( PRESENT ( PRINT           ) ) QPRINT    = PRINT
      IF ( PRESENT ( CONVERGENCE     ) ) TOLERANCE = CONVERGENCE

      ! . Initialize the convergence flag.
      QCONVERGE = .FALSE.

      ! . Initialize the window free energies.
      F = 0.0_DP

      ! . Calculate the umbrella potential factors.
      DO I = 1,NBINS
         UMBPOT(1:NWINDOWS,I) = 0.5_DP * FORCE(1:NWINDOWS) * ( COORDINATE(I) - REFERENCE(1:NWINDOWS) )**2
      END DO
      UMBPOTE = EXP ( - UMBPOT / RT )

      ! . Write out a header.
      IF ( QPRINT ) THEN
	 CALL PRINT_TABLE_OPTIONS ( COLUMNS = 2, HEADER_COLOR = "#AAAA00", PAGEWIDTH = 30, VARIABLEWIDTHS = (/ 10, 20 /) )
	 CALL PRINT_TABLE_START
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "WHAM Equation Solver", COLSPAN = 2, HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Iteration",   HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Convergence", HEADER = .TRUE. )
      END IF

      ! . Start the loop over iterations.
      DO ITER = 1,MAXIT

         ! . Save F.
         FOLD = F

         ! . Calculate the full distribution function.
         DO I = 1,NBINS
            RHO(I) = RAWDATA(I) / SUM ( REAL ( FREQUENCY(1:NWINDOWS), DP ) * &
                     EXP ( - ( UMBPOT(1:NWINDOWS,I) - F(1:NWINDOWS) ) / RT ) )
         END DO

         ! . Redetermine the window free energies (using a simple integration formula).
         F(1:NWINDOWS) = - RT * LOG ( MATMUL ( UMBPOTE(1:NWINDOWS,1:NBINS), RHO(1:NBINS) ) )

         ! . Calculate the maximum difference in the free energies of successive cycles.
         DIFF = MAXVAL ( ABS ( F - FOLD ) )

         ! . Do some printing.
         IF ( QPRINT ) THEN
            WRITE ( PRINT_LINE, "(I20)"    ) ITER ; CALL PRINT_TABLE_ELEMENT
            WRITE ( PRINT_LINE, "(F20.10)" ) DIFF ; CALL PRINT_TABLE_ELEMENT
         END IF

         ! . Check for convergence.
         IF ( DIFF < TOLERANCE ) THEN

            ! . Normalize RHO.
            RHO = RHO / SUM ( RHO )

            ! . Calculate the PMF.
            DO I = 1,NBINS
               IF ( RHO(I) > 0.0_DP ) THEN
                  PMF(I) = - RT * LOG ( RHO(I) )
               ELSE
                  PMF(I) = HUGE ( 0.0_DP )
               END IF
            END DO

            ! . Adjust the PMF so that its lowest value is zero.
            PMF = PMF - MINVAL ( PMF )

            ! . Exit.
            QCONVERGE = .TRUE. ; EXIT

         END IF

      END DO

      ! . Finish writing.
      IF ( QPRINT ) CALL PRINT_TABLE_STOP

      END SUBROUTINE WHAM_SOLVE

   END SUBROUTINE WHAM_ANALYSE

END MODULE WHAM
