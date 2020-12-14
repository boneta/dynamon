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
!==============================================================================
!                  The Conjugate Gradient Minimization Module
!==============================================================================
!
! . Routines:
!
!   CONJUGATE_GRADIENT_MINIMIZE           Minimize a multidimension function.
!
! . Notes:
!
!   This code is a mish-mash of various conjugate gradient procedures. It is
!   not divided into separate subroutines because the SGI compiler has problems
!   when the passed procedure FCALC is called from within an internal
!   subroutine.
!
!==============================================================================
MODULE CONJUGATE_GRADIENT

! . Utility data structure declarations.
USE DEFINITIONS, ONLY : DP
!USE PRINTING,    ONLY : PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
!                        PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT,              &
!			PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CONJUGATE_GRADIENT_MINIMIZE ( FCALC, X, STATUS, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, GRADIENT_TOLERANCE )
   !-----------------------------------------------------------------------------------------------------------------------

   ! . External subroutine declarations.
   INTERFACE
      SUBROUTINE FCALC ( X, F, G )
         USE DEFINITIONS, ONLY : DP
         REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)  :: X
         REAL ( KIND = DP ),               INTENT(OUT) :: F
         REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: G
      END SUBROUTINE FCALC
   END INTERFACE

   ! . Essential arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: X

   ! . Optional arguments.
   INTEGER,            INTENT(IN),  OPTIONAL :: PRINT_FREQUENCY, STEP_NUMBER
   INTEGER,            INTENT(OUT), OPTIONAL :: STATUS
   REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: GRADIENT_TOLERANCE, STEP_SIZE

   ! . Local scalars.
   INTEGER            :: IFAIL, NPRINT, NSTEP, NVAR
   REAL ( KIND = DP ) :: GTOLERANCE, STPSIZ, TOLGRD

   !---------------------------------------------------------------------------
   ! . Minimization Variables.
   !---------------------------------------------------------------------------
   ! . Local parameters.
   INTEGER, PARAMETER :: MSEARCH = 5, NTRIES = 2

   ! . Local scalars.
   INTEGER            :: IRETRY, MBEST, MCALLS, MSTART, NBEST, NRESTART, NSEARCH
   LOGICAL            :: QBETA, QGAMMA, QTERMINATE
   REAL ( KIND = DP ) :: BETA, BETFAC, DDBEST, DDFIRST, DDNEW, DDSPLINE, DFGUESS, F, F_BEST, F_DELTA, F_FIRST, &
                         GAMFAC, GAMMA, GSPLINE, G2, G2_BEST, STBEST, STCURR, STMAX, TEMP

   ! . Local arrays.
! -- this is too much for the heap... :(
!   REAL ( KIND = DP ), DIMENSION(1:SIZE(X)) :: G, G_BEST, G_FIRST, G_RESTART, X_BEST, X_RESTART, X_SEARCH
   REAL ( KIND = DP ), DIMENSION(:), ALLOCATABLE :: G, G_BEST, G_FIRST, G_RESTART, X_BEST, X_RESTART, X_SEARCH

   !---------------------------------------------------------------------------
   ! . Initialization.
   !---------------------------------------------------------------------------
   ! . Initialize the IFAIL counter.
   IFAIL = 0

   ! . Find the number of variables.
   NVAR = SIZE ( X )

allocate( g(1:nvar), g_best(1:nvar), g_first(1:nvar), g_restart(1:nvar), x_best(1:nvar), x_restart(1:nvar), x_search(1:nvar) )

   ! . Skip the calculation if there are no variables.
   IF ( NVAR <= 0 ) RETURN

   ! . Initialize some algorithm options.
   NPRINT = 0
   NSTEP  = 0
   STPSIZ = 1.0E-3_DP
   TOLGRD = 1.0E-3_DP

   ! . Assign the input parameters.
   IF ( PRESENT ( GRADIENT_TOLERANCE ) ) TOLGRD = GRADIENT_TOLERANCE
   IF ( PRESENT ( PRINT_FREQUENCY    ) ) NPRINT = PRINT_FREQUENCY
   IF ( PRESENT ( STEP_NUMBER        ) ) NSTEP  = STEP_NUMBER
   IF ( PRESENT ( STEP_SIZE          ) ) STPSIZ = STEP_SIZE

   ! . Check the input parameters.
   IF ( ( NPRINT < 0 ) .OR. ( NPRINT > NSTEP ) ) NPRINT = 0

   ! . Calculate the gradient tolerance.
   GTOLERANCE = REAL ( NVAR, DP ) * TOLGRD * TOLGRD

   ! . Print out the control parameters.
   IF ( NPRINT > 0 ) THEN

      ! . Print out the options.
!      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#55FF55", VARIABLEWIDTH = 16 )
!      CALL PRINT_SUMMARY_START ( "Conjugate-Gradient Minimization Calculation" )
!      WRITE ( PRINT_LINE, "(I16)"    ) NVAR   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Variables"  )
!      WRITE ( PRINT_LINE, "(I16)"    ) NSTEP  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations" )
!      WRITE ( PRINT_LINE, "(I16)"    ) NPRINT ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency"      )
!      WRITE ( PRINT_LINE, "(G16.10)" ) TOLGRD ; CALL PRINT_SUMMARY_ELEMENT ( "Gradient Tolerance"   )
!      WRITE ( PRINT_LINE, "(G16.10)" ) STPSIZ ; CALL PRINT_SUMMARY_ELEMENT ( "Step Size"            )
!      CALL PRINT_SUMMARY_STOP

      ! . Print out the header for the iterations.
!      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 3, HEADER_COLOR = "#55FF55", PAGEWIDTH = 60, VARIABLEWIDTHS = (/ 10, 20, 20 /) )
!      CALL PRINT_TABLE_START
!      CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",         HEADER = .TRUE. )
!      CALL PRINT_TABLE_ELEMENT ( TEXT = "Function",     HEADER = .TRUE. )
!      CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS Gradient", HEADER = .TRUE. )

write( 6, "(80('-'))" )
write( 6, "(a)" ) "Conjugate-Gradient Minimization Calculation"
write( 6, "(80('-'))" )
write( 6, "(a,i16)" ) "Number of Variables:   ", nvar
write( 6, "(a,i16)" ) "Number of Iterations:  ", nstep
write( 6, "(a,i16)" ) "Print Frequency:       ", nprint
write( 6, "(a,g16.10)" ) "Gradient Tolerance:    ", tolgrd
write( 6, "(a,g16.10)" ) "Step Size:             ", stpsiz
write( 6, "(50('-'))" )

   END IF

   ! . Check the step number.
   IF ( NSTEP <= 0 ) RETURN

   !---------------------------------------------------------------------------
   ! . Set up the minimization.
   !---------------------------------------------------------------------------
   ! . Initialize some counters.
   MBEST    = 1 ! The best point.
   MCALLS   = 1 ! The number of function evaluations.
   NBEST    = 0 ! The most recent iteration that decreases F.
   NRESTART = 0 ! The most recent restart iteration (unless doing steepest descents when it is zero).
   NSEARCH  = 0 ! The number of search directions.

   ! . Evaluate the function and its derivatives at the initial point.
   CALL FCALC ( X, F, G )

   ! . Calculate the square of G.
   G2 = DOT_PRODUCT ( G, G )

   ! . Do some printing.
   IF ( NPRINT == 1 ) THEN
!      WRITE ( PRINT_LINE, "(I10)"   ) MCALLS                          ; CALL PRINT_TABLE_ELEMENT
!      WRITE ( PRINT_LINE, "(F20.8)" ) F                               ; CALL PRINT_TABLE_ELEMENT
!      WRITE ( PRINT_LINE, "(F20.8)" ) SQRT ( G2 / REAL ( NVAR, DP ) ) ; CALL PRINT_TABLE_ELEMENT
write( 6, "(i10,f20.8,f20.8)" ) mcalls, f, sqrt ( g2 / real ( nvar, dp ) )
   END IF

   ! . Save the best point.
   F_BEST  = F
   G_BEST  = G
   G2_BEST = G2
   X_BEST  = X

   ! . Check for convergence.
   IF ( G2 <= GTOLERANCE ) THEN
      GO TO 9999
   ! . Check the number of steps.
   ELSE IF ( NSTEP == 1 ) THEN
      IFAIL = 131
      GO TO 9999
   END IF

   ! . Set the initial search direction (steepest descent) and directional derivative.
   X_SEARCH = -G
   DDNEW    = -G2

   ! . Set the estimated size of the reduction in the function value per iteration.
   DFGUESS = STPSIZ

   ! . Set the best step size.
   STBEST = STPSIZ / G2_BEST

   !---------------------------------------------------------------------------
   ! . Perform the loop over search directions.
   !---------------------------------------------------------------------------
   ! . Start an infinite loop.
   DO

      ! . Increment the loop counter.
      NSEARCH = NSEARCH + 1

      ! . Store the initial function value and gradient.
      F_FIRST = F
      G_FIRST = G

      ! . Calculate the initial directional derivative.
      DDFIRST = DOT_PRODUCT ( G, X_SEARCH )

      ! . Check to see if the search direction is uphill.
      IF ( DDFIRST >= 0.0_DP ) THEN
	 IFAIL = 130 ; EXIT
	 EXIT
      END IF

      ! . Do some more initialization.
      IRETRY = -1
      MSTART = MCALLS
      DDBEST = DDFIRST
      STMAX  = -1.0_DP

      ! . Set the step length and the step length to the best point.
      STCURR = MIN ( STBEST, ABS ( DFGUESS / DDFIRST ) )
      STBEST = 0.0_DP

      !------------------------------------------------------------------------
      ! . Perform the line search loop.
      !------------------------------------------------------------------------
      ! . Start a second infinite loop.
      DO

         ! . Initialize the QBETA flag.
	 QBETA = .TRUE.

	 ! . Determine the new point.
	 X = X_BEST + STCURR * X_SEARCH

	 ! . Proceed for a non-zero change.
	 IF ( MAXVAL ( ABS ( X - X_BEST ) ) > 0.0_DP ) THEN

	    ! . Increment the number of function evaluations.
	    MCALLS = MCALLS + 1

	    ! . Evaluate the function and its derivatives.
	    CALL FCALC ( X, F, G )

            ! . Determine the new directional derivative and the square of G.
	    DDNEW = DOT_PRODUCT ( G, X_SEARCH )
	    G2    = DOT_PRODUCT ( G, G )

	    ! . Do some printing.
	    IF ( NPRINT > 0 ) THEN
	       IF ( MOD ( MCALLS, NPRINT ) == 0 ) THEN
!		  WRITE ( PRINT_LINE, "(I10)"   ) MCALLS                          ; CALL PRINT_TABLE_ELEMENT
!		  WRITE ( PRINT_LINE, "(F20.8)" ) F                               ; CALL PRINT_TABLE_ELEMENT
!		  WRITE ( PRINT_LINE, "(F20.8)" ) SQRT ( G2 / REAL ( NVAR, DP ) ) ; CALL PRINT_TABLE_ELEMENT
write( 6, "(i10,f20.8,f20.8)" ) mcalls, f, sqrt ( g2 / real ( nvar, dp ) )
               END IF
	    END IF

	    ! . Determine the difference in function values.
	    F_DELTA = F - F_BEST

	    ! . The function value has decreased or remained the same.
	    IF ( F_DELTA <= 0.0_DP ) THEN

               ! . Save the best point.
               IF ( ( F_DELTA < 0.0_DP ) .OR. ( DDNEW / DDBEST >= -1.0_DP ) ) THEN
		  MBEST   = MCALLS
		  F_BEST  = F
		  G_BEST  = G
		  G2_BEST = G2
		  X_BEST  = X
               END IF

               ! . Check for convergence.
               IF ( G2 <= GTOLERANCE ) GO TO 9999

	    END IF

	    ! . Check the number of function evaluations.
	    IF ( MCALLS >= NSTEP ) THEN
	       IFAIL = 131 ; GO TO 9999
	    END IF

	    ! . Determine the second derivative of the spline at STBEST using only the gradient.
	    TEMP     = ( F_DELTA + F_DELTA ) / STCURR - DDNEW - DDBEST
	    DDSPLINE = ( DDNEW - DDBEST ) / STCURR

            ! . Reset DDBEST, STBEST and STMAX.
	    IF ( MCALLS > MBEST ) THEN
	       STMAX  = STBEST + STCURR
	    ELSE
               IF ( DDBEST * DDNEW <= 0.0_DP ) STMAX = STBEST
               STBEST = STBEST + STCURR
               DDBEST = DDNEW
               STCURR = -STCURR
	    END IF

            ! . Determine the second derivative of the spline at STBEST for a non-zero energy change.
	    IF ( F_DELTA /= 0.0_DP ) DDSPLINE = DDSPLINE + ( TEMP + TEMP ) / STCURR

            ! . The directional derivative is non-zero.
	    IF ( DDBEST /= 0.0_DP ) THEN

               ! . Force at least two steps in the line search.
	       IF ( MCALLS <= MSTART + 1 ) THEN
        	  QBETA = .FALSE.

               ! . The line search has not yet reached convergence.
	       ELSE IF ( ABS ( DDBEST / DDFIRST ) > 0.2_DP ) THEN

		  ! . Abandon the line search if it has involved too many function evaluations.
		  IF ( MCALLS >= ( MBEST + MSEARCH ) ) THEN
        	     IFAIL = 129 ; GO TO 9999
		  ! . Continue the line search.
		  ELSE
		     QBETA = .FALSE.
		  END IF
               END IF
	    END IF

	 ! . STCURR is effectively zero.
	 ELSE

	    ! . Check to see if the line search is to be abandoned.
	    IF ( ( MCALLS > MSTART + 1 ) .OR. ( ABS ( DDBEST / DDFIRST ) > 0.2_DP ) ) THEN
	       IFAIL = 129 ; GO TO 9999
	    END IF

	 END IF

         ! . Apply the BETA test.
	 IF ( QBETA ) THEN

            ! . Initialize the line search termination flag.
	    QTERMINATE = .FALSE.

	    ! . Restore the values of F, G and X for the best point.
	    IF ( MCALLS /= MBEST ) THEN
	       F = F_BEST
	       G = G_BEST
	       X = X_BEST
	    END IF

	    ! . Calculate the value of BETA for the new search direction.
	    BETFAC = DOT_PRODUCT ( G, G_FIRST )
	    BETA   = ( G2_BEST - BETFAC ) / ( DDBEST - DDFIRST )

            ! . Test for termination of the line search.
	    IF ( ABS ( BETA * DDBEST ) <= 0.2_DP * G2_BEST ) THEN
	       QTERMINATE = .TRUE.

            ! . The BETA test has failed.
	    ELSE

               ! . Increment the retry counter.
	       IRETRY = IRETRY + 1

               ! . Terminate the line search.
	       IF ( IRETRY > 0 ) THEN
		  QTERMINATE = .TRUE.

               ! . Continue the line search for one more step.
	       ELSE

        	  ! . Abandon the line search if it has involved too many function evaluations.
		  IF ( MCALLS >= ( MBEST + MSEARCH ) ) THEN
        	     IFAIL = 129 ; GO TO 9999
		  END IF

	       END IF
	    END IF

            ! . Check to see if the line search is to be terminated.
	    IF ( QTERMINATE ) EXIT

         END IF

         ! . Determine a new step length.
	 STCURR = 0.5_DP * ( STMAX - STBEST )
	 IF ( STMAX < -0.5_DP ) STCURR = 9.0_DP * STBEST
	 GSPLINE = DDBEST + STCURR * DDSPLINE
	 IF ( DDBEST * GSPLINE < 0.0_DP ) STCURR = STCURR * DDBEST / ( DDBEST - GSPLINE )

      END DO

      !------------------------------------------------------------------------
      ! . Determine a new search direction.
      !------------------------------------------------------------------------
      ! . See if the line search has resulted in a reduction in energy.
      IF ( F < F_FIRST ) NBEST = NSEARCH

      ! . Exit if too many attempts have been made to reduce the function value.
      IF ( NSEARCH >= ( NBEST + NTRIES ) ) THEN
	 IFAIL = 132 ; EXIT
      END IF

      ! . Set DFGUESS to the predicted reduction in F for the next iteration.
      DFGUESS = STBEST * DDFIRST

      ! . Start with a new direction.
      IF ( IRETRY > 0 ) THEN

         ! . Use a steepest descent direction.
         NRESTART =  0
         X_SEARCH = -G

      ! . Modify the present search direction.
      ELSE

         ! . Initialize QGAMMA.
	 QGAMMA = .FALSE.

         ! . Determine whether GAMMA is to be used in the new search direction.
	 IF ( ( NRESTART /= 0 ) .AND. ( ( NSEARCH - NRESTART ) < NVAR ) .AND. ( ABS ( BETFAC ) < 0.2_DP * G2_BEST ) ) THEN

	    ! . Calculate GAMMA.
	    GAMMA = DOT_PRODUCT ( G, G_RESTART ) / GAMFAC

	    ! . Use GAMMA if the new search direction is sufficiently downhill.
	    IF ( ABS ( BETA * DDBEST + GAMMA * DOT_PRODUCT ( G, X_RESTART ) ) < 0.2_DP * G2_BEST ) THEN
	       QGAMMA   = .TRUE.
               X_SEARCH = -G + BETA * X_SEARCH + GAMMA * X_RESTART
	    END IF

	 END IF

	 ! . GAMMA was not used.
	 IF ( .NOT. QGAMMA ) THEN

            ! . Reset the restart variables.
	    NRESTART  = NSEARCH
	    GAMFAC    = DDBEST - DDFIRST
	    G_RESTART = G - G_FIRST
	    X_RESTART = X_SEARCH

            ! . Determine the new search direction.
	    X_SEARCH  = -G + BETA * X_SEARCH

	 END IF

      END IF

   ! . End of the line search loop.
   END DO

   ! . End of the minimization.
   9999 CONTINUE

   ! . Restore the best point if the minimization has failed.
   IF ( IFAIL /= 0 ) X = X_BEST

   !---------------------------------------------------------------------------
   ! . Finish up.
   !---------------------------------------------------------------------------
   ! . Do some more printing.
   IF ( NPRINT > 0 ) THEN

      ! . Print out the terminator.
!      CALL PRINT_TABLE_STOP
write( 6, "(50('-'))" )

      ! . Print out the status.
      SELECT CASE ( IFAIL )
!      CASE (   0 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Gradient tolerance reached."      )
!      CASE ( 129 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Line search abandoned."           )
!      CASE ( 130 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Search direction uphill."         )
!      CASE ( 131 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Too many steps."                  )
!      CASE ( 132 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Unable to reduce function value." )
      CASE (   0 ) ; write( 6, "(/a/)" ) "Minimization Status: Gradient tolerance reached."
      CASE ( 129 ) ; write( 6, "(/a/)" ) "Minimization Status: Line search abandoned."
      CASE ( 130 ) ; write( 6, "(/a/)" ) "Minimization Status: Search direction uphill."
      CASE ( 131 ) ; write( 6, "(/a/)" ) "Minimization Status: Too many steps."
      CASE ( 132 ) ; write( 6, "(/a/)" ) "Minimization Status: Unable to reduce function value."
      END SELECT

   END IF

   ! . Assign a value to STATUS if necessary.
   IF ( PRESENT ( STATUS ) ) STATUS = IFAIL

deallocate( g, g_best, g_first, g_restart, x_best, x_restart, x_search )

   END SUBROUTINE CONJUGATE_GRADIENT_MINIMIZE

END MODULE CONJUGATE_GRADIENT

