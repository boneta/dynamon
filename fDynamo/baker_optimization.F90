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
!                         The Baker Optimization Module
!===============================================================================
!
! . Subroutines:
!
!   BAKER_SEARCH                    Search for a stationary point.
!
!===============================================================================
MODULE BAKER_OPTIMIZATION

! . Module declarations.
USE DEFINITIONS,     ONLY : DP
USE DIAGONALIZATION, ONLY : SYMMETRIC_UPPER
!USE PRINTING,        ONLY : PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
!                            PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT,              &
!			    PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP

IMPLICIT NONE
PRIVATE
PUBLIC :: BAKER_SEARCH

!==============================================================================
CONTAINS
!==============================================================================

   !--------------------------------------------------------------------------------------------------------
   SUBROUTINE BAKER_SEARCH ( FCALC, X, STATUS, FOLLOW_MODE, FOLLOW_VARIABLE, PRINT_FREQUENCY, STEP_NUMBER, &
                                       LOCATE_SADDLE, USE_NR_STEP, GRADIENT_TOLERANCE, MAXIMUM_EIGENVALUE, &
                                       MAXIMUM_STEP, MINIMUM_EIGENVALUE )
   !--------------------------------------------------------------------------------------------------------

   ! . External subroutine declarations.
   INTERFACE
      SUBROUTINE FCALC ( X, F, G, H )
         USE DEFINITIONS, ONLY : DP
         REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)  :: X
         REAL ( KIND = DP ),               INTENT(OUT) :: F
         REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: G
         REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: H
      END SUBROUTINE FCALC
   END INTERFACE

   ! . Essential arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: X

   ! . Optional arguments.
   INTEGER,            INTENT(IN),  OPTIONAL :: FOLLOW_MODE, FOLLOW_VARIABLE, PRINT_FREQUENCY, STEP_NUMBER
   INTEGER,            INTENT(OUT), OPTIONAL :: STATUS
   LOGICAL,            INTENT(IN),  OPTIONAL :: LOCATE_SADDLE, USE_NR_STEP
   REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: GRADIENT_TOLERANCE, MAXIMUM_EIGENVALUE, MAXIMUM_STEP, MINIMUM_EIGENVALUE

   ! . Local scalars.
   INTEGER            :: CURMOD, CURVAR, IFAIL, ISTEP, NEGREQ, NPRINT, NSTEP, NUMNEG, NVAR
   LOGICAL            :: QNR
   REAL ( KIND = DP ) :: EIGMAX, EIGMIN, F, GRMS, MAXSTP, TOLGRD, TMP

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: DX, EIGVAL, G, GX, H, MODVEC
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: EIGVEC

   ! . Initialize the IFAIL counter.
   IFAIL = 0

   ! . Find the number of variables.
   NVAR = SIZE ( X )

   ! . Skip the calculation if there are no variables.
   IF ( NVAR <= 0 ) RETURN

   ! . Initialize some algorithm options.
   CURMOD = 0
   CURVAR = 0
   NEGREQ = 0
   NPRINT = 0
   NSTEP  = 0
   QNR    = .FALSE.
   EIGMAX = 1.0E+5_DP
   EIGMIN = 1.0E-3_DP
   MAXSTP = 0.3_DP
   TOLGRD = 1.0E-3_DP

   ! . Find the type of optimization.
   IF ( PRESENT ( LOCATE_SADDLE ) ) THEN
      IF ( LOCATE_SADDLE ) NEGREQ = 1
   END IF

   ! . Assign the remaining input parameters.
   IF ( PRESENT ( FOLLOW_MODE        ) ) CURMOD = FOLLOW_MODE
   IF ( PRESENT ( FOLLOW_VARIABLE    ) ) CURVAR = FOLLOW_VARIABLE
   IF ( PRESENT ( GRADIENT_TOLERANCE ) ) TOLGRD = GRADIENT_TOLERANCE
   IF ( PRESENT ( MAXIMUM_EIGENVALUE ) ) EIGMAX = MAXIMUM_EIGENVALUE
   IF ( PRESENT ( MAXIMUM_STEP       ) ) MAXSTP = MAXIMUM_STEP
   IF ( PRESENT ( MINIMUM_EIGENVALUE ) ) EIGMIN = MINIMUM_EIGENVALUE
   IF ( PRESENT ( PRINT_FREQUENCY    ) ) NPRINT = PRINT_FREQUENCY
   IF ( PRESENT ( STEP_NUMBER        ) ) NSTEP  = STEP_NUMBER
   IF ( PRESENT ( USE_NR_STEP        ) ) QNR    = USE_NR_STEP

   ! . Check the input parameters for a minimum search.
   IF ( NEGREQ == 0 ) THEN

      ! . Reset CURMOD and CURVAR.
      CURMOD = 0
      CURVAR = 0

   ! . Check the input parameters for a saddle point search.
   ELSE

      ! . Check CURVAR.
      IF ( CURVAR > NVAR ) CURVAR = 0

      ! . If a variable is being followed set CURMOD to a value greater than 1.
      IF ( CURVAR > 0 ) CURMOD = 2

      ! . Check CURMOD.
      IF ( ( CURMOD < 1 ) .OR. ( CURMOD > NVAR ) ) CURMOD = 1

   END IF

   ! . Check the print frequency.
   IF ( ( NPRINT < 0 ) .OR. ( NPRINT > NSTEP ) ) NPRINT = 0

   ! . Print out the control parameters.
   IF ( NPRINT > 0 ) THEN

      ! . Print out the options.
 !     CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#55FF55", VARIABLEWIDTH = 16 )
write( 6, "(80('-'))")
      IF ( NEGREQ > 0 ) THEN
!         CALL PRINT_SUMMARY_START ( "Baker Search for Saddle Point" )
write( 6, "(a)" ) "Baker Search for Saddle Point"
      ELSE
!         CALL PRINT_SUMMARY_START ( "Baker Search for Minimum" )
write( 6, "(a)" ) "Baker Search for Minimum"
      END IF
write( 6, "(80('-'))")
!      WRITE ( PRINT_LINE, "(I16)"    ) NVAR   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Variables"  )
!      WRITE ( PRINT_LINE, "(I16)"    ) NSTEP  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations" )
!      WRITE ( PRINT_LINE, "(I16)"    ) NPRINT ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency"      )
!      WRITE ( PRINT_LINE, "(L16)"    ) QNR    ; CALL PRINT_SUMMARY_ELEMENT ( "Newton-Raphson Step"  )
!      WRITE ( PRINT_LINE, "(I16)"    ) CURMOD ; CALL PRINT_SUMMARY_ELEMENT ( "Follow Mode"          )
!      WRITE ( PRINT_LINE, "(I16)"    ) CURVAR ; CALL PRINT_SUMMARY_ELEMENT ( "Follow Variable"      )
!      WRITE ( PRINT_LINE, "(G16.10)" ) MAXSTP ; CALL PRINT_SUMMARY_ELEMENT ( "Maximum Step"         )
!      WRITE ( PRINT_LINE, "(G16.10)" ) TOLGRD ; CALL PRINT_SUMMARY_ELEMENT ( "Gradient Tolerance"   )
!      WRITE ( PRINT_LINE, "(G16.10)" ) EIGMAX ; CALL PRINT_SUMMARY_ELEMENT ( "Maximum Eigenvalue"   )
!      WRITE ( PRINT_LINE, "(G16.10)" ) EIGMIN ; CALL PRINT_SUMMARY_ELEMENT ( "Minimum Eigenvalue"   )
!      CALL PRINT_SUMMARY_STOP
write( 6, "(a,i16)" ) "Number of Variables:  ", nvar
write( 6, "(a,i16)" ) "Number of Iterations: ", nstep
write( 6, "(a,i16)" ) "Print Frequency:      ", nprint
write( 6, "(a,L16)" ) "Newton-Raphson Step:  ", qnr
write( 6, "(a,i16)" ) "Follow Mode:          ", curmod
write( 6, "(a,i16)" ) "Follow Variable:      ", curvar
write( 6, "(a,g16.10)" ) "Maximum Step:         ", maxstp
write( 6, "(a,g16.10)" ) "Gradient Tolerance:   ", tolgrd
write( 6, "(a,g16.10)" ) "Maximum Eigenvalue:   ", eigmax
write( 6, "(a,g16.10)" ) "Minimum Eigenvalue:   ", eigmin
write( 6, "(80('-'))")

      ! . Print out the header for the iterations.
!      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#55FF55", PAGEWIDTH = 80, VARIABLEWIDTHS = (/ 10, 20, 20, 10, 20 /) )
!      CALL PRINT_TABLE_START
!      CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",         HEADER = .TRUE. )
!      CALL PRINT_TABLE_ELEMENT ( TEXT = "Function",     HEADER = .TRUE. )
!      CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS Gradient", HEADER = .TRUE. )
!      CALL PRINT_TABLE_ELEMENT ( TEXT = "Order",        HEADER = .TRUE. )
!      CALL PRINT_TABLE_ELEMENT ( TEXT = "Eigenvalue",   HEADER = .TRUE. )
write( 6, "(a10,a20,a20,a10,a20)" ) "Step", "Function", "RMS Gradient", "Order", "Eigenvalue"
write( 6, "(80('-'))")

   END IF

   !---------------------------------------------------------------------------
   ! . Perform the optimization.
   !---------------------------------------------------------------------------
   ! . Allocate some temporary arrays.
   ALLOCATE ( DX(1:NVAR), EIGVAL(1:NVAR), EIGVEC(1:NVAR,1:NVAR), G(1:NVAR), GX(1:NVAR), &
              H(1:(NVAR*(NVAR+1))/2), MODVEC(1:NVAR) )

   ! . Initialize the step vector.
   DX = 0.0_DP

   ! . Loop over the steps.
   DO ISTEP = 1,NSTEP

      ! . Form the new variable vector.
      X = X + DX

      ! . Calculate the function and its derivatives.
      CALL FCALC ( X, F, G, H )

      ! . Diagonalize the Hessian matrix.
      CALL SYMMETRIC_UPPER ( H, EIGVAL, EIGVEC )

      ! . Check the magnitudes of the Hessian eigenvalues.
      WHERE ( ABS ( EIGVAL ) < EIGMIN ) EIGVAL = SIGN ( EIGMIN, EIGVAL )
      WHERE ( ABS ( EIGVAL ) > EIGMAX ) EIGVAL = SIGN ( EIGMAX, EIGVAL )

      ! . Determine the number of negative eigenvalues.
      NUMNEG = COUNT ( EIGVAL < 0.0_DP )
      TMP = EIGVAL(1)

      ! . Transform the gradient vector to the local Hessian modes.
      GX = MATMUL ( G, EIGVEC )

      ! . Calculate the step vector.
!write(87,"(a)") "coor"
!write(87,"(f20.10)") x
!write(87,"(a)") "grad"
!write(87,"(f20.10)") g
!write(87,"(a)") "values"
!write(87,"(f20.10)") eigval
!write(87,"(a)") "grad2hess"
!write(87,"(f20.10)") gx
      CALL GET_STEP
!write(87,"(a)") "curmod"
!write(87,"(f20.10)") modvec
!write(87,"(a)") "step"
!write(87,"(f20.10)") dx
      IF ( IFAIL /= 0 ) GO TO 10

      ! . Find the RMS gradient.
      GRMS = SQRT ( DOT_PRODUCT ( G, G ) / REAL ( NVAR, DP ) )

      ! . Do some printing.
      IF ( NPRINT > 0 ) THEN
         IF ( MOD ( ISTEP, NPRINT ) == 0 ) THEN
!            WRITE ( PRINT_LINE, "(I10)"   ) ISTEP  ; CALL PRINT_TABLE_ELEMENT
!            WRITE ( PRINT_LINE, "(F20.8)" ) F      ; CALL PRINT_TABLE_ELEMENT
!            WRITE ( PRINT_LINE, "(F20.8)" ) GRMS   ; CALL PRINT_TABLE_ELEMENT
!            WRITE ( PRINT_LINE, "(I10)"   ) NUMNEG ; CALL PRINT_TABLE_ELEMENT
!            WRITE ( PRINT_LINE, "(F20.8)" ) TMP    ; CALL PRINT_TABLE_ELEMENT
write( 6, "(i10,f20.8,f20.9,2i5,f20.8)" ) istep, f, grms, numneg, curmod, tmp
         END IF
      END IF

      ! . Check for convergence.
      IF ( ( GRMS < TOLGRD ) .AND. ( NUMNEG == NEGREQ ) ) GO TO 10

   END DO

   ! . There have been too many iterations.
   IFAIL = 132

   ! . Error or normal termination.
   10 CONTINUE

   ! . Do some more printing.
   IF ( NPRINT > 0 ) THEN

      ! . Print out the terminator.
!      CALL PRINT_TABLE_STOP
write( 6, "(80('-'))")

      ! . Print out the status.
      SELECT CASE ( IFAIL )
!      CASE (   0 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Baker Search Status: Gradient tolerance reached."       )
!      CASE ( 129 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Baker Search Status: Too many LAMBDA iterations."       )
!      CASE ( 130 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Baker Search Status: Error in determination of LAMBDA." )
!      CASE ( 131 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Baker Search Status: Step size too small."              )
!      CASE ( 132 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Baker Search Status: Too many steps."                   )
      CASE (   0 ) ; write( 6, "(/a/)" ) "Baker Search Status: Gradient tolerance reached."
      CASE ( 129 ) ; write( 6, "(/a/)" ) "Baker Search Status: Too many LAMBDA iterations."
      CASE ( 130 ) ; write( 6, "(/a/)" ) "Baker Search Status: Error in determination of LAMBDA."
      CASE ( 131 ) ; write( 6, "(/a/)" ) "Baker Search Status: Step size too small."
      CASE ( 132 ) ; write( 6, "(/a/)" ) "Baker Search Status: Too many steps."
      END SELECT

   END IF

   ! . Assign a value to STATUS if necessary.
   IF ( PRESENT ( STATUS ) ) STATUS = IFAIL

   ! . Deallocate all temporary space.
   DEALLOCATE ( DX, EIGVAL, EIGVEC, G, GX, H, MODVEC )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !------------------
      SUBROUTINE GET_STEP
      !------------------

      ! . Local parameters.
      INTEGER,            PARAMETER :: MAXIT = 999
      REAL ( KIND = DP ), PARAMETER :: LARGE = 1.0E+6_DP, MINSTP = 1.0E-1_DP, STEP = 50.0_DP, TOL1 = 1.0E-4_DP, TOL2 = 1.0E-8_DP

      ! . Local scalars.
      INTEGER            :: I, IT, LOWER
      REAL ( KIND = DP ) :: BIGELE, DXSIZE, LAMBDA, LAMBDA0, LAMBDA1, LAMBDA2, OVERLAP, SCALE, TEMP

      ! . Initialize the step vector.
      DX = 0.0_DP

      ! . Do a N-R step if the number of negative eigenvalues is OK and QNEWTON is on.
      IF ( QNR .AND. ( NUMNEG == NEGREQ ) ) THEN

         ! . Calculate the Newton-Raphson step.
         DX = DX + MATMUL ( EIGVEC, ( GX / EIGVAL ) )

      ! . Take a P-RFO step for a TS search and an RFO step for a minimum search.
      ELSE

         ! . A saddle point is being searched for.
         IF ( NEGREQ > 0 ) THEN

            ! . A specific mode is active.
            IF ( CURMOD > 1 ) THEN

               ! . On the first step determine which mode to follow.
               IF ( ISTEP == 1 ) THEN

                  ! . Find the mode with the largest magnitude for a particular variable.
                  IF ( CURVAR > 0 ) THEN

                     ! . Find the eigenvector with the largest magnitude in the given direction.
                     CURMOD = 1
                     BIGELE = ABS ( EIGVEC(CURVAR,1) )
                     DO I = 2,NVAR
                        IF ( ABS ( EIGVEC(CURVAR,I) ) > BIGELE ) THEN
                           BIGELE = ABS ( EIGVEC(CURVAR,I) )
                           CURMOD = I
                        END IF
                     END DO

                  END IF

               ! . Determine which eigenvector has the greatest overlap with the existing mode being followed.
               ELSE

                  CURMOD = 1
                  BIGELE = ABS ( DOT_PRODUCT ( EIGVEC(1:NVAR,1), MODVEC ) )
                  DO I = 2,NVAR
                     OVERLAP = ABS ( DOT_PRODUCT ( EIGVEC(1:NVAR,I), MODVEC ) )
                     IF ( OVERLAP > BIGELE ) THEN
                        BIGELE = OVERLAP
                        CURMOD = I
                     END IF
                  END DO

               END IF

               ! . Save the eigenvector for the mode being followed.
               MODVEC = EIGVEC(1:NVAR,CURMOD)

            END IF

            ! . The mode scaling factor is normal.
            IF ( ABS ( GX(CURMOD) ) > TOL1 ) THEN

               ! . Calculate LAMBDA0.
               LAMBDA0 = 0.5_DP * ( EIGVAL(CURMOD) + SQRT ( EIGVAL(CURMOD)**2 + 4.0_DP * GX(CURMOD)**2 ) )

               ! . Calculate the mode scaling factor.
               SCALE = GX(CURMOD) / ( LAMBDA0 - EIGVAL(CURMOD) )

            ! . The mode scaling factor is not normal.
            ELSE

               ! . Use a Newton-Raphson step if the curvature is OK.
               IF ( NUMNEG == NEGREQ ) THEN
                  SCALE = - GX(CURMOD) / EIGVAL(CURMOD)
               ! . Use a minimum step length.
               ELSE
                  SCALE = MINSTP
               END IF

            END IF

            ! . Calculate the step for the maximization.
            DX = SCALE * EIGVEC(1:NVAR,CURMOD)

         END IF


         ! . Find the step for the minimization.
         IF ( ( NVAR - NEGREQ ) > 0 ) THEN

            ! . Find the starting mode.
            IF ( CURMOD == 1 ) THEN
               LOWER = 2
            ELSE
               LOWER = 1
            END IF

            ! . Initialize LAMBDA.
            LAMBDA = 0.0_DP
            IF ( EIGVAL(LOWER) < 0.0_DP ) THEN
               LAMBDA  = EIGVAL(LOWER) - STEP
               LAMBDA1 = EIGVAL(LOWER)
               LAMBDA2 = - LARGE
            END IF

            ! . Start of the iterations.
            DO IT = 1,MAXIT

               TEMP = 0.0_DP
               DO I = 1,NVAR
                  IF ( I == CURMOD ) CYCLE
                  TEMP = TEMP + ( GX(I) * GX(I) ) / ( LAMBDA - EIGVAL(I) )
               END DO

               ! . Check for LAMBDA convergence.
               IF ( ABS ( LAMBDA - TEMP ) < TOL2 ) GO TO 10

               ! . Straightforward iterations.
               IF ( EIGVAL(LOWER) > 0.0_DP ) THEN
                  LAMBDA = TEMP
               ! . Reduce the range for the search.
               ELSE
                  IF ( TEMP < LAMBDA ) LAMBDA1 = LAMBDA
                  IF ( TEMP > LAMBDA ) LAMBDA2 = LAMBDA
                  IF ( LAMBDA2 > - LARGE ) THEN
                     LAMBDA = 0.5_DP * ( LAMBDA1 + LAMBDA2 )
                  ELSE IF ( LAMBDA2 == - LARGE ) THEN
                     LAMBDA = LAMBDA - STEP
                  END IF
               END IF

            END DO

            ! . There have been too many iterations.
            IFAIL = 129

            ! . Normal termination.
            10 CONTINUE

            ! . Check the value of LAMBDA.
            IF ( ( LAMBDA > EIGVAL(LOWER) ) .OR. ( EIGVAL(LOWER) > 0.0_DP .AND. LAMBDA > 0.0_DP ) ) THEN
               IFAIL = 130
            END IF

!write(87,"(a)") "lambda"
!write(87,"(f20.10)") lambda
            ! . Modify the CURMOD eigenvalues and eigenvectors.
	    IF ( CURMOD > 0 ) THEN
	       EIGVAL(CURMOD)        = LAMBDA - 1.0_DP
	       EIGVEC(1:NVAR,CURMOD) = 0.0_DP
	    END IF

            ! . Calculate the step for the minimization.
	    DX = DX + MATMUL ( EIGVEC, ( GX / ( LAMBDA - EIGVAL ) ) )

         END IF
      END IF

      ! . Calculate the step size.
      DXSIZE = SQRT ( DOT_PRODUCT ( DX, DX ) )

      ! . The step size is very small.
      IF ( DXSIZE < TOL2 ) THEN
         IFAIL = 131
      ! . Reduce the step size.
      ELSE IF ( DXSIZE > MAXSTP ) THEN
         DX = ( MAXSTP / DXSIZE ) * DX
      END IF


      END SUBROUTINE GET_STEP

   END SUBROUTINE BAKER_SEARCH

END MODULE BAKER_OPTIMIZATION
