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
!                          The BFGS Optimization Module
!===============================================================================
!
! . Subroutines:
!
!   OPTIMIZE_BFGS_CALCULATE                Search for a minimum.
!
! . Notes:
!
!   This subroutine follows closely the CADPAC version of the BFGS algorithm.
!
!===============================================================================
MODULE OPTIMIZE_BFGS

! . Module declarations.
USE CONSTANTS,             ONLY : AU_TO_KJ, BOHRS_TO_ANGSTROMS
USE DEFINITIONS,           ONLY : DP
USE DIAGONALIZATION,       ONLY : SYMMETRIC_UPPER
USE IO_UNITS,              ONLY : OUTPUT

USE ATOMS,                 ONLY : ATMCRD, ATMFIX, NATOMS, NFREE
USE NORMAL_MODE_UTILITIES, ONLY : PROJECT_ROTATION_TRANSLATION
USE POTENTIAL_ENERGY,      ONLY : ATMDER, ETOTAL, GRADIENT

IMPLICIT NONE
PRIVATE
PUBLIC :: OPTIMIZE_BFGS_CALCULATE, &
          MAXSEARCH, NPRINT, STEPMAX, TOLGRD
SAVE

! . Parameters.
REAL ( KIND = DP ), PARAMETER :: ALPHA0 = 1.0_DP, ALPHL = 1.0E-3_DP, ALPHM = 1.1_DP, BROY = 0.0_DP, &
                                 ESMALL = 1.0E-3_DP * AU_TO_KJ, &
                                 STEPL  = 5.0E-3_DP * BOHRS_TO_ANGSTROMS

! . Options.
INTEGER            :: MAXSEARCH = 50, NPRINT = 1
REAL ( KIND = DP ) :: EIGMIN    = 0.005_DP  * AU_TO_KJ / BOHRS_TO_ANGSTROMS**2, &
                      STEPMAX   = 0.2_DP    *            BOHRS_TO_ANGSTROMS,    &
                      THRES10   = 0.0005_DP * AU_TO_KJ / BOHRS_TO_ANGSTROMS**2, &
                      THRES20   = 5.0E-5_DP * AU_TO_KJ / BOHRS_TO_ANGSTROMS**2, &
                      TOLGRD    = 0.5_DP

!===============================================================================
CONTAINS
!===============================================================================

   !---------------------------------------------
   SUBROUTINE OPTIMIZE_BFGS_CALCULATE ( HESSIAN )
   !---------------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: HESSIAN

   ! . Scalars.
   INTEGER            :: ECALLS, IFAIL, NRESET, NSEARCH, NVAR, NVARTR
   LOGICAL            :: QRESET
   REAL ( KIND = DP ) :: ALPHA, ELAST, ENEW, GRMS, STEPMA

   ! . Arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: GDEL, GLAST, GNEW, XDEL, XLAST, XLAST2, XNEW
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: HINV, HINV0

   ! . Initialize some scalars.
   ALPHA   = ALPHA0
   ECALLS  = 0
   IFAIL   = 0
   NRESET  = 0
   NSEARCH = 0
   QRESET  = .FALSE.

   ! . Check the print parameter.
   IF ( ( NPRINT < 0 ) .OR. ( NPRINT > MAXSEARCH ) ) NPRINT = 1

   ! . Determine the number of variables.
   NVAR = 3 * NFREE
   IF ( NVAR <= 0 ) RETURN
   NVARTR = ( NVAR * ( NVAR + 1 ) ) / 2

   ! . Allocate space.
   ALLOCATE ( GDEL(1:NVAR), GLAST(1:NVAR), GNEW(1:NVAR), HINV(1:NVAR,1:NVAR), HINV0(1:NVAR,1:NVAR), XDEL(1:NVAR), XLAST(1:NVAR), &
              XLAST2(1:NVAR), XNEW (1:NVAR))

   ! . Fill XNEW.
   CALL VARIABLES_FILL ( XNEW, ATMCRD )

   ! . Get the starting function value.
   CALL EGCALC ( XNEW, ENEW, GNEW )

   ! . Determine GRMS.
   GRMS = SQRT ( DOT_PRODUCT ( GNEW, GNEW ) / NVAR )

   ! . Save the variables, etc.
   ELAST  = ENEW
   GLAST  = GNEW
   XLAST  = XNEW
   XLAST2 = HUGE ( 0.0_DP )

   ! . Setup the inverse Hessian.
   CALL SETUP_INVERSE_HESSIAN

   ! . Set HINV.
   HINV = HINV0

   ! . Print out some data.
   IF ( NPRINT > 0 ) THEN

      ! . Print out the header.
      WRITE ( OUTPUT, "(1X)" )
      WRITE ( OUTPUT, "(25('-'),A,24('-'))" ) " BFGS Minimization Calculation "

      ! . Print out the options.
      WRITE ( OUTPUT, "(A,I16,2X,A,I16)"    ) "Number of Variables  = ", NVAR,     &
                                              "Number of Searches   = ", MAXSEARCH
      WRITE ( OUTPUT, "(A,I16,2X,A,G16.10)" ) "Print Frequency      = ", NPRINT,   &
                                              "Gradient Tolerance   = ", TOLGRD
      WRITE ( OUTPUT, "(A,G16.10)" )          "Maximum Step Size    = ", STEPMAX

      ! . Print out the terminator.
      WRITE ( OUTPUT, "(80('-'))" )

      ! . Print out the header for the iterations.
      WRITE ( OUTPUT, "(1X)" )
      WRITE ( OUTPUT, "(60('-'))" )
      WRITE ( OUTPUT, "(A,2X,A,11X,A,10X,A)" ) "Search", "Calls", "Function", "RMS Gradient"
      WRITE ( OUTPUT, "(60('-'))" )

      ! . Write out the starting data.
      WRITE ( OUTPUT, "(2I10,2F20.8)" ) 0, ECALLS, ELAST, GRMS

   END IF

   ! . Start the loop over searches.
   DO

      ! . Check for convergence.
      IF ( GRMS < TOLGRD ) EXIT

      ! . Check the number of steps.
      IF ( NSEARCH > MAXSEARCH ) THEN
         IFAIL = 131
         EXIT
      END IF

      ! . Get the new search direction.
      XDEL = - MATMUL ( HINV, GNEW )

      ! . Perform a linear search.
      CALL LINE_SEARCH
      IF ( IFAIL /= 0 ) EXIT

      ! . Reset the inverse Hessian back to its starting value if necessary.
      IF ( QRESET ) THEN
         HINV = HINV0
         CYCLE
      END IF

      ! . Update the inverse Hessian.
      CALL UPDATE_BFGS

   ! . End of the loop.
   END DO

   ! . Finish up.
   ! . Save the best coordinates.
   CALL VARIABLES_EMPTY ( XLAST, ATMCRD )

   ! . Deallocate space.
   DEALLOCATE ( GDEL, GLAST, GNEW, HINV, HINV0, XDEL, XLAST, XLAST2, XNEW )

   ! . Do some more printing.
   IF ( NPRINT > 0 ) THEN

      ! . Print out the terminator.
      WRITE ( OUTPUT, "(60('-'))" )

      ! . Print out the status.
      WRITE ( OUTPUT, "(1X)" )
      SELECT CASE ( IFAIL )
      CASE (   0 ) ; WRITE ( OUTPUT, "(A)" ) "Minimization Status: Gradient tolerance reached."
      CASE ( 129 ) ; WRITE ( OUTPUT, "(A)" ) "Minimization Status: Too many successive increases in energy."
      CASE ( 130 ) ; WRITE ( OUTPUT, "(A)" ) "Minimization Status: Too many Hessian resets."
      CASE ( 131 ) ; WRITE ( OUTPUT, "(A)" ) "Minimization Status: Maximum number of iterations reached."
      END SELECT

   END IF

   !============================================================================
   CONTAINS
   !============================================================================

      !---------------------
      SUBROUTINE LINE_SEARCH
      !---------------------

      ! . Scalars.
      INTEGER            :: NRISES
      REAL ( KIND = DP ) :: ALPH, ALPH0, GS0, SLAST, SMAX

      ! . The inverse Hessian is not to be reset.
      IF ( .NOT. QRESET ) THEN

         ! . Set the value of ALPHA.
         IF ( NSEARCH > 0 )THEN
            ALPHA = ALPHM
         ELSE
            ALPHA = ALPHA0
         END IF

         ! . Determine some constants.
         SMAX = ALPHA * SQRT ( DOT_PRODUCT ( XDEL, XDEL ) )
         IF ( NSEARCH > 0 ) THEN
            SLAST = 4.0_DP * SQRT ( SUM ( ( XLAST - XLAST2 )**2 ) )
            SLAST = MAX ( SLAST, STEPL  )
            SLAST = MIN ( SLAST, STEPMA )
         ELSE
            SLAST = STEPMA
         END IF
         IF ( SMAX  > SLAST ) ALPHA = ALPHA * ( SLAST / SMAX )
         IF ( ALPHA < ALPHL ) ALPHA = ALPHL

      END IF

      ! . Initialize some other variables.
      NRISES = 0
      QRESET = .FALSE.
      ALPH   = ALPHA

      ! . Search along XDEL.
   10 CONTINUE
      XNEW = XLAST + ALPH * XDEL

      ! . Evaluate the function at the new point.
      CALL EGCALC ( XNEW, ENEW, GNEW )

      ! . Check to see if the energy goes up.
      IF ( ( ENEW - ELAST ) > ESMALL ) THEN

         ! . Determine ALPH and ALPH0.
         GS0   = - DOT_PRODUCT ( GLAST, XDEL )
         ALPH0 = ALPH
         ALPH  = ( GS0 * ALPH0**2 ) / ( 2.0_DP * ( GS0 * ALPH0 + ENEW - ELAST ) )
         IF ( ALPH <= ALPH0 ) THEN
            IF ( ALPH > ALPHL ) THEN 
               NRISES = NRISES + 1
               IF ( NRISES > 4 ) THEN
                  IFAIL = 129
                  RETURN
               ELSE
                  GO TO 10
               END IF
            END IF
         END IF

         ! . The inverse Hessian needs to be reset.
         QRESET = .TRUE.
         NRESET = NRESET + 1
         IF ( NRESET > 3 ) IFAIL = 130
         RETURN

      END IF

      ! . The search has been successful.
      NSEARCH = NSEARCH + 1
      GRMS    = SQRT ( DOT_PRODUCT ( GNEW, GNEW ) / NVAR )

      ! . Write out some information.
      IF ( NPRINT > 0 ) THEN
         IF ( MOD ( NSEARCH, NPRINT ) == 0 ) WRITE ( OUTPUT, "(2I10,2F20.8)" ) NSEARCH, ECALLS, ELAST, GRMS
      END IF

      ! . Calculate XDEL and GDEL.
      XDEL = ALPH * XDEL
      GDEL = GNEW - GLAST

      ! . Save the data about the last point.
      XLAST2 = XLAST
      ELAST  = ENEW
      GLAST  = GNEW
      XLAST  = XNEW

      END SUBROUTINE LINE_SEARCH

      !----------------------------
      SUBROUTINE EGCALC ( X, E, G )
      !----------------------------

      ! . Argument declarations.
      REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)  :: X
      REAL ( KIND = DP ),               INTENT(OUT) :: E
      REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: G

      ! . Update the number of energy calls.
      ECALLS = ECALLS + 1

      ! . Move X to the coordinates.
      CALL VARIABLES_EMPTY ( X, ATMCRD )

      ! . Calculate the energy and its first derivatives.
      CALL GRADIENT ( PRINT = .FALSE. ) ; E = ETOTAL ; CALL VARIABLES_FILL ( G, ATMDER )

      END SUBROUTINE EGCALC

      !-------------------------------
      SUBROUTINE SETUP_INVERSE_HESSIAN
      !-------------------------------

      ! . Scalars.
      INTEGER            :: I, J, NZ, NZERO
      REAL ( KIND = DP ) :: SHIFT, THRES1, THRES2

      ! . Arrays.
      REAL ( KIND = DP ), DIMENSION(1:NVAR)        :: EIGENVALUES
      REAL ( KIND = DP ), DIMENSION(1:NVARTR)      :: H
      REAL ( KIND = DP ), DIMENSION(1:NVAR,1:NVAR) :: EIGENVECTORS

      ! . Set some constants.
      IF ( EIGMIN == 0.0_DP ) THEN
         THRES1 = THRES10
      ELSE
         THRES1 = EIGMIN
      END IF
      THRES2 = THRES20

      ! . Get the starting Hessian and the starting step length.
      ! . A starting Hessian was passed in.
      IF ( PRESENT ( HESSIAN ) ) THEN
          H      = HESSIAN
          STEPMA = 1.5_DP * STEPMAX
      ! . Set the starting Hessian to the identity matrix.
      ELSE
          H = 0.0_DP
          DO I = 1,NVAR
             H((I*(I+1))/2) = 1.0_DP
          END DO
          STEPMA = STEPMAX
      END IF

      ! . Project out rotations and translations.
      CALL PROJECT_ROTATION_TRANSLATION ( H, .FALSE. )

      ! . Diagonalize the matrix.
      CALL SYMMETRIC_UPPER ( H, EIGENVALUES, EIGENVECTORS )

      ! . Make small eigenvalues zero.
      WHERE ( ABS ( EIGENVALUES ) < THRES2 ) EIGENVALUES = 0.0_DP

      ! . Determine the number of zero eigenvalues.
      NZERO = COUNT ( EIGENVALUES == 0.0_DP )

      ! . Compress the eigenvalue and eigenvector arrays to remove zero eigenvalues.
      IF ( NZERO > 0 ) THEN
         NZ = 0
         DO I = 1,NVAR
            IF ( EIGENVALUES(I) /= 0.0_DP ) THEN
               NZ = NZ + 1
               EIGENVALUES(NZ)    = EIGENVALUES(I)
               EIGENVECTORS(:,NZ) = EIGENVECTORS(:,I)
            END IF
         END DO
      END IF

      ! . Shift the eigenvalues.
      SHIFT = 0.0_DP
      IF ( EIGENVALUES(1) < THRES1 ) SHIFT = THRES1 - EIGENVALUES(1)
      EIGENVALUES(1:NVAR-NZERO) = EIGENVALUES(1:NVAR-NZERO) + SHIFT

      ! . Invert the eigenvalues.
      WHERE ( EIGENVALUES(1:NVAR-NZERO) /= 0.0_DP ) EIGENVALUES(1:NVAR-NZERO) = 1.0_DP / EIGENVALUES(1:NVAR-NZERO)

      ! . Construct the starting inverse Hessian.
      DO I = 1,NVAR
         DO J = 1,I
            HINV0(I,J) = SUM ( EIGENVALUES(1:NVAR-NZERO) * EIGENVECTORS(I,1:NVAR-NZERO) * EIGENVECTORS(J,1:NVAR-NZERO) )
            HINV0(J,I) = HINV0(I,J)
         END DO
      END DO

      END SUBROUTINE SETUP_INVERSE_HESSIAN

      !---------------------
      SUBROUTINE UPDATE_BFGS
      !---------------------

      ! . Scalars.
      INTEGER            :: I, J
      REAL ( KIND = DP ) :: A, AB, B, BFG, DFP

      ! . Arrays.
      REAL ( KIND = DP ), DIMENSION(1:NVAR) :: Z

      ! . Get some intermediate factors.
      Z  = MATMUL ( HINV, GDEL )
      A  = DOT_PRODUCT ( XDEL, GDEL )
      B  = DOT_PRODUCT ( Z,    GDEL )
      AB = ( A + B ) / ( A * A )

      ! . Do the Broyden rank 2 update.
      DO I = 1,NVAR
         DO J = 1,I
            DFP =      XDEL(I) * XDEL(J) / A  -             Z(I)           * Z(J)   / B
            BFG = AB * XDEL(I) * XDEL(J)      - ( XDEL(I) * Z(J) + XDEL(J) * Z(I) ) / A
            HINV(I,J) = HINV(I,J) + BROY * DFP  + ( 1.0_DP - BROY ) * BFG
            HINV(J,I) = HINV(I,J)
         END DO
      END DO

      END SUBROUTINE UPDATE_BFGS

      !---------------------------------------
      SUBROUTINE VARIABLES_EMPTY ( X, ATMDAT )
      !---------------------------------------

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(:),   INTENT(IN)  :: X
      REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: ATMDAT

      ! . Local scalars.
      INTEGER :: IATOM, II

      ! . Loop over the free atoms.
      II = -3
      DO IATOM = 1,NATOMS
         IF ( .NOT. ATMFIX(IATOM) ) THEN
            II = II + 3
            ATMDAT(1:3,IATOM) = X(II+1:II+3)
         END IF
      END DO

      END SUBROUTINE VARIABLES_EMPTY

      !--------------------------------------
      SUBROUTINE VARIABLES_FILL ( X, ATMDAT )
      !--------------------------------------

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT) :: X
      REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)  :: ATMDAT

      ! . Local scalars.
      INTEGER :: IATOM, II

      ! . Loop over the free atoms.
      II = -3
      DO IATOM = 1,NATOMS
         IF ( .NOT. ATMFIX(IATOM) ) THEN
            II = II + 3
            X(II+1:II+3) = ATMDAT(1:3,IATOM)
         END IF
      END DO

      END SUBROUTINE VARIABLES_FILL

   END SUBROUTINE OPTIMIZE_BFGS_CALCULATE

END MODULE OPTIMIZE_BFGS
