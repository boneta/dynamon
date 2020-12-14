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
!
! . Contributed by Ignacio Galvan.
!
!===============================================================================
MODULE LBFGS

USE DEFINITIONS,      ONLY : DP
USE ATOMS,            ONLY : ATMCRD, ATMFIX, NATOMS, NFREE
USE POTENTIAL_ENERGY, ONLY : ETOTAL, ATMDER, GRADIENT
USE PRINTING

IMPLICIT NONE
PRIVATE
PUBLIC :: LBFGS_INITIALIZE, LBFGS_DATA, LBFGS_STEP, OPTIMIZE_LBFGS

REAL ( KIND = DP ), DIMENSION(:),   ALLOCATABLE :: X, G, RHO
REAL ( KIND = DP ), DIMENSION(:,:), ALLOCATABLE :: DX, DG
REAL ( KIND = DP ) :: HINIT

INTEGER :: N, M, ITER

!======================================================================

CONTAINS

!======================================================================

SUBROUTINE LBFGS_INITIALIZE ( NDIM, NSTEP, DELETE )

  INTEGER, INTENT(IN) :: NDIM, NSTEP
  LOGICAL, INTENT(IN), OPTIONAL :: DELETE

  LOGICAL :: DEL

  N = NDIM
  M = NSTEP
  ITER = 0

  DEL = .FALSE.
  IF ( PRESENT(DELETE) ) DEL = DELETE
  IF ( ALLOCATED(X)    ) DEALLOCATE( X )
  IF ( ALLOCATED(G)    ) DEALLOCATE( G )
  IF ( ALLOCATED(RHO)  ) DEALLOCATE( RHO )
  IF ( ALLOCATED(DX)   ) DEALLOCATE( DX )
  IF ( ALLOCATED(DG)   ) DEALLOCATE( DG )

  IF (.NOT. DEL) &
    ALLOCATE( X(N), G(N), RHO(M), DX(M,N), DG(M,N) )

  HINIT = 1.0_DP

END SUBROUTINE LBFGS_INITIALIZE

!======================================================================

SUBROUTINE LBFGS_DATA ( COORD, GRAD )

  REAL ( KIND = DP ), DIMENSION(N), INTENT(IN) :: COORD, GRAD

  REAL ( KIND = DP ) :: YS, YY
  INTEGER :: I

  IF ( ITER > M ) THEN
    DX = CSHIFT ( DX, 1 )
    DG = CSHIFT ( DG, 1 )
    RHO = CSHIFT ( RHO, 1 )
  END IF

  IF ( ITER > 0 ) THEN
    I = MIN(ITER,M)
    DX(I,:) = COORD - X
    DG(I,:) = GRAD - G
    YS = DOT_PRODUCT( DX(I,:), DG(I,:) )
    YY = DOT_PRODUCT( DG(I,:), DG(I,:) )
    RHO(I) = 1.0_DP / YS
    HINIT = YS / YY
  END IF
  X = COORD
  G = GRAD

END SUBROUTINE LBFGS_DATA

!======================================================================

FUNCTION LBFGS_STEP ( )

  REAL ( KIND = DP ), DIMENSION(N) :: LBFGS_STEP

  REAL ( KIND = DP ), DIMENSION(M) :: AUX
  INTEGER :: I

  LBFGS_STEP = -G

  DO I = MIN(ITER,M), 1, -1
    AUX(I) = RHO(I) * DOT_PRODUCT( DX(I,:), LBFGS_STEP )
    LBFGS_STEP = LBFGS_STEP - AUX(I) * DG(I,:)
  END DO

  LBFGS_STEP = HINIT*LBFGS_STEP

  DO I = 1, MIN(ITER,M)
    AUX(I) = AUX(I) - RHO(I) * DOT_PRODUCT( DG(I,:), LBFGS_STEP )
    LBFGS_STEP = LBFGS_STEP + AUX(I) * DX(I,:)
  END DO

  ITER = ITER + 1

  IF ( ITER == 1 ) LBFGS_STEP = LBFGS_STEP / DOT_PRODUCT( G, G )

END FUNCTION LBFGS_STEP

!======================================================================

SUBROUTINE OPTIMIZE_LBFGS ( NUPD, STATUS, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, GRADIENT_TOLERANCE )

  ! . Optional scalar arguments.
  INTEGER,            INTENT(IN),  OPTIONAL :: PRINT_FREQUENCY, STEP_NUMBER, NUPD
  INTEGER,            INTENT(OUT), OPTIONAL :: STATUS
  REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: GRADIENT_TOLERANCE, STEP_SIZE

  ! . Local allocatable arrays.
  REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: COORD, GRAD, STEP

  ! . Local scalars.
  REAL ( KIND = DP ) :: E, G2, STEP_SCALE, TOLGRD, STPSIZ, GTOLERANCE
  INTEGER :: I, NPRINT, NSTEP, NVAR

  M = 4
  NPRINT = 0
  NSTEP  = 0
  STPSIZ = 1.0E-1_DP
  TOLGRD = 1.0E-3_DP
  NVAR = 3*NFREE

  IF ( PRESENT ( NUPD               ) ) M = NUPD
  IF ( PRESENT ( GRADIENT_TOLERANCE ) ) TOLGRD = GRADIENT_TOLERANCE
  IF ( PRESENT ( PRINT_FREQUENCY    ) ) NPRINT = PRINT_FREQUENCY
  IF ( PRESENT ( STEP_NUMBER        ) ) NSTEP  = STEP_NUMBER
  IF ( PRESENT ( STEP_SIZE          ) ) STPSIZ = STEP_SIZE

  GTOLERANCE = REAL ( NVAR, DP ) * TOLGRD * TOLGRD

  ! . Initialize the STATUS counter if necessary.
  IF ( PRESENT(STATUS) ) STATUS=0

  ! . Check the number of free atoms.
  IF ( NFREE < 0 ) RETURN

  ! . Check the input parameters.
  IF ( ( NPRINT < 0 ) .OR. ( NPRINT > NSTEP ) ) NPRINT = 0

  ! . Print out the control parameters.
  IF ( NPRINT > 0 ) THEN

    ! . Print out the options.
    CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#55FF55", VARIABLEWIDTH = 16 )
    CALL PRINT_SUMMARY_START ( "L-BFGS Minimization Calculation" )
    WRITE ( PRINT_LINE, "(I16)"    ) NVAR   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Variables"  )
    WRITE ( PRINT_LINE, "(I16)"    ) NSTEP  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations" )
    WRITE ( PRINT_LINE, "(I16)"    ) NPRINT ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency"      )
    WRITE ( PRINT_LINE, "(G16.10)" ) TOLGRD ; CALL PRINT_SUMMARY_ELEMENT ( "Gradient Tolerance"   )
    WRITE ( PRINT_LINE, "(G16.10)" ) STPSIZ ; CALL PRINT_SUMMARY_ELEMENT ( "Step Size"            )
    WRITE ( PRINT_LINE, "(I16)"    ) M      ; CALL PRINT_SUMMARY_ELEMENT ( "L-BFGS Update Steps"  )
    CALL PRINT_SUMMARY_STOP

    ! . Print out the header for the iterations.
    CALL PRINT_TABLE_OPTIONS ( COLUMNS = 3, HEADER_COLOR = "#55FF55", PAGEWIDTH = 60, VARIABLEWIDTHS = (/ 10, 20, 20 /) )
    CALL PRINT_TABLE_START
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",         HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Function",     HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS Gradient", HEADER = .TRUE. )

  END IF

  ! . Allocate space for the variables.
  ALLOCATE ( COORD(1:NVAR), GRAD(1:NVAR), STEP(1:NVAR) )

  ! . Move the coordinates to X.
  CALL VARIABLES_FILL ( COORD, ATMCRD )

  ! . Initialize the L-BFGS minimizer.
  CALL LBFGS_INITIALIZE ( NVAR, M )

  ! . Get the initial gradient.
  CALL VARIABLES_EMPTY ( COORD, ATMCRD )
  CALL GRADIENT ( PRINT = .FALSE. ) ; E = ETOTAL
  CALL VARIABLES_FILL ( GRAD, ATMDER )
  G2 = DOT_PRODUCT( GRAD, GRAD )

  ! . Do some printing.
  IF ( NPRINT == 1 ) THEN
    WRITE ( PRINT_LINE, "(I10)"   ) 0                               ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F20.8)" ) E                               ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F20.8)" ) SQRT ( G2 / REAL ( NVAR, DP ) ) ; CALL PRINT_TABLE_ELEMENT
  END IF

  DO I = 1, NSTEP
    IF ( G2 <= GTOLERANCE ) EXIT

    ! . Get the L-BFGS step.
    CALL LBFGS_DATA ( COORD, -GRAD )
    STEP = LBFGS_STEP ( )

    ! . Scale down the step if necessary.
    STEP_SCALE = SQRT( DOT_PRODUCT( STEP, STEP ) / REAL( NVAR, DP ) )
    STEP_SCALE = MIN( STPSIZ/STEP_SCALE, 1.0_DP )
    STEP = STEP_SCALE*STEP

    ! . Move the coordinates.
    COORD = COORD + STEP

    ! . Get the gradient
    CALL VARIABLES_EMPTY ( COORD, ATMCRD )
    CALL GRADIENT ( PRINT = .FALSE. ) ; E = ETOTAL
    CALL VARIABLES_FILL ( GRAD, ATMDER )
    G2 = DOT_PRODUCT( GRAD, GRAD )

    ! . Do some printing.
    IF ( NPRINT > 0 ) THEN
      IF ( MOD ( I, NPRINT ) == 0 ) THEN
        WRITE ( PRINT_LINE, "(I10)"   ) I                               ; CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F20.8)" ) E                               ; CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F20.8)" ) SQRT ( G2 / REAL ( NVAR, DP ) ) ; CALL PRINT_TABLE_ELEMENT
      END IF
    END IF
  END DO

  ! . Do some more printing.
  IF ( NPRINT > 0 ) THEN
    ! . Print out the terminator.
    CALL PRINT_TABLE_STOP

    ! . Print out the status.
    IF ( G2 <= GTOLERANCE ) THEN
      CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Gradient tolerance reached." )
    ELSE
      CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Too many steps."             )
    END IF
  END IF

  ! . Save the optimized coordinates.
  CALL VARIABLES_EMPTY ( COORD, ATMCRD )

  ! . Deallocate the temporary array.
  DEALLOCATE ( COORD, GRAD, STEP )

  CONTAINS

  SUBROUTINE VARIABLES_EMPTY ( COORD, ATMDAT )

  ! . Array arguments.
  REAL ( KIND = DP ), DIMENSION(:),   INTENT(IN)  :: COORD
  REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: ATMDAT

  ! . Local scalars.
  INTEGER :: IATOM, II

  ! . Loop over the free atoms.
  II = -3
  DO IATOM = 1,NATOMS
    IF ( .NOT. ATMFIX(IATOM) ) THEN
      II = II + 3
      ATMDAT(1:3,IATOM) = COORD(II+1:II+3)
    END IF
  END DO

  END SUBROUTINE VARIABLES_EMPTY

  SUBROUTINE VARIABLES_FILL ( COORD, ATMDAT )

  ! . Array arguments.
  REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT) :: COORD
  REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)  :: ATMDAT

  ! . Local scalars.
  INTEGER :: IATOM, II

  ! . Loop over the free atoms.
  II = -3
  DO IATOM = 1,NATOMS
    IF ( .NOT. ATMFIX(IATOM) ) THEN
      II = II + 3
      COORD(II+1:II+3) = ATMDAT(1:3,IATOM)
    END IF
  END DO

  END SUBROUTINE VARIABLES_FILL

END SUBROUTINE OPTIMIZE_LBFGS

!======================================================================

END MODULE LBFGS
