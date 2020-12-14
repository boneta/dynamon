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
!                          The NEB-SPLINE Module
!===============================================================================
!
! . Subroutines:
!
!   NEB_OPTIMIZE_SPLINE    Optimize a chain trajectory.
!   VERLET_MINIMIZE_MAX    Minimize the forces using a damped verlet algorithm
!                          moving only one image at a time.
!   LBFGS_MINIMIZE         Minimize the forces using the L-BFGS algorithm
!   GET_LINEAR_CONSTRAINTS Impose Ekart conditions for translations and rotations.
!   PROJECT                Project the velocity onto the force vector.
!   TANGENT                Calculates the tangent at a chain point.
!   PARGRAD                Calculates the gradient parallel to the tangent.
!   INTERPOLATE            Does a polynomial interpolation of the images.
!   FORTOT                 Calculates the force acting on the images.
!   BUILD_SPLINE           Generates a interpolating spline for the images.
!   REDISTRIBUTE           Redistributes the images evenly along the spline.
!   DISTANCE               Calculate the distance between a pair of images.
!   CRD_TO_X               Copy ATMCRD to the variable array X.
!   X_TO_CRD               Copy X to ATMCRD.
!   CUBIC_SPINE            Calculates the spline coefficients.
!   INTEGRATE_SPLINE       Calculate the length of a spline by integration.
!   FIND_POSITION          Find a point in a spline at a given length.
!   SPLINE_PROP            Calculate geometrical properties for a spline.
!   FIND_MAX               Find a possible maximum in a spline segment.
!   GROWING_STRING         Generates an initial guess path.
!   CHAIN_INTERPOLATE      Generates a new chain as a spline interpolation.
!===============================================================================
MODULE NEB_SPLINE

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_LINE, PRINT_LINEBREAK, PRINT_PARAGRAPH, PRINT_PARAGRAPH_START, &
                        PRINT_PARAGRAPH_STOP, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS,  &
                        PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT,        &
                        PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP, PRINT_TEXT
USE CONSTANTS

USE ATOMS,              ONLY : ATMCRD, ATMFIX, NATOMS, NATOMSQM, NFIXED, NFREE
USE DCD_IO
USE POTENTIAL_ENERGY,   ONLY : ATMDER, ENERGY, GRADIENT, ETOTAL, EQM
USE STRING,             ONLY : ENCODE_INTEGER
USE MOPAC_DENSITY,      ONLY : DENSITY_READ, DENSITY_WRITE, DENSITY_GUESS
USE LINEAR_ALGEBRA,     ONLY : NORM
USE SUPERIMPOSE,        ONLY : SUPERIMPOSE_QUATERNION
USE LBFGS

IMPLICIT NONE
PRIVATE
PUBLIC :: NEB_OPTIMIZE_SPLINE, GROWING_STRING, CHAIN_INTERPOLATE
#ifndef PGPC
SAVE
#endif

! . Module scalars.
INTEGER            :: NCHAIN, NCONSTRAINTS, N3

! . Module arrays.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)     :: D1, ECHAIN, ECHAINQM,  XBEG, XEND, &
                                                     GBEG, GEND
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:)   :: G_CONSTRAINT
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:,:) :: SPLINE

!==============================================================================
CONTAINS
!==============================================================================


!-------------------------------------------------------------------------------------------------
SUBROUTINE NEB_OPTIMIZE_SPLINE (FILE_IN, FILE_OUT, PRINT_FREQUENCY, STEP_NUMBER,        &
                  GRADIENT_TOLERANCE, STEP_SIZE, NBINS, SAVE_DENSITIES, READ_DENSITIES, &
                  TEMP, METHOD, OPTIM_METHOD, OPTIM_STEP, MAX_FILE, INTERP_FILE )
!-------------------------------------------------------------------------------------------------

  ! . Scalar arguments.
  CHARACTER ( LEN = * ), INTENT(IN) :: FILE_IN, FILE_OUT
  ! . Optional scalar arguments for the CHAIN energy.
  LOGICAL, INTENT(IN), OPTIONAL :: SAVE_DENSITIES, READ_DENSITIES

  ! . Optional scalar arguments for the minimization routine.
  CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: METHOD, OPTIM_METHOD
  INTEGER,               INTENT(IN), OPTIONAL :: PRINT_FREQUENCY, STEP_NUMBER, OPTIM_STEP
  REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: GRADIENT_TOLERANCE, STEP_SIZE, TEMP

  ! . Optional scalar argument for the final interpolation.
  CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: MAX_FILE, INTERP_FILE
  INTEGER,               INTENT(IN), OPTIONAL :: NBINS

  ! . Local scalars.
  INTEGER            :: I, II, NPOINTS, METH, OPMETH, ENEVAL, MINI_ITER
  REAL ( KIND = DP ) :: EREF, T, BETA
  LOGICAL :: SAVEDE, READDE
  CHARACTER (LEN = 256 ) :: MAXFILE, INTERPFILE

  ! . Local arrays.
  REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)     :: G, X, PFLUX, SPLINE_LENGTH, SPLINE_CURV
  REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:)   :: TMPCRD, ENERGY_SPLINE
  CHARACTER (LEN = 3 ) :: ISTRING

  ! . Local types.
  TYPE(DCD_TYPE) :: DCD_IN, DCD_OUT

  ! . There are no free atoms.
  IF ( NFREE <= 0 ) RETURN

  ! . Initialize the CHAIN energy parameters.
  IF ( PRESENT ( SAVE_DENSITIES) ) THEN
     SAVEDE = SAVE_DENSITIES
  ELSE IF (NATOMSQM>0) THEN
    SAVEDE = .TRUE.
  ELSE
    SAVEDE = .FALSE.
  END IF

  IF ( PRESENT ( READ_DENSITIES) ) THEN
    READDE = READ_DENSITIES
  ELSE
    READDE = .FALSE.
  END IF

  T = 0.0_DP
  METH = 0
  IF ( PRESENT(TEMP) ) THEN
    T = TEMP
    BETA = 1.0_DP / R / TEMP
    METH = 1
    IF ( T < 1.0_DP ) METH = 0
    IF ( PRESENT(METHOD) ) THEN
      SELECT CASE ( METHOD )
      CASE ( "DIFF" )
        METH = 1
      CASE ( "INT" )
        METH = 2
      CASE ( "INT2" )
        METH = 3
      CASE DEFAULT
        CALL PRINT_ERROR ( "NEB_OPTIMIZE_SPLINE", "METHOD can only be INT, INT2 or DIFF." )
      END SELECT
    END IF
  END IF

  ! . Default is LBFGS
  OPMETH = 2
  IF ( PRESENT(OPTIM_METHOD) ) THEN
    SELECT CASE ( OPTIM_METHOD )
    CASE ( "VERLET" )
      OPMETH = 1
    CASE ( "LBFGS" )
      OPMETH = 2
    CASE DEFAULT
      CALL PRINT_ERROR ( "NEB_OPTIMIZE_SPLINE", "Unrecognized optimization method." )
    END SELECT
  END IF

  IF ( PRESENT(OPTIM_STEP) ) THEN
    MINI_ITER = OPTIM_STEP
  ELSE
    MINI_ITER = 20
  END IF

  MAXFILE = ''
  IF ( PRESENT(MAX_FILE) ) MAXFILE = MAX_FILE
  INTERPFILE = ''
  IF ( PRESENT(INTERP_FILE) ) INTERPFILE = INTERP_FILE

  IF ((SAVEDE .OR. READDE) .AND. NATOMSQM==0) CALL PRINT_ERROR("NEB_OPTIMIZE_SPLINE", &
   "If there are no QM atoms, the density cannot be saved or read.")

  ! . Activate the input trajectory.
  CALL DCD_INITIALIZE ( DCD_IN )
  CALL DCD_ACTIVATE_READ ( FILE_IN, DCD_IN )

  ! . Get the number of points in the chain.
  NPOINTS = DCD_IN%NFRAMES

  ! . Print out the options.
  CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#8888FF", VARIABLEWIDTH = 16 )
  CALL PRINT_SUMMARY_START ( "Nudged-Elastic-Band Spline Path Calculation" )
  WRITE ( PRINT_LINE, "(I16)"    ) NPOINTS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Structures" )
  WRITE ( PRINT_LINE, "(G16.10)" ) T       ; CALL PRINT_SUMMARY_ELEMENT ( "Temperature"          )
  CALL PRINT_SUMMARY_STOP

  ! . Write out the trajectory names.
  CALL PRINT_PARAGRAPH_START
  CALL PRINT_TEXT ( TEXT = "Chain Optimization:  Input Trajectory = "//TRIM ( FILE_IN  ) )
  CALL PRINT_LINEBREAK
  CALL PRINT_TEXT ( TEXT = "                    Output Trajectory = "//TRIM ( FILE_OUT ) )
  CALL PRINT_PARAGRAPH_STOP

  ! . Calculate some counters.
  NCHAIN = NPOINTS - 2
  N3     = 3 * NFREE

  ! . Allocate some temporary arrays.
  ALLOCATE ( D1(1:NCHAIN+1), ECHAIN(1:NCHAIN+2), ECHAINQM(1:NCHAIN+2), G(1:NCHAIN*N3), &
      G_CONSTRAINT(1:N3,6), TMPCRD(1:3,1:NATOMS), X(1:NCHAIN*N3), &
      XBEG(1:N3), XEND(1:N3), GBEG(1:N3), GEND(1:N3), PFLUX(1:NCHAIN+1), &
      SPLINE(1:4,1:NCHAIN+1,1:N3), SPLINE_LENGTH(1:NCHAIN+1), SPLINE_CURV(1:NCHAIN+2), &
      ENERGY_SPLINE(1:4,1:NCHAIN+2) )

  ! . Save the existing coordinates.
  TMPCRD = ATMCRD

  ! . Read in the first structure from the trajectory, calculate its energy and save its coordinates.
  ! . Note that if there are free atoms their coordinates are ALWAYS present in ATMCRD.
  CALL DCD_READ ( DCD_IN, ATMCRD )
  CALL DENSITY_GUESS
  CALL GRADIENT ( PRINT = .FALSE. ) ; ECHAIN(1) = ETOTAL; ECHAINQM(1)=EQM
  CALL CRD_TO_X ( ATMCRD, XBEG )
  CALL CRD_TO_X (ATMDER, GBEG)

  ! . Read in the intermediate structures from the trajectory.
  DO I = 1,NCHAIN
    CALL DCD_READ ( DCD_IN, ATMCRD )
    II = ( I - 1 ) * N3
    CALL CRD_TO_X ( ATMCRD, X(II+1:II+N3) )
  END DO

  ! . Read in the last structure from the trajectory, calculate its energy and save its coordinates.
  CALL DCD_READ ( DCD_IN, ATMCRD )
  CALL DENSITY_GUESS
  CALL GRADIENT ( PRINT = .FALSE. ) ; ECHAIN(NCHAIN+2) = ETOTAL; ECHAINQM(NCHAIN+2)=EQM
  CALL CRD_TO_X ( ATMCRD, XEND )
  CALL CRD_TO_X ( ATMDER, GEND )

  EREF = MAX(ECHAIN(1), ECHAIN(NCHAIN+2))

  ! . Deactivate the trajectory.
  CALL DCD_DEACTIVATE ( DCD_IN  )

  ! . Calculate the linear constraints using the central structure.
  II = ( ( NCHAIN + 1 ) / 2 - 1 ) * N3
  CALL GET_LINEAR_CONSTRAINTS ( X(II+1:II+N3) )

  ! . Perform the minimization.
  SELECT CASE ( OPMETH )
  CASE (1)
    CALL VERLET_MINIMIZE_MAX (X, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, GRADIENT_TOLERANCE)
  CASE (2)
    CALL LBFGS_MINIMIZE (X, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, GRADIENT_TOLERANCE)
  END SELECT

  ! . Activate the output trajectory.
  CALL DCD_INITIALIZE ( DCD_OUT )
  CALL DCD_ACTIVATE_WRITE ( FILE_OUT, DCD_OUT, "CORD", NATOMS, NFIXED, NPOINTS, &
                            QFIX = ATMFIX, PRINT=.FALSE. )

  ! . Write out the first structure of the chain to the output trajectory.
  CALL X_TO_CRD ( XBEG, ATMCRD )
  CALL DCD_WRITE ( DCD_OUT, ATMCRD )

  ! . Write out the intermediate structures of the chain.
  DO I = 1,NCHAIN
     II = ( I - 1 ) * N3
     CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
     CALL DCD_WRITE ( DCD_OUT, ATMCRD )
  END DO

  ! . Write out the last structure of the chain.
  CALL X_TO_CRD ( XEND, ATMCRD )
  CALL DCD_WRITE ( DCD_OUT, ATMCRD )

  ! . Deactivate the trajectory.
  CALL DCD_DEACTIVATE ( DCD_OUT )

  ! . Calculate and print the flux.
  IF ( METH > 0 ) THEN
    CALL FLUX ( PFLUX )
    CALL PRINT_TABLE_OPTIONS ( COLUMNS = 2, HEADER_COLOR = "#55FFFF", PAGEWIDTH = 40, &
                               VARIABLEWIDTHS = (/ 8, 16 /) )
    CALL PRINT_TABLE_START
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Flux", COLSPAN = 2, HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Structures",        HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Partial 1/Flux",    HEADER = .TRUE. )
    DO I = 1, NCHAIN+1
      WRITE ( PRINT_LINE, "(I3,', ',I3)"    ) I, I+1   ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(E16.9)"         ) PFLUX(I) ; CALL PRINT_TABLE_ELEMENT
    END DO
    CALL PRINT_TABLE_ELEMENT ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(A,E16.9)" ) "Flux proportional to ", 1.0_DP / SUM(PFLUX)
    CALL PRINT_TABLE_ELEMENT( COLSPAN = 2 )
    CALL PRINT_TABLE_STOP
  END IF

  ! . Write out the header.
  CALL PRINT_TABLE_OPTIONS ( COLUMNS = 6, HEADER_COLOR = "#8888FF", VARIABLEWIDTHS = (/ 10, 15, 15, 15, 10, 12 /) )
  CALL PRINT_TABLE_START
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Chain Structures", COLSPAN = 6, HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Structure",        HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Energy",           HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "QM Energy",        HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Gradient",         HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Curvature",        HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Dist(i,i+1)",      HEADER = .TRUE. )

  ! . Loop over the structures.
  DO I = 1,NCHAIN+2
    WRITE ( PRINT_LINE, "(I10)"   ) I              ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) ECHAIN(I)      ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) ECHAINQM(I)    ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) PARGRAD(I)     ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F10.5)" ) SPLINE_CURV(I) ; CALL PRINT_TABLE_ELEMENT
    IF ( I == ( NCHAIN + 2 ) ) THEN
      WRITE ( PRINT_LINE, "(F12.5)" ) !0.0
    ELSE
      WRITE ( PRINT_LINE, "(F12.5)" ) D1(I)
    END IF
    CALL PRINT_TABLE_ELEMENT
  END DO
  CALL PRINT_TABLE_STOP

  ! . Do the simplex interpolation.
  CALL INTERPOLATE (NBINS)
  ! . Reset the coordinates.
  ATMCRD = TMPCRD

  ! . Deallocate the temporary arrays.

  DEALLOCATE ( D1, ECHAIN, ECHAINQM, G, G_CONSTRAINT, TMPCRD, X, XBEG, XEND, GBEG, GEND, PFLUX, &
               SPLINE, SPLINE_LENGTH, SPLINE_CURV, ENERGY_SPLINE )


  !===========================================================================
  CONTAINS
  !===========================================================================

  !--------------------------------------
  SUBROUTINE GET_LINEAR_CONSTRAINTS ( X )
  !--------------------------------------

  ! . Array arguments.
  REAL ( KIND = DP ), DIMENSION(1:N3), INTENT(IN) :: X

  ! . Local scalars.
  INTEGER            :: I, II, IVEC, JVEC
  REAL ( KIND = DP ) :: DOTFAC

  ! . Initialize NCONSTRAINTS.
  NCONSTRAINTS = 0

  ! . Return if there are fixed atoms.
  IF ( NFIXED > 0 ) RETURN

  ! . Initialize the number of constraint vectors.
  NCONSTRAINTS = 6

  ! . Initialize the vectors.
  G_CONSTRAINT = 0.0_DP

  ! . Loop over the atoms.
  DO I = 1,NATOMS

    ! . Calculate the array index.
    II = 3 * ( I - 1 )

    ! . The translations.
    G_CONSTRAINT(II+1,1) = 1.0_DP
    G_CONSTRAINT(II+2,2) = 1.0_DP
    G_CONSTRAINT(II+3,3) = 1.0_DP

    ! . The rotations.
    G_CONSTRAINT(II+2,4) =   X(II+3)
    G_CONSTRAINT(II+3,4) = - X(II+2)
    G_CONSTRAINT(II+1,5) = - X(II+3)
    G_CONSTRAINT(II+3,5) =   X(II+1)
    G_CONSTRAINT(II+1,6) =   X(II+2)
    G_CONSTRAINT(II+2,6) = - X(II+1)

  END DO

  ! . Orthonormalize the vectors.
  DO IVEC = 1,NCONSTRAINTS

    ! . Subtract the contributions from vectors of lower index.
    DO JVEC = 1,(IVEC-1)
      DOTFAC = DOT_PRODUCT ( G_CONSTRAINT(1:N3,IVEC), G_CONSTRAINT(1:N3,JVEC) )
      G_CONSTRAINT(1:N3,IVEC) = G_CONSTRAINT(1:N3,IVEC) - DOTFAC * G_CONSTRAINT(1:N3,JVEC)
    END DO

    ! . Normalize the vector.
    DOTFAC = NORM ( G_CONSTRAINT(1:N3,IVEC))
    G_CONSTRAINT(1:N3,IVEC) = G_CONSTRAINT(1:N3,IVEC) / DOTFAC

  END DO

  END SUBROUTINE GET_LINEAR_CONSTRAINTS

  !-------------------------------------------------------------------
  SUBROUTINE CONSTRAIN(VECTOR)
  !-------------------------------------------------------------------

  REAL ( KIND = DP ), DIMENSION(N3) :: VECTOR
  INTEGER :: I

  IF ( NCONSTRAINTS > 0 ) THEN
    DO I = 1,NCONSTRAINTS
      VECTOR = VECTOR - DOT_PRODUCT ( G_CONSTRAINT(1:N3,I), VECTOR ) * G_CONSTRAINT(1:N3,I)
    END DO
  END IF

  END SUBROUTINE CONSTRAIN

  !-------------------------------------------------------------------
  SUBROUTINE VERLET_MINIMIZE_MAX (X, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, FORTOL,  &
                              & STATUS)
  !-------------------------------------------------------------------

  ! . Array arguments.
  REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: X

  ! . Optional arguments.
  INTEGER,            INTENT(IN),  OPTIONAL :: PRINT_FREQUENCY, STEP_NUMBER
  INTEGER,            INTENT(OUT), OPTIONAL :: STATUS
  REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: FORTOL, STEP_SIZE

  ! . Local scalars.
  INTEGER            :: I, II, J, IFAIL, NPRINT, NSTEP, NVAR
  REAL ( KIND = DP ) :: GTOLERANCE, STPSIZ, TOLGRD

  ! . Local arrays
  REAL (KIND = DP), DIMENSION (1:SIZE(X)) :: VEL, FORCENEW, FORCE
  REAL (KIND = DP), DIMENSION (1:NCHAIN) :: MODFORCE
  INTEGER :: IMOVE

  !---------------------------------------------------------------------------
  ! . Initialization.
  !---------------------------------------------------------------------------
  ! . Initialize the IFAIL counter.
  IFAIL = 0

  ! . Find the number of variables.
  NVAR = SIZE ( X )

  ! . Skip the calculation if there are no variables.
  IF ( NVAR <= 0 ) RETURN

  ! . Initialize some algorithm options.
  NPRINT = 0
  NSTEP  = 0
  STPSIZ = 1.0E-3_DP
  TOLGRD = 2.0E-2_DP

  ! . Assign the input parameters.
  IF ( PRESENT ( FORTOL             ) ) TOLGRD = FORTOL
  IF ( PRESENT ( PRINT_FREQUENCY    ) ) NPRINT = PRINT_FREQUENCY
  IF ( PRESENT ( STEP_NUMBER        ) ) NSTEP  = STEP_NUMBER
  IF ( PRESENT ( STEP_SIZE          ) ) STPSIZ = STEP_SIZE

  ! . Check the input parameters.
  IF ( ( NPRINT < 0 ) .OR. ( NPRINT > NSTEP ) ) NPRINT = 0

  ! . Calculate the gradient tolerance.
  GTOLERANCE = REAL ( N3, DP ) * TOLGRD * TOLGRD

  ! . Print out the control parameters.
  IF ( NPRINT > 0 ) THEN

    ! . Print out the options.
    CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#55FF55", VARIABLEWIDTH = 16 )
    CALL PRINT_SUMMARY_START ( "Verlet Minimization Calculation" )
    WRITE ( PRINT_LINE, "(I16)"    ) NVAR   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Variables"  )
    WRITE ( PRINT_LINE, "(I16)"    ) NSTEP  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations" )
    WRITE ( PRINT_LINE, "(G16.10)" ) TOLGRD ; CALL PRINT_SUMMARY_ELEMENT ( "Force Tolerance"      )
    WRITE ( PRINT_LINE, "(I16)"    ) NPRINT ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency"      )
    WRITE ( PRINT_LINE, "(G16.10)" ) STPSIZ ; CALL PRINT_SUMMARY_ELEMENT ( "Step Size"            )
    CALL PRINT_SUMMARY_STOP

    ! . Print out the header for the iterations.
    CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#55FF55", PAGEWIDTH = 80, &
          & VARIABLEWIDTHS = (/ 6, 12, 15, 13, 13 /) )
    CALL PRINT_TABLE_START
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",           HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Grad. eval.",    HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Max. RMS Force", HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Path length",    HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Max. curv.",     HEADER = .TRUE. )

  END IF

  ! .  Calculate the initial energies.
  DO I = 1,NCHAIN
    II = ( I - 1 ) * N3
    CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
    IF (READDE) THEN
      CALL ENCODE_INTEGER(I,ISTRING)
      CALL DENSITY_READ ("dens"//FILE_IN//TRIM(ISTRING), PRINT=.FALSE.)
    ELSE
      IF ( SAVEDE .AND. I > 1 ) THEN
        CALL ENCODE_INTEGER(I-1,ISTRING)
        CALL DENSITY_READ ("dens"//FILE_OUT//TRIM(ISTRING), PRINT=.FALSE.)
      ELSE
        CALL DENSITY_GUESS
      END IF
    END IF
    CALL ENERGY (PRINT=.FALSE.)
    ECHAIN(I+1)=ETOTAL
    ECHAINQM(I+1)=EQM
    IF (SAVEDE) THEN
      CALL ENCODE_INTEGER(I,ISTRING)
      CALL DENSITY_WRITE ("dens"//FILE_OUT//TRIM(ISTRING), PRINT=.FALSE.)
    END IF
  END DO

  ! . Check the step number.
  IF ( NSTEP <= 0 ) RETURN

  !---------------------------------------------------------------------------
  ! . Set up the minimization.
  !---------------------------------------------------------------------------

  ENEVAL = 0
  CALL BUILD_SPLINE( )
  VEL = 0.0_DP
  DO I = 1, NCHAIN
    CALL GET_GRADIENT(X, NVAR, I)
  END DO
  DO I = 1, NCHAIN
    CALL FORTOT(FORCE, I, MODF2=MODFORCE(I))
  END DO
  FORCENEW = FORCE
  ! . Do some printing.
  IF ( NPRINT > 0 ) THEN
    WRITE ( PRINT_LINE, "(I6)"   ) 1
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(I12)"   ) ENEVAL
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.8)" ) AMUA2PS2_TO_KJMOL*SQRT(MAXVAL(MODFORCE)/REAL(N3,DP))
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F11.8)" ) SUM( SPLINE_LENGTH ) / SQRT( REAL(N3) )
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F11.5)" ) MAXVAL( SPLINE_CURV )
    CALL PRINT_TABLE_ELEMENT
  END IF
  ! . Check for convergence.
  IF (AMUA2PS2_TO_KJMOL**2*MAXVAL(MODFORCE) <= GTOLERANCE) THEN
    GO TO 9999
  ! . Check the number of steps.
  ELSE IF ( NSTEP == 1 ) THEN
    IFAIL = 131
    GO TO 9999
  END IF

  !---------------------------------------------------------------------------
  ! . Perform the minimization.
  !---------------------------------------------------------------------------

  DO I=2, NSTEP
    FORCE = FORCENEW
    IMOVE = SUM(MAXLOC(MODFORCE))
    II = (IMOVE - 1) *N3
    X(II+1:II+N3) = X(II+1:II+N3) + STPSIZ*VEL(II+1:II+N3) + STPSIZ**2 * &
                    FORCE(II+1:II+N3) / 2

    CALL GET_GRADIENT(X, NVAR, IMOVE)
    CALL FORTOT( FORCE, IMOVE, MODF2=MODFORCE(IMOVE) )
    CALL BUILD_SPLINE( )
    IF ( MAXVAL(SPLINE_LENGTH) > 1.5_DP * MINVAL(SPLINE_LENGTH) ) THEN
      CALL REDISTRIBUTE( )
      DO J = 1, NCHAIN
        CALL GET_GRADIENT(X, NVAR, J)
      END DO
      DO J = 1, NCHAIN
        CALL FORTOT(FORCE, J, MODF2=MODFORCE(J))
      END DO
    ELSE
      DO J = IMOVE-1, IMOVE+1
        IF (J < 1 .OR. J > NCHAIN) CYCLE
        CALL FORTOT(FORCE, J, MODF2=MODFORCE(J))
      END DO
    END IF

    FORCENEW = FORCE
    VEL(II+1:II+N3) = VEL(II+1:II+N3) + STPSIZ / 2 * (FORCE(II+1:II+N3) +&
                      FORCENEW(II+1:II+N3))
    CALL PROJECT (VEL, FORCENEW)

    ! . Do some printing.
    IF ( NPRINT > 0 ) THEN
      IF ( MOD ( I, NPRINT ) == 0 ) THEN
        WRITE ( PRINT_LINE, "(I6)"    ) I
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(I12)"   ) ENEVAL
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F15.8)" ) AMUA2PS2_TO_KJMOL*SQRT(MAXVAL(MODFORCE)/REAL(N3,DP))
        CALL PRINT_TABLE_ELEMENT
        WRITE ( print_line, "(F11.8)" ) SUM( SPLINE_LENGTH ) / SQRT( REAL(N3) )
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F11.5)" ) MAXVAL( SPLINE_CURV )
        CALL PRINT_TABLE_ELEMENT

        ! . Activate the output trajectory.
        CALL DCD_INITIALIZE ( DCD_OUT )
        CALL DCD_ACTIVATE_WRITE ( FILE_OUT, DCD_OUT, "CORD", NATOMS, NFIXED, NPOINTS, &
                                  QFIX = ATMFIX, PRINT=.FALSE. )

        ! . Write out the first structure of the chain to the output trajectory.
        CALL X_TO_CRD ( XBEG, ATMCRD )
        CALL DCD_WRITE ( DCD_OUT, ATMCRD )

        ! . Write out the intermediate structures of the chain.
        DO J = 1,NCHAIN
          II = ( J - 1 ) * N3
          CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
          CALL DCD_WRITE ( DCD_OUT, ATMCRD )
        END DO

        ! . Write out the last structure of the chain.
        CALL X_TO_CRD ( XEND, ATMCRD )
        CALL DCD_WRITE ( DCD_OUT, ATMCRD )

        ! . Deactivate the trajectory.
        CALL DCD_DEACTIVATE ( DCD_OUT )
      END IF
    END IF
    IF (AMUA2PS2_TO_KJMOL**2*MAXVAL(MODFORCE) <= GTOLERANCE) EXIT
  END DO
  IF ( I > NSTEP ) THEN
    IFAIL=131
  ELSE
    IFAIL=0
  END IF

  ! . End of the minimization.
  9999 CONTINUE

  !---------------------------------------------------------------------------
  ! . Finish up.
  !---------------------------------------------------------------------------

  ! . Calculate the distances between neighbouring structures.
  D1 = SPLINE_LENGTH / SQRT( REAL(N3) )

  ! . Do some more printing.
  IF ( NPRINT > 0 ) THEN

    ! . Print out the terminator.
    CALL PRINT_TABLE_STOP

    ! . Print out the status.
    SELECT CASE ( IFAIL )
    CASE (   0 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Force tolerance reached.")
    CASE ( 131 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Too many steps."         )
    END SELECT

  END IF

  ! . Assign a value to STATUS if necessary.
  IF ( PRESENT ( STATUS ) ) STATUS = IFAIL

  END SUBROUTINE VERLET_MINIMIZE_MAX


  !-------------------------------------------------------------------
  SUBROUTINE LBFGS_MINIMIZE (X, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, FORTOL,  &
                             STATUS)
  !-------------------------------------------------------------------
  ! . Array arguments.
  REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: X

  ! . Optional arguments.
  INTEGER,            INTENT(IN),  OPTIONAL :: PRINT_FREQUENCY, STEP_NUMBER
  INTEGER,            INTENT(OUT), OPTIONAL :: STATUS
  REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: FORTOL, STEP_SIZE

  ! . Local scalars.
  INTEGER            :: IFAIL, NVAR, NSTEP, NPRINT, I, II, J, K, IMOVE
  REAL ( KIND = DP ) :: GTOLERANCE, STPSIZ, TOLGRD, STPSCL, FCONV

  ! . Local arrays
  REAL (KIND = DP), DIMENSION (1:SIZE(X)) :: FORCE, STEP
  REAL (KIND = DP), DIMENSION (1:NCHAIN)  :: MODFORCE

  !---------------------------------------------------------------------------
  ! . Initialization.
  !---------------------------------------------------------------------------

  ! . Initialize the IFAIL counter.
  IFAIL = 0

  ! . Find the number of variables.
  NVAR = SIZE ( X )

  ! . Skip the calculation if there are no variables.
  IF ( NVAR <= 0 ) RETURN

  ! . Initialize some algorithm options.
  NPRINT = 0
  NSTEP  = 0
  STPSIZ = 1.0_DP / SQRT( REAL(MINI_ITER, DP) )
  TOLGRD = 2.0E-2_DP

  ! . Assign the input parameters.
  IF ( PRESENT ( FORTOL          ) ) TOLGRD = FORTOL
  IF ( PRESENT ( PRINT_FREQUENCY ) ) NPRINT = PRINT_FREQUENCY
  IF ( PRESENT ( STEP_NUMBER     ) ) NSTEP  = STEP_NUMBER
  IF ( PRESENT ( STEP_SIZE       ) ) STPSIZ = STEP_SIZE / SQRT( REAL(MINI_ITER, DP) )

  ! . Check the input parameters.
  IF ( ( NPRINT < 0 ) .OR. ( NPRINT > NSTEP ) ) NPRINT = 0

  ! . Calculate the gradient tolerance.
  GTOLERANCE = REAL( N3, DP ) * TOLGRD * TOLGRD

  ! . Print out the control parameters.
  IF ( NPRINT > 0 ) THEN

    ! . Print out the options.
    CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#55FF55", VARIABLEWIDTH = 16 )
    CALL PRINT_SUMMARY_START ( "L-BFGS Minimization Calculation" )
    WRITE ( PRINT_LINE, "(I16)"    ) NVAR   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Variables"  )
    WRITE ( PRINT_LINE, "(I16)"    ) NSTEP  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations" )
    WRITE ( PRINT_LINE, "(G16.10)" ) TOLGRD ; CALL PRINT_SUMMARY_ELEMENT ( "Force Tolerance"      )
    WRITE ( PRINT_LINE, "(I16)"    ) NPRINT ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency"      )
    WRITE ( PRINT_LINE, "(G16.10)" ) STPSIZ ; CALL PRINT_SUMMARY_ELEMENT ( "Max. Step Size"       )
    CALL PRINT_SUMMARY_STOP

    ! . Print out the header for the iterations.
    CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#55FF55", PAGEWIDTH = 80, &
          & VARIABLEWIDTHS = (/ 6, 12, 15, 13, 13 /) )
    CALL PRINT_TABLE_START
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",           HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Grad. eval.",    HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Max. RMS Force", HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Path length",    HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Max. curv.",     HEADER = .TRUE. )

  END IF

  ! .  Calculate the initial energies.
  DO I = 1,NCHAIN
    II = ( I - 1 ) * N3
    CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
    IF (READDE) THEN
      CALL ENCODE_INTEGER(I,ISTRING)
      CALL DENSITY_READ ("dens"//FILE_IN//TRIM(ISTRING), PRINT=.FALSE.)
    ELSE
      IF ( SAVEDE .AND. I > 1 ) THEN
        CALL ENCODE_INTEGER(I-1,ISTRING)
        CALL DENSITY_READ ("dens"//FILE_OUT//TRIM(ISTRING), PRINT=.FALSE.)
      ELSE
        CALL DENSITY_GUESS
      END IF
    END IF
    CALL ENERGY (PRINT=.FALSE.)
    ECHAIN(I+1)=ETOTAL
    ECHAINQM(I+1)=EQM
    IF (SAVEDE) THEN
      CALL ENCODE_INTEGER(I,ISTRING)
      CALL DENSITY_WRITE ("dens"//FILE_OUT//TRIM(ISTRING), PRINT=.FALSE.)
    END IF
  END DO

  ! . Check the step number.
  IF ( NSTEP <= 0 ) RETURN

  !---------------------------------------------------------------------------
  ! . Set up the minimization.
  !---------------------------------------------------------------------------

  ENEVAL = 0
  CALL BUILD_SPLINE( )

  DO J = 1, NCHAIN
    CALL GET_GRADIENT(X, NVAR, J)
  END DO
  DO J = 1, NCHAIN
    CALL FORTOT(FORCE, J, MODF2=MODFORCE(J))
  END DO
  ! . Do some printing.
  IF ( NPRINT > 0 ) THEN
    WRITE ( PRINT_LINE, "(I6)"   ) 1
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(I12)"   ) ENEVAL
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.8)" ) AMUA2PS2_TO_KJMOL*SQRT(MAXVAL(MODFORCE)/REAL(N3,DP))
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F11.8)" ) SUM( SPLINE_LENGTH ) / SQRT( REAL(N3) )
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F11.5)" ) MAXVAL( SPLINE_CURV )
    CALL PRINT_TABLE_ELEMENT
  END IF

  !---------------------------------------------------------------------------
  ! . Perform the minimization.
  !---------------------------------------------------------------------------

  DO I = 2, NSTEP

    ! . Find the image with largest force and setup the LBFGS minimizer
    IMOVE = SUM( MAXLOC( MODFORCE ) )
    FCONV = MODFORCE(IMOVE) / REAL( MAX(10,MINI_ITER/2), DP )
    II = ( IMOVE - 1 ) * N3
    CALL LBFGS_INITIALIZE( N3, 4 )

    ! . Optimize this image until convergence
    DO K = 1, MINI_ITER
      CALL LBFGS_DATA( X(II+1:II+N3), -AMUA2PS2_TO_KJMOL*FORCE(II+1:II+N3) )
      STEP(II+1:II+N3) = LBFGS_STEP( )

      STPSCL = NORM( STEP(II+1:II+N3) )
      STPSCL = MIN( STPSIZ*MAXVAL(SPLINE_LENGTH) / STPSCL, 1.0_DP )
      STEP(II+1:II+N3) = STPSCL * STEP(II+1:II+N3)
      X(II+1:II+N3) = X(II+1:II+N3) + STEP(II+1:II+N3)

      CALL GET_GRADIENT(X, NVAR, IMOVE)
      CALL FORTOT(FORCE, IMOVE, MODF2=MODFORCE(IMOVE))

      IF ( MODFORCE(IMOVE) < FCONV) EXIT
    END DO

    CALL LBFGS_INITIALIZE( 0, 0, DELETE=.TRUE. )

    ! . Create an interpolating spline and redistribute if necessary
    CALL BUILD_SPLINE( )
    IF ( MAXVAL(SPLINE_LENGTH) > 1.5_DP * MINVAL(SPLINE_LENGTH) ) THEN
      CALL REDISTRIBUTE( )
      DO J = 1, NCHAIN
        CALL GET_GRADIENT(X, NVAR, J)
      END DO
      DO J = 1, NCHAIN
        CALL FORTOT(FORCE, J, MODF2=MODFORCE(J))
      END DO
      CALL BUILD_SPLINE( ONLY_ENERGY=.TRUE. )
    ELSE
      DO J = IMOVE-1, IMOVE+1
        IF (J < 1 .OR. J > NCHAIN) CYCLE
        CALL FORTOT(FORCE, J, MODF2=MODFORCE(J))
      END DO
    END IF

    ! . Print some data
    IF ( NPRINT > 0 ) THEN
      IF ( MOD ( I, NPRINT ) == 0 ) THEN
        WRITE ( PRINT_LINE, "(I6)"   ) I
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(I12)"   ) ENEVAL
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F15.8)" ) AMUA2PS2_TO_KJMOL*SQRT(MAXVAL(MODFORCE)/REAL(N3,DP))
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F11.8)" ) SUM( SPLINE_LENGTH ) / SQRT( REAL(N3) )
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F11.5)" ) MAXVAL( SPLINE_CURV )
        CALL PRINT_TABLE_ELEMENT

        ! . Activate the output trajectory.
        CALL DCD_INITIALIZE ( DCD_OUT )
        CALL DCD_ACTIVATE_WRITE ( FILE_OUT, DCD_OUT, "CORD", NATOMS, NFIXED, NPOINTS, &
                                  QFIX = ATMFIX, PRINT=.FALSE. )

        ! . Write out the first structure of the chain to the output trajectory.
        CALL X_TO_CRD ( XBEG, ATMCRD )
        CALL DCD_WRITE ( DCD_OUT, ATMCRD )

        ! . Write out the intermediate structures of the chain.
        DO J = 1,NCHAIN
           II = ( J - 1 ) * N3
           CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
           CALL DCD_WRITE ( DCD_OUT, ATMCRD )
        END DO

        ! . Write out the last structure of the chain.
        CALL X_TO_CRD ( XEND, ATMCRD )
        CALL DCD_WRITE ( DCD_OUT, ATMCRD )

        ! . Deactivate the trajectory.
        CALL DCD_DEACTIVATE ( DCD_OUT )
      END IF
    END IF

    IF (AMUA2PS2_TO_KJMOL**2*MAXVAL(MODFORCE) <= GTOLERANCE) EXIT
  END DO

  IF (I>NSTEP) THEN
    IFAIL=131
  ELSE
    IFAIL=0
  END IF

  !---------------------------------------------------------------------------
  ! . Finish up.
  !---------------------------------------------------------------------------

  ! . Calculate the distances between neighbouring structures.
  D1 = SPLINE_LENGTH / SQRT( REAL(N3) )

  ! . Do some more printing.
  IF ( NPRINT > 0 ) THEN

     ! . Print out the terminator.
     CALL PRINT_TABLE_STOP

     ! . Print out the status.
     SELECT CASE ( IFAIL )
     CASE (   0 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Force tolerance reached.")
     CASE ( 131 ) ; CALL PRINT_PARAGRAPH ( TEXT = "Minimization Status: Too many steps."         )
     END SELECT

  END IF

  ! . Assign a value to STATUS if necessary.
  IF ( PRESENT ( STATUS ) ) STATUS = IFAIL

  END SUBROUTINE LBFGS_MINIMIZE

  !---------------------------------------------------------------------
  SUBROUTINE PROJECT (VEL, FOR)
  !-------------------------------------------------------------------

  ! . Array arguments.
  REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: VEL
  REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: FOR

  ! . Local scalars.
  REAL ( KIND = DP ) :: VELFOR
  INTEGER            :: I

  ! . Project the full velocity on the full force.
!  VELFOR = DOT_PRODUCT (VEL, FOR)
!  IF (VELFOR > 0.0_DP) THEN
!    VEL = VELFOR * FOR / DOT_PRODUCT (FOR, FOR)
!  ELSE
!    VEL = 0.0_DP
!  END IF

  ! . Project each structure velocity on each structure force.
!  DO I = 1, SIZE(VEL), N3
!    VELFOR = DOT_PRODUCT (VEL(I:I+N3-1), FOR(I:I+N3-1))
!    IF (VELFOR > 0.0) THEN
!      VEL(I:I+N3-1) = VELFOR * FOR(I:I+N3-1) / DOT_PRODUCT (FOR(I:I+N3-1), FOR(I:I+N3-1))
!    ELSE
!      VEL(I:I+N3-1) = 0.0_DP
!    END IF
!  END DO

  ! . Project each atom velocity on each atom force.
  DO I = 1, SIZE(VEL), 3
    VELFOR = DOT_PRODUCT (VEL(I:I+2), FOR(I:I+2))
    IF (VELFOR > 0.0) THEN
      VEL(I:I+2) = VELFOR * FOR(I:I+2) / DOT_PRODUCT (FOR(I:I+2), FOR(I:I+2))
    ELSE
      VEL(I:I+2) = 0.0_DP
    END IF
  END DO

  ! . Project each velocity component on each force component.
!  WHERE (VEL*FOR < 0.0_DP) VEL=0.0_DP

  END SUBROUTINE PROJECT

  !---------------------------------------------------------------------
  SUBROUTINE TANGENT(TANG, DELTARA, DELTARB, I)
  !---------------------------------------------------------------------
  INTEGER, INTENT(IN) :: I
  REAL (KIND = DP), DIMENSION (N3), INTENT (OUT) :: TANG, DELTARA, DELTARB
  INTEGER :: II, IB, IA, J
  REAL (KIND = DP) :: VB, VX, VA

  IF ((I > 1) .AND. (I < NCHAIN+2)) THEN
    J = I-1
    IF (J==1) THEN
      DELTARB = DISTANCE ( XBEG, X(1:N3) )
      DELTARA = DISTANCE ( X(1:N3), X(N3+1:2*N3) )
    ELSE IF (J==NCHAIN) THEN
      II = ( NCHAIN - 1 ) * N3
      IB = ( NCHAIN - 2 ) * N3
      DELTARB = DISTANCE ( X(IB+1:IB+N3), X(II+1:II+N3) )
      DELTARA = DISTANCE ( X(II+1:II+N3), XEND )
    ELSE
      IB = ( J - 2 ) * N3
      II = ( J - 1 ) * N3
      IA = ( J     ) * N3
      DELTARB = DISTANCE ( X(IB+1:IB+N3), X(II+1:II+N3) )
      DELTARA = DISTANCE ( X(II+1:II+N3), X(IA+1:IA+N3) )
    END IF
    VB = ECHAIN(I - 1)
    VX = ECHAIN(I)
    VA = ECHAIN(I + 1)
    IF ((VA==VX) .AND. (VX==VB)) THEN
      TANG = DELTARA+DELTARB
    ELSE IF ((VA>=VX) .AND. (VX>=VB)) THEN
      TANG = DELTARA
    ELSE IF ((VB>=VX) .AND. (VX>=VA)) THEN
      TANG = DELTARB
    ELSE IF ((VA<VX) .AND. (VX>VB)) THEN
      TANG = DELTARA*(VX-VB)+DELTARB*(VX-VA)
    ELSE IF ((VA>VX) .AND. (VX<VB)) THEN
      TANG = DELTARA*(VA-VX)+DELTARB*(VB-VX)
    ELSE
      CALL PRINT_ERROR ("TANGENT", "Tangent error. Probably due to NAN.")
      STOP
    END IF
  ELSE
    IF (I==1)        TANG = X(1:N3)-XBEG
    IF (I==NCHAIN+2) TANG = XEND - X((NCHAIN-1)*N3+1:NCHAIN*N3)
  END IF

  CALL CONSTRAIN ( TANG )

  TANG = TANG/NORM(TANG)

  END SUBROUTINE TANGENT

  !---------------------------------------------------------------------
  FUNCTION PARGRAD(I)
  !---------------------------------------------------------------------

  ! . Scalar arguments
  INTEGER, INTENT(IN):: I

  ! . Function value.
  REAL (KIND = DP) :: PARGRAD

  ! . Local scalars.
  INTEGER :: II

  ! . Local arrays.
  REAL (KIND = DP), DIMENSION (N3) :: TANG

  ! . Calculate the tangent to each image.
  CALL SPLINE_PROP( I, TANGENT=TANG )

  ! . Calculate parallel force
  IF ( I == 1 ) THEN
    PARGRAD = DOT_PRODUCT(GBEG, TANG)
  ELSE IF ( I == NCHAIN+2 ) THEN
    PARGRAD = DOT_PRODUCT(GEND, TANG)
  ELSE
    II = ( I - 2 ) * N3
    PARGRAD = DOT_PRODUCT(G(II+1:II+N3), TANG)
  END IF

  END FUNCTION PARGRAD

  !---------------------------------------------------------------------
  SUBROUTINE FLUX(PARTIAL)
  !---------------------------------------------------------------------

  REAL( KIND = DP ), DIMENSION(NCHAIN+1), INTENT(OUT):: PARTIAL

  REAL( KIND = DP ), DIMENSION(8) :: AX, AW
  REAL( KIND = DP )               :: S, F, ENERG, EBASE
  INTEGER                         :: I, J

  ! . Setup constants for Gaussian quadrature
  AX(1) = 0.183434642495650_DP; AX(8) = -AX(1)
  AX(2) = 0.525532409916329_DP; AX(7) = -AX(2)
  AX(3) = 0.796666477413627_DP; AX(6) = -AX(3)
  AX(4) = 0.960289856497536_DP; AX(5) = -AX(4)
  AW(1) = 0.362683783378362_DP; AW(8) =  AW(1)
  AW(2) = 0.313706645877887_DP; AW(7) =  AW(2)
  AW(3) = 0.222381034453374_DP; AW(6) =  AW(3)
  AW(4) = 0.101228536290376_DP; AW(5) =  AW(4)

  EBASE = MIN( ECHAIN(1), ECHAIN(NCHAIN+2) )

  ! . Integrate the flux in this spline segment
  DO I = 1, NCHAIN+1
    PARTIAL(I) = 0.0_DP
    F = 0.5_DP * SPLINE_LENGTH(I)
    DO J = 1, SIZE(AX)
      S = 0.5_DP * ( AX(J) + 1.0_DP )
      ENERG = ENERGY_SPLINE(1,I)       + ENERGY_SPLINE(2,I) * S     + &
              ENERGY_SPLINE(3,I) * S*S + ENERGY_SPLINE(4,I) * S*S*S
      PARTIAL(I) = PARTIAL(I) + F * AW(J) * EXP( BETA * ( ENERG - EBASE ) )
    END DO
  END DO

  END SUBROUTINE FLUX

  !---------------------------------------------------------------------
  SUBROUTINE INTERPOLATE (NBININ)
  !---------------------------------------------------------------------

  ! . Input scalars.
  INTEGER,            INTENT(IN),  OPTIONAL :: NBININ

  ! . Local arrays.
  REAL ( KIND = DP ), DIMENSION(NCHAIN+1,N3) :: COORD

  ! . Local scalars.
  INTEGER :: NBINS, IREF, NMAX
  REAL (KIND = DP) :: STEP, POS, RELPOS, ESTIME, SUMREF, FACTOR
  LOGICAL :: FOUND

  ! . Local types.
  TYPE ( DCD_TYPE ) :: DCD_OUT

  ! .  Check if there is an input and if its value is reasonable.
  NBINS = 0
  IF ( PRESENT ( NBININ  ) ) THEN
   NBINS = NBININ
  ELSE
   NBINS = 0
  END IF
  IF (NBINS==0) THEN
    CALL PRINT_PARAGRAPH (Text="No interpolation requested.")
  ELSE IF (NBINS<NPOINTS) THEN
    CALL PRINT_PARAGRAPH (Text="NBINS should be greater than the number of chain structures.")
  ELSE

    ! . Write out the header.
    CALL PRINT_TABLE_OPTIONS ( COLUMNS = 4, HEADER_COLOR = "#8888FF", VARIABLEWIDTHS = (/ 10, 15, 15, 15 /) )
    CALL PRINT_TABLE_START
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Polinomial Interpolation", COLSPAN = 4, HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Structure",                HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Distance",                 HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Interpolated points",      HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Original points",          HEADER = .TRUE. )

    ! . Calculate the distance step.
    STEP = SUM(SPLINE_LENGTH) / (NBINS-1)
    FACTOR = 1.0_DP / SQRT( REAL(N3) )

    IF (TRIM(INTERPFILE) /= '') OPEN(78,FILE=TRIM(INTERPFILE))
    IREF = 1
    SUMREF = 0.0_DP
    ESTIME = ENERGY_SPLINE(1,IREF)
    WRITE ( PRINT_LINE, "(I10)"   ) IREF            ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) SUMREF * FACTOR ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) ESTIME          ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) ECHAIN(IREF)    ; CALL PRINT_TABLE_ELEMENT
    IF (TRIM(INTERPFILE) /= '') WRITE(78,"(I10,3F15.5)") IREF, SUMREF*FACTOR, ESTIME, ECHAIN(IREF)

    ! . Loop over the structures.
    DO I = 1, NBINS-2
      ! . Calculate the coordinate and the pair of images it is limited by.
      POS = I * STEP
      DO WHILE ( SUMREF + SPLINE_LENGTH(IREF) < POS )
        SUMREF = SUMREF + SPLINE_LENGTH(IREF)
        IREF = IREF + 1
        ESTIME = ENERGY_SPLINE(1,IREF)
        WRITE ( PRINT_LINE, "(I10)"   ) IREF            ; CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F15.5)" ) SUMREF * FACTOR ; CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F15.5)" ) ESTIME          ; CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F15.5)" ) ECHAIN(IREF)    ; CALL PRINT_TABLE_ELEMENT
        IF (TRIM(INTERPFILE) /= '') WRITE(78,"(I10,3F15.5)") IREF, SUMREF*FACTOR, ESTIME, ECHAIN(IREF)
      END DO

      RELPOS = ( POS - SUMREF ) / SPLINE_LENGTH(IREF)

      ESTIME = ENERGY_SPLINE(1,IREF) + &
               ENERGY_SPLINE(2,IREF) * RELPOS + &
               ENERGY_SPLINE(3,IREF) * RELPOS**2 + &
               ENERGY_SPLINE(4,IREF) * RELPOS**3
      WRITE ( PRINT_LINE, "(A10)"   ) "*"          ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F15.5)" ) POS * FACTOR ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F15.5)" ) ESTIME       ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(A15)"   ) ""           ; CALL PRINT_TABLE_ELEMENT
      IF (TRIM(INTERPFILE) /= '') WRITE(78,"(A10,2F15.5)") "*", POS*FACTOR, ESTIME
    END DO

    SUMREF = SUM(SPLINE_LENGTH(1:IREF))
    ESTIME = ENERGY_SPLINE(1,IREF) + ENERGY_SPLINE(2,IREF) + &
             ENERGY_SPLINE(3,IREF) + ENERGY_SPLINE(4,IREF)
    WRITE ( PRINT_LINE, "(I10)"   ) IREF+1          ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) SUMREF * FACTOR ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) ESTIME          ; CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.5)" ) ECHAIN(IREF+1)  ; CALL PRINT_TABLE_ELEMENT
    IF (TRIM(INTERPFILE) /= '') WRITE(78,"(I10,3F15.5)") IREF+1, SUMREF*FACTOR, ESTIME, ECHAIN(IREF+1)

    CALL PRINT_TABLE_STOP
    IF (TRIM(INTERPFILE) /= '') CLOSE(78)

  END IF

  ! . Find the maxima in the energy interpolation
  CALL PRINT_TABLE_OPTIONS ( COLUMNS = 2, HEADER_COLOR = "#8888FF", PAGEWIDTH = 60, VARIABLEWIDTHS = (/ 15, 15 /) )
  CALL PRINT_TABLE_START
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Maxima found", COLSPAN = 2, HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Distance",                  HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Interpolated energy",       HEADER = .TRUE. )
  FACTOR = 1.0_DP / SQRT( REAL(N3) )
  NMAX = 0
  DO I = 1, NCHAIN+1
    CALL FIND_MAX( ENERGY_SPLINE(:,I), FOUND, POS, ESTIME )
    IF ( FOUND ) THEN
      NMAX = NMAX+1
      COORD(NMAX,:) = SPLINE(1,I,:)           + SPLINE(2,I,:) * POS + &
                      SPLINE(3,I,:) * POS*POS + SPLINE(4,I,:) * POS*POS*POS
      POS = SUM( SPLINE_LENGTH(1:I-1) ) + POS * SPLINE_LENGTH(I)
      WRITE ( PRINT_LINE, "(F15.5)" ) POS * FACTOR ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F15.5)" ) ESTIME       ; CALL PRINT_TABLE_ELEMENT
    END IF
  END DO
  CALL PRINT_TABLE_STOP

  ! . Store the maxima structures if asked for
  IF ( TRIM(MAXFILE) /= '' ) THEN
    CALL DCD_INITIALIZE ( DCD_OUT )
    CALL DCD_ACTIVATE_WRITE ( MAXFILE, DCD_OUT, "CORD", NATOMS, NFIXED, NMAX, QFIX = ATMFIX )
    DO I = 1, NMAX
      CALL X_TO_CRD ( COORD(I,:), ATMCRD )
      CALL DCD_WRITE ( DCD_OUT, ATMCRD )
    END DO
    CALL DCD_DEACTIVATE ( DCD_OUT )
  END IF

  END SUBROUTINE INTERPOLATE

  !---------------------------------------------------------------------
  SUBROUTINE GET_GRADIENT(X, N, ISTRUCT)
  !---------------------------------------------------------------------

  INTEGER, INTENT(IN)                          :: N, ISTRUCT
  REAL ( KIND = DP ), DIMENSION(N), INTENT(IN) :: X

  INTEGER :: II

  II = ( ISTRUCT - 1 ) * N3

  CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )

  IF (SAVEDE) THEN
    CALL ENCODE_INTEGER(ISTRUCT,ISTRING)
    CALL DENSITY_READ ("dens"//FILE_OUT//TRIM(ISTRING), PRINT=.FALSE.)
  ELSE
    CALL DENSITY_GUESS
  END IF
  CALL GRADIENT ( PRINT = .FALSE. )
  IF (SAVEDE) THEN
    CALL ENCODE_INTEGER(ISTRUCT,ISTRING)
    CALL DENSITY_WRITE ("dens"//FILE_OUT//TRIM(ISTRING), PRINT=.FALSE.)
  END IF

  ECHAIN(ISTRUCT+1)   = ETOTAL
  ECHAINQM(ISTRUCT+1) = EQM

  CALL CRD_TO_X ( ATMDER, G(II+1:II+N3) )

  ENEVAL = ENEVAL + 1

  END SUBROUTINE GET_GRADIENT

  !---------------------------------------------------------------------
  SUBROUTINE FORTOT(FORCE, ISTRUCT, MODF2)
  !---------------------------------------------------------------------

  ! . Scalar arguments
  INTEGER, INTENT(IN) :: ISTRUCT
  REAL ( KIND = DP ), INTENT(OUT), OPTIONAL :: MODF2

  ! . Array arguments
  REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: FORCE

  ! . Local arrays
  REAL (KIND = DP), DIMENSION(N3) :: DELTARB, DELTARA, TANG, GRAD, NORMVEC

  ! . Local scalars.
  REAL ( KIND = DP) :: VA, VX, VB, BETAU
  INTEGER :: II

  IF ( ( ISTRUCT < 1 ) .OR. ( ISTRUCT > NCHAIN ) ) RETURN

  IF ( ( ISTRUCT == 1 ) .OR. ( ISTRUCT == NCHAIN ) ) THEN
    II = ( ISTRUCT - 1 ) * N3
    FORCE(II+1:II+N3) = 0.0_DP
    IF ( PRESENT(MODF2) ) MODF2 = 0.0_DP
  END IF

  VB = ECHAIN(ISTRUCT)
  VX = ECHAIN(ISTRUCT+1)
  VA = ECHAIN(ISTRUCT+2)

  ! . Calculate the index into the coordinate array.
  II = ( ISTRUCT - 1 ) * N3

  ! . Calculate the tangent to each image.
  CALL TANGENT ( TANG, DELTARA, DELTARB, ISTRUCT+1 )
  !CALL SPLINE_PROP ( ISTRUCT+1, TANGENT=TANG )

  ! . Calculate the gradient.
  SELECT CASE (METH)
  CASE ( 0 )
    ! . Temperature = 0 K
    GRAD = G(II+1:II+N3)
  CASE ( 1 )
    ! . Differential method
    NORMVEC = DELTARA - DELTARB
    !CALL SPLINE_PROP( ISTRUCT+1, NORMAL=NORMVEC )
    NORMVEC = NORMVEC - DOT_PRODUCT(NORMVEC, TANG) * TANG
    ! . Set the normal to zero in the linear case
    IF ( NORM(NORMVEC) > 1.0E-7_DP ) THEN
      NORMVEC = NORMVEC/NORM(NORMVEC)
    ELSE
      NORMVEC = 0.0_DP
    END IF
    GRAD = G(II+1:II+N3) - SPLINE_CURV(ISTRUCT+1)/BETA * NORMVEC
  CASE ( 2 )
    ! . Integral method
    GRAD = ( (NORM(DELTARB)+NORM(DELTARA)) * BETA * G(II+1:II+N3) + &
             DELTARB/NORM(DELTARB) * ( 1.0_DP + EXP(BETA*(VB-VX)) ) - &
             DELTARA/NORM(DELTARA) * ( 1.0_DP + EXP(BETA*(VA-VX)) ) ) * &
           0.5_DP * EXP(BETA*VX)
  CASE ( 3 )
    ! . Integral method #2 (old)
    ! (where does this come from?)
    BETAU = BETA * ( VA - VB )
    IF ( ABS(BETAU) < 40.0_DP ) THEN
      GRAD = G(II+1:II+N3) + 2.0_DP/BETA * &
             ( DELTARB/NORM(DELTARB) - EXP(BETAU) * DELTARA/NORM(DELTARA) ) / &
             ( NORM(DELTARB) - EXP(BETAU) * NORM(DELTARA) )
    ELSE IF ( BETAU > 0.0_DP ) THEN
      GRAD = G(II+1:II+N3) - 2.0_DP/BETA * DELTARA/NORM(DELTARA)**2
    ELSE IF ( BETAU < 0.0_DP ) THEN
      GRAD = G(II+1:II+N3) + 2.0_DP/BETA * DELTARB/NORM(DELTARB)**2
    END IF
  END SELECT

  ! . Project out the rigid body motions from the force.
  CALL CONSTRAIN ( GRAD )

  ! . Calculate (perpendicular) force and change to atomic units.
  FORCE (II+1:II+N3) = -KJMOL_TO_AMUA2PS2 * ( GRAD - DOT_PRODUCT(GRAD, TANG) * TANG )

  IF ( PRESENT(MODF2) ) THEN
    MODF2 = DOT_PRODUCT( FORCE(II+1:II+N3), FORCE(II+1:II+N3) )
  END IF

  END SUBROUTINE FORTOT

  !---------------------------------------------------------------------
  SUBROUTINE BUILD_SPLINE ( ONLY_ENERGY )
  !---------------------------------------------------------------------

  ! . Scalar arguments.
  LOGICAL, INTENT(IN), OPTIONAL               :: ONLY_ENERGY

  ! . Local arrays.
  REAL( KIND = DP ), DIMENSION((NCHAIN+2)*N3) :: COORDS

  ! . Local scalars.
  REAL( KIND = DP )                           :: FA,FB
  INTEGER                                     :: I, NP
  LOGICAL                                     :: QOE

  ! . Setup initial values.
  NP = NCHAIN+2
  QOE = .FALSE.
  IF ( PRESENT(ONLY_ENERGY) ) QOE = ONLY_ENERGY

  ! . Calculate the spline if needed.
  IF ( .NOT. QOE ) THEN

    COORDS(1:N3)                          = XBEG
    COORDS(N3+1:(NCHAIN+1)*N3)            = X
    COORDS((NCHAIN+1)*N3+1:(NCHAIN+2)*N3) = XEND

    ! . Calculate a interpolating spline.
    CALL CUBIC_SPLINE ( COORDS )

    ! . Calculate segment lengths and curvatures.
    DO I = 1, NP-1
      CALL INTEGRATE_SPLINE( I, 0.0_DP, 1.0_DP, SPLINE_LENGTH(I) )
      CALL SPLINE_PROP( I, CURVATURE=SPLINE_CURV(I) )
    END DO
    CALL SPLINE_PROP( NP, CURVATURE=SPLINE_CURV(NP) )

  END IF

  ! . Calculate a spline for the energy profile.
  FB = PARGRAD(1)
  DO I = 1, NP-1
    FA = FB
    FB = PARGRAD(I+1)
    ENERGY_SPLINE(1,I) = ECHAIN(I)
    ENERGY_SPLINE(2,I) = FA * SPLINE_LENGTH(I)
    ENERGY_SPLINE(3,I) = 3.0_DP * ( ECHAIN(I+1) - ECHAIN(I) ) - ( 2.0_DP * FA + FB ) * SPLINE_LENGTH(I)
    ENERGY_SPLINE(4,I) = 2.0_DP * ( ECHAIN(I) - ECHAIN(I+1) ) + ( FA + FB ) * SPLINE_LENGTH(I)
  END DO

  END SUBROUTINE BUILD_SPLINE

  !---------------------------------------------------------------------
  SUBROUTINE REDISTRIBUTE ( )
  !---------------------------------------------------------------------

  ! . Local arrays
  REAL( KIND = DP ), DIMENSION(0:NCHAIN+1) :: POS

  ! . Local scalars
  REAL( KIND = DP )                        :: TOTAL_LENGTH, LENGTH, S
  INTEGER                                  :: I, II, IND

  TOTAL_LENGTH = SUM( SPLINE_LENGTH )

  ! . Get the desired locatin in the spline for each image.
  DO I = 0, NCHAIN+1
    POS(I) = I * TOTAL_LENGTH / REAL( NCHAIN+1, DP )
  END DO

  ! . Calculate a interpolated structure for each image.
  TOTAL_LENGTH = 0.0_DP
  IND = 1
  DO I = 1, NCHAIN
    DO WHILE ( TOTAL_LENGTH + SPLINE_LENGTH(IND) <= POS(I) )
      TOTAL_LENGTH = TOTAL_LENGTH + SPLINE_LENGTH(IND)
      IND = IND + 1
    END DO
    LENGTH = POS(I) - TOTAL_LENGTH

    CALL FIND_POSITION ( IND, SPLINE_LENGTH(IND), LENGTH, S )

    II = ( I - 1 ) * N3
    X(II+1:II+N3) = SPLINE(1,IND,:)       + SPLINE(2,IND,:) * S + &
                    SPLINE(3,IND,:) * S*S + SPLINE(4,IND,:) * S*S*S
  END DO

  CALL BUILD_SPLINE( )

  END SUBROUTINE REDISTRIBUTE

END SUBROUTINE NEB_OPTIMIZE_SPLINE

!---------------------------
FUNCTION DISTANCE ( X1, X2 )
!---------------------------

! . Array arguments.
REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: X1, X2

! . Function declarations.
REAL ( KIND = DP ), DIMENSION (SIZE(X1)) :: DISTANCE

! . Calculate the distance.
DISTANCE = ( X2 - X1 )

END FUNCTION DISTANCE

!-----------------------------
SUBROUTINE CRD_TO_X ( CRD, X )
!-----------------------------

! . Array arguments.
REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(IN)  :: CRD
REAL ( KIND = DP ), DIMENSION(1:N3),         INTENT(OUT) :: X

! . Local scalars.
INTEGER :: IATOM, II

! . Copy coordinates for free atoms only.
II = -3
DO IATOM = 1,NATOMS
  IF ( .NOT. ATMFIX(IATOM) ) THEN
    II = II + 3
    X(II+1:II+3) = CRD(1:3,IATOM)
  END IF
END DO

END SUBROUTINE CRD_TO_X

!-----------------------------
SUBROUTINE X_TO_CRD ( X, CRD )
!-----------------------------

! . Array arguments.
REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(OUT) :: CRD
REAL ( KIND = DP ), DIMENSION(1:N3),         INTENT(IN)  :: X

! . Local scalars.
INTEGER :: IATOM, II

! . Copy coordinates for free atoms only.
II = -3
DO IATOM = 1,NATOMS
  IF ( .NOT. ATMFIX(IATOM) ) THEN
    II = II + 3
    CRD(1:3,IATOM) = X(II+1:II+3)
  END IF
END DO

END SUBROUTINE X_TO_CRD

!---------------------------------------------------------
SUBROUTINE CUBIC_SPLINE ( COORD )
!---------------------------------------------------------

! . Array arguments
REAL ( KIND = DP ), DIMENSION(N3*(NCHAIN+2)), INTENT(IN) :: COORD

! . Local scalars
INTEGER :: I, II, JJ, NPT

! . Local arrays
REAL ( KIND = DP ) :: GAM(NCHAIN+2), DERIV(NCHAIN+2,N3)

NPT = NCHAIN+2

! . Setup and solve the equation system to get the coefficients
GAM(1) = 0.5_DP
DO I = 2, NPT-1
  GAM(I) = 1.0_DP / ( 4.0_DP - GAM(I-1) )
END DO
GAM(NPT) = 1.0_DP / ( 2.0_DP - GAM(I-1) )

DERIV(1,:) = 3.0_DP * ( COORD(N3+1:2*N3) - COORD(1:N3) )
DO I = 2, NPT-1
  II = ( I     ) * N3
  JJ = ( I - 2 ) * N3
  DERIV(I,:) = 3.0_DP * ( COORD(II+1:II+N3) - COORD(JJ+1:JJ+N3) )
END DO
JJ = ( NPT - 2 ) * N3
DERIV(NPT,:) = 3.0_DP * ( COORD(JJ+N3+1:JJ+2*N3) - COORD(JJ+1:JJ+N3) )

DERIV(1,:) = GAM(1) * DERIV(1,:)
DO I = 2, NPT
  DERIV(I,:) = GAM(I) * (DERIV(I,:) - DERIV(I-1,:) )
END DO
DO I = NPT-1, 1, -1
  DERIV(I,:) = DERIV(I,:) - GAM(I)*DERIV(I+1,:)
END DO

! . Calculate the spline coefficients
DO I = 1, NPT-1
  II = ( I - 1 ) * N3
  JJ = ( I     ) * N3
  SPLINE(1,I,:) = COORD(II+1:II+N3)
  SPLINE(2,I,:) = DERIV(I,:)
  SPLINE(3,I,:) = 3.0_DP * ( COORD(JJ+1:JJ+N3) - COORD(II+1:II+N3) ) - 2.0_DP * DERIV(I,:) - DERIV(I+1,:)
  SPLINE(4,I,:) = 2.0_DP * ( COORD(II+1:II+N3) - COORD(JJ+1:JJ+N3) ) +  DERIV(I,:) + DERIV(I+1,:)
END DO

END SUBROUTINE CUBIC_SPLINE

!----------------------------------------------------------------
SUBROUTINE INTEGRATE_SPLINE ( SEG, S_INI, S_FIN, LENGTH )
!----------------------------------------------------------------

! . Scalar arguments
INTEGER, INTENT(IN)             :: SEG
REAL ( KIND = DP ), INTENT(IN)  :: S_INI, S_FIN
REAL ( KIND = DP ), INTENT(OUT) :: LENGTH

! . Local scalars
INTEGER :: I
REAL ( KIND = DP ) :: F, S, SPD

! . Local arrays
REAL( KIND = DP ), DIMENSION(8) :: AX, AW

! . Set up constants for Gaussian quadrature
AX(1) = 0.183434642495650_DP; AX(8) = -AX(1)
AX(2) = 0.525532409916329_DP; AX(7) = -AX(2)
AX(3) = 0.796666477413627_DP; AX(6) = -AX(3)
AX(4) = 0.960289856497536_DP; AX(5) = -AX(4)
AW(1) = 0.362683783378362_DP; AW(8) =  AW(1)
AW(2) = 0.313706645877887_DP; AW(7) =  AW(2)
AW(3) = 0.222381034453374_DP; AW(6) =  AW(3)
AW(4) = 0.101228536290376_DP; AW(5) =  AW(4)

! . Integrate the "speed" between S_INI and S_FIN
F = 0.5_DP * ( S_FIN - S_INI )
LENGTH = 0.0_DP
DO I = 1, SIZE(AX)
  S = F * AX(I) + 0.5_DP * ( S_INI + S_FIN )
  CALL SPLINE_PROP( SEG, S, SPEED=SPD )
  LENGTH = LENGTH + F * AW(I) * SPD
END DO

END SUBROUTINE INTEGRATE_SPLINE

!---------------------------------------------------
SUBROUTINE FIND_POSITION (  SEG, S_LEN, DIST, POS )
!---------------------------------------------------

! . Scalar arguments
INTEGER, INTENT(IN)             :: SEG
REAL ( KIND = DP ), INTENT(IN)  :: S_LEN, DIST
REAL ( KIND = DP ), INTENT(OUT) :: POS

! . Local scalars
REAL ( KIND = DP ) :: A, B, LA, LB, LENGTH

! . Set up initial values
A = 0.0_DP ; LA = -DIST
B = 1.0_DP ; LB = S_LEN - DIST
POS = B

! . Find the root with the secant method
DO WHILE ( ABS( LB/S_LEN ) > 1.0E-4_DP )
  POS = B -  LB * ( B - A ) / ( LB - LA )
  CALL INTEGRATE_SPLINE( SEG, 0.0_DP, POS, LENGTH )
  A = B   ; LA = LB
  B = POS ; LB = LENGTH - DIST
END DO

END SUBROUTINE FIND_POSITION

!-------------------------------------------------------------------
SUBROUTINE SPLINE_PROP ( SEG, S, SPEED, TANGENT, NORMAL, CURVATURE )
!-------------------------------------------------------------------

! . Scalar arguments
INTEGER, INTENT(IN) :: SEG
REAL ( KIND = DP), INTENT(IN), OPTIONAL  :: S
REAL ( KIND = DP), INTENT(OUT), OPTIONAL :: SPEED, CURVATURE

! . Array arguments
REAL ( KIND = DP), DIMENSION(N3), INTENT(OUT), OPTIONAL :: TANGENT, NORMAL

! . Local scalars
REAL ( KIND = DP) :: INVSPEED, DERS, SK

! . Local arrays
REAL ( KIND = DP), DIMENSION(N3) :: DER1, DER2, VT, VN

! . Calculate the first derivative, the speed and the tangent
IF ( PRESENT(S) ) THEN
  DER1 = SPLINE(2,SEG,:) + 2.0_DP * SPLINE(3,SEG,:) * S + 3.0_DP * SPLINE(4,SEG,:) * S*S
ELSE
  IF ( SEG < NCHAIN+2 ) THEN
    DER1 = SPLINE(2,SEG,:)
  ELSE
    DER1 = SPLINE(2,SEG-1,:) + 2.0_DP * SPLINE(3,SEG-1,:) + 3.0_DP * SPLINE(4,SEG-1,:)
  END IF
END IF

IF ( PRESENT(SPEED) ) SPEED = NORM(DER1)
INVSPEED = 1.0_DP/NORM(DER1)
VT = INVSPEED * DER1
IF ( PRESENT(TANGENT) ) TANGENT = VT

! . Calculate the second derivative, the curvature and the normal
IF ( PRESENT(NORMAL) .OR. PRESENT(CURVATURE) ) THEN

  IF ( PRESENT(S) ) THEN
    DER2 = 2.0_DP * SPLINE(3,SEG,:) + 6.0_DP * SPLINE(4,SEG,:) * S
    DERS = 2.0_DP * SUM(        SPLINE(2,SEG,:)*SPLINE(3,SEG,:) + &
                         2.0_DP*SPLINE(3,SEG,:)**2 * S + &
                         3.0_DP*SPLINE(2,SEG,:)*SPLINE(4,SEG,:) * S + &
                         9.0_DP*SPLINE(3,SEG,:)*SPLINE(4,SEG,:) * S*S + &
                         9.0_DP*SPLINE(4,SEG,:)**2 * S*S*S ) * INVSPEED
  ELSE
    IF ( SEG < NCHAIN+2 ) THEN
      DER2 = 2.0_DP * SPLINE(3,SEG,:)
      DERS = 2.0_DP * SUM( SPLINE(2,SEG,:)*SPLINE(3,SEG,:) ) * INVSPEED
    ELSE
      DER2 = 2.0_DP * SPLINE(3,SEG-1,:) + 6.0_DP * SPLINE(4,SEG-1,:)
      DERS = 2.0_DP * SUM(        SPLINE(2,SEG-1,:)*SPLINE(3,SEG-1,:) + &
                           2.0_DP*SPLINE(3,SEG-1,:)**2 + &
                           3.0_DP*SPLINE(2,SEG-1,:)*SPLINE(4,SEG-1,:) + &
                           9.0_DP*SPLINE(3,SEG-1,:)*SPLINE(4,SEG-1,:) + &
                           9.0_DP*SPLINE(4,SEG-1,:)**2 ) * INVSPEED
    END IF
  END IF

  VN = INVSPEED * ( DER2 - DERS*VT )
  SK = NORM(VN)
  IF ( PRESENT(NORMAL) ) THEN
    IF ( SK > 1.0E-7_DP ) THEN
      NORMAL = VN / SK
    ELSE
      NORMAL = 0.0_DP
    END IF
  END IF
  IF ( PRESENT(CURVATURE) ) CURVATURE = INVSPEED * SK

END IF

END SUBROUTINE SPLINE_PROP

!---------------------------------------
SUBROUTINE FIND_MAX( COEF, FOUND, X, Y )
!---------------------------------------

! . Array arguments
REAL ( KIND = DP ), DIMENSION(4), INTENT(IN) :: COEF

! . Scalar arguments
LOGICAL, INTENT(OUT)            :: FOUND
REAL ( KIND = DP ), INTENT(OUT) :: X, Y

! . Local scalars
REAL ( KIND = DP ) :: C, D, DISCR

X = 0.0_DP
FOUND = .FALSE.
C = COEF(3)
D = 3.0_DP*COEF(4)

! . Return if the polynomial has no extrema
DISCR = C*C - D*COEF(2)
IF ( DISCR <= 0.0_DP ) RETURN

! . Calculate the maximum
X = ( -C - SQRT(DISCR) ) / D

! . Calculate the energy if the root is not off-limits
IF ( X > 0.0_DP .AND. X < 1.0_DP ) THEN
  Y = COEF(1) + COEF(2)*X + COEF(3)*X*X + COEF(4)*X*X*X
  FOUND = .TRUE.
END IF

END SUBROUTINE FIND_MAX

!-----------------------------------------------------------------------------------
SUBROUTINE GROWING_STRING( FILE_OUT, NIMG, COORD_INI, COORD_FIN, STEP_FACTOR )
!-----------------------------------------------------------------------------------

CHARACTER ( LEN = * ), INTENT(IN) :: FILE_OUT
INTEGER, INTENT(IN)               :: NIMG
REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: COORD_INI, COORD_FIN
REAL ( KIND = DP ), OPTIONAL, INTENT(IN) :: STEP_FACTOR

REAL ( KIND = DP ), DIMENSION(:), ALLOCATABLE :: XINI, XFIN, X, G, TANG, STEP
REAL ( KIND = DP ) :: GTOLERANCE, STPSIZ, FCT
REAL ( KIND = DP ), PARAMETER :: SAME = 1.0E-10_DP
CHARACTER (LEN = 3 ) :: ISTRING
LOGICAL :: SAVEDE
INTEGER :: I, J
TYPE ( DCD_TYPE ) :: DCD_OUT

FCT = 1.0_DP
IF ( PRESENT(STEP_FACTOR) ) FCT = STEP_FACTOR

IF ( NATOMSQM > 0 ) THEN
  SAVEDE = .TRUE.
ELSE
  SAVEDE = .FALSE.
END IF

N3 = 3*NFREE
ALLOCATE ( XINI(N3), XFIN(N3), X(N3), G(N3), TANG(N3), STEP(N3) )

IF ( NFIXED > 0 ) THEN
  DO I = 1,NATOMS
    IF ( ATMFIX(I) ) THEN
      IF ( MAXVAL ( ABS ( COORD_INI(1:3,I) - COORD_FIN(1:3,I) ) ) > SAME ) THEN
        CALL PRINT_ERROR ( "GROWING_STRING", "Fixed atom coordinate mismatch in reactant and product structures.", I )
      END IF
    END IF
  END DO
ELSE
  CALL SUPERIMPOSE_QUATERNION ( COORD_INI, COORD_FIN, PRINT = .FALSE. )
END IF

GTOLERANCE = REAL(N3) * ( 0.1_DP / AMUA2PS2_TO_KJMOL )**2
STPSIZ = SQRT( REAL(N3) ) * 0.01_DP * FCT

CALL DCD_INITIALIZE ( DCD_OUT )
CALL DCD_ACTIVATE_WRITE ( FILE_OUT, DCD_OUT, "CORD", NATOMS, NFIXED, NIMG, QFIX = ATMFIX )

CALL CRD_TO_X ( COORD_INI, X )
CALL CRD_TO_X ( COORD_FIN, XFIN )

CALL DCD_WRITE ( DCD_OUT, COORD_INI )
ATMCRD = COORD_INI
CALL DENSITY_GUESS
CALL GRADIENT ( PRINT = .TRUE. )
IF ( SAVEDE ) THEN
  CALL ENCODE_INTEGER(NIMG-1,ISTRING)
  CALL DENSITY_WRITE ("dens"//TRIM(FILE_OUT)//TRIM(ISTRING), PRINT=.TRUE.)
END IF

DO I = NIMG-1, 2, -1
  XINI = X
  TANG = XFIN-XINI
  X = XINI + 1.0_DP/REAL(I,DP) * TANG
  TANG = TANG / NORM(TANG)

  CALL LBFGS_INITIALIZE( N3, 4 )

  ! . Optimize the current image
  CALL X_TO_CRD ( X, ATMCRD )
  CALL DENSITY_GUESS
  CALL GRADIENT ( PRINT = .FALSE. )
  CALL CRD_TO_X ( ATMDER, G )
  G = G - DOT_PRODUCT(G, TANG) * TANG
  CALL LBFGS_INITIALIZE( N3, 4 )
  DO J = 1, 1000
    CALL LBFGS_DATA( X, -AMUA2PS2_TO_KJMOL*G )
    STEP = LBFGS_STEP( )
    STEP = MIN( STPSIZ/NORM(STEP), 1.0_DP ) * STEP
    X = X + STEP
    CALL X_TO_CRD ( X, ATMCRD )

    IF ( SAVEDE .AND. I < NIMG-1 ) THEN
      CALL ENCODE_INTEGER(I,ISTRING)
      CALL DENSITY_READ ("dens"//TRIM(FILE_OUT)//TRIM(ISTRING), PRINT=.FALSE.)
    ELSE
      CALL DENSITY_GUESS
    END IF
    CALL GRADIENT ( PRINT = .FALSE. )
    IF ( SAVEDE ) THEN
      CALL ENCODE_INTEGER(I-1,ISTRING)
      CALL DENSITY_WRITE ("dens"//TRIM(FILE_OUT)//TRIM(ISTRING), PRINT=.FALSE.)
    END IF

    CALL CRD_TO_X ( ATMDER, G )
    G = G - DOT_PRODUCT(G, TANG) * TANG
    IF ( DOT_PRODUCT(G,G) < GTOLERANCE ) EXIT
  END DO
  CALL LBFGS_INITIALIZE( 0, 0, DELETE=.TRUE. )

  CALL X_TO_CRD ( X, ATMCRD )
  CALL DCD_WRITE ( DCD_OUT, ATMCRD )
END DO

CALL DCD_WRITE ( DCD_OUT, COORD_FIN )
ATMCRD = COORD_FIN
CALL DENSITY_GUESS
CALL GRADIENT ( PRINT = .TRUE. )
IF ( SAVEDE ) THEN
  CALL ENCODE_INTEGER(0,ISTRING)
  CALL DENSITY_WRITE ("dens"//TRIM(FILE_OUT)//TRIM(ISTRING), PRINT=.TRUE.)
END IF

CALL DCD_DEACTIVATE ( DCD_OUT )

DEALLOCATE ( XINI, XFIN, X, G, TANG, STEP )

END SUBROUTINE GROWING_STRING

!-----------------------------------------------------------------------------------
SUBROUTINE CHAIN_INTERPOLATE( FILE_IN, FILE_OUT, NIMG, IM_INI, IM_FIN, KEEP_IMAGES )
!-----------------------------------------------------------------------------------

! . Scalar arguments
CHARACTER ( LEN = * ), INTENT(IN) :: FILE_IN, FILE_OUT
INTEGER, INTENT(IN)               :: NIMG
INTEGER, INTENT(IN), OPTIONAL     :: IM_INI, IM_FIN
LOGICAL, INTENT(IN), OPTIONAL     :: KEEP_IMAGES

! . Local arrays
REAL ( KIND = DP ), DIMENSION(:), ALLOCATABLE :: X, LENGTH, POS, COORD

! . Local scalars
INTEGER            :: NUM, I, J, II, IND, INI, FIN
LOGICAL            :: KEEP
REAL ( KIND = DP ) :: TOTAL_LENGTH, S, S_LEN

! . Local types
TYPE ( DCD_TYPE ) :: DCD_IN, DCD_OUT

! . Initialization and checks
KEEP = .FALSE.
IF ( PRESENT(KEEP_IMAGES) ) KEEP = KEEP_IMAGES

CALL DCD_INITIALIZE ( DCD_IN )
CALL DCD_ACTIVATE_READ ( FILE_IN, DCD_IN )

NUM = DCD_IN%NFRAMES
N3 = 3*NFREE

INI = 1
FIN = NUM
IF ( PRESENT(IM_INI) ) INI = IM_INI
IF ( PRESENT(IM_FIN) ) FIN = IM_FIN

IF ( INI > NUM .OR. FIN > NUM ) THEN
  CALL PRINT_ERROR ( "CHAIN_INTERPOLATE", "Wrong initial or final image.")
END IF

NUM = FIN - INI + 1
NCHAIN = NUM -2

IF ( KEEP .AND. ( NIMG < NUM  .OR. MOD(NIMG-NUM,NUM-1) /= 0 ) ) THEN
  CALL PRINT_ERROR ( "CHAIN_INTERPOLATE", "Wrong number of images to add if keeping the original ones.")
END IF

ALLOCATE ( X(NUM*N3), SPLINE(4,NUM,N3), LENGTH(NUM-1), POS(NIMG), COORD(N3) )

! . Read the coordinates from the input file
DO I = 1,INI-1
  CALL DCD_READ ( DCD_IN, ATMCRD )
END DO
DO I = 1,NUM
  II = ( I - 1 ) * N3
  CALL DCD_READ ( DCD_IN, ATMCRD )
  CALL CRD_TO_X ( ATMCRD, X(II+1:II+N3) )
END DO
CALL DCD_DEACTIVATE ( DCD_IN  )

! . Generate a cubic spline through all the points
CALL CUBIC_SPLINE( X )

! . Calcutate the length
DO I = 1, NUM-1
  CALL INTEGRATE_SPLINE( I, 0.0_DP, 1.0_DP, LENGTH(I) )
END DO
TOTAL_LENGTH = SUM( LENGTH )

! . Activate the output trajectory
CALL DCD_INITIALIZE ( DCD_OUT )
CALL DCD_ACTIVATE_WRITE ( FILE_OUT, DCD_OUT, "CORD", NATOMS, NFIXED, NIMG, QFIX = ATMFIX )

IF ( KEEP ) THEN

  ! . Calculate the number of interpolated images between every two original ones
  IND = (NIMG-NUM) / (NUM-1)

  DO I = 1, NUM-1
    II = ( I - 1 ) * N3

    ! . Write the first image in the segment
    CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
    CALL DCD_WRITE ( DCD_OUT, ATMCRD )

    ! . Write the interpolated images
    S_LEN = LENGTH(I) / REAL( IND+1, DP )
    DO J = 1, IND
      CALL FIND_POSITION ( I, LENGTH(I), J*S_LEN, S )
      COORD = SPLINE(1,I,:)       + SPLINE(2,I,:) * S + &
              SPLINE(3,I,:) * S*S + SPLINE(4,I,:) * S*S*S
      CALL X_TO_CRD ( COORD, ATMCRD )
      CALL DCD_WRITE ( DCD_OUT, ATMCRD )
    END DO

  END DO

  ! . Write the last image
  II = ( NUM - 1 ) * N3
  CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
  CALL DCD_WRITE ( DCD_OUT, ATMCRD )

ELSE

  ! . Calculate the position of the desired images
  DO I = 1, NIMG
    POS(I) = (I-1) * TOTAL_LENGTH / REAL( NIMG-1, DP )
  END DO

  ! . Write the first image
  CALL X_TO_CRD ( X(1:N3), ATMCRD )
  CALL DCD_WRITE ( DCD_OUT, ATMCRD )

  ! . Write the interpolated images
  TOTAL_LENGTH = 0.0_DP
  IND = 1
  DO I = 2, NIMG-1
    DO WHILE ( TOTAL_LENGTH + LENGTH(IND) <= POS(I) )
      TOTAL_LENGTH = TOTAL_LENGTH + LENGTH(IND)
      IND = IND + 1
    END DO
    S_LEN = POS(I) - TOTAL_LENGTH
    CALL FIND_POSITION ( IND, LENGTH(IND), S_LEN, S )
    COORD = SPLINE(1,IND,:)       + SPLINE(2,IND,:) * S + &
            SPLINE(3,IND,:) * S*S + SPLINE(4,IND,:) * S*S*S
    CALL X_TO_CRD ( COORD, ATMCRD )
    CALL DCD_WRITE ( DCD_OUT, ATMCRD )
  END DO

  ! . Write the last image
  II = ( NUM - 1 ) * N3
  CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
  CALL DCD_WRITE ( DCD_OUT, ATMCRD )

END IF

CALL DCD_DEACTIVATE ( DCD_OUT )
DEALLOCATE ( X, SPLINE, LENGTH, POS, COORD )

END SUBROUTINE CHAIN_INTERPOLATE

END MODULE NEB_SPLINE
