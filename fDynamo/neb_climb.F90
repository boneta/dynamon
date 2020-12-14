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
!                             The NEB Module
!===============================================================================
!
! . Subroutines:
!
!   NEB_OPTIMIZE_CLIMB     Optimize a chain trajectory.
!   VERLET_MINIMIZE_MAX    Minimize the forces using a damped verlet algorithm 
!                          moving only one image at a time. 
!   LBFGS_MINIMIZE         Minimize the forces using the L-BFGS algorithm
!   GET_LINEAR_CONSTAINTS  Impose Ekart conditions for translations and rotations.
!   PROJECT                Project the velocity onto the force vector.
!   TANGENT                Calculates the tangent at a chain point.
!   PARGRAD                Calculates the gradient parallel to the tangent.
!   CURVA                  Calculates the curvature of the final chain.
!   INTERPOLATE            Does a polynomial interpolation of the images.
!   FORTOT                 Calculates the forces acting on the NEB.
!                          The forces are a perpendicular compontent (FORPER) and
!                          a parallel component (FORPAR).
!   DIST_ARRAY             Create an array of all the distances between images.
!   DISTANCE               Calculate the distance between a pair of images.
!   CRD_TO_X               Copy ATMCRD to the variable array X.
!   X_TO_CRD               Copy X to ATMCRD.
!===============================================================================
MODULE NEB_CLIMB

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_LINE, PRINT_LINEBREAK, PRINT_PARAGRAPH, PRINT_PARAGRAPH_START, &
                        PRINT_PARAGRAPH_STOP, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS,               &
                        PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT,                     &
                        PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP, PRINT_TEXT
USE CONSTANTS

USE ATOMS,              ONLY : ATMCRD, ATMFIX, NATOMS, NATOMSQM, NFIXED, NFREE
USE DCD_IO
USE POTENTIAL_ENERGY,   ONLY : ATMDER, ENERGY, GRADIENT, ETOTAL, EQM
USE STRING,             ONLY : ENCODE_INTEGER
USE MOPAC_DENSITY,      ONLY : DENSITY_READ, DENSITY_WRITE, DENSITY_GUESS
USE LINEAR_ALGEBRA,     ONLY : NORM
USE LBFGS

IMPLICIT NONE
PRIVATE
PUBLIC :: NEB_OPTIMIZE_CLIMB
#ifndef PGPC
SAVE
#endif

! . Module scalars.
INTEGER            :: NCHAIN, NCONSTRAINTS, N3

! . Module arrays.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: D1, ECHAIN, ECHAINQM,  XBEG, XEND, &
                                                   GBEG, GEND, W
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: G_CONSTRAINT

!==============================================================================
CONTAINS
!==============================================================================


!-------------------------------------------------------------------------------------------------
SUBROUTINE NEB_OPTIMIZE_CLIMB (FILE_IN, FILE_OUT,  PRINT_FREQUENCY, STEP_NUMBER,   &
                  GRADIENT_TOLERANCE, STEP_SIZE,KFOR, NBINS,   SAVE_DENSITIES, &
                  READ_DENSITIES, TEMP, METHOD, KMIN, THRESHOLD, OPTIM_METHOD )
!-------------------------------------------------------------------------------------------------
  
  ! . Scalar arguments.
  CHARACTER ( LEN = * ), INTENT(IN) :: FILE_IN, FILE_OUT
  ! . Optional scalar arguments for the CHAIN energy.
  REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: KFOR, KMIN, THRESHOLD
  LOGICAL, INTENT(IN), OPTIONAL :: SAVE_DENSITIES, READ_DENSITIES
     
  ! . Optional scalar arguments for the minimization routine.
  CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: METHOD, OPTIM_METHOD
  INTEGER,               INTENT(IN), OPTIONAL :: PRINT_FREQUENCY, STEP_NUMBER
  REAL ( KIND = DP ),    INTENT(IN), OPTIONAL :: GRADIENT_TOLERANCE, STEP_SIZE, TEMP
  
  ! . Optional scalar argument for the final interpolation.
  INTEGER,            INTENT(IN), OPTIONAL :: NBINS
  
  ! . Local scalars.
  INTEGER            :: I, II, NPOINTS, METH, OPMETH
  REAL ( KIND = DP ) :: KFAC, DELK, EREF, ETHR, T, BETA
  LOGICAL :: SAVEDE, READDE
      
  ! . Local arrays.
  REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: G, X, PFLUX
  REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: TMPCRD
  CHARACTER (LEN = 3 ) :: ISTRING

  ! . Local types.
  TYPE(DCD_TYPE) :: DCD_IN, DCD_OUT

  ! . There are no free atoms.
  IF ( NFREE <= 0 ) RETURN   

  ! . Initialize the CHAIN energy parameters.
  IF ( PRESENT ( KFOR  ) ) THEN
   KFAC = KFOR
  ELSE
   KFAC = 1000.0_DP
  END IF
  IF ( PRESENT ( KMIN  ) ) THEN
   DELK = KFAC-KMIN
  ELSE
   DELK = 0.0_DP
  END IF
  IF ( PRESENT ( THRESHOLD ) ) THEN
    ETHR = THRESHOLD
  ELSE
    ETHR = -1.0_DP
  END IF
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
        CALL PRINT_ERROR ( "NEB_OPTIMIZE_CLIMB", "METHOD can only be INT, INT2 or DIFF." )
      END SELECT
    END IF
  END IF

  OPMETH = 1
  IF ( PRESENT(OPTIM_METHOD) ) THEN
    SELECT CASE ( OPTIM_METHOD )
    CASE ( "VERLET" )
      OPMETH = 1
    CASE ( "LBFGS" )
      OPMETH = 2
    CASE DEFAULT
      CALL PRINT_ERROR ( "NEB_OPTIMIZE_CLIMB", "Unrecognized optimization method." )
    END SELECT
  END IF

  IF ((SAVEDE .OR. READDE) .AND. NATOMSQM==0) CALL PRINT_ERROR("NEB_OPTIMIZE_CLIMB", &
   "If there are no QM atoms, the density cannot be saved or read.")

  ! . Activate the input trajectory.
  CALL DCD_INITIALIZE ( DCD_IN )
  CALL DCD_ACTIVATE_READ ( FILE_IN, DCD_IN )
  
  ! . Get the number of points in the chain.
  NPOINTS = DCD_IN%NFRAMES
  
  ! . Print out the options.
  CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#8888FF", VARIABLEWIDTH = 16 )
  CALL PRINT_SUMMARY_START ( "Nudged-Elastic-Band Path Calculation" )
  WRITE ( PRINT_LINE, "(I16)"    ) NPOINTS   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Structures" )
  WRITE ( PRINT_LINE, "(G16.10)" ) KFAC      ; CALL PRINT_SUMMARY_ELEMENT ( "Force Constant"       )
  WRITE ( PRINT_LINE, "(G16.10)" ) KFAC-DELK ; CALL PRINT_SUMMARY_ELEMENT ( "Min. Force Const."    )
  WRITE ( PRINT_LINE, "(G16.10)" ) ETHR      ; CALL PRINT_SUMMARY_ELEMENT ( "Energy Threshold"     )
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
      XBEG(1:N3), XEND(1:N3), GBEG(1:N3), GEND(1:N3), W(1:N3), PFLUX(1:NCHAIN+1) )

  ! . Define the weights
  W = 1.0_DP

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
     WRITE ( PRINT_LINE, "(I10)"   ) I           ; CALL PRINT_TABLE_ELEMENT
     WRITE ( PRINT_LINE, "(F15.5)" ) ECHAIN(I)   ; CALL PRINT_TABLE_ELEMENT
     WRITE ( PRINT_LINE, "(F15.5)" ) ECHAINQM(I) ; CALL PRINT_TABLE_ELEMENT
     WRITE ( PRINT_LINE, "(F15.5)" ) PARGRAD(I-1); CALL PRINT_TABLE_ELEMENT
     IF ( I == ( NCHAIN + 2 ) .OR. I==1 ) THEN
       WRITE ( PRINT_LINE, "(F10.3)" ) !0.0 
     ELSE      
       WRITE ( PRINT_LINE, "(F10.5)" ) CURVA(I-1)
     END IF
     CALL PRINT_TABLE_ELEMENT
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
    
  DEALLOCATE ( D1, ECHAIN, ECHAINQM, G, G_CONSTRAINT, TMPCRD, X, XBEG, XEND, GBEG, GEND, W, PFLUX )

  
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
   REAL ( KIND = DP ) :: GTOLERANCE, STPSIZ, TOLGRD, MODF2, MEAN_DIST

 
   ! . Local arrays
   REAL (KIND = DP), DIMENSION (1:SIZE(X)) :: VEL, FORCENEW, FORCE, FORPER, FORPAR
   REAL (KIND = DP), DIMENSION (1:NCHAIN) :: FORMOD
   INTEGER, DIMENSION (1) :: IMAX
   
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
   TOLGRD = 1.0E-3_DP

   ! . Assign the input parameters.
   IF ( PRESENT ( FORTOL             ) ) TOLGRD = FORTOL
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
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#55FF55", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Verlet Minimization Calculation" )
      WRITE ( PRINT_LINE, "(I16)"    ) NVAR   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Variables"  )
      WRITE ( PRINT_LINE, "(I16)"    ) NSTEP  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations" )
      WRITE ( PRINT_LINE, "(I16)"    ) NPRINT ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency"      )
      WRITE ( PRINT_LINE, "(G16.10)" ) TOLGRD ; CALL PRINT_SUMMARY_ELEMENT ( "Force Tolerance"      )
      WRITE ( PRINT_LINE, "(G16.10)" ) T      ; CALL PRINT_SUMMARY_ELEMENT ( "Temperature"          )
      WRITE ( PRINT_LINE, "(G16.10)" ) STPSIZ ; CALL PRINT_SUMMARY_ELEMENT ( "Step Size"            )
      CALL PRINT_SUMMARY_STOP

      ! . Print out the header for the iterations.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#55FF55", PAGEWIDTH = 80, &
            & VARIABLEWIDTHS = (/ 6, 15, 15, 15, 11 /) )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",         HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS Force",    HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS F perp.",    HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Mean Dist.",    HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Dist. std. dev.",    HEADER = .TRUE. )

   END IF
   
   ! .  Calculate the initial energies.
   DO I = 1,NCHAIN
     II = ( I - 1 ) * N3
     CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
     IF (READDE) THEN
       CALL ENCODE_INTEGER(I,ISTRING)
       CALL DENSITY_READ ("dens"//FILE_IN//TRIM(ISTRING), PRINT=.FALSE.)
     ELSE
       CALL DENSITY_GUESS
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

   VEL = 0.0_DP
   DO I = 1, NCHAIN
     CALL FORTOT(X, KFAC, NVAR, FORPER, FORPAR , I)
     II = (I - 1) * N3
     FORCENEW(II+1:II+N3) = FORPER(II+1:II+N3) + FORPAR(II+1:II+N3)
     FORMOD(I) = SUM((FORCENEW(II+1:II+N3))**2)
   END DO
   IMAX = MAXLOC(FORMOD)
   ! . Calculate the square of the force.
   MODF2= AMUA2PS2_TO_KJMOL**2 * DOT_PRODUCT ( FORCENEW, FORCENEW )
   D1 = DIST_ARRAY()/ SQRT( REAL (N3) )
   MEAN_DIST = SUM ( D1(1:NCHAIN+1) ) / REAL ( (NCHAIN + 1 ), DP )
   ! . Do some printing.
   IF ( NPRINT > 0 ) THEN
      WRITE ( PRINT_LINE, "(I6)"   ) 1
      CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F15.8)" ) SQRT ( MODF2 / REAL ( NVAR, DP ) )
      CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F15.8)" ) AMUA2PS2_TO_KJMOL*SQRT(SUM (FORPER**2)/ REAL ( NVAR, DP ) )
      CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F11.8)" ) MEAN_DIST
      CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F11.8)" ) SQRT(SUM ( (D1(1:NCHAIN+1)-MEAN_DIST)**2 )&
                                     / REAL (  NCHAIN + 1 , DP ))
      CALL PRINT_TABLE_ELEMENT
   END IF
   ! . Check for convergence.
   IF ( MODF2 <= GTOLERANCE ) THEN
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
     II = (IMAX(1) - 1) *N3
     X(II+1:II+N3) = X(II+1:II+N3) + STPSIZ*VEL(II+1:II+N3) + STPSIZ**2 * & 
                    FORCE(II+1:II+N3) / 2
     CALL FORTOT(X, KFAC, NVAR, FORPER, FORPAR , IMAX(1))
     FORCENEW = FORPAR + FORPER
     VEL(II+1:II+N3) = VEL(II+1:II+N3) + STPSIZ / 2 * (FORCE(II+1:II+N3) +&
                       FORCENEW(II+1:II+N3))
     CALL PROJECT (VEL, FORCENEW)
     DO J = 1, NCHAIN
       II = (J - 1) * N3
       FORMOD(J) = SUM((FORCENEW(II+1:II+N3))**2)
     END DO
     IMAX = MAXLOC(FORMOD)
     MODF2 = AMUA2PS2_TO_KJMOL**2 * DOT_PRODUCT (FORCENEW, FORCENEW)
     ! . Do some printing.
     IF ( NPRINT > 0 ) THEN
       IF ( MOD ( I, NPRINT ) == 0 ) THEN
         D1 = DIST_ARRAY()/ SQRT( REAL (N3) )
         MEAN_DIST = SUM ( D1(1:NCHAIN+1) ) / REAL ( (NCHAIN + 1 ), DP )  
         WRITE ( PRINT_LINE, "(I6)"   ) I
         CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F15.8)" ) SQRT ( MODF2 / REAL ( NVAR, DP ) )
         CALL PRINT_TABLE_ELEMENT
         WRITE ( print_line, "(F15.8)" ) AMUA2PS2_TO_KJMOL* &
              SQRT(dot_product (forper, forper)/ REAL ( NVAR, DP ) )
         CALL PRINT_TABLE_ELEMENT
         !write ( print_line, "(F15.8)" ) SQRT(ABS(dot_product (forpar, forper))/ REAL ( NVAR, DP ) )
         WRITE ( print_line, "(F11.8)" ) MEAN_DIST
         CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(F11.8)" ) SQRT(SUM ( (D1(1:NCHAIN+1)-MEAN_DIST)**2 &
                                         ) / REAL (  NCHAIN + 1 , DP ))
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
     IF ( MODF2 <= GTOLERANCE ) EXIT
   END DO
   IF (I>NSTEP) THEN
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
  D1 = DIST_ARRAY () / SQRT ( REAL ( N3, DP ) )
  
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
  INTEGER            :: IFAIL, NVAR, NSTEP, NPRINT, I, II, J, IMOVE
  REAL ( KIND = DP ) :: GTOLERANCE, STPSIZ, TOLGRD, MODF2, MEAN_DIST, STPSCL, MAXSTEP

  ! . Local arrays
  REAL (KIND = DP), DIMENSION (1:SIZE(X)) :: FORCE, FORPER, FORPAR, STEP

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
  STPSIZ = 1.0_DP
  TOLGRD = 1.0E-3_DP

  ! . Assign the input parameters.
  IF ( PRESENT ( FORTOL             ) ) TOLGRD = FORTOL
  IF ( PRESENT ( PRINT_FREQUENCY    ) ) NPRINT = PRINT_FREQUENCY
  IF ( PRESENT ( STEP_NUMBER        ) ) NSTEP  = STEP_NUMBER
  IF ( PRESENT ( STEP_SIZE          ) ) STPSIZ = STEP_SIZE

  ! . Check the input parameters.
  IF ( ( NPRINT < 0 ) .OR. ( NPRINT > NSTEP ) ) NPRINT = 0

  ! . Calculate the gradient tolerance.
  GTOLERANCE = REAL( NVAR, DP ) * TOLGRD * TOLGRD

  ! . Print out the control parameters.
  IF ( NPRINT > 0 ) THEN

    ! . Print out the options.
    CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#55FF55", VARIABLEWIDTH = 16 )
    CALL PRINT_SUMMARY_START ( "L-BFGS Minimization Calculation" )
    WRITE ( PRINT_LINE, "(I16)"    ) NVAR   ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Variables"  )
    WRITE ( PRINT_LINE, "(I16)"    ) NSTEP  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations" )
    WRITE ( PRINT_LINE, "(I16)"    ) NPRINT ; CALL PRINT_SUMMARY_ELEMENT ( "Print Frequency"      )
    WRITE ( PRINT_LINE, "(G16.10)" ) TOLGRD ; CALL PRINT_SUMMARY_ELEMENT ( "Force Tolerance"      )
    WRITE ( PRINT_LINE, "(G16.10)" ) T      ; CALL PRINT_SUMMARY_ELEMENT ( "Temperature"          )
    WRITE ( PRINT_LINE, "(G16.10)" ) STPSIZ ; CALL PRINT_SUMMARY_ELEMENT ( "Max. Step Size"       )
    CALL PRINT_SUMMARY_STOP

    ! . Print out the header for the iterations.
    CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#55FF55", PAGEWIDTH = 80, &
          & VARIABLEWIDTHS = (/ 6, 15, 15, 15, 11 /) )
    CALL PRINT_TABLE_START
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",         HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS Force",    HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS F perp.",    HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Mean Dist.",    HEADER = .TRUE. )
    CALL PRINT_TABLE_ELEMENT ( TEXT = "Dist. std. dev.",    HEADER = .TRUE. )

  END IF

  ! .  Calculate the initial energies.
  DO I = 1,NCHAIN
    II = ( I - 1 ) * N3
    CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )
    IF (READDE) THEN
      CALL ENCODE_INTEGER(I,ISTRING)
      CALL DENSITY_READ ("dens"//FILE_IN//TRIM(ISTRING), PRINT=.FALSE.)
    ELSE
      CALL DENSITY_GUESS
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

  DO J = 1, NCHAIN
    CALL FORTOT(X, KFAC, NVAR, FORPER, FORPAR, J)
    II = (J - 1) * N3
    FORCE(II+1:II+N3) = FORPER(II+1:II+N3) + FORPAR(II+1:II+N3)
  END DO
  MODF2 = AMUA2PS2_TO_KJMOL**2 * DOT_PRODUCT (FORCE, FORCE)
  D1 = DIST_ARRAY()/ SQRT( REAL (N3) )
  MEAN_DIST = SUM ( D1(1:NCHAIN+1) ) / REAL ( (NCHAIN + 1 ), DP )
  ! . Do some printing.
  IF ( NPRINT > 0 ) THEN
    WRITE ( PRINT_LINE, "(I6)"   ) 1
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.8)" ) SQRT( MODF2 / REAL(NVAR,DP) )
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F15.8)" ) AMUA2PS2_TO_KJMOL*SQRT(SUM (FORPER**2)/ REAL ( NVAR, DP ) )
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F11.8)" ) MEAN_DIST
    CALL PRINT_TABLE_ELEMENT
    WRITE ( PRINT_LINE, "(F11.8)" ) SQRT(SUM ( (D1(1:NCHAIN+1)-MEAN_DIST)**2 )&
                                    / REAL (  NCHAIN + 1 , DP ))
    CALL PRINT_TABLE_ELEMENT
  END IF

  CALL LBFGS_INITIALIZE( NVAR, 4 )

  DO I = 2,NSTEP
    IF (MODF2 <= GTOLERANCE) EXIT

    ! . This conversion seems to do more damage than good
    !CALL LBFGS_DATA( X, -FORCE/KJMOL_TO_AMUA2PS2 )
    CALL LBFGS_DATA( X, -FORCE )
    STEP = LBFGS_STEP( )

    MAXSTEP = 0.0_DP
    DO J = 1, NCHAIN
      II = ( J - 1 ) * N3
      STPSCL = SQRT( DOT_PRODUCT( STEP(II+1:II+N3), STEP(II+1:II+N3) ) / REAL( N3, DP ) )
      IF ( STPSCL > MAXSTEP ) THEN
        MAXSTEP = STPSCL
        IMOVE = J
      END IF
      STPSCL = MIN( STPSIZ/STPSCL, 1.0_DP )
      STEP(II+1:II+N3) = STPSCL * STEP(II+1:II+N3)
    END DO
    !MAXSTEP = SQRT( DOT_PRODUCT( STEP, STEP ) / REAL( N3, DP ) )
    !STPSCL = MIN( STPSIZ/MAXSTEP, 1.0_DP )
    !STEP = STPSCL * STEP

    ! . Move all images at once
    X = X + STEP
    DO J = 1, NCHAIN
      CALL FORTOT(X, KFAC, NVAR, FORPER, FORPAR, J)
      II = (J - 1) * N3
      FORCE(II+1:II+N3) = FORPER(II+1:II+N3) + FORPAR(II+1:II+N3)
    END DO

    ! . Move just one image
    !II = ( IMOVE - 1 ) * N3
    !X(II+1:II+N3) = X(II+1:II+N3) + STEP(II+1:II+N3)
    !CALL FORTOT(X, KFAC, NVAR, FORPER, FORPAR, IMOVE)
    !FORCE(II+1:II+N3) = FORPER(II+1:II+N3) + FORPAR(II+1:II+N3)

    MODF2 = AMUA2PS2_TO_KJMOL**2 * DOT_PRODUCT (FORCE, FORCE)

    IF ( NPRINT > 0 ) THEN
      IF ( MOD ( I, NPRINT ) == 0 ) THEN
        D1 = DIST_ARRAY()/ SQRT( REAL (N3) )
        MEAN_DIST = SUM ( D1(1:NCHAIN+1) ) / REAL ( (NCHAIN + 1 ), DP )  
        WRITE ( PRINT_LINE, "(I6)"   ) I
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F15.8)" ) SQRT( MODF2 / REAL(NVAR,DP) )
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F15.8)" ) AMUA2PS2_TO_KJMOL* &
             SQRT(DOT_PRODUCT (FORPER, FORPER)/ REAL ( NVAR, DP ) )
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F11.8)" ) MEAN_DIST
        CALL PRINT_TABLE_ELEMENT
        WRITE ( PRINT_LINE, "(F11.8)" ) SQRT(SUM ( (D1(1:NCHAIN+1)-MEAN_DIST)**2 &
                                        ) / REAL (  NCHAIN + 1 , DP ))
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

    IF (MODF2 <= GTOLERANCE) EXIT
  END DO

  CALL LBFGS_INITIALIZE(0,0,DELETE=.TRUE.)

  IF (I>NSTEP) THEN
    IFAIL=131
  ELSE
    IFAIL=0
  END IF

  !---------------------------------------------------------------------------
  ! . Finish up.
  !---------------------------------------------------------------------------

  ! . Calculate the distances between neighbouring structures.
  D1 = DIST_ARRAY () / SQRT ( REAL ( N3, DP ) )

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
!    VELFOR = DOT_PRODUCT (VEL, FOR)
!    IF (VELFOR > 0.0_DP) THEN
!      VEL = VELFOR * FOR / DOT_PRODUCT (FOR, FOR)
!    ELSE
!      VEL = 0.0_DP
!    END IF

    ! . Project each structure velocity on each structure force.
!    DO I = 1, SIZE(VEL), N3
!     VELFOR = DOT_PRODUCT (VEL(I:I+N3-1), FOR(I:I+N3-1))
!      IF (VELFOR > 0.0) THEN
!       VEL(I:I+N3-1) = VELFOR * FOR(I:I+N3-1) / DOT_PRODUCT (FOR(I:I+N3-1), FOR(I:I+N3-1))
!      ELSE
!        VEL(I:I+N3-1) = 0.0_DP
!      END IF
!    END DO

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
!    WHERE (VEL*FOR < 0.0_DP) VEL=0.0_DP

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
  REAL (KIND = DP), DIMENSION (N3) :: TANG, DELTARA, DELTARB
  

  ! . Calculate the tangent to each image.
  CALL TANGENT ( TANG, DELTARA, DELTARB, I+1 )
  ! . Calculate the index into the coordinate array.
  II = ( I - 1 ) * N3
  ! . Calculate parallel force
  IF (I==0) THEN
    PARGRAD = DOT_PRODUCT(GBEG, TANG)
  ELSE IF (I==NCHAIN+1) THEN
    PARGRAD = DOT_PRODUCT(GEND, TANG)
  ELSE 
    PARGRAD = DOT_PRODUCT(G(II+1:II+N3), TANG)
  END IF
  PARGRAD = PARGRAD * SQRT (REAL (N3, DP))

  END FUNCTION PARGRAD

  !---------------------------------------------------------------------      
  FUNCTION CURVA(ISTRUCT)
  !---------------------------------------------------------------------

  ! . Scalar arguments
  INTEGER, INTENT(IN):: ISTRUCT
  
  ! . Function value.
  REAL (KIND = DP) :: CURVA
  
  ! . Local scalars.
  INTEGER :: I
  
  ! . Local arrays.
  REAL (KIND = DP), DIMENSION (-1:1,N3) :: TANGI
  REAL (KIND = DP), DIMENSION (N3) :: DELTARA, DELTARB
  REAL (KIND = DP) :: DB, DA
  
  IF (ISTRUCT<1 .OR. ISTRUCT>NCHAIN) THEN
    CURVA = 0.0_DP
    RETURN
  END IF

  DO I=-1,1
    CALL TANGENT ( TANGI(I,:), DELTARA, DELTARB, ISTRUCT+I+1 )
    IF (I == 0) THEN
      DA = NORM(DELTARA)
      DB = NORM(DELTARB)
    END IF
  END DO

  CURVA = ACOS(DOT_PRODUCT(TANGI(-1,:), TANGI(1,:)))/ (DB+DA)
 
  END FUNCTION CURVA  

  !---------------------------------------------------------------------      
  SUBROUTINE FLUX(PARTIAL)
  !---------------------------------------------------------------------

  ! . Array arguments
  REAL ( KIND = DP ), DIMENSION(NCHAIN+1), INTENT(OUT):: PARTIAL

  ! . Local scalars
  REAL ( KIND = DP ) :: EBASE
  INTEGER :: I

  EBASE = MIN( ECHAIN(1), ECHAIN(NCHAIN+2) )

  DO I = 1, NCHAIN+1
    PARTIAL(I) = 0.5_DP * D1(I) * ( EXP( BETA * (ECHAIN(I  )-EBASE) ) + &
                                    EXP( BETA * (ECHAIN(I+1)-EBASE) ) )
  END DO

  END SUBROUTINE FLUX

  !---------------------------------------------------------------------      
  SUBROUTINE INTERPOLATE (NBININ)
  !---------------------------------------------------------------------
  
  ! . Input scalars.
  INTEGER,            INTENT(IN),  OPTIONAL :: NBININ
  
  ! . Local scalars.
  INTEGER :: NBINS, IREF
  REAL (KIND = DP) :: STEP, POS, A, B, C, D, ESTIME, SUMREF
  LOGICAL :: RECALC
  
  ! .  Check if there is an input and if its value is reasonable.
  NBINS = 0
  IF ( PRESENT ( NBININ  ) ) THEN
   NBINS = NBININ
  ELSE
   NBINS = 0
  END IF
  IF (NBINS==0) THEN
    CALL PRINT_PARAGRAPH (Text="No interpolation requested.")
    RETURN
  ELSE IF (NBINS<NPOINTS) THEN
    CALL PRINT_ERROR ("INTERPOLATE", " NBINS has to be greater than the number of chain structures.")
  END IF

  ! . Calcualte the distance step.
  STEP = SUM(D1)/(NBINS-1)

  ! . Write out the header.     
  CALL PRINT_TABLE_OPTIONS ( COLUMNS = 3, HEADER_COLOR = "#8888FF", VARIABLEWIDTHS = (/ 15, 15, 15 /) )
  CALL PRINT_TABLE_START
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Polinomial Interpolation", COLSPAN = 3, HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Distance",   HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Interpolated points",      HEADER = .TRUE. )
  CALL PRINT_TABLE_ELEMENT ( TEXT = "Original points",      HEADER = .TRUE. )
  
  ! . Loop over the structures.
  IREF = 1
  RECALC = .TRUE.
  SUMREF=0.0_DP
  DO I = 2,NBINS-1
     ! . Calculate the coordinate and the pair of images it is limited by.
     POS = (I-1)*STEP
     DO WHILE (POS > SUM(D1(1:IREF)))
       IREF = IREF + 1
       RECALC = .TRUE.
       SUMREF = SUM(D1(1:IREF)) - D1(IREF)
     END DO
     
     IF (RECALC) THEN
       A = 2*(ECHAIN(IREF)-ECHAIN(IREF+1))/D1(IREF)**3+ (PARGRAD(IREF-1)+ &
           PARGRAD(IREF))/D1(IREF)**2 
       B = 3*(ECHAIN(IREF+1)-ECHAIN(IREF))/D1(IREF)**2 -(2*PARGRAD(IREF-1)+ &
           PARGRAD(IREF))/D1(IREF)
       C = PARGRAD(IREF-1)
       D = ECHAIN(IREF)
       RECALC = .FALSE.
       ESTIME = A*(SUM(D1(1:IREF-1))-SUMREF)**3+B*(SUM(D1(1:IREF-1))-SUMREF)**2&
                +C*(SUM(D1(1:IREF-1))-SUMREF)+D
       WRITE ( PRINT_LINE, "(F15.5)"   ) SUM(D1(1:IREF-1))
                                                     CALL PRINT_TABLE_ELEMENT
       WRITE ( PRINT_LINE, "(F15.5)" ) ESTIME      ; CALL PRINT_TABLE_ELEMENT
       WRITE ( PRINT_LINE, "(F15.5)" ) ECHAIN(IREF); CALL PRINT_TABLE_ELEMENT
   END IF
     ESTIME = A*(POS-SUMREF)**3+B*(POS-SUMREF)**2+C*(POS-SUMREF)+D
     WRITE ( PRINT_LINE, "(F15.5)" ) POS           ; CALL PRINT_TABLE_ELEMENT
     WRITE ( PRINT_LINE, "(F15.5)" ) ESTIME        ; CALL PRINT_TABLE_ELEMENT
     WRITE ( PRINT_LINE, "(A15)" )   ""            ; CALL PRINT_TABLE_ELEMENT
  
  END DO
  WRITE ( PRINT_LINE, "(F15.5)"   ) SUM(D1(1:IREF)) 
                                                  CALL PRINT_TABLE_ELEMENT
  WRITE ( PRINT_LINE, "(F15.5)" ) ECHAIN(IREF+1); CALL PRINT_TABLE_ELEMENT
  WRITE ( PRINT_LINE, "(F15.5)" ) ECHAIN(IREF+1); CALL PRINT_TABLE_ELEMENT
  
  CALL PRINT_TABLE_STOP
 
  END SUBROUTINE INTERPOLATE  
   
  !---------------------------------------------------------------------      
  SUBROUTINE FORTOT(X, K, N, FORPER, FORPAR, ISTRUCT)
  !---------------------------------------------------------------------

  ! . Scalar arguments
  REAL ( KIND = DP ), INTENT(IN):: K
  INTEGER, INTENT(IN):: N, ISTRUCT
  
  ! . Array arguments
  REAL ( KIND = DP ), DIMENSION(N), INTENT(INOUT) :: FORPER, FORPAR
  REAL ( KIND = DP ), DIMENSION(N), INTENT(IN) :: X

  ! . Local arrays
  REAL (KIND = DP), DIMENSION(N3) :: DELTARB, DELTARA, TANG, FORCE, NORMAL
  
  ! . Local scalars.
  REAL ( KIND = DP) :: VA, VX, VB, KB, KA, EMAX, BETAU
  INTEGER :: II, I

  ! . Calculate the index into the coordinate array.
  II = ( ISTRUCT - 1 ) * N3
  ! . Copy the variables to the coordinate array.
  CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )

  ! . Calculate the energy and derivatives.
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

  ! . Save the energy.
  ECHAIN(ISTRUCT+1) = ETOTAL
  ECHAINQM(ISTRUCT+1) = EQM

  ! . Copy over the derivatives.
  CALL CRD_TO_X ( ATMDER, G(II+1:II+N3) )

  DO I = ISTRUCT-1, ISTRUCT+1
    IF ((I<1) .OR. (I>NCHAIN)) CYCLE 

    ! . Calculate the tangent to each image.
    CALL TANGENT ( TANG, DELTARA, DELTARB, I+1 )

    VB = ECHAIN(I)
    VX = ECHAIN(I+1)
    VA = ECHAIN(I+2)

    ! . Calculate the index into the coordinate array.
    II = ( I - 1 ) * N3

    ! . Calculate the normal.
    NORMAL = DELTARB - DELTARA
    NORMAL = NORMAL - DOT_PRODUCT(NORMAL, TANG) * TANG
    CALL CONSTRAIN ( NORMAL )
    ! . Set the normal to zero in the linear case
    IF ( NORM(NORMAL) > 1.0E-7_DP ) THEN
      NORMAL = NORMAL/NORM(NORMAL)
    ELSE
      NORMAL = 0.0_DP
    END IF

    ! . Calculate the force.
    SELECT CASE (METH)
    CASE ( 0 )
      ! . Temperature = 0 K
      FORCE = G(II+1:II+N3)
    CASE ( 1 )
      ! . Differential method
      FORCE = G(II+1:II+N3) + CURVA(I)/BETA * NORMAL
    CASE ( 2 )
      ! . Integral method
      FORCE = ( (NORM(DELTARB)+NORM(DELTARA)) * BETA * G(II+1:II+N3) + &
                DELTARB/NORM(DELTARB) * ( 1.0_DP + EXP(BETA*(VB-VX)) ) - &
                DELTARA/NORM(DELTARA) * ( 1.0_DP + EXP(BETA*(VA-VX)) ) ) * &
              0.5_DP * EXP(BETA*VX)
    CASE ( 3 )
      ! . Integral method #2 (old)
      ! (where does this come from?)
      BETAU = BETA * ( VA - VB )
      IF ( ABS(BETAU) < 40.0_DP ) THEN
        FORCE = G(II+1:II+N3) + 2.0_DP/BETA * &
                ( DELTARB/NORM(DELTARB) - EXP(BETAU) * DELTARA/NORM(DELTARA) ) / &
                ( NORM(DELTARB) - EXP(BETAU) * NORM(DELTARA) )
      ELSE IF ( BETAU > 0.0_DP ) THEN
        FORCE = G(II+1:II+N3) - 2.0_DP/BETA * DELTARA/NORM(DELTARA)**2
      ELSE IF ( BETAU < 0.0_DP ) THEN
        FORCE = G(II+1:II+N3) + 2.0_DP/BETA * DELTARB/NORM(DELTARB)**2
      END IF
    END SELECT

    ! . Project out the rigid body motions from the force.
    CALL CONSTRAIN ( FORCE )

    ! . Calculate parallel force
    IF ((VA<VX) .AND. (VX>VB) .AND. (MIN(VX-VA,VX-VB) > ETHR) .AND. (ETHR >= 0.0_DP) ) THEN
      FORPAR (II+1:II+N3) = DOT_PRODUCT(FORCE,TANG)*TANG
    ELSE IF ((VA>VX) .AND. (VX<VB) .AND. (MIN(VA-VX,VB-VX) > ETHR) .AND. (ETHR >= 0.0_DP) ) THEN
      FORPAR (II+1:II+N3) = -DOT_PRODUCT(FORCE,TANG)*TANG
    ELSE
      IF ( DELK > 0.0_DP ) THEN
        EMAX = MAXVAL(ECHAIN)
        KB = K - DELK * MIN( (EMAX - MAX( ECHAIN(I  ), ECHAIN(I+1) ) )/(EMAX-EREF), 1.0_DP )
        KA = K - DELK * MIN( (EMAX - MAX( ECHAIN(I+1), ECHAIN(I+2) ) )/(EMAX-EREF), 1.0_DP )
      ELSE
        KA = K; KB = K
      END IF
      FORPAR (II+1:II+N3) = TANG * (KA*NORM(DELTARA) - KB*NORM(DELTARB))
    END IF

    ! . Calculate perpendicular force
    FORPER (II+1:II+N3) = -( FORCE - DOT_PRODUCT(FORCE, TANG) * TANG )

    !. Change to atomic units.
    FORPAR(II+1:II+N3) =  KJMOL_TO_AMUA2PS2 * FORPAR(II+1:II+N3)
    FORPER(II+1:II+N3) =  KJMOL_TO_AMUA2PS2 * FORPER(II+1:II+N3)
  END DO

  END SUBROUTINE FORTOT
    
  !---------------------------
  FUNCTION DIST_ARRAY ()
  !---------------------------
  
  ! . Array arguments
  REAL (KIND = DP), DIMENSION (NCHAIN+1) :: DIST_ARRAY
  
  ! . Local scalar
  INTEGER :: I, II, JJ
  
  ! . Calculate the distances between neighbouring structures.
  DIST_ARRAY(1) = NORM ( DISTANCE ( XBEG, X(1:N3) ) )
  DO I = 2,NCHAIN
     II = ( I - 2 ) * N3
     JJ = II + N3
     DIST_ARRAY(I) = NORM ( DISTANCE ( X(II+1:II+N3), X(JJ+1:JJ+N3) ) )
  END DO
  II = ( NCHAIN - 1 ) * N3
  DIST_ARRAY(NCHAIN+1) = NORM ( DISTANCE ( X(II+1:II+N3), XEND ) )
  
  END FUNCTION DIST_ARRAY    
  
END SUBROUTINE NEB_OPTIMIZE_CLIMB


!---------------------------
FUNCTION DISTANCE ( X1, X2 )
!---------------------------

! . Array arguments.
REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: X1, X2

! . Function declarations.
REAL ( KIND = DP ), DIMENSION (SIZE(X1)) :: DISTANCE

! . Calculate the distance.
DISTANCE = W * ( X2 - X1 )

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

END MODULE NEB_CLIMB
