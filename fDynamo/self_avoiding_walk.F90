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
!                         The Self-Avoiding Walk Module
!===============================================================================
!
! . Subroutines:
!
!   CHAIN_EXPAND                    Expand a chain trajectory.
!   CHAIN_INITIALIZE                Initialize a chain trajectory.
!   CHAIN_OPTIMIZE                  Optimize a chain trajectory.
!   CHAIN_TEST                      Test the CHAIN_OPTIMIZE function, EGCALC.
!   EGCALC                          Calculate the CHAIN_OPTIMIZE function.
!
!   CRD_TO_X                        Copy ATMCRD to the variable array X.
!   X_TO_CRD                        Copy X to ATMCRD.
!
! . Notes:
!
!   It may be necessary to add an extra command to save the density matrices
!   of each structure during the optimization.
!
!===============================================================================
MODULE SELF_AVOIDING_WALK

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_ERROR, PRINT_LINE, PRINT_LINEBREAK, PRINT_PARAGRAPH, PRINT_PARAGRAPH_START, &
                        PRINT_PARAGRAPH_STOP, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS,               &
			PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT,                     &
			PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP, PRINT_TEXT

USE CONJUGATE_GRADIENT, ONLY : CONJUGATE_GRADIENT_MINIMIZE

USE ATOMS,              ONLY : ATMCRD, ATMFIX, NATOMS, NFIXED, NFREE
USE DCD_IO
USE POTENTIAL_ENERGY,   ONLY : ATMDER, ENERGY, GRADIENT, ETOTAL
USE SUPERIMPOSE,        ONLY : SUPERIMPOSE_QUATERNION

IMPLICIT NONE
PRIVATE
PUBLIC :: CHAIN_EXPAND, CHAIN_INITIALIZE, CHAIN_OPTIMIZE
#ifndef PGPC
SAVE
#endif

! . Module scalars.
INTEGER            :: NCHAIN, NCONSTRAINTS, N3
REAL ( KIND = DP ) :: GAMFAC, LAMFAC, RHOFAC

! . Module arrays.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: D1, D2, ECHAIN, XBEG, XEND
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: G_CONSTRAINT

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------------------------------
   SUBROUTINE CHAIN_EXPAND ( FILE_IN, FILE_OUT, NINSERT )
   !-----------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE_IN, FILE_OUT
   INTEGER,               INTENT(IN) :: NINSERT

   ! . Local scalars.
   INTEGER :: IADD, IOLD, NPTNEW, NPTOLD

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: CRDBEG, CRDEND, DXSTEP

   ! . Local types.
   TYPE(DCD_TYPE) :: DCD_IN, DCD_OUT

   ! . There are no free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Check the number of points to insert between each point on the trajectory.
   IF ( NINSERT < 1 ) CALL PRINT_ERROR ( "CHAIN_EXPAND", "Invalid number of insertions." )

   ! . Activate the input trajectory.
   CALL DCD_INITIALIZE ( DCD_IN )
   CALL DCD_ACTIVATE_READ ( FILE_IN, DCD_IN )

   ! . Get the number of frames on the trajectory.
   NPTOLD = DCD_IN%NFRAMES

   ! . Calculate and check the number of new points on the trajectory.
   NPTNEW = NPTOLD + ( NPTOLD - 1 ) * NINSERT
   IF ( MOD ( NPTNEW, 2 ) == 0 ) CALL PRINT_ERROR ( "CHAIN_EXPAND", "The number of trajectory points must be odd." )
 
   ! . Activate the output trajectory.
   CALL DCD_INITIALIZE ( DCD_OUT )
   CALL DCD_ACTIVATE_WRITE ( FILE_OUT, DCD_OUT, "CORD", NATOMS, NFIXED, NPTNEW, QFIX = ATMFIX )

   ! . Read in and write out the first structure.
   CALL DCD_READ  ( DCD_IN,  CRDBEG )
   CALL DCD_WRITE ( DCD_OUT, CRDBEG )

   ! . Loop over the points in the old trajectory.
   DO IOLD = 1,(NPTOLD-1)

      ! . Read the structure.
      CALL DCD_READ ( DCD_IN, CRDEND )

      ! . Orient the two structures with respect to each other (if all atoms are free).
      IF ( NFIXED <= 0 ) CALL SUPERIMPOSE_QUATERNION ( CRDBEG, CRDEND, PRINT = .FALSE. )

      ! . Calculate the step between structures.
! . Error corrected by Ignacio.
      DXSTEP = ( CRDEND - CRDBEG ) / REAL ( ( NINSERT + 1 ), DP )

      ! . Form the intermediate structures.
      DO IADD = 1,NINSERT
         CRDBEG = CRDBEG + DXSTEP
         CALL DCD_WRITE ( DCD_OUT, CRDBEG )
      END DO

      ! . Write out the final structure.
      CALL DCD_WRITE ( DCD_OUT, CRDEND )

      ! . Reset CRDBEG.
      CRDBEG = CRDEND

   END DO

   ! . Deactivate the trajectories.
   CALL DCD_DEACTIVATE ( DCD_IN  )
   CALL DCD_DEACTIVATE ( DCD_OUT )

   ! . Write out some information.
   CALL PRINT_PARAGRAPH_START
   WRITE ( PRINT_LINE, "(A,I6,A)"   ) "Chain Expansion:", NPTOLD, " structures read from  "//TRIM ( FILE_IN  )
   CALL PRINT_TEXT
   CALL PRINT_LINEBREAK
   WRITE ( PRINT_LINE, "(17X,I6,A)" )                     NPTNEW, " structures written to "//TRIM ( FILE_OUT )
   CALL PRINT_TEXT
   CALL PRINT_PARAGRAPH_STOP

   END SUBROUTINE CHAIN_EXPAND

   !-----------------------------------------------------------------
   SUBROUTINE CHAIN_INITIALIZE ( FILE, NPOINTS, REACTANTS, PRODUCTS )
   !-----------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE
   INTEGER,               INTENT(IN) :: NPOINTS

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(IN) :: PRODUCTS, REACTANTS

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: SAME = 1.0E-10_DP

   ! . Local scalars.
   INTEGER :: I

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: CRDEND, CRDXYZ, DXSTEP

   ! . Local types.
   TYPE(DCD_TYPE) :: TRAJECTORY

   ! . There are no free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Check the number of points on the trajectory.
   IF ( NPOINTS < 3 ) THEN
      CALL PRINT_ERROR ( "CHAIN_INITIALIZE", "Too few points in trajectory." )
   ELSE IF ( MOD ( NPOINTS, 2 ) == 0 ) THEN
      CALL PRINT_ERROR ( "CHAIN_INITIALIZE", "The number of trajectory points must be odd." )
   END IF

   ! . Activate the trajectory.
   CALL DCD_INITIALIZE ( TRAJECTORY )
   CALL DCD_ACTIVATE_WRITE ( FILE, TRAJECTORY, "CORD", NATOMS, NFIXED, NPOINTS, QFIX = ATMFIX )

   ! . Copy the reactant and product coordinates.
   CRDXYZ = REACTANTS
   CRDEND = PRODUCTS

   ! . There are fixed atoms.
   IF ( NFIXED > 0 ) THEN

      ! . Check that the coordinates of the fixed atoms in both structures are the same.
      DO I = 1,NATOMS
         IF ( ATMFIX(I) ) THEN
            IF ( MAXVAL ( ABS ( CRDXYZ(1:3,I) - CRDEND(1:3,I) ) ) > SAME ) THEN
               CALL PRINT_ERROR ( "CHAIN_INITIALIZE", "Fixed atom coordinate mismatch in reactant and product structures.", I )
            END IF
         END IF
      END DO

   ! . There are no fixed atoms.
   ELSE

      ! . Orient the structure in PRODUCTS with respect to that in REACTANTS.
      CALL SUPERIMPOSE_QUATERNION ( CRDXYZ, CRDEND, PRINT = .FALSE. )

   END IF

   ! . Calculate the step size between coordinate sets.
   DXSTEP = ( CRDEND - CRDXYZ ) / REAL ( ( NPOINTS - 1 ), DP )

   ! . Write out the reactant coordinates to the trajectory.
   CALL DCD_WRITE ( TRAJECTORY, REACTANTS )

   ! . Calculate and write out the intermediate structures.
   DO I = 1,(NPOINTS-2)
      CRDXYZ = CRDXYZ + DXSTEP
      CALL DCD_WRITE ( TRAJECTORY, CRDXYZ )
   END DO

   ! . Write out the product coordinates.
   CALL DCD_WRITE ( TRAJECTORY, CRDEND )

   ! . Deactivate the trajectory.
   CALL DCD_DEACTIVATE ( TRAJECTORY )

   ! . Write out some information.
   WRITE ( PRINT_LINE, "(A,I6,A)" ) "Chain Initialization: ", NPOINTS, " structures written to "//TRIM ( FILE )
   CALL PRINT_PARAGRAPH

   END SUBROUTINE CHAIN_INITIALIZE

   !-------------------------------------------------------------------------------------------------
   SUBROUTINE CHAIN_OPTIMIZE ( FILE_IN, FILE_OUT, GAMMA, LAMBDA, RHO, PRINT_FREQUENCY, STEP_NUMBER, &
                                                                      GRADIENT_TOLERANCE, STEP_SIZE )
   !-------------------------------------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE_IN, FILE_OUT

   ! . Optional scalar arguments for the CHAIN energy.
   REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: GAMMA, LAMBDA, RHO

   ! . Optional scalar arguments for the minimization routine.
   INTEGER,            INTENT(IN),  OPTIONAL :: PRINT_FREQUENCY, STEP_NUMBER
   REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: GRADIENT_TOLERANCE, STEP_SIZE

   ! . Local scalars.
   INTEGER            :: I, IFLAG, II, NPOINTS
   REAL ( KIND = DP ) :: DIST1, DIST2, F

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: G, X
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: TMPCRD

   ! . Local types.
   TYPE(DCD_TYPE) :: DCD_IN, DCD_OUT

   ! . There are no free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Initialize the CHAIN energy parameters.
   GAMFAC =  100.0_DP
   LAMFAC =    2.0_DP
   RHOFAC = 5000.0_DP

   ! . Check the input options for the CHAIN energy.
   IF ( PRESENT ( GAMMA  ) ) GAMFAC = GAMMA
   IF ( PRESENT ( LAMBDA ) ) LAMFAC = LAMBDA
   IF ( PRESENT ( RHO    ) ) RHOFAC = RHO

   ! . Activate the input trajectory.
   CALL DCD_INITIALIZE ( DCD_IN )
   CALL DCD_ACTIVATE_READ ( FILE_IN, DCD_IN )

   ! . Get the number of points in the chain.
   NPOINTS = DCD_IN%NFRAMES

   ! . Check the data on the trajectory.
   IF ( MOD ( NPOINTS, 2 ) == 0 ) CALL PRINT_ERROR ( "CHAIN_OPTIMIZE", "The number of trajectory points must be odd." )

   ! . Activate the output trajectory.
   CALL DCD_INITIALIZE ( DCD_OUT )
   CALL DCD_ACTIVATE_WRITE ( FILE_OUT, DCD_OUT, "CORD", NATOMS, NFIXED, NPOINTS, QFIX = ATMFIX )

   ! . Print out the options.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#8888FF", VARIABLEWIDTH = 16 )
   CALL PRINT_SUMMARY_START ( "Self-Avoiding Walk Path Calculation" )
   WRITE ( PRINT_LINE, "(I16)"    ) NPOINTS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Structures" )
   WRITE ( PRINT_LINE, "(G16.10)" ) GAMFAC  ; CALL PRINT_SUMMARY_ELEMENT ( "Gamma"                )
   WRITE ( PRINT_LINE, "(G16.10)" ) LAMFAC  ; CALL PRINT_SUMMARY_ELEMENT ( "Lambda"               )
   WRITE ( PRINT_LINE, "(G16.10)" ) RHOFAC  ; CALL PRINT_SUMMARY_ELEMENT ( "Rho"                  )
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
   ALLOCATE ( D1(1:NCHAIN+1), D2(1:NCHAIN), ECHAIN(1:NCHAIN+2), G(1:NCHAIN*N3), &
              G_CONSTRAINT(1:N3,6), TMPCRD(1:3,1:NATOMS), X(1:NCHAIN*N3), XBEG(1:N3), XEND(1:N3) )

   ! . Save the existing coordinates.
   TMPCRD = ATMCRD

   ! . Read in the first structure from the trajectory, calculate its energy and save its coordinates.
   ! . Note that if there are free atoms their coordinates are ALWAYS present in ATMCRD.
   CALL DCD_READ ( DCD_IN, ATMCRD )
   CALL ENERGY ( PRINT = .FALSE. ) ; ECHAIN(1) = ETOTAL
   CALL CRD_TO_X ( ATMCRD, XBEG )

   ! . Read in the intermediate structures from the trajectory.
   DO I = 1,NCHAIN
      CALL DCD_READ ( DCD_IN, ATMCRD )
      II = ( I - 1 ) * N3
      CALL CRD_TO_X ( ATMCRD, X(II+1:II+N3) )
   END DO

   ! . Read in the last structure from the trajectory, calculate its energy and save its coordinates.
   CALL DCD_READ ( DCD_IN, ATMCRD )
   CALL ENERGY ( PRINT = .FALSE. ) ; ECHAIN(NCHAIN+2) = ETOTAL
   CALL CRD_TO_X ( ATMCRD, XEND )
   
   ! . Calculate the linear constraints using the central structure.
   II = ( ( NCHAIN + 1 ) / 2 - 1 ) * N3
   CALL GET_LINEAR_CONSTRAINTS ( X(II+1:II+N3) )

   ! . Perform the minimization.
   CALL CONJUGATE_GRADIENT_MINIMIZE ( EGCALC, X, IFLAG, PRINT_FREQUENCY, STEP_NUMBER, STEP_SIZE, GRADIENT_TOLERANCE )

   ! . Calculate the energies of the best structures.
   CALL EGCALC ( X, F, G )

   ! . Scale the distances.
   D1 = D1 / SQRT ( REAL ( N3, DP ) )
   D2 = SQRT ( D2 / REAL ( N3, DP ) )

   ! . Write out the header.
   CALL PRINT_TABLE_OPTIONS ( COLUMNS = 4, HEADER_COLOR = "#8888FF", VARIABLEWIDTHS = (/ 20, 20, 20, 20 /) )
   CALL PRINT_TABLE_START
   CALL PRINT_TABLE_ELEMENT ( TEXT = "Chain Structures", COLSPAN = 4, HEADER = .TRUE. )
   CALL PRINT_TABLE_ELEMENT ( TEXT = "Structure",   HEADER = .TRUE. )
   CALL PRINT_TABLE_ELEMENT ( TEXT = "Energy",      HEADER = .TRUE. )
   CALL PRINT_TABLE_ELEMENT ( TEXT = "Dist(i,i+1)", HEADER = .TRUE. )
   CALL PRINT_TABLE_ELEMENT ( TEXT = "Dist(i,i+2)", HEADER = .TRUE. )

   ! . Loop over the structures.
   DO I = 1,NCHAIN+2
      IF ( I == ( NCHAIN + 1 ) ) THEN
         DIST1 = D1(I)  ; DIST2 = 0.0_DP
      ELSE IF ( I == ( NCHAIN + 2 ) ) THEN
         DIST1 = 0.0_DP ; DIST2 = 0.0_DP
      ELSE
         DIST1 = D1(I)  ; DIST2 = D2(I)
      END IF
      WRITE ( PRINT_LINE, "(I20)"   ) I         ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F20.5)" ) ECHAIN(I) ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F20.5)" ) DIST1     ; CALL PRINT_TABLE_ELEMENT
      WRITE ( PRINT_LINE, "(F20.5)" ) DIST2     ; CALL PRINT_TABLE_ELEMENT
   END DO

   ! . Write out the terminator.
   CALL PRINT_TABLE_STOP

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

   ! . Deactivate the trajectories.
   CALL DCD_DEACTIVATE ( DCD_IN  )
   CALL DCD_DEACTIVATE ( DCD_OUT )

   ! . Reset the coordinates.
   ATMCRD = TMPCRD

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( D1, D2, ECHAIN, G, G_CONSTRAINT, TMPCRD, X, XBEG, XEND )

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
         DOTFAC = SQRT ( DOT_PRODUCT ( G_CONSTRAINT(1:N3,IVEC), G_CONSTRAINT(1:N3,IVEC) ) )
         G_CONSTRAINT(1:N3,IVEC) = G_CONSTRAINT(1:N3,IVEC) / DOTFAC

      END DO

      END SUBROUTINE GET_LINEAR_CONSTRAINTS

   END SUBROUTINE CHAIN_OPTIMIZE

   !----------------------------
   SUBROUTINE EGCALC ( X, F, G )
   !----------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT) :: F

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)  :: X
   REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: G

   ! . Local scalars.
   INTEGER            :: I, II, J, JJ
   REAL ( KIND = DP ) :: DFAC, D1AVE, D1FAC

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:NCHAIN)     :: D2FAC
   REAL ( KIND = DP ), DIMENSION(1:N3)         :: DX

   ! . Initialize the function value and its derivatives.
   F              = 0.0_DP
   G(1:NCHAIN*N3) = 0.0_DP

   ! . Loop over the structures in the chain.
   DO I = 1,NCHAIN

      ! . Calculate the index into the coordinate array.
      II = ( I - 1 ) * N3

      ! . Copy the variables to the coordinate array.
      CALL X_TO_CRD ( X(II+1:II+N3), ATMCRD )

      ! . Calculate the energy and derivatives.
      CALL GRADIENT ( PRINT = .FALSE. )

      ! . Save the energy.
      ECHAIN(I+1) = ETOTAL

      ! . Add in the energy of the structure to the function value.
      F = F + ETOTAL

      ! . Copy over the derivatives.
      CALL CRD_TO_X ( ATMDER, G(II+1:II+N3) )

   END DO

   ! . Calculate the distances between neighbouring structures.
   D1(1) = DISTANCE ( XBEG, X(1:N3) )
   DO I = 2,NCHAIN
      II = ( I - 2 ) * N3
      JJ = II + N3
      D1(I) = DISTANCE ( X(II+1:II+N3), X(JJ+1:JJ+N3) )
   END DO
   II = ( NCHAIN - 1 ) * N3
   D1(NCHAIN+1) = DISTANCE ( X(II+1:II+N3), XEND )

   ! . Calculate the average neighbouring distance.
   D1AVE = SUM ( D1(1:NCHAIN+1) ) / REAL ( ( NCHAIN + 1 ), DP )

   ! . Calculate the CHAIN terms.
   IF ( GAMFAC /= 0.0_DP ) THEN

      ! . Calculate the CHAIN constraint contribution to the function value.
      F = F + GAMFAC * SUM ( ( D1(1:NCHAIN+1) - D1AVE )**2 )

      ! . Calculate the CHAIN constraint gradient terms.
      D1FAC   = 2.0_DP * GAMFAC * ( 1.0_DP - D1AVE / D1(1) )
      G(1:N3) = G(1:N3) - D1FAC * ( XBEG - X(1:N3) )
      DO I = 2,NCHAIN
         II    = ( I - 2 ) * N3
         JJ    = II + N3
         D1FAC = 2.0_DP * GAMFAC * ( 1.0_DP - D1AVE / D1(I) )
         DX    = D1FAC * ( X(II+1:II+N3) - X(JJ+1:JJ+N3) )
         G(II+1:II+N3) = G(II+1:II+N3) + DX
         G(JJ+1:JJ+N3) = G(JJ+1:JJ+N3) - DX
      END DO
      II    = ( NCHAIN - 1 ) * N3
      D1FAC = 2.0_DP * GAMFAC * ( 1.0_DP - D1AVE / D1(NCHAIN+1) )
      G(II+1:II+N3) = G(II+1:II+N3) + D1FAC * ( X(II+1:II+N3) - XEND )

   END IF

   ! . Calculate the distances squared between next nearest neighbours.
   II = N3
   D2(1) = DISTANCE_SQUARED ( XBEG, X(II+1:II+N3) )
   DO I = 2,NCHAIN-1
      II = ( I - 2 ) * N3
      JJ = II + 2 * N3
      D2(I) = DISTANCE_SQUARED ( X(II+1:II+N3), X(JJ+1:JJ+N3) )
   END DO
   II = ( NCHAIN - 2 ) * N3
   D2(NCHAIN) = DISTANCE_SQUARED ( X(II+1:II+N3), XEND )

   ! . Calculate the REPULSION terms.
   IF ( RHOFAC /= 0.0_DP ) THEN

      ! . Calculate some intermediate factors for the REPULSION constraint calculation.
      D2FAC(1:NCHAIN) = RHOFAC * EXP ( - ( LAMFAC / ( D1AVE * D1AVE ) ) * D2(1:NCHAIN) )

      ! . Calculate the REPULSION constraint contribution to the function value.
      F = F + SUM ( D2FAC(1:NCHAIN) ) / LAMFAC

      ! . Scale D2FAC for the gradient calculation.
      D2FAC(1:NCHAIN) = - 2.0_DP * D2FAC(1:NCHAIN) / ( D1AVE * D1AVE )

      ! . Calculate the REPULSION constraint gradient terms.
      II = N3
      G(II+1:II+N3) = G(II+1:II+N3) - D2FAC(1) * ( XBEG - X(II+1:II+N3) )
      DO I = 2,NCHAIN-1
         II = ( I - 2 ) * N3
         JJ = II + 2 * N3
         DX = D2FAC(I) * ( X(II+1:II+N3) - X(JJ+1:JJ+N3) )
         G(II+1:II+N3) = G(II+1:II+N3) + DX
         G(JJ+1:JJ+N3) = G(JJ+1:JJ+N3) - DX
      END DO
      II = ( NCHAIN - 2 ) * N3
      G(II+1:II+N3) = G(II+1:II+N3) + D2FAC(NCHAIN) * ( X(II+1:II+N3) - XEND )

      ! . Calculate the factor for the calculation of the remaining gradient terms.
      DFAC = - DOT_PRODUCT ( D2FAC(1:NCHAIN), D2(1:NCHAIN) ) / ( D1AVE * REAL ( ( NCHAIN + 1 ), DP ) )

      ! . Calculate the REPULSION constraint gradient terms due to the average distance factor.
      G(1:N3) = G(1:N3) - ( DFAC / D1(1) ) * ( XBEG - X(1:N3) )
      DO I = 2,NCHAIN
         II = ( I - 2 ) * N3
         JJ = II + N3
         DX = ( DFAC / D1(I) ) * ( X(II+1:II+N3) - X(JJ+1:JJ+N3) )
         G(II+1:II+N3) = G(II+1:II+N3) + DX
         G(JJ+1:JJ+N3) = G(JJ+1:JJ+N3) - DX
      END DO
      II = ( NCHAIN - 1 ) * N3
      G(II+1:II+N3) = G(II+1:II+N3) + ( DFAC / D1(NCHAIN+1) ) * ( X(II+1:II+N3) - XEND )

   END IF

   ! . Project out the rigid body motions from the gradient.
   IF ( NCONSTRAINTS > 0 ) THEN
      DO I = 1,NCHAIN
         II = ( I - 1 ) * N3
         DO J = 1,NCONSTRAINTS
            G(II+1:II+N3) = G(II+1:II+N3) - DOT_PRODUCT ( G_CONSTRAINT(1:N3,J), G(II+1:II+N3) ) * G_CONSTRAINT(1:N3,J)
         END DO
      END DO
   END IF

   !===========================================================================
   CONTAINS
   !===========================================================================

      !---------------------------
      FUNCTION DISTANCE ( X1, X2 )
      !---------------------------

      ! . Function declarations.
      REAL ( KIND = DP ) :: DISTANCE

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: X1, X2

      ! . Calculate the distance.
      DX = X1 - X2
      DISTANCE = SQRT ( DOT_PRODUCT ( DX, DX ) )

      END FUNCTION DISTANCE

      !-----------------------------------
      FUNCTION DISTANCE_SQUARED ( X1, X2 )
      !-----------------------------------

      ! . Function declarations.
      REAL ( KIND = DP ) :: DISTANCE_SQUARED

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: X1, X2

      ! . Calculate the distance.
      DX = X1 - X2
      DISTANCE_SQUARED = DOT_PRODUCT ( DX, DX )

      END FUNCTION DISTANCE_SQUARED

   END SUBROUTINE EGCALC

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

END MODULE SELF_AVOIDING_WALK
