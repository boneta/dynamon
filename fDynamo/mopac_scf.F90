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
!                             The Quantum SCF Module
!===============================================================================
!
! . Scalars:
!
!   ACURCY		            the accuracy tolerance.
!   DAMPF		            the damping factor.
!   DAMP0		            the zero damping factor.
!   IEXTRP		            the extrapolation frequency.
!   MAXIT		            the maximum number of SCF iterations.
!   NDIIS                           the number of DIIS matrices.
!   OSHIFT                          the open-shell level-shift parameter.
!   QFORCEUHF                       the flag to always force a UHF calculation.
!   VSHIFT		            the virtual-orbital level-shift parameter.
!
! . Principal subroutines:
!
!   MOPAC_SCF_CALCULATE             perform an SCF calculation.
!   MOPAC_SCF_INITIALIZE            initialize the SCF options.
!   MOPAC_SCF_OPTIONS               set the SCF options.
!
! . Private functions and subroutines:
!
!   DAVIDSON_DAMPING                damp a Fock matrix.
!   DIIS                            perform a DIIS extrapolation.
!   GET_DAMP_FACTOR                 determine the damp factor for an iteration.
!   SPIN_VALUES                     calculate Sz and S^2 for a wavefunction.
!   
!   SYMMETRIC_MATRIX_DIAGONAL       get the diagonal of a symmetric matrix.
!   SYMMETRIC_MATRIX_DOT_PRODUCT    the scalar product of two matrices.
!   SYMMETRIC_MATRIX_MULTIPLY_2     multiply two matrices.
!
!===============================================================================
MODULE MOPAC_SCF

! . Module declarations.
USE DEFINITIONS,       ONLY : DP
USE PRINTING,          ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
                              PRINT_SUMMARY_START, PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT,          &
			      PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP,              &
                              PRINT_PARAGRAPH
USE DIAGONALIZATION,   ONLY : SYMMETRIC_UPPER
USE LINEAR_EQUATIONS,  ONLY : LINEAR_EQUATIONS_SOLVE
USE MOPAC_DATA,        ONLY : BFINDEX, DENMATA, DENMATB, DENMAT, HCORE, NALPHA, NBASIS, NBASTR, NBETA, NPIBEADS
USE MOPAC_DENSITY,     ONLY : DENSITY_CALCULATE
USE MOPAC_FOCK_MATRIX, ONLY : FOCK_MATRIX_RESTRICTED, FOCK_MATRIX_UNRESTRICTED
USE MOPAC_INTEGRALS,   ONLY : INTEGRALS_HCORE

IMPLICIT NONE
PRIVATE
PUBLIC :: MOPAC_SCF_CALCULATE, MOPAC_SCF_OPTIONS, QSCFBOMB, MOPAC_SCF_INITIALIZE
#ifndef PGPC
SAVE
#endif

! . Module scalars.
INTEGER            :: IEXTRP    = 15,        MAXIT = 100,    NDIIS = 19,     NDIISP = 20 ! = NDIIS + 1
LOGICAL            :: QFORCEUHF = .FALSE., QSCFBOMB = .TRUE.
REAL ( KIND = DP ) :: ACURCY    = 1.0E-8_DP, DAMPF = 0.0_DP, DAMP0 = 0.0_DP, SHIFTO = 0.0_DP, SHIFTV = 0.0_DP

! . DIIS variables.
INTEGER                                         :: MATNUM
INTEGER,            ALLOCATABLE, DIMENSION(:)   :: MATIND
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: BCOEFF

!===============================================================================
CONTAINS
!===============================================================================

   !---------------------------------------------
   SUBROUTINE MOPAC_SCF_CALCULATE ( EHF, QPRINT )
   !---------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN)  :: QPRINT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NPIBEADS), INTENT(OUT) :: EHF

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: CVGTOL = 0.005_DP, DMPTOL = 1.0_DP, SHFTOL = 0.2_DP

   ! . Other local scalars.
   INTEGER            :: FLAG, IBEAD, ICALL, IDAMP, ITER
   LOGICAL            :: CVGED, DMPFLG, LEVSHF, QBETA, QOPEN, QUHF, XTPFLG
   REAL ( KIND = DP ) :: DAMP, DEHF, DEHFOLD, DIFF, DMPSAV, EHF0, OCCNUM, OSHIFT, VSHIFT

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: FOCK1, FOCK2, ORBENE1, ORBENE2, TMPMAT
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: ERRMAT1, ERRMAT2, FCKDMP1, FCKDMP2, FCKOLD1, FCKOLD2, &
                                                      MORBS1, MORBS2

   !----------------------------------------------------------------------------
   ! . Set up the calculation.
   !----------------------------------------------------------------------------
   ! . Initialize the energy.
   EHF = 0.0_DP

   ! . Check the number of basis functions and electrons.
   IF ( ( NBASIS <= 0 ) .OR. ( ( NALPHA <= 0 ) .AND. ( NBETA <= 0 ) ) ) RETURN

   ! . Set the UHF flag.
   QUHF = QFORCEUHF .OR. ( NALPHA /= NBETA )

   ! . Set the flag for the existence of BETA orbitals.
   QBETA = QUHF .AND. ( NBETA > 0 )

   ! . Set the flag for the existence of open shell orbitals.
   QOPEN = ( NALPHA /= NBETA )

   ! . Check to see if the density matrices exist and are in the correct place.
   ! . Spin-unrestricted calculation.
   IF ( QUHF ) THEN

      ! . Check QFORCEUHF.
      IF ( QFORCEUHF ) THEN

         ! . Allocate DENMATA and DENMATB if they do not exist.
	 IF ( .NOT. ALLOCATED ( DENMATA ) ) THEN
	    ALLOCATE ( DENMATA(1:NBASTR,1:NPIBEADS) ) ; DENMATA = 0.5_DP * DENMAT
	 END IF
	 IF ( .NOT. ALLOCATED ( DENMATB ) ) THEN
	    ALLOCATE ( DENMATB(1:NBASTR,1:NPIBEADS) ) ; DENMATB = 0.5_DP * DENMAT
	 END IF

      END IF

      ! . There is nothing to do here except zero DENMATB if NBETA == 0.
      IF ( NBETA == 0 ) THEN
         DENMAT  = DENMATA
         DENMATB = 0.0_DP
      END IF

   ! . Spin-restricted calculation.
   ELSE

      ! . Copy DENMAT to DENMATA.
      ALLOCATE ( DENMATA(1:NBASTR,1:NPIBEADS) ) ; DENMATA = DENMAT ; DEALLOCATE ( DENMAT )

   END IF

   ! . Reallocate the DIIS arrays.
   IF ( ALLOCATED ( MATIND ) ) DEALLOCATE ( MATIND )
   IF ( ALLOCATED ( BCOEFF ) ) DEALLOCATE ( BCOEFF )
   ALLOCATE ( BCOEFF(1:NDIIS,1:NDIIS), MATIND(1:NDIIS) )

   !----------------------------------------------------------------------------
   ! . Loop over the SCFs.
   !----------------------------------------------------------------------------
   ! . Loop over the number of PI beads.
   DO IBEAD = 1,NPIBEADS

      ! . Set up the one-electron matrix.
      CALL INTEGRALS_HCORE  ( IBEAD )

      ! . Initialize some scalars.
      ICALL = 0
      IDAMP = 0
      DEHF  = 0.0_DP
      DIFF  = 0.0_DP
      EHF0  = 0.0_DP

      ! . Set the initial values for the damping and level shift parameters.
      DAMP   = DAMPF
      OSHIFT = SHIFTO
      VSHIFT = SHIFTV

      ! . Set the level shifting flag.
      IF ( QOPEN ) THEN
	 LEVSHF = ( OSHIFT > 0.0_DP ) .OR. ( VSHIFT > 0.0_DP )
      ELSE
	 LEVSHF = ( VSHIFT > 0.0_DP )
      END IF

      ! . Set the damping flag.
      DMPFLG = ( DAMP > DMPTOL )

      ! . Set the extrapolation flag.
      XTPFLG = ( .NOT. DMPFLG ) .AND. ( .NOT. LEVSHF )

      ! . Set the occupation number for the calculation of the density matrix.
      IF ( QUHF ) THEN
	 OCCNUM = 1.0_DP
      ELSE
	 OCCNUM = 2.0_DP
      END IF

      !----------------------------------------------------------------------------
      ! . Allocate temporary arrays.
      !----------------------------------------------------------------------------
      ! . Allocate the arrays for orbital set 1.
      ALLOCATE ( FOCK1(1:NBASTR), MORBS1(1:NBASIS,1:NBASIS), ORBENE1(1:NBASIS), &
                      ERRMAT1(1:NBASTR,1:NDIIS), FCKDMP1(1:NBASTR,1:2), FCKOLD1(1:NBASTR,1:NDIIS), STAT = FLAG )

      ! . Allocate the arrays for orbital set 2.
      IF ( QUHF ) THEN
	 ALLOCATE ( FOCK2(1:NBASTR), MORBS2(1:NBASIS,1:NBASIS), ORBENE2(1:NBASIS), &
                      ERRMAT2(1:NBASTR,1:NDIIS), FCKDMP2(1:NBASTR,1:2), FCKOLD2(1:NBASTR,1:NDIIS), STAT = FLAG )
      END IF

      ! . Allocate extra temporary arrays.
      ALLOCATE ( TMPMAT(1:NBASTR), STAT = FLAG )

      ! . Check FLAG.
      IF ( FLAG /= 0 ) CALL PRINT_ERROR ( "MOPAC_SCF_CALCULATE", "Unable to allocate sufficient array space." )

      !----------------------------------------------------------------------------
      ! . Do some printing.
      !----------------------------------------------------------------------------
      ! . Write out the header information.
      IF ( QPRINT ) THEN

	 ! . The header.
         CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800" )
	 IF ( QUHF ) THEN
	    CALL PRINT_SUMMARY_START ( "Spin-Unrestricted SCF Calculation" )
	 ELSE
	    CALL PRINT_SUMMARY_START ( "Spin-Restricted SCF Calculation" )
	 END IF

	 ! . Write out the SCF parameters.
	 WRITE ( PRINT_LINE, "(I14)"   ) MAXIT     ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations"   )
	 WRITE ( PRINT_LINE, "(G14.8)" ) ACURCY    ; CALL PRINT_SUMMARY_ELEMENT ( "Accuracy Tolerance"     )
	 WRITE ( PRINT_LINE, "(L14)"   ) QFORCEUHF ; CALL PRINT_SUMMARY_ELEMENT ( "Force UHF"              )
	 WRITE ( PRINT_LINE, "(I14)"   ) IEXTRP    ; CALL PRINT_SUMMARY_ELEMENT ( "Extrapolation Freq."    )
	 WRITE ( PRINT_LINE, "(G14.8)" ) DAMPF     ; CALL PRINT_SUMMARY_ELEMENT ( "Damping Factor"         )
	 WRITE ( PRINT_LINE, "(G14.8)" ) DAMP0     ; CALL PRINT_SUMMARY_ELEMENT ( "Zero Damping Factor"    )
	 WRITE ( PRINT_LINE, "(G14.8)" ) OSHIFT    ; CALL PRINT_SUMMARY_ELEMENT ( "Open-Shell Level Shift" )
	 WRITE ( PRINT_LINE, "(G14.8)" ) VSHIFT    ; CALL PRINT_SUMMARY_ELEMENT ( "Virtual Level Shift"    )
	 WRITE ( PRINT_LINE, "(I14)"   ) NDIIS     ; CALL PRINT_SUMMARY_ELEMENT ( "No. of DIIS Matrices"   )
	 CALL PRINT_SUMMARY_STOP

	 ! . Write out the column headers for the cycle information.
	 CALL PRINT_TABLE_OPTIONS ( COLUMNS = 5, HEADER_COLOR = "#FF8800", VARIABLEWIDTHS = (/ 8, 18, 18, 18, 18 /) )
	 CALL PRINT_TABLE_START
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Cycle",          HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Energy",         HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Energy Change",  HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Density Change", HEADER = .TRUE. )
	 CALL PRINT_TABLE_ELEMENT ( TEXT = "Damp/Shift",     HEADER = .TRUE. )

      END IF

    ! . Perform 0SCF calculation
    IF( MAXIT <= 0 ) THEN
! #######################################################################################################################
	 ! . Build the two-electron part of the Fock matrix.
	 IF ( QUHF ) THEN
            CALL FOCK_MATRIX_UNRESTRICTED ( IBEAD, DENMAT(:,IBEAD), DENMATA(:,IBEAD), DENMATB(:,IBEAD), FOCK1, FOCK2 )
	 ELSE
            CALL FOCK_MATRIX_RESTRICTED ( IBEAD, DENMATA(:,IBEAD), FOCK1 )
	 END IF

	 ! . Produce the full FOCK matrices by adding in the one-electron matrix.
	 FOCK1 = FOCK1 + HCORE
	 IF ( QBETA ) FOCK2 = FOCK2 + HCORE

	 ! . Calculate the HF energy.
	 EHF0        = EHF(IBEAD)
	 EHF(IBEAD)  = 0.5_DP * ( SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA(:,IBEAD), HCORE ) + &
                                  SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA(:,IBEAD), FOCK1 ) )
	 IF ( QBETA ) THEN
            EHF(IBEAD) = EHF(IBEAD) + 0.5_DP * ( SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATB(:,IBEAD), HCORE ) + &
	                                         SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATB(:,IBEAD), FOCK2 ) )
	 END IF

	 ! . Calculate the energy differences.
	 DEHFOLD = DEHF
	 DEHF    = EHF(IBEAD) - EHF0

	 ! . Print out data for the current iteration.
	 IF ( QPRINT ) THEN
	    WRITE ( PRINT_LINE, "(I8)"    ) -1         ; CALL PRINT_TABLE_ELEMENT
	    WRITE ( PRINT_LINE, "(G18.8)" ) EHF(IBEAD) ; CALL PRINT_TABLE_ELEMENT
	    WRITE ( PRINT_LINE, "(G18.8)" ) DEHF       ; CALL PRINT_TABLE_ELEMENT
	    WRITE ( PRINT_LINE, "(G18.8)" ) DIFF       ; CALL PRINT_TABLE_ELEMENT
            IF ( DMPFLG ) THEN
	       WRITE ( PRINT_LINE, "(G18.8)" ) DAMP
	    ELSE IF ( LEVSHF ) THEN
	       WRITE ( PRINT_LINE, "(G18.8)" ) MAX ( OSHIFT, VSHIFT )
	    ELSE
	       WRITE ( PRINT_LINE, "(G18.8)" ) 0.0_DP
	    END IF
	    CALL PRINT_TABLE_ELEMENT
	 END IF

	 ! . Set the convergence flag.
	 CVGED = .TRUE.
! #######################################################################################################################
    END IF

      !----------------------------------------------------------------------------
      ! . Do the SCF.
      !----------------------------------------------------------------------------
      ! . Start the SCF procedure.
      DO ITER = 1,MAXIT

	 ! . Build the two-electron part of the Fock matrix.
	 IF ( QUHF ) THEN
            CALL FOCK_MATRIX_UNRESTRICTED ( IBEAD, DENMAT(:,IBEAD), DENMATA(:,IBEAD), DENMATB(:,IBEAD), FOCK1, FOCK2 )
	 ELSE
            CALL FOCK_MATRIX_RESTRICTED ( IBEAD, DENMATA(:,IBEAD), FOCK1 )
	 END IF

	 ! . Produce the full FOCK matrices by adding in the one-electron matrix.
	 FOCK1 = FOCK1 + HCORE
	 IF ( QBETA ) FOCK2 = FOCK2 + HCORE

	 ! . Calculate the HF energy.
	 EHF0        = EHF(IBEAD)
	 EHF(IBEAD)  = 0.5_DP * ( SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA(:,IBEAD), HCORE ) + &
                                  SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATA(:,IBEAD), FOCK1 ) )
	 IF ( QBETA ) THEN
            EHF(IBEAD) = EHF(IBEAD) + 0.5_DP * ( SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATB(:,IBEAD), HCORE ) + &
	                                         SYMMETRIC_MATRIX_DOT_PRODUCT ( DENMATB(:,IBEAD), FOCK2 ) )
	 END IF

	 ! . Calculate the energy differences.
	 DEHFOLD = DEHF
	 DEHF    = EHF(IBEAD) - EHF0

	 ! . Print out data for the current iteration.
	 IF ( QPRINT ) THEN
	    WRITE ( PRINT_LINE, "(I8)"    ) ITER       ; CALL PRINT_TABLE_ELEMENT
	    WRITE ( PRINT_LINE, "(G18.8)" ) EHF(IBEAD) ; CALL PRINT_TABLE_ELEMENT
	    WRITE ( PRINT_LINE, "(G18.8)" ) DEHF       ; CALL PRINT_TABLE_ELEMENT
	    WRITE ( PRINT_LINE, "(G18.8)" ) DIFF       ; CALL PRINT_TABLE_ELEMENT
            IF ( DMPFLG ) THEN
	       WRITE ( PRINT_LINE, "(G18.8)" ) DAMP
	    ELSE IF ( LEVSHF ) THEN
	       WRITE ( PRINT_LINE, "(G18.8)" ) MAX ( OSHIFT, VSHIFT )
	    ELSE
	       WRITE ( PRINT_LINE, "(G18.8)" ) 0.0_DP
	    END IF
	    CALL PRINT_TABLE_ELEMENT
	 END IF

	 ! . Set the convergence flag.
	 CVGED = ( DEHF < CVGTOL ) .AND. ( DIFF < ACURCY ) .AND. ( ITER > 1 )

	 ! . Exit if there is convergence.
	 IF ( CVGED ) EXIT

	 !-------------------------------------------------------------------------
	 ! . Modify the Fock matrices to aid convergence.
	 !-------------------------------------------------------------------------
	 ! . Check to see if extrapolation is to be used.
	 IF ( DMPFLG .OR. LEVSHF ) THEN
            XTPFLG = ( ( ABS ( DEHF ) < CVGTOL ) .AND. ( DAMP < DMPTOL ) .AND. ( MAX ( OSHIFT, VSHIFT ) < SHFTOL ) ) .OR. &
                     ( ( ABS ( DEHF ) < CVGTOL ) .AND. ( DIFF < ACURCY ) )
            IF ( XTPFLG ) THEN
               DMPFLG = .FALSE.
               LEVSHF = .FALSE.
	       OSHIFT = 0.0_DP
               VSHIFT = 0.0_DP
            END IF
	 ! . Check to see if damping is to be switched back on.
	 ELSE
            IF ( DEHF >= CVGTOL ) THEN
               DAMP   = DMPTOL
               DMPFLG = .TRUE.
	       LEVSHF = .FALSE.
	       IDAMP  = 0
               XTPFLG = .FALSE.
            END IF
	 END IF

	 ! . Level shifting after first iteration only.
	 IF ( LEVSHF .AND. ( ITER >  1 ) ) THEN

            ! . There are open shell orbitals.
	    IF ( QOPEN ) THEN
               CALL DENSITY_CALCULATE ( TMPMAT, MORBS1(1:NBASIS,NBETA+1:NALPHA), NALPHA-NBETA, 1.0_DP )
               FOCK1 = FOCK1 + OSHIFT * TMPMAT
            END IF

            ! . Modify FOCK for the virtual orbitals.
            CALL DENSITY_CALCULATE ( TMPMAT, MORBS1(1:NBASIS,NALPHA+1:NBASIS), NBASIS-NALPHA, 1.0_DP )
            FOCK1 = FOCK1 + VSHIFT * TMPMAT

            ! . Beta orbitals.
	    IF ( QBETA ) THEN
               CALL DENSITY_CALCULATE ( TMPMAT, MORBS2(1:NBASIS,NALPHA+1:NBASIS), NBASIS-NBETA, 1.0_DP )
               FOCK2 = FOCK2 + VSHIFT * TMPMAT
            END IF

	 ! . Damping.
	 ELSE IF ( DMPFLG ) THEN

            ! . Get the damp factor after the second iteration.
	    IF ( ITER > 2 ) CALL GET_DAMP_FACTOR ( ITER, DEHF, DEHFOLD, DAMP )

            ! . Do the damping.
            IDAMP = IDAMP + 1
            CALL DAVIDSON_DAMPING ( IDAMP, IEXTRP, DAMP, DAMP0, DMPSAV, FOCK1, FCKDMP1(:,1), FCKDMP1(:,2) )
	    IF ( QBETA ) CALL DAVIDSON_DAMPING ( IDAMP, IEXTRP, DAMP, DAMP0, DMPSAV, FOCK2, FCKDMP2(:,1), FCKDMP2(:,2) )

	 ! . Extrapolation.
	 ELSE IF ( XTPFLG ) THEN
            IF ( QBETA ) THEN
               CALL DIIS ( ICALL, DENMATA(:,IBEAD), ERRMAT1, FOCK1, FCKOLD1, &
	                	  DENMATB(:,IBEAD), ERRMAT2, FOCK2, FCKOLD2  )
	    ELSE
               CALL DIIS ( ICALL, DENMATA(:,IBEAD), ERRMAT1, FOCK1, FCKOLD1 )
	    END IF
	 END IF

	 !-------------------------------------------------------------------------
	 ! . First orbital set.
	 !-------------------------------------------------------------------------
         ! . Diagonalize the Fock matrix.
         CALL SYMMETRIC_UPPER ( FOCK1, ORBENE1, MORBS1 )

	 ! . Remove the level shift from the virtual orbital energies.
	 IF ( LEVSHF .AND. ( ITER > 1 ) ) THEN
	    IF ( QOPEN ) ORBENE1(NBETA+1:NALPHA) = ORBENE1(NBETA+1:NALPHA) - OSHIFT
            ORBENE1(NALPHA+1:NBASIS) = ORBENE1(NALPHA+1:NBASIS) - VSHIFT
	 END IF

	 ! . Save the old density matrix temporarily.
	 TMPMAT = DENMATA(:,IBEAD)

	 ! . Form the density matrix.
	 CALL DENSITY_CALCULATE ( DENMATA(:,IBEAD), MORBS1(1:NBASIS,1:NALPHA), NALPHA, OCCNUM )

	 ! . Check to see if the density matrix has converged.
	 DIFF = MAXVAL ( ABS ( DENMATA(:,IBEAD) - TMPMAT ) )

         ! . Copy DENMATA to DENMAT.
	 IF ( QUHF ) DENMAT(:,IBEAD) = DENMATA(:,IBEAD)

	 !-------------------------------------------------------------------------
	 ! . Second orbital set.
	 !-------------------------------------------------------------------------
	 IF ( QBETA ) THEN

            ! . Diagonalize the Fock matrix.
            CALL SYMMETRIC_UPPER ( FOCK2, ORBENE2, MORBS2 )

	    ! . Remove the level shift from the virtual orbital energies.
	    IF ( LEVSHF .AND. ( ITER > 1 ) ) ORBENE2(NBETA+1:NBASIS) = ORBENE2(NBETA+1:NBASIS) - VSHIFT

            ! . Save the old density matrix temporarily.
            TMPMAT = DENMATB(:,IBEAD)

	    ! . Form the beta-density matrix.
            CALL DENSITY_CALCULATE ( DENMATB(:,IBEAD), MORBS2(1:NBASIS,1:NBETA), NBETA, 1.0_DP )

	    ! . Check to see if the density matrix has converged.
	    DIFF = MAX ( DIFF, MAXVAL ( ABS ( DENMATB(:,IBEAD) - TMPMAT ) ) )

            ! . Add in DENMATB to DENMAT.
	    DENMAT(:,IBEAD) = DENMAT(:,IBEAD) + DENMATB(:,IBEAD)

	 END IF

	 ! . Check on convergence.
	 IF ( CVGED ) EXIT

      END DO

      !----------------------------------------------------------------------------
      ! . Finish up.
      !----------------------------------------------------------------------------
      ! . Finish up the table.
      IF ( QPRINT ) CALL PRINT_TABLE_STOP

      ! . Fail if the SCF has not converged.
      IF ( .NOT. CVGED ) THEN
          IF ( QSCFBOMB ) THEN
              CALL PRINT_ERROR ( "MOPAC_SCF", "Excessive number of SCF iterations." )
          ELSE
              CALL PRINT_PARAGRAPH ( TEXT = "Excessive number of SCF iterations." )
          END IF
      END IF

      ! . Deallocate most of the temporary arrays.
      DEALLOCATE ( FOCK1, ERRMAT1, FCKDMP1, FCKOLD1 )
      IF ( QBETA ) DEALLOCATE ( FOCK2, ERRMAT2, FCKDMP2, FCKOLD2 )

      ! . Calculate and print the spin values for a UHF wavefunction.
      IF ( QPRINT .AND. QUHF ) CALL SPIN_VALUES ( DENMATA(:,IBEAD), DENMATB(:,IBEAD) )

      ! . Deallocate the remaining temporary arrays.
      DEALLOCATE ( MORBS1, ORBENE1, TMPMAT )
      IF ( QBETA ) DEALLOCATE ( MORBS2, ORBENE2 )

   END DO

   ! . Make sure the density matrices are in the correct place for a spin-restricted calculation.
   IF ( .NOT. QUHF ) THEN

      ! . Copy DENMATA to DENMAT.
      ALLOCATE ( DENMAT(1:NBASTR,1:NPIBEADS) ) ; DENMAT = DENMATA ; DEALLOCATE ( DENMATA )

   END IF

   END SUBROUTINE MOPAC_SCF_CALCULATE

   !------------------------------
   SUBROUTINE MOPAC_SCF_INITIALIZE
   !------------------------------

   ! . Initialize the SCF parameters.
   IEXTRP    =  15
   MAXIT     = 100
   NDIIS     =  19
   NDIISP    = NDIIS + 1

   QFORCEUHF = .FALSE.

   ACURCY    = 1.0E-8_DP
   DAMPF     = 0.0_DP
   DAMP0     = 0.0_DP
   SHIFTO    = 0.0_DP
   SHIFTV    = 0.0_DP

   END SUBROUTINE MOPAC_SCF_INITIALIZE

   !-----------------------------------------------------------------------------------------------
   SUBROUTINE MOPAC_SCF_OPTIONS ( ACCURACY, DAMP_FACTOR, DIIS_MATRICES, EXTRAPOLATION, FORCE_UHF, &
                                  ITERATIONS, OPENSHELL_SHIFT, VIRTUAL_SHIFT, ZERO_DAMP_FACTOR,   &
				  PRINT )
   !-----------------------------------------------------------------------------------------------

   ! . Optional argument declarations.
   INTEGER,            INTENT(IN), OPTIONAL :: DIIS_MATRICES, EXTRAPOLATION, ITERATIONS
   LOGICAL,            INTENT(IN), OPTIONAL :: FORCE_UHF, PRINT
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: ACCURACY, DAMP_FACTOR, OPENSHELL_SHIFT, VIRTUAL_SHIFT, ZERO_DAMP_FACTOR

   ! . Local scalars.
   LOGICAL :: QPRINT

   ! . Check each of the arguments.
   IF ( PRESENT ( ACCURACY         ) ) ACURCY	 = ACCURACY
   IF ( PRESENT ( DAMP_FACTOR      ) ) DAMPF	 = DAMP_FACTOR
   IF ( PRESENT ( DIIS_MATRICES    ) ) NDIIS	 = DIIS_MATRICES
   IF ( PRESENT ( EXTRAPOLATION    ) ) IEXTRP	 = EXTRAPOLATION
   IF ( PRESENT ( FORCE_UHF        ) ) QFORCEUHF  = FORCE_UHF
   IF ( PRESENT ( ITERATIONS       ) ) MAXIT	 = ITERATIONS
   IF ( PRESENT ( OPENSHELL_SHIFT  ) ) SHIFTO	 = OPENSHELL_SHIFT
   IF ( PRESENT ( VIRTUAL_SHIFT    ) ) SHIFTV	 = VIRTUAL_SHIFT
   IF ( PRESENT ( ZERO_DAMP_FACTOR ) ) DAMP0	 = ZERO_DAMP_FACTOR   
 
   ! . Set NDIISP.
   NDIISP = NDIIS + 1
 
   ! . Check the PRINT argument.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Check and reset the values of the level shift parameters.
   IF ( ( SHIFTO >  0.0_DP ) .AND. ( SHIFTV <= 0.0_DP ) ) SHIFTV = SHIFTO
   IF ( ( SHIFTO <= 0.0_DP ) .AND. ( SHIFTV >  0.0_DP ) ) SHIFTO = SHIFTV

   ! . Print the parameters.
   IF ( QPRINT ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800" )
      CALL PRINT_SUMMARY_START ( "Quantum SCF Options" )
      WRITE ( PRINT_LINE, "(I14)"   ) MAXIT     ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Iterations"   )
      WRITE ( PRINT_LINE, "(G14.8)" ) ACURCY    ; CALL PRINT_SUMMARY_ELEMENT ( "Accuracy Tolerance"     )
      WRITE ( PRINT_LINE, "(L14)"   ) QFORCEUHF ; CALL PRINT_SUMMARY_ELEMENT ( "Force UHF"              )
      WRITE ( PRINT_LINE, "(I14)"   ) IEXTRP    ; CALL PRINT_SUMMARY_ELEMENT ( "Extrapolation Freq."    )
      WRITE ( PRINT_LINE, "(G14.8)" ) DAMPF     ; CALL PRINT_SUMMARY_ELEMENT ( "Damping Factor"         )
      WRITE ( PRINT_LINE, "(G14.8)" ) DAMP0     ; CALL PRINT_SUMMARY_ELEMENT ( "Zero Damping Factor"    )
      WRITE ( PRINT_LINE, "(G14.8)" ) SHIFTO    ; CALL PRINT_SUMMARY_ELEMENT ( "Open-Shell Level Shift" )
      WRITE ( PRINT_LINE, "(G14.8)" ) SHIFTV    ; CALL PRINT_SUMMARY_ELEMENT ( "Virtual Level Shift"    )
      WRITE ( PRINT_LINE, "(I14)"   ) NDIIS     ; CALL PRINT_SUMMARY_ELEMENT ( "No. of DIIS Matrices"   )
      CALL PRINT_SUMMARY_STOP
   END IF

   END SUBROUTINE MOPAC_SCF_OPTIONS

!-------------------------------------------------------------------------------
! . Private subroutines.
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   SUBROUTINE DAVIDSON_DAMPING ( ITER, IEXTRP, DAMP, DAMP0, DMPSAV, D0, D1, D2 )
   !----------------------------------------------------------------------------

   ! . Scalar argument declarations.
   INTEGER,            INTENT(IN)    :: IEXTRP, ITER
   REAL ( KIND = DP ), INTENT(IN)    :: DAMP, DAMP0
   REAL ( KIND = DP ), INTENT(INOUT) :: DMPSAV

   ! . Array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:NBASTR), INTENT(INOUT) :: D0, D1, D2

   ! . Local scalars.
   INTEGER            :: N
   REAL ( KIND = DP ) :: CUTOFF, DUM11, DUM21, DUM22, EPS

   ! . Parameter declarations.
   REAL ( KIND = DP ), PARAMETER :: TOL = 0.01_DP

   ! . On the first iteration save the first matrix only.
   IF ( ITER <= 1 ) THEN
      D1 = D0
   ELSE
      ! . Perform the damping without extrapolation for ITER > 1.
      IF ( IEXTRP <= 0 ) THEN
	 D0 = ( D0 + DAMP * D1 ) / ( 1.0_DP + DAMP )
	 D1 = D0
      ELSE
	 ! . Perform the damping for ITER = 2.
	 IF ( ITER == 2 ) THEN
            D0 = ( D0 + DAMP * D1 ) / ( 1.0_DP + DAMP )
            D2 = D1
            D1 = D0
            DMPSAV = DAMP
	 ELSE
            ! . Check whether any extrapolation needs to be performed.
            N = ITER - ( ITER / IEXTRP ) * IEXTRP
            ! . Perform damping without extrapolation for ITER > 2.
            IF ( N /= 0 ) THEN
               D0 = ( D0  + DAMP * D1 ) / ( 1.0_DP + DAMP )
               D2 = D1
               D1 = D0
               DMPSAV = DAMP
            ELSE
               ! . Perform damping and extrapolation for ITER > 2.
               CUTOFF = 0.5_DP - DAMP0
               IF ( CUTOFF >= 1.0_DP ) CUTOFF = 0.99_DP
               IF ( CUTOFF <= 0.0_DP ) CUTOFF = TOL
               D0    = ( D0 + DMPSAV * D1 ) / ( 1.0_DP + DMPSAV )
               DUM11 = DOT_PRODUCT ( ( D0 - D1 ), ( D0 - D1 ) )
               DUM21 = DOT_PRODUCT ( ( D1 - D2 ), ( D0 - D1 ) )
               DUM22 = DOT_PRODUCT ( ( D1 - D2 ), ( D1 - D2 ) )
               EPS   = ( DUM11 - DUM21 ) / ( DUM21 - DUM22 )
               IF ( DUM21 * DUM21 < 0.5_DP * DUM11 * DUM22 ) EPS = 0.0_DP
               IF ( EPS > CUTOFF ) EPS = CUTOFF
               D0 = ( D0 - EPS * D1 ) / ( 1.0_DP - EPS )
               D2 = D1
               D1 = D0
               DMPSAV = -EPS * ( 1.0_DP + DMPSAV ) + DMPSAV
            END IF
	 END IF
      END IF
   END IF

   END SUBROUTINE DAVIDSON_DAMPING

   !--------------------------------------------------------------
   SUBROUTINE DIIS ( NCALLS, DENMATA, ERRMATA, FCKNEWA, FCKOLDA, &
                             DENMATB, ERRMATB, FCKNEWB, FCKOLDB  )
   !--------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(INOUT) :: NCALLS

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NBASTR),         INTENT(IN)    :: DENMATA
   REAL ( KIND = DP ), DIMENSION(1:NBASTR),         INTENT(INOUT) :: FCKNEWA
   REAL ( KIND = DP ), DIMENSION(1:NBASTR,1:NDIIS), INTENT(INOUT) :: ERRMATA, FCKOLDA

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NBASTR),         INTENT(IN),    OPTIONAL :: DENMATB
   REAL ( KIND = DP ), DIMENSION(1:NBASTR),         INTENT(INOUT), OPTIONAL :: FCKNEWB
   REAL ( KIND = DP ), DIMENSION(1:NBASTR,1:NDIIS), INTENT(INOUT), OPTIONAL :: ERRMATB, FCKOLDB

   ! . Local scalars.
   INTEGER :: I, IERR, II, IMAT, J, JJ, JMAT, NDIM, NEWEST, START, STOP
   LOGICAL :: QBAD

   ! . Local linear equation arrays.
   REAL ( KIND = DP ), DIMENSION(1:NDIISP,1:NDIISP) :: A
   REAL ( KIND = DP ), DIMENSION(1:NDIISP)          :: B

   ! . Increment NCALLS.
   NCALLS = NCALLS + 1

   ! . Initialization.
   IF ( NCALLS == 1 ) THEN
      ! . Initialize the matrix indexing variables.
      MATNUM = 0
      DO I = 1,NDIIS
         MATIND(I) = I
      END DO
   END IF

   ! . Update the indexing variables.
   IF ( MATNUM < NDIIS ) THEN
      MATNUM = MATNUM + 1
      NEWEST = MATIND(MATNUM)
   ELSE
      MATNUM = NDIIS
      NEWEST = MATIND(1)
      DO I = 1,(NDIIS-1)
         MATIND(I) = MATIND(I+1)
      END DO
      MATIND(NDIIS) = NEWEST
   END IF

   ! . Alpha orbitals or the full set of orbitals.
   ! . Copy the current Fock matrix to the storage array.
   FCKOLDA(1:NBASTR,NEWEST) = FCKNEWA

   ! . Calculate ( F*P - P*F ).
   CALL SYMMETRIC_MATRIX_MULTIPLY_2 ( FCKNEWA, DENMATA, ERRMATA(1,NEWEST),   1.0_DP, .TRUE.  )
   CALL SYMMETRIC_MATRIX_MULTIPLY_2 ( DENMATA, FCKNEWA, ERRMATA(1,NEWEST), - 1.0_DP, .FALSE. )

   ! . Evaluate the scalar products with the new error vector.
   DO I = 1,MATNUM
      IMAT = MATIND(I)
      BCOEFF(IMAT,NEWEST) = 2.0_DP * DOT_PRODUCT ( ERRMATA(:,IMAT), ERRMATA(:,NEWEST) )
      BCOEFF(NEWEST,IMAT) = BCOEFF(IMAT,NEWEST)
   END DO

   ! . Beta orbitals.
   IF ( PRESENT ( DENMATB ) ) THEN

      ! . Copy the current Fock matrix to the storage array.
      FCKOLDB(1:NBASTR,NEWEST) = FCKNEWB

      ! . Calculate ( F*P - P*F ).
      CALL SYMMETRIC_MATRIX_MULTIPLY_2 ( FCKNEWB, DENMATB, ERRMATB(1,NEWEST),   1.0_DP, .TRUE.  )
      CALL SYMMETRIC_MATRIX_MULTIPLY_2 ( DENMATB, FCKNEWB, ERRMATB(1,NEWEST), - 1.0_DP, .FALSE. )

      ! . Evaluate the scalar products with the new error vector and add into the existing matrix.
      ! . Note that a simple dot product is OK here as the error matrices are antisymmetric and
      ! . so the diagonal elements are zero.
      DO I = 1,MATNUM
         IMAT = MATIND(I)
         BCOEFF(IMAT,NEWEST) = BCOEFF(IMAT,NEWEST) + 2.0_DP * DOT_PRODUCT ( ERRMATB(:,IMAT), ERRMATB(:,NEWEST) )
         BCOEFF(NEWEST,IMAT) = BCOEFF(IMAT,NEWEST)
      END DO

   END IF

   ! . Restore FCKNEWA.
   FCKNEWA = FCKOLDA(1:NBASTR,NEWEST)

   ! . Stop if there are too few matrices.
   IF ( MATNUM <= 1 ) GO TO 9999

   ! . Set up and solve the linear equations.
   ! . Initialization.
   START = 1
   STOP  = MATNUM

   ! . Top of loop.
   10 CONTINUE
   NDIM = STOP - START + 2

   ! . Zero the matrix A and the RHS B.
   A = 0.0_DP
   B = 0.0_DP

   ! . Copy the B-coefficients to the matrix A.
   II = 0
   DO I = START,STOP
      II   = II + 1
      IMAT = MATIND(I)
      JJ   = 0
      DO J = START,STOP
         JJ   = JJ + 1
         JMAT = MATIND(J)
         A(II,JJ) = BCOEFF(IMAT,JMAT)
      END DO
   END DO

   ! . Fill the remaining elements of A and B.
   DO I = 1,(NDIM-1)
      A(NDIM,I) = - 1.0_DP
      A(I,NDIM) = - 1.0_DP
   END DO
   B(NDIM) = - 1.0_DP

   ! . Solve the matrix equation.
   IERR = 0
   CALL LINEAR_EQUATIONS_SOLVE ( NDIM, A(1:NDIISP,1:NDIM), NDIISP, B(1:NDIM), IERR )
   QBAD = ( IERR /= 0 )

   ! . Check for an ill-conditioned solution.
   IF ( QBAD ) THEN

      ! . Repeat the step.
      START  = START + 1
      MATNUM = STOP - START + 1
      IF ( MATNUM <= 1 ) GO TO 20
      GO TO 10

   END IF

   ! . Check that the indexing arrays are valid.
   20 CONTINUE
   IF ( START > 1 ) THEN
      ! . Reset the indexing array.
      II = 0
      DO I = START,STOP
         II = II + 1
         MATIND(II) = MATIND(I)
      END DO
      DO I = (MATNUM+1),NDIIS
         MATIND(I) = MATIND(I-1)+1
         IF ( MATIND(I) > NDIIS ) MATIND(I) = MATIND(I) - NDIIS
      END DO
      IF ( MATNUM <= 1 ) GOTO 9999
   END IF

   ! . Calculate the revised Fock matrices.
   FCKNEWA = 0.0_DP
   DO I = 1,MATNUM
      FCKNEWA = FCKNEWA + B(I) * FCKOLDA(:,MATIND(I))
   END DO
   IF ( PRESENT ( DENMATB ) ) THEN
      FCKNEWB = 0.0_DP
      DO I = 1,MATNUM
         FCKNEWB = FCKNEWB + B(I) * FCKOLDB(:,MATIND(I))
      END DO
   END IF

   ! . Exit due to too few matrices.
   9999 CONTINUE

   END SUBROUTINE DIIS

   !-------------------------------------------------
   SUBROUTINE GET_DAMP_FACTOR ( ITER, DE, DEP, DAMP )
   !-------------------------------------------------

   !. Argument declarations.
   INTEGER,            INTENT(IN)    :: ITER
   REAL ( KIND = DP ), INTENT(IN)    :: DE, DEP
   REAL ( KIND = DP ), INTENT(INOUT) :: DAMP

   ! . Parameter declarations.
   REAL ( KIND = DP ), PARAMETER :: DMPMAX = 256.0_DP, FAC = 16.0_DP

   ! . Local scalars.
   REAL ( KIND = DP ), SAVE :: DEAVG

   ! . Dtermine the average energy.
   SELECT CASE ( ITER )
   CASE ( 1 )   ; DEAVG = 0.0_DP
   CASE ( 2 )   ; DEAVG = ABS ( DE )
   CASE DEFAULT ; DEAVG = ( ABS ( DE ) + ABS ( DEP ) + ( 0.2_DP * DEAVG ) ) / 2.2_DP
   END SELECT

   ! . Branch on the five cases to be distinguished.
   IF ( DE > 0.0_DP ) THEN

      ! . DE > 0.0 , DEP > 0.0.
      IF ( DEP > 0.0_DP ) THEN
         DAMP = 4.0_DP * MAX ( DAMP, DEAVG )
         IF ( DE >= ( 4.0_DP * DEP ) ) THEN
            DAMP = FAC * MAX ( DAMP, DEAVG )
         ELSE IF ( DE <= ( 0.25_DP * DEP ) ) THEN
            DAMP = DAMP / FAC
         ELSE
            DAMP = ( DE / DEP )**2 * MAX ( DAMP, DEAVG )
         END IF

      ! . DE > 0.0 , DEP < 0.0.
      ELSE
         DAMP = 4.0_DP * MAX ( DAMP, DEAVG )
         IF ( DE > 0.5_DP * DEAVG ) DAMP = DAMP * FAC
         IF ( ( DE - DEP ) < ( 0.2_DP * DEAVG ) ) THEN
            DAMP = DAMP / FAC
         END IF
      END IF

   ELSE

      ! . DE < 0.0 , DEP > 0.0.
      IF ( DEP > 0.0_DP ) THEN
         DAMP = 4.0_DP * MAX ( DAMP, DEAVG )
         IF ( -DE > DEAVG ) DAMP = DAMP * FAC
         IF ( ( -DE + DEP ) < DEAVG ) THEN
            DAMP = DAMP / FAC
         END IF

      ! . DE < 0.0 , DEP < 0.0.
      ELSE

         ! . DE > DEP.
         IF ( DE > DEP ) THEN
            IF ( DE <= ( 0.25_DP * DEP ) ) THEN
               DAMP = ( DE / DEP )**2 * MAX ( DAMP, DEAVG )
            ELSE
               DAMP = DAMP / FAC
            END IF

         ! . DE < DEP.
         ELSE
            IF ( ABS ( DE ) >= ( 2.0_DP * DEAVG ) ) THEN
               DAMP = FAC * MAX ( DAMP, DEAVG )
            ELSE
               IF ( ABS ( DE ) <= ( 0.5_DP * DEAVG ) ) DAMP = DAMP / FAC
            END IF
         END IF
      END IF
   END IF

   ! . Limit how big DAMP can get.
   DAMP = MIN ( DAMP, DMPMAX )

   END SUBROUTINE GET_DAMP_FACTOR

   !----------------------------------------
   SUBROUTINE SPIN_VALUES ( DENALP, DENBET )
   !----------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: DENALP, DENBET

   ! . Local scalars.
   REAL ( KIND = DP ) :: SZ, S2

   ! . Calculate Sz.
   SZ = 0.5_DP * REAL ( NALPHA - NBETA, DP )

   ! . Calculate S2.
   S2 = SZ**2 + 0.5_DP * REAL ( NALPHA + NBETA, DP ) - SYMMETRIC_MATRIX_DOT_PRODUCT ( DENALP, DENBET )

   ! . Write out the spin quantum numbers.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800" )
   CALL PRINT_SUMMARY_START ( "Wavefunction Spin Values" )
   WRITE ( PRINT_LINE, "(G14.8)" ) SZ ; CALL PRINT_SUMMARY_ELEMENT ( "<Sz>"  )
   WRITE ( PRINT_LINE, "(G14.8)" ) S2 ; CALL PRINT_SUMMARY_ELEMENT ( "<S^2>" )
   CALL PRINT_SUMMARY_STOP

   END SUBROUTINE SPIN_VALUES

!-------------------------------------------------------------------------------
! . Symmetric matrix procedures.
!-------------------------------------------------------------------------------

   !-----------------------------------------------------
   FUNCTION SYMMETRIC_MATRIX_COLUMN ( COLUMN, N, MATRIX )
   !-----------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: COLUMN, N

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: MATRIX

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:N) :: SYMMETRIC_MATRIX_COLUMN

   ! . Local scalars.
   INTEGER :: I

   ! . I <= COLUMN.
   SYMMETRIC_MATRIX_COLUMN(1:COLUMN) = MATRIX(BFINDEX(COLUMN)+1:BFINDEX(COLUMN)+COLUMN)

   ! . I > COLUMN.
   DO I = (COLUMN+1),N
      SYMMETRIC_MATRIX_COLUMN(I) = MATRIX(BFINDEX(I)+COLUMN)
   END DO

   END FUNCTION SYMMETRIC_MATRIX_COLUMN

   !-----------------------------------------------
   FUNCTION SYMMETRIC_MATRIX_DIAGONAL ( N, MATRIX )
   !-----------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: N

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: MATRIX

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:N) :: SYMMETRIC_MATRIX_DIAGONAL

   ! . Local scalars.
   INTEGER :: I

   ! . Loop over the elements in the matrix.
   DO I = 1,N
      SYMMETRIC_MATRIX_DIAGONAL(I) = MATRIX(BFINDEX(I)+I)
   END DO

   END FUNCTION SYMMETRIC_MATRIX_DIAGONAL

   !---------------------------------------------
   FUNCTION SYMMETRIC_MATRIX_DOT_PRODUCT ( A, B )
   !---------------------------------------------

   ! . Function declaration.
   REAL ( KIND = DP ) :: SYMMETRIC_MATRIX_DOT_PRODUCT

   ! . Array argument declarations.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: A, B

   ! . Calculate the dot product.
   SYMMETRIC_MATRIX_DOT_PRODUCT = 2.0_DP * DOT_PRODUCT ( A, B ) - &
                                           DOT_PRODUCT ( SYMMETRIC_MATRIX_DIAGONAL ( NBASIS, A ), &
                                                         SYMMETRIC_MATRIX_DIAGONAL ( NBASIS, B )  )

   END FUNCTION SYMMETRIC_MATRIX_DOT_PRODUCT

   !----------------------------------------------------------------------
   SUBROUTINE SYMMETRIC_MATRIX_MULTIPLY_2 ( A, B, C, FACTOR, QINITIALIZE )
   !----------------------------------------------------------------------

   ! . Calculate C = C + FACTOR * A * B.

   ! . Scalar argument declarations.
   LOGICAL,            INTENT(IN) :: QINITIALIZE
   REAL ( KIND = DP ), INTENT(IN) :: FACTOR

   ! . Array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:NBASTR), INTENT(IN)    :: A, B
   REAL ( KIND = DP ), DIMENSION(1:NBASTR), INTENT(INOUT) :: C

   ! . Local scalars.
   INTEGER :: I, IK, K

   ! . Initialize C if necessary.
   IF ( QINITIALIZE ) C = 0.0_DP

   ! . Check FACTOR.
   IF ( FACTOR == 0.0_DP ) RETURN

   ! . Initialize the matrix index.
   IK = 0

   ! . Loop over the rows of A.
   DO I = 1,NBASIS

      ! . Loop over the columns of B.
      DO K = 1,I
         IK = IK + 1
         C(IK) = C(IK) + FACTOR * DOT_PRODUCT ( SYMMETRIC_MATRIX_COLUMN ( I, NBASIS, A ), &
	                                        SYMMETRIC_MATRIX_COLUMN ( K, NBASIS, B )  )
      END DO
   END DO

   END SUBROUTINE SYMMETRIC_MATRIX_MULTIPLY_2

END MODULE MOPAC_SCF
