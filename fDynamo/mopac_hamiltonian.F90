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
!                           The MOPAC Hamiltonian Module
!===============================================================================
!
! . Subroutines:
!
!   MOPAC_SETUP        specify the MOPAC Hamiltonian and set up the calculation.
!
!===============================================================================
MODULE MOPAC_HAMILTONIAN

! . Utility data structure declarations.
USE CONSTANTS,        ONLY : AMU_TO_KG, EV_TO_KJ, KBOLTZ, KCAL_TO_KJ, NAVOGADRO, PI, PLANCK
USE DEFINITIONS,      ONLY : DP
USE ELEMENTS,         ONLY : NELEMENTS
USE PRINTING,         ONLY : PRINT_ERROR, PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, &
                             PRINT_SUMMARY_START, PRINT_SUMMARY_STOP

USE ATOMS,            ONLY : ATMMAS, ATMNUM, ATMQMI, NATOMS, NATOMSMM, NATOMSQM
USE GAUSSIAN_BASIS
USE MM_TERMS,         ONLY : ANGLES, ATMCHG, ATMCHG14, ATMEPS, ATMEPS14, ATME14I, ATME14J, BONDS, DIHEDRALS, IMPROPERS, &
                             NANGLES, NBONDS, NDIHEDRALS, NIMPROPERS
USE MOPAC_DATA
USE MOPAC_DENSITY,    ONLY : DENSITY_GUESS
USE MOPAC_PARAMETERS

IMPLICIT NONE
PRIVATE
PUBLIC :: MOPAC_SETUP

!===============================================================================
CONTAINS
!===============================================================================

   !--------------------------------------------------------------------------------------------------------------
   SUBROUTINE MOPAC_SETUP ( METHOD, CHARGE, MULTIPLICITY, SELECTION, LINK_ATOM_DISTANCES, PIATOMS, PITEMPERATURE )
   !--------------------------------------------------------------------------------------------------------------

   ! . Argument declarations.
   CHARACTER ( LEN = * ), INTENT(IN) :: METHOD

   ! . Optional scalar argument declarations.
   INTEGER, INTENT(IN),            OPTIONAL :: CHARGE, MULTIPLICITY
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: PITEMPERATURE

   ! . Optional array argument declarations.
   INTEGER,            DIMENSION(:,:),      INTENT(IN), OPTIONAL :: PIATOMS
   LOGICAL,            DIMENSION(1:NATOMS), INTENT(IN), OPTIONAL :: SELECTION
   REAL ( KIND = DP ), DIMENSION(:),        INTENT(IN), OPTIONAL :: LINK_ATOM_DISTANCES

   ! . Local scalars.
   INTEGER            :: I, IATOM, IBEAD, IBOUNDARY, IFIRST, ILAST, IPIATOM, IQM, J, K, KK, L, LL, &
                         NBOUNDARY, NHEAVY, NI, NLIGHT, NPIHEAVY, NPILIGHT, NUMORB, N0HEAVY, N0PIHEAVY
   LOGICAL            :: QPIATOM
   REAL ( KIND = DP ) :: EAT, ELECS, PITEMP

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: JIND4, KIND4
   INTEGER,              DIMENSION(1:NATOMS)  :: QATMNUM
   LOGICAL,              DIMENSION(1:NATOMS)  :: QMFLAG

   ! . Local link atom arrays.
   INTEGER,            ALLOCATABLE, DIMENSION(:) :: PARTNER
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: LENGTHS

   !----------------------------------------------------------------------------
   ! . Initialization.
   !----------------------------------------------------------------------------
   ! . Initialize the MOPAC data structure.
   CALL MOPAC_DATA_INITIALIZE

   ! . Set the Hamiltonian option.
   HAMILTONIAN = METHOD
   IF ( ( HAMILTONIAN /= "AM1 " ) .AND. ( HAMILTONIAN /= "MNDO" ) .AND. &
        ( HAMILTONIAN /= "PDDG" ) .AND. ( HAMILTONIAN /= "PM3 " ) .AND. &
        ( HAMILTONIAN /= "RM1 " ) ) THEN
      CALL PRINT_ERROR ( "MOPAC_SETUP", "Unknown MOPAC Hamiltonian." )
   END IF

   ! . Initialize the MOPAC parameter data with the Hamiltonian option.
   CALL MOPAC_PARAMETERS_INITIALIZE ( HAMILTONIAN )

   ! . Check the optional CHARGE and MULTIPLICITY arguments.
   IF ( PRESENT ( CHARGE       ) ) TOTCHG = CHARGE
   IF ( PRESENT ( MULTIPLICITY ) ) MULTIP = MULTIPLICITY

   ! . There is a QM atom selection array.
   IF ( PRESENT ( SELECTION ) ) THEN
      QMFLAG = SELECTION
   ! . All atoms are QM by default.
   ELSE
      QMFLAG = .TRUE.
   END IF

   ! . Set the NATOMSMM and NATOMSQM counters.
   NATOMSMM = COUNT ( .NOT. QMFLAG )
   NATOMSQM = COUNT (       QMFLAG )

   ! . Initialize ATMQMI.
   ATMQMI = 0

   ! . Check the number of QM atoms.
   IF ( NATOMSQM <= 0 ) RETURN

   !----------------------------------------------------------------------------
   ! . Get the atomic numbers of the quantum atoms.
   !----------------------------------------------------------------------------
   ! . Identify the types of the quantum atoms.
   CALL CHECK_QUANTUM_ATOMS

   ! . Increment the NATOMSQM counter by NBOUNDARY.
   ! . Note that NATOMSMM remains the same. 
   NATOMSQM = NATOMSQM + NBOUNDARY

   !----------------------------------------------------------------------------
   ! . Check for PI atoms.
   !----------------------------------------------------------------------------
   ! . The PI atoms argument is present.
   IF ( PRESENT ( PIATOMS ) ) THEN

      ! . Get the number of polymers and the number of atoms/polymer.
      NPIATOMS = SIZE ( PIATOMS, 1 )
      NPIBEADS = SIZE ( PIATOMS, 2 )

      ! . There are PI atoms.
      IF ( ( NPIBEADS > 1 ) .AND. ( NPIATOMS > 0 ) ) THEN

         ! . Allocate the PILIST array.
         ALLOCATE ( PILIST(1:NPIATOMS,1:NPIBEADS) )

         ! . Save PIATOMS.
         PILIST = PIATOMS

         ! . Check to see that all PI atoms are quantum.
         DO IBEAD = 1,NPIBEADS
            DO IATOM = 1,NPIATOMS
               IF ( .NOT. QMFLAG(PILIST(IATOM,IBEAD)) ) THEN
                  CALL PRINT_ERROR ( "MOPAC_HAMILTONIAN", "A PI atom is not quantum.", PILIST(IATOM,IBEAD) )
               END IF
            END DO
         END DO

         ! . Check the mass and number of corresponding atoms in each bead.
         DO IATOM = 1,NPIATOMS
            DO IBEAD = 2,NPIBEADS
               IF ( (  ATMMAS(PILIST(IATOM,1)) /=  ATMMAS(PILIST(IATOM,IBEAD)) ) .OR. &
                    ( QATMNUM(PILIST(IATOM,1)) /= QATMNUM(PILIST(IATOM,IBEAD)) ) ) THEN
                  CALL PRINT_ERROR ( "MOPAC_HAMILTONIAN", &
		                "Data mismatch between equivalent PI atoms in different polymers.", IATOM )
               END IF
            END DO
         END DO

         ! . Check to see if an atom occurs twice.
         DO IBEAD = 1,NPIBEADS
            DO IATOM = 1,NPIATOMS
               IF ( COUNT ( ( PILIST - PILIST(IATOM,IBEAD) ) == 0 ) /= 1 ) THEN
                  CALL PRINT_ERROR ( "MOPAC_HAMILTONIAN", "An atom occurs more than once in the definition of the PI atoms." )
               END IF
            END DO
         END DO

         ! . Check the PI temperature.
         IF ( PRESENT ( PITEMPERATURE ) ) THEN
            PITEMP = PITEMPERATURE
         ELSE
            PITEMP = 300.0_DP
         END IF

         ! . Calculate the force constant (kJ mol^-1 A^-2).
         PIFC = 0.5_DP * REAL ( NPIBEADS, DP ) * ( 1.0E-3_DP * AMU_TO_KG * NAVOGADRO ) * &
                ( ( 1.0E-10_DP * KBOLTZ * PITEMP ) / ( PLANCK / ( 2.0_DP * PI ) ) )**2 

      ! . There are no PI atoms.
      ELSE
         NPIATOMS = 0
         NPIBEADS = 1
      END IF

   END IF

   !----------------------------------------------------------------------------
   ! . Get the QM data array sizes and counters.
   !----------------------------------------------------------------------------
   ! . Initialize some accumulators.
   ATHEAT = 0.0_DP
   EAT    = 0.0_DP
   ELECS  = 0.0_DP

   ! . Initialise the atom counters.
   N0HEAVY   = 0
   N0PIHEAVY = 0
   NHEAVY    = 0
   NLIGHT    = 0
   NPIHEAVY  = 0
   NPILIGHT  = 0
   QPIATOM   = .FALSE.

   ! . Loop over the atoms.
   DO IATOM = 1,NATOMS

      ! . The atom is not QM.
      IF ( QATMNUM(IATOM) <= 0 ) CYCLE

      ! . There are PI atoms.
      IF ( NPIATOMS > 0 ) THEN

         ! . Skip PI atoms other than those in the first polymer.
         IF ( ANY ( PILIST(:,2:NPIBEADS) == IATOM ) ) CYCLE

         ! . Set the PI atom flag for atoms in the first polymer.
         QPIATOM = ANY ( PILIST(:,1) == IATOM )

      END IF

      ! . Check for an invalid atomic number.
      NI = QATMNUM(IATOM)
      IF ( ( NI <= 0 ) .OR. ( NI > NELEMENTS ) ) THEN
         CALL PRINT_ERROR ( "MOPAC_SETUP", "An atom has an invalid atomic number.", NI )
      END IF

      ! . Check on parameter availability.
      IF ( .NOT. SEPAR(NI) ) THEN
         CALL PRINT_ERROR ( "MOPAC_SETUP", "Parameters unavailable for an element.", NI )
      END IF

      ! . Accumulate some variables.
      ATHEAT = ATHEAT + EHEAT(NI)
      EAT    = EAT    + EISOL(NI)
      ELECS  = ELECS  +  CORE(NI)
      ! . Determine the orbital data.
      NUMORB = NATORB(NI)
!write ( 6, "(3i6,2f25.15)" ) iatom, ni, numorb, core(ni), elecs 

      ! . Check for D-elements.
      IF ( NUMORB > 4 ) THEN
         CALL PRINT_ERROR ( "MOPAC_SETUP", "Only S and SP elements allowed." )
      END IF

      ! . Increment the number of elements counters.
      IF ( NUMORB == 0 ) THEN
          N0HEAVY = N0HEAVY + 1
          IF ( QPIATOM ) N0PIHEAVY = N0PIHEAVY + 1
      ELSE IF ( NUMORB == 1 ) THEN
         NLIGHT = NLIGHT + 1
         IF ( QPIATOM ) NPILIGHT = NPILIGHT + 1
      ELSE
         NHEAVY = NHEAVY + 1
         IF ( QPIATOM ) NPIHEAVY = NPIHEAVY + 1
      END IF

   END DO

   ! . Find the total number of electrons.
   NELEC = - TOTCHG + NINT ( ELECS )

   ! . Find the number of alpha and beta electrons.
   NALPHA = ( NELEC + MULTIP - 1 ) / 2
   NBETA  = ( NELEC - MULTIP + 1 ) / 2
!write ( 6, "(a,5i6,f25.15)" ) "na, nb, ne = ", nalpha, nbeta, nelec, multip, totchg, elecs
   ! . Check for inconsistencies between the multiplicity and the number of electrons.
   IF ( ( NALPHA + NBETA ) /= NELEC ) CALL PRINT_ERROR ( "MOPAC_SETUP", "Impossible multiplicity for the system." )

   ! . Determine the heat of formation counter (in kJ/mole).
   ATHEAT = ( KCAL_TO_KJ * ATHEAT ) - ( EV_TO_KJ * EAT )

   ! . Determine the number of orbitals and some orbital related quantities.
   NBASIS = ( 4 * NHEAVY ) + NLIGHT
   NBASTR = ( NBASIS * ( NBASIS + 1 ) ) / 2

   ! . Determine the number of two-electron integrals.
   N2ELEC = (50*NHEAVY*(NHEAVY-1)) + (10*NHEAVY*NLIGHT) + (NLIGHT*(NLIGHT-1))/2

   ! . Set the PI basis function quantities.
   NPIBASIS = ( 4 * NPIHEAVY ) + NPILIGHT
   NPIHELE  = NBASTR - ( ( NBASIS - NPIBASIS ) * ( NBASIS - NPIBASIS + 1 ) ) / 2
   NPI2ELEC = (  50 * ( NPIHEAVY * ( NPIHEAVY - 1        ) ) ) + & ! PI heavy/PI heavy
              ( 100 *   NPIHEAVY * ( NHEAVY   - NPIHEAVY ) )   + & ! PI heavy/other heavy
              ( 10  *   NPIHEAVY *   NPILIGHT                ) + & ! PI heavy/PI light
              ( 10  *   NPIHEAVY * ( NLIGHT   - NPILIGHT ) )   + & ! PI heavy/other light
              ( 10  *   NPILIGHT * ( NHEAVY   - NPIHEAVY ) )   + & ! PI light/other heavy
              (         NPILIGHT * ( NLIGHT   - NPILIGHT ) )   + & ! PI light/other light
              (       ( NPILIGHT * ( NPILIGHT - 1 ) ) / 2  )       ! PI light/PI light

   ! . Reset N2ELEC.
   N2ELEC = N2ELEC - NPI2ELEC

   !----------------------------------------------------------------------------
   ! . Fill the QM data arrays.
   !----------------------------------------------------------------------------
   ! . Allocate QMATOM and the BETA and USPD arrays for the orbitals.
   ALLOCATE ( BETA(1:NBASIS), QMATOM(1:NATOMSQM), USPD(1:NBASIS) )

   ! . Initialize the basis function index counters.
   IBOUNDARY = 0
   IFIRST    = 1
   ILAST     = 0
   IQM       = 0

   ! . Fill the orbital data arrays for non-PI atoms.
   DO IATOM = 1,NATOMS

      ! . The atom is not QM.
      IF ( QATMNUM(IATOM) <= 0 ) CYCLE

      ! . Skip PI atoms for the moment.
      IF ( NPIATOMS > 0 ) THEN
         IF ( ANY ( PILIST == IATOM ) ) CYCLE
      END IF

      ! . Get the atomic number and the number of orbitals.
      NI     = QATMNUM(IATOM)
      NUMORB = NATORB(NI)

      ! . Determine the orbital indexing arrays.
      ILAST = IFIRST + NUMORB - 1

      ! . Increment the QM atom counter.
      IQM = IQM + 1

      ! . Set ATMQMI.
      ATMQMI(IATOM) = IQM

      ! . Fill QMATOM.
      QMATOM(IQM)%ATOM   = IATOM
      QMATOM(IQM)%BFIRST = IFIRST
      QMATOM(IQM)%BLAST  = ILAST
      QMATOM(IQM)%COPY   = 0
      QMATOM(IQM)%NUMBER = NI

      ! . The atom is not a boundary atom.
      IF ( QMFLAG(IATOM) ) THEN
	 QMATOM(IQM)%LENGTH    = 0.0_DP
	 QMATOM(IQM)%PARTNER   = 0
	 QMATOM(IQM)%QBOUNDARY = .FALSE.
      ! . The atom is a boundary atom.
      ELSE
         IBOUNDARY = IBOUNDARY + 1
         QMATOM(IQM)%LENGTH    = LENGTHS(IBOUNDARY)
	 QMATOM(IQM)%PARTNER   = PARTNER(IBOUNDARY)
	 QMATOM(IQM)%QBOUNDARY = .TRUE.
      END IF

      IF ( NUMORB > 0 ) THEN
      ! . Do the s orbital terms.
      BETA(IFIRST) = BETAS(NI)
      USPD(IFIRST) = USS(NI)

      ! . Do the p orbital terms.
      IF ( IFIRST < ILAST ) THEN
         BETA((IFIRST+1):ILAST) = BETAP(NI)
         USPD((IFIRST+1):ILAST) = UPP(NI)
      END IF

      ! . Update IFIRST.
      IFIRST = ILAST + 1
      END IF

   END DO

   ! . Deallocate the boundary atom arrays.
   IF ( NBOUNDARY > 0 ) DEALLOCATE ( LENGTHS, PARTNER )

   ! . Fill the orbital data arrays for PI atoms.
   DO IPIATOM = 1,NPIATOMS

      ! . Get the atom number of the atom in the first polymer.
      IATOM = PILIST(IPIATOM,1)

      ! . Get the atomic number and the number of orbitals.
      NI     = QATMNUM(IATOM)
      NUMORB = NATORB(NI)

      ! . Determine the orbital indexing arrays.
      ILAST = IFIRST + NUMORB - 1

      ! . Loop over all the beads.
      DO IBEAD = 1,NPIBEADS

         ! . Increment the QM atom counter.
         IQM = IQM + 1

         ! . Set ATMQMI.
         ATMQMI(PILIST(IPIATOM,IBEAD)) = IQM

         ! . Fill QMATOM.
         QMATOM(IQM)%ATOM   = PILIST(IPIATOM,IBEAD)
         QMATOM(IQM)%BFIRST = IFIRST
         QMATOM(IQM)%BLAST  = ILAST
         QMATOM(IQM)%COPY   = IBEAD
         QMATOM(IQM)%NUMBER = NI

         ! . Initialize all the link atom data for PI atoms.
	 QMATOM(IQM)%LENGTH    = 0.0_DP
	 QMATOM(IQM)%PARTNER   = 0
	 QMATOM(IQM)%QBOUNDARY = .FALSE.

      END DO

      ! . Do the s orbital terms.
      BETA(IFIRST) = BETAS(NI)
      USPD(IFIRST) = USS(NI)

      ! . Do the p orbital terms.
      IF ( IFIRST < ILAST ) THEN
         BETA((IFIRST+1):ILAST) = BETAP(NI)
         USPD((IFIRST+1):ILAST) = UPP(NI)
      END IF

      ! . Update IFIRST.
      IFIRST = ILAST + 1

   END DO

   ! . Set up the basis set.
   CALL GAUSSIAN_SETUP ( NHEAVY, NLIGHT, NPIHEAVY, NPILIGHT )

   ! . Set up the triangle indexing array.
   ALLOCATE ( BFINDEX(1:NBASIS) )
   DO I = 1,NBASIS
      BFINDEX(I) = ( I * ( I - 1 ) ) / 2
   END DO

   ! . Allocate space for the JIND4 and KIND4 arrays.
   ALLOCATE ( JIND4(1:4,1:4,1:4,1:4,1:2), KIND4(1:4,1:4,1:4,1:4,1:2) )

   ! . Set up the gather-scatter arrays used by the Fock matrix builder to index the semiempirical integrals.
   KK = 0
   LL = 0
   DO I = 1,4
      DO J = 1,I
         LL = LL + 1
         JIND4(1,1,I,J,2) = LL
         JIND4(1,1,J,I,2) = LL
         JIND4(I,J,1,1,2) = LL
         JIND4(J,I,1,1,2) = LL
         KIND4(I,1,1,J,2) = LL
         KIND4(J,1,1,I,2) = LL
         KIND4(I,1,J,1,2) = LL
         KIND4(J,1,I,1,2) = LL
         DO K = 1,4
            DO L = 1,K
               KK = KK + 1
               KIND4(K,I,J,L,1) = KK
               KIND4(K,J,I,L,1) = KK
               KIND4(L,I,J,K,1) = KK
               KIND4(L,J,I,K,1) = KK
               JIND4(K,L,I,J,1) = KK
               JIND4(L,K,I,J,1) = KK
               JIND4(K,L,J,I,1) = KK
               JIND4(L,K,J,I,1) = KK
            END DO
         END DO
      END DO
   END DO

   ! . Save the JIND2, JIND3 and KIND2 arrays.
   JIND2 = RESHAPE ( JIND4, SHAPE ( JIND2 ) )
   JIND3 = RESHAPE ( JIND4, SHAPE ( JIND3 ) )
   KIND2 = RESHAPE ( KIND4, SHAPE ( KIND2 ) )

   ! . Deallocate the JIND4 and KIND4 arrays.
   DEALLOCATE ( JIND4, KIND4 )

   ! . Write out some information.
   CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800" )
   CALL PRINT_SUMMARY_START ( "MOPAC Semi-Empirical QM Data" )
   CALL PRINT_SUMMARY_ELEMENT ( "Hamiltonian", TEXT = HAMILTONIAN )
   WRITE ( PRINT_LINE, "(I14)"   ) NATOMSQM  ; CALL PRINT_SUMMARY_ELEMENT ( "Number of QM Atoms"     )
   WRITE ( PRINT_LINE, "(I14)"   ) TOTCHG    ; CALL PRINT_SUMMARY_ELEMENT ( "System Charge"          )
   WRITE ( PRINT_LINE, "(I14)"   ) MULTIP    ; CALL PRINT_SUMMARY_ELEMENT ( "System Multiplicity"    )
   WRITE ( PRINT_LINE, "(I14)"   ) NBASIS    ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Orbitals"     )
   WRITE ( PRINT_LINE, "(I14)"   ) NALPHA    ; CALL PRINT_SUMMARY_ELEMENT ( "Occupied Orbitals"      )
   WRITE ( PRINT_LINE, "(I14)"   ) NALPHA    ; CALL PRINT_SUMMARY_ELEMENT ( "Number Alpha Electrons" )
   WRITE ( PRINT_LINE, "(I14)"   ) NBETA     ; CALL PRINT_SUMMARY_ELEMENT ( "Number Beta Electrons"  )
   WRITE ( PRINT_LINE, "(I14)"   ) NSHELL    ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Shells"       )
   IF ( NSHELL > 0 ) THEN
       WRITE ( PRINT_LINE, "(I14)"   ) ( KSTART(NSHELL) + KNG(NSHELL) - 1 )
   ELSE
       WRITE ( PRINT_LINE, "(I14)"   ) 0
   END IF
   CALL PRINT_SUMMARY_ELEMENT ( "Number of Primitives" )
   WRITE ( PRINT_LINE, "(I14)"   ) N2ELEC    ; CALL PRINT_SUMMARY_ELEMENT ( "Number of 2e Integrals" )
   WRITE ( PRINT_LINE, "(F14.3)" ) ATHEAT    ; CALL PRINT_SUMMARY_ELEMENT ( "Energy Base Line"       )
   WRITE ( PRINT_LINE, "(I14)"   ) NBOUNDARY ; CALL PRINT_SUMMARY_ELEMENT ( "No. of Boundary Atoms"  )
   WRITE ( PRINT_LINE, "(I14)"   ) NBASTR    ; CALL PRINT_SUMMARY_ELEMENT ( "Matrix Size"            )
   CALL PRINT_SUMMARY_STOP

   ! . Write out some information for the PI atoms.
   IF ( NPIATOMS > 0 ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800" )
      CALL PRINT_SUMMARY_START ( "Path Integral Atom Data" )
      WRITE ( PRINT_LINE, "(I14)"   ) NPIATOMS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of PI Atoms"     )
      WRITE ( PRINT_LINE, "(I14)"   ) NPIBEADS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of PI Polymers"  )
      WRITE ( PRINT_LINE, "(F14.4)" ) PITEMP   ; CALL PRINT_SUMMARY_ELEMENT ( "Temperature"            )
      WRITE ( PRINT_LINE, "(F14.4)" ) PIFC     ; CALL PRINT_SUMMARY_ELEMENT ( "Force Constant"         )
      WRITE ( PRINT_LINE, "(I14)"   ) NPIBASIS ; CALL PRINT_SUMMARY_ELEMENT ( "Number of PI Orbitals"  )
      WRITE ( PRINT_LINE, "(I14)"   ) NPI2ELEC ; CALL PRINT_SUMMARY_ELEMENT ( "Number of 2e Integrals" )
      CALL PRINT_SUMMARY_STOP
   END IF

   ! . Modify the MM terms data structure.
   CALL REMOVE_MM_TERMS

   ! . Set up the density and molecular orbital arrays.
   CALL DENSITY_GUESS

   !============================================================================
   CONTAINS
   !============================================================================

      !-----------------------------
      SUBROUTINE CHECK_QUANTUM_ATOMS
      !-----------------------------

      ! . Local scalars.
      INTEGER :: I, IBOND, IBOUNDARY,  J

      ! . Initialization.
      NBOUNDARY = 0
      QATMNUM   = 0

      ! . Loop over the atoms.
      DO I = 1,NATOMS

         ! . The atom is a quantum atom.
	 IF ( QMFLAG(I) ) QATMNUM(I) = ATMNUM(I)

      END DO

      ! . Check the number of bonds.
      IF ( NBONDS > 0 ) THEN

	 ! . Loop over the bonds.
	 DO IBOND = 1,NBONDS

            ! . Get the atoms in the bond.
            I = BONDS(IBOND)%I
            J = BONDS(IBOND)%J

            ! . There are bonds between QM and MM atoms so make the MM atom quantum.
            IF ( QMFLAG(I) .AND. ( .NOT. QMFLAG(J) ) ) THEN
               NBOUNDARY  = NBOUNDARY + 1
	       QATMNUM(J) = -I
	    ELSE IF ( QMFLAG(J) .AND. ( .NOT. QMFLAG(I) ) ) THEN
               NBOUNDARY  = NBOUNDARY + 1
	       QATMNUM(I) = -J
            END IF

	 END DO

      END IF

      ! . Process any boundary atom data.
      IF ( NBOUNDARY > 0 ) THEN

         ! . Allocate the boundary atom arrays.
	 ALLOCATE ( LENGTHS(1:NBOUNDARY), PARTNER(1:NBOUNDARY) )

         ! . Fill the partners array.
	 IBOUNDARY = 0
	 DO I = 1,NATOMS

            ! . The atom is a boundary atom.
	    IF ( QATMNUM(I) < 0 ) THEN

               ! . Fill PARTNER and reset QATMNUM.
	       IBOUNDARY = IBOUNDARY + 1
	       PARTNER(IBOUNDARY) = ABS ( QATMNUM(I) )
	       QATMNUM(I)         = 1

	    END IF
	 END DO

         ! . Check for the presence of the LINK_ATOM_DISTANCES argument.
	 IF ( PRESENT ( LINK_ATOM_DISTANCES ) ) THEN

            ! . The array is the correct size.
	    IF ( SIZE ( LINK_ATOM_DISTANCES ) == NBOUNDARY ) THEN
               LENGTHS = LINK_ATOM_DISTANCES
	    ! . The array is of the wrong size.
	    ELSE
	       CALL PRINT_ERROR ( "MOPAC_SETUP", "LINK_ATOM_DISTANCES argument has the wrong size.", NBOUNDARY )
	    END IF

         ! . Use default link atom - QM atom distances.
	 ELSE
	    LENGTHS = 1.0_DP
	 END IF

      ! . There are no boundary atoms.
      ELSE

         ! . Check for the presence of the LINK_ATOM_DISTANCES argument.
	 IF ( PRESENT ( LINK_ATOM_DISTANCES ) ) THEN
	    CALL PRINT_ERROR ( "MOPAC_SETUP", "LINK_ATOM_DISTANCES argument present but there are no boundary atoms." )
	 END IF

      END IF

      END SUBROUTINE CHECK_QUANTUM_ATOMS

      !-----------------------------
      FUNCTION QSCALE ( I, J, K, L )
      !-----------------------------

      ! . Function declarations.
      LOGICAL :: QSCALE

      ! . Scalar arguments.
      INTEGER, INTENT(IN)           :: I, J
      INTEGER, INTENT(IN), OPTIONAL :: K, L

      ! . Local arrays.
      INTEGER, DIMENSION(1:4) :: COPY

      ! . Initialization.
      COPY = 0

      ! . Get the copies to which the path integral atoms belong.
      IF ( QMFLAG(I) ) COPY(1) = QMATOM(ATMQMI(I))%COPY
      IF ( QMFLAG(J) ) COPY(2) = QMATOM(ATMQMI(J))%COPY
      IF ( PRESENT ( K ) ) THEN
         IF ( QMFLAG(K) ) COPY(3) = QMATOM(ATMQMI(K))%COPY
      END IF
      IF ( PRESENT ( L ) ) THEN
         IF ( QMFLAG(L) ) COPY(4) = QMATOM(ATMQMI(L))%COPY
      END IF

      ! . Set QSCALE.
      QSCALE = ANY ( COPY > 0 )

      END FUNCTION QSCALE

      !----------------------------
      FUNCTION QSKIP ( I, J, K, L )
      !----------------------------

      ! . Function declarations.
      LOGICAL :: QSKIP

      ! . Scalar arguments.
      INTEGER, INTENT(IN)           :: I, J
      INTEGER, INTENT(IN), OPTIONAL :: K, L

      ! . Local scalars.
      INTEGER :: MAXCOPY

      ! . Local arrays.
      INTEGER, DIMENSION(1:4) :: COPY

      ! . Check to see if all the atoms are QM.
      QSKIP = ( QATMNUM(I) > 0 ) .AND. ( QATMNUM(J) > 0 )
      IF ( PRESENT ( K ) ) QSKIP = QSKIP .AND. ( QATMNUM(K) > 0 )
      IF ( PRESENT ( L ) ) QSKIP = QSKIP .AND. ( QATMNUM(L) > 0 )

      ! . Return at this stage if all atoms are QM (including boundary atoms).
      IF ( QSKIP ) RETURN

      ! . Initialization.
      COPY = 0

      ! . Get the copies to which the path integral atoms belong.
      IF ( QMFLAG(I) ) COPY(1) = QMATOM(ATMQMI(I))%COPY
      IF ( QMFLAG(J) ) COPY(2) = QMATOM(ATMQMI(J))%COPY
      IF ( PRESENT ( K ) ) THEN
         IF ( QMFLAG(K) ) COPY(3) = QMATOM(ATMQMI(K))%COPY
      END IF
      IF ( PRESENT ( L ) ) THEN
         IF ( QMFLAG(L) ) COPY(4) = QMATOM(ATMQMI(L))%COPY
      END IF

      ! . Get the maximum value of COPY.
      MAXCOPY = MAXVAL ( COPY )

      ! . Adjust COPY.
      WHERE ( COPY == 0 ) COPY = COPY + MAXCOPY

      ! . Set QSKIP.
      QSKIP = ( MINVAL ( COPY ) /= MAXCOPY )

      END FUNCTION QSKIP

      !-------------------------
      SUBROUTINE REMOVE_MM_TERMS
      !-------------------------

      ! . Local scalars.
      INTEGER :: I, IINT, J, N, NANG, NBND, NDIH, NE14, NIMP

      ! . Local arrays.
      INTEGER, ALLOCATABLE, DIMENSION(:)        :: EX14I, EX14J

      ! . Check the number of MM terms.
      IF ( ( NANGLES <= 0 ) .AND. ( NBONDS <= 0 ) .AND. ( NDIHEDRALS <= 0 ) .AND. ( NIMPROPERS <= 0 ) ) RETURN

      ! . Loop over the angles.
      N = 0
      DO I = 1,NANGLES

         ! . Skip angles between QM atoms or those involving PI atoms from different copies.
         IF ( QSKIP ( ANGLES(I)%I, ANGLES(I)%J, ANGLES(I)%K ) ) CYCLE

         ! . Increment the angle number.
         N = N + 1

         ! . Keep the angle.
         ANGLES(N) = ANGLES(I)

         ! . Check to see if the MM term is to be scaled.
	 IF ( QSCALE ( ANGLES(I)%I, ANGLES(I)%J, ANGLES(I)%K ) ) ANGLES(N)%FC = ANGLES(N)%FC / REAL ( NPIBEADS, DP )

      END DO
      NANG    = NANGLES - N
      NANGLES = N

      ! . Loop over the bonds.
      N = 0
      DO I = 1,NBONDS

         ! . Skip bonds between QM atoms or those involving PI atoms from different copies.
         IF ( QSKIP ( BONDS(I)%I, BONDS(I)%J ) ) THEN

            ! . Only skip bonds between quantum atoms but not those between boundary and quantum atoms.
	    IF ( QMFLAG(BONDS(I)%I) .AND. QMFLAG(BONDS(I)%J) ) CYCLE

	 END IF

         ! . Increment the bond number.
         N = N + 1

         ! . Keep the bond.
         BONDS(N) = BONDS(I)

         ! . Check to see if the MM term is to be scaled.
	 IF ( QSCALE ( BONDS(I)%I, BONDS(I)%J ) ) BONDS(N)%FC = BONDS(N)%FC / REAL ( NPIBEADS, DP )

      END DO
      NBND   = NBONDS - N
      NBONDS = N

      ! . Loop over the dihedrals.
      N = 0
      DO I = 1,NDIHEDRALS

         ! . Skip dihedrals between QM atoms or those involving PI atoms from different copies.
         IF ( QSKIP ( DIHEDRALS(I)%I, DIHEDRALS(I)%J, DIHEDRALS(I)%K, DIHEDRALS(I)%L ) ) CYCLE

         ! . Increment the dihedral number.
         N = N + 1

         ! . Keep the dihedral.
         DIHEDRALS(N) = DIHEDRALS(I)

         ! . Check to see if the MM term is to be scaled.
	 IF ( QSCALE ( DIHEDRALS(I)%I, DIHEDRALS(I)%J, DIHEDRALS(I)%K, DIHEDRALS(I)%L ) ) &
	                                               DIHEDRALS(N)%FC = DIHEDRALS(N)%FC / REAL ( NPIBEADS, DP )

      END DO
      NDIH       = NDIHEDRALS - N
      NDIHEDRALS = N

      ! . Loop over the impropers.
      N = 0
      DO I = 1,NIMPROPERS

         ! . Skip impropers between QM atoms or those involving PI atoms from different copies.
         IF ( QSKIP ( IMPROPERS(I)%I, IMPROPERS(I)%J, IMPROPERS(I)%K, IMPROPERS(I)%L ) ) CYCLE

         ! . Increment the improper number.
         N = N + 1

         ! . Keep the improper.
         IMPROPERS(N) = IMPROPERS(I)

         ! . Check to see if the MM term is to be scaled.
	 IF ( QSCALE ( IMPROPERS(I)%I, IMPROPERS(I)%J, IMPROPERS(I)%K, IMPROPERS(I)%L ) ) &
	                                               IMPROPERS(N)%FC = IMPROPERS(N)%FC / REAL ( NPIBEADS, DP )

      END DO
      NIMP       = NIMPROPERS - N
      NIMPROPERS = N

      ! . Remove exclusions between QM atoms.
      N = 0

      ! . Initialize the temporary exclusion arrays.
      ALLOCATE ( EX14I(1:NATOMS+1), EX14J(1:SIZE(ATME14J)) )

      ! . Loop over atoms.
      DO I = 1,NATOMS

         ! . Set the new index.
         EX14I(I) = N

         ! . Loop over the exclusions.
         DO IINT = (ATME14I(I)+1),ATME14I(I+1)

            ! . Get the second atom of the exclusion.
            J = ATME14J(IINT)

            ! . Skip exclusions between QM atoms but keep them if either is a boundary atom. This is needed
	    ! . because the 1-4 lists are used to calculate the MM/MM and QM/MM 1-4 interactions.
            IF ( QMFLAG(I) .AND. QMFLAG(J) ) CYCLE

            ! . Copy the exclusion.
            N = N + 1
            EX14J(N) = J

         END DO
      END DO

      ! . Set the last index.
      EX14I(NATOMS+1) = N

      ! . Save the number of exclusions removed.
      NE14 = SIZE ( ATME14J ) - N

      ! . Reallocate the exclusion arrays.
      DEALLOCATE ( ATME14I, ATME14J )
      ALLOCATE ( ATME14I(1:NATOMS+1), ATME14J(1:N) )

      ! . Save the exclusions.
      ATME14I(1:NATOMS+1) = EX14I(1:NATOMS+1)
      ATME14J(1:N)        = EX14J(1:N)

      ! . Deallocate the temporary arrays.
      DEALLOCATE ( EX14I, EX14J )

      ! . Zero the charges on the QM atoms except for the boundary atoms.
      WHERE ( QMFLAG ) ATMCHG   = 0.0_DP
      WHERE ( QMFLAG ) ATMCHG14 = 0.0_DP

      ! . Scale the LJ parameters on PI atoms by 1/NPIBEADS.
      IF ( NPIATOMS > 0 ) THEN
         DO I = 1,NPIBEADS
            DO J = 1,NPIATOMS
               ATMEPS(PILIST(J,I))   = ATMEPS(PILIST(J,I))   / REAL ( NPIBEADS, DP )
               ATMEPS14(PILIST(J,I)) = ATMEPS14(PILIST(J,I)) / REAL ( NPIBEADS, DP )
            END DO
         END DO
      END IF

      ! . Write out a summary of the MM terms.
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#FF8800" )
      CALL PRINT_SUMMARY_START ( "MM Terms Removed" )
      WRITE ( PRINT_LINE, "(I14)" ) NBND ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Bonds"     )
      WRITE ( PRINT_LINE, "(I14)" ) NANG ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Angles"    )
      WRITE ( PRINT_LINE, "(I14)" ) NDIH ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Dihedrals" )
      WRITE ( PRINT_LINE, "(I14)" ) NIMP ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Impropers" )
      WRITE ( PRINT_LINE, "(I14)" ) NE14 ; CALL PRINT_SUMMARY_ELEMENT ( "Number of 1-4 Excl." )
      WRITE ( PRINT_LINE, "(F14.4)" ) SUM( ATMCHG(1:NATOMS) ); CALL PRINT_SUMMARY_ELEMENT ( "Total MM Charge " )
      CALL PRINT_SUMMARY_STOP

      END SUBROUTINE REMOVE_MM_TERMS

   END SUBROUTINE MOPAC_SETUP

END MODULE MOPAC_HAMILTONIAN

 
