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
!                           The MOPAC Gradients Module
!===============================================================================
!
! . Subroutines.
!
!   GRADIENTS                  Calculate the MOPAC gradients.
!
!===============================================================================
MODULE MOPAC_GRADIENTS

! . Module declarations.
USE CONSTANTS,          ONLY : ANGSTROMS_TO_BOHRS, AU_TO_EV, EV_TO_KJ
USE DEFINITIONS,        ONLY : DP

USE ATOMS,              ONLY : ATMCRD, ATMFIX, ATMQMI, NATOMSMM, NATOMSQM
USE ENERGY_NON_BONDING, ONLY : NBLIST_TYPE, NBLISTQM_FIRST, NNBLISTQM
USE GAUSSIAN_BASIS,     ONLY : GAUSSIAN_OVERLAP_DERIVATIVES, KATOM, KLOC, NPISHELL, NSHELL, PIKATOM
USE MOPAC_DATA,         ONLY : BETA, DENMAT, DENMATA, DENMATB, HAMILTONIAN, MOPAC_DATA_ATOM_POSITION, &
                               MOPAC_DATA_LA_GRADIENT, NBASIS, NBASTR, NPIBASIS, NPIBEADS, NPIHELE, QMATOM
USE MOPAC_INTEGRALS,    ONLY : DXE1BA, DXE1BB, DYE1BA, DYE1BB, DZE1BA, DZE1BB, DXQM, DYQM, DZQM, REPP
USE MOPAC_PARAMETERS,   ONLY : AD, ALP, AM, AQ, CORE, DD, FN1, FN2, FN3, NATORB, PDDGE, PDDGC, QQ

IMPLICIT NONE
PRIVATE
PUBLIC :: GRADIENTS

!===============================================================================
CONTAINS
!===============================================================================

   !----------------------------------------
   SUBROUTINE GRADIENTS ( VIRIAL, GRADIENT )
   !----------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT) :: GRADIENT

   ! . Local scalars.
   LOGICAL :: QAM1PM3

   ! . Set the AM1/PM3 flag.
   QAM1PM3 = ( HAMILTONIAN == "AM1 " ) .OR. ( HAMILTONIAN == "PDDG" ) .OR. ( HAMILTONIAN == "PM3 " ) .OR. ( HAMILTONIAN == "RM1 " )

   ! . Calculate the purely QM derivatives.
   CALL DERIVATIVES_QM

   ! . Calculate the QM/MM derivatives.
   IF ( NATOMSMM > 0 ) CALL DERIVATIVES_QM_MM

   !============================================================================
   CONTAINS
   !============================================================================

      !------------------------
      SUBROUTINE DERIVATIVES_QM
      !------------------------

      ! . Local scalars.
      INTEGER            :: I, IATOM, IBASIS, ICOPY, IFIRST, IJ, ILAST, IQM, ISHELL, ISTART,    &
                            ISTOP, J, JATOM, JBASIS, JCOPY, JFIRST, JLAST, JQM, JSHELL, JSTART, &
                            JSTOP, K, L, NNOTRI, P
      LOGICAL            :: QUHF
      REAL ( KIND = DP ) :: BETAIJ, PFAC

      ! . Local arrays.
      INTEGER, DIMENSION(1:2)                         :: NDI
      REAL ( KIND = DP ), DIMENSION(1:3,1:2)          :: CDI
      REAL ( KIND = DP ), DIMENSION(1:3)              :: ENG, MODGVEC
      REAL ( KIND = DP ), DIMENSION(1:171)            :: PDIT
      REAL ( KIND = DP ), DIMENSION(1:NBASTR-NPIHELE) :: DENTMP
      REAL ( KIND = DP ), DIMENSION(1:171,1:NPIBEADS) :: PDALP, PDBET, PDIA
      REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMSQM)   :: GRADQM

      ! . Local allocatable arrays.
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   ::   SX,   SY,   SZ
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: PISX, PISY, PISZ

      ! . Initialize the QM gradient.
      GRADQM = 0.0_DP

      !-------------------------------------------------------------------------
      ! . Calculate the derivatives of the overlap terms.
      !-------------------------------------------------------------------------
      ! . Initialize some counters.
      NNOTRI = NBASTR - NPIHELE

      ! . Allocate the overlap derivative matrix arrays.
      ALLOCATE ( PISX(1:NPIHELE,1:NPIBEADS), PISY(1:NPIHELE,1:NPIBEADS), PISZ(1:NPIHELE,1:NPIBEADS), &
                 SX(1:NBASTR-NPIHELE), SY(1:NBASTR-NPIHELE), SZ(1:NBASTR-NPIHELE) )

      ! . Calculate the derivatives of the overlap matrix.
      CALL GAUSSIAN_OVERLAP_DERIVATIVES ( SX, SY, SZ, PISX, PISY, PISZ )

      ! . Scale the derivatives by the BETA factors.
      I = 0
      DO IBASIS = 1,(NBASIS-NPIBASIS)
         DO JBASIS = 1,IBASIS
            BETAIJ = 0.5_DP * ( BETA(IBASIS) + BETA(JBASIS) )
            I = I + 1
            SX(I) = BETAIJ * SX(I)
            SY(I) = BETAIJ * SY(I)
            SZ(I) = BETAIJ * SZ(I)
         END DO
      END DO
      I = I - NNOTRI
      DO IBASIS = (NBASIS-NPIBASIS+1),NBASIS
         DO JBASIS = 1,IBASIS
            BETAIJ = 0.5_DP * ( BETA(IBASIS) + BETA(JBASIS) )
            I = I + 1
            PISX(I,1:NPIBEADS) = BETAIJ * PISX(I,1:NPIBEADS)
            PISY(I,1:NPIBEADS) = BETAIJ * PISY(I,1:NPIBEADS)
            PISZ(I,1:NPIBEADS) = BETAIJ * PISZ(I,1:NPIBEADS)
         END DO
      END DO

      ! . Scale the matrices by the appropriate factors.
      PISX = 2.0_DP * ANGSTROMS_TO_BOHRS * EV_TO_KJ * PISX / REAL ( NPIBEADS, DP )
      PISY = 2.0_DP * ANGSTROMS_TO_BOHRS * EV_TO_KJ * PISY / REAL ( NPIBEADS, DP )
      PISZ = 2.0_DP * ANGSTROMS_TO_BOHRS * EV_TO_KJ * PISZ / REAL ( NPIBEADS, DP )
        SX = 2.0_DP * ANGSTROMS_TO_BOHRS * EV_TO_KJ *   SX / REAL ( NPIBEADS, DP )
        SY = 2.0_DP * ANGSTROMS_TO_BOHRS * EV_TO_KJ *   SY / REAL ( NPIBEADS, DP )
        SZ = 2.0_DP * ANGSTROMS_TO_BOHRS * EV_TO_KJ *   SZ / REAL ( NPIBEADS, DP )

      ! . There are PI atoms.
      IF ( NPIBEADS > 1 ) THEN

         ! . Save the non-PI atom elements of the first density matrix.
         DENTMP(1:NBASTR-NPIHELE) = DENMAT(1:NBASTR-NPIHELE,1)

         ! . Fill the elements of the density matrix with the sum of all the matrices.
         DO I = 1,NBASTR-NPIHELE
            DENMAT(I,1) = SUM ( DENMAT(I,1:NPIBEADS) )
         END DO
         
      END IF

      ! . Add the overlap terms into the derivatives.
      I = 0
      DO ISHELL = 1,(NSHELL-NPISHELL)
         IQM    = KATOM(ISHELL)
         ISTART = KLOC(ISHELL)
         IF ( ISHELL /= NSHELL ) THEN
            ISTOP  = KLOC(ISHELL+1) - 1
         ELSE
            ISTOP = NBASIS
         END IF
         DO IBASIS = ISTART,ISTOP
            DO JSHELL = 1,ISHELL
               JQM    = KATOM(JSHELL)
               JSTART = KLOC(JSHELL)
               IF ( JSHELL /= NSHELL ) THEN
                  JSTOP  = MIN ( IBASIS, (KLOC(JSHELL+1)-1) )
               ELSE
                  JSTOP  = MIN ( IBASIS, NBASIS )
               END IF

               ! . Calculate the derivatives.
               DO JBASIS = JSTART,JSTOP
                  I = I + 1
                  PFAC = DENMAT(I,1)
                  GRADQM(1,IQM) = GRADQM(1,IQM) + PFAC * SX(I)
                  GRADQM(2,IQM) = GRADQM(2,IQM) + PFAC * SY(I)
                  GRADQM(3,IQM) = GRADQM(3,IQM) + PFAC * SZ(I)
                  GRADQM(1,JQM) = GRADQM(1,JQM) - PFAC * SX(I)
                  GRADQM(2,JQM) = GRADQM(2,JQM) - PFAC * SY(I)
                  GRADQM(3,JQM) = GRADQM(3,JQM) - PFAC * SZ(I)
               END DO
            END DO
         END DO
      END DO

      ! . There are PI atoms.
      IF ( NPIBEADS > 1 ) THEN

         ! . Restore the density matrix.
         DENMAT(1:NBASTR-NPIHELE,1) = DENTMP(1:NBASTR-NPIHELE)

         ! . Loop over the polymers.
         DO P = 1,NPIBEADS

            ! . Make sure KATOM is up to date.
            KATOM(NSHELL-NPISHELL+1:NSHELL) = PIKATOM(1:NPISHELL,P)

            ! . Loop over the PI atom shells.
            I = 0
            DO ISHELL = (NSHELL-NPISHELL+1),NSHELL
               IQM    = KATOM(ISHELL)
               ISTART = KLOC(ISHELL)
               IF ( ISHELL /= NSHELL ) THEN
                  ISTOP = KLOC(ISHELL+1) - 1
               ELSE
                  ISTOP = NBASIS
               END IF
               DO IBASIS = ISTART,ISTOP
                  DO JSHELL = 1,ISHELL
                     JQM    = KATOM(JSHELL)
                     JSTART = KLOC(JSHELL)
                     IF ( JSHELL /= NSHELL ) THEN
                        JSTOP  = MIN ( IBASIS, (KLOC(JSHELL+1)-1) )
                     ELSE
                        JSTOP  = MIN ( IBASIS, NBASIS )
                     END IF

                     ! . Calculate the derivatives.
                     DO JBASIS = JSTART,JSTOP
                        I = I + 1
                        PFAC = DENMAT(I+NNOTRI,P)
                        GRADQM(1,IQM) = GRADQM(1,IQM) + PFAC * PISX(I,P)
                        GRADQM(2,IQM) = GRADQM(2,IQM) + PFAC * PISY(I,P)
                        GRADQM(3,IQM) = GRADQM(3,IQM) + PFAC * PISZ(I,P)
                        GRADQM(1,JQM) = GRADQM(1,JQM) - PFAC * PISX(I,P)
                        GRADQM(2,JQM) = GRADQM(2,JQM) - PFAC * PISY(I,P)
                        GRADQM(3,JQM) = GRADQM(3,JQM) - PFAC * PISZ(I,P)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF

      ! . Deallocate some temporary space.
      DEALLOCATE ( PISX, PISY, PISZ, SX, SY, SZ )

      !-------------------------------------------------------------------------
      ! . Calculate the derivatives of the two-electron terms.
      !-------------------------------------------------------------------------
      ! . Set the UHF flag.
      QUHF = ALLOCATED ( DENMATA ) .AND. ALLOCATED ( DENMATB )

      ! . Outer loop over atoms.
      DO IQM = 1,NATOMSQM

         ! . Get some information about the atom.
         IATOM  = QMATOM(IQM)%ATOM
         ICOPY  = QMATOM(IQM)%COPY
         IFIRST = QMATOM(IQM)%BFIRST
         ILAST  = QMATOM(IQM)%BLAST
         NDI(2) = QMATOM(IQM)%NUMBER
         CDI(1:3,2) = MOPAC_DATA_ATOM_POSITION ( IQM, ATMCRD )

         ! . Inner loop over atoms.
         DO JQM = 1,(IQM-1)

            ! . Get some information about the atom.
            JATOM  = QMATOM(JQM)%ATOM
            JCOPY  = QMATOM(JQM)%COPY
            JFIRST = QMATOM(JQM)%BFIRST
            JLAST  = QMATOM(JQM)%BLAST
            NDI(1) = QMATOM(JQM)%NUMBER
            CDI(1:3,1) = MOPAC_DATA_ATOM_POSITION ( JQM, ATMCRD )

            ! . Do nothing if both atoms are fixed.
            IF ( ATMFIX(IATOM) .AND. ATMFIX(JATOM) ) CYCLE

            ! . Two non-PI atoms but for the case where other PI atoms are present.
            IF ( ( ICOPY == 0 ) .AND. ( JCOPY == 0 ) .AND. ( NPIBEADS > 1 ) ) THEN

               ! . Form the closed shell diatomic matrices.
               IJ = 0
               DO I = JFIRST,JLAST
                  K = I*(I-1)/2+JFIRST-1
                  DO J = JFIRST,I
                     IJ = IJ+1
                     K  = K+1
                     PDIA(IJ,1:NPIBEADS) = DENMAT(K,1:NPIBEADS)
                     PDIT(IJ)            = SUM ( DENMAT(K,1:NPIBEADS) )
                  END DO
               END DO

               ! . Get the second atom - first atom intersection.
               DO I = IFIRST,ILAST
                  L = I*(I-1)/2
                  K = L+JFIRST-1
                  DO J = JFIRST,JLAST
                     IJ = IJ+1
                     K  = K+1
                     PDIA(IJ,1:NPIBEADS) = DENMAT(K,1:NPIBEADS)
                     PDIT(IJ)            = SUM ( DENMAT(K,1:NPIBEADS) )
                  END DO
                  K = L+IFIRST-1
                  DO L = IFIRST,I
                     IJ = IJ+1
                     K  = K+1
                     PDIA(IJ,1:NPIBEADS) = DENMAT(K,1:NPIBEADS)
                     PDIT(IJ)            = SUM ( DENMAT(K,1:NPIBEADS) )
                  END DO
               END DO

               ! . Form the open shell diatomic matrices.
	       ! . Spin-unrestricted calculation.
	       IF ( QUHF ) THEN

        	  IJ = 0
        	  DO I = JFIRST,JLAST
                     K = I*(I-1)/2+JFIRST-1
                     DO J = JFIRST,I
                	IJ = IJ+1
                	K  = K+1
                	PDALP(IJ,1:NPIBEADS) = DENMATA(K,1:NPIBEADS)
                	PDBET(IJ,1:NPIBEADS) = DENMATB(K,1:NPIBEADS)
                     END DO
        	  END DO

        	  ! . Get the second atom - first atom intersection.
        	  DO I = IFIRST,ILAST
                     L = I*(I-1)/2
                     K = L+JFIRST-1
                     DO J = JFIRST,JLAST
                	IJ = IJ+1
                	K  = K+1
                	PDALP(IJ,1:NPIBEADS) = DENMATA(K,1:NPIBEADS)
                	PDBET(IJ,1:NPIBEADS) = DENMATB(K,1:NPIBEADS)
                     END DO
                     K = L+IFIRST-1
                     DO L = IFIRST,I
                	IJ = IJ+1
                	K  = K+1
                	PDALP(IJ,1:NPIBEADS) = DENMATA(K,1:NPIBEADS)
                	PDBET(IJ,1:NPIBEADS) = DENMATB(K,1:NPIBEADS)
                     END DO
        	  END DO

               ! . Spin-restricted calculation.
               ELSE
	          PDALP(1:IJ,1:NPIBEADS) = 0.5_DP * PDIA(1:IJ,1:NPIBEADS)
		  PDBET(1:IJ,1:NPIBEADS) = 0.5_DP * PDIA(1:IJ,1:NPIBEADS)
	       END IF

               ! . Calculate the derivatives.
               CALL ANALYT ( PDIT, PDIA, PDALP, PDBET, NPIBEADS, CDI, NDI, JFIRST, JLAST, IFIRST, ILAST, ENG, QAM1PM3 )

            ! . No PI atoms or a PI atom with a non-PI atom or two non-PI atoms in same polymer.
            ELSE IF ( ( ( ICOPY == 0 ) .AND. ( JCOPY == 0 ) .AND. ( NPIBEADS == 1 ) ) .OR. &
                      ( ( ICOPY == 0 ) .OR.  ( JCOPY == 0 ) .OR.  ( ICOPY == JCOPY ) ) ) THEN

               ! . Get the density matrix index.
               P = MAX ( ICOPY, JCOPY, 1 )

               ! . Form the closed shell diatomic matrices.
               IJ = 0
               DO I = JFIRST,JLAST
                  K = I*(I-1)/2+JFIRST-1
                  DO J = JFIRST,I
                     IJ = IJ+1
                     K  = K+1
                     PDIA(IJ,1) = DENMAT(K,P)
                     PDIT(IJ)   = DENMAT(K,P)
                  END DO
               END DO

               ! . Get the second atom - first atom intersection.
               DO I = IFIRST,ILAST
                  L = I*(I-1)/2
                  K = L+JFIRST-1
                  DO J = JFIRST,JLAST
                     IJ = IJ+1
                     K  = K+1
                     PDIA(IJ,1) = DENMAT(K,P)
                     PDIT(IJ)   = DENMAT(K,P)
                  END DO
                  K = L+IFIRST-1
                  DO L = IFIRST,I
                     IJ = IJ+1
                     K  = K+1
                     PDIA(IJ,1) = DENMAT(K,P)
                     PDIT(IJ)	= DENMAT(K,P)
                  END DO
               END DO

               ! . Form the open shell diatomic matrices.
	       ! . Spin-unrestricted calculation.
	       IF ( QUHF ) THEN

        	  IJ = 0
        	  DO I = JFIRST,JLAST
                     K = I*(I-1)/2+JFIRST-1
                     DO J = JFIRST,I
                	IJ = IJ+1
                	K  = K+1
                	PDALP(IJ,1) = DENMATA(K,P)
                	PDBET(IJ,1) = DENMATB(K,P)
                     END DO
        	  END DO

        	  ! . Get the second atom - first atom intersection.
        	  DO I = IFIRST,ILAST
                     L = I*(I-1)/2
                     K = L+JFIRST-1
                     DO J = JFIRST,JLAST
                	IJ = IJ+1
                	K  = K+1
                	PDALP(IJ,1) = DENMATA(K,P)
                	PDBET(IJ,1) = DENMATB(K,P)
                     END DO
                     K = L+IFIRST-1
                     DO L = IFIRST,I
                	IJ = IJ+1
                	K  = K+1
                	PDALP(IJ,1) = DENMATA(K,P)
                	PDBET(IJ,1) = DENMATB(K,P)
                     END DO
        	  END DO

               ! . Spin-restricted calculation.
               ELSE
	          PDALP(1:IJ,1) = 0.5_DP * PDIA(1:IJ,1)
		  PDBET(1:IJ,1) = 0.5_DP * PDIA(1:IJ,1)
	       END IF

               ! . Calculate the derivatives.
               CALL ANALYT ( PDIT, PDIA(:,1:1), PDALP(:,1:1), PDBET(:,1:1), 1, CDI, NDI, JFIRST, JLAST, IFIRST, ILAST, &
	                                                                                                  ENG, QAM1PM3 )

            ! . No other atom combination is appropriate.
            ELSE
               CYCLE
            END IF

            ! . Add in the contribution to the derivatives.
            GRADQM(1:3,IQM) = GRADQM(1:3,IQM) + ENG(1:3) / REAL ( NPIBEADS, DP )
            GRADQM(1:3,JQM) = GRADQM(1:3,JQM) - ENG(1:3) / REAL ( NPIBEADS, DP )

         END DO
      END DO

      ! . Put the gradients into the proper array and calculate the virial.
      ! . Loop over the QM atoms.
      DO IQM = 1,NATOMSQM

         ! . Get the atom index of the atom.
	 IATOM = QMATOM(IQM)%ATOM

         ! . The atom is a boundary atom.
	 IF ( QMATOM(IQM)%QBOUNDARY ) THEN

            ! . Get the modified gradient vector.
	    MODGVEC = MOPAC_DATA_LA_GRADIENT ( IQM, ATMCRD, GRADQM(1:3,IQM) )

            ! . Add in the contributions to the gradients.
	    GRADIENT(1:3,IATOM)               = GRADIENT(1:3,IATOM)                                 + MODGVEC
	    GRADIENT(1:3,QMATOM(IQM)%PARTNER) = GRADIENT(1:3,QMATOM(IQM)%PARTNER) + GRADQM(1:3,IQM) - MODGVEC

            ! . Add in the contribution to the virial (shifting it to the partner atom).
	    VIRIAL = VIRIAL + DOT_PRODUCT ( ATMCRD(1:3,QMATOM(IQM)%PARTNER), GRADQM(1:3,IQM) )

         ! . The atom is a normal QM atom.
	 ELSE

            ! . Add in the contribution to the gradient.
	    GRADIENT(1:3,IATOM) = GRADIENT(1:3,IATOM) + GRADQM(1:3,IQM)

            ! . Add in the contribution to the virial. This formula is OK as the minimum image convention is not being used.
	    VIRIAL = VIRIAL + DOT_PRODUCT ( ATMCRD(1:3,IATOM), GRADQM(1:3,IQM) )

	 END IF

      END DO

      END SUBROUTINE DERIVATIVES_QM

      !---------------------------
      SUBROUTINE DERIVATIVES_QM_MM
      !---------------------------

      ! . Local scalars.
      INTEGER            :: I, IATOM, II, IINT, IN, IREC, J, JATOM, JJ, MATOM, N1, NINTQ, P, QATOM, QINDX
      REAL ( KIND = DP ) :: PDEN

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3)              :: GVEC, MODGVEC
      REAL ( KIND = DP ), DIMENSION(1:10)             :: DENFAC
      REAL ( KIND = DP ), DIMENSION(1:NBASTR-NPIHELE) :: DENTMP

      ! . Local types.
      TYPE(NBLIST_TYPE), POINTER :: NBCURRENT, NBNEXT

      ! . There are PI atoms.
      IF ( NPIBEADS > 1 ) THEN

         ! . Save the non-PI atom elements of the first density matrix.
         DENTMP(1:NBASTR-NPIHELE) = DENMAT(1:NBASTR-NPIHELE,1)

         ! . Fill the elements of the density matrix with the sum of all the matrices.
         DO I = 1,NBASTR-NPIHELE
            DENMAT(I,1) = SUM ( DENMAT(I,1:NPIBEADS) )
         END DO
         
      END IF

      ! . Initialization.
      DENFAC     = 2.0_DP
      DENFAC(1)  = 1.0_DP
      DENFAC(3)  = 1.0_DP
      DENFAC(6)  = 1.0_DP
      DENFAC(10) = 1.0_DP

      ! . Scale DENFAC by the number of PI polymers.
      DENFAC = DENFAC / REAL ( NPIBEADS, DP )

      ! . Initialize the number of interactions counter.
      IINT  = 0
      NINTQ = 0

      ! . Loop over the records in the non-bonding interaction lists.
      DO IREC = 1,NNBLISTQM

	 ! . Get the next list.
	 IF ( IREC == 1 ) THEN
            NBCURRENT => NBLISTQM_FIRST
            NBNEXT    => NBLISTQM_FIRST%NEXT_LIST
	 ELSE
            NBCURRENT => NBNEXT
            NBNEXT    => NBNEXT%NEXT_LIST
	 END IF

         ! . Get the first atom of the interaction.
	 IATOM = NBCURRENT%ATOM

         ! . Loop over the interactions for the atom.
         DO IN = 1,SIZE(NBCURRENT%INTERACTIONS)

            ! . Get the second atom of the interaction.
            JATOM = NBCURRENT%INTERACTIONS(IN)

            ! . Increment the total number of interactions counter.
	    IINT = IINT + 1

            ! . Determine which is the QM atom. A negative J indicates J is quantum.
            IF ( JATOM > 0 ) THEN
               QATOM =   IATOM
               MATOM =   JATOM
            ELSE
               QATOM = - JATOM
               MATOM =   IATOM
            END IF

            ! . Get some information for the QM atom.
	    QINDX = ATMQMI(QATOM)
            N1    = QMATOM(QINDX)%BFIRST
            P     = QMATOM(QINDX)%COPY
            IF ( P == 0 ) P = 1

            ! . Evaluate the derivatives involving the first E1B integral.
            II   = ( N1 * ( N1 + 1 ) ) / 2
            PDEN = DENFAC(1) * DENMAT(II,P)

            ! . Determine the gradient vector.
	    GVEC(1) = - PDEN * DXE1BA(IINT)
	    GVEC(2) = - PDEN * DYE1BA(IINT)
	    GVEC(3) = - PDEN * DZE1BA(IINT)

            ! . Add in the contribution to the MM atom gradients.
            GRADIENT(1:3,MATOM) = GRADIENT(1:3,MATOM) - GVEC

            ! . Add in the contribution to the QM atom gradients.
	    IF ( QMATOM(QINDX)%QBOUNDARY ) THEN
	       MODGVEC = MOPAC_DATA_LA_GRADIENT ( QINDX, ATMCRD, GVEC )
               GRADIENT(1:3,QATOM)                 = GRADIENT(1:3,QATOM)                        + MODGVEC
               GRADIENT(1:3,QMATOM(QINDX)%PARTNER) = GRADIENT(1:3,QMATOM(QINDX)%PARTNER) + GVEC - MODGVEC
	    ELSE
               GRADIENT(1:3,QATOM) = GRADIENT(1:3,QATOM) + GVEC
	    END IF

            ! . Calculate the contribution to the virial.
	    VIRIAL = VIRIAL + GVEC(1) * DXQM(IINT) + GVEC(2) * DYQM(IINT) + GVEC(3) * DZQM(IINT)

            ! . Non-hydrogen atoms.
            IF ( NATORB(QMATOM(QINDX)%NUMBER) > 1 ) THEN

               !. Increment the interaction number.
               NINTQ = NINTQ + 1

               ! . Evaluate the derivatives for the remaining E1B integrals.
               JJ = 1
               DO I = (N1 + 1),QMATOM(QINDX)%BLAST
                  II = (I * (I - 1)) / 2 + N1 - 1
                  DO J = N1,I
                     II = II + 1
                     JJ = JJ + 1
                     PDEN = DENFAC(JJ) * DENMAT(II,P)

        	     ! . Determine the gradient vector.
		     GVEC(1) = - PDEN * DXE1BB(JJ,NINTQ)
		     GVEC(2) = - PDEN * DYE1BB(JJ,NINTQ)
		     GVEC(3) = - PDEN * DZE1BB(JJ,NINTQ)

        	     ! . Add in the contribution to the MM atom gradients.
        	     GRADIENT(1:3,MATOM) = GRADIENT(1:3,MATOM) - GVEC

        	     ! . Add in the contribution to the QM atom gradients.
		     IF ( QMATOM(QINDX)%QBOUNDARY ) THEN
			MODGVEC = MOPAC_DATA_LA_GRADIENT ( QINDX, ATMCRD, GVEC )
        		GRADIENT(1:3,QATOM)                 = GRADIENT(1:3,QATOM)                        + MODGVEC
        		GRADIENT(1:3,QMATOM(QINDX)%PARTNER) = GRADIENT(1:3,QMATOM(QINDX)%PARTNER) + GVEC - MODGVEC
		     ELSE
        		GRADIENT(1:3,QATOM) = GRADIENT(1:3,QATOM) + GVEC
		     END IF

        	     ! . Calculate the contribution to the virial.
		     VIRIAL = VIRIAL + GVEC(1) * DXQM(IINT) + GVEC(2) * DYQM(IINT) + GVEC(3) * DZQM(IINT)

                  END DO
               END DO
            END IF
         END DO
      END DO

      ! . Restore the density matrix if there are PI atoms.
      IF ( NPIBEADS > 1 ) DENMAT(1:NBASTR-NPIHELE,1) = DENTMP(1:NBASTR-NPIHELE)

      ! . Free the QM/MM interaction derivative arrays.
      IF ( ALLOCATED ( DXE1BA ) ) DEALLOCATE ( DXE1BA )
      IF ( ALLOCATED ( DYE1BA ) ) DEALLOCATE ( DYE1BA )
      IF ( ALLOCATED ( DZE1BA ) ) DEALLOCATE ( DZE1BA )
      IF ( ALLOCATED ( DXE1BB ) ) DEALLOCATE ( DXE1BB )
      IF ( ALLOCATED ( DYE1BB ) ) DEALLOCATE ( DYE1BB )
      IF ( ALLOCATED ( DZE1BB ) ) DEALLOCATE ( DZE1BB )
      IF ( ALLOCATED ( DXQM   ) ) DEALLOCATE ( DXQM   )
      IF ( ALLOCATED ( DYQM   ) ) DEALLOCATE ( DYQM   )
      IF ( ALLOCATED ( DZQM   ) ) DEALLOCATE ( DZQM   )

      END SUBROUTINE DERIVATIVES_QM_MM

   END SUBROUTINE GRADIENTS

   !------------------------------------------------------------------------------------------------------------
   SUBROUTINE ANALYT ( PSUMT, PSUMA, PALPHA, PBETA, NP, COORD, NAT, JFIRST, JLAST, IFIRST, ILAST, ENG, QAM1PM3 )
   !------------------------------------------------------------------------------------------------------------

   ! . Scalar argument declarations.
   INTEGER, INTENT(IN) :: IFIRST, ILAST, JFIRST, JLAST, NP
   LOGICAL, INTENT(IN) :: QAM1PM3

   ! . Array argument declarations.
   INTEGER,            DIMENSION(1:2),        INTENT(IN)  :: NAT
   REAL ( KIND = DP ), DIMENSION(1:3,1:2),    INTENT(IN)  :: COORD
   REAL ( KIND = DP ), DIMENSION(1:3),        INTENT(OUT) :: ENG
   REAL ( KIND = DP ), DIMENSION(1:171),      INTENT(IN)  :: PSUMT
   REAL ( KIND = DP ), DIMENSION(1:171,1:NP), INTENT(IN)  :: PALPHA, PBETA, PSUMA

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: PDDG_EXPONENT = 10.0_DP

   ! . Local scalars.
   INTEGER            :: I, IA, ID, IG, IK, ISP, IX, J, JA, JD, JK, K, KK, KL, L, LL, M, MK, ML, MN, N, NI, NJ, NK, NL
   REAL ( KIND = DP ) :: AA, ANAM1, AX, BB, DFAC, DEL1, DM, F3, PFAC, RIJ, RR, R0, R2, TERMAA, TERMAB, TERMNC, ZAF, ZBF

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:4,1:2) :: CORINT
   REAL ( KIND = DP ), DIMENSION(1:22)    :: DG, G
   REAL ( KIND = DP ), DIMENSION(1:100)   :: DR
   REAL ( KIND = DP ), DIMENSION(1:3)     :: EAA, EAB, ENUC

   ! . Initialize the local derivative arrays.
   EAA  = 0.0_DP
   EAB  = 0.0_DP
   ENG  = 0.0_DP
   ENUC = 0.0_DP

   ! . Initialize some counters.
   JD = JLAST - JFIRST + 1
   JA = 1
   ID = ILAST - IFIRST + 1 + JD
   IA = JD + 1

   ! . Get the interatomic distance.
   I  = 2
   NI = NAT(I)
   J  = 1
   NJ = NAT(J)
   R2 = SUM ( ( COORD(1:3,I) - COORD(1:3,J) )**2 )
   RIJ = SQRT ( R2 )
   R0  = ANGSTROMS_TO_BOHRS * RIJ
   RR  = R0 * R0

   ! . Calculate the two-electron integrals.
   CALL REPP ( NI, NJ, R0, G, CORINT )

   ! . Loop over the Cartesian components.
   DO IX = 1,3
      DEL1   = COORD(IX,I) - COORD(IX,J)
      TERMAA = 0.0_DP
      TERMAB = 0.0_DP
      ISP = 0
      CALL DELRI  ( DG, NI, NJ, R0, DEL1 )
      CALL DELMOL ( COORD, I, J, NI, NJ, IA, ID, JA, JD, IX, RIJ, DEL1, ISP, DG, DR, G )

      ! . Calculate the first derivative of the nuclear repulsion term.
      ! . Skip invalid terms.
      IF ( (RIJ.LT.1.0_DP) .AND. (NATORB(NI)*NATORB(NJ).EQ.0) ) THEN
         TERMNC = 0.0_DP
         GO TO 10
      END IF

      ! . Treat N-H and O-H interactions.
      IF ( (NI.EQ.1) .AND. ( (NJ.EQ.7) .OR. (NJ.EQ.8) ) ) THEN
         F3 = 1.0_DP + EXP(-ALP(1)*RIJ)+RIJ*EXP(-ALP(NJ)*RIJ)
         DFAC = DG(1)*F3-G(1)*(DEL1/RIJ)*(ALP(1)*EXP(-ALP(1)*RIJ)+(ALP(NJ)*RIJ-1.0_DP)*EXP(-ALP(NJ)*RIJ))
         ELSE IF ( ((NI.EQ.7) .OR. (NI.EQ.8)) .AND. (NJ.EQ.1) ) THEN
            F3 = 1.0_DP+EXP(-ALP(1)*RIJ)+RIJ*EXP(-ALP(NI)*RIJ)
            DFAC = DG(1)*F3-G(1)*(DEL1/RIJ)*(ALP(1)*EXP(-ALP(1)*RIJ)+(ALP(NI)*RIJ-1.0_DP)*EXP(-ALP(NI)*RIJ))
      ! . Treat other interactions.
      ELSE
         F3 = 1.0_DP+EXP(-ALP(NI)*RIJ)+EXP(-ALP(NJ)*RIJ)
         DFAC = DG(1)*F3-G(1)*(DEL1/RIJ)*(ALP(NI)*EXP(-ALP(NI)*RIJ)+ALP(NJ)*EXP(-ALP(NJ)*RIJ))
      END IF
      TERMNC = CORE(NI) * CORE(NJ) * DFAC

      ! . Calculate the extra core-core repulsion terms required for an AM1 calculation.
      IF ( QAM1PM3 ) THEN
         ANAM1 = 0.0_DP
         DO IG = 1,4
            IF ( ABS ( FN1(NI,IG) ) .GT. 0.0_DP ) THEN
               ANAM1 = ANAM1+FN1(NI,IG)*(1.0_DP/(RIJ*RIJ) + 2.0_DP*FN2(NI,IG)*(RIJ-FN3(NI,IG))/RIJ)* &
                                               EXP(MAX(-30.0_DP,-FN2(NI,IG)*(RIJ-FN3(NI,IG))**2))
            END IF
            IF ( ABS ( FN1(NJ,IG) ) .GT. 0.0_DP ) THEN
               ANAM1 = ANAM1+FN1(NJ,IG)*(1.0_DP/(RIJ*RIJ) + 2.0_DP*FN2(NJ,IG)*(RIJ-FN3(NJ,IG))/RIJ)* &
                                               EXP(MAX(-30.0_DP,-FN2(NJ,IG)*(RIJ-FN3(NJ,IG))**2))
            END IF
         END DO
         ANAM1  = ANAM1 * CORE(NI) * CORE(NJ)
         TERMNC = TERMNC - ANAM1 * DEL1 / RIJ
      END IF

      ! . PDDG specific terms.
      IF ( HAMILTONIAN == "PDDG" ) THEN
        ANAM1 = 0.0_DP
        ZAF = CORE(NI)/(CORE(NI) + CORE(NJ))
        ZBF = CORE(NJ)/(CORE(NI) + CORE(NJ))
        DO IK = 1,2
           DO JK = 1,2
              DM = RIJ - PDDGE(NI,IK) - PDDGE(NJ,JK)
              AX = PDDG_EXPONENT * DM * DM
!write ( 6, "(2i10,3f20.10)" ) ik, jk, dm, ax, (ZAF * PDDGC(NI,IK) + ZBF * PDDGC(NJ,JK))
              ANAM1 = ANAM1 + (ZAF * PDDGC(NI,IK) + ZBF * PDDGC(NJ,JK)) * &
                      2.0_DP * PDDG_EXPONENT * DM * EXP( -AX ) 
            END DO
         END DO
         TERMNC = TERMNC - ANAM1 * DEL1 / RIJ
      END IF

      10 CONTINUE

      ! . Calculate the core-electron attraction derivatives.
      ! . Atom core I affecting AOs on J.
      ISP = 0
      DO M = JA,JD
         BB = 1.0_DP
         DO N = M,JD
            MN = M+N*(N-1)/2
            ISP = ISP+1
            TERMAB = TERMAB-BB*CORE(NI)*PSUMT(MN)*DR(ISP)
            BB = 2.0_DP
         END DO
      END DO

      ! . Atom core J affecting AOs on I.
      K = MAX(JD-JA+1,1)
      K = (K*(K+1))/2
      ISP = -K+1
      DO M = IA,ID
         BB = 1.0_DP
         DO N = M,ID
            MN = M+N*(N-1)/2
            ISP = ISP+K
            TERMAB = TERMAB-BB*CORE(NJ)*PSUMT(MN)*DR(ISP)
            BB = 2.0_DP
         END DO
      END DO

      ! . Calculate the two electron terms.
      ISP = 0
      DO K = IA,ID
         AA = 1.0_DP
         KK = (K*(K-1))/2
         DO L = K,ID
            LL = (L*(L-1))/2
            DO M = JA,JD
               BB = 1.0_DP
               DO N = M,JD
                  ISP = ISP+1
                  KL = K+LL
                  MN = M+N*(N-1)/2

                  ! . Coulomb term.
                  PFAC   = SUM ( PSUMA(KL,1:NP) * PSUMA(MN,1:NP) )
                  TERMAA = TERMAA + AA * BB * PFAC * DR(ISP)
                  MK = M+KK
                  NK = N+KK
                  ML = M+LL
                  NL = N+LL

                  ! . Exchange term.
                  PFAC   = SUM ( PALPHA(MK,1:NP) * PALPHA(NL,1:NP) + PALPHA(NK,1:NP) * PALPHA(ML,1:NP) + &
		                 PBETA (MK,1:NP) * PBETA (NL,1:NP) + PBETA (NK,1:NP) * PBETA (ML,1:NP)   )
                  TERMAA = TERMAA - 0.5_DP * AA * BB * DR(ISP) * PFAC
                  BB = 2.0_DP

               END DO
            END DO
            AA = 2.0_DP
         END DO
      END DO

      ! . Accumulate the results.
      EAA(IX)  = TERMAA
      EAB(IX)  = TERMAB
      ENUC(IX) = TERMNC

   END DO

   ! . Sum the total contributions to the derivatives.
   ENG = EV_TO_KJ * ( EAA + EAB + REAL ( NP, DP ) * ENUC )

   !============================================================================
   CONTAINS
   !============================================================================

      !---------------------------------------------------------------------------------------
      SUBROUTINE DELMOL ( COORD, I, J, NI, NJ, IA, ID, JA, JD, IX, RIJ, DEL1, ISP, DG, DR, G )
      !---------------------------------------------------------------------------------------

      ! . Scalar argument declarations.
      INTEGER,            INTENT(IN)    :: I, IA, ID, IX, J, JA, JD, NI, NJ
      INTEGER,            INTENT(INOUT) :: ISP
      REAL ( KIND = DP ), INTENT(IN)    :: DEL1, RIJ

      ! . Array argument declarations.
      REAL ( KIND = DP ), DIMENSION(1:3,1:2), INTENT(IN)    :: COORD
      REAL ( KIND = DP ), DIMENSION(1:22),    INTENT(IN)    :: DG, G
      REAL ( KIND = DP ), DIMENSION(1:100),   INTENT(INOUT) :: DR

      ! . Local scalars.
      INTEGER            :: IB, JB, K, KK, L, LL, M, MM, N, NN
      REAL ( KIND = DP ) :: TEMP1, TEMP2

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: TDX, TDY, TDZ, TX, TY, TZ

      IF ( ( NI.GT.1 ) .OR. ( NJ.GT.1) ) THEN
         CALL ROTAT ( COORD, I, J, IX, RIJ, DEL1, 2, TDX, TDY, TDZ, TX, TY, TZ )
      ENDIF
      IB=MAX(IA,ID)
      JB=MAX(JA,JD)
      DO K=IA,IB
         KK=K-IA
         DO L=K,IB
            LL=L-IA
            DO M=JA,JB
               MM=M-JA
               DO N=M,JB
                  NN=N-JA
                  ISP=ISP+1
                  IF(NN.EQ.0)THEN
                     IF(LL.EQ.0) THEN
                        ! . (SS/SS)
                        DR(ISP)=DG(1)
                     ELSEIF(KK.EQ.0) THEN
                        ! . (SP/SS)
                        DR(ISP)=DG(2)*TX(LL)+G(2)*TDX(LL)
                     ELSE
                        ! . (PP/SS)
                        DR(ISP)=DG(3)*TX(KK)*TX(LL)+G(3)*(TDX(KK)*TX(LL)+TX(KK)*TDX(LL)) &
                               +DG(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))+G(4)*(TDY(KK)*TY(LL)+TY(KK)*TDY(LL) &
                               +TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))
                     ENDIF
                  ELSEIF(MM.EQ.0) THEN
                     IF(LL.EQ.0) THEN
                        ! . (SS/SP)
                        DR(ISP)=DG(5)*TX(NN)+G(5)*TDX(NN)
                     ELSEIF(KK.EQ.0) THEN
                        ! . (SP/SP)
                        DR(ISP)=DG(6)*TX(LL)*TX(NN)+G(6)*(TDX(LL)*TX(NN)+TX(LL)*TDX(NN)) &
                               +DG(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+G(7)*(TDY(LL)*TY(NN)+TY(LL)*TDY(NN) &
                               +TDZ(LL)*TZ(NN)+TZ(LL)*TDZ(NN))
                     ELSE
                        ! . (PP/SP)
                        DR(ISP)=DG(8)*TX(KK)*TX(LL)*TX(NN)+G(8)*(TDX(KK)*TX(LL)*TX(NN)+TX(KK)*TDX(LL)*TX(NN) &
                   +TX(KK)*TX(LL)*TDX(NN))+DG(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN) &
             +G(9)*((TDY(KK)*TY(LL)+TY(KK)*TDY(LL)+TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))*TX(NN) &
                   +(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TDX(NN))+DG(10)*(TX(KK)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                     +TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN)))+G(10)*(TDX(KK)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                    +TDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(KK)*(TDY(LL)*TY(NN)+TY(LL)*TDY(NN) &
                            +TDZ(LL)*TZ(NN)+TZ(LL)*TDZ(NN))+TX(LL)*(TDY(KK)*TY(NN)+TY(KK)*TDY(NN) &
                            +TDZ(KK)*TZ(NN)+TZ(KK)*TDZ(NN)))
                     ENDIF
                  ELSEIF(LL.EQ.0) THEN
                     ! . (SS/PP)
                     DR(ISP)=DG(11)*TX(MM)*TX(NN)+G(11)*(TDX(MM)*TX(NN)+TX(MM)*TDX(NN)) &
             +DG(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+G(12)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN) &
                    +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN))
                  ELSEIF(KK.EQ.0) THEN
                     ! . (SP/PP)
                     DR(ISP)=DG(13)*TX(LL)*TX(MM)*TX(NN)+G(13)*(TDX(LL)*TX(MM)*TX(NN)+TX(LL)*TDX(MM)*TX(NN) &
                    +TX(LL)*TX(MM)*TDX(NN))+DG(14)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN)) &
             +G(14)*(TDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(LL)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN) &
                            +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN)))+DG(15)*(TY(LL)*(TY(MM)*TX(NN)+TY(NN)*TX(MM)) &
                     +TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM)))+G(15)*(TDY(LL)*(TY(MM)*TX(NN)+TY(NN)*TX(MM)) &
                    +TDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))+TY(LL)*(TDY(MM)*TX(NN)+TY(MM)*TDX(NN) &
                            +TDY(NN)*TX(MM)+TY(NN)*TDX(MM))+TZ(LL)*(TDZ(MM)*TX(NN)+TZ(MM)*TDX(NN) &
                            +TDZ(NN)*TX(MM)+TZ(NN)*TDX(MM)))
                  ELSE
                     ! . (PP/PP)
                     DR(ISP)=DG(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)+G(16)*(TDX(KK)*TX(LL)*TX(MM)*TX(NN) &
                    +TX(KK)*TDX(LL)*TX(MM)*TX(NN)+TX(KK)*TX(LL)*TDX(MM)*TX(NN) &
                    +TX(KK)*TX(LL)*TX(MM)*TDX(NN))+DG(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN) &
             +G(17)*((TDY(KK)*TY(LL)+TY(KK)*TDY(LL)+TDZ(KK)*TZ(LL)+TZ(KK)*TDZ(LL))*TX(MM)*TX(NN) &
                    +(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*(TDX(MM)*TX(NN)+TX(MM)*TDX(NN))) &
             +DG(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+G(18)*((TDX(KK)*TX(LL)+TX(KK)*TDX(LL)) &
                       *(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(KK)*TX(LL)*(TDY(MM)*TY(NN)+TY(MM)*TDY(NN) &
                                   +TDZ(MM)*TZ(NN)+TZ(MM)*TDZ(NN)))
                     DR(ISP)=DR(ISP)+DG(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN)) &
             +G(19)*(TDY(KK)*TY(LL)*TY(MM)*TY(NN)+TY(KK)*TDY(LL)*TY(MM)*TY(NN)+TY(KK)*TY(LL)*TDY(MM)*TY(NN) &
                       +TY(KK)*TY(LL)*TY(MM)*TDY(NN)+TDZ(KK)*TZ(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TDZ(LL)*TZ(MM)*TZ(NN) &
                       +TZ(KK)*TZ(LL)*TDZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TDZ(NN)) &
             +DG(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)*TZ(MM))) &
                        +TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
                     ! . To avoid compiler difficulties the expression is divided.
                     TEMP1 = TDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))  &
                            +TDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))  &
                            +TX(KK)*(TDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TDX(NN)*(TY(LL)*TY(MM)+TZ(LL)*TZ(MM))) &
                            +TX(LL)*(TDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TDX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))
                     TEMP2 = TX(KK)*(TX(MM)*(TDY(LL)*TY(NN)+TY(LL)*TDY(NN)+TDZ(LL)*TZ(NN)+TZ(LL)*TDZ(NN)) &
                            +TX(NN)*(TDY(LL)*TY(MM)+TY(LL)*TDY(MM)+TDZ(LL)*TZ(MM)+TZ(LL)*TDZ(MM))) &
                            +TX(LL)*(TX(MM)*(TDY(KK)*TY(NN)+TY(KK)*TDY(NN)+TDZ(KK)*TZ(NN)+TZ(KK)*TDZ(NN)) &
                            +TX(NN)*(TDY(KK)*TY(MM)+TY(KK)*TDY(MM)+TDZ(KK)*TZ(MM)+TZ(KK)*TDZ(MM)))
                     DR(ISP)=DR(ISP)+G(20)*(TEMP1+TEMP2)
                     DR(ISP)=DR(ISP)+DG(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TY(MM)*TY(NN)) &
             +G(21)*(TDY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TY(KK)*TDY(LL)*TZ(MM)*TZ(NN)+TY(KK)*TY(LL)*TDZ(MM)*TZ(NN) &
                       +TY(KK)*TY(LL)*TZ(MM)*TDZ(NN)+TDZ(KK)*TZ(LL)*TY(MM)*TY(NN)+TZ(KK)*TDZ(LL)*TY(MM)*TY(NN) &
                       +TZ(KK)*TZ(LL)*TDY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TY(MM)*TDY(NN))
                     DR(ISP)=DR(ISP)+DG(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN)) &
             +G(22)*((TDY(KK)*TZ(LL)+TY(KK)*TDZ(LL)+TDZ(KK)*TY(LL)+TZ(KK)*TDY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN)) &
                       +(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TDY(MM)*TZ(NN)+TY(MM)*TDZ(NN)+TDZ(MM)*TY(NN)+TZ(MM)*TDY(NN)))
                  END IF
               END DO
            END DO
         END DO
      END DO

      END SUBROUTINE DELMOL

      !----------------------------------------
      SUBROUTINE DELRI ( DG, NI, NJ, RR, DEL1 )
      !----------------------------------------

      ! . Scalar argument declarations.
      INTEGER,            INTENT(IN) :: NI, NJ
      REAL ( KIND = DP ), INTENT(IN) :: DEL1, RR

      ! . Array argument declarations.
      REAL ( KIND = DP ), DIMENSION(1:22), INTENT(OUT) :: DG

      ! . Local scalars.
      REAL ( KIND = DP ) :: ADD, ADE, ADQ, AED, AEE, AEQ, AQD, AQE, AQQ, DA, DB, &
                            DXDX, DXQXZ, DZDZ, DZE, DZQXX, DZQZZ, EDZ, EE,       &
                            EQXX, EQZZ, QA, QB, QXXDZ, QXXE, QXXQXX, QXXQYY,     &
                            QXXQZZ, QXYQXY, QXZDX, QXZQXZ, QZZDZ, QZZE, QZZQXX,  &
                            QZZQZZ, TERM

      TERM = ( ANGSTROMS_TO_BOHRS * ANGSTROMS_TO_BOHRS * AU_TO_EV * DEL1 ) / RR
      DA = DD(NI)
      DB = DD(NJ)
      QA = QQ(NI)
      QB = QQ(NJ)
      ! . Hydrogen - hydrogen.
      AEE=0.25_DP*(1.0_DP/AM(NI)+1.0_DP/AM(NJ))**2
      EE    =-RR/(SQRT(RR**2+AEE))**3
      DG(1)=TERM*EE
      IF(NATORB(NI).LE.2.AND.NATORB(NJ).LE.2) RETURN
      IF(NATORB(NI).LE.2) GO TO 10
      ! . Heavy atom - hydrogen.
      ADE=0.25_DP*(1.0_DP/AD(NI)+1.0_DP/AM(NJ))**2
      AQE=0.25_DP*(1.0_DP/AQ(NI)+1.0_DP/AM(NJ))**2
      DZE   = (RR+DA)/(SQRT((RR+DA)**2+ADE))**3-(RR-DA)/(SQRT((RR-DA)**2+ADE))**3
      QZZE  =-(RR+2.0_DP*QA)/(SQRT((RR+2.0_DP*QA)**2+AQE))**3-(RR-2.0_DP*QA)/(SQRT((RR-2.0_DP*QA)**2+AQE))**3 &
             +(2.0_DP*RR)/(SQRT(RR**2+AQE))**3
      QXXE  =-(2.0_DP*RR)/(SQRT(RR**2+4.0_DP*QA**2+AQE))**3+(2.0_DP*RR)/(SQRT(RR**2+AQE))**3
      DG(2)=-(TERM*DZE)/2.0_DP
      DG(3)=TERM*(EE+QZZE/4.0_DP)
      DG(4)=TERM*(EE+QXXE/4.0_DP)
      IF(NATORB(NJ).LE.2) RETURN
      ! . Hydrogen - Heavy atom.
   10 AED=0.25_DP*(1.0_DP/AM(NI)+1.0_DP/AD(NJ))**2
      AEQ=0.25_DP*(1.0_DP/AM(NI)+1.0_DP/AQ(NJ))**2
      EDZ   = (RR-DB)/(SQRT((RR-DB)**2+AED))**3-(RR+DB)/(SQRT((RR+DB)**2+AED))**3
      EQZZ  =-(RR-2.0_DP*QB)/(SQRT((RR-2.0_DP*QB)**2+AEQ))**3-(RR+2.0_DP*QB)/(SQRT((RR+2.0_DP*QB)**2+AEQ))**3 &
             +(2.0_DP*RR)/(SQRT(RR**2+AEQ))**3
      EQXX  =-(2.0_DP*RR)/(SQRT(RR**2+4.0_DP*QB**2+AEQ))**3+(2.0_DP*RR)/(SQRT(RR**2+AEQ))**3
      DG(5)=-(TERM*EDZ)/2.0_DP
      DG(11)=TERM*(EE+EQZZ/4.0_DP)
      DG(12)=TERM*(EE+EQXX/4.0_DP)
      IF(NATORB(NI).LE.2) RETURN
      ! . Heavy atom - Heavy atom.
      ADD=0.25_DP*(1.0_DP/AD(NI)+1.0_DP/AD(NJ))**2
      ADQ=0.25_DP*(1.0_DP/AD(NI)+1.0_DP/AQ(NJ))**2
      AQD=0.25_DP*(1.0_DP/AQ(NI)+1.0_DP/AD(NJ))**2
      AQQ=0.25_DP*(1.0_DP/AQ(NI)+1.0_DP/AQ(NJ))**2
      DXDX  =-(2.0_DP*RR)/(SQRT(RR**2+(DA-DB)**2+ADD))**3+(2.0_DP*RR)/(SQRT(RR**2+(DA+DB)**2+ADD))**3
      DZDZ  =-(RR+DA-DB)/(SQRT((RR+DA-DB)**2+ADD))**3-(RR-DA+DB)/(SQRT((RR-DA+DB)**2+ADD))**3 &
             +(RR-DA-DB)/(SQRT((RR-DA-DB)**2+ADD))**3+(RR+DA+DB)/(SQRT((RR+DA+DB)**2+ADD))**3
      DZQXX = 2.0_DP*(RR+DA)/(SQRT((RR+DA)**2+4.0_DP*QB**2+ADQ))**3-2.0_DP*(RR-DA)/(SQRT((RR-DA)**2+4.0_DP*QB**2+ADQ))**3 &
             -2.0_DP*(RR+DA)/(SQRT((RR+DA)**2+ADQ))**3+2.0_DP*(RR-DA)/(SQRT((RR-DA)**2+ADQ))**3
      QXXDZ = 2.0_DP*(RR-DB)/(SQRT((RR-DB)**2+4.0_DP*QA**2+AQD))**3-2.0_DP*(RR+DB)/(SQRT((RR+DB)**2+4.0_DP*QA**2+AQD))**3 &
             -2.0_DP*(RR-DB)/(SQRT((RR-DB)**2+AQD))**3+2.0_DP*(RR+DB)/(SQRT((RR+DB)**2+AQD))**3
      DZQZZ = (RR+DA-2.0_DP*QB)/(SQRT((RR+DA-2.0_DP*QB)**2+ADQ))**3-(RR-DA-2.0_DP*QB)/(SQRT((RR-DA-2.0_DP*QB)**2+ADQ))**3 &
             +(RR+DA+2.0_DP*QB)/(SQRT((RR+DA+2.0_DP*QB)**2+ADQ))**3-(RR-DA+2.0_DP*QB)/(SQRT((RR-DA+2.0_DP*QB)**2+ADQ))**3 &
             +2.0_DP*(RR-DA)/(SQRT((RR-DA)**2+ADQ))**3-2.0_DP*(RR+DA)/(SQRT((RR+DA)**2+ADQ))**3
      QZZDZ = (RR+2.0_DP*QA-DB)/(SQRT((RR+2.0_DP*QA-DB)**2+AQD))**3-(RR+2.0_DP*QA+DB)/(SQRT((RR+2.0_DP*QA+DB)**2+AQD))**3 &
             +(RR-2.0_DP*QA-DB)/(SQRT((RR-2.0_DP*QA-DB)**2+AQD))**3-(RR-2.0_DP*QA+DB)/(SQRT((RR-2.0_DP*QA+DB)**2+AQD))**3 &
             -2.0_DP*(RR-DB)/(SQRT((RR-DB)**2+AQD))**3+2.0_DP*(RR+DB)/(SQRT((RR+DB)**2+AQD))**3
      QXXQXX=-(2.0_DP*RR)/(SQRT(RR**2+4.0_DP*(QA-QB)**2+AQQ))**3-(2.0_DP*RR)/(SQRT(RR**2+4.0_DP*(QA+QB)**2+AQQ))**3 &
             +(4.0_DP*RR)/(SQRT(RR**2+4.0_DP*QA**2+AQQ))**3+(4.0_DP*RR)/(SQRT(RR**2+4.0_DP*QB**2+AQQ))**3 &
             -(4.0_DP*RR)/(SQRT(RR**2+AQQ))**3
      QXXQYY=-(4.0_DP*RR)/(SQRT(RR**2+4.0_DP*QA**2+4.0_DP*QB**2+AQQ))**3+(4.0_DP*RR)/(SQRT(RR**2+4.0_DP*QA**2+AQQ))**3 &
             +(4.0_DP*RR)/(SQRT(RR**2+4.0_DP*QB**2+AQQ))**3-(4.0_DP*RR)/(SQRT(RR**2+AQQ))**3
      QXXQZZ= &
           -2.0_DP*(RR-2.0_DP*QB)/(SQRT((RR-2.0_DP*QB)**2+4.0_DP*QA**2+AQQ))**3-&
            2.0_DP*(RR+2.0_DP*QB)/(SQRT((RR+2.0_DP*QB)**2+4.0_DP*QA**2+AQQ))**3 &
             +2.0_DP*(RR-2.0_DP*QB)/(SQRT((RR-2.0_DP*QB)**2+AQQ))**3+2.0_DP*(RR+2.0_DP*QB)/(SQRT((RR+2.0_DP*QB)**2+AQQ))**3 &
             +(4.0_DP*RR)/(SQRT(RR**2+4.0_DP*QA**2+AQQ))**3-(4.0_DP*RR)/(SQRT(RR**2+AQQ))**3
      QZZQXX= &
           -2.0_DP*(RR+2.0_DP*QA)/(SQRT((RR+2.0_DP*QA)**2+4.0_DP*QB**2+AQQ))**3-&
            2.0_DP*(RR-2.0_DP*QA)/(SQRT((RR-2.0_DP*QA)**2+4.0_DP*QB**2+AQQ))**3 &
             +2.0_DP*(RR+2.0_DP*QA)/(SQRT((RR+2.0_DP*QA)**2+AQQ))**3+2.0_DP*(RR-2.0_DP*QA)/(SQRT((RR-2.0_DP*QA)**2+AQQ))**3 &
             +(4.0_DP*RR)/(SQRT(RR**2+4.0_DP*QB**2+AQQ))**3-(4.0_DP*RR)/(SQRT(RR**2+AQQ))**3
      QZZQZZ= &
           -(RR+2.0_DP*QA-2.0_DP*QB)/(SQRT((RR+2.0_DP*QA-2.0_DP*QB)**2+AQQ))**3-&
            (RR+2.0_DP*QA+2.0_DP*QB)/(SQRT((RR+2.0_DP*QA+2.0_DP*QB)**2+AQQ))**3 &
           -(RR-2.0_DP*QA-2.0_DP*QB)/(SQRT((RR-2.0_DP*QA-2.0_DP*QB)**2+AQQ))**3-&
            (RR-2.0_DP*QA+2.0_DP*QB)/(SQRT((RR-2.0_DP*QA+2.0_DP*QB)**2+AQQ))**3 &
             +2.0_DP*(RR-2.0_DP*QA)/(SQRT((RR-2.0_DP*QA)**2+AQQ))**3+2.0_DP*(RR+2.0_DP*QA)/(SQRT((RR+2.0_DP*QA)**2+AQQ))**3 &
             +2.0_DP*(RR-2.0_DP*QB)/(SQRT((RR-2.0_DP*QB)**2+AQQ))**3+2.0_DP*(RR+2.0_DP*QB)/(SQRT((RR+2.0_DP*QB)**2+AQQ))**3 &
             -(4.0_DP*RR)/(SQRT(RR**2+AQQ))**3
      DXQXZ = 2.0_DP*(RR-QB)/(SQRT((RR-QB)**2+(DA-QB)**2+ADQ))**3-2.0_DP*(RR+QB)/(SQRT((RR+QB)**2+(DA-QB)**2+ADQ))**3 &
             -2.0_DP*(RR-QB)/(SQRT((RR-QB)**2+(DA+QB)**2+ADQ))**3+2.0_DP*(RR+QB)/(SQRT((RR+QB)**2+(DA+QB)**2+ADQ))**3
      QXZDX = 2.0_DP*(RR+QA)/(SQRT((RR+QA)**2+(QA-DB)**2+AQD))**3-2.0_DP*(RR-QA)/(SQRT((RR-QA)**2+(QA-DB)**2+AQD))**3 &
             -2.0_DP*(RR+QA)/(SQRT((RR+QA)**2+(QA+DB)**2+AQD))**3+2.0_DP*(RR-QA)/(SQRT((RR-QA)**2+(QA+DB)**2+AQD))**3
      QXYQXY=-(4.0_DP*RR)/(SQRT(RR**2+2.0_DP*(QA-QB)**2+AQQ))**3-(4.0_DP*RR)/(SQRT(RR**2+2.0_DP*(QA+QB)**2+AQQ))**3 &
             +(8.0_DP*RR)/(SQRT(RR**2+2.0_DP*(QA**2+QB**2)+AQQ))**3
      QXZQXZ=-2.0_DP*(RR+QA-QB)/(SQRT((RR+QA-QB)**2+(QA-QB)**2+AQQ))**3+2.0_DP*(RR+QA+QB)/(SQRT((RR+QA+QB)**2+(QA-QB)**2+AQQ))**3 &
             +2.0_DP*(RR-QA-QB)/(SQRT((RR-QA-QB)**2+(QA-QB)**2+AQQ))**3-2.0_DP*(RR-QA+QB)/(SQRT((RR-QA+QB)**2+(QA-QB)**2+AQQ))**3 &
             +2.0_DP*(RR+QA-QB)/(SQRT((RR+QA-QB)**2+(QA+QB)**2+AQQ))**3-2.0_DP*(RR+QA+QB)/(SQRT((RR+QA+QB)**2+(QA+QB)**2+AQQ))**3 &
             -2.0_DP*(RR-QA-QB)/(SQRT((RR-QA-QB)**2+(QA+QB)**2+AQQ))**3+2.0_DP*(RR-QA+QB)/(SQRT((RR-QA+QB)**2+(QA+QB)**2+AQQ))**3
      DG(6)=(TERM*DZDZ)/4.0_DP
      DG(7)=(TERM*DXDX)/4.0_DP
      DG(8)=-TERM*(EDZ/2.0_DP+QZZDZ/8.0_DP)
      DG(9)=-TERM*(EDZ/2.0_DP+QXXDZ/8.0_DP)
      DG(10)=-(TERM*QXZDX)/8.0_DP
      DG(13)=-TERM*(DZE/2.0_DP+DZQZZ/8.0_DP)
      DG(14)=-TERM*(DZE/2.0_DP+DZQXX/8.0_DP)
      DG(15)=-(TERM*DXQXZ)/8.0_DP
      DG(16)=TERM*(EE+EQZZ/4.0_DP+QZZE/4.0_DP+QZZQZZ/16.0_DP)
      DG(17)=TERM*(EE+EQZZ/4.0_DP+QXXE/4.0_DP+QXXQZZ/16.0_DP)
      DG(18)=TERM*(EE+EQXX/4.0_DP+QZZE/4.0_DP+QZZQXX/16.0_DP)
      DG(19)=TERM*(EE+EQXX/4.0_DP+QXXE/4.0_DP+QXXQXX/16.0_DP)
      DG(20)=(TERM*QXZQXZ)/16.0_DP
      DG(21)=TERM*(EE+EQXX/4.0_DP+QXXE/4.0_DP+QXXQYY/16.0_DP)
      DG(22)=TERM*(QXXQXX-QXXQYY)/32.0_DP

      END SUBROUTINE DELRI

      !------------------------------------------------------------------------------
      SUBROUTINE ROTAT ( COORD, I, J, IX, RIJ, DEL1, IDX, TDX, TDY, TDZ, TX, TY, TZ )
      !------------------------------------------------------------------------------

      ! . Scalar argument declarations.
      INTEGER,            INTENT(IN) :: I, IDX, IX, J
      REAL ( KIND = DP ), INTENT(IN) :: RIJ, DEL1

      ! . Array argument declarations.
      REAL ( KIND = DP ), DIMENSION(1:3,1:2), INTENT(IN)  :: COORD
      REAL ( KIND = DP ), DIMENSION(1:3),     INTENT(OUT) :: TDX, TDY, TDZ, TX, TY, TZ

      ! . Local variables.
      INTEGER            :: IJK
      REAL ( KIND = DP ) :: RXY, RYZ, RZX, TERM, XD, YD, ZD

      XD=COORD(1,I)-COORD(1,J)
      YD=COORD(2,I)-COORD(2,J)
      ZD=COORD(3,I)-COORD(3,J)
      RXY=SQRT(XD*XD+YD*YD)
      RYZ=SQRT(YD*YD+ZD*ZD)
      RZX=SQRT(ZD*ZD+XD*XD)
      DO IJK=1,3
         TX(IJK)=0.0_DP
         TY(IJK)=0.0_DP
         TZ(IJK)=0.0_DP
         TDX(IJK)=0.0_DP
         TDY(IJK)=0.0_DP
         TDZ(IJK)=0.0_DP
      END DO
      IF(RXY.LT.1.0E-4_DP) THEN
         !   MOLECULAR Z AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(3)=1.0_DP
         IF(ZD.LT.0.0_DP) TX(3)=-1.0_DP
         TY(2)=1.0_DP
         TZ(1)=TX(3)
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.1) TDX(1)=1.0_DP/RIJ
         IF(IX.EQ.2) TDX(2)=1.0_DP/RIJ
         IF(IX.EQ.1) TDZ(3)=-1.0_DP/RIJ
         IF(IX.EQ.2) TDY(3)=-TX(3)/RIJ
      ELSEIF(RYZ.LT.1.0E-4_DP) THEN
         !   MOLECULAR X AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=1.0_DP
         IF(XD.LT.0.0_DP) TX(1)=-1.0_DP
         TY(2)=TX(1)
         TZ(3)=1.0_DP
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.2) TDX(2)=1.0_DP/RIJ
         IF(IX.EQ.3) TDX(3)=1.0_DP/RIJ
         IF(IX.EQ.2) TDY(1)=-1.0_DP/RIJ
         IF(IX.EQ.3) TDZ(1)=-TX(1)/RIJ
      ELSEIF(RZX.LT.1.0E-4_DP) THEN
         !   MOLECULAR Y AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(2)=1.0_DP
         IF(YD.LT.0.0_DP) TX(2)=-1.0_DP
         TY(1)=-TX(2)
         TZ(3)=1.0_DP
         IF(IDX.EQ.1) RETURN
         IF(IX.EQ.1) TDX(1)=1.0_DP/RIJ
         IF(IX.EQ.3) TDX(3)=1.0_DP/RIJ
         IF(IX.EQ.1) TDY(2)=1.0_DP/RIJ
         IF(IX.EQ.3) TDZ(2)=-TX(2)/RIJ
      ELSE
         TX(1)=XD/RIJ
         TX(2)=YD/RIJ
         TX(3)=ZD/RIJ
         TZ(3)=RXY/RIJ
         TY(1)=-TX(2)*SIGN(+1.0_DP,TX(1))/TZ(3)
         TY(2)=ABS(TX(1)/TZ(3))
         TY(3)=0.0_DP
         TZ(1)=-TX(1)*TX(3)/TZ(3)
         TZ(2)=-TX(2)*TX(3)/TZ(3)
         IF(IDX.EQ.1) RETURN
         TERM=DEL1/(RIJ*RIJ)
         IF(IX.EQ.1)THEN
            TDX(1)=1.0_DP/RIJ-TX(1)*TERM
            TDX(2)=-TX(2)*TERM
            TDX(3)=-TX(3)*TERM
            TDZ(3)=TX(1)/RXY-TZ(3)*TERM
         ELSEIF(IX.EQ.2) THEN
            TDX(1)=-TX(1)*TERM
            TDX(2)=1.0_DP/RIJ-TX(2)*TERM
            TDX(3)=-TX(3)*TERM
            TDZ(3)=TX(2)/RXY-TZ(3)*TERM
         ELSEIF(IX.EQ.3)THEN
            TDX(1)=-TX(1)*TERM
            TDX(2)=-TX(2)*TERM
            TDX(3)=1.0_DP/RIJ-TX(3)*TERM
            TDZ(3)=-TZ(3)*TERM
         ENDIF
         TDY(1)=-TDX(2)/TZ(3)+TX(2)*TDZ(3)/TZ(3)**2
         IF(TX(1).LT.0.0_DP) TDY(1)=-TDY(1)
         TDY(2)=TDX(1)/TZ(3)-TX(1)*TDZ(3)/TZ(3)**2
         IF(TX(1).LT.0.0_DP) TDY(2)=-TDY(2)
         TDY(3)=0.0_DP
         TDZ(1)=-TX(3)*TDX(1)/TZ(3)-TX(1)*TDX(3)/TZ(3)+TX(1)*TX(3)*TDZ(3)/TZ(3)**2
         TDZ(2)=-TX(3)*TDX(2)/TZ(3)-TX(2)*TDX(3)/TZ(3)+TX(2)*TX(3)*TDZ(3)/TZ(3)**2
      END IF

      END SUBROUTINE ROTAT

   END SUBROUTINE ANALYT

END MODULE MOPAC_GRADIENTS
