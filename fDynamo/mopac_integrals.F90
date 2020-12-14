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
!                           The MOPAC Integrals Module
!===============================================================================
!
! . Subroutines:
!
!   INTEGRALS                    Calculate the MOPAC integrals.
!   INTEGRALS_HCORE              Calculate HCORE for the nth QM calculation.
!   REPP                         Needed for the gradient calculation.
!
! . SEINTC subroutines (public):
!
!   SEINTC_INTEGRALAB
!   SEINTC_INTEGRALB
!   SEINTC_SETUP
!
! . SEINTC subroutines (private):
!
!   SEINTC_TNAB
!   SEINTC_TNB
!
! . QM/MM arrays:
!
!   DXE1BA, DYE1BA, DZE1BA,        Arrays used for storing integrals for
!   DXE1BB, DYE1BB, DZE1BB,        calculation of the QM/MM derivatives.
!   DXQM,   DYQM,   DZQM
!
!===============================================================================
MODULE MOPAC_INTEGRALS

! . Module declarations.
USE CONSTANTS,          ONLY : ANGSTROMS_TO_BOHRS, AU_TO_EV, EV_TO_KJ
USE DEFINITIONS,        ONLY : DP
USE PRINTING,           ONLY : PRINT_ERROR

USE ATOMS,              ONLY : ATMCRD, ATMQMI, NATOMS, NATOMSMM, NATOMSQM
USE ENERGY_NON_BONDING, ONLY : CUT_OFF, CUT_ON, NBLIST_TYPE, NBLISTQM_FIRST, NNBLISTQM, NQMINT, QIMAGE
USE GAUSSIAN_BASIS,     ONLY : GAUSSIAN_OVERLAP
USE MM_TERMS,           ONLY : ATMCHG
USE MOPAC_DATA,         ONLY : BETA, HAMILTONIAN, HCORE, HSTORE, MOPAC_DATA_ATOM_POSITION, MOPAC_DATA_LA_GRADIENT, &
                               NBASIS, NBASTR, NPIATOMS, NPIBASIS, NPIBEADS, NPIHELE, NPI2ELEC, N2ELEC, PIHCORE,   &
			       PIHDIAG, PISETEI, QMATOM, SETEI, USPD
USE MOPAC_PARAMETERS,   ONLY : AD, ALP, AM, AQ, CORE, DD, FN1, FN2, FN3, NATORB, PDDGE, PDDGC, QQ
USE SYMMETRY,           ONLY : BOXL

IMPLICIT NONE
#ifdef	HACK_LAMB
PUBLIC  :: MOPAC_LAMBDA_VALUE, MOPAC_LAMBDA_DELTA
#else
PUBLIC
#endif
PRIVATE :: CUTOFFB, CUTONB, GAMMAB, R2OFFB, R2ONB, T0FAC, T2FAC, T4FAC, T6FAC
#ifndef PGPC
SAVE
#endif

#ifdef	HACK_LAMB
REAL( KIND = DP ) :: MOPAC_LAMBDA_VALUE = 1.0_DP
REAL( KIND = DP ) :: MOPAC_LAMBDA_DELTA = 0.0_DP
#endif

! . Module scalars (for the SEINTC subroutines).
REAL ( KIND = DP ) :: CUTOFFB = 0.0_DP, CUTONB = 0.0_DP, GAMMAB = 1.0_DP, R2OFFB = 0.0_DP, R2ONB = 0.0_DP, &
                      T0FAC   = 0.0_DP, T2FAC  = 0.0_DP, T4FAC  = 0.0_DP, T6FAC  = 0.0_DP

! . Module arrays.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: DXE1BA, DYE1BA, DZE1BA, DXQM, DYQM, DZQM
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: DXE1BB, DYE1BB, DZE1BB

!===============================================================================
CONTAINS
!===============================================================================

   !--------------------------------------------------
   SUBROUTINE INTEGRALS ( ENUCLEAR, VIRIAL, GRADIENT )
   !--------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(INOUT) :: VIRIAL

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT)             :: ENUCLEAR
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: GRADIENT

   ! . Local scalars.
   INTEGER            :: FLAG, I, IBASIS, ICOPY, IFIRST, IINT, INDTEI, ILAST, IQM, IREC, J, JBASIS, JCOPY, &
                         JFIRST, JLAST, JQM, NI, NJ, NNOTRI, NQMINTH, P, QINDX
   LOGICAL            :: QAM1PM3, QGRADIENT
   REAL ( KIND = DP ) :: BETAIJ, ENUC, ENUCLR

   ! . Local arrays.
   INTEGER,            DIMENSION(1:NPIBEADS) :: PIINDTEI
   REAL ( KIND = DP ), DIMENSION(1:NPIBEADS) :: PIENUCLR
   REAL ( KIND = DP ), DIMENSION(1:3)        :: RI, RJ
   REAL ( KIND = DP ), DIMENSION(1:10)       :: E1B, E2A

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: OVERLAP
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: PIOVERLAP

   ! . Local types.
   TYPE(NBLIST_TYPE), POINTER :: NBCURRENT, NBNEXT

   !----------------------------------------------------------------------------
   ! . Setup.
   !----------------------------------------------------------------------------
   ! . Set the AM1/PM3 flag.
   QAM1PM3 = ( HAMILTONIAN == "AM1 " ) .OR. ( HAMILTONIAN == "PDDG " ) .OR. ( HAMILTONIAN == "PM3 " ) .OR. ( HAMILTONIAN == "RM1 " )

   ! . Set the derivative flag.
   QGRADIENT = PRESENT ( GRADIENT )

   ! . Allocate the integral arrays.
   IF ( ALLOCATED ( HCORE   ) ) DEALLOCATE ( HCORE   )
   IF ( ALLOCATED ( HSTORE  ) ) DEALLOCATE ( HSTORE  )
   IF ( ALLOCATED ( PIHCORE ) ) DEALLOCATE ( PIHCORE )
   IF ( ALLOCATED ( PIHDIAG ) ) DEALLOCATE ( PIHDIAG )
   IF ( ALLOCATED ( PISETEI ) ) DEALLOCATE ( PISETEI )
   IF ( ALLOCATED ( SETEI   ) ) DEALLOCATE ( SETEI   )
   ALLOCATE ( HCORE(1:NBASTR), OVERLAP(1:NBASTR-NPIHELE), PIHCORE(1:NPIHELE,1:NPIBEADS),       &
              PIOVERLAP(1:NPIHELE,1:NPIBEADS), PISETEI(1:NPI2ELEC,1:NPIBEADS), SETEI(1:N2ELEC) )

   ! . Allocate PIHDIAG.
   IF ( NPIBEADS > 1 ) THEN
      ALLOCATE ( HSTORE(1:NBASTR-NPIHELE), PIHDIAG(1:10,1:(NATOMSQM-(NPIATOMS*NPIBEADS)),1:NPIBEADS) )
   ELSE
      ALLOCATE ( HSTORE(1:1), PIHDIAG(1:1,1:1,1:1) )
   END IF

   !----------------------------------------------------------------------------
   ! . Start the calculation of HCORE.
   !----------------------------------------------------------------------------
   ! . Initialize some counters.
   NNOTRI = NBASTR - NPIHELE

   ! . Initialize the one-electron matrix.
   HCORE   = 0.0_DP
   PIHCORE = 0.0_DP
   PIHDIAG = 0.0_DP

   ! . Fill the diagonal elements of the one-electron matrix with the kinetic energy terms.
   I = 0
   DO IBASIS = 1,(NBASIS-NPIBASIS)
      I = I + IBASIS
      HCORE(I) = USPD(IBASIS)
   END DO
   I = I - NNOTRI
   DO IBASIS = (NBASIS-NPIBASIS+1),NBASIS
      I = I + IBASIS
      PIHCORE(I,1:NPIBEADS) = USPD(IBASIS)
   END DO

   ! . Calculate the overlap integrals using a Gaussian (STO-6G) approximation.
   CALL GAUSSIAN_OVERLAP ( OVERLAP, PIOVERLAP )

   ! . Multiply the atom-atom overlaps by (Bi+Bj)/2 but zero the same atom terms.
   ! . Note that the overlaps on the same atom are zero except for the diagonal
   ! . terms.
   I = 0
   DO IBASIS = 1,(NBASIS-NPIBASIS)
      ! . Calculate the off-diagonal terms.
      DO JBASIS = 1,(IBASIS-1)
         BETAIJ = 0.5_DP * ( BETA(IBASIS) + BETA(JBASIS) )
         I = I + 1
         OVERLAP(I) = BETAIJ * OVERLAP(I)
      END DO
      ! . Zero the diagonal terms.
      I = I + 1
      OVERLAP(I) = 0.0_DP
   END DO
   I = I - NNOTRI
   DO IBASIS = (NBASIS-NPIBASIS+1),NBASIS
      ! . Calculate the off-diagonal terms.
      DO JBASIS = 1,(IBASIS-1)
         BETAIJ = 0.5_DP * ( BETA(IBASIS) + BETA(JBASIS) )
         I = I + 1
         PIOVERLAP(I,1:NPIBEADS) = BETAIJ * PIOVERLAP(I,1:NPIBEADS)
      END DO
      ! . Zero the diagonal terms.
      I = I + 1
      PIOVERLAP(I,1:NPIBEADS) = 0.0_DP
   END DO

   ! . Add the overlap terms into the HCORE matrix.
     HCORE(1:NBASTR-NPIHELE) =   HCORE(1:NBASTR-NPIHELE) +   OVERLAP(:)
   PIHCORE                   = PIHCORE                   + PIOVERLAP

   ! . Deallocate the overlap array.
   DEALLOCATE ( OVERLAP, PIOVERLAP )

   !----------------------------------------------------------------------------
   ! . Calculate the TEI contributions.
   !----------------------------------------------------------------------------
   ! . Initialize the TEI counters.
   INDTEI   = 1
   PIINDTEI = 1

   ! . Initialize the nuclear energies.
   ENUCLR   = 0.0_DP
   PIENUCLR = 0.0_DP

   ! . Outer loop over QM atoms.
   DO IQM = 1,NATOMSQM

      ! . Get some information for the atom.
      ICOPY  = QMATOM(IQM)%COPY
      IFIRST = QMATOM(IQM)%BFIRST
      ILAST  = QMATOM(IQM)%BLAST
      NI     = QMATOM(IQM)%NUMBER
      RI     = MOPAC_DATA_ATOM_POSITION ( IQM, ATMCRD )

      ! . Inner loop over QM atoms.
      DO JQM = 1,(IQM-1)

         ! . Get some information for the atom.
         JCOPY  = QMATOM(JQM)%COPY
         JFIRST = QMATOM(JQM)%BFIRST
         JLAST  = QMATOM(JQM)%BLAST
         NJ     = QMATOM(JQM)%NUMBER
         RJ     = MOPAC_DATA_ATOM_POSITION ( JQM, ATMCRD )

         ! . Both atoms are non-PI.
         IF ( ( ICOPY == 0 ) .AND. ( JCOPY == 0 ) ) THEN

            ! . Calculate the two-electron integrals for the atom pair.
            CALL ROTATEI ( NI, NJ, RI, RJ, SETEI(INDTEI:), INDTEI, E1B, E2A, ENUC, QAM1PM3 )
!if ( abs ( enuc ) > 1000000.0 ) then
!write ( 6, "(4i10,f25.15)" ) qmatom(iqm)%atom, qmatom(jqm)%atom, iqm, jqm, ENUC
!end if
            ! . Sum the core-core energy term into the total.
            ENUCLR = ENUCLR + ENUC

            ! . Add in the electron-nuclear attraction terms for atom I.
            IF ( NATORB(NI) > 0 ) THEN
            J = 0
            DO IBASIS = IFIRST,ILAST
               I = (IBASIS*(IBASIS-1))/2 + IFIRST-1
               DO JBASIS = IFIRST,IBASIS
                  I = I + 1
                  J = J + 1
                  HCORE(I) = HCORE(I) + E1B(J)
               END DO
            END DO
            END IF

            ! . Add in the electron-nuclear attraction terms for atom J.
            IF ( NATORB(NJ) > 0 ) THEN
            J = 0
            DO IBASIS = JFIRST,JLAST
               I = (IBASIS*(IBASIS-1))/2 + JFIRST-1
               DO JBASIS = JFIRST,IBASIS
                  I = I + 1
                  J = J + 1
                  HCORE(I) = HCORE(I) + E2A(J)
               END DO
            END DO
            END IF

         ! . One atom is not-PI or both are PI and in the same polymer.
         ELSE IF ( ( ICOPY == 0 ) .OR. ( JCOPY == 0 ) .OR. ( ICOPY == JCOPY ) ) THEN

            ! . Get the appropriate block number.
            P = MAX ( ICOPY, JCOPY )

            ! . Calculate the two-electron integrals for the atom pair.
            CALL ROTATEI ( NI, NJ, RI, RJ, PISETEI(PIINDTEI(P):,P), PIINDTEI(P), E1B, E2A, ENUC, QAM1PM3 )

            ! . Sum the core-core energy term into the total.
            PIENUCLR(P) = PIENUCLR(P) + ENUC

            ! . Add in the electron-nuclear attraction terms for atom I.
            IF ( ICOPY == 0 ) THEN
               J = 0
               DO IBASIS = IFIRST,ILAST
                  DO JBASIS = IFIRST,IBASIS
                     J = J + 1
                     PIHDIAG(J,IQM,P) = PIHDIAG(J,IQM,P) + E1B(J)
                  END DO
               END DO
            ELSE
               J = 0
               DO IBASIS = IFIRST,ILAST
                  I = (IBASIS*(IBASIS-1))/2 + IFIRST-1 - NNOTRI
                  DO JBASIS = IFIRST,IBASIS
                     I = I + 1
                     J = J + 1
                     PIHCORE(I,P) = PIHCORE(I,P) + E1B(J)
                  END DO
               END DO
            END IF

            ! . Add in the electron-nuclear attraction terms for atom J.
            IF ( JCOPY == 0 ) THEN
               J = 0
               DO IBASIS = JFIRST,JLAST
                  DO JBASIS = JFIRST,IBASIS
                     J = J + 1
                     PIHDIAG(J,JQM,P) = PIHDIAG(J,JQM,P) + E2A(J)
                  END DO
               END DO
            ELSE
               J = 0
               DO IBASIS = JFIRST,JLAST
                  I = (IBASIS*(IBASIS-1))/2 + JFIRST-1 - NNOTRI
                  DO JBASIS = JFIRST,IBASIS
                     I = I + 1
                     J = J + 1
                     PIHCORE(I,P) = PIHCORE(I,P) + E2A(J)
                  END DO
               END DO
            END IF

         END IF
      END DO
   END DO

   !----------------------------------------------------------------------------
   ! . Calculate the one-electron QM/MM electrostatic interaction integrals.
   !----------------------------------------------------------------------------
   IF ( NATOMSMM > 0 ) THEN

      ! . Re-allocate the arrays for a derivative calculation.
      IF ( QGRADIENT ) THEN

         ! . Initialize the number of QM/MM interactions not involving a QM hydrogen.
	 NQMINTH = 0

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
	    I = NBCURRENT%ATOM

            ! . Loop over the interactions for the atom.
            DO IINT = 1,SIZE(NBCURRENT%INTERACTIONS)

               ! . Get the second atom of the interaction.
               J = NBCURRENT%INTERACTIONS(IINT)

               ! . Determine which is the QM atom. A negative J indicates J is quantum.
               IF ( J > 0 ) THEN
		  QINDX = ATMQMI(I)
               ELSE
		  QINDX = ATMQMI(-J)
               END IF

               ! . Increment the number of non-hydrogen interactions.
               IF ( NATORB(QMATOM(QINDX)%NUMBER) > 1 ) NQMINTH = NQMINTH + 1

	    END DO
	 END DO

         ! . Allocate the derivative integral arrays.
         ALLOCATE (   DXQM(1:NQMINT),         DYQM(1:NQMINT),         DZQM(1:NQMINT),                   &
	            DXE1BA(1:NQMINT),       DYE1BA(1:NQMINT),       DZE1BA(1:NQMINT),                   &
                    DXE1BB(2:10,1:NQMINTH), DYE1BB(2:10,1:NQMINTH), DZE1BB(2:10,1:NQMINTH), STAT = FLAG )
         IF ( FLAG /= 0 ) CALL PRINT_ERROR ( "INTEGRALS", "Unable to allocate MOPAC derivative integral arrays." )

      END IF

      ! . Calculate the interactions.
      CALL SEINTC

   END IF

   !----------------------------------------------------------------------------
   ! . Finish up.
   !----------------------------------------------------------------------------
   ! . Save the nuclear energies.
   ENUCLEAR = ENUCLR + PIENUCLR

   ! . Save the invariant part of the one-electron matrix.
   IF ( NPIBEADS > 1 ) HSTORE = HCORE(1:NBASTR-NPIHELE)

   !===============================================================================
   CONTAINS
   !===============================================================================

      !----------------
      SUBROUTINE SEINTC
      !----------------

      ! . Local variables.
      INTEGER            :: I, IATOM, II, IINT, IN, IREC, J, JATOM, JJ, N1, N2, MATOM, NINTQ, NQM, NTERM, P, QATOM, QINDX
      REAL ( KIND = DP ) :: ALPHA, CGQM, D1, D2, D4, F1, F2, F3, GAMMA, RHO0, RHO1, RHO2, R2OFF, R2ON, SWITCH, XQ, YQ, ZQ

      ! . Local integral scalars.
      REAL ( KIND = DP ) :: SS,   SPZ,   PZPZ,   PPPP,  &
                            DSS,  DSPZ,  DPZPZ,  DPPPP, &
                            SFCT1M,   SFCT1P,  SFCT20,  SFCT2M,  SFCT2P,  SFCT4P, &
                            DSFCT1M, DSFCT1P, DSFCT20, DSFCT2M, DSFCT2P, DSFCT4P

      ! . Local intermediate scalars.
      REAL ( KIND = DP ) :: CGL, DXENUC, DYENUC, DZENUC, DXPPPP, DYPPPP, DZPPPP, &
                            DXPZPZ, DYPZPZ, DZPZPZ, DXSPZ, DYSPZ, DZSPZ,         &
                            DXSS, DYSS, DZSS, ENUC, RQM, RQM2, RQMB, RQMI,       &
                            TEMP1, TEMP2, TEMP3, TEMP4, TEMP5,                   &
                            XN, YN, ZN, XN2, YN2, ZN2, XQM, YQM, ZQM

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3)  :: GVEC, MODGVEC, RI, RQ
      REAL ( KIND = DP ), DIMENSION(1:10) :: E1B

      ! . Local types.
      TYPE(NBLIST_TYPE), POINTER :: NBCURRENT, NBNEXT

      ! . QM/MM parameters.
      REAL ( KIND = DP ), PARAMETER :: ALPMM = 5.0_DP, RHO0MM = 0.0_DP

      ! . Calculate some factors.
      R2OFF = CUT_OFF**2
      R2ON  = CUT_ON**2
      GAMMA = ( R2OFF - R2ON )**3
      CALL SEINTC_SETUP

      ! . Initialization.
      IF ( QGRADIENT ) THEN
         DXE1BA = 0.0_DP ; DXE1BB = 0.0_DP ; DXQM = 0.0_DP
         DYE1BA = 0.0_DP ; DYE1BB = 0.0_DP ; DYQM = 0.0_DP
         DZE1BA = 0.0_DP ; DZE1BB = 0.0_DP ; DZQM = 0.0_DP
      END IF

      ! . Initialize the number of interactions counters.
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

         ! . Get the position vector for IATOM if it is a quantum or a boundary atom.
	 IF ( ATMQMI(IATOM) > 0 ) RI = MOPAC_DATA_ATOM_POSITION ( ATMQMI(IATOM), ATMCRD )

         ! . Loop over the interactions for the atom.
         DO IN = 1,SIZE(NBCURRENT%INTERACTIONS)

            ! . Get the second atom of the interaction.
            JATOM = NBCURRENT%INTERACTIONS(IN)

            ! . Increment the total interactions counter.
	    IINT = IINT + 1

            ! . Determine which is the QM atom. A negative J indicates J is quantum.
            IF ( JATOM > 0 ) THEN
               QATOM =   IATOM
               MATOM =   JATOM
	       QINDX =   ATMQMI(QATOM)
	       RQ    =   RI
            ELSE
               QATOM = - JATOM
               MATOM =   IATOM
	       QINDX =   ATMQMI(QATOM)
               RQ    =   MOPAC_DATA_ATOM_POSITION ( QINDX, ATMCRD )
            END IF

            ! . Define some constants for the QM atom.
            NQM   = QMATOM(QINDX)%NUMBER
            ALPHA = ALP(NQM)
            CGQM  = CORE(NQM)
            N1    = QMATOM(QINDX)%BFIRST
            N2    = QMATOM(QINDX)%BLAST
            P     = QMATOM(QINDX)%COPY
            XQ    = RQ(1)
            YQ    = RQ(2)
            ZQ    = RQ(3)
            RHO0  = ( 0.5_DP / AM(NQM) + RHO0MM )**2

            IF ( QAM1PM3 ) NTERM = COUNT ( ABS ( FN1(NQM,1:4) ) > 0.0_DP )

            ! . Get the MM atom charge.
#ifdef	HACK_LAMB
            CGL = ATMCHG(MATOM) * MOPAC_LAMBDA_VALUE
#else
            CGL = ATMCHG(MATOM)
#endif

            ! . Define the interatomic distance vector.
            XQM = ( XQ - ATMCRD(1,MATOM) )
            YQM = ( YQ - ATMCRD(2,MATOM) )
            ZQM = ( ZQ - ATMCRD(3,MATOM) )
            IF ( QIMAGE ) THEN
               XQM = XQM - BOXL(1) * ANINT ( XQM / BOXL(1), DP )
               YQM = YQM - BOXL(2) * ANINT ( YQM / BOXL(2), DP )
               ZQM = ZQM - BOXL(3) * ANINT ( ZQM / BOXL(3), DP )
            END IF
            RQM2 = XQM * XQM + YQM * YQM + ZQM * ZQM

            ! . Increment the interaction number.
            IF ( NATORB(NQM) > 1 ) NINTQ = NINTQ + 1

            ! . Check the interaction distance.
            IF ( RQM2 > R2OFF ) CYCLE

            ! . Calculate some more distance dependent terms.
            RQM  = SQRT ( RQM2 )
            RQMI = 1.0_DP / RQM
            XN   = RQMI * XQM
            YN   = RQMI * YQM
            ZN   = RQMI * ZQM

            ! . Calculate the switching function.
            IF ( RQM2 > R2ON ) THEN
               SWITCH = ( R2OFF - RQM2 )**2 * ( R2OFF + 2.0_DP * RQM2 - 3.0_DP * R2ON ) / GAMMA
            ELSE
               SWITCH = 1.0_DP
            END IF

            ! . QM atoms with one orbital.
            IF ( NATORB(NQM) == 1 ) THEN

               ! . Calculate some factors.
               RQMB = ANGSTROMS_TO_BOHRS * RQM

               ! . Calculate the SS integral and its derivative.
               CALL SEINTC_INTEGRALB ( RQMB, RHO0, SWITCH, SS, DSS )

               ! . Calculate the integrals.
               E1B(1)    = AU_TO_EV * SS
               E1B(2:10) = 0.0_DP

               ! . The derivatives.
               IF ( QGRADIENT ) THEN

                  ! . Scale DSS.
                  DSS = ANGSTROMS_TO_BOHRS * DSS

                  ! . The derivatives.
                  DXE1BA(IINT) = XN * DSS
                  DYE1BA(IINT) = YN * DSS
                  DZE1BA(IINT) = ZN * DSS

               END IF

            ! . QM atoms with more than one orbital.
            ELSE

               ! . Calculate some intermediate quantities.
               RHO1 = ( 0.5_DP / AD(NQM) + RHO0MM )**2
               RHO2 = ( 0.5_DP / AQ(NQM) + RHO0MM )**2
               D1   = DD(NQM)
               D2   = 2.0_DP * QQ(NQM)
               D4   = D2 * D2
               XN2  = XN * XN
               YN2  = YN * YN
               ZN2  = ZN * ZN

               ! . Calculate some factors.
               RQMB = ANGSTROMS_TO_BOHRS * RQM

               ! . Calculate the SS integral and its derivative.
               CALL SEINTC_INTEGRALB ( RQMB, RHO0, SWITCH, SS, DSS )

               ! . Calculate the SPZ integral and its derivative.
               CALL SEINTC_INTEGRALAB ( RQMB,  D1, RHO1, SWITCH, SFCT1P, DSFCT1P )
               CALL SEINTC_INTEGRALAB ( RQMB, -D1, RHO1, SWITCH, SFCT1M, DSFCT1M )
               SPZ  = 0.5_DP * (  SFCT1P -  SFCT1M )
               DSPZ = 0.5_DP * ( DSFCT1P - DSFCT1M )

               ! . Calculate the PZPZ integral and its derivative.
               CALL SEINTC_INTEGRALB  ( RQMB,      RHO2, SWITCH, SFCT20, DSFCT20 )
               CALL SEINTC_INTEGRALAB ( RQMB,  D2, RHO2, SWITCH, SFCT2P, DSFCT2P )
               CALL SEINTC_INTEGRALAB ( RQMB, -D2, RHO2, SWITCH, SFCT2M, DSFCT2M )
               PZPZ  =  SS + 0.25_DP * (  SFCT2P +  SFCT2M ) - 0.5_DP *  SFCT20
               DPZPZ = DSS + 0.25_DP * ( DSFCT2P + DSFCT2M ) - 0.5_DP * DSFCT20

               ! . Calculate the PPPP integral and its derivative.
               CALL SEINTC_INTEGRALB ( RQMB, ( D4 + RHO2 ), SWITCH, SFCT4P, DSFCT4P )
               PPPP  =  SS + 0.5_DP * (  SFCT4P -  SFCT20 )
               DPPPP = DSS + 0.5_DP * ( DSFCT4P - DSFCT20 )

               ! . Fill arrays with the calculated integrals.
               E1B(1)  = SS
               E1B(2)  = XN * SPZ
               E1B(3)  = XN2 * PZPZ + ( YN2 + ZN2 ) * PPPP
               E1B(4)  = YN * SPZ
               E1B(5)  = XN * YN * ( PZPZ - PPPP )
               E1B(6)  = YN2 * PZPZ + ( XN2 + ZN2 ) * PPPP
               E1B(7)  = ZN * SPZ
               E1B(8)  = XN * ZN * ( PZPZ - PPPP )
               E1B(9)  = YN * ZN * ( PZPZ - PPPP )
               E1B(10) = ZN2 * PZPZ + ( XN2 + YN2 ) * PPPP

               ! . Scale the integrals.
               E1B(1)    = AU_TO_EV * E1B(1)
               E1B(2:10) = AU_TO_EV * CGL * E1B(2:10)

               ! . Calculate some intermediate arrays for the derivatives.
               IF ( QGRADIENT ) THEN

                  ! . Scale the derivative integrals.
                  DSS   = ANGSTROMS_TO_BOHRS * DSS
                  DSPZ  = ANGSTROMS_TO_BOHRS * DSPZ
                  DPZPZ = ANGSTROMS_TO_BOHRS * DPZPZ
                  DPPPP = ANGSTROMS_TO_BOHRS * DPPPP

                  ! . Calculate the vector components of the derivatives.
                  DXSS  = XN * DSS
                  DYSS  = YN * DSS
                  DZSS  = ZN * DSS

                  DXSPZ = XN * DSPZ
                  DYSPZ = YN * DSPZ
                  DZSPZ = ZN * DSPZ

                  DXPZPZ = XN * DPZPZ
                  DYPZPZ = YN * DPZPZ
                  DZPZPZ = ZN * DPZPZ

                  DXPPPP = XN * DPPPP
                  DYPPPP = YN * DPPPP
                  DZPPPP = ZN * DPPPP

                  ! . Calculate the derivatives.
                  TEMP1 = PZPZ - PPPP
                  TEMP2 = YN2 + ZN2
                  TEMP3 = XN2 + ZN2
                  TEMP4 = RQMI * SPZ
                  TEMP5 = XN2 + YN2

                  DXE1BA(IINT) = DXSS
                  DYE1BA(IINT) = DYSS
                  DZE1BA(IINT) = DZSS

                  DXE1BB(2,NINTQ) = TEMP2 * TEMP4 + XN * DXSPZ
                  DYE1BB(2,NINTQ) = -XN * (TEMP4 * YN - DYSPZ)
                  DZE1BB(2,NINTQ) = -XN * (TEMP4 * ZN - DZSPZ)

                  DXE1BB(3,NINTQ) = XN2 * DXPZPZ + TEMP2 * DXPPPP + 2.0_DP * XN * RQMI * TEMP2 * TEMP1
                  DYE1BB(3,NINTQ) = XN2 * DYPZPZ + TEMP2 * DYPPPP - 2.0_DP * YN * RQMI * XN2 * TEMP1
                  DZE1BB(3,NINTQ) = XN2 * DZPZPZ + TEMP2 * DZPPPP - 2.0_DP * ZN * RQMI * XN2 * TEMP1

                  DXE1BB(4,NINTQ) = -YN * (XN * TEMP4 - DXSPZ)
                  DYE1BB(4,NINTQ) = TEMP4 * TEMP3 + YN * DYSPZ
                  DZE1BB(4,NINTQ) = -YN * (ZN * TEMP4 - DZSPZ)

                  DXE1BB(5,NINTQ) = YN * (RQMI * TEMP1 * (TEMP2 - XN2) + XN * (DXPZPZ - DXPPPP))
                  DYE1BB(5,NINTQ) = XN * (RQMI * TEMP1 * (TEMP3 - YN2) + YN * (DYPZPZ - DYPPPP))
                  DZE1BB(5,NINTQ) = XN * YN * ((DZPZPZ - DZPPPP) - 2.0_DP * ZN * RQMI * TEMP1)

                  DXE1BB(6,NINTQ) = YN2 * DXPZPZ + TEMP3 * DXPPPP - 2.0_DP * XN * YN2 * RQMI * TEMP1
                  DYE1BB(6,NINTQ) = YN2 * DYPZPZ + TEMP3 * (DYPPPP + 2.0_DP * YN * RQMI * TEMP1)
                  DZE1BB(6,NINTQ) = YN2 * DZPZPZ + TEMP3 * DZPPPP - 2.0_DP * ZN * YN2 * RQMI * TEMP1

                  DXE1BB(7,NINTQ) = -ZN * (XN * TEMP4 - DXSPZ)
                  DYE1BB(7,NINTQ) = -ZN * (YN * TEMP4 - DYSPZ)
                  DZE1BB(7,NINTQ) = TEMP5 * TEMP4 + ZN * DZSPZ

                  DXE1BB(8,NINTQ) = ZN * (RQMI * TEMP1 * (TEMP2 - XN2) + XN * (DXPZPZ - DXPPPP))
                  DYE1BB(8,NINTQ) = XN * ZN * ((DYPZPZ - DYPPPP) - 2.0_DP * YN * RQMI * TEMP1)
                  DZE1BB(8,NINTQ) = XN * (RQMI * TEMP1 * (TEMP5 - ZN2) + ZN * (DZPZPZ - DZPPPP))

                  DXE1BB(9,NINTQ) = YN * ZN * ((DXPZPZ - DXPPPP) - 2.0_DP * XN * RQMI * TEMP1)
                  DYE1BB(9,NINTQ) = ZN * (RQMI * TEMP1 * (TEMP3 - YN2) + YN * (DYPZPZ - DYPPPP))
                  DZE1BB(9,NINTQ) = YN * (RQMI * TEMP1 * (TEMP5 - ZN2) + ZN * (DZPZPZ - DZPPPP))

                  DXE1BB(10,NINTQ) = ZN2 * DXPZPZ + TEMP5 * DXPPPP - 2.0_DP * XN * ZN2 * RQMI * TEMP1
                  DYE1BB(10,NINTQ) = ZN2 * DYPZPZ + TEMP5 * DYPPPP - 2.0_DP * YN * ZN2 * RQMI * TEMP1
                  DZE1BB(10,NINTQ) = ZN2 * DZPZPZ + TEMP5 * (DZPPPP + 2.0_DP * ZN * RQMI * TEMP1)

                  ! . Multiply the derivative integrals by CGL.
                  DXE1BB(2:10,NINTQ) = CGL * DXE1BB(2:10,NINTQ)
                  DYE1BB(2:10,NINTQ) = CGL * DYE1BB(2:10,NINTQ)
                  DZE1BB(2:10,NINTQ) = CGL * DZE1BB(2:10,NINTQ)
               END IF
            END IF

            ! . Calculate the nuclear/nuclear terms.
            TEMP1 = E1B(1)
            TEMP2 = CGQM * CGL
            TEMP3 = ABS ( TEMP2 ) * EXP ( - ALPHA * RQM )
            TEMP4 = ABS ( TEMP2 ) * EXP ( - ALPMM * RQM )
            ENUC  = TEMP1 * ( TEMP2 + TEMP3 + TEMP4 )

            IF ( QGRADIENT ) THEN
               DXENUC = AU_TO_EV * ( TEMP2 + TEMP3 + TEMP4 ) * DXE1BA(IINT) - XN * TEMP1 * ( ALPHA * TEMP3 + ALPMM * TEMP4 )
               DYENUC = AU_TO_EV * ( TEMP2 + TEMP3 + TEMP4 ) * DYE1BA(IINT) - YN * TEMP1 * ( ALPHA * TEMP3 + ALPMM * TEMP4 )
               DZENUC = AU_TO_EV * ( TEMP2 + TEMP3 + TEMP4 ) * DZE1BA(IINT) - ZN * TEMP1 * ( ALPHA * TEMP3 + ALPMM * TEMP4 )
            END IF

            IF ( QAM1PM3 ) THEN
               DO I = 1,NTERM
                  F1    = FN1(NQM,I)
                  F2    = FN2(NQM,I)
                  F3    = FN3(NQM,I)
                  TEMP1 = EXP ( MAX ( -30.0_DP, - F2 * ( RQM - F3 )**2 ) )
                  ENUC  = ENUC + CGQM * CGL * RQMI * F1 * TEMP1
                  IF ( QGRADIENT ) THEN
                     TEMP2 = CGQM * CGL / RQM2
                     DXENUC = DXENUC - F1 * TEMP1 * TEMP2 * XN * ( 1.0_DP + 2.0_DP * F2 * RQM * ( RQM - F3 ) )
                     DYENUC = DYENUC - F1 * TEMP1 * TEMP2 * YN * ( 1.0_DP + 2.0_DP * F2 * RQM * ( RQM - F3 ) )
                     DZENUC = DZENUC - F1 * TEMP1 * TEMP2 * ZN * ( 1.0_DP + 2.0_DP * F2 * RQM * ( RQM - F3 ) )
                  END IF
               END DO
            END IF
!if ( abs ( enuc ) > 1000000.0 ) then
!write ( 6, "(2I10,2F25.15)" ) IATOM, JATOM, RQM, ENUC
!end if
            ! . Scale the E1B(1) integral by CGL now ENUC is finished.
            E1B(1) = CGL * E1B(1)
            IF ( QGRADIENT ) THEN
               DXE1BA(IINT) = CGL * DXE1BA(IINT)
               DYE1BA(IINT) = CGL * DYE1BA(IINT)
               DZE1BA(IINT) = CGL * DZE1BA(IINT)
            END IF

            ! . Accumulate ENUCLR.
            IF ( P == 0 ) THEN
               ENUCLR = ENUCLR + ENUC
            ELSE
               PIENUCLR(P) = PIENUCLR(P) + ENUC
            END IF

            ! . Add the integrals into the one-electron matrix.
            IF ( P == 0 ) THEN
               JJ = 0
               DO I = N1,N2
                  II = (I * (I - 1)) / 2 + N1 - 1
                  DO J = N1,I
                     II = II + 1
                     JJ = JJ + 1
                     HCORE(II) = HCORE(II) - E1B(JJ)
                  END DO
               END DO
            ELSE
               JJ = 0
               DO I = N1,N2
                  II = (I * (I - 1)) / 2 + N1 - 1 - NNOTRI
                  DO J = N1,I
                     II = II + 1
                     JJ = JJ + 1
                     PIHCORE(II,P) = PIHCORE(II,P) - E1B(JJ)
                  END DO
               END DO
            END IF

            ! . Put the nuclear/nuclear terms into the derivative arrays.
            IF ( QGRADIENT ) THEN

               ! . Calculate some scale factors.
               TEMP1 = EV_TO_KJ
               IF ( P > 0 ) TEMP1 = TEMP1 / REAL ( NPIBEADS, DP )

               ! . Determine the gradient vector.
	       GVEC(1) = TEMP1 * DXENUC
	       GVEC(2) = TEMP1 * DYENUC
	       GVEC(3) = TEMP1 * DZENUC

               ! . Add in the contribution to the MM atom gradients.
               GRADIENT(1:3,MATOM) = GRADIENT(1:3,MATOM) - GVEC

               ! . Add in the contribution to the QM atom gradients.
	       IF ( QMATOM(QINDX)%QBOUNDARY ) THEN

                  ! . Calculate the gradients.
	          MODGVEC = MOPAC_DATA_LA_GRADIENT ( QINDX, ATMCRD, GVEC )
                  GRADIENT(1:3,QATOM)                 = GRADIENT(1:3,QATOM)                        + MODGVEC
                  GRADIENT(1:3,QMATOM(QINDX)%PARTNER) = GRADIENT(1:3,QMATOM(QINDX)%PARTNER) + GVEC - MODGVEC

                  ! . Shift the interatomic vector to the partner atom.
		  XQM = XQM + ( ATMCRD(1,QMATOM(QINDX)%PARTNER) - XQ )
		  YQM = YQM + ( ATMCRD(2,QMATOM(QINDX)%PARTNER) - YQ )
		  ZQM = ZQM + ( ATMCRD(3,QMATOM(QINDX)%PARTNER) - ZQ )

	       ELSE
                  GRADIENT(1:3,QATOM) = GRADIENT(1:3,QATOM) + GVEC
	       END IF

               ! . Determine the contribution to the virial.
	       VIRIAL = VIRIAL + GVEC(1) * XQM + GVEC(2) * YQM + GVEC(3) * ZQM

               ! . Save the displacement vector.
	       DXQM(IINT) = XQM
	       DYQM(IINT) = YQM
	       DZQM(IINT) = ZQM

            END IF
         END DO
      END DO

      ! . Scale the derivatives correctly.
      IF ( QGRADIENT ) THEN
         DXE1BA = AU_TO_EV * EV_TO_KJ * DXE1BA
         DYE1BA = AU_TO_EV * EV_TO_KJ * DYE1BA
         DZE1BA = AU_TO_EV * EV_TO_KJ * DZE1BA
         DXE1BB = AU_TO_EV * EV_TO_KJ * DXE1BB
         DYE1BB = AU_TO_EV * EV_TO_KJ * DYE1BB
         DZE1BB = AU_TO_EV * EV_TO_KJ * DZE1BB
      END IF

      END SUBROUTINE SEINTC

   END SUBROUTINE INTEGRALS

   !-----------------------------------
   SUBROUTINE INTEGRALS_HCORE ( ICALC )
   !-----------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: ICALC

   ! . Local scalars.
   INTEGER :: I, IBASIS, IFIRST, ILAST, IQM, J, JBASIS

   ! . There are no PI atoms (HCORE is complete).
   IF ( NPIBASIS <= 0 ) RETURN

   ! . Copy over HSTORE.
   HCORE(1:NBASTR-NPIHELE) = HSTORE

   ! . Fill HCORE with the appropriate missing PI atom elements.
   HCORE(NBASTR-NPIHELE+1:NBASTR) = PIHCORE(1:NPIHELE,ICALC)

   ! . Complement the diagonal elements of the non-PI atoms.
   DO IQM = 1,(NATOMSQM-(NPIATOMS*NPIBEADS))

      ! . Get some information for the atom.
      IFIRST = QMATOM(IQM)%BFIRST
      ILAST  = QMATOM(IQM)%BLAST

      ! . Loop over the elements.
      J = 0
      DO IBASIS = IFIRST,ILAST
         I = (IBASIS*(IBASIS-1))/2 + IFIRST-1
         DO JBASIS = IFIRST,IBASIS
            I = I + 1
            J = J + 1
            HCORE(I) = HCORE(I) + PIHDIAG(J,IQM,ICALC)
         END DO
      END DO
   END DO

   END SUBROUTINE INTEGRALS_HCORE

   !---------------------------------------------------
   SUBROUTINE MOPAC_ESP ( METHOD, POINTS, DENMAT, PHI )
   !---------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: METHOD

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(IN)  :: DENMAT
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)  :: POINTS
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT) :: PHI

   ! . Local scalars.
   INTEGER :: NPOINTS
   LOGICAL :: QAM1PM3

   ! . Local allocatable arrays.
!   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: ONEMAT

   ! . Set the AM1/PM3 flag.
   QAM1PM3 = ( HAMILTONIAN == "AM1 " ) .OR. ( HAMILTONIAN == "PDDG " ) .OR. ( HAMILTONIAN == "PM3 " ) .OR. ( HAMILTONIAN == "RM1 " )

   ! . Allocate space.
!   ALLOCATE ( ONEMAT(1:NBASTR) )

   ! . Initialization.
   NPOINTS = SIZE ( POINTS, 2 )
   PHI     = 0.0_DP

   ! . Do the calculation.
   IF ( METHOD == "MNDO" ) THEN
      CALL MNDOESP
   ELSE
      CALL PRINT_ERROR ( "MOPAC_ESP", "Unknown ESP method." )
   END IF

   ! . Finish up.
!   DEALLOCATE ( ONEMAT )

   !===============================================================================
   CONTAINS
   !===============================================================================

      !-----------------
      SUBROUTINE MNDOESP
      !-----------------

      ! . Local variables.
      INTEGER            :: I, II, J, JJ, MATOM, NQM, N1, N2, NTERM, QATOM, QINDX
      REAL ( KIND = DP ) :: ALPHA, CGQM, D1, D2, D4, FACT, F1, F2, F3, PHILOCAL, RHO0, RHO1, RHO2, XQ, YQ, ZQ

      ! . Local integral scalars.
      REAL ( KIND = DP ) :: SS, SPZ, PZPZ, PPPP, DSS, DSPZ, DPZPZ, DPPPP, &
                            SFCT1M,   SFCT1P,  SFCT20,  SFCT2M,  SFCT2P,  SFCT4P, &
                            DSFCT1M, DSFCT1P, DSFCT20, DSFCT2M, DSFCT2P, DSFCT4P

      ! . Local intermediate scalars.
      REAL ( KIND = DP ) :: ENUC, RQM, RQM2, RQMB, RQMI, TEMP1, TEMP2, TEMP3, TEMP4, XN, YN, ZN, XN2, YN2, ZN2, XQM, YQM, ZQM

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3)  :: RQ
      REAL ( KIND = DP ), DIMENSION(1:10) :: E1B

      ! . QM/MM parameters.
      REAL ( KIND = DP ), PARAMETER :: ALPMM = 5.0_DP, RHO0MM = 0.0_DP

      ! . Loop over the points.
      DO MATOM = 1,NPOINTS

         ! . Initialization.
         PHILOCAL = 0.0_DP

         ! . Loop over QM atoms.
         DO QINDX = 1,NATOMSQM

            ! . Get some information for the atom.
            QATOM = QMATOM(QINDX)%ATOM
            RQ    = MOPAC_DATA_ATOM_POSITION ( QINDX, ATMCRD )

            ! . Define some constants for the QM atom.
            NQM   = QMATOM(QINDX)%NUMBER
            ALPHA = ALP(NQM)
            CGQM  = CORE(NQM)
            N1    = QMATOM(QINDX)%BFIRST
            N2    = QMATOM(QINDX)%BLAST
            XQ    = RQ(1)
            YQ    = RQ(2)
            ZQ    = RQ(3)
            RHO0  = ( 0.5_DP / AM(NQM) + RHO0MM )**2

            IF ( QAM1PM3 ) NTERM = COUNT ( ABS ( FN1(NQM,1:4) ) > 0.0_DP )

            ! . Define the interatomic distance vector.
            XQM = ( XQ - POINTS(1,MATOM) )
            YQM = ( YQ - POINTS(2,MATOM) )
            ZQM = ( ZQ - POINTS(3,MATOM) )
            RQM2 = XQM * XQM + YQM * YQM + ZQM * ZQM

            ! . Calculate some more distance dependent terms.
            RQM  = SQRT ( RQM2 )
            RQMI = 1.0_DP / RQM
            XN   = RQMI * XQM
            YN   = RQMI * YQM
            ZN   = RQMI * ZQM

            ! . QM atoms with one orbital.
            IF ( NATORB(NQM) == 1 ) THEN

               ! . Calculate some factors.
               RQMB = ANGSTROMS_TO_BOHRS * RQM

               ! . Calculate the SS integral and its derivative.
               CALL MNDOESP_INTEGRALB ( RQMB, RHO0, SS, DSS )

               ! . Calculate the integrals.
               E1B(1)    = AU_TO_EV * SS
               E1B(2:10) = 0.0_DP

            ! . QM atoms with more than one orbital.
            ELSE

               ! . Calculate some intermediate quantities.
               RHO1 = ( 0.5_DP / AD(NQM) + RHO0MM )**2
               RHO2 = ( 0.5_DP / AQ(NQM) + RHO0MM )**2
               D1   = DD(NQM)
               D2   = 2.0_DP * QQ(NQM)
               D4   = D2 * D2
               XN2  = XN * XN
               YN2  = YN * YN
               ZN2  = ZN * ZN

               ! . Calculate some factors.
               RQMB = ANGSTROMS_TO_BOHRS * RQM

               ! . Calculate the SS integral and its derivative.
               CALL MNDOESP_INTEGRALB ( RQMB, RHO0, SS, DSS )

               ! . Calculate the SPZ integral and its derivative.
               CALL MNDOESP_INTEGRALAB ( RQMB,  D1, RHO1, SFCT1P, DSFCT1P )
               CALL MNDOESP_INTEGRALAB ( RQMB, -D1, RHO1, SFCT1M, DSFCT1M )
               SPZ  = 0.5_DP * (  SFCT1P -  SFCT1M )
               DSPZ = 0.5_DP * ( DSFCT1P - DSFCT1M )

               ! . Calculate the PZPZ integral and its derivative.
               CALL MNDOESP_INTEGRALB  ( RQMB,      RHO2, SFCT20, DSFCT20 )
               CALL MNDOESP_INTEGRALAB ( RQMB,  D2, RHO2, SFCT2P, DSFCT2P )
               CALL MNDOESP_INTEGRALAB ( RQMB, -D2, RHO2, SFCT2M, DSFCT2M )
               PZPZ  =  SS + 0.25_DP * (  SFCT2P +  SFCT2M ) - 0.5_DP *  SFCT20
               DPZPZ = DSS + 0.25_DP * ( DSFCT2P + DSFCT2M ) - 0.5_DP * DSFCT20

               ! . Calculate the PPPP integral and its derivative.
               CALL MNDOESP_INTEGRALB ( RQMB, ( D4 + RHO2 ), SFCT4P, DSFCT4P )
               PPPP  =  SS + 0.5_DP * (  SFCT4P -  SFCT20 )
               DPPPP = DSS + 0.5_DP * ( DSFCT4P - DSFCT20 )

               ! . Fill arrays with the calculated integrals.
               E1B(1)  = SS
               E1B(2)  = XN * SPZ
               E1B(3)  = XN2 * PZPZ + ( YN2 + ZN2 ) * PPPP
               E1B(4)  = YN * SPZ
               E1B(5)  = XN * YN * ( PZPZ - PPPP )
               E1B(6)  = YN2 * PZPZ + ( XN2 + ZN2 ) * PPPP
               E1B(7)  = ZN * SPZ
               E1B(8)  = XN * ZN * ( PZPZ - PPPP )
               E1B(9)  = YN * ZN * ( PZPZ - PPPP )
               E1B(10) = ZN2 * PZPZ + ( XN2 + YN2 ) * PPPP

               ! . Scale the integrals.
               E1B(1)    = AU_TO_EV * E1B(1)
               E1B(2:10) = AU_TO_EV * E1B(2:10)

            END IF

            ! . Calculate the nuclear/nuclear terms.
            TEMP1 = E1B(1)
            TEMP2 = CGQM
            TEMP3 = ABS ( TEMP2 ) * EXP ( - ALPHA * RQM )
            TEMP4 = ABS ( TEMP2 ) * EXP ( - ALPMM * RQM )
            ENUC  = TEMP1 * ( TEMP2 + TEMP3 + TEMP4 )
            IF ( QAM1PM3 ) THEN
               DO I = 1,NTERM
                  F1    = FN1(NQM,I)
                  F2    = FN2(NQM,I)
                  F3    = FN3(NQM,I)
                  TEMP1 = EXP ( MAX ( -30.0_DP, - F2 * ( RQM - F3 )**2 ) )
                  ENUC  = ENUC + CGQM * RQMI * F1 * TEMP1
               END DO
            END IF

            ! . Accumulate phi.
            PHILOCAL = PHILOCAL + ENUC

            ! . Add the integrals into the one-electron matrix.
            JJ = 0
            DO I = N1,N2
               II = (I * (I - 1)) / 2 + N1 - 1
               DO J = N1,I
                  II = II + 1
                  JJ = JJ + 1
                  IF ( I == J ) THEN
                      FACT = 1.0_DP
                  ELSE
                      FACT = 2.0_DP
                  END IF
                  PHILOCAL = PHILOCAL - FACT * DENMAT(II) * E1B(JJ)
               END DO
            END DO

         ! . QM atoms.
         END DO

         ! . Set the ESP.
         PHI(MATOM) = PHILOCAL 

      ! . Points.
      END DO

      END SUBROUTINE MNDOESP

      !-----------------------------------------------
      SUBROUTINE MNDOESP_INTEGRALAB ( R, A, B, F, DF )
      !-----------------------------------------------

      ! . Scalar arguments.
      REAL ( KIND = DP ), INTENT(IN)  :: R, A, B
      REAL ( KIND = DP ), INTENT(OUT) :: F, DF

      ! . Calculate F and DF.
      F  = 1.0_DP / SQRT ( ( R + A )**2 + B )
      DF = - ( R + A ) * F * F * F

      END SUBROUTINE MNDOESP_INTEGRALAB

      !-------------------------------------------
      SUBROUTINE MNDOESP_INTEGRALB ( R, B, F, DF )
      !-------------------------------------------

      ! . Scalar arguments.
      REAL ( KIND = DP ), INTENT(IN)  :: R, B
      REAL ( KIND = DP ), INTENT(OUT) :: F, DF

      ! . Calculate F and DF.
      F  = 1.0_DP / SQRT ( R**2 + B )
      DF = - R * F * F * F

      END SUBROUTINE MNDOESP_INTEGRALB

   END SUBROUTINE MOPAC_ESP

   !------------------------------------------
   SUBROUTINE REPP ( NI, NJ, RIJ, RI, CORINT )
   !------------------------------------------

   !     REPP calculates the two-electron integrals and electron-nuclear
   !     attraction terms in the local frame between the atoms I and J.
   !
   !     Input:   NI, NJ        the atomic numbers.
   !              RIJ           interatomic distance in atomic units.
   !              AM, AD, AQ    arrays of two-electron one-centre integrals.
   !              CORE          array on atomic nuclear charges.
   !              DD            array of dipole charge separations.
   !              QQ            array of quadrupole charge separations.
   !
   !    Output:   RI            two-electron repulsion integrals.
   !              CORINT          array of electron-nuclear attraction integrals.
   !

   ! . Scalar argument declarations.
   INTEGER,            INTENT(IN) :: NI, NJ
   REAL ( KIND = DP ), INTENT(IN) :: RIJ

   ! . Array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:4,1:2), INTENT(OUT) :: CORINT
   REAL ( KIND = DP ), DIMENSION(1:22),    INTENT(OUT) :: RI

   ! . Local scalars.
   REAL ( KIND = DP ) :: ADD, ADE, ADQ, AED, AEE, AEQ, AQD, AQE, AQQ, DA, DB,  &
                         DXDX, DXQXZ, DZDZ, DZE, DZQXX, DZQZZ, EDZ, EE, EQXX,  &
                         EQZZ, QA, QB, QXXDZ, QXXE, QXXQXX, QXXQYY, QXXQZZ,    &
                         QXYQXY, QXZQXZ, QXZDX, QZZDZ, QZZE, QZZQXX, QZZQZZ, R

   ! . Local parameter declarations.
   REAL ( KIND = DP ), PARAMETER :: EIGHT = 8.0_DP, FD = 4.0_DP, FOUR = 4.0_DP, OD = 1.0_DP, &
                                    PT25  = 0.25_DP, SIXTN = 16.0_DP, TD = 2.0_DP

   ! . Initialize some variables.
   R      = RIJ
   CORINT = 0.0_DP
   RI     = 0.0_DP

   ! . Define the charge separations.
   DA = DD(NI)
   DB = DD(NJ)
   QA = QQ(NI)
   QB = QQ(NJ)

   ! . Hydrogen - Hydrogen terms.
   AEE = 0.25_DP * ( OD/AM(NI) + OD/AM(NJ) )**2
   EE  = OD / SQRT ( R**2+AEE )
   RI(1) = EE * AU_TO_EV
   CORINT(1,1) = CORE(NJ) * RI(1)
   CORINT(1,2) = CORE(NI) * RI(1)
   IF ( (NATORB(NI) < 3) .AND. (NATORB(NJ) < 3) ) RETURN
   IF ( (NATORB(NI) < 3) ) GOTO 10

   ! . Heavy atom - Hydrogen terms.
   ADE=0.25_DP*(OD/AD(NI)+OD/AM(NJ))**2
   AQE=0.25_DP*(OD/AQ(NI)+OD/AM(NJ))**2
   DZE=-OD/SQRT((R+DA)**2+ADE)+OD/SQRT((R-DA)**2+ADE)
   QZZE=OD/SQRT((R-TD*QA)**2+AQE)-TD/SQRT(R**2+AQE)+OD/SQRT((R+TD*QA)**2+AQE)
   QXXE=TD/SQRT(R**2+FD*QA**2+AQE)-TD/SQRT(R**2+AQE)
   DZE=DZE/2.0_DP
   QXXE=QXXE/FOUR
   QZZE=QZZE/FOUR
   RI(2)=-DZE
   RI(3)=EE+QZZE
   RI(4)=EE+QXXE
   IF ( NATORB(NJ) < 3 ) GOTO 20
   ! . Hydrogen - Heavy atom terms.
   10 CONTINUE
   AED=0.25_DP*(OD/AM(NI)+OD/AD(NJ))**2
   AEQ=0.25_DP*(OD/AM(NI)+OD/AQ(NJ))**2
   EDZ=-OD/SQRT((R-DB)**2+AED)+OD/SQRT((R+DB)**2+AED)
   EQZZ=OD/SQRT((R-TD*QB)**2+AEQ)-TD/SQRT(R**2+AEQ)+OD/SQRT((R+TD*QB)**2+AEQ)
   EQXX=TD/SQRT(R**2+FD*QB**2+AEQ)-TD/SQRT(R**2+AEQ)
   EDZ=EDZ/2.0_DP
   EQXX=EQXX/FOUR
   EQZZ=EQZZ/FOUR
   RI(5)=-EDZ
   RI(11)=EE+EQZZ
   RI(12)=EE+EQXX
   IF ( NATORB(NI) < 3 ) GOTO 20
   ! . Heavy atom - Heavy atom terms.
   ADD=0.25_DP*(OD/AD(NI)+OD/AD(NJ))**2
   ADQ=0.25_DP*(OD/AD(NI)+OD/AQ(NJ))**2
   AQD=0.25_DP*(OD/AQ(NI)+OD/AD(NJ))**2
   AQQ=0.25_DP*(OD/AQ(NI)+OD/AQ(NJ))**2
   DXDX=TD/SQRT(R**2+(DA-DB)**2+ADD)-TD/SQRT(R**2+(DA+DB)**2+ADD)
   DZDZ=OD/SQRT((R+DA-DB)**2+ADD)+OD/SQRT((R-DA+DB)**2+ADD)-OD/SQRT((R-DA-DB)**2+ADD)-OD/SQRT((R+DA+DB)**2+ADD)
   DZQXX=-TD/SQRT((R+DA)**2+FD*QB**2+ADQ)+TD/SQRT((R-DA)**2+FD*QB**2+ADQ)+TD/SQRT((R+DA)**2+ADQ)-TD/SQRT((R-DA)**2+ADQ)
   QXXDZ=-TD/SQRT((R-DB)**2+FD*QA**2+AQD)+TD/SQRT((R+DB)**2+FD*QA**2+AQD)+TD/SQRT((R-DB)**2+AQD)-TD/SQRT((R+DB)**2+AQD)
   DZQZZ=-OD/SQRT((R+DA-TD*QB)**2+ADQ)+OD/SQRT((R-DA-TD*QB)**2+ADQ)-OD/SQRT((R+DA+TD*QB)**2+ADQ)+OD/SQRT((R-DA+TD*QB) &
          **2+ADQ)-TD/SQRT((R-DA)**2+ADQ)+TD/SQRT((R+DA)**2+ADQ)
   QZZDZ=-OD/SQRT((R+TD*QA-DB)**2+AQD)+OD/SQRT((R+TD*QA+DB)**2+AQD)-OD/SQRT((R-TD*QA-DB)**2+AQD)+OD/SQRT((R-2.0_DP &
          *QA+DB)**2+AQD)+TD/SQRT((R-DB)**2+AQD)-TD/SQRT((R+DB)**2+AQD)
   QXXQXX=TD/SQRT(R**2+FD*(QA-QB)**2+AQQ)+TD/SQRT(R**2+FD*(QA+QB)**2+ &
          AQQ)-FD/SQRT(R**2+FD*QA**2+AQQ)-FD/SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT(R**2+AQQ)
   QXXQYY=FD/SQRT(R**2+FD*QA**2+FD*QB**2+AQQ)-FD/SQRT(R**2+FD*QA**2+AQQ)-FD/SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT(R**2+AQQ)
   QXXQZZ=TD/SQRT((R-TD*QB)**2+FD*QA**2+AQQ)+TD/SQRT((R+TD*QB)**2+FD* &
          QA**2+AQQ)-TD/SQRT((R-TD*QB)**2+AQQ)-TD/SQRT((R+TD*QB)**2+AQQ)-FD/ &
          SQRT(R**2+FD*QA**2+AQQ)+FD/SQRT(R**2+AQQ)
   QZZQXX=TD/SQRT((R+TD*QA)**2+FD*QB**2+AQQ)+TD/SQRT((R-TD*QA)**2+FD* &
          QB**2+AQQ)-TD/SQRT((R+TD*QA)**2+AQQ)-TD/SQRT((R-TD*QA)**2+AQQ)-FD/ &
          SQRT(R**2+FD*QB**2+AQQ)+FD/SQRT(R**2+AQQ)
   QZZQZZ=OD/SQRT((R+TD*QA-TD*QB)**2+AQQ)+OD/SQRT((R+TD*QA+TD*QB)**2+ &
          AQQ)+OD/SQRT((R-TD*QA-TD*QB)**2+AQQ)+OD/SQRT((R-TD*QA+TD*QB)**2+AQQ) &
          -TD/SQRT((R-TD*QA)**2+AQQ)-TD/SQRT((R+TD*QA)**2+AQQ)-TD/SQRT((R- &
          TD*QB)**2+AQQ)-TD/SQRT((R+TD*QB)**2+AQQ)+FD/SQRT(R**2+AQQ)
   DXQXZ=-TD/SQRT((R-QB)**2+(DA-QB)**2+ADQ)+TD/SQRT((R+QB)**2+(DA-QB) &
         **2+ADQ)+TD/SQRT((R-QB)**2+(DA+QB)**2+ADQ)-TD/SQRT((R+QB)**2+(DA+QB)**2+ADQ)
   QXZDX=-TD/SQRT((R+QA)**2+(QA-DB)**2+AQD)+TD/SQRT((R-QA)**2+(QA-DB) &
         **2+AQD)+TD/SQRT((R+QA)**2+(QA+DB)**2+AQD)-TD/SQRT((R-QA)**2+(QA+DB)**2+AQD)
   QXYQXY=FD/SQRT(R**2+TD*(QA-QB)**2+AQQ)+FD/SQRT(R**2+TD*(QA+QB)**2+AQQ)-EIGHT/SQRT(R**2+TD*(QA**2+QB**2)+AQQ)
   QXZQXZ=TD/SQRT((R+QA-QB)**2+(QA-QB)**2+AQQ)-TD/SQRT((R+QA+QB)**2+( & 
          QA-QB)**2+AQQ)-TD/SQRT((R-QA-QB)**2+(QA-QB)**2+AQQ)+TD/SQRT((R-QA+  &
          QB)**2+(QA-QB)**2+AQQ)-TD/SQRT((R+QA-QB)**2+(QA+QB)**2+AQQ)+TD/SQRT &
          ((R+QA+QB)**2+(QA+QB)**2+AQQ)+TD/SQRT((R-QA-QB)**2+(QA+QB)**2+AQQ &
          )-TD/SQRT((R-QA+QB)**2+(QA+QB)**2+AQQ)
   DXDX=DXDX/FOUR
   DZDZ=DZDZ/FOUR
   DZQXX=DZQXX/EIGHT
   QXXDZ=QXXDZ/EIGHT
   DZQZZ=DZQZZ/EIGHT
   QZZDZ=QZZDZ/EIGHT
   DXQXZ=DXQXZ/EIGHT
   QXZDX=QXZDX/EIGHT
   QXXQXX=QXXQXX/16.0_DP
   QXXQYY=QXXQYY/16.0_DP
   QXXQZZ=QXXQZZ/16.0_DP
   QZZQXX=QZZQXX/16.0_DP
   QZZQZZ=QZZQZZ/16.0_DP
   QXZQXZ=QXZQXZ/16.0_DP
   QXYQXY=QXYQXY/16.0_DP
   RI(6)=DZDZ
   RI(7)=DXDX
   RI(8)=-EDZ-QZZDZ
   RI(9)=-EDZ-QXXDZ
   RI(10)=-QXZDX
   RI(13)=-DZE-DZQZZ
   RI(14)=-DZE-DZQXX
   RI(15)=-DXQXZ
   RI(16)=EE+EQZZ+QZZE+QZZQZZ
   RI(17)=EE+EQZZ+QXXE+QXXQZZ
   RI(18)=EE+EQXX+QZZE+QZZQXX
   RI(19)=EE+EQXX+QXXE+QXXQXX
   RI(20)=QXZQXZ
   RI(21)=EE+EQXX+QXXE+QXXQYY
   RI(22)=0.5_DP*(QXXQXX-QXXQYY)

   20 CONTINUE

   ! . Convert the integrals to eV.
   RI(2:22) = RI(2:22) * AU_TO_EV

   ! . Calculate the electron-nuclear attraction integrals.
   CORINT(2,1) = CORE(NJ) * RI(2)
   CORINT(3,1) = CORE(NJ) * RI(3)
   CORINT(4,1) = CORE(NJ) * RI(4)
   CORINT(2,2) = CORE(NI) * RI(5)
   CORINT(3,2) = CORE(NI) * RI(11)
   CORINT(4,2) = CORE(NI) * RI(12)

   END SUBROUTINE REPP

   !-------------------------------------------------------------------
   SUBROUTINE ROTATEI ( NI, NJ, XI, XJ, W, KR, E1B, E2A, ENUC, QAM1PM3 )
   !-------------------------------------------------------------------

   !
   !     ROTATE calculates the two-electron integrals and the electron-nuclear
   !     attraction integrals. The integrals are evaluated first in the local
   !     frame by REPP and then transformed to the molecular frame by ROTATE.
   !
   !     Input:     NI, NJ      atomic numbers of first and second atoms.
   !                XI, XJ      coordinates of atoms (in atomic units).
   !
   !     Output:    W           array of two-electron integrals.
   !                E1B         attraction integrals for electron of atom NI
   !                            and nucleus of NJ.
   !                E2A         attraction integrals for electron of atom NJ
   !                            and nucleus of NI.
   !                ENUC        nuclear-nuclear repulsion term.
   !
   !
   !     The integrals in the local frame are stored in RI and they have the
   !     order :
   !
   !         (ss/ss) = 1,     (so/ss) = 2,   (oo/ss) = 3,   (pp/ss) = 4,
   !         (ss/os) = 5,     (so/so) = 6,   (sp/sp) = 7,   (oo/so) = 8,
   !         (pp/so) = 9,     (po/sp) = 10,  (ss/oo) = 11,  (ss/pp) = 12,
   !         (so/oo) = 13,    (so/pp) = 14,  (sp/op) = 15,  (oo/oo) = 16,
   !         (pp/oo) = 17,    (oo/pp) = 18,  (pp/pp) = 19,  (po/po) = 20,
   !                       (pp/p*p*) = 21,  (p*p/p*p) = 22
   !
   !     where o is a p-sigma orbital and p and p* are p-pi orbitals.
   !
   !     The nuclear attraction integrals are stored in CORINT(kl/ij) in the
   !     order :
   !              (ss/) = 1,  (so/) = 2, (oo/) = 3, (pp/) = 4
   !
   !     where ij = 1 if the orbitals are on atom NI and ij = 2 if they are
   !     on atom NJ.
   !

   ! . Scalar argument declarations.
   INTEGER,            INTENT(IN)    :: NI, NJ
   INTEGER,            INTENT(INOUT) :: KR
   LOGICAL,            INTENT(IN)    :: QAM1PM3
   REAL ( KIND = DP ), INTENT(OUT)   :: ENUC

   ! . Array argument declarations.
   REAL ( KIND = DP ), DIMENSION(1:10), INTENT(OUT)   :: E1B, E2A
   REAL ( KIND = DP ), DIMENSION(1:3),  INTENT(IN)    :: XI, XJ
   REAL ( KIND = DP ), DIMENSION(:),    INTENT(INOUT) :: W

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: PDDG_EXPONENT = 10.0_DP

   ! . Local scalars.
   INTEGER            :: I, IB, IG, II, IJ, IK, J, JB, JJ, JK, K, KI, KK, L, LL, NT
   LOGICAL            :: SI, SK
   REAL ( KIND = DP ) :: A, AX, D, GAM, RIJ, SCALE, ZAF, ZBF

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:4,1:2) :: CORINT
   REAL ( KIND = DP ), DIMENSION(1:22)    :: RI
   REAL ( KIND = DP ), DIMENSION(1:3)     :: X, Y, Z

   ! . Equivalence variables for CORINT.
   REAL ( KIND = DP ) :: CSS1, CSP1, CPPS1, CPPP1, CSS2, CSP2, CPPS2, CPPP2

   ! . Equivalence variables for X, Y and Z.
   REAL ( KIND = DP ) :: X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3

   ! . Local parameter declarations.
   REAL ( KIND = DP ), PARAMETER :: APROX1 = 1.0_DP - 1.0E-8_DP

   ! . Find the internuclear distance in bohrs.
   X   = ANGSTROMS_TO_BOHRS * ( XI - XJ )
   RIJ = SQRT ( DOT_PRODUCT ( X, X ) )

   ! . Calculate the integrals in the local frame.
   CALL REPP ( NI, NJ, RIJ, RI, CORINT )

   ! . Save the first integral for the electron-core attractions later.
   GAM = RI(1)

   ! . The repulsion integrals in the molecular frame are stored in the
   ! . order in which they are used : (IJ/KL) where J .LE. I and L .LE. K
   ! . and L varies most rapidly and I least rapidly.
   A = 1.0_DP / RIJ
   X = A * X
   Z(3)=0.0_DP
   IF ( ABS ( X(3) ) > APROX1 ) THEN
      X(3)=SIGN(1.0_DP,X(3))
      Y(1)=0.0_DP
      Y(2)=1.0_DP
      Y(3)=0.0_DP
      Z(1)=1.0_DP
      Z(2)=0.0_DP
   ELSE
      Z(3)=SQRT(1.0_DP-X(3)**2)
      A=1.0_DP/Z(3)
      Y(1)=-A*X(2)*SIGN(1.0_DP,X(1))
      Y(2)=ABS(A*X(1))
      Y(3)=0.0_DP
      Z(1)=-A*X(1)*X(3)
      Z(2)=-A*X(2)*X(3)
   END IF

   ! . Replace the equivalence statement for X, Y and Z by specific assignments.
   X1 = X(1)
   X2 = X(2)
   X3 = X(3)
   Y1 = Y(1)
   Y2 = Y(2)
   Y3 = Y(3)
   Z1 = Z(1)
   Z2 = Z(2)
   Z3 = Z(3)

   IB=MIN(NATORB(NI),4)
   JB=MIN(NATORB(NJ),4)
   KI=0
   DO I=1,IB
      SI=I.EQ.1
      II=I-1
      DO J=1,I
         JJ=J-1
         IJ=0
         IF (JJ.EQ.0) IJ=-1
         IF (SI) IJ=+1
         DO K=1,JB
            KK=K-1
            SK=KK.GT.0
            DO L=1,K
               KI=KI+1
               IF (SK) GOTO 180
               ! . Integral is of the type (IJ/SS).
               IF ( IJ < 0 ) THEN
                  ! . (PS/SS).
                  W(KI)=RI(2)*X(II)
               ELSE IF ( IJ == 0 ) THEN
                  ! . (PP/SS).
                  W(KI)=RI(3)*X(II)*X(JJ)+RI(4)*(Y(II)*Y(JJ)+Z(II)*Z(JJ))
               ELSE
                  ! . (SS/SS).
                  W(KI)=RI(1)
               END IF
               GO TO 260
180             LL=L-1
               IF (LL.GT.0) GOTO 220
               ! . Integral is of the type (IJ/PS).
               IF ( IJ < 0 ) THEN
                  ! . (PS/PS).
                  W(KI)=RI(6)*X(II)*X(KK)+RI(7)*(Y(II)*Y(KK)+Z(II)*Z(KK))
               ELSE IF ( IJ == 0 ) THEN
                  ! . (PP/PS).
                  W(KI)=X(KK)*(RI(8)*X(II)*X(JJ)+RI(9)*(Y(II)*Y(JJ)+Z(II)*Z(JJ))) &
                        +RI(10)*(X(II)*(Y(JJ)*Y(KK)+Z(JJ)*Z(KK))+X(JJ)*(Y(II)*Y(KK)+Z(II)*Z(KK)))
               ELSE
                  ! . (SS/PS).
                  W(KI)=RI(5)*X(KK)
               END IF
               GO TO 260
               ! . Integral is of the type (IJ/PP).
220             CONTINUE
               IF ( IJ < 0 ) THEN
                  ! . (PS/PP).
                  W(KI)=X(II)*(RI(13)*X(KK)*X(LL)+RI(14)*(Y(KK)*Y(LL)+Z(KK)*Z(LL) &
                        ))+RI(15)*(Y(II)*(Y(KK)*X(LL)+Y(LL)*X(KK))+Z(II)*(Z(KK)*X(LL)+Z(LL)*X(KK)))
                ELSE IF ( IJ == 0 ) THEN
                  ! . (PP/PP).
                  W(KI)=(RI(16)*X(II)*X(JJ)+RI(17)*(Y(II)*Y(JJ)+Z(II)*Z(JJ)))*X(KK)*X(LL)+   &
                        RI(18)*X(II)*X(JJ)*(Y(KK)*Y(LL)+Z(KK)*Z(LL))+                        &
                        RI(19)*(Y(II)*Y(JJ)*Y(KK)*Y(LL)+Z(II)*Z(JJ)*Z(KK)*Z(LL))+RI(20)      &
                        *(X(II)*(   X(KK)*(Y(JJ)*Y(LL)+Z(JJ)*Z(LL))+X(LL)*(Y(JJ)*Y(KK)+Z(JJ) &
                        *Z(KK))   )+X(JJ)*(X(KK)*(Y(II)*Y(LL)+Z(II)*Z(LL))+X(LL)*(Y(II)*     &
                        Y(KK)+Z(II)*Z(KK))))+RI(21)*(Y(II)*Y(JJ)*Z(KK)*Z(LL)+Z(II)*Z(JJ)     &
                        *Y(KK)*Y(LL))+RI(22)*(Y(II)*Z(JJ)+Z(II)*Y(JJ))*(Y(KK)*Z(LL)+Z(KK)*Y(LL))
               ELSE
                  ! . (SS/PP).
                  W(KI)=RI(11)*X(KK)*X(LL)+RI(12)*(Y(KK)*Y(LL)+Z(KK)*Z(LL))
               END IF
260             CONTINUE
            END DO
         END DO
      END DO
   END DO
280 CONTINUE

   ! . Replace the equivalence statement for CORINT by specific assignments.
   CSS1  = CORINT(1,1)
   CSP1  = CORINT(2,1)
   CPPS1 = CORINT(3,1)
   CPPP1 = CORINT(4,1)
   CSS2  = CORINT(1,2)
   CSP2  = CORINT(2,2)
   CPPS2 = CORINT(3,2)
   CPPP2 = CORINT(4,2)

   ! . Compute the electron-core attraction terms.
   E1B = 0.0_DP
   E2A = 0.0_DP
   IF ( NATORB(NI) > 0 ) THEN
   E1B(1)=-CSS1
   IF(NI.GT.1) THEN
      E1B(2) = -CSP1 *X1
      E1B(3) = -CPPS1*X1**2-CPPP1*(Y1**2+Z1**2)
      E1B(4) = -CSP1 *X2
      E1B(5) = -CPPS1*X1*X2-CPPP1*(Y1*Y2+Z1*Z2)
      E1B(6) = -CPPS1*X2*X2-CPPP1*(Y2*Y2+Z2*Z2)
      E1B(7) = -CSP1 *X3
      E1B(8) = -CPPS1*X1*X3-CPPP1*(Y1*Y3+Z1*Z3)
      E1B(9) = -CPPS1*X2*X3-CPPP1*(Y2*Y3+Z2*Z3)
      E1B(10)= -CPPS1*X3*X3-CPPP1*(Y3*Y3+Z3*Z3)
   END IF
   END IF
   IF ( NATORB(NJ) > 0 ) THEN
   E2A(1)=-CSS2
   IF(NJ.GT.1) THEN
      E2A(2) = -CSP2 *X1
      E2A(3) = -CPPS2*X1**2-CPPP2*(Y1**2+Z1**2)
      E2A(4) = -CSP2 *X2
      E2A(5) = -CPPS2*X1*X2-CPPP2*(Y1*Y2+Z1*Z2)
      E2A(6) = -CPPS2*X2*X2-CPPP2*(Y2*Y2+Z2*Z2)
      E2A(7) = -CSP2 *X3
      E2A(8) = -CPPS2*X1*X3-CPPP2*(Y1*Y3+Z1*Z3)
      E2A(9) = -CPPS2*X2*X3-CPPP2*(Y2*Y3+Z2*Z3)
      E2A(10)= -CPPS2*X3*X3-CPPP2*(Y3*Y3+Z3*Z3)
   END IF
   END IF
   ! . Compute the core-core repulsion terms.
      IF(ABS(CORE(NI)).GT.20.AND.ABS(CORE(NJ)).GT.20) THEN
         ENUC=0.0_DP
         RETURN
      ELSE IF (RIJ.LT.1.0_DP.AND.NATORB(NI)*NATORB(NJ).EQ.0) THEN
         ENUC=0.0_DP
         RETURN
      END IF
      ! . Convert RIJ to Angstroms.
      RIJ   = RIJ / ANGSTROMS_TO_BOHRS
      SCALE = EXP(-ALP(NI)*RIJ)+EXP(-ALP(NJ)*RIJ)
      NT=NI+NJ
      IF(NT.EQ.8.OR.NT.EQ.9) THEN
         IF(NI.EQ.7.OR.NI.EQ.8) SCALE=SCALE+(RIJ-1.0_DP)*EXP(-ALP(NI)*RIJ)
         IF(NJ.EQ.7.OR.NJ.EQ.8) SCALE=SCALE+(RIJ-1.0_DP)*EXP(-ALP(NJ)*RIJ)
      END IF
      ENUC = CORE(NI)*CORE(NJ)*GAM
      SCALE=ABS(SCALE*ENUC)
      ! . Compute the AM1/PM3 specific terms.
      IF ( QAM1PM3 ) THEN
         DO IG = 1,4
            IF ( ABS ( FN1(NI,IG) ) > 0.0_DP ) THEN
               SCALE = SCALE + (CORE(NI)*CORE(NJ)/RIJ)*FN1(NI,IG)*EXP ( MAX ( -30.0_DP, -FN2(NI,IG)*(RIJ-FN3(NI,IG))**2 ) )
            END IF
            IF ( ABS ( FN1(NJ,IG) ) > 0.0_DP )  THEN
               SCALE = SCALE + (CORE(NI)*CORE(NJ)/RIJ)*FN1(NJ,IG)*EXP ( MAX ( -30.0_DP, -FN2(NJ,IG)*(RIJ-FN3(NJ,IG))**2 ) )
            END IF
         END DO
      END IF

      ! . PDDG specific terms.
      IF ( HAMILTONIAN == "PDDG" ) THEN 
          ZAF = CORE(NI)/(CORE(NI) + CORE(NJ))
          ZBF = CORE(NJ)/(CORE(NI) + CORE(NJ))
          DO IK= 1,2   
             DO JK = 1,2
                D = RIJ - PDDGE(NI,IK) - PDDGE(NJ,JK)
                AX = PDDG_EXPONENT * D * D
!write ( 6, "(2i10,3f20.10)" ) ik, jk, d, ax, (ZAF * PDDGC(NI,IK) + ZBF * PDDGC(NJ,JK))
                SCALE = SCALE + (ZAF * PDDGC(NI,IK) + ZBF * PDDGC(NJ,JK)) * EXP(-AX)
             END DO
          END DO
      END IF

      ENUC = ENUC + SCALE

   ! . Increment the index to the two-electron integral array.
   KR = KR + KI

   END SUBROUTINE ROTATEI

!------------------------------------------------------------------------
! . SEINTC Subroutines.
!------------------------------------------------------------------------

   !------------------------------------------------------
   SUBROUTINE SEINTC_INTEGRALAB ( R, A, B, SWITCH, F, DF )
   !------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN)  :: R, A, B, SWITCH
   REAL ( KIND = DP ), INTENT(OUT) :: F, DF

   ! . Local scalars.
   REAL ( KIND = DP ) :: IOFF, ION, IMODOFF, IMODON, IMODR, IR

   ! . Calculate IMODOFF.
   CALL SEINTC_TNAB ( CUTOFFB, A, B, IMODOFF, IOFF )

   ! . The interaction is close range.
   IF ( R <= CUTONB ) THEN

      ! . Calculate IR.
      IR = 1.0_DP / SQRT ( ( R + A )**2 + B )

      ! . Calculate IMODON.
      CALL SEINTC_TNAB ( CUTONB, A, B, IMODON, ION )

      ! . Calculate the potential.
      F = IR - ( IMODOFF - IMODON ) - ION

   ! . The interaction is in the switching region.
   ELSE

      ! . Calculate IMODR.
      CALL SEINTC_TNAB ( R, A, B, IMODR, IR )

      ! . Calculate the potential.
      F = IMODR - IMODOFF

   END IF

   ! . Calculate the derivative.
   DF = - SWITCH * ( R + A ) * IR * IR * IR

   END SUBROUTINE SEINTC_INTEGRALAB

   !--------------------------------------------------
   SUBROUTINE SEINTC_INTEGRALB ( R, B, SWITCH, F, DF )
   !--------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN)  :: R, B, SWITCH
   REAL ( KIND = DP ), INTENT(OUT) :: F, DF

   ! . Local scalars.
   REAL ( KIND = DP ) :: IOFF, ION, IMODOFF, IMODON, IMODR, IR

   ! . Calculate IMODOFF.
   CALL SEINTC_TNB ( CUTOFFB, B, IMODOFF, IOFF )

   ! . The interaction is close range.
   IF ( R <= CUTONB ) THEN

      ! . Calculate IR.
      IR = 1.0_DP / SQRT ( R**2 + B )

      ! . Calculate IMODON.
      CALL SEINTC_TNB ( CUTONB, B, IMODON, ION )

      ! . Calculate the potential.
      F = IR - ( IMODOFF - IMODON ) - ION

   ! . The interaction is in the switching region.
   ELSE

      ! . Calculate IMODR.
      CALL SEINTC_TNB ( R, B, IMODR, IR )

      ! . Calculate the potential.
      F = IMODR - IMODOFF

   END IF

   ! . Calculate the derivative.
   DF = - SWITCH * R * IR * IR * IR

   END SUBROUTINE SEINTC_INTEGRALB

   !----------------------
   SUBROUTINE SEINTC_SETUP
   !----------------------

   ! . Calculate some cutoff factors.
   CUTOFFB = ANGSTROMS_TO_BOHRS * CUT_OFF
   CUTONB  = ANGSTROMS_TO_BOHRS * CUT_ON
   R2OFFB  = CUTOFFB**2
   R2ONB   = CUTONB**2
   GAMMAB  = ( R2OFFB - R2ONB )**3

   ! . Calculate the TN prefactors.
   T6FAC =   2.0_DP / GAMMAB
   T4FAC = - 3.0_DP * ( R2OFFB + R2ONB ) / GAMMAB
   T2FAC =   6.0_DP * R2OFFB * R2ONB / GAMMAB
   T0FAC =   R2OFFB * R2OFFB * ( R2OFFB - 3.0_DP * R2ONB ) / GAMMAB

   END SUBROUTINE SEINTC_SETUP

   !-------------------------------------------
   SUBROUTINE SEINTC_TNAB ( R, A, B, IMOD, T0 )
   !-------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN)  :: R, A, B
   REAL ( KIND = DP ), INTENT(OUT) :: IMOD, T0

   ! . Local scalars.
   REAL ( KIND = DP ) :: A2, B2, FACT, LOGFAC, R2, RA, RASQB, T2, T4, T6

   ! . Calculate some intermediate factors for the TN calculations.
   A2     = A * A
   B2     = B * B
   R2     = R * R
   RA     = R + A
   RASQB  = RA * RA + B
   FACT   = SQRT ( RASQB )
   LOGFAC = LOG ( RA + FACT )

   ! . T0.
   T0 = 1.0_DP / FACT

   ! . T2.
   T2 = - ( 2.0_DP * A2 + 2.0_DP * B + 4.0_DP * A * R + R2 ) * T0 + 2.0_DP * A * LOGFAC

   ! . T4.
   T4 = ( - 22.0_DP * A2 * A2 - 14.0_DP * A2 * B + 8.0_DP * B2 - 34.0_DP * A2 * A * R + 26.0_DP * A * B * R &
                     - 6.0_DP * A2 * R2 + 4.0_DP * B * R2 + 2.0_DP * A * R2 * R - R2 * R2 ) * T0 / 3.0_DP + &
                                                        2.0_DP * A * ( 2.0_DP * A2 - 3.0_DP * B ) * LOGFAC

   ! . T6.
   T6 = ( - 274.0_DP * A2 * A2 * A2 + 333.0_DP * A2 * A2 * B + 543.0_DP * A2 * B2 - 64.0_DP * B2 * B -                 &
            394.0_DP * A2 * A2 * A * R + 1207.0_DP * A2 * A * B * R - 289.0_DP * A * B2 * R - 60.0_DP * A2 * A2 * R2 + &
            223.0_DP * A2 * B * R2 - 32.0_DP * B2 * R2 + 20.0_DP * A2 * A * R2 * R - 43.0_DP * A * B * R2 * R -        &
            10.0_DP * A2 * R2 * R2 + 8.0_DP * B * R2 * R2 + 6.0_DP * A * R2 * R2 * R - 4.0_DP * R2 * R2 * R2 ) * T0 /  &
            20.0_DP + 3.0_DP * A * ( 8.0_DP * A2 * A2 - 40.0_DP * A2 * B + 15.0_DP * B2 ) * LOGFAC / 4.0_DP

   ! . Calculate IMOD.
   IMOD = T0FAC * T0 + T2FAC * T2 + T4FAC * T4 + T6FAC * T6

   END SUBROUTINE SEINTC_TNAB

   !---------------------------------------
   SUBROUTINE SEINTC_TNB ( R, B, IMOD, T0 )
   !---------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN)  :: R, B
   REAL ( KIND = DP ), INTENT(OUT) :: IMOD, T0

   ! . Local scalars.
   REAL ( KIND = DP ) :: B2, FACT, R2, RSQB, T2, T4, T6

   ! . Calculate some intermediate factors for the TN calculations.
   B2   = B * B
   R2   = R * R
   RSQB = R * R + B
   FACT = SQRT ( RSQB )

   ! . T0.
   T0 = 1.0_DP / FACT

   ! . T2.
   T2 = - ( 2.0_DP * B + R2 ) * T0

   ! . T4.
   T4 = ( 8.0_DP * B2 + 4.0_DP * B * R2 - R2 * R2 ) * T0 / 3.0_DP

   ! . T6.
   T6 = ( - 64.0_DP * B2 * B - 32.0_DP * B2 * R2 + 8.0_DP * B * R2 * R2 - 4.0_DP * R2 * R2 * R2 ) * T0 / 20.0_DP

   ! . Calculate IMOD.
   IMOD = T0FAC * T0 + T2FAC * T2 + T4FAC * T4 + T6FAC * T6

   END SUBROUTINE SEINTC_TNB

END MODULE MOPAC_INTEGRALS

