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
!                          The MOPAC Fock Matrix Module
!===============================================================================
!
! . Subroutines:
!
!   FOCK_MATRIX_RESTRICTED        Calculate a spin-restricted Fock matrix.
!   FOCK_MATRIX_UNRESTRICTED      Calculate a spin-unrestricted Fock matrix.
!
!===============================================================================
MODULE MOPAC_FOCK_MATRIX

! . Module declarations.
USE DEFINITIONS,      ONLY : DP

USE ATOMS,            ONLY : NATOMSQM
USE MOPAC_DATA,       ONLY : JIND2, JIND3, KIND2, NBASTR, NPIBEADS, PISETEI, QMATOM, SETEI
USE MOPAC_PARAMETERS, ONLY : GPP, GP2, GSP, GSS, HSP, NATORB

IMPLICIT NONE
PRIVATE
PUBLIC :: FOCK_MATRIX_RESTRICTED, FOCK_MATRIX_UNRESTRICTED

!===============================================================================
CONTAINS
!===============================================================================

   !----------------------------------------------------------
   SUBROUTINE FOCK_MATRIX_RESTRICTED ( NCALC, DENCLS, FCKCLS )
   !----------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: NCALC

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NBASTR), INTENT(IN)  :: DENCLS
   REAL ( KIND = DP ), DIMENSION(1:NBASTR), INTENT(OUT) :: FCKCLS

   ! . Local scalars.
   INTEGER            :: I, IAB, ICOPY, IFIRST, ILAST, INDTEI, IP, IQM, I1, I2, J, JAB, JCOPY, JFIRST, &
                         JJ, JLAST, JP, JQM, KA, KFI, KFJ, L, M, NI, NJ, P, PIINDTEI
   REAL ( KIND = DP ) :: PTPOP

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)     :: FDIAG
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:,:) :: PDIAG

   !----------------------------------------------------------------------------
   ! . Calculate the one-centre part of the Fock matrix.
   !----------------------------------------------------------------------------
   ! . Initialize the Fock matrix.
   FCKCLS = 0.0_DP

   ! . Allocate the temporary arrays.
   ALLOCATE ( FDIAG(10*NATOMSQM), PDIAG(1:4,1:4,1:NATOMSQM) )

   ! . Set the Fock matrix counter.
   IF ( NPIBEADS == 1 ) THEN
      P = 0
   ELSE
      P = NCALC
   END IF

   ! . Compute the diagonal terms for the Fock matrix involving two-electron
   ! . integrals on the same atom (these are stored and are not calculated).
   DO IQM = 1,NATOMSQM

      ! . Get some information for the atom.
      ICOPY  = QMATOM(IQM)%COPY
      IFIRST = QMATOM(IQM)%BFIRST
      ILAST  = QMATOM(IQM)%BLAST
      NI     = QMATOM(IQM)%NUMBER
      IF ( NATORB(NI) == 0 ) CYCLE

      ! . Skip atoms of the wrong index.
      IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

      ! . Calculate the density matrix contributions.
      PTPOP  = 0.0_DP

      ! . Hydrogen.
      IF ( NATORB(NI) <= 1 ) THEN
         CONTINUE
      ! . SP elements.
      ELSE
         PTPOP = DENCLS((ILAST*(ILAST+1))/2) + DENCLS(((ILAST-1)*(ILAST))/2) + DENCLS(((ILAST-2)*(ILAST-1))/2)
      END IF

      ! . ss diagonal terms.
      KA = (IFIRST*(IFIRST+1))/2
      FCKCLS(KA) = FCKCLS(KA) + ( 0.5_DP * DENCLS(KA) * GSS(NI) ) + PTPOP * ( GSP(NI) - 0.5_DP * HSP(NI) )
      IF ( NI < 3 ) CYCLE

      L = KA
      DO J = (IFIRST+1),ILAST
         M = L+IFIRST
         L = L+J
         ! . pp diagonal terms.
         FCKCLS(L) = FCKCLS(L) + DENCLS(KA) * (GSP(NI)-0.5_DP*HSP(NI)) + 0.5_DP * DENCLS(L) * GPP(NI) + &
                     ( PTPOP - DENCLS(L) ) * GP2(NI) - 0.25_DP * (PTPOP-DENCLS(L)) * (GPP(NI)-GP2(NI))
         ! . sp terms.
         FCKCLS(M) = FCKCLS(M) + ( 2.0_DP * DENCLS(M) * HSP(NI) ) - 0.5_DP * DENCLS(M) * ( HSP(NI)+GSP(NI) )
      END DO

      ! . pp* terms.
      DO J = (IFIRST+1),(ILAST-1)
         DO L = (J+1),ILAST
            M = (L*(L-1))/2+J
            FCKCLS(M) = FCKCLS(M) + DENCLS(M) * (GPP(NI)-GP2(NI)) - 0.25_DP * DENCLS(M) * (GPP(NI)+GP2(NI))
         END DO
      END DO
   END DO

   !----------------------------------------------------------------------------
   ! . Calculate the two-centre part of the Fock matrix.
   !----------------------------------------------------------------------------
   ! . Initialize some counters.
   INDTEI   = 1
   KFI      = 1
   PIINDTEI = 1
   FDIAG    = 0.0_DP
   PDIAG    = 0.0_DP

   ! . Outer loop over QM atoms.
   DO IQM = 1,NATOMSQM

      ! . Get some information for the atom.
      ICOPY  = QMATOM(IQM)%COPY
      IFIRST = QMATOM(IQM)%BFIRST
      ILAST  = QMATOM(IQM)%BLAST
      IAB    = ILAST - IFIRST + 1
      NI     = QMATOM(IQM)%NUMBER
      IF ( NATORB(NI) == 0 ) CYCLE

      ! . Skip atoms of the wrong index.
      IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

      IF ( IAB == 1 ) THEN
         PDIAG(1,1,IQM) = DENCLS((IFIRST*(IFIRST+1))/2)
      ELSE
         I1 = ((IFIRST-1)*(IFIRST-2))/2
         IP = 0
         DO I = IFIRST,ILAST
            IP = IP + 1
            I1 = I1 + I - 1
            I2 = I1 + IFIRST - 1
            JP = 0
            DO J = IFIRST,I
               JP = JP + 1
               I2 = I2 + 1
               PDIAG(JP,IP,IQM) = DENCLS(I2)
               PDIAG(IP,JP,IQM) = DENCLS(I2)
            END DO
         END DO
      END IF
      KFJ = 1

      ! . Inner loop over QM atoms.
      DO JQM = 1,(IQM-1)

         ! . Get some information for the atom.
         JCOPY  = QMATOM(JQM)%COPY
         JFIRST = QMATOM(JQM)%BFIRST
         JLAST  = QMATOM(JQM)%BLAST
         JAB    = JLAST - JFIRST + 1
         NJ     = QMATOM(JQM)%NUMBER
         IF ( NATORB(NJ) == 0 ) CYCLE

         ! . Skip atoms of the wrong index.
         IF ( ( JCOPY /= 0 ) .AND. ( JCOPY /= P ) ) CYCLE

         ! . Both atoms are non-PI atoms.
         IF ( ( ICOPY == 0 ) .AND. ( JCOPY == 0 ) ) THEN
            CALL SEFOCK ( FDIAG(KFI:), FDIAG(KFJ:),   SETEI(  INDTEI:),     INDTEI, PDIAG(:,:,IQM), PDIAG(:,:,JQM) )
         ! . One or both atoms are PI atoms.
         ELSE
            CALL SEFOCK ( FDIAG(KFI:), FDIAG(KFJ:), PISETEI(PIINDTEI:,P), PIINDTEI, PDIAG(:,:,IQM), PDIAG(:,:,JQM) )
         END IF

         ! . Increment the counter for FDIAG.
         KFJ = KFJ + (JAB*(JAB+1))/2

      END DO

      ! . Increment the counter for FDIAG.
      KFI = KFI + (IAB*(IAB+1))/2

   END DO

   ! . Add in the diagonal parts of the Fock matrix.
   JJ = 0
   DO IQM = 1,NATOMSQM

      ! . Get some information for the atom.
      ICOPY  = QMATOM(IQM)%COPY
      IFIRST = QMATOM(IQM)%BFIRST
      ILAST  = QMATOM(IQM)%BLAST
      I1     = ((IFIRST-1)*(IFIRST-2))/2
      NI     = QMATOM(IQM)%NUMBER
      IF ( NATORB(NI) == 0 ) CYCLE

      ! . Skip atoms of the wrong index.
      IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

      ! . Loop over the atom basis functions.
      DO I = IFIRST,ILAST
         I1 = I1 + I - 1
         I2 = I1 + IFIRST - 1
         DO J = IFIRST,I
            I2 = I2 + 1
            JJ = JJ + 1
            FCKCLS(I2) = FCKCLS(I2) + FDIAG(JJ)
         END DO
      END DO
   END DO

   ! . Deallocate space.
   DEALLOCATE ( FDIAG, PDIAG )

   !============================================================================
   CONTAINS
   !============================================================================

      !----------------------------------------------------
      SUBROUTINE SEFOCK ( FII, FJJ, SETEI, KR, P2II, P2JJ )
      !----------------------------------------------------

      !
      !     SEFCC1 constructs the two-centre two-electron part of the Fock
      !     matrix given :
      !
      !     F - the partial packed Fock matrix, FII - the atom ii part of the
      !     Fock matrix, FJJ - the atom jj part of the Fock matrix, W - the
      !     two-centre two-electron integrals and P - the density matrix.
      !
      !     Atom indices are IFIRST to ILAST for atom I and JFIRST to JLAST
      !     for atom J.
      !

      ! . Scalar argument declarations.
      INTEGER, INTENT(INOUT) :: KR

      ! . Array argument declarations.
      REAL ( KIND = DP ), DIMENSION(:),    INTENT(INOUT) :: FII,  FJJ
      REAL ( KIND = DP ), DIMENSION(1:16), INTENT(IN)    :: P2II, P2JJ
      REAL ( KIND = DP ), DIMENSION(:),    INTENT(IN)    :: SETEI

      ! . Local scalars.
      INTEGER            :: I, IAC, IAC2, II, IJ, IJAC, J, JAC, JAC2, JJ, K, KK, KL, NTYPE
      REAL ( KIND = DP ) :: SUM

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:16) :: PAIJ

      ! . Establish NTYPE. If 1 both atoms in the diatomic have 4 orbitals,
      ! . if 2 one atom has 4 orbitals and the other 1 and if 3 both atoms
      ! . have single orbitals.
      IAC = ILAST - IFIRST + 1
      JAC = JLAST - JFIRST + 1
      IF ( ( IAC + JAC ) == 8 ) THEN
         NTYPE = 1
      ELSE IF ( ( IAC + JAC ) == 5 ) THEN
         NTYPE = 2
      ELSE
         IJ = (IFIRST*(IFIRST-1))/2 + JFIRST
         FCKCLS(IJ) = FCKCLS(IJ) - ( 0.5_DP * SETEI(1) * DENCLS(IJ) )
         FII(1)   = FII(1) + ( SETEI(1) * P2JJ(1) )
         FJJ(1)   = FJJ(1) + ( SETEI(1) * P2II(1) )
         KR = KR + 1
         RETURN
      END IF

      ! . Unpack off-diagonal spin density matrix elements.
      KK = ((IFIRST-1)*(IFIRST-2))/2
      DO I = IFIRST,ILAST
         KK = KK + I - 1
         JJ = ( I - IFIRST ) * JAC
         K  = KK + JFIRST - 1
         IF ( JAC == 1 ) THEN
            PAIJ(JJ+1) = 0.5_DP * DENCLS(K+1)
         ELSE
            DO J = JFIRST,JLAST
               JJ = JJ + 1
               K  = K  + 1
               PAIJ(JJ) = 0.5_DP * DENCLS(K)
            END DO
         ENDIF
      END DO

      ! . Construct the off-diagonal exchange contributions.
      IJAC = IAC * JAC
      II   = ((IFIRST-1)*(IFIRST-2))/2
      DO I = 1,IAC
         II = II + IFIRST + I - 2
         IJ = II + JFIRST - 1
         DO J = 1,JAC
            IJ  = IJ + 1
            SUM = 0.0_DP
            DO KL = 1,IJAC
               SUM = SUM + ( SETEI( KIND2(KL,I,J,NTYPE) ) * PAIJ(KL) )
            END DO
            FCKCLS(IJ) = FCKCLS(IJ) - SUM
         END DO
      END DO

      ! . Construct the on-diagonal Coulombic contributions.
      ! . Atom I affected by atom J.
      JAC2 = JAC * JAC
      JJ   = 0
      DO I = 1,IAC
         DO J = 1,I
            JJ = JJ + 1
            IF ( JAC2 == 1 ) THEN
               FII(JJ) = FII(JJ) + ( SETEI( JIND2(1,I,J,2) ) * P2JJ(1) )
            ELSE
               SUM = 0.0_DP
               DO KL = 1,JAC2
                  SUM = SUM + ( SETEI( JIND2(KL,I,J,NTYPE) ) * P2JJ(KL) )
               END DO
               FII(JJ) = FII(JJ) + SUM
            END IF
         END DO
      END DO

      ! . Atom J affected by atom I.
      IAC2 = IAC * IAC
      JJ   = 0
      DO I = 1,JAC
         DO J = 1,I
            JJ = JJ + 1
            IF ( IAC2 == 1 ) THEN
               FJJ(JJ) = FJJ(JJ) + ( SETEI( JIND3(I,J,1,2) ) * P2II(1) )
            ELSE
               SUM = 0.0_DP
               DO KL = 1,IAC2
                  SUM = SUM + ( SETEI( JIND3(I,J,KL,NTYPE) ) * P2II(KL) )
               END DO
               FJJ(JJ) = FJJ(JJ) + SUM
            END IF
         END DO
      END DO

      ! . Increment KR counter to line up W matrix for the next pair of atoms.
      KR = KR + ( ( IAC*(IAC+1) * JAC*(JAC+1) ) / 4 )

      END SUBROUTINE SEFOCK

   END SUBROUTINE FOCK_MATRIX_RESTRICTED

   !------------------------------------------------------------------------------------
   SUBROUTINE FOCK_MATRIX_UNRESTRICTED ( NCALC, DENTOT, DENALP, DENBET, FCKALP, FCKBET )
   !------------------------------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: NCALC

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:NBASTR), INTENT(IN)  :: DENALP, DENBET, DENTOT
   REAL ( KIND = DP ), DIMENSION(1:NBASTR), INTENT(OUT) :: FCKALP, FCKBET

   ! . Local scalars.
   INTEGER            :: I, IAB, ICOPY, IFIRST, ILAST, INDTEI, IP, IQM, I1, I2, J, JAB, JCOPY, JFIRST, &
                         JJ, JLAST, JP, JQM, KA, KFI, KFJ, L, M, NI, NJ, P, PIINDTEI
   REAL ( KIND = DP ) :: PAPOP, PBPOP, PTPOP

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)     :: FDIAG
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:,:) :: PDIAG

   !----------------------------------------------------------------------------
   ! . Calculate the one-centre part of the Fock matrix.
   !----------------------------------------------------------------------------
   ! . Initialize the Fock matrix.
   FCKALP = 0.0_DP
   FCKBET = 0.0_DP

   ! . Allocate the temporary arrays.
   ALLOCATE ( FDIAG(10*NATOMSQM), PDIAG(1:4,1:4,1:NATOMSQM) )

   ! . Set the Fock matrix counter.
   IF ( NPIBEADS == 1 ) THEN
      P = 0
   ELSE
      P = NCALC
   END IF

   ! . Compute the diagonal terms for the Fock matrix involving two-electron
   ! . integrals on the same atom (these are stored and are not calculated).
   DO IQM = 1,NATOMSQM

      ! . Get some information for the atom.
      ICOPY  = QMATOM(IQM)%COPY
      IFIRST = QMATOM(IQM)%BFIRST
      ILAST  = QMATOM(IQM)%BLAST
      NI     = QMATOM(IQM)%NUMBER

      ! . Skip atoms of the wrong index.
      IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

      ! . Calculate the density matrix contributions.
      PAPOP = 0.0_DP ; PBPOP = 0.0_DP ; PTPOP = 0.0_DP

      ! . Hydrogen.
      IF ( NATORB(NI) <= 1 ) THEN
         CONTINUE
      ! . SP elements.
      ELSE
         PAPOP = DENALP((ILAST*(ILAST+1))/2) + DENALP(((ILAST-1)*(ILAST))/2) + DENALP(((ILAST-2)*(ILAST-1))/2)
         PBPOP = DENBET((ILAST*(ILAST+1))/2) + DENBET(((ILAST-1)*(ILAST))/2) + DENBET(((ILAST-2)*(ILAST-1))/2)
         PTPOP = DENTOT((ILAST*(ILAST+1))/2) + DENTOT(((ILAST-1)*(ILAST))/2) + DENTOT(((ILAST-2)*(ILAST-1))/2)
      END IF

      ! . ss diagonal terms.
      KA = (IFIRST*(IFIRST+1))/2
      FCKALP(KA) = FCKALP(KA) + DENBET(KA) * GSS(NI) + PTPOP * GSP(NI) - PAPOP * HSP(NI)
      FCKBET(KA) = FCKBET(KA) + DENALP(KA) * GSS(NI) + PTPOP * GSP(NI) - PBPOP * HSP(NI)
      IF ( NI < 3 ) CYCLE

      L = KA
      DO J = (IFIRST+1),ILAST
         M = L+IFIRST
         L = L+J
         ! . pp diagonal terms.
         FCKALP(L) = FCKALP(L) + DENTOT(KA) * GSP(NI) - DENALP(KA) * HSP(NI) + &
                                 DENBET(L)  * GPP(NI) + ( PTPOP - DENTOT(L) ) * GP2(NI) - &
                                 0.5_DP * ( PAPOP - DENALP(L) ) * ( GPP(NI) - GP2(NI) )
         FCKBET(L) = FCKBET(L) + DENTOT(KA) * GSP(NI) - DENBET(KA) * HSP(NI) + &
                                 DENALP(L)  * GPP(NI) + ( PTPOP - DENTOT(L) ) * GP2(NI) - &
                                 0.5_DP * ( PBPOP - DENBET(L) ) * ( GPP(NI) - GP2(NI) )
	 ! . sp terms.
         FCKALP(M) = FCKALP(M) + ( 2.0_DP * DENTOT(M) * HSP(NI) ) - DENALP(M) * ( HSP(NI)+GSP(NI) )
         FCKBET(M) = FCKBET(M) + ( 2.0_DP * DENTOT(M) * HSP(NI) ) - DENBET(M) * ( HSP(NI)+GSP(NI) )
      END DO

      ! . pp* terms.
      DO J = (IFIRST+1),(ILAST-1)
         DO L = (J+1),ILAST
            M = (L*(L-1))/2+J
            FCKALP(M) = FCKALP(M) + DENTOT(M) * (GPP(NI)-GP2(NI)) - 0.5_DP * DENALP(M) * (GPP(NI)+GP2(NI))
            FCKBET(M) = FCKBET(M) + DENTOT(M) * (GPP(NI)-GP2(NI)) - 0.5_DP * DENBET(M) * (GPP(NI)+GP2(NI))
         END DO
      END DO
   END DO

   !----------------------------------------------------------------------------
   ! . Calculate the two-centre part of the Fock matrix.
   !----------------------------------------------------------------------------
   ! . Initialize some counters.
   INDTEI   = 1
   KFI      = 1
   PIINDTEI = 1
   FDIAG    = 0.0_DP
   PDIAG    = 0.0_DP

   ! . Outer loop over QM atoms.
   DO IQM = 1,NATOMSQM

      ! . Get some information for the atom.
      ICOPY  = QMATOM(IQM)%COPY
      IFIRST = QMATOM(IQM)%BFIRST
      ILAST  = QMATOM(IQM)%BLAST
      IAB    = ILAST - IFIRST + 1
      NI     = QMATOM(IQM)%NUMBER

      ! . Skip atoms of the wrong index.
      IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

      IF ( IAB == 1 ) THEN
         PDIAG(1,1,IQM) = DENTOT((IFIRST*(IFIRST+1))/2)
      ELSE
         I1 = ((IFIRST-1)*(IFIRST-2))/2
         IP = 0
         DO I = IFIRST,ILAST
            IP = IP + 1
            I1 = I1 + I - 1
            I2 = I1 + IFIRST - 1
            JP = 0
            DO J = IFIRST,I
               JP = JP + 1
               I2 = I2 + 1
               PDIAG(JP,IP,IQM) = DENTOT(I2)
               PDIAG(IP,JP,IQM) = DENTOT(I2)
            END DO
         END DO
      END IF
      KFJ = 1

      ! . Inner loop over QM atoms.
      DO JQM = 1,(IQM-1)

         ! . Get some information for the atom.
         JCOPY  = QMATOM(JQM)%COPY
         JFIRST = QMATOM(JQM)%BFIRST
         JLAST  = QMATOM(JQM)%BLAST
         JAB    = JLAST - JFIRST + 1
         NJ     = QMATOM(JQM)%NUMBER

         ! . Skip atoms of the wrong index.
         IF ( ( JCOPY /= 0 ) .AND. ( JCOPY /= P ) ) CYCLE

         ! . Both atoms are non-PI atoms.
         IF ( ( ICOPY == 0 ) .AND. ( JCOPY == 0 ) ) THEN
            CALL SEFOCK ( FDIAG(KFI:), FDIAG(KFJ:),   SETEI(  INDTEI:),     INDTEI, PDIAG(:,:,IQM), PDIAG(:,:,JQM) )
         ! . One or both atoms are PI atoms.
         ELSE
            CALL SEFOCK ( FDIAG(KFI:), FDIAG(KFJ:), PISETEI(PIINDTEI:,P), PIINDTEI, PDIAG(:,:,IQM), PDIAG(:,:,JQM) )
         END IF

         ! . Increment the counter for FDIAG.
         KFJ = KFJ + (JAB*(JAB+1))/2

      END DO

      ! . Increment the counter for FDIAG.
      KFI = KFI + (IAB*(IAB+1))/2

   END DO

   ! . Add in the diagonal parts of the Fock matrix.
   JJ = 0
   DO IQM = 1,NATOMSQM

      ! . Get some information for the atom.
      ICOPY  = QMATOM(IQM)%COPY
      IFIRST = QMATOM(IQM)%BFIRST
      ILAST  = QMATOM(IQM)%BLAST
      I1     = ((IFIRST-1)*(IFIRST-2))/2
      NI     = QMATOM(IQM)%NUMBER

      ! . Skip atoms of the wrong index.
      IF ( ( ICOPY /= 0 ) .AND. ( ICOPY /= P ) ) CYCLE

      ! . Loop over the atom basis functions.
      DO I = IFIRST,ILAST
         I1 = I1 + I - 1
         I2 = I1 + IFIRST - 1
         DO J = IFIRST,I
            I2 = I2 + 1
            JJ = JJ + 1
            FCKALP(I2) = FCKALP(I2) + FDIAG(JJ)
            FCKBET(I2) = FCKBET(I2) + FDIAG(JJ)
         END DO
      END DO
   END DO

   ! . Deallocate space.
   DEALLOCATE ( FDIAG, PDIAG )

   !============================================================================
   CONTAINS
   !============================================================================

      !----------------------------------------------------
      SUBROUTINE SEFOCK ( FII, FJJ, SETEI, KR, P2II, P2JJ )
      !----------------------------------------------------

      !
      !     SEFCC1 constructs the two-centre two-electron part of the Fock
      !     matrix given :
      !
      !     F - the partial packed Fock matrix, FII - the atom ii part of the
      !     Fock matrix, FJJ - the atom jj part of the Fock matrix, W - the
      !     two-centre two-electron integrals and P - the density matrix.
      !
      !     Atom indices are IFIRST to ILAST for atom I and JFIRST to JLAST
      !     for atom J.
      !

      ! . Scalar argument declarations.
      INTEGER, INTENT(INOUT) :: KR

      ! . Array argument declarations.
      REAL ( KIND = DP ), DIMENSION(:),     INTENT(INOUT) :: FII,  FJJ
      REAL ( KIND = DP ), DIMENSION(1:16),  INTENT(IN)    :: P2II, P2JJ
      REAL ( KIND = DP ), DIMENSION(:),     INTENT(IN)    :: SETEI

      ! . Local scalars.
      INTEGER            :: I, IAC, IAC2, II, IJ, IJAC, J, JAC, JAC2, JJ, K, KK, KL, NTYPE
      REAL ( KIND = DP ) :: SUM, SUMA, SUMB

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:16) :: PAIJ, PBIJ

      ! . Establish NTYPE. If 1 both atoms in the diatomic have 4 orbitals,
      ! . if 2 one atom has 4 orbitals and the other 1 and if 3 both atoms
      ! . have single orbitals.
      IAC = ILAST - IFIRST + 1
      JAC = JLAST - JFIRST + 1
      IF ( ( IAC + JAC ) == 8 ) THEN
         NTYPE = 1
      ELSE IF ( ( IAC + JAC ) == 5 ) THEN
         NTYPE = 2
      ELSE
         IJ = (IFIRST*(IFIRST-1))/2 + JFIRST
         FCKALP(IJ) = FCKALP(IJ) - ( SETEI(1) * DENALP(IJ) )
         FCKBET(IJ) = FCKBET(IJ) - ( SETEI(1) * DENBET(IJ) )
         FII(1)     = FII(1)     + ( SETEI(1) * P2JJ(1)    )
         FJJ(1)     = FJJ(1)     + ( SETEI(1) * P2II(1)    )
         KR = KR + 1
         RETURN
      END IF

      ! . Unpack off-diagonal spin density matrix elements.
      KK = ((IFIRST-1)*(IFIRST-2))/2
      DO I = IFIRST,ILAST
         KK = KK + I - 1
         JJ = ( I - IFIRST ) * JAC
         K  = KK + JFIRST - 1
         IF ( JAC == 1 ) THEN
            PAIJ(JJ+1) = DENALP(K+1)
            PBIJ(JJ+1) = DENBET(K+1)
         ELSE
            DO J = JFIRST,JLAST
               JJ = JJ + 1
               K  = K  + 1
               PAIJ(JJ) = DENALP(K)
               PBIJ(JJ) = DENBET(K)
            END DO
         ENDIF
      END DO

      ! . Construct the off-diagonal exchange contributions.
      IJAC = IAC * JAC
      II   = ((IFIRST-1)*(IFIRST-2))/2
      DO I = 1,IAC
         II = II + IFIRST + I - 2
         IJ = II + JFIRST - 1
         DO J = 1,JAC
            IJ   = IJ + 1
            SUMA = 0.0_DP
            SUMB = 0.0_DP
            DO KL = 1,IJAC
               SUMA = SUMA + ( SETEI( KIND2(KL,I,J,NTYPE) ) * PAIJ(KL) )
               SUMB = SUMB + ( SETEI( KIND2(KL,I,J,NTYPE) ) * PBIJ(KL) )
            END DO
            FCKALP(IJ) = FCKALP(IJ) - SUMA
            FCKBET(IJ) = FCKBET(IJ) - SUMB
         END DO
      END DO

      ! . Construct the on-diagonal Coulombic contributions.
      ! . Atom I affected by atom J.
      JAC2 = JAC * JAC
      JJ   = 0
      DO I = 1,IAC
         DO J = 1,I
            JJ = JJ + 1
            IF ( JAC2 == 1 ) THEN
               FII(JJ) = FII(JJ) + ( SETEI( JIND2(1,I,J,2) ) * P2JJ(1) )
            ELSE
               SUM = 0.0_DP
               DO KL = 1,JAC2
                  SUM = SUM + ( SETEI( JIND2(KL,I,J,NTYPE) ) * P2JJ(KL) )
               END DO
               FII(JJ) = FII(JJ) + SUM
            END IF
         END DO
      END DO

      ! . Atom J affected by atom I.
      IAC2 = IAC * IAC
      JJ   = 0
      DO I = 1,JAC
         DO J = 1,I
            JJ = JJ + 1
            IF ( IAC2 == 1 ) THEN
               FJJ(JJ) = FJJ(JJ) + ( SETEI( JIND3(I,J,1,2) ) * P2II(1) )
            ELSE
               SUM = 0.0_DP
               DO KL = 1,IAC2
                  SUM = SUM + ( SETEI( JIND3(I,J,KL,NTYPE) ) * P2II(KL) )
               END DO
               FJJ(JJ) = FJJ(JJ) + SUM
            END IF
         END DO
      END DO

      ! . Increment KR counter to line up W matrix for the next pair of atoms.
      KR = KR + ( ( IAC*(IAC+1) * JAC*(JAC+1) ) / 4 )

      END SUBROUTINE SEFOCK

   END SUBROUTINE FOCK_MATRIX_UNRESTRICTED

END MODULE MOPAC_FOCK_MATRIX
