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
!                             The MOPAC Data Module
!===============================================================================
!
! . System options:
!
!   TOTCHG                   the total charge on the system.
!   MULTIP                   the multiplicity of the system.
!   NALPHA                   the number of occupied alpha orbitals.
!   NBETA                    the number of occupied alpha orbitals.
!   NELEC                    the number of electrons in the system.
!
!   NBASIS                   the number of basis functions.
!   NBASTR                   ( NBASIS * ( NBASIS + 1 ) ) / 2.
!
!   QUSEBAQM                 the BA(QM)/MM non-bonding interaction flag.
!
! . Approximation options:
!
!   HAMILTONIAN              the semiempirical Hamiltonian to use.
!
! . AM1/MNDO/PDDG/PM3/RM1 data:
!
!   ATHEAT                   the energy baseline (kJ mol^-1).
!   N2ELEC                   the number of two-electron integrals.
!
! . Atom data:
!
!   QMATOM                   data about the QM atoms.
!
! . Integer arrays:
!
!   BFINDEX                  the pair index array (IA(I) = I(I-1)/2).
!
! . Real arrays:
!
!   BETA                     the parameters for the two-centre, one-electron
!                            resonance integrals in a MNDO/AM1 calculation.
!   DENMAT                   the density matrix.
!   DENMATA                  the alpha density matrix.
!   DENMATB                  the beta density matrix.
!   HCORE                    the one-electron matrix.
!   OVERLAP                  the overlap matrix.
!   SETEI                    the two-electron integrals.
!   USPD                     the kinetic energy integrals.
!
! . Indexing arrays for the semi-empirical two-electron integrals:
!
!   JIND2, JIND3,            arrays used to index the semiempirical two-electron
!   KIND2                    integrals. J holds the gather indices for the
!                            Coulomb integrals and K holds those for the
!                            exchange integrals. The last index refers to the
!                            number of AOs on the atoms - 1 is for systems with
!                            AOs on each atom and 2 is for systems with 1 and 4
!                            or 4 and 1 AOs on each atom.
!
! . Functions:
!
!   MOPAC_DATA_ATOM_POSITION return the position of a quantum atom.
!   MOPAC_DATA_LA_GRADIENT   return the modified gradient vector for a link atom.
!
! . Subroutines:
!
!   MOPAC_DATA_INITIALIZE    initialize the MOPAC data structure.
!
! . Path integral scalars:
!
!   NPIATOMS                 the number of path integral atoms.
!   NPIBASIS                 the number of basis functions for a single
!                            PI polymer.
!   NPIBEADS                 the number of polymers (always >= 1).
!   NPIHELE                  the number of PI atom one-electron terms/polymer.
!   NPI2ELEC                 the number of PI atom TEIs/polymer.
!
!   PIFC                     the PI harmonic energy force constant prefactor.
!
! . Path integral arrays:
!
!   HSTORE                   the invariant part of the one-electron matrix.
!   PIHCORE                  the PI atom one-electron matrix terms.
!   PIHDIAG                  the non-PI atom diagonal one-electron matrix elements.
!   PILIST                   the PI atoms for each polymer.
!   PISETEI                  the PI atom two-electron integrals.
!
!===============================================================================
MODULE MOPAC_DATA

! . Module declarations.
USE DEFINITIONS,    ONLY : DP
USE LINEAR_ALGEBRA, ONLY : NORMALIZE

USE ATOMS,          ONLY : NATOMSQM

IMPLICIT NONE
PUBLIC
#ifndef PGPC
SAVE
#endif

! . Type definitions.
TYPE QMATOM_TYPE
   INTEGER            :: ATOM, BFIRST, BLAST, COPY, NUMBER, PARTNER
   LOGICAL            :: QBOUNDARY
   REAL ( KIND = DP ) :: LENGTH
END TYPE QMATOM_TYPE

! . QM atom data.
TYPE(QMATOM_TYPE), ALLOCATABLE, DIMENSION(:) :: QMATOM

! . General scalar data.
CHARACTER ( LEN = 4 ) :: HAMILTONIAN
INTEGER               :: MULTIP   = 1, NALPHA = 0, NBETA = 0, NELEC = 0, N2ELEC = 0, TOTCHG = 0
LOGICAL               :: QQMMMORIGINAL = .TRUE., QUSEBAQM = .FALSE., REALDENSITY = .FALSE.
REAL ( KIND = DP )    :: ATHEAT   = 0.0_DP

! . Scalar basis function data.
INTEGER :: NBASIS = 0, NBASTR = 0

! . General array data.
INTEGER,            ALLOCATABLE, DIMENSION(:)   :: BFINDEX
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: BETA, HCORE, OVERLAP, SETEI, USPD
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: DENMAT, DENMATA, DENMATB

! . Index arrays.
INTEGER, DIMENSION(1:16,1:4,1:4,1:2) :: JIND2
INTEGER, DIMENSION(1:4,1:4,1:16,1:2) :: JIND3
INTEGER, DIMENSION(1:16,1:4,1:4,1:2) :: KIND2

! . Path integral scalars.
INTEGER            :: NPIATOMS = 0, NPIBASIS = 0, NPIBEADS = 1, NPIHELE = 0, NPI2ELEC = 0
REAL ( KIND = DP ) :: PIFC = 0.0_DP

! . Path integral arrays.
INTEGER,            ALLOCATABLE, DIMENSION(:,:)   :: PILIST
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)     :: HSTORE
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:)   :: PIHCORE, PISETEI
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:,:) :: PIHDIAG

!===============================================================================
CONTAINS
!===============================================================================

   !-------------------------------
   SUBROUTINE MOPAC_DATA_INITIALIZE
   !-------------------------------

   ! . Scalar data.
   HAMILTONIAN   = "AM1 "
   MULTIP        = 1
   NBASIS        = 0
   NBASTR        = 0
   NALPHA        = 0
   NBETA         = 0
   NELEC         = 0
   N2ELEC        = 0
   TOTCHG        = 0
   QQMMMORIGINAL = .TRUE.
   QUSEBAQM      = .FALSE.
   REALDENSITY   = .FALSE.
   ATHEAT        = 0.0_DP

   ! . Array data.
   IF ( ALLOCATED ( BFINDEX ) ) DEALLOCATE ( BFINDEX )
   IF ( ALLOCATED ( BETA    ) ) DEALLOCATE ( BETA    )
   IF ( ALLOCATED ( DENMAT  ) ) DEALLOCATE ( DENMAT  )
   IF ( ALLOCATED ( DENMATA ) ) DEALLOCATE ( DENMATA )
   IF ( ALLOCATED ( DENMATB ) ) DEALLOCATE ( DENMATB )
   IF ( ALLOCATED ( HCORE   ) ) DEALLOCATE ( HCORE   )
   IF ( ALLOCATED ( OVERLAP ) ) DEALLOCATE ( OVERLAP )
   IF ( ALLOCATED ( SETEI   ) ) DEALLOCATE ( SETEI   )
   IF ( ALLOCATED ( USPD    ) ) DEALLOCATE ( USPD    )
   IF ( ALLOCATED ( QMATOM  ) ) DEALLOCATE ( QMATOM  )

   ! . Path integral scalars.
   NPIATOMS = 0
   NPIBASIS = 0
   NPIBEADS = 1
   NPIHELE  = 0
   NPI2ELEC = 0
   PIFC     = 0.0_DP

   ! . Path integral arrays.
   IF ( ALLOCATED ( HSTORE  ) ) DEALLOCATE ( HSTORE  )
   IF ( ALLOCATED ( PILIST  ) ) DEALLOCATE ( PILIST  )
   IF ( ALLOCATED ( PIHCORE ) ) DEALLOCATE ( PIHCORE )
   IF ( ALLOCATED ( PIHDIAG ) ) DEALLOCATE ( PIHDIAG )
   IF ( ALLOCATED ( PISETEI ) ) DEALLOCATE ( PISETEI )

   END SUBROUTINE MOPAC_DATA_INITIALIZE

   !-------------------------------------------------------
   FUNCTION MOPAC_DATA_ATOM_POSITION ( QATOM, COORDINATES )
   !-------------------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3) :: MOPAC_DATA_ATOM_POSITION

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: QATOM

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . The QM atom index is valid.
   IF ( ( QATOM >= 1 ) .AND. ( QATOM <= NATOMSQM ) ) THEN

      ! . The atom is a boundary atom.
      IF ( QMATOM(QATOM)%QBOUNDARY ) THEN

         ! . Calculate the link atom position.
         MOPAC_DATA_ATOM_POSITION = COORDINATES(1:3,QMATOM(QATOM)%PARTNER) + QMATOM(QATOM)%LENGTH *  &
	                                        NORMALIZE ( COORDINATES(1:3,QMATOM(QATOM)%ATOM) -  &
						            COORDINATES(1:3,QMATOM(QATOM)%PARTNER) )

      ! . The atom is a normal atom.
      ELSE
         MOPAC_DATA_ATOM_POSITION = COORDINATES(1:3,QMATOM(QATOM)%ATOM)
      END IF

   ! . The QM atom index is invalid.
   ELSE

      ! . Initialize MOPAC_DATA_ATOM_POSITION.
      MOPAC_DATA_ATOM_POSITION = HUGE ( 0.0_DP )

   END IF

   END FUNCTION MOPAC_DATA_ATOM_POSITION

   !---------------------------------------------------------------
   FUNCTION MOPAC_DATA_LA_GRADIENT ( QATOM, COORDINATES, GRADIENT )
   !---------------------------------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3) :: MOPAC_DATA_LA_GRADIENT

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: QATOM

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3), INTENT(IN) :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN) :: COORDINATES

   ! . Local scalars.
   REAL ( KIND = DP ) :: RMQ

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DRMQ

   ! . The QM atom index is valid.
   IF ( ( QATOM >= 1 ) .AND. ( QATOM <= NATOMSQM ) ) THEN

      ! . The atom is a boundary atom.
      IF ( QMATOM(QATOM)%QBOUNDARY ) THEN

         ! . Get the interatomic vector.
	 DRMQ = COORDINATES(1:3,QMATOM(QATOM)%ATOM) - COORDINATES(1:3,QMATOM(QATOM)%PARTNER)

         ! . Find the vector's size.
	 RMQ  = SQRT ( DOT_PRODUCT ( DRMQ, DRMQ ) )

         ! . Normalize the vector.
	 DRMQ = DRMQ / RMQ

         ! . Calculate the modified gradient vector.
         MOPAC_DATA_LA_GRADIENT = ( QMATOM(QATOM)%LENGTH / RMQ ) * ( GRADIENT - DOT_PRODUCT ( DRMQ, GRADIENT ) * DRMQ )

      END IF
   END IF

   END FUNCTION MOPAC_DATA_LA_GRADIENT

END MODULE MOPAC_DATA

