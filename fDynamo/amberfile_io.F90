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
!                           The AMBER FILE I/O Module
!===============================================================================
!
! . Subroutines:
!
!   AMBERCRD_WRITE       Write out coordinates in AMBER CRD format.
!   AMBERTOP_WRITE       Write out system information in AMBER TOP format.
!
!===============================================================================
!------------------
MODULE AMBERFILE_IO
!------------------

! . Module declarations.
USE CONSTANTS,    ONLY : KCAL_TO_KJ
USE DEFINITIONS,  ONLY : DP
USE FILES,        ONLY : NEXT_UNIT
USE IO_UNITS,     ONLY : OUTPUT
USE PRINTING,     ONLY : PRINT_ERROR, PRINT_PARAGRAPH

USE ATOMS
USE CONNECTIVITY
USE MM_TERMS
USE SEQUENCE
USE SYMMETRY

IMPLICIT NONE
PRIVATE
PUBLIC :: AMBERCRD_WRITE, AMBERTOP_WRITE

!==============================================================================
CONTAINS
!==============================================================================

    !--------------------------------------------------------------------------
    SUBROUTINE AMBERCRD_WRITE ( FILE, DATA, INCLUDESYMMETRY, SELECTION, TITLE )
    !--------------------------------------------------------------------------

    ! . Optional scalar arguments.
    CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE, TITLE
    LOGICAL,               INTENT(IN), OPTIONAL :: INCLUDESYMMETRY

    ! . Optional array arguments.
    LOGICAL,            DIMENSION(1:NATOMS),     INTENT(IN), OPTIONAL :: SELECTION
    REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(IN), OPTIONAL :: DATA

    ! . Local scalars.
    CHARACTER ( LEN = 16 ) :: TAG
    INTEGER                :: I, IOSTAT, J, LTAG, N, NATM, UNIT
    LOGICAL                :: QINCLUDESYMMETRY

    ! . Local arrays.
    REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: ATMDAT

    ! . Find the size of the data.
    IF ( PRESENT ( SELECTION ) ) THEN
       NATM = COUNT ( SELECTION )
    ELSE
       NATM = NATOMS
    END IF
    ALLOCATE ( ATMDAT(1:3*NATM) )

    ! . Fill ATMDAT.
    IF ( PRESENT ( DATA ) ) THEN
        TAG = "Data"
        IF ( PRESENT ( SELECTION ) ) THEN
            N = 0
            DO I = 1,NATOMS
                IF ( SELECTION(I) ) THEN
                    DO J = 1,3
                       N = N + 1
                       ATMDAT(N) = DATA(J,I)
                    END DO
                END IF
            END DO
        ELSE
            N = 0
            DO I = 1,NATOMS
                DO J = 1,3
                   N = N + 1
                   ATMDAT(N) = DATA(J,I)
                END DO
            END DO
        END IF
    ! . By default write out the contents of ATMCRD.
    ELSE
        TAG = "Coordinates"
        IF ( PRESENT ( SELECTION ) ) THEN
            N = 0
            DO I = 1,NATOMS
                IF ( SELECTION(I) ) THEN
                    DO J = 1,3
                       N = N + 1
                       ATMDAT(N) = ATMCRD(J,I)
                    END DO
                END IF
            END DO
        ELSE
            N = 0
            DO I = 1,NATOMS
                DO J = 1,3
                   N = N + 1
                   ATMDAT(N) = ATMCRD(J,I)
                END DO
            END DO
        END IF
    END IF

    ! . Get the length of the tag.
    LTAG = LEN_TRIM ( TAG )

    ! . Check to see if the FILE argument is present.
    IF ( PRESENT ( FILE ) ) THEN

       ! . Get the next unit number.
       UNIT = NEXT_UNIT ( )

       ! . Open the file.
       OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = IOSTAT )

       ! . If there has been an error return.
       IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "AMBERCRD_WRITE", "I/O Error.", IOSTAT )

       ! . Do some printing.
       CALL PRINT_PARAGRAPH ( TEXT = TAG(1:LTAG)//" written to " // FILE )

    ! . The unit for writing is the output stream.
    ELSE

       ! . Assign the unit number.
       UNIT = OUTPUT

       ! . Do some printing.
       CALL PRINT_PARAGRAPH ( TEXT = TAG(1:LTAG)//" written to the output stream." )
       CALL PRINT_BLANKLINE
       CALL PRINT_NOFORMAT_START

       ! . Write out the separator.
       WRITE ( UNIT, "('!',69('='))" )

    END IF

    ! . Write out a title.
    IF ( PRESENT ( TITLE ) ) THEN
       WRITE ( UNIT, "(A)" ) TITLE
    ELSE
       WRITE ( UNIT, "(A)" ) "AMBER Coordinate File."
    END IF

    ! . Write out the number of atoms to the file.
    WRITE ( UNIT, "(I6)" ) NATM

    ! . Write out coordinates for the selected atoms.
    WRITE ( UNIT, "(6F12.7)" ) ATMDAT
    DEALLOCATE ( ATMDAT )

    ! . Get symmetry options.
    IF ( PRESENT ( INCLUDESYMMETRY ) ) THEN
        QINCLUDESYMMETRY = INCLUDESYMMETRY
    ELSE
        QINCLUDESYMMETRY = .TRUE.
    END IF

    ! . Write out symmetry.
    IF ( QBOX .AND. QINCLUDESYMMETRY ) THEN
        WRITE ( UNIT, "(3F12.7)" ) BOXL(1), BOXL(2), BOXL(3)
    END IF

    ! . Finish off.
    IF ( UNIT == OUTPUT ) THEN

       ! . Write out the terminator.
       WRITE ( UNIT, "('!',69('='))" )

       ! . Finish noformatting.
       CALL PRINT_NOFORMAT_STOP

    ! . Close the file.
    ELSE
       CLOSE ( UNIT )
    END IF

    END SUBROUTINE AMBERCRD_WRITE

    !----------------------------------------------------------------------------
    SUBROUTINE AMBERTOP_WRITE ( FILE, INCLUDESYMMETRY, LASTSOLUTERESIDUE, TITLE )
    !----------------------------------------------------------------------------

    ! . The 1-4 scale factors need to be set manually (they are usually 0.5 for both electrostatics and Lennard-Jones).

    ! . Optional scalar arguments.
    CHARACTER ( LEN = * ), INTENT(IN)           :: FILE
    CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: TITLE
    INTEGER,               INTENT(IN), OPTIONAL :: LASTSOLUTERESIDUE
    LOGICAL,               INTENT(IN), OPTIONAL :: INCLUDESYMMETRY

    ! . Local parameters.
    REAL(DP), PARAMETER :: DIFFTOL = 1.0E-6_DP, ELECTRONTOKCAL = 18.2223_DP

    ! . Local scalars.
    INTEGER :: I, IOSTAT, UNIT
    LOGICAL :: QINCLUDESYMMETRY

    ! . AMBER pointers.
    INTEGER :: NLJ, NBONH, MBONA, NTHETH, MTHETA, NPHIH, MPHIA, NHPARM, NPARM,        &
               NNB, NBONA, NTHETA, NPHIA, NUMBND, NUMANG, NPTRA, NATYP, NPHB,         &
               IFPERT, NBPER, NGPER, NDPER, MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP, &
               NUMEXT

    ! .Other AMBER scalars.
    INTEGER :: IPTRES, NSPM, NSPSOL

    ! . Bond, angle, dihedral/improper arrays.
    INTEGER,  ALLOCATABLE, DIMENSION(:)   :: AHINDICES, AQINDICES, BHINDICES, BQINDICES, NSP, THINDICES, TQINDICES
    REAL(DP), ALLOCATABLE, DIMENSION(:)   :: AEQ, AFC, BEQ, BFC, TFC, TPERIOD, TPHASE

    ! . Exclusion arrays.
    INTEGER,  ALLOCATABLE, DIMENSION(:)   :: EXCLI, EXCLJ
    
    ! . LJ arrays.
    INTEGER,  ALLOCATABLE, DIMENSION(:)   :: LJINDICES
    INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: LJTABLE
    REAL(DP), ALLOCATABLE, DIMENSION(:)   :: LJACOEF, LJBCOEF

    !---------------------------------------------------------------------------
    ! . Initialization.
    !---------------------------------------------------------------------------
    ! . No QM atoms.
    IF ( ( NATOMS <= 0 ) .OR. ( NATOMSQM > 0 ) ) RETURN

    ! . Get the next unit number.
    UNIT = NEXT_UNIT ( )

    ! . Open the file.
    OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = IOSTAT )

    ! . If there has been an error return.
    IF ( IOSTAT /= 0 ) CALL PRINT_ERROR ( "AMBERTOP_WRITE", "I/O Error.", IOSTAT )

    ! . Do some printing.
    CALL PRINT_PARAGRAPH ( TEXT = "AMBER system written to " // FILE )

    ! . Get symmetry options.
    IF ( PRESENT ( LASTSOLUTERESIDUE ) ) THEN
        IPTRES = MIN ( MAX ( LASTSOLUTERESIDUE, 0 ), NRESID )
    ELSE
        IPTRES = NRESID
    END IF
    IF ( PRESENT ( INCLUDESYMMETRY ) ) THEN
        QINCLUDESYMMETRY = INCLUDESYMMETRY
    ELSE
        QINCLUDESYMMETRY = .TRUE.
    END IF

    !---------------------------------------------------------------------------
    ! . Get all data.
    !---------------------------------------------------------------------------
    ! . Sequence information.
    CALL SEQUENCE_UNKNOWN

    ! . Number of atom types.
    CALL FINDATOMTYPES

    ! . Find bond, angle and dihedral/improper types.
    CALL FINDBONDTYPES
    CALL FINDANGLETYPES
    CALL FINDTORSIONTYPES

    ! . Find the exclusions.
    CALL FINDEXCLUSIONS

    ! . Find the Lennard-Jones types.
    CALL FINDLJTYPES

    ! . Get the size of the largest residue.
    NMXRS = 0
    DO I = 1,NRESID
        NMXRS = MAX ( NMXRS, RESIND(I+1) - RESIND(I) )
    END DO

    ! . Zero counters.
    NHPARM = 0
    NPARM  = 0
    NPHB   = 0
    IFPERT = 0
    NBPER  = 0
    NGPER  = 0
    NDPER  = 0
    MBPER  = 0
    MGPER  = 0
    MDPER  = 0
    IFBOX  = 0
    IFCAP  = 0
    NUMEXT = 0

    ! . Find symmetry information.
    CALL FINDSYMMETRYINFORMATION

    !---------------------------------------------------------------------------
    ! . Writing - header.
    !---------------------------------------------------------------------------
    ! . Version.
    WRITE ( UNIT, "(A)" ) "%VERSION  VERSION_STAMP = V0001.000  DATE = __/__/__  __:__:__" 

    ! . Title.
    WRITE ( UNIT, "(A)" ) "%FLAG TITLE"
    WRITE ( UNIT, "(A)" ) "%FORMAT(20A4)"
    IF ( PRESENT ( TITLE ) ) THEN
       WRITE ( UNIT, "(A)" ) TITLE
    ELSE
       WRITE ( UNIT, "(A)" ) "AMBER Topology File."
    END IF

    ! . Pointers.
    WRITE ( UNIT, "(A)" ) "%FLAG POINTERS"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) NATOMS, NLJ, NBONH, MBONA, NTHETH, MTHETA, NPHIH, MPHIA, NHPARM, NPARM, &
                             NNB, NRESID, NBONA, NTHETA, NPHIA, NUMBND, NUMANG, NPTRA, NATYP, NPHB,  &
                             IFPERT, NBPER, NGPER, NDPER, MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP,  &
                             NUMEXT

    !---------------------------------------------------------------------------
    ! . Writing - sections.
    !---------------------------------------------------------------------------
    WRITE ( UNIT, "(A)" ) "%FLAG ATOM_NAME"
    WRITE ( UNIT, "(A)" ) "%FORMAT(20A4)"
    WRITE ( UNIT, "(20A4)" ) ATMNAM

    WRITE ( UNIT, "(A)" ) "%FLAG CHARGE"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) ATMCHG * ELECTRONTOKCAL

    WRITE ( UNIT, "(A)" ) "%FLAG MASS"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) ATMMAS

    WRITE ( UNIT, "(A)" ) "%FLAG ATOM_TYPE_INDEX"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) LJINDICES(1:NATOMS)

    WRITE ( UNIT, "(A)" ) "%FLAG NUMBER_EXCLUDED_ATOMS"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) EXCLI(1:NATOMS)

    WRITE ( UNIT, "(A)" ) "%FLAG NONBONDED_PARM_INDEX"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) LJTABLE

    WRITE ( UNIT, "(A)" ) "%FLAG RESIDUE_LABEL"
    WRITE ( UNIT, "(A)" ) "%FORMAT(20A4)"
    WRITE ( UNIT, "(20A4)" ) RESNAM(1:NRESID)

    WRITE ( UNIT, "(A)" ) "%FLAG RESIDUE_POINTER"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) RESIND(1:NRESID) + 1

    WRITE ( UNIT, "(A)" ) "%FLAG BOND_FORCE_CONSTANT"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) BFC(1:NUMBND)

    WRITE ( UNIT, "(A)" ) "%FLAG BOND_EQUIL_VALUE"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) BEQ(1:NUMBND)

    WRITE ( UNIT, "(A)" ) "%FLAG ANGLE_FORCE_CONSTANT"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) AFC(1:NUMANG)

    WRITE ( UNIT, "(A)" ) "%FLAG ANGLE_EQUIL_VALUE"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) AEQ(1:NUMANG)

    WRITE ( UNIT, "(A)" ) "%FLAG DIHEDRAL_FORCE_CONSTANT"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) TFC(1:NPTRA)

    WRITE ( UNIT, "(A)" ) "%FLAG DIHEDRAL_PERIODICITY"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) TPERIOD(1:NPTRA)

    WRITE ( UNIT, "(A)" ) "%FLAG DIHEDRAL_PHASE"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) TPHASE(1:NPTRA)

    WRITE ( UNIT, "(A)" ) "%FLAG SOLTY"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) ( 0.0_DP, I = 1,NATYP )

    WRITE ( UNIT, "(A)" ) "%FLAG LENNARD_JONES_ACOEF"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) LJACOEF

    WRITE ( UNIT, "(A)" ) "%FLAG LENNARD_JONES_BCOEF"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) LJBCOEF

    WRITE ( UNIT, "(A)" ) "%FLAG BONDS_INC_HYDROGEN"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) BHINDICES(1:3*NBONH)

    WRITE ( UNIT, "(A)" ) "%FLAG BONDS_WITHOUT_HYDROGEN"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) BQINDICES(1:3*MBONA)

    WRITE ( UNIT, "(A)" ) "%FLAG ANGLES_INC_HYDROGEN"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) AHINDICES(1:4*NTHETH)

    WRITE ( UNIT, "(A)" ) "%FLAG ANGLES_WITHOUT_HYDROGEN"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) AQINDICES(1:4*MTHETA)

    WRITE ( UNIT, "(A)" ) "%FLAG DIHEDRALS_INC_HYDROGEN"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) THINDICES(1:5*NPHIH)

    WRITE ( UNIT, "(A)" ) "%FLAG DIHEDRALS_WITHOUT_HYDROGEN"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) TQINDICES(1:5*MPHIA)

    WRITE ( UNIT, "(A)" ) "%FLAG EXCLUDED_ATOMS_LIST"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) EXCLJ(1:NNB)

    WRITE ( UNIT, "(A)" ) "%FLAG HBOND_ACOEF"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) ( 0.0_DP, I = 1,NPHB )

    WRITE ( UNIT, "(A)" ) "%FLAG HBOND_BCOEF"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) ( 0.0_DP, I = 1,NPHB )

    WRITE ( UNIT, "(A)" ) "%FLAG HBCUT"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) ( 0.0_DP, I = 1,NPHB )

    WRITE ( UNIT, "(A)" ) "%FLAG AMBER_ATOM_TYPE"
    WRITE ( UNIT, "(A)" ) "%FORMAT(20A4)"
    WRITE ( UNIT, "(20A4)" ) ATMTYP(1:NATOMS)

    WRITE ( UNIT, "(A)" ) "%FLAG TREE_CHAIN_CLASSIFICATION"
    WRITE ( UNIT, "(A)" ) "%FORMAT(20A4)"
    WRITE ( UNIT, "(20A4)" ) ( "    ", I = 1,NATOMS )

    WRITE ( UNIT, "(A)" ) "%FLAG JOIN_ARRAY"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) ( 0, I = 1,NATOMS )

    WRITE ( UNIT, "(A)" ) "%FLAG IROTAT"
    WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
    WRITE ( UNIT, "(10I8)" ) ( 0, I = 1,NATOMS )

    ! . Solvent information.
    IF ( IFBOX > 0 ) THEN

        WRITE ( UNIT, "(A)" ) "%FLAG SOLVENT_POINTERS"                                                      
        WRITE ( UNIT, "(A)" ) "%FORMAT(3I8)"
        WRITE ( UNIT, "(3I8)" ) IPTRES, NSPM, NSPSOL

        WRITE ( UNIT, "(A)" ) "%FLAG ATOMS_PER_MOLECULE"
        WRITE ( UNIT, "(A)" ) "%FORMAT(10I8)"
        WRITE ( UNIT, "(10I8)" ) ( NSP(I), I = 1,NSPM )

        WRITE ( UNIT, "(A)" ) "%FLAG BOX_DIMENSIONS"
        WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
        WRITE ( UNIT, "(5E16.4)" ) 90.0_DP, BOXL(1), BOXL(2), BOXL(3)

    END IF

    WRITE ( UNIT, "(A)" ) "%FLAG RADII"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) ( 0.0_DP, I = 1,NATOMS )

    WRITE ( UNIT, "(A)" ) "%FLAG SCREEN"
    WRITE ( UNIT, "(A)" ) "%FORMAT(5E16.8)"
    WRITE ( UNIT, "(5E16.8)" ) ( 0.0_DP, I = 1,NATOMS )

    !---------------------------------------------------------------------------
    ! . Termination.
    !---------------------------------------------------------------------------
    ! . Close the file.
    CLOSE ( UNIT )

    ! . Deallocate space.
    DEALLOCATE ( AEQ, AFC, AHINDICES, AQINDICES )
    DEALLOCATE ( BEQ, BFC, BHINDICES, BQINDICES )
    DEALLOCATE ( TFC, TPERIOD, TPHASE, THINDICES, TQINDICES )
    DEALLOCATE ( EXCLI, EXCLJ )
    DEALLOCATE ( LJACOEF, LJBCOEF, LJINDICES, LJTABLE )
    IF ( IFBOX > 0 ) DEALLOCATE ( NSP )

    !===========================================================================
    CONTAINS
    !===========================================================================

        !------------------------
        SUBROUTINE FINDANGLETYPES
        !------------------------

        ! . Local variables.
        LOGICAL  :: QAH
        INTEGER  :: AINDEX, I, T

        ! . Create the type and index arrays - slow version.
        ALLOCATE ( AEQ(1:NANGLES), AFC(1:NANGLES), AHINDICES(1:4*NANGLES), AQINDICES(1:4*NANGLES) )
        MTHETA = 0
        NTHETH = 0
        NUMANG = 0
        DO I = 1,NANGLES
            ! . Hydrogen-containing.
            QAH = ( ( ATMNUM(ANGLES(I)%I) == 1 ) .OR. ( ATMNUM(ANGLES(I)%J) == 1 ) .OR. ( ATMNUM(ANGLES(I)%K) == 1 ) )
            ! . Types.
            DO T = 1,NUMANG
                IF ( ( ABS ( ANGLES(I)%EQ - AEQ(T) ) < DIFFTOL ) .AND. ( ABS ( ANGLES(I)%FC - AFC(T) ) < DIFFTOL ) ) THEN
                    AINDEX = T
                    GO TO 10
                END IF
            END DO
            NUMANG = NUMANG + 1
            AEQ(NUMANG) = ANGLES(I)%EQ
            AFC(NUMANG) = ANGLES(I)%FC
            AINDEX = NUMANG
            10 CONTINUE
            ! . Index arrays.
            IF ( QAH ) THEN
                AHINDICES(4*NTHETH+1) = 3 * ( ANGLES(I)%I - 1 )
                AHINDICES(4*NTHETH+2) = 3 * ( ANGLES(I)%J - 1 )
                AHINDICES(4*NTHETH+3) = 3 * ( ANGLES(I)%K - 1 )
                AHINDICES(4*NTHETH+4) = AINDEX
                NTHETH = NTHETH + 1
            ELSE
                AQINDICES(4*MTHETA+1) = 3 * ( ANGLES(I)%I - 1 )
                AQINDICES(4*MTHETA+2) = 3 * ( ANGLES(I)%J - 1 )
                AQINDICES(4*MTHETA+3) = 3 * ( ANGLES(I)%K - 1 )
                AQINDICES(4*MTHETA+4) = AINDEX
                MTHETA = MTHETA + 1
            END IF
        END DO

        ! . Convert AFC to kcal.
        AFC(1:NUMANG) = AFC(1:NUMANG) / KCAL_TO_KJ

        ! . Set some counters.
        NTHETA = MTHETA

        END SUBROUTINE FINDANGLETYPES

        !-----------------------
        SUBROUTINE FINDATOMTYPES
        !-----------------------

        ! . Local variables.
        INTEGER :: I, T
        CHARACTER ( LEN = ATOM_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: ATYPES

        ! . Find the number of atom types - slow version.
        ALLOCATE ( ATYPES(1:NATOMS) )
        NATYP = 0
        DO I = 1,NATOMS
            DO T = 1,NATYP
                IF ( ATMTYP(I) == ATYPES(T) ) GO TO 10
             END DO
            NATYP = NATYP + 1
            ATYPES(NATYP) = ATMTYP(I)
           10 CONTINUE
        END DO
        DEALLOCATE ( ATYPES )

        END SUBROUTINE FINDATOMTYPES

        !-----------------------
        SUBROUTINE FINDBONDTYPES
        !-----------------------

        ! . Local variables.
        LOGICAL  :: QBH
        INTEGER  :: BINDEX, I, T

        ! . Create the type and index arrays - slow version.
        ALLOCATE ( BEQ(1:NBONDS), BFC(1:NBONDS), BHINDICES(1:3*NBONDS), BQINDICES(1:3*NBONDS) )
        MBONA  = 0
        NBONH  = 0
        NUMBND = 0
        DO I = 1,NBONDS
            ! . Hydrogen-containing.
            QBH = ( ( ATMNUM(BONDS(I)%I) == 1 ) .OR. ( ATMNUM(BONDS(I)%J) == 1 ) )
            ! . Types.
            DO T = 1,NUMBND
                IF ( ( ABS ( BONDS(I)%EQ - BEQ(T) ) < DIFFTOL ) .AND. ( ABS ( BONDS(I)%FC - BFC(T) ) < DIFFTOL ) ) THEN
                    BINDEX = T
                    GO TO 10
                END IF
            END DO
            NUMBND = NUMBND + 1
            BEQ(NUMBND) = BONDS(I)%EQ
            BFC(NUMBND) = BONDS(I)%FC
            BINDEX = NUMBND
            10 CONTINUE
            ! . Index arrays.
            IF ( QBH ) THEN
                BHINDICES(3*NBONH+1) = 3 * ( BONDS(I)%I - 1 )
                BHINDICES(3*NBONH+2) = 3 * ( BONDS(I)%J - 1 )
                BHINDICES(3*NBONH+3) = BINDEX
                NBONH = NBONH + 1
            ELSE
                BQINDICES(3*MBONA+1) = 3 * ( BONDS(I)%I - 1 )
                BQINDICES(3*MBONA+2) = 3 * ( BONDS(I)%J - 1 )
                BQINDICES(3*MBONA+3) = BINDEX
                MBONA = MBONA + 1
            END IF
        END DO

        ! . Convert BFC to kcal.
        BFC(1:NUMBND) = BFC(1:NUMBND) / KCAL_TO_KJ

        ! . Set some counters.
        NBONA = MBONA

        END SUBROUTINE FINDBONDTYPES

        !------------------------
        SUBROUTINE FINDEXCLUSIONS
        !------------------------

        ! . Local variables.
        INTEGER :: I, IEX, N

        ! . Allocate space.
        ALLOCATE ( EXCLI(1:NATOMS), EXCLJ(1:SIZE(ATMEXCJ)+NATOMS) ) ! . Treat case where an atom (in the limit all atoms) have no exclusions.

        ! . Find exclusion arrays.
        NNB = 0
        DO I = 1,NATOMS

            ! . Number of exclusions for the atom.
            N = ATMEXCI(I+1) - ATMEXCI(I)

            ! . Zero exclusions - add 1 zero exclusion.
            IF ( N <= 0 ) THEN
                NNB = NNB + 1
                EXCLI(I)   = 1
                EXCLJ(NNB) = 0
            ! . Non-zero exclusions.
            ELSE
                EXCLI(I) = N
                DO IEX = ATMEXCI(I)+1,ATMEXCI(I+1)
                    NNB = NNB + 1
                    EXCLJ(NNB) = ATMEXCJ(IEX)
                END DO
            END IF
        END DO

        END SUBROUTINE FINDEXCLUSIONS

        !---------------------
        SUBROUTINE FINDLJTYPES
        !---------------------

        ! . Local variables.
        INTEGER  :: I, IJ, J, T
        REAL(DP) :: EIJ, SIJ6
        REAL(DP), ALLOCATABLE, DIMENSION(:) :: LJEPS, LJSIG

        ! . Create the LJ type array - slow version.
        ALLOCATE ( LJEPS(1:NATOMS), LJSIG(1:NATOMS), LJINDICES(1:NATOMS) )
        NLJ = 0
        DO I = 1,NATOMS
            DO T = 1,NLJ
                IF ( ( ABS ( ATMEPS(I) - LJEPS(T) ) < DIFFTOL ) .AND. ( ABS ( ATMSIG(I) - LJSIG(T) ) < DIFFTOL ) ) THEN
                    LJINDICES(I) = T
                    GO TO 10
                END IF
            END DO
            NLJ = NLJ + 1
            LJEPS(NLJ)   = ATMEPS(I)
            LJSIG(NLJ)   = ATMSIG(I)
            LJINDICES(I) = NLJ
            10 CONTINUE
        END DO

        ! . Allocate space.
        ALLOCATE ( LJACOEF(1:(NLJ*(NLJ+1))/2), LJBCOEF(1:(NLJ*(NLJ+1))/2), LJTABLE(1:NLJ,1:NLJ) )

        ! . Create the LJACOEF and LJBCOEF arrays.
        IJ = 0
        DO I = 1,NLJ
            DO J = 1,I
                IJ = IJ + 1
                SIJ6 = ( LJSIG(I) * LJSIG(J) )**6
                EIJ  =   LJEPS(I) * LJEPS(J) * SIJ6
                LJACOEF(IJ) = EIJ * SIJ6
                LJBCOEF(IJ) = EIJ
            END DO
        END DO
        LJACOEF = LJACOEF / KCAL_TO_KJ
        LJBCOEF = LJBCOEF / KCAL_TO_KJ

        ! . Create LJTABLE.
        IJ = 0
        DO I = 1,NLJ
            DO J = 1,(I-1)
                IJ = IJ + 1
                LJTABLE(I,J) = IJ
                LJTABLE(J,I) = IJ
            END DO
            IJ = IJ + 1
            LJTABLE(I,I) = IJ
        END DO

        ! . Finish up.
        DEALLOCATE ( LJEPS, LJSIG )

        END SUBROUTINE FINDLJTYPES

        !---------------------------------
        SUBROUTINE FINDSYMMETRYINFORMATION
        !---------------------------------

        ! . AMBER assumes that atoms in molecules are contiguous with solute appearing before solvent.
        ! . It also appears that IPTRES is not used and NSPSOL is used little except in ptraj. However,
        ! . NSP is necessary.

        ! . Scalars.
        INTEGER :: I, N
        LOGICAL :: QOK

        ! . Arrays.
        INTEGER, ALLOCATABLE, DIMENSION(:) :: INDICES

        ! . Check for symmetry.
        IF ( QBOX .AND. QINCLUDESYMMETRY ) THEN

            ! . Find isolates and use them for molecules.
            ALLOCATE ( INDICES(1:NATOMS) )
            CALL MAKEISOLATEARRAY ( INDICES )

            ! . Check that indices only increase.
            QOK = .TRUE.
            DO I = 2,NATOMS
                IF ( INDICES(I) < INDICES(I-1) ) THEN
                    QOK = .FALSE.
                    EXIT
                END IF
            END DO

            ! . Everything OK.
            IF ( QOK ) THEN

                ! . Find the number of molecules.
                NSPM = MAXVAL ( INDICES )

                ! . Allocate and fill the molecule array.
                ALLOCATE ( NSP(1:NSPM) ) ; NSP = 0
                DO I = 1,NATOMS
                    N = INDICES(I)
                    NSP(N) = NSP(N) + 1
                END DO

                ! . Set the IFBOX flag.
                IFBOX = 1

                ! . Set NSPSOL.
                IF ( IPTRES == 0 ) THEN
                    NSPSOL = 1
                ELSE IF ( IPTRES == NRESID ) THEN
                    NSPSOL = 0
                ELSE
                    N = RESIND(IPTRES+1)
                    NSPSOL = INDICES(N+1)
                END IF

            ! . There is a problem.
            ELSE
                WRITE ( OUTPUT, "(A)" ) "WARNING> Symmetry switched off as there are molecules with non-contiguous atoms."
            END IF

            ! . Finish up.
            DEALLOCATE ( INDICES )

        END IF

        END SUBROUTINE FINDSYMMETRYINFORMATION

        !--------------------------
        SUBROUTINE FINDTORSIONTYPES
        !--------------------------

        ! . This subroutine may not be fully general as far as treatment of 1-4s is concerned.
        ! . A problem is the presence of NULL dihedrals - AMBER includes them, DYNAMO does not.
        ! . Impropers are assumed not to contribute to 1-4s.

        ! . Local variables.
        LOGICAL  :: QH, QP
        INTEGER  :: I, I1, I4, N, NDTEMP, NNULL, T, TINDEX, X, Y
        INTEGER, ALLOCATABLE, DIMENSION(:)            :: IPERIOD, TMPE14I, TMPE14J
        INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: TMPBND
        INTEGER,              DIMENSION(:,:), POINTER :: TMPANG, TMPDIH
 
        ! . Initialize the pointers.
        NULLIFY ( TMPANG, TMPDIH )

        ! . Find the number of unique dihedrals (i.e. excluding multiple terms).
        NDTEMP = 1
        DO I = 2,NDIHEDRALS
            IF ( ( DIHEDRALS(I)%I /= DIHEDRALS(I-1)%I ) .OR. ( DIHEDRALS(I)%J /= DIHEDRALS(I-1)%J ) .OR. &
                 ( DIHEDRALS(I)%K /= DIHEDRALS(I-1)%K ) .OR. ( DIHEDRALS(I)%L /= DIHEDRALS(I-1)%L ) ) THEN
                NDTEMP = NDTEMP + 1
            END IF
        END DO

        ! . Determine the number of NULL dihedral terms (as the difference between 1,4 exclusions and the number of unique dihedrals).
        NNULL = MAX ( 0, SIZE ( ATME14J ) - NDTEMP )

        ! . Create the type and index arrays - slow version.
        N = NDIHEDRALS + NIMPROPERS + NNULL
        ALLOCATE ( IPERIOD(1:N), TFC(1:N), TPERIOD(1:N), TPHASE(1:N), THINDICES(1:5*N), TQINDICES(1:5*N) )

        ! . Copy the exclusion array.
        ALLOCATE ( TMPE14I(1:SIZE(ATME14I)), TMPE14J(1:SIZE(ATME14J)) )
        TMPE14I = ATME14I
        TMPE14J = ATME14J

        ! . Recreate the dihedral list - assuming that all bonds are still present.
        ALLOCATE ( TMPBND(1:2,1:NBONDS) )
        DO I = 1,NBONDS
           TMPBND(1,I) = BONDS(I)%I
           TMPBND(2,I) = BONDS(I)%J
        END DO
        CALL ORDER_BONDS            (         TMPBND )
        CALL CONNECTIVITY_ANGLES    ( TMPANG, TMPBND )
        CALL CONNECTIVITY_DIHEDRALS ( TMPDIH, TMPBND, TMPANG )
        NDTEMP = SIZE ( TMPDIH, 2 )

        ! . Initialization.
        MPHIA = 0
        NPHIH = 0
        NPTRA = 0

        ! . Proper dihedrals.
        DO I = 1,NDIHEDRALS
            ! . Hydrogen-containing.
            QH = ( ( ATMNUM(DIHEDRALS(I)%I) == 1 ) .OR. ( ATMNUM(DIHEDRALS(I)%J) == 1 ) .OR. &
                   ( ATMNUM(DIHEDRALS(I)%K) == 1 ) .OR. ( ATMNUM(DIHEDRALS(I)%L) == 1 ) )

            ! . Check for a previous 1-4 interaction and flag the current one.
            ! . Do not include the 1-4 interaction if it isn't found in the list.
            QP = .FALSE.
            I1 = MIN ( DIHEDRALS(I)%I, DIHEDRALS(I)%L )
            I4 = MAX ( DIHEDRALS(I)%I, DIHEDRALS(I)%L )
            DO X = (TMPE14I(I1)+1),TMPE14I(I1+1)
                Y = TMPE14J(X)
                ! . Already occurred.
                IF ( Y == -I4 ) THEN
                    QP = .TRUE.
                    GO TO 5
                ! . Should be flagged.
                ELSE IF ( Y == I4 ) THEN
                    TMPE14J(X) = - Y
                    GO TO 5
                END IF
            END DO

            ! . Skip interactions not found in the list.
            QP = .TRUE.
            5 CONTINUE

            ! . Types.
            DO T = 1,NPTRA
                IF ( ( ABS ( DIHEDRALS(I)%FC - TFC(T) )       < DIFFTOL ) .AND. &
                     ( ABS ( DIHEDRALS(I)%PHASE - TPHASE(T) ) < DIFFTOL ) .AND. &
                           ( DIHEDRALS(I)%PERIOD          == IPERIOD(T) ) ) THEN
                    TINDEX = T
                    GO TO 10
                END IF
            END DO
            NPTRA = NPTRA + 1
            TFC(NPTRA)     = DIHEDRALS(I)%FC
            TPHASE(NPTRA)  = DIHEDRALS(I)%PHASE
            IPERIOD(NPTRA) = DIHEDRALS(I)%PERIOD
            TINDEX = NPTRA
            10 CONTINUE
            ! . Index arrays.
            IF ( QH ) THEN
                THINDICES(5*NPHIH+1) = 3 * ( DIHEDRALS(I)%I - 1 )
                THINDICES(5*NPHIH+2) = 3 * ( DIHEDRALS(I)%J - 1 )
                IF ( QP ) THEN
                    THINDICES(5*NPHIH+3) = - 3 * ( DIHEDRALS(I)%K - 1 )
                ELSE
                    THINDICES(5*NPHIH+3) =   3 * ( DIHEDRALS(I)%K - 1 )
                END IF
                THINDICES(5*NPHIH+4) = 3 * ( DIHEDRALS(I)%L - 1 )
                THINDICES(5*NPHIH+5) = TINDEX
                NPHIH = NPHIH + 1
            ELSE
                TQINDICES(5*MPHIA+1) = 3 * ( DIHEDRALS(I)%I - 1 )
                TQINDICES(5*MPHIA+2) = 3 * ( DIHEDRALS(I)%J - 1 )
                IF ( QP ) THEN
                    TQINDICES(5*MPHIA+3) = - 3 * ( DIHEDRALS(I)%K - 1 )
                ELSE
                    TQINDICES(5*MPHIA+3) =   3 * ( DIHEDRALS(I)%K - 1 )
                END IF
                TQINDICES(5*MPHIA+4) = 3 * ( DIHEDRALS(I)%L - 1 )
                TQINDICES(5*MPHIA+5) = TINDEX
                MPHIA = MPHIA + 1
            END IF
        END DO

        ! . NULL dihedrals.
        IF ( NNULL > 0 ) THEN

            ! . Add a NULL type.
            NPTRA = NPTRA + 1
            TFC(NPTRA)     = 0.0_DP
            TPHASE(NPTRA)  = 0.0_DP
            IPERIOD(NPTRA) = 2
            TINDEX = NPTRA

            ! . Loop over full dihedral list and write out NULL dihedrals.
            DO I = 1,NDTEMP

                ! . Hydrogen-containing.
                QH = ( ( ATMNUM(TMPDIH(1,I)) == 1 ) .OR. ( ATMNUM(TMPDIH(2,I)) == 1 ) .OR. &
                       ( ATMNUM(TMPDIH(3,I)) == 1 ) .OR. ( ATMNUM(TMPDIH(4,I)) == 1 ) )

                ! . Check for a previous 1-4 interaction or flag a match.
                QP = .FALSE.
                I1 = MIN ( TMPDIH(1,I), TMPDIH(4,I) )
                I4 = MAX ( TMPDIH(1,I), TMPDIH(4,I) )
                DO X = (TMPE14I(I1)+1),TMPE14I(I1+1)
                    Y = TMPE14J(X)
                    ! . Already occurred.
                    IF ( Y == -I4 ) THEN
                        QP = .TRUE.
                        GO TO 20
                    ! . Should be flagged and included.
                    ELSE IF ( Y == I4 ) THEN
                        TMPE14J(X) = - Y
                        GO TO 20
                    END IF
                END DO

                ! . Skip any 1-4 interactions that are not found.
                QP = .TRUE.
                20 CONTINUE

                ! . Only write out those that haven't already occurred.
                IF ( .NOT. QP ) THEN
                    IF ( QH ) THEN
                        THINDICES(5*NPHIH+1) = 3 * ( TMPDIH(1,I) - 1 )
                        THINDICES(5*NPHIH+2) = 3 * ( TMPDIH(2,I) - 1 )
                        THINDICES(5*NPHIH+3) = 3 * ( TMPDIH(3,I) - 1 )
                        THINDICES(5*NPHIH+4) = 3 * ( TMPDIH(4,I) - 1 )
                        THINDICES(5*NPHIH+5) = TINDEX
                        NPHIH = NPHIH + 1
                    ELSE
                        TQINDICES(5*MPHIA+1) = 3 * ( TMPDIH(1,I) - 1 )
                        TQINDICES(5*MPHIA+2) = 3 * ( TMPDIH(2,I) - 1 )
                        TQINDICES(5*MPHIA+3) = 3 * ( TMPDIH(3,I) - 1 )
                        TQINDICES(5*MPHIA+4) = 3 * ( TMPDIH(4,I) - 1 )
                        TQINDICES(5*MPHIA+5) = TINDEX
                        MPHIA = MPHIA + 1
                    END IF
                END IF
            END DO

        END IF

        ! . Improper dihedrals - 3rd and 4th indices always negative.
        DO I = 1,NIMPROPERS
            ! . Hydrogen-containing.
            QH = ( ( ATMNUM(IMPROPERS(I)%I) == 1 ) .OR. ( ATMNUM(IMPROPERS(I)%J) == 1 ) .OR. &
                   ( ATMNUM(IMPROPERS(I)%K) == 1 ) .OR. ( ATMNUM(IMPROPERS(I)%L) == 1 ) )
            ! . Types.
            DO T = 1,NPTRA
                IF ( ( ABS ( IMPROPERS(I)%FC - TFC(T) )       < DIFFTOL ) .AND. &
                     ( ABS ( IMPROPERS(I)%PHASE - TPHASE(T) ) < DIFFTOL ) .AND. &
                           ( IMPROPERS(I)%PERIOD          == IPERIOD(T) ) ) THEN
                    TINDEX = T
                    GO TO 30
                END IF
            END DO
            NPTRA = NPTRA + 1
            TFC(NPTRA)     = IMPROPERS(I)%FC
            TPHASE(NPTRA)  = IMPROPERS(I)%PHASE
            IPERIOD(NPTRA) = IMPROPERS(I)%PERIOD
            TINDEX = NPTRA
            30 CONTINUE
            ! . Index arrays.
            IF ( QH ) THEN
                THINDICES(5*NPHIH+1) =   3 * ( IMPROPERS(I)%I - 1 )
                THINDICES(5*NPHIH+2) =   3 * ( IMPROPERS(I)%J - 1 )
                THINDICES(5*NPHIH+3) = - 3 * ( IMPROPERS(I)%K - 1 )
                THINDICES(5*NPHIH+4) = - 3 * ( IMPROPERS(I)%L - 1 )
                THINDICES(5*NPHIH+5) = TINDEX
                NPHIH = NPHIH + 1
            ELSE
                TQINDICES(5*MPHIA+1) =   3 * ( IMPROPERS(I)%I - 1 )
                TQINDICES(5*MPHIA+2) =   3 * ( IMPROPERS(I)%J - 1 )
                TQINDICES(5*MPHIA+3) = - 3 * ( IMPROPERS(I)%K - 1 )
                TQINDICES(5*MPHIA+4) = - 3 * ( IMPROPERS(I)%L - 1 )
                TQINDICES(5*MPHIA+5) = TINDEX
                MPHIA = MPHIA + 1
            END IF
        END DO

        ! . Convert AFC to kcal.
        TFC(1:NPTRA) = TFC(1:NPTRA) / KCAL_TO_KJ

        ! . Fill TPERIOD.
        TPERIOD = REAL ( IPERIOD, DP )

        ! . Convert TPHASE.
        TPHASE = ACOS ( TPHASE )

        ! . Set some counters.
        NPHIA = MPHIA

        ! . Deallocate space.
        DEALLOCATE ( IPERIOD, TMPANG, TMPBND, TMPDIH, TMPE14I, TMPE14J )

        END SUBROUTINE FINDTORSIONTYPES

        !--------------------------------------
        SUBROUTINE MAKEISOLATEARRAY ( INDICES )
        !--------------------------------------

        ! . Make the isolate indices.

        ! . Arguments.
        INTEGER, DIMENSION(1:NATOMS), INTENT(OUT) :: INDICES

        ! . Scalars.
        INTEGER :: C, I, J, M, N, NISOLATE, S

        ! . Arrays.
        INTEGER, ALLOCATABLE, DIMENSION(:) :: CONNECTIONSI, CONNECTIONSJ, STACK

        ! . Allocate space and initialize.
        ALLOCATE ( CONNECTIONSI(1:NATOMS+1), CONNECTIONSJ(1:2*NBONDS), STACK(1:NATOMS) )
        CONNECTIONSI = 0 ; CONNECTIONSJ = 0 ; STACK = 0

        ! . Find the number of connections per atom.
        DO N = 1,NBONDS
            I = BONDS(N)%I
            J = BONDS(N)%J
            STACK(I) = STACK(I) + 1
            STACK(J) = STACK(J) + 1
        END DO

        ! . Generate the index array.
        N = 0
        DO I = 1,NATOMS
            N = N + STACK(I)
            CONNECTIONSI(I+1) = N
        END DO

        ! . Fill CONNECTIONSJ.
        STACK = 0
        DO N = 1,NBONDS
            I = BONDS(N)%I
            J = BONDS(N)%J
            STACK(I) = STACK(I) + 1
            STACK(J) = STACK(J) + 1
            CONNECTIONSJ(CONNECTIONSI(I)+STACK(I)) = J
            CONNECTIONSJ(CONNECTIONSI(J)+STACK(J)) = I
        END DO

        ! . Initialization.
        INDICES  = 0
        NISOLATE = 0
        STACK    = 0

        ! . Loop over all atoms.
        DO S = 1,NATOMS
            IF ( INDICES(S) <= 0 ) THEN

                ! . Start the new isolate.
                NISOLATE   = NISOLATE + 1
                INDICES(S) = NISOLATE

                ! . Put the point on the stack.
                N = 1
                STACK(N) = S

                ! . Loop over all atoms in the stack.
                M = 0
                DO WHILE ( M < N )
                    M = M + 1
                    I = STACK(M)

                    ! . Loop over connections for the atom.
                    DO C = (CONNECTIONSI(I)+1),CONNECTIONSI(I+1)
                        J = CONNECTIONSJ(C)
                        IF ( INDICES(J) <= 0 ) THEN
                            INDICES(J) = NISOLATE
                            N = N + 1
                            STACK(N) = J
                        END IF
                    END DO
                END DO
            END IF
        END DO

        ! . Finish up.
        DEALLOCATE ( CONNECTIONSI, CONNECTIONSJ, STACK )
        IF ( ANY ( INDICES <= 0 ) ) CALL PRINT_ERROR ( "AMBERCRD_WRITE", "Invalid molecule index generation." )

        END SUBROUTINE MAKEISOLATEARRAY

    END SUBROUTINE AMBERTOP_WRITE

END MODULE AMBERFILE_IO
