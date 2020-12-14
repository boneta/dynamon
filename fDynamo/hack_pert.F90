MODULE PERT

#ifdef   HACK_PERT

! .Module declarations.
USE DEFINITIONS,        ONLY : DP
USE FILES,              ONLY : NEXT_UNIT
USE ATOMS,              ONLY : NATOMS
USE PRINTING,           ONLY : PRINT_ERROR, PRINT_LINE, PRINT_TEXT

IMPLICIT NONE
PUBLIC
SAVE

REAL( KIND=DP ) :: PRV_PERT_EL, PRV_PERT_LJ, PRV_PERT_QM
REAL( KIND=DP ) :: CUR_PERT_EL, CUR_PERT_LJ, CUR_PERT_QM
REAL( KIND=DP ) :: NXT_PERT_EL, NXT_PERT_LJ, NXT_PERT_QM
REAL( KIND=DP ) :: LAMBDA_VALUE = 1._DP, LAMBDA_DELTA = 0._DP
LOGICAL :: PERT_WRITE, PERT_CALC_EL = .FALSE., PERT_CALC_LJ = .FALSE.
INTEGER :: PERT_UNIT
LOGICAL, DIMENSION(:), ALLOCATABLE :: PERT_SELE

!===============================================================================
CONTAINS
!===============================================================================

   !----------------------------
   SUBROUTINE PERT_WRITING_START
   !----------------------------

   PERT_WRITE = .TRUE.

   IF( PERT_UNIT > -1 ) THEN
       WRITE( PERT_UNIT, "(A,2F20.10)" ) "#", LAMBDA_VALUE, LAMBDA_DELTA
   ELSE
       PERT_WRITE = .FALSE.
   END IF

   END SUBROUTINE PERT_WRITING_START


   !-------------------------
   SUBROUTINE PERT_WRITING_STOP
   !-------------------------

   IF( PERT_WRITE .AND. PERT_UNIT > -1 ) THEN
       CALL FLUSH( PERT_UNIT )
       CLOSE( PERT_UNIT )
       PERT_UNIT  = -1
       PERT_WRITE = .FALSE.
   END IF

   END SUBROUTINE PERT_WRITING_STOP


   !------------------------------------------------------------------
   SUBROUTINE PERT_DEFINE( METHOD, SELECTION, L_VALUE, L_DELTA, FILE )
   !------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN=* ), INTENT(IN) :: METHOD
   LOGICAL, DIMENSION(1:NATOMS) :: SELECTION
   REAL( KIND=DP ), INTENT(IN) :: L_VALUE, L_DELTA
   CHARACTER( LEN=* ), INTENT(IN), OPTIONAL :: FILE

   IF( ALLOCATED( PERT_SELE ) ) DEALLOCATE( PERT_SELE )
   ALLOCATE( PERT_SELE(1:NATOMS) )

   PERT_SELE    = SELECTION
   PERT_UNIT    = -1
   PERT_WRITE   = .FALSE.
   LAMBDA_VALUE = L_VALUE
   LAMBDA_DELTA = L_DELTA

   PERT_CALC_EL = ( TRIM( METHOD ) == "ELECTROSTATIC" )
   PERT_CALC_LJ = ( TRIM( METHOD ) == "LENNARD-JONES" )

   IF( .NOT. ( PERT_CALC_EL .OR. PERT_CALC_LJ ) ) STOP "* Perturbation: ELECTROSTATIC or LENNARD-JONES"

   IF( PRESENT( FILE ) ) THEN
       PERT_UNIT  = NEXT_UNIT()
       OPEN( UNIT = PERT_UNIT, FILE = TRIM( FILE ), ACTION = "WRITE", FORM = "FORMATTED" )
   END IF

   WRITE( PRINT_LINE, "(A)" ) ""
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A)" ) "----- Free Energy Perturbation -------------------------------------------------"
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,A14)" )  "Method:     ", TRIM( METHOD )
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,F14.6)" ) "Lambda value:", LAMBDA_VALUE
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A,F14.6)" ) "Lambda delta:", LAMBDA_DELTA
   CALL PRINT_TEXT
   IF( PERT_UNIT > -1 ) THEN
       WRITE( PRINT_LINE, "(A,A14)" )  "File:       ", TRIM( FILE )
       CALL PRINT_TEXT
   END IF
   WRITE( PRINT_LINE, "(A)" )  "--------------------------------------------------------------------------------"
   CALL PRINT_TEXT

   END SUBROUTINE PERT_DEFINE

#endif

END MODULE PERT
