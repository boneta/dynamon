MODULE ZBND

#ifdef	HACK_ZBND

! .Module declarations.
USE CONSTANTS,          ONLY : TO_DEGREES
USE DEFINITIONS,        ONLY : DP
USE FILES,              ONLY : NEXT_UNIT
USE ATOMS,              ONLY : ATMCRD, NATOMS, symbol, atmnum
USE PRINTING,           ONLY : PRINT_ERROR, PRINT_LINE, PRINT_TEXT

IMPLICIT NONE
PRIVATE
PUBLIC :: ZBND_DEFINE, CALC_ZBND, ZBND_WRITING_START, &
          ZBND_WRITING_STOP, ZBND_CLEAN, zbnd_mcb, zbnd_micb
SAVE

REAL( KIND=DP ) :: Z_EQ, Z_FC, X_FC, Y_FC
REAL( KIND=DP ), DIMENSION(1:3) :: ZBND_CEN
REAL( KIND=DP ), DIMENSION(1:3,1:3) :: ZBND_MCB, ZBND_MICB
INTEGER, DIMENSION(:), ALLOCATABLE :: ZBND_IDX
REAL( KIND=DP ), DIMENSION(:,:), ALLOCATABLE :: ZBND_CRD, ZBND_XXX
LOGICAL :: ZBND_WRITE = .FALSE., ZBND_CALC = .FALSE.
INTEGER :: ZBND_UNIT = -1, ZBND_NAT, ZBND_WHO, ZBND_OHW

!===============================================================================
CONTAINS
!===============================================================================

   !--------------------------------------
   SUBROUTINE CALC_ZBND( ECORR, GRADIENT )
   !--------------------------------------
   REAL ( KIND = DP ), INTENT(INOUT) :: ECORR
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT), OPTIONAL :: GRADIENT

   INTEGER :: I
   REAL( KIND = DP ) :: SZ, ZZ, DZ, DX, DY
   REAL( KIND = DP ), DIMENSION(1:3) :: GZ, GXY

   ECORR = 0.0_DP
   IF( .NOT. ZBND_CALC ) RETURN

!write(900,"(i6/)") zbnd_nat
   DO I = 1, ZBND_NAT
       ZBND_XXX(1:3,I) = MATMUL( ZBND_MCB(1:3,1:3), ATMCRD(1:3,ZBND_IDX(I)) - ZBND_CEN(1:3) )
!write(900,"(a2,3f20.10)") symbol(atmnum(zbnd_idx(i))), zbnd_xxx(1:3,i)
   END DO

   ZZ = DABS( ZBND_XXX(3,ZBND_WHO) )
   DZ = ZZ - Z_EQ
   ECORR = 0.5_DP * Z_FC * DZ * DZ
   IF( PRESENT( GRADIENT ) ) THEN
       SZ = DSIGN( 1._DP, ZBND_XXX(3,ZBND_WHO) )
       GZ(1:3) = MATMUL( ZBND_MICB(1:3,1:3), (/ .0_DP, .0_DP, Z_FC * DZ * SZ /) )
!       IF( X_FC == .0_DP .AND. Y_FC == .0_DP ) THEN
           DO I = 1, ZBND_NAT
               GRADIENT(1:3,ZBND_IDX(I)) = GRADIENT(1:3,ZBND_IDX(I)) + GZ(1:3)
           END DO
!       ELSE
!           GRADIENT(1:3,ZBND_OHW) = GRADIENT(1:3,ZBND_OHW) + GZ(1:3)
!       END IF
   END IF

   DO I = 1, ZBND_NAT
       DX = ZBND_XXX(1,I) - ZBND_CRD(1,I)
       ECORR = ECORR + 0.5_DP * X_FC * DX * DX
       DY = ZBND_XXX(2,I) - ZBND_CRD(2,I)
       ECORR = ECORR + 0.5_DP * Y_FC * DY * DY
       IF( PRESENT( GRADIENT ) ) THEN
           GXY(1:3) = MATMUL( ZBND_MICB(1:3,1:3), (/ X_FC * DX, Y_FC * DY, .0_DP /) )
           GRADIENT(1:3,ZBND_IDX(I)) = GRADIENT(1:3,ZBND_IDX(I)) + GXY(1:3)
       END IF
   END DO

   IF( ZBND_WRITE ) THEN
       WRITE( ZBND_UNIT, "(F20.10)" ) ZZ
       CALL FLUSH( ZBND_UNIT )
   END IF
   END SUBROUTINE CALC_ZBND


   !----------------------------------
   SUBROUTINE ZBND_WRITING_START
   !----------------------------------
   IF( ZBND_UNIT > -1 ) THEN
       ZBND_WRITE = .TRUE.
       WRITE( ZBND_UNIT, "(2F20.10)" ) Z_FC, Z_EQ
   END IF
   END SUBROUTINE ZBND_WRITING_START


   !----------------------------------
   SUBROUTINE ZBND_WRITING_STOP
   !----------------------------------
   IF( ZBND_WRITE ) THEN
       CALL FLUSH( ZBND_UNIT )
       CLOSE( ZBND_UNIT )
       ZBND_UNIT  = -1
       ZBND_WRITE = .FALSE.
   END IF
   END SUBROUTINE ZBND_WRITING_STOP


   !-------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ZBND_DEFINE ( ORIGIN, DIRECTION, SELECTION, Z_ATOM, REFERENCE, Z_CONSTANT, X_CONSTANT, Y_CONSTANT, FILE )
   !-------------------------------------------------------------------------------------------------------------------
   REAL( KIND=DP ), DIMENSION(1:3), INTENT(IN) :: ORIGIN, DIRECTION
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN) :: SELECTION
   INTEGER, INTENT(IN) :: Z_ATOM
   REAL( KIND=DP ), INTENT(IN) :: REFERENCE, Z_CONSTANT, X_CONSTANT, Y_CONSTANT
   CHARACTER( LEN=* ), INTENT(IN), OPTIONAL :: FILE

   INTEGER :: I, J
   REAL( KIND=DP ) :: TMP
   REAL( KIND=DP ), DIMENSION(1:3) :: V1, V2, V3

   ZBND_CEN(1:3) = ORIGIN(1:3)

   Z_EQ = REFERENCE
   Z_FC = Z_CONSTANT
   X_FC = X_CONSTANT
   Y_FC = Y_CONSTANT

   ZBND_OHW = Z_ATOM
   ZBND_WHO = -1
   ZBND_NAT = COUNT( SELECTION )
   ALLOCATE( ZBND_IDX(1:ZBND_NAT), ZBND_CRD(1:3,1:ZBND_NAT), ZBND_XXX(1:3,1:ZBND_NAT) )
   J = 1
   DO I = 1, NATOMS
       IF( SELECTION(I) ) THEN
           ZBND_IDX(J) = I
           IF( I == Z_ATOM ) ZBND_WHO = J
           J = J + 1
       END IF
   END DO

   IF( PRESENT( FILE ) ) THEN
       ZBND_UNIT = NEXT_UNIT()
       OPEN( UNIT = ZBND_UNIT, FILE = TRIM( FILE ), ACTION = "WRITE", FORM = "FORMATTED" )
   END IF

   ZBND_CALC  = .TRUE.

   V1(1:3) = DIRECTION(1:3) / DSQRT( DOT_PRODUCT( DIRECTION(1:3), DIRECTION(1:3) ) )
   IF( V1(3) /= .0_DP ) THEN
       V2(1:3) = (/ 1.0_DP, 1.0_DP, - ( V1(1) + V1(2) ) / V1(3) /)
   ELSE
       IF( V1(2) /= .0_DP ) THEN
           V2(1:3) = (/ 1._DP, - V1(1) / V1(2), .0_DP /)
       ELSE
           V2(1:3) = (/ .0_DP, 1._DP, .0_DP /)
       END IF
   END IF
   V2(1:3) = V2(1:3) / DSQRT( DOT_PRODUCT( V2(1:3), V2(1:3) ) )
   V3(1) = V1(2) * V2(3) - V1(3) * V2(2)
   V3(2) = V1(3) * V2(1) - V1(1) * V2(3)
   V3(3) = V1(1) * V2(2) - V1(2) * V2(1)
   V3(1:3) = V3(1:3) / DSQRT( DOT_PRODUCT( V3(1:3), V3(1:3) ) )

   ZBND_MCB(1:3,1) = (/ V3(1), V2(1), v1(1) /)
   ZBND_MCB(1:3,2) = (/ V3(2), V2(2), v1(2) /)
   ZBND_MCB(1:3,3) = (/ V3(3), V2(3), v1(3) /)
!   ZBND_MCB(1:3,1) = V3(1:3)
!   ZBND_MCB(1:3,2) = V2(1:3)
!   ZBND_MCB(1:3,3) = V1(1:3)

! - should be orthogonal
!   ZBND_MICB(1,1) = ZBND_MCB(2,2) * ZBND_MCB(3,3) - ZBND_MCB(2,3) * ZBND_MCB(3,2)
!   ZBND_MICB(2,1) = ZBND_MCB(2,3) * ZBND_MCB(3,1) - ZBND_MCB(2,1) * ZBND_MCB(3,3)
!   ZBND_MICB(3,1) = ZBND_MCB(2,1) * ZBND_MCB(3,2) - ZBND_MCB(2,2) * ZBND_MCB(3,1)
!   ZBND_MICB(1,2) = ZBND_MCB(1,3) * ZBND_MCB(3,2) - ZBND_MCB(1,2) * ZBND_MCB(3,3)
!   ZBND_MICB(2,2) = ZBND_MCB(1,1) * ZBND_MCB(3,3) - ZBND_MCB(1,3) * ZBND_MCB(3,1)
!   ZBND_MICB(3,2) = ZBND_MCB(1,2) * ZBND_MCB(3,1) - ZBND_MCB(1,1) * ZBND_MCB(3,2)
!   ZBND_MICB(1,3) = ZBND_MCB(1,2) * ZBND_MCB(2,3) - ZBND_MCB(1,3) * ZBND_MCB(2,2)
!   ZBND_MICB(2,3) = ZBND_MCB(1,3) * ZBND_MCB(2,1) - ZBND_MCB(1,1) * ZBND_MCB(2,3)
!   ZBND_MICB(3,3) = ZBND_MCB(1,1) * ZBND_MCB(2,2) - ZBND_MCB(1,2) * ZBND_MCB(2,1)
!   TMP = ZBND_MCB(1,1) * ( ZBND_MCB(2,2) * ZBND_MCB(3,3) - ZBND_MCB(2,3) * ZBND_MCB(3,2) ) - &
!         ZBND_MCB(1,2) * ( ZBND_MCB(2,1) * ZBND_MCB(3,3) - ZBND_MCB(3,1) * ZBND_MCB(2,3) ) + &
!         ZBND_MCB(1,3) * ( ZBND_MCB(2,1) * ZBND_MCB(3,2) - ZBND_MCB(3,1) * ZBND_MCB(2,2) )
!   ZBND_MICB = ZBND_MICB / TMP
   ZBND_MICB = TRANSPOSE( ZBND_MCB )

   DO I = 1, ZBND_NAT
       ZBND_CRD(1:3,I) = MATMUL( ZBND_MCB(1:3,1:3), ATMCRD(1:3,ZBND_IDX(I)) - ZBND_CEN(1:3) )
   END DO

   WRITE( PRINT_LINE, "(A)" ) "- Z-Binding Constraint Setup -----------------------------------------------"
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A16,3F12.6)" ) "Origin:", ZBND_CEN
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A16,3F12.6)" ) "Direction:", V1
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A16,2I20)" ) "Z-Atom:", ZBND_OHW, ZBND_WHO
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A16,F20.10)" ) "Reference:", Z_EQ
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A16,F20.10)" ) "Force Constant:", Z_FC
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A16,F20.10)" ) "X-Force Const.:", X_FC
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A16,F20.10)" ) "Y-Force Const.:", Y_FC
   CALL PRINT_TEXT
   WRITE( PRINT_LINE, "(A16,I20)" ) "Selected Atoms:", ZBND_NAT
   CALL PRINT_TEXT
   IF( ZBND_UNIT > -1 ) THEN
       WRITE( PRINT_LINE, "(A16,A)" ) "File:", TRIM( FILE )
       CALL PRINT_TEXT
   END IF
write(*,"(a)") "----------------------------------------------------------------------------"
write(*,"(a)")                 "graphics top color 4"
write(*,"(a,3f8.3,a,3f8.3,a)") "graphics top cylinder {",zbnd_cen,"} {",zbnd_cen+v1,"} radius .1 resolution 10 filled yes"
write(*,"(a,3f8.3,a,3f8.3,a)") "graphics top cone     {",zbnd_cen+v1,"} {",zbnd_cen+v1+v1,"} radius .2 resolution 10"
   WRITE( PRINT_LINE, "(A)" ) "----------------------------------------------------------------------------"
   CALL PRINT_TEXT
   END SUBROUTINE ZBND_DEFINE


   SUBROUTINE ZBND_CLEAN
   IF( ALLOCATED( ZBND_IDX ) )  DEALLOCATE( ZBND_IDX )
   IF( ALLOCATED( ZBND_CRD ) )  DEALLOCATE( ZBND_CRD )
   IF( ALLOCATED( ZBND_XXX ) )  DEALLOCATE( ZBND_XXX )
   END SUBROUTINE ZBND_CLEAN

#endif

END MODULE ZBND
