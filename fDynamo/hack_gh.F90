MODULE GH

#ifdef	HACK_GH

! . Module declarations.
USE DEFINITIONS,      ONLY : DP
USE FILES,            ONLY : NEXT_UNIT
USE ATOMS,            ONLY : ATMCRD, NATOMS, ATMMAS
USE VELOCITY,         ONLY : ATMVEL
USE POTENTIAL_ENERGY, ONLY : ATMDER, ATMDER_CONSTRAINT

! -- Module for the GH calculation of a transfer coordinate

IMPLICIT NONE
PRIVATE
PUBLIC :: GH_FORCES, GH_VV_INIT, GH_VV_A, GH_VV_B, GH_VV_C, GH_DEFINE
SAVE

LOGICAL :: GH_ON = .FALSE.
INTEGER :: A_ATOM, D_ATOM, T_ATOM
INTEGER       :: ITERSHAKE, uno, dos, tres, GH_UNIT = -1
REAL(KIND=DP) :: M1X, M1Y, M1Z, M2X, M2Y, M2Z, RAB, RBC, F1, F2, F11
REAL(KIND=DP) :: D1, D2, D3, D4, D5, D6, D7, D8, D9
REAL(KIND=DP) :: B1, B2, B3, B4, B5, B6, B7, B8, B9, SRC !,FATO
REAL(KIND=DP) :: KILO1, KILO2, KILO3, TMP, vr_x, vr_y, vr_z, vt_x, vt_y, vt_z, ve_x, ve_y, ve_z
REAL(KIND=DP) :: RC_OLD, RC_NEW, DELTARC, FACMIO, XACR2
REAL(KIND=DP), DIMENSION(3,3) :: REFCRD, REFCRD2, TMPCRD, U, U_o, U_t, TMPCRD2, TMPCRD3 !, FATM

CONTAINS

SUBROUTINE GH_FORCES
   IMPLICIT NONE

   IF( .NOT. GH_ON ) RETURN

   M1X = .0_DP; M1Y = .0_DP; M1Z = .0_DP
   M2X = .0_DP; M2Y = .0_DP; M2Z = .0_DP

   RAB = .0_DP; RBC = .0_DP

   F11 = .0_DP; F1 = .0_DP; F2 = .0_DP

   M1X = - ATMCRD (1,T_ATOM) + ATMCRD(1,A_ATOM)
   M1Y = - ATMCRD (2,T_ATOM) + ATMCRD(2,A_ATOM)
   M1Z = - ATMCRD (3,T_ATOM) + ATMCRD(3,A_ATOM)

   M2X = - ATMCRD (1,T_ATOM) + ATMCRD(1,D_ATOM)
   M2Y = - ATMCRD (2,T_ATOM) + ATMCRD(2,D_ATOM)
   M2Z = - ATMCRD (3,T_ATOM) + ATMCRD(3,D_ATOM)
              
   RAB = DSQRT ( M1X**2 + M1Y**2 + M1Z **2)
   RBC = DSQRT ( M2X**2 + M2Y**2 + M2Z **2)

   B1 = (M1X/RAB)
   B2 = (M1Y/RAB)
   B3 = (M1Z/RAB)
   B4 = (M2X/RBC - M1X/RAB)
   B5 = (M2Y/RBC - M1Y/RAB)
   B6 = (M2Z/RBC - M1Z/RAB)
   B7 = - (M2X/RBC)
   B8 = - (M2Y/RBC)
   B9 = - (M2Z/RBC)

  ATMDER = ATMDER - ATMDER_CONSTRAINT

   F2 = ( B1*ATMDER(1,A_ATOM)/ATMMAS(A_ATOM) + B2*ATMDER(2,A_ATOM)/ATMMAS(A_ATOM) + B3*ATMDER(3,A_ATOM)/ATMMAS(A_ATOM) +&
		B4*ATMDER(1,T_ATOM)/ATMMAS(T_ATOM) + B5*ATMDER(2,T_ATOM)/ATMMAS(T_ATOM) + B6*ATMDER(3,T_ATOM)/ATMMAS(T_ATOM) +&
		B7*ATMDER(1,D_ATOM)/ATMMAS(D_ATOM) + B8*ATMDER(2,D_ATOM)/ATMMAS(D_ATOM) + B9*ATMDER(3,D_ATOM)/ATMMAS(D_ATOM) )

  ATMDER = ATMDER + ATMDER_CONSTRAINT
              
   F11= ( (B1**2)/ATMMAS(A_ATOM) + (B2**2)/ATMMAS(A_ATOM) + (B3**2)/ATMMAS(A_ATOM) +& 
		(B4**2)/ATMMAS(T_ATOM) + (B5**2)/ATMMAS(T_ATOM) + (B6**2)/ATMMAS(T_ATOM) +&
		(B7**2)/ATMMAS(D_ATOM) + (B8**2)/ATMMAS(D_ATOM) + (B9**2)/ATMMAS(D_ATOM) )
              
   F1 = 1._DP / F11

   SRC = B1*ATMVEL(1,A_ATOM) + B2*ATMVEL(2,A_ATOM) + B3*ATMVEL(3,A_ATOM) + &
		B4*ATMVEL(1,T_ATOM) + B5*ATMVEL(2,T_ATOM) + B6*ATMVEL(3,T_ATOM) + &
		B7*ATMVEL(1,D_ATOM) + B8*ATMVEL(2,D_ATOM) + B9*ATMVEL(3,D_ATOM) 

   WRITE( GH_UNIT,"(F12.6,F12.6,F12.6,F12.6)") (RBC-RAB), F1*F2, F1, SRC

END SUBROUTINE GH_FORCES

SUBROUTINE GH_VV_INIT( FACV, FACR2 )
   IMPLICIT NONE
   REAL( KIND=DP ), INTENT(IN) :: FACV, FACR2

   IF( .NOT. GH_ON ) RETURN

   FACMIO = FACV / FACR2
   XACR2 = FACR2
   ATMVEL(1:3,A_ATOM) = 0.0_DP
   ATMVEL(1:3,D_ATOM) = 0.0_DP
   ATMVEL(1:3,T_ATOM) = 0.0_DP

END SUBROUTINE GH_VV_INIT

SUBROUTINE GH_VV_A
   IMPLICIT NONE

   IF( .NOT. GH_ON ) RETURN

   M1X = .0_DP; M1Y = .0_DP; M1Z = .0_DP
   M2X = .0_DP; M2Y = .0_DP; M2Z = .0_DP

   RAB = .0_DP; RBC = .0_DP

   F11 = .0_DP; F1 = .0_DP; F2 = .0_DP

   REFCRD (1,1) = ATMCRD (1,A_ATOM)
   REFCRD (2,1) = ATMCRD (2,A_ATOM)
   REFCRD (3,1) = ATMCRD (3,A_ATOM)
   REFCRD (1,2) = ATMCRD (1,T_ATOM)
   REFCRD (2,2) = ATMCRD (2,T_ATOM)
   REFCRD (3,2) = ATMCRD (3,T_ATOM)
   REFCRD (1,3) = ATMCRD (1,D_ATOM)
   REFCRD (2,3) = ATMCRD (2,D_ATOM)
   REFCRD (3,3) = ATMCRD (3,D_ATOM)

END SUBROUTINE GH_VV_A

SUBROUTINE GH_VV_B( ISTEP )
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: ISTEP

   IF( .NOT. GH_ON ) RETURN

   itershake = 0

   vr_x = REFCRD (1,1) - REFCRD(1,3)
   vr_y = REFCRD (2,1) - REFCRD(2,3)
   vr_z = REFCRD (3,1) - REFCRD(3,3)

   Tmp = DSQRT ( vr_x**2 + vr_y**2 + vr_z**2 )

   vr_x = vr_x / Tmp
   vr_y = vr_y / Tmp
   vr_z = vr_z / Tmp

   vt_x = REFCRD(1,2) - REFCRD(1,3)
   vt_y = REFCRD(2,2) - REFCRD(2,3) 
   vt_z = REFCRD(3,2) - REFCRD(3,3)

   Tmp = DSQRT ( vt_x**2 + vt_y**2 + vt_z**2 )

   vt_x = vt_x / Tmp
   vt_y = vt_y / Tmp
   vt_z = vt_z / Tmp

   U_o = 0.0_DP

   U_o (1,1) = vr_y * vt_z  - vr_z * vt_y
   U_o (1,2) = - vr_x * vt_z  + vt_x * vr_z
   U_o (1,3) = vr_x * vt_y  - vr_y * vt_x

   Tmp = DSQRT ( U_o (1,1)**2 + U_o (1,2)**2 + U_o (1,3)**2 )

   U_o (1,1) = U_o (1,1) / Tmp
   U_o (1,2) = U_o (1,2) / Tmp
   U_o (1,3) = U_o (1,3) / Tmp
   Tmp = DACOS ( U_o(1,1)*vr_x + U_o(1,2)*vr_y + U_o(1,3)*vr_z )

   U_o (2,1) = vr_y * U_o (1,3) - vr_z * U_o (1,2)
   U_o (2,2) = - vr_x * U_o (1,3) + vr_z * U_o (1,1)
   U_o (2,3) = vr_x * U_o (1,2) - vr_y * U_o (1,1)

   Tmp = DSQRT ( U_o (2,1)**2 + U_o (2,2)**2 + U_o (2,3)**2 )

   U_o (2,1) = U_o (2,1) / Tmp
   U_o (2,2) = U_o (2,2) / Tmp
   U_o (2,3) = U_o (2,3) / Tmp

   Tmp = DACOS ( U_o(1,1)*U_o(2,1) + U_o(1,2)*U_o(2,2) + U_o(1,3)*U_o(2,3) )

   U_o (3,1) =  vr_x
   U_o (3,2) =  vr_y
   U_o (3,3) =  vr_z

   Tmp = DSQRT ( U_o (3,1)**2 + U_o (3,2)**2 + U_o (3,3)**2 )

   U_o (3,1) = U_o (3,1) / Tmp
   U_o (3,2) = U_o (3,2) / Tmp
   U_o (3,3) = U_o (3,3) / Tmp

   vr_x = ATMCRD (1,A_ATOM) - ATMCRD(1,D_ATOM)
   vr_y = ATMCRD (2,A_ATOM) - ATMCRD(2,D_ATOM)
   vr_z = ATMCRD (3,A_ATOM) - ATMCRD(3,D_ATOM)

   Tmp = DSQRT ( vr_x**2 + vr_y**2 + vr_z**2 )

   vr_x = vr_x / Tmp
   vr_y = vr_y / Tmp
   vr_z = vr_z / Tmp

   vt_x = ATMCRD(1,T_ATOM)-ATMCRD(1,D_ATOM)
   vt_y = ATMCRD(2,T_ATOM)-ATMCRD(2,D_ATOM) 
   vt_z = ATMCRD(3,T_ATOM)-ATMCRD(3,D_ATOM)

   Tmp = DSQRT ( vt_x**2 + vt_y**2 + vt_z**2 )

   vt_x = vt_x / Tmp
   vt_y = vt_y / Tmp
   vt_z = vt_z / Tmp

   U (1,1) = vr_y * vt_z  - vr_z * vt_y
   U (1,2) = - vr_x * vt_z  + vr_z * vt_x
   U (1,3) = vr_x * vt_y  - vr_y * vt_x

   Tmp = DSQRT ( U (1,1)**2 + U (1,2)**2 + U (1,3)**2 )

   U (1,1) = U (1,1) / Tmp
   U (1,2) = U (1,2) / Tmp
   U (1,3) = U (1,3) / Tmp


   U (2,1) = vr_y * U (1,3) - vr_z * U (1,2)
   U (2,2) = - vr_x * U (1,3) + vr_z * U (1,1)
   U (2,3) = vr_x * U (1,2) - vr_y * U (1,1)

   Tmp = DSQRT ( U (2,1)**2 + U (2,2)**2 + U (2,3)**2 )

   U (2,1) = U (2,1) / Tmp
   U (2,2) = U (2,2) / Tmp
   U (2,3) = U (2,3) / Tmp

   U (3,1) =  vr_x
   U (3,2) =  vr_y
   U (3,3) =  vr_z

   DO uno = 1, 3
    DO dos = 1, 3
       U_t (uno,dos) = U (dos,uno)
    END DO
   END DO


   KILO1 = REFCRD (1,3) - ATMCRD (1,D_ATOM)
   KILO2 = REFCRD (2,3) - ATMCRD (2,D_ATOM)
   KILO3 = REFCRD (3,3) - ATMCRD (3,D_ATOM)

   REFCRD2 (1,1) = REFCRD (1,1) - KILO1
   REFCRD2 (2,1) = REFCRD (2,1) - KILO2
   REFCRD2 (3,1) = REFCRD (3,1) - KILO3
   REFCRD2 (1,2) = REFCRD (1,2) - KILO1
   REFCRD2 (2,2) = REFCRD (2,2) - KILO2
   REFCRD2 (3,2) = REFCRD (3,2) - KILO3
   REFCRD2 (1,3) = REFCRD (1,3) - KILO1
   REFCRD2 (2,3) = REFCRD (2,3) - KILO2
   REFCRD2 (3,3) = REFCRD (3,3) - KILO3
   
   TMPCRD = .0_DP
   TMPCRD2= .0_DP

   DO uno=1,3
      DO dos=1,3
         DO tres=1,3
	    TMPCRD(uno,dos) = TMPCRD (uno,dos) + (U_o (uno,tres) * REFCRD2 (tres,dos))
         END DO
      END DO
   END DO

   vr_x = TMPCRD (1,1) - TMPCRD(1,3)
   vr_y = TMPCRD (2,1) - TMPCRD(2,3)
   vr_z = TMPCRD (3,1) - TMPCRD(3,3)

   Tmp = DSQRT ( vr_x**2 + vr_y**2 + vr_z**2 )

   vr_x = vr_x / Tmp
   vr_y = vr_y / Tmp
   vr_z = vr_z / Tmp

   DO uno=1,3
      DO dos=1,3
         DO tres=1,3
	    TMPCRD2(uno,dos) = TMPCRD2 (uno,dos) + (U_t (uno,tres) * TMPCRD (tres,dos))
         END DO
      END DO
   END DO

   REFCRD2 (1,1) = TMPCRD2 (1,1) + KILO1
   REFCRD2 (2,1) = TMPCRD2 (2,1) + KILO2
   REFCRD2 (3,1) = TMPCRD2 (3,1) + KILO3
   REFCRD2 (1,2) = TMPCRD2 (1,2) + KILO1
   REFCRD2 (2,2) = TMPCRD2 (2,2) + KILO2
   REFCRD2 (3,2) = TMPCRD2 (3,2) + KILO3
   REFCRD2 (1,3) = TMPCRD2 (1,3) + KILO1
   REFCRD2 (2,3) = TMPCRD2 (2,3) + KILO2
   REFCRD2 (3,3) = TMPCRD2 (3,3) + KILO3
 
   vr_x = REFCRD2 (1,1) - REFCRD2(1,3)
   vr_y = REFCRD2 (2,1) - REFCRD2(2,3)
   vr_z = REFCRD2 (3,1) - REFCRD2(3,3)

   Tmp = DSQRT ( vr_x**2 + vr_y**2 + vr_z**2 )

   vr_x = vr_x / Tmp
   vr_y = vr_y / Tmp
   vr_z = vr_z / Tmp

   ve_x = ATMCRD (1,A_ATOM) - ATMCRD(1,D_ATOM)
   ve_y = ATMCRD (2,A_ATOM) - ATMCRD(2,D_ATOM)
   ve_z = ATMCRD (3,A_ATOM) - ATMCRD(3,D_ATOM)

   Tmp = DSQRT ( ve_x**2 + ve_y**2 + ve_z**2 )

   ve_x = ve_x / Tmp
   ve_y = ve_y / Tmp
   ve_z = ve_z / Tmp

   M1X = REFCRD2 (1,2) - REFCRD2(1,1)
   M1Y = REFCRD2 (2,2) - REFCRD2(2,1)
   M1Z = REFCRD2 (3,2) - REFCRD2(3,1)

   M2X = REFCRD2 (1,2) - REFCRD2(1,3)
   M2Y = REFCRD2 (2,2) - REFCRD2(2,3)
   M2Z = REFCRD2 (3,2) - REFCRD2(3,3)
              
   RAB = DSQRT ( M1X**2 + M1Y**2 + M1Z **2)
   RBC = DSQRT ( M2X**2 + M2Y**2 + M2Z **2)

   RC_OLD = RBC - RAB

   B1 = (M1X/RAB)
   B2 = (M1Y/RAB)
   B3 = (M1Z/RAB)
   B4 = (M2X/RBC - M1X/RAB)
   B5 = (M2Y/RBC - M1Y/RAB)
   B6 = (M2Z/RBC - M1Z/RAB)
   B7 = - (M2X/RBC)
   B8 = - (M2Y/RBC)
   B9 = - (M2Z/RBC)

   F11= ( ((B1**2)/ATMMAS(A_ATOM)) + ((B2**2)/ATMMAS(A_ATOM)) + ((B3**2)/ATMMAS(A_ATOM)) +& 
          ((B4**2)/ATMMAS(T_ATOM)) + ((B5**2)/ATMMAS(T_ATOM)) + ((B6**2)/ATMMAS(T_ATOM)) +&
	  ((B7**2)/ATMMAS(D_ATOM)) + ((B8**2)/ATMMAS(D_ATOM)) + ((B9**2)/ATMMAS(D_ATOM)) )
              
   F1 = 1._DP / F11

   DELTARC = B1 * (ATMCRD (1,A_ATOM) - REFCRD2 (1,1)) + B2 * (ATMCRD(2,A_ATOM) - REFCRD2 (2,1)) + &
             B3 * (ATMCRD(3,A_ATOM) - REFCRD2(3,1)) + B4 * (ATMCRD (1,T_ATOM) - REFCRD2 (1,2)) + &
             B5 * (ATMCRD(2,T_ATOM) - REFCRD2 (2,2)) + B6 * (ATMCRD(3,T_ATOM) - REFCRD2(3,2)) +&
             B7 * (ATMCRD (1,D_ATOM) - REFCRD2 (1,3)) + B8 * (ATMCRD(2,D_ATOM) - REFCRD2 (2,3)) + &
             B9 * (ATMCRD(3,D_ATOM) - REFCRD2(3,3)) 

   ITERSHAKE = ITERSHAKE + 1

   M1X = - ATMCRD (1,T_ATOM) + ATMCRD(1,A_ATOM)
   M1Y = - ATMCRD (2,T_ATOM) + ATMCRD(2,A_ATOM)
   M1Z = - ATMCRD (3,T_ATOM) + ATMCRD(3,A_ATOM)

   M2X = - ATMCRD (1,T_ATOM) + ATMCRD(1,D_ATOM)
   M2Y = - ATMCRD (2,T_ATOM) + ATMCRD(2,D_ATOM)
   M2Z = - ATMCRD (3,T_ATOM) + ATMCRD(3,D_ATOM)

   RAB = DSQRT ( M1X**2 + M1Y**2 + M1Z **2)
   RBC = DSQRT ( M2X**2 + M2Y**2 + M2Z **2)

   RC_NEW = RBC - RAB

   ATMCRD(1,A_ATOM) = ATMCRD(1,A_ATOM) - (DELTARC * F1 * B1)/ ATMMAS (A_ATOM)
   ATMCRD(2,A_ATOM) = ATMCRD(2,A_ATOM) - (DELTARC * F1 * B2)/ ATMMAS (A_ATOM)
   ATMCRD(3,A_ATOM) = ATMCRD(3,A_ATOM) - (DELTARC * F1 * B3)/ ATMMAS (A_ATOM)
   ATMCRD(1,T_ATOM) = ATMCRD(1,T_ATOM) - (DELTARC * F1 * B4)/ ATMMAS (T_ATOM)
   ATMCRD(2,T_ATOM) = ATMCRD(2,T_ATOM) - (DELTARC * F1 * B5)/ ATMMAS (T_ATOM)
   ATMCRD(3,T_ATOM) = ATMCRD(3,T_ATOM) - (DELTARC * F1 * B6)/ ATMMAS (T_ATOM)
   ATMCRD(1,D_ATOM) = ATMCRD(1,D_ATOM) - (DELTARC * F1 * B7)/ ATMMAS (D_ATOM)
   ATMCRD(2,D_ATOM) = ATMCRD(2,D_ATOM) - (DELTARC * F1 * B8)/ ATMMAS (D_ATOM)
   ATMCRD(3,D_ATOM) = ATMCRD(3,D_ATOM) - (DELTARC * F1 * B9)/ ATMMAS (D_ATOM)

!   FATM(1,A_ATOM) =  (DELTARC * F1 * B1**2)/ (XACR2*ATMMAS (A_ATOM))
!   FATM(2,A_ATOM) =  (DELTARC * F1 * B2**2)/ (XACR2*ATMMAS (A_ATOM))
!   FATM(3,A_ATOM) =  (DELTARC * F1 * B3**2)/ (XACR2*ATMMAS (A_ATOM))
!   FATM(1,T_ATOM) =  (DELTARC * F1 * B4**2)/ (XACR2*ATMMAS (T_ATOM))
!   FATM(2,T_ATOM) =  (DELTARC * F1 * B5**2)/ (XACR2*ATMMAS (T_ATOM))
!   FATM(3,T_ATOM) =  (DELTARC * F1 * B6**2)/ (XACR2*ATMMAS (T_ATOM))
!   FATM(1,D_ATOM) =  (DELTARC * F1 * B7**2)/ (XACR2*ATMMAS (D_ATOM))
!   FATM(2,D_ATOM) =  (DELTARC * F1 * B8**2)/ (XACR2*ATMMAS (D_ATOM))
!   FATM(3,D_ATOM) =  (DELTARC * F1 * B9**2)/ (XACR2*ATMMAS (D_ATOM))
!
!   FATO = FATM(1,A_ATOM)+FATM(2,A_ATOM)+FATM(3,A_ATOM)+FATM(1,T_ATOM)+FATM(2,T_ATOM)+&
!          FATM(3,T_ATOM)+FATM(1,D_ATOM)+FATM(2,D_ATOM)+FATM(3,D_ATOM)

   IF (ITERSHAKE .GT. 1000) THEN
      WRITE (*,*) 'Excesive number of iterations in shake '
      stop
   end if

   IF( ISTEP .NE. 1 ) THEN
       ATMVEL(1,A_ATOM) = ATMVEL(1,A_ATOM) - FACMIO*((DELTARC * F1 * B1)/ ATMMAS (A_ATOM))
       ATMVEL(2,A_ATOM) = ATMVEL(2,A_ATOM) - FACMIO*((DELTARC * F1 * B2)/ ATMMAS (A_ATOM))
       ATMVEL(3,A_ATOM) = ATMVEL(3,A_ATOM) - FACMIO*((DELTARC * F1 * B3)/ ATMMAS (A_ATOM))
       ATMVEL(1,T_ATOM) = ATMVEL(1,T_ATOM) - FACMIO*((DELTARC * F1 * B4)/ ATMMAS (T_ATOM))
       ATMVEL(2,T_ATOM) = ATMVEL(2,T_ATOM) - FACMIO*((DELTARC * F1 * B5)/ ATMMAS (T_ATOM))
       ATMVEL(3,T_ATOM) = ATMVEL(3,T_ATOM) - FACMIO*((DELTARC * F1 * B6)/ ATMMAS (T_ATOM))
       ATMVEL(1,D_ATOM) = ATMVEL(1,D_ATOM) - FACMIO*((DELTARC * F1 * B7)/ ATMMAS (D_ATOM))
       ATMVEL(2,D_ATOM) = ATMVEL(2,D_ATOM) - FACMIO*((DELTARC * F1 * B8)/ ATMMAS (D_ATOM))
       ATMVEL(3,D_ATOM) = ATMVEL(3,D_ATOM) - FACMIO*((DELTARC * F1 * B9)/ ATMMAS (D_ATOM))
   END IF
   
END SUBROUTINE GH_VV_B


SUBROUTINE GH_VV_C( ATMTMP )
   IMPLICIT NONE
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT) :: ATMTMP

   IF( .NOT. GH_ON ) RETURN

   ATMTMP(1,A_ATOM) = ATMTMP(1,A_ATOM) - FACMIO*((DELTARC * F1 * B1)/ ATMMAS (A_ATOM))
   ATMTMP(2,A_ATOM) = ATMTMP(2,A_ATOM) - FACMIO*((DELTARC * F1 * B2)/ ATMMAS (A_ATOM))
   ATMTMP(3,A_ATOM) = ATMTMP(3,A_ATOM) - FACMIO*((DELTARC * F1 * B3)/ ATMMAS (A_ATOM))
   ATMTMP(1,T_ATOM) = ATMTMP(1,T_ATOM) - FACMIO*((DELTARC * F1 * B4)/ ATMMAS (T_ATOM))
   ATMTMP(2,T_ATOM) = ATMTMP(2,T_ATOM) - FACMIO*((DELTARC * F1 * B5)/ ATMMAS (T_ATOM))
   ATMTMP(3,T_ATOM) = ATMTMP(3,T_ATOM) - FACMIO*((DELTARC * F1 * B6)/ ATMMAS (T_ATOM))
   ATMTMP(1,D_ATOM) = ATMTMP(1,D_ATOM) - FACMIO*((DELTARC * F1 * B7)/ ATMMAS (D_ATOM))
   ATMTMP(2,D_ATOM) = ATMTMP(2,D_ATOM) - FACMIO*((DELTARC * F1 * B8)/ ATMMAS (D_ATOM))
   ATMTMP(3,D_ATOM) = ATMTMP(3,D_ATOM) - FACMIO*((DELTARC * F1 * B9)/ ATMMAS (D_ATOM))

!   write (89,*) FATO*F1

END SUBROUTINE GH_VV_C


SUBROUTINE GH_DEFINE( ACCEPTOR, DONNOR, TRANSFERRED, FILE )
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: ACCEPTOR, DONNOR, TRANSFERRED
   CHARACTER( LEN=* ), INTENT(IN) :: FILE

   ! ----------------------------------------------------------------------
   ! Modifica los numeros segun tu sistema. 
   ! ----------------------------------------------------------------------
   ! A_ATOM = atomo aceptor
   ! D_ATOM = atomo dador
   ! T_ATOM = atomo que se transfiere

   A_ATOM = ACCEPTOR
   D_ATOM = DONNOR
   T_ATOM = TRANSFERRED

   GH_ON = .TRUE.

   GH_UNIT = NEXT_UNIT()
   OPEN( UNIT = GH_UNIT, FILE = TRIM( FILE ), ACTION = "WRITE", FORM = "FORMATTED" )

END SUBROUTINE GH_DEFINE

#endif

END 
