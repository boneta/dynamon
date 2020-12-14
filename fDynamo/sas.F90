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
!                     Solvent Accessible Surface Area Module
!===============================================================================
!
! . Subroutines:
!
!   SURFACE_CALCULATE         Calculate the surface area and its derivatives.
!
! . Notes:
!
!   Calculate the area of a selection of atoms and its derivatives
!   with respect to the atom coordinates. The routine is a modification
!   of Richmond's original program for DYNAMO.
!
!   There is only one input argument. It is optional and labelled and
!   is the probe radius. It takes the default of 1.4 Angstroms.
!
! . Notes :
!
!   MAXARC is the maximum number of partial arcs for the IRth sphere.
!   MAXEPT is the maximum number of overlap end points for a circle of
!          intersection.
!   MAXOVR is the maximum number of overlapping spheres.
!
!   The routine calculates the areas of the atoms and the derivatives
!   of the dot product of the atom areas and their weights.
!
!===============================================================================
MODULE SOLVENT_ACCESSIBLE_SURFACE

! . Module declarations.
USE CONSTANTS,    ONLY : PI
USE DEFINITIONS,  ONLY : DP
USE PRINTING,     ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH

USE ATOMS,        ONLY : ATMCRD, NATOMS

IMPLICIT NONE
PRIVATE
PUBLIC :: SURFACE_CALCULATE

! . Module parameters.
REAL ( KIND = DP ), PARAMETER :: DEFAULT_PROBE_RADIUS = 1.4_DP

!===============================================================================
CONTAINS
!===============================================================================

   !-------------------------------------------------------------------------------------
   SUBROUTINE SURFACE_CALCULATE ( RADIUS, PROBE_RADIUS, SELECTION, WEIGHTS, TOTAL_AREA, &
                                                 AREA_DERIVATIVES, ATOM_SURFACES, PRINT )
   !-------------------------------------------------------------------------------------

   ! . Scalar arguments.
   LOGICAL,            INTENT(IN),  OPTIONAL :: PRINT
   REAL ( KIND = DP ), INTENT(IN),  OPTIONAL :: PROBE_RADIUS
   REAL ( KIND = DP ), INTENT(OUT), OPTIONAL :: TOTAL_AREA

   ! . Array arguments.
   LOGICAL,            DIMENSION(1:NATOMS),     INTENT(IN),  OPTIONAL :: SELECTION
   REAL ( KIND = DP ), DIMENSION(1:NATOMS),     INTENT(IN)            :: RADIUS
   REAL ( KIND = DP ), DIMENSION(1:NATOMS),     INTENT(IN),  OPTIONAL :: WEIGHTS
   REAL ( KIND = DP ), DIMENSION(1:NATOMS),     INTENT(OUT), OPTIONAL :: ATOM_SURFACES
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(OUT), OPTIONAL :: AREA_DERIVATIVES

   ! . Local parameters.
   INTEGER, PARAMETER :: MAXARC = 201, MAXEPT = 300, MAXOVR = 200

   ! . Local scalars.
   INTEGER            :: I, IB, II, IN, IO, IR, J, JB, K, KARC, K1, L, M, MI, N
   LOGICAL            :: ISKIPS, LTOP, QPRINT
   REAL ( KIND = DP ) :: ARCLEN, ARCSUM, AREA, AXX, AXY, AXZ, AYX, AYY,       &
                         AZX, AZY, AZZ, BGL, BK,                              &
                         BSQK, BSQL, CCSQ, CC, DAX, DAY, DAZ, DEAL, DECL, DK, &
                         DSQL, DTKAL, DTKCL, DTLAL, DTLCL, EXANG,             &
                         FACA, FACB, FACC, GACA, GACB, GK, GL,                &
                         P, PID2, PIX2, PIX4, PROBE,                          &
                         RCN, RIK, RISQK, RISQL, RMINUS, RPLUS, RR, RRSQ,     &
                         RRX2, S, SIG, SIGSQ,                                 &
                         T, TB, TD, TF, THE, THERK, TI, TK1, TK2, TOTSUR, TR, &
                         TT, TX, TXB, TXK, TXL, TXR, TY, TYB, TYK, TYL, TYR,  &
                         TZ, TZK, TZL, T1, T2, UXL, UYL, UZL, V,              &
                         WGT, WXL, WXLSQ, XR, XYSQ, YR, ZR


   ! . Atoms arrays.
   LOGICAL,            DIMENSION(1:NATOMS)     :: QSELCT
   REAL ( KIND = DP ), DIMENSION(1:NATOMS)     :: ATMRAD, ATMWGT
   REAL ( KIND = DP ), DIMENSION(1:NATOMS)     :: ATMSAS
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: ATMSAD

   ! . Arc arrays.
   INTEGER, DIMENSION(1:MAXARC) :: KENT, KOUT

   ! . End point arrays.
   INTEGER,            DIMENSION(1:MAXEPT) :: LT
   REAL ( KIND = DP ), DIMENSION(1:MAXEPT) :: ARCF, ARCI, EX

   ! . Overlap arrays.
   INTEGER,            DIMENSION(1:MAXOVR) :: IDER, INTAG, INTAG1, ITAG, SYDER
   LOGICAL,            DIMENSION(1:MAXOVR) :: ISKIP
   REAL ( KIND = DP ), DIMENSION(1:MAXOVR) :: B, BG, BSQ, BSQ1, B1, DSQ, DSQ1, GR, &
                                              RI, RISQ, THER, UX,  UY,  UZ,        &
                                              XC,  YC,  ZC, XC1, YC1, ZC1

   !----------------------------------------------------------------------------
   ! . Initialization.
   !----------------------------------------------------------------------------
   ! . Check the number of atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Get the probe radius.
   IF ( PRESENT ( PROBE_RADIUS ) ) THEN
      PROBE = PROBE_RADIUS
   ELSE
      PROBE = DEFAULT_PROBE_RADIUS
   END IF

   ! . Set the atom radii.
   ATMRAD = RADIUS + PROBE

   ! . Get the atom selection.
   IF ( PRESENT ( SELECTION ) ) THEN
      QSELCT = SELECTION
   ELSE
      QSELCT = .TRUE.
   END IF

   ! . Get the atom weights.
   IF ( PRESENT ( WEIGHTS ) ) THEN
      ATMWGT = WEIGHTS
   ELSE
      ATMWGT = 1.0_DP
   END IF

   !----------------------------------------------------------------------------
   ! . Do the calculation.
   !----------------------------------------------------------------------------
   ! . Initialization.

   ! . Constants.
   PIX2  = 2.0_DP * PI
   PIX4  = 4.0_DP * PI
   PID2  = PI / 2.0_DP

   ! . Overlap significance (also used to test if the spheres are colinear).
   SIG   = 1.0E-5_DP
   SIGSQ = SIG * SIG

   ! . Overlap arrays.
   IDER  = 0
   SYDER = 0

   ! . The surface areas and their derivatives.
   ATMSAD = 0.0_DP
   ATMSAS = 0.0_DP

   ! . Calculate the areas.

   ! . Outer loop over atoms.
   DO IR = 1,NATOMS

      ! . Check to see if the atom has been selected.
         IF ( .NOT. QSELCT(IR) ) GO TO 160
      ! . Local initialization.
      ! . Integers.
         IB = 0
         IO = 1
         JB = 0
      ! . Reals.
         AREA   = 0.0_DP
         ARCLEN = 0.0_DP
         EXANG  = 0.0_DP
      ! . Atom radius.
         RR   = ATMRAD(IR)
         RRSQ = RR * RR
         RRX2 = 2.0_DP * RR
      ! . Atom coordinates.
         XR   = ATMCRD(1,IR)
         YR   = ATMCRD(2,IR)
         ZR   = ATMCRD(3,IR)
      ! . Atom weight.
         WGT  = ATMWGT(IR)
      ! . Inner loop over atoms.
         DO IN = 1,NATOMS
      ! . Check to see if the atom has been selected.
            IF ( .NOT. QSELCT(IN) ) GO TO 10
      ! . Check to see if the IN sphere is next to the IR sphere.
            RPLUS = RR + ATMRAD(IN)
            TX = ATMCRD(1,IN) - XR
            IF ( ABS ( TX ) >= RPLUS ) GO TO 10
            TY = ATMCRD(2,IN) - YR
            IF ( ABS ( TY ) >= RPLUS ) GO TO 10
            TZ = ATMCRD(3,IN) - ZR
            IF ( ABS ( TZ ) >= RPLUS ) GO TO 10
            IF ( IN == IR ) GO TO 10
      ! . Check for overlap of spheres by testing centre to centre distance
      ! . against the sum and the difference of the radii.
            XYSQ = TX**2 + TY**2
            IF ( XYSQ < SIGSQ ) THEN
               TX = SIG
               TY = 0.0_DP
               XYSQ = SIGSQ
            END IF
            CCSQ = XYSQ + TZ**2
            CC   = SQRT ( CCSQ )
            IF ( ( RPLUS - CC ) <= SIG ) GO TO 10
            RMINUS = RR - ATMRAD(IN)
      ! . Check to see if the IR sphere is completely buried.
            IF ( ( CC - ABS ( RMINUS ) ) <= SIG ) THEN
               IF ( RMINUS <= 0.0_DP ) GO TO 150
      ! . Calculate the overlap parameters.
            ELSE
               XC1(IO) = TX
               YC1(IO) = TY
               ZC1(IO) = TZ
               DSQ1(IO) = XYSQ
               BSQ1(IO) = CCSQ
               B1(IO) = CC
               GR(IO) = ( CCSQ + RPLUS * RMINUS ) / ( RRX2 * B1(IO) )
               INTAG1(IO) = IN
               IO = IO + 1
               IF ( IO > MAXOVR ) CALL PRINT_ERROR ( "SURFACE_CALCULATE", "Number of overlaps overflow." )
            END IF
      ! . End of inner loop over atoms.
   10       CONTINUE
         ENDDO
      ! . Adjust the number of overlaps.
         IO = IO - 1
      ! . Zero overlaps.
         IF ( IO == 0 ) THEN
            AREA = PIX4
            GO TO 140
      ! . One overlap.
         ELSE IF ( IO == 1 ) THEN
      ! . Assign some variables.
      ! . Integers.
            INTAG(1) = INTAG1(1)
            K        = 1
      ! . Reals.
            BK   = B1(1)
            BSQK = BSQ1(1)
            TXK  = XC1(1)
            TYK  = YC1(1)
            TZK  = ZC1(1)
      ! . Code from later. Included here to avoid an unpredictable jump into a DO loop.
            ARCSUM = PIX2
            IB = IB + 1
            ARCLEN = ARCLEN + ( GR(K) * ARCSUM )
            IN = INTAG(K)
            T1 = ARCSUM*RRSQ* ( BSQK - RRSQ+ATMRAD(IN)**2 ) / ( RRX2 * BSQK * BK )
      ! . Add contributions into the derivatives.
            ATMSAD(1,IR) = ATMSAD(1,IR) - ( WGT * TXK * T1 )
            ATMSAD(2,IR) = ATMSAD(2,IR) - ( WGT * TYK * T1 )
            ATMSAD(3,IR) = ATMSAD(3,IR) - ( WGT * TZK * T1 )
            ATMSAD(1,IN) = ATMSAD(1,IN) + ( WGT * TXK * T1 )
            ATMSAD(2,IN) = ATMSAD(2,IN) + ( WGT * TYK * T1 )
            ATMSAD(3,IN) = ATMSAD(3,IN) + ( WGT * TZK * T1 )
            GO TO 130
         END IF
      ! . Multiple overlaps.
      ! . Sort the IN spheres by their degree of overlap with the IR sphere.
         CALL SORTAG ( GR, IO, ITAG )
      ! . Loop over the overlaps.
         DO L = 1,IO
            K  = ITAG(L)
            IN = INTAG1(K)
            INTAG(L) = IN
            XC(L) = XC1(K)
            YC(L) = YC1(K)
            ZC(L) = ZC1(K)
            DSQ(L) = DSQ1(K)
            B(L)   = B1(K)
            BSQ(L) = BSQ1(K)
            ISKIP(L) = .FALSE.
         END DO
         DO L = 1,IO
            GL = GR(L) * RR
            BG(L)   = B(L) * GL
            RISQ(L) = RRSQ - GL**2
            RI(L)   = SQRT ( RISQ(L) )
            ! . Radius of the INth circle on the surface of the sphere.
            THER(L) = PID2 - ASIN ( GR(L) )
         END DO
      ! . Find the boundary of the inaccessible area on the IRth sphere.
      ! . Outer loop over overlaps.
         DO K = 1,(IO-1)
            ! . Check to see if the overlap is to be skipped.
            IF ( ISKIP(K) ) GO TO 30
            ! . Local assignments.
            TXK = XC(K)
            TYK = YC(K)
            TZK = ZC(K)
            BK  = B(K)
            THERK = THER(K)
            K1 = K + 1
            ! . Inner loop over overlaps.
            DO L = (K+1),IO
               ! . Check to see if the overlap is to be skipped.
               IF ( ISKIP(L) ) GO TO 20
               ! . Is L circle intersecting K circle?
               ! . Distance between circle centers and the sum of radii.
               CC = ACOS ( ( TXK*XC(L) + TYK*YC(L) + TZK*ZC(L) ) / ( BK*B(L) ) )
               TD = THERK + THER(L)
               ! . Circles enclose separate regions?
               IF ( CC >= TD ) GO TO 20
               ! . Circle L completely inside circle K?
               IF ( ( CC + THER(L) ) >= THER(K) ) THEN
                  ! . Circles essentially parallel?
                  IF ( CC <= SIG ) THEN
                     ISKIP(L) = .TRUE.
                  ELSE
                  ! . IR sphere completely buried?
                     IF ( ( PIX2 - CC ) <= TD ) GO TO 150
                  END IF
               ELSE
                  ISKIP(L) = .TRUE.
               ENDIF
   20          CONTINUE
            END DO
            ! . End of outer loop over overlaps.
   30       CONTINUE
         END DO
         ! . Find the T value of the circle intersections.
         ! . Outer loop over overlaps.
         DO K = 1,IO
            ! . Check to see if the overlap is to be skipped.
            IF ( ISKIP(K) ) GO TO 90
            ISKIPS   = ISKIP(K)
            ISKIP(K) = .TRUE.
            KARC  = 0
            LTOP  = .FALSE.
            TXK   = XC(K)
            TYK   = YC(K)
            TZK   = ZC(K)
            DK    = SQRT(DSQ(K))
            BSQK  = BSQ(K)
            BK    = B(K)
            GK    = GR(K)*RR
            RISQK = RISQ(K)
            RIK   = RI(K)
            THERK = THER(K)
            ! . Rotation matrix elements
            T1  = TZK / ( BK * DK )
            AXX = TXK * T1
            AXY = TYK * T1
            AXZ = DK  / BK
            AYX = TYK / DK
            AYY = TXK / DK
            AZX = TXK / BK
            AZY = TYK / BK
            AZZ = TZK / BK
            ! . Inner loop over overlaps.
            DO L = 1,IO
               ! . Check to see if the overlap is to be skipped.
               IF ( ISKIP(L) ) GO TO 40
               TXL = XC(L)
               TYL = YC(L)
               TZL = ZC(L)
               ! . Rotate spheres so K vector colinear with z-axis.
               UXL  = TXL*AXX + TYL*AXY - TZL*AXZ
               UYL  = TYL*AYY - TXL*AYX
               UZL  = TXL*AZX + TYL*AZY + TZL*AZZ
               IF ( ACOS ( UZL / B(L) ) >= ( THERK + THER(L) ) ) GO TO 40
               GL   = GR(L) * RR
               DSQL = UXL**2 + UYL**2
               TB   = UZL*GK - BG(L)
               TXB  = UXL*TB
               TYB  = UYL*TB
               TD   = RIK*DSQL
               TR   = SQRT ( MAX ( ( RISQK*DSQL - TB**2 ), 0.0_DP ) )
               TXR  = UXL*TR
               TYR  = UYL*TR
               ! . T values of intersection for the Kth circle.
               TB  = ( TXB + TYR ) / TD
               IF ( ABS ( TB ) > 1.0_DP ) TB = SIGN ( 1.0_DP, TB )
               TK1 = ACOS ( TB )
               IF ( ( TYB - TXR ) < 0.0_DP ) TK1 = PIX2 - TK1
               TB  = ( TXB - TYR ) / TD
               IF ( ABS ( TB ) > 1.0_DP ) TB = SIGN ( 1.0_DP, TB )
               TK2 = ACOS ( TB )
               IF ( ( TYB + TXR ) < 0.0_DP ) TK2 = PIX2 - TK2
               THE = - ACOS ( ( RRSQ*UZL - GK*BG(L) ) / ( RIK*RI(L)*B(L) ) )
               ! . Is TK1 entry or exit point?  check T = 0 point.
               ! . TI is exit point, TF is entry point.
               IF ( ( ( ACOS ( ( UZL*GK - UXL*RIK ) / ( B(L)*RR ) ) - THER(L) ) * ( TK2 - TK1 ) ) > 0.0_DP ) THEN
                  TI = TK1
                  TF = TK2
               ELSE
                  TI = TK2
                  TF = TK1
               ENDIF
               KARC = KARC + 1
               IF ( KARC >= MAXEPT ) CALL PRINT_ERROR ( "SURFACE_CALCULATE", "Overlap end point overflow." )
               IF ( TF <= TI ) THEN
                  ARCF(KARC) = TF
                  ARCI(KARC) = 0.0_DP
                  TF = PIX2
                  LT(KARC) = L
                  EX(KARC) = THE
                  LTOP = .TRUE.
                  KARC = KARC + 1
               END IF
               ARCF(KARC) = TF
               ARCI(KARC) = TI
               LT(KARC) = L
               EX(KARC) = THE
               UX(L) = UXL
               UY(L) = UYL
               UZ(L) = UZL
               ! . End of inner loop over overlaps.
   40          CONTINUE
            END DO
            ISKIP(K) = ISKIPS
            ! . Special case: K circle without intersections?
            IF ( KARC <= 0 ) GO TO 70
            ! . General case: sum up arclength and set connectivity code.
            CALL SORTAG ( ARCI, KARC, ITAG )
            ARCSUM = ARCI(1)
            MI = ITAG(1)
            T  = ARCF(MI)
            N  = MI
            IF ( KARC == 1 ) GO TO 50
         ! . Loop over overlap end points.
            DO J = 2,KARC
               M = ITAG(J)
               IF ( T < ARCI(J) ) THEN
                  ARCSUM = ARCSUM + ARCI(J) - T
                  EXANG  = EXANG + EX(N)
                  JB = JB + 1
                  IF ( JB >= MAXARC ) CALL PRINT_ERROR ( "SURFACE_CALCULATE", "Partial arc overflow." )
                  L = LT(N)
                  IDER(L)  = IDER(L)  + 1
                  SYDER(L) = SYDER(L) + 1
                  KENT(JB) = L*1024 + K
                  L = LT(M)
                  IDER(L)  = IDER(L)  + 1
                  SYDER(L) = SYDER(L) - 1
                  KOUT(JB) = K*1024 + L
               END IF
               TT = ARCF(M)
               IF ( TT >= T ) THEN
                  T = TT
                  N = M
               END IF
            ! . End of loop over overlap end points.
            END DO
   50       ARCSUM = ARCSUM + PIX2 - T
            IF ( .NOT. LTOP ) THEN
               EXANG = EXANG + EX(N)
               JB = JB + 1
               L  = LT(N)
               IDER(L)  = IDER(L)  + 1
               SYDER(L) = SYDER(L) + 1
               KENT(JB) = L*1024 + K
               L  = LT(MI)
               IDER(L)  = IDER(L)  + 1
               SYDER(L) = SYDER(L) - 1
               KOUT(JB) = K*1024 + L
            END IF
            ! . Calculate derivatives.
            DO L = 1,IO
               IF ( IDER(L) == 0 ) GOTO 60
               RCN   = IDER(L) * RRSQ
               IDER(L) = 0
               UZL   = UZ(L)
               GL    = GR(L)*RR
               BGL   = BG(L)
               BSQL  = BSQ(L)
               RISQL = RISQ(L)
               WXLSQ = BSQL - UZL**2
               WXL   = SQRT ( WXLSQ )
               P     = BGL - GK*UZL
               V     = RISQK*WXLSQ - P**2
               V     = SQRT ( MAX ( V, 1.0E-10_DP ) )
               T1    = RR * ( GK * ( BGL - BSQL ) + UZL * ( BGL - RRSQ ) ) / ( V * RISQL * BSQL )
               DEAL  =  - WXL*T1
               DECL  =  - UZL*T1 - RR/V
               DTKAL = ( WXLSQ - P ) / ( WXL * V )
               DTKCL = ( UZL - GK ) / V
               S     = GK*B(L) - GL*UZL
               T1    = 2.0_DP * GK - UZL
               T2    = RRSQ - BGL
               DTLAL = - ( RISQL*WXLSQ*B(L)*T1 - S * ( WXLSQ*T2 + RISQL*BSQL ) ) / ( RISQL*WXL*BSQL*V )
               DTLCL = - ( RISQL*B(L) * ( UZL*T1 - BGL ) - UZL*T2*S ) / ( RISQL*BSQL*V )
               GACA  = RCN * ( DEAL - ( GK*DTKAL - GL*DTLAL ) / RR ) / WXL
               GACB  = GK - UZL * GL / B(L)
               GACB  = GACB * SYDER(L) * RR / WXLSQ
               SYDER(L) = 0
               FACA = UX(L)*GACA - UY(L)*GACB
               FACB = UY(L)*GACA + UX(L)*GACB
               FACC = RCN * ( DECL - ( GK*DTKCL - GL*DTLCL ) / RR )
               DAX  = AXX*FACA - AYX*FACB + AZX*FACC
               DAY  = AXY*FACA + AYY*FACB + AZY*FACC
               DAZ  = AZZ*FACC - AXZ*FACA
               IN   = INTAG(L)
               ATMSAD(1,IR) = ATMSAD(1,IR) + ( WGT * DAX )
               ATMSAD(2,IR) = ATMSAD(2,IR) + ( WGT * DAY )
               ATMSAD(3,IR) = ATMSAD(3,IR) + ( WGT * DAZ )
               ATMSAD(1,IN) = ATMSAD(1,IN) - ( WGT * DAX )
               ATMSAD(2,IN) = ATMSAD(2,IN) - ( WGT * DAY )
               ATMSAD(3,IN) = ATMSAD(3,IN) - ( WGT * DAZ )
   ! . End of derivative loop.
   60          CONTINUE
            END DO
            GO TO 80
   70       ARCSUM = PIX2
            IB = IB + 1
   80       ARCLEN = ARCLEN + ( GR(K) * ARCSUM )
            IN = INTAG(K)
            T1 = ARCSUM*RRSQ* ( BSQK - RRSQ+ATMRAD(IN)**2 ) / ( RRX2 * BSQK * BK )
            ATMSAD(1,IR) = ATMSAD(1,IR) - ( WGT * TXK * T1 )
            ATMSAD(2,IR) = ATMSAD(2,IR) - ( WGT * TYK * T1 )
            ATMSAD(3,IR) = ATMSAD(3,IR) - ( WGT * TZK * T1 )
            ATMSAD(1,IN) = ATMSAD(1,IN) + ( WGT * TXK * T1 )
            ATMSAD(2,IN) = ATMSAD(2,IN) + ( WGT * TYK * T1 )
            ATMSAD(3,IN) = ATMSAD(3,IN) + ( WGT * TZK * T1 )
   ! . End of outer loop over overlaps.
   90       CONTINUE
         END DO
         IF ( ARCLEN == 0.0_DP ) GO TO 150
         IF ( JB == 0 ) GO TO 130
! . Find the number of independent boundaries.
         J = 0
         DO K = 1,JB
            IF ( KOUT(K) == 0 ) GO TO 120
            I = K
   100      N = KOUT(I)
            KOUT(I) = 0
            J = J + 1
            DO II = 1,JB
               IF ( N /= KENT(II) ) GO TO 110
               IF ( II == K ) THEN
                  IB = IB + 1
                  IF ( J == JB ) GO TO 130
                  GO TO 120
               ELSE
                  I = II
                  GO TO 100
               ENDIF
   ! . End of inner loop over boundaries.
   110         CONTINUE
            END DO
   ! . End of outer loop over boundaries.
   120      CONTINUE
         END DO
         IB = IB + 1
         CALL PRINT_ERROR ( "SURFACE_CALCULATE", "Connectivity error." )
   130   AREA = IB * PIX2 + EXANG + ARCLEN
         AREA = MOD ( AREA, PIX4 )
   140   AREA = AREA * RRSQ
   ! . Assign the area.
   150   CONTINUE
         ATMSAS(IR) = AREA
   ! . End of outer loop over atoms.
   160 CONTINUE
   END DO

   ! . Scale the atom areas by the weights.
   ATMSAS = ATMSAS * ATMWGT

   ! . Calculate the weighted total area.
   TOTSUR = SUM ( ATMSAS )

   !----------------------------------------------------------------------------
   ! . Finish up.
   !----------------------------------------------------------------------------
   ! . Set the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .FALSE.
   END IF

   ! . Do some printing.
   IF ( QPRINT ) THEN
      WRITE ( PRINT_LINE, "(A,F16.2,A)" ) "Total solvent accessible surface area = ", TOTSUR, " A^2."
      CALL PRINT_PARAGRAPH
   END IF

   ! . Fill the output variables.
   IF ( PRESENT ( TOTAL_AREA       ) ) TOTAL_AREA       = TOTSUR
   IF ( PRESENT ( AREA_DERIVATIVES ) ) AREA_DERIVATIVES = ATMSAD
   IF ( PRESENT ( ATOM_SURFACES    ) ) ATOM_SURFACES    = ATMSAS

   END SUBROUTINE SURFACE_CALCULATE

   !------------------------------
   SUBROUTINE SORTAG ( A, N, TAG )
   !------------------------------

   ! . Argument declarations.
   INTEGER,                            INTENT(IN)  :: N
   INTEGER,            DIMENSION(1:N), INTENT(OUT) :: TAG(N)
   REAL ( KIND = DP ), DIMENSION(1:N), INTENT(OUT) :: A(N)

   ! . Local scalars.
   INTEGER  I, IJ, J, K, L, M, TG
   REAL ( KIND = DP )   T, TT

   ! . Local arrays.
   INTEGER  IL(16), IU(16)

   ! . Initialisation.
   DO I = 1,N
      TAG(I) = I
   END DO
   I = 1
   J = N
   M = 1

   ! . Top of loop.
   10 IF ( I >= J ) GOTO 70
   20 K  = I
      IJ = ( J + I ) / 2
      T  = A(IJ)
      IF ( A(I) > T ) THEN
         A(IJ)  =  A(I)
         A(I)   = T
         T      = A(IJ)
         TG     = TAG(IJ)
         TAG(IJ) = TAG(I)
         TAG(I)  = TG
      ENDIF
      L = J
      IF ( A(J) >= T ) GOTO 40
      A(IJ)   = A(J)
      A(J)    = T
      T       = A(IJ)
      TG      = TAG(IJ)
      TAG(IJ) = TAG(J)
      TAG(J)  = TG
      IF ( A(I) <= T ) GOTO 40
      A(IJ)   = A(I)
      A(I)    = T
      T       = A(IJ)
      TG      = TAG(IJ)
      TAG(IJ) = TAG(I)
      TAG(I)  = TG
      GOTO 40
   30 A(L)   = A(K)
      A(K)   = TT
      TG     = TAG(L)
      TAG(L) = TAG(K)
      TAG(K) = TG
   40 L = L - 1
      IF ( A(L) > T ) GOTO 40
      TT = A(L)
   50 K  = K + 1
      IF ( A(K) < T ) GOTO 50
      IF ( K <= L) GOTO 30
      IF ( ( L - I ) <= ( J - K ) ) GOTO 60
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GOTO 80
   60 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GOTO 80
   70 M = M - 1
   ! . Check for exit.
      IF ( M == 0 ) RETURN
      I = IL(M)
      J = IU(M)
   80 IF ( ( J - I ) >= 1 ) GOTO 20
      IF ( I == 1 ) GOTO 10
      I = I - 1
   90 I = I + 1
      IF ( I == J ) GOTO 70
      T = A(I+1)
      IF ( A(I) <= T ) GOTO 90
      TG = TAG(I+1)
      K  = I
  100 A(K+1)   = A(K)
      TAG(K+1) = TAG(K)
      K = K - 1
      IF ( T < A(K) ) GOTO 100
      A(K+1)   = T
      TAG(K+1) = TG
      GO TO 90

   END SUBROUTINE SORTAG

END MODULE SOLVENT_ACCESSIBLE_SURFACE
