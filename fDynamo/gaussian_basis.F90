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
!                           The Gaussian Basis Module
!===============================================================================
!
! . The Gaussian Basis Function Data Structure:
!
!     BFLAB         the basis function names.
!
!     EX            holds the Gaussian exponents.
!     CS            holds the s-coefficients for the primitive shells.
!     CP            holds the p-coefficients for the primitive shells.
!     CD            holds the d-coefficients for the primitive shells.
!     CF            holds the f-coefficients for the primitive shells.
!     CG            holds the g-coefficients for the primitive shells.
!
!     KATOM         the QM numbering of the atom with the primitive shell.
!     KLOC          contains the starting atomic orbital basis function
!                   number for the primitive shells.
!     KMIN, KMAX    the starting and finishing number of the basis function
!                   types within each shell.
!     KNG           contains the number of primitive Gaussians in the
!                   primitive shells.
!     KSTART        the index to the first primitive Gaussian in each shell.
!     KTYPE         contains the maximum angular quantum number for the
!                   primitive shells.
!     PIKATOM       the KATOM array for shells on PI atoms.
!     NPISHELL      the number of primitive shells on PI atoms.
!     NSHELL        the number of primitive shells.
!
! . Quadrature Parameters:
!
!      H, W         used by GAUSSIAN_OVERLAP and GAUSSIAN_OVERLAP_DERIVATIVES.
!
! . Subroutines:
!
!      GAUSSIAN_INITIALIZE            initialize the basis set data structure.
!      GAUSSIAN_OVERLAP               calculate the overlap matrix, S.
!      GAUSSIAN_OVERLAP_DERIVATIVES   calculate the derivatives of S.
!      GAUSSIAN_SETUP                 set up the basis functions.
!
! . Notes:
!
!   In the case that PI atoms are present, the basis functions for only a
!   single copy of the PI atoms are stored. They are automatically the last
!   basis functions in the list.
!
!===============================================================================
MODULE GAUSSIAN_BASIS

! . Module declarations.
USE CONSTANTS,        ONLY : ANGSTROMS_TO_BOHRS
USE DEFINITIONS,      ONLY : DP
USE ELEMENTS,         ONLY : SYMBOL
USE PRINTING,         ONLY : PRINT_ERROR
USE STRING,           ONLY : ENCODE_INTEGER

USE ATOMS,            ONLY : ATMCRD, ATMNUM, ATMQMI, NATOMS, NATOMSQM
USE MOPAC_DATA,       ONLY : BFINDEX, MOPAC_DATA_ATOM_POSITION, NBASIS, NBASTR, NPIATOMS, NPIBEADS, PILIST, QMATOM
USE MOPAC_PARAMETERS, ONLY : NATORB, ZP, ZS

IMPLICIT NONE
PUBLIC
#ifndef PGPC
SAVE
#endif

! . Some real constants used in the integrals routines.
REAL ( KIND = DP ), PARAMETER :: PI32 = 5.56832799683170_DP, RLN10 = 2.30258_DP, SQRT3 = 1.73205080756888_DP

! . Scalar basis function data.
INTEGER :: NPISHELL = 0, NSHELL = 0

! . Array basis function data.
CHARACTER ( LEN = 12 ), ALLOCATABLE, DIMENSION(:)   :: BFLAB
INTEGER,                ALLOCATABLE, DIMENSION(:)   :: KSTART, KATOM, KTYPE, KNG, KLOC, KMIN, KMAX
INTEGER,                ALLOCATABLE, DIMENSION(:,:) :: PIKATOM
REAL ( KIND = DP ),     ALLOCATABLE, DIMENSION(:)   :: EX, CS, CP, CD, CF, CG

! . Quadrature parameters.
REAL ( KIND = DP ), DIMENSION(1:21), PARAMETER :: H =  (/ 0.0_DP,               -0.707106781186548_DP, &
                                                          0.707106781186548_DP, -1.22474487139159_DP,  &
                                                          0.0_DP,                1.22474487139159_DP,  &
                                                         -1.65068012388578_DP,  -0.524647623275290_DP, &
                                                          0.524647623275290_DP,  1.65068012388578_DP,  &
                                                         -2.02018287045609_DP,  -0.958572464613819_DP, &
                                                          0.0_DP,                0.958572464613819_DP, &
                                                          2.02018287045609_DP,  -2.350604973674_DP,    &
                                                         -1.335849074014_DP,    -0.436077411928_DP,    &
                                                          0.436077411928_DP,     1.335849074014_DP,    &
                                                          2.350604973674_DP /)

REAL ( KIND = DP ), DIMENSION(1:21), PARAMETER :: W =  (/ 1.77245385090552_DP,   0.8862269254528_DP,    &
                                                          0.8862269254528_DP,    0.2954089751509_DP,    &
                                                          1.181635900604_DP,     0.2954089751509_DP,    &
                                                          8.131283544725E-02_DP, 8.049140900055E-01_DP, &
                                                          8.049140900055E-01_DP, 8.131283544725E-02_DP, &
                                                          1.995324205905E-02_DP, 3.936193231522E-01_DP, &
                                                          9.453087204829E-01_DP, 3.936193231522E-01_DP, &
                                                          1.995324205905E-02_DP, 4.530009905509E-03_DP, &
                                                          1.570673203229E-01_DP, 7.246295952244E-01_DP, &
                                                          7.246295952244E-01_DP, 1.570673203229E-01_DP, &
                                                          4.530009905509E-03_DP /)

!===============================================================================
CONTAINS
!===============================================================================

   !-----------------------------
   SUBROUTINE GAUSSIAN_INITIALIZE
   !-----------------------------

   ! . Scalar basis function data.
   NPISHELL = 0
   NSHELL   = 0 

   ! . Array basis function data.
   IF ( ALLOCATED ( BFLAB   ) ) DEALLOCATE ( BFLAB   )
   IF ( ALLOCATED ( KSTART  ) ) DEALLOCATE ( KSTART  )
   IF ( ALLOCATED ( KATOM   ) ) DEALLOCATE ( KATOM   )
   IF ( ALLOCATED ( KTYPE   ) ) DEALLOCATE ( KTYPE   )
   IF ( ALLOCATED ( KNG     ) ) DEALLOCATE ( KNG     )
   IF ( ALLOCATED ( KLOC    ) ) DEALLOCATE ( KLOC    )
   IF ( ALLOCATED ( KMIN    ) ) DEALLOCATE ( KMIN    )
   IF ( ALLOCATED ( KMAX    ) ) DEALLOCATE ( KMAX    )
   IF ( ALLOCATED ( PIKATOM ) ) DEALLOCATE ( PIKATOM )
   IF ( ALLOCATED ( EX      ) ) DEALLOCATE ( EX      )
   IF ( ALLOCATED ( CS      ) ) DEALLOCATE ( CS      )
   IF ( ALLOCATED ( CP      ) ) DEALLOCATE ( CP      )
   IF ( ALLOCATED ( CD      ) ) DEALLOCATE ( CD      )
   IF ( ALLOCATED ( CF      ) ) DEALLOCATE ( CF      )
   IF ( ALLOCATED ( CG      ) ) DEALLOCATE ( CG      )

   END SUBROUTINE GAUSSIAN_INITIALIZE

   !-------------------------------------------------
   SUBROUTINE GAUSSIAN_OVERLAP ( OVERLAP, PIOVERLAP )
   !-------------------------------------------------

   ! . Argument declarations.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT) ::   OVERLAP
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: PIOVERLAP

   ! . ONE_ELECTRON.DAT declaration.
   INTEGER            :: NI, NJ
   REAL ( KIND = DP ) :: XINT, YINT, ZINT, T, X0, Y0, Z0, XI, YI, ZI, XJ, YJ, ZJ

   ! . Local scalars.
   INTEGER :: IBEAD, NNOTRI

   ! . Calculate the overlap.
   CALL GET_OVERLAP ( OVERLAP, 1, NSHELL - NPISHELL, 0 )

   ! . There are PI atoms.
   IF ( NPIBEADS > 1 ) THEN

      ! . Set the element number counter.
      NNOTRI = SIZE ( OVERLAP )

      ! . Loop over the polymers.
      DO IBEAD = 1,NPIBEADS

         ! . Make sure KATOM is up to date.
         KATOM(NSHELL-NPISHELL+1:NSHELL) = PIKATOM(1:NPISHELL,IBEAD)

         ! . Get the overlap terms.
         CALL GET_OVERLAP ( PIOVERLAP(:,IBEAD), NSHELL - NPISHELL + 1, NSHELL, NNOTRI )
      END DO

   END IF

   !===============================================================================
   CONTAINS
   !===============================================================================

      !--------------------------------------------------------
      SUBROUTINE GET_OVERLAP ( OVERLAP, NSTART, NSTOP, INDINC )
      !--------------------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(IN) :: INDINC, NSTART, NSTOP

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: OVERLAP

      ! . Local scalars.
      INTEGER            :: I, IG, II, IJ, IN, ITOL, I1, I2, J, JG, JGMAX, JJ, JN, &
                            J1, J2, LI, LIT, LJ, LJT, LOCI, LOCJ, MAX, MAXI,       &
                            MAXJ, MINI, MINJ, N, NN, NX, NY, NZ
      LOGICAL            :: DOUBLE, IANDJ
      REAL ( KIND = DP ) :: AA, AA1, AI, AJ, ARRI, AX, AXI, AY, AYI, AZ, AZI,   &
                            CDI, CDJ, CPI, CPJ, CSI, CSJ, DUM, DUM1, DUM2, FAC, &
                            RR, TOL, T1, T2, YZ
   
      ! . Local arrays.
      INTEGER, DIMENSION(1:10) :: IX, IY, IZ, JX, JY, JZ
      INTEGER, DIMENSION(1:36) :: IJX, IJY, IJZ
      REAL ( KIND = DP ), DIMENSION(1:3)  :: RI, RJ
      REAL ( KIND = DP ), DIMENSION(1:36) :: DIJ, S
      REAL ( KIND = DP ), DIMENSION(1:27) :: XIN, YIN, ZIN

      ! . Local data declarations.
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0/
      DATA IX / 1, 4, 1, 1, 7, 1, 1, 4, 4, 1/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1/
      DATA IY / 1, 1, 4, 1, 1, 7, 1, 4, 1, 4/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1/
      DATA IZ / 1, 1, 1, 4, 1, 1, 7, 1, 4, 4/

      ! . Set up some constants.
      ITOL = 15
      TOL  = RLN10 * REAL ( ITOL, DP )

      ! . Initialise the one-electron matrix.
      OVERLAP = 0.0_DP

      ! . Calculate the S matrix.
      ! . ISHELL.
      DO 9000 II = NSTART,NSTOP

      ! . Get information about the shell.
      RI = MOPAC_DATA_ATOM_POSITION ( KATOM(II), ATMCRD ) * ANGSTROMS_TO_BOHRS
      XI = RI(1)
      YI = RI(2)
      ZI = RI(3)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI

      ! . JSHELL.
      DO 8000 JJ = 1,II

      ! . Get information about the shell.
      RJ = MOPAC_DATA_ATOM_POSITION ( KATOM(JJ), ATMCRD ) * ANGSTROMS_TO_BOHRS
      XJ = RJ(1)
      YJ = RJ(2)
      ZJ = RJ(3)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IANDJ=II.EQ.JJ

      ! . Prepare the indices for pairs of (I,J) functions.
      IJ=0
      MAX=MAXJ
      DO I=MINI,MAXI
         NX=IX(I)
         NY=IY(I)
         NZ=IZ(I)
         IF(IANDJ) MAX=I
         DO J=MINJ,MAX
            IJ=IJ+1
            IJX(IJ)=NX+JX(J)
            IJY(IJ)=NY+JY(J)
            IJZ(IJ)=NZ+JZ(J)
         END DO
      END DO

      DO 60 I=1,IJ
      S(I)=0.0_DP
      60 CONTINUE
      ! . I PRIMITIVE.
      JGMAX=J2
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      ! . J PRIMTIVE.
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=1.0_DP/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC=DEXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
      ! . Density factor.
      DOUBLE=IANDJ.AND.IG.NE.JG
      MAX=MAXJ
      NN=0
      DO I=MINI,MAXI
         GO TO (70,80,110,110,90,110,110,100,110,110),I
         70 DUM1=CSI*FAC
         GO TO 110
         80 DUM1=CPI*FAC
         GO TO 110
         90 DUM1=CDI*FAC
         GO TO 110
         100 DUM1=DUM1*SQRT3
         110 IF(IANDJ) MAX=I
         DO J=MINJ,MAX
            GO TO (120,130,160,160,140,160,160,150,160,160),J
            120 DUM2=DUM1*CSJ
            IF(.NOT.DOUBLE) GO TO 160
            IF(I.GT.1) GO TO 125
            DUM2=DUM2+DUM2
            GO TO 160
            125 DUM2=DUM2+CSI*CPJ*FAC
            GO TO 160
            130 DUM2=DUM1*CPJ
            IF(DOUBLE) DUM2=DUM2+DUM2
            GO TO 160
            140 DUM2=DUM1*CDJ
            IF(DOUBLE) DUM2=DUM2+DUM2
            GO TO 160
            150 DUM2=DUM2*SQRT3
            160 NN=NN+1
            DIJ(NN)=DUM2
         END DO
      END DO
      ! . Overlap.
      T= SQRT(AA1)
      T1=-2.0_DP*AJ*AJ*T
      T2=-0.5_DP*T
      X0=AX
      Y0=AY
      Z0=AZ
      IN=-3
      DO I=1,LIT
         IN=IN+3
         NI=I
         DO J=1,LJT
            JN=IN+J
            NJ=J
            CALL STVXYZ
            XIN(JN)=XINT*T
            YIN(JN)=YINT*T
            ZIN(JN)=ZINT*T
            NJ=J+2
            CALL STVXYZ
            XIN(JN+9)=XINT*T1
            YIN(JN+9)=YINT*T1
            ZIN(JN+9)=ZINT*T1
            NJ=J-2
            IF(NJ.GT.0) GO TO 320
            XINT=0.0_DP
            YINT=0.0_DP
            ZINT=0.0_DP
            GO TO 330
            320 CALL STVXYZ
            330 N=(J-1)*(J-2)
            DUM= FLOAT(N)*T2
            XIN(JN+18)=XINT*DUM
            YIN(JN+18)=YINT*DUM
            ZIN(JN+18)=ZINT*DUM
         END DO
      END DO
      DO 350 I=1,IJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      YZ=YIN(NY)*ZIN(NZ)
      DUM  = YZ*XIN(NX)
      DUM1 = (XIN(NX+9)+XIN(NX+18))*YZ+(YIN(NY+9)+YIN(NY+18))*XIN(NX)*ZIN(NZ)+(ZIN(NZ+9)+ZIN(NZ+18))*XIN(NX)*YIN(NY)
      S(I)=S(I)+DIJ(I)*DUM
       350 CONTINUE
      6000 CONTINUE
      7000 CONTINUE

      ! . Set up the overlap matrix.
      MAX=MAXJ
      NN=0
      DO I=MINI,MAXI
         LI=LOCI+I
         IN=BFINDEX(LI) - INDINC
         IF(IANDJ) MAX=I
         DO J=MINJ,MAX
            LJ=LOCJ+J
            JN=LJ+IN
            NN=NN+1
            OVERLAP(JN) = S(NN)
         END DO
      END DO

         8000 CONTINUE
      9000 CONTINUE

      END SUBROUTINE GET_OVERLAP

      !----------------
      SUBROUTINE STVXYZ
      !----------------

      ! . Local scalars.
      INTEGER            :: I, IMAX, IMIN, NPTS
      REAL ( KIND = DP ) :: AX, AY, AZ, BX, BY, BZ, DUM, PTX, PTY, PTZ, PX, PY, PZ

      ! . Local arrays.
      INTEGER, DIMENSION(1:6) :: MIN, MAX

      ! . Local data declarations.
      DATA MIN /1, 2, 4, 7, 11, 16/
      DATA MAX /1, 3, 6, 10, 15, 21/

      XINT=0.0_DP
      YINT=0.0_DP
      ZINT=0.0_DP
      NPTS=(NI+NJ-2)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO I=IMIN,IMAX
         DUM=W(I)
         PX=DUM
         PY=DUM
         PZ=DUM
         DUM=H(I)*T
         PTX=DUM+X0
         PTY=DUM+Y0
         PTZ=DUM+Z0
         AX=PTX-XI
         AY=PTY-YI
         AZ=PTZ-ZI
         BX=PTX-XJ
         BY=PTY-YJ
         BZ=PTZ-ZJ
         GO TO (5,4,3,2,1),NI
         1 PX=PX*AX
         PY=PY*AY
         PZ=PZ*AZ
         2 PX=PX*AX
         PY=PY*AY
         PZ=PZ*AZ
         3 PX=PX*AX
         PY=PY*AY
         PZ=PZ*AZ
         4 PX=PX*AX
         PY=PY*AY
         PZ=PZ*AZ
         5 GO TO (12,11,10,9,8,7,6),NJ
         6 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         7 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         8 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         9 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         10 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         11 PX=PX*BX
         PY=PY*BY
         PZ=PZ*BZ
         12 CONTINUE
         XINT=XINT+PX
         YINT=YINT+PY
         ZINT=ZINT+PZ
      END DO

      END SUBROUTINE STVXYZ

   END SUBROUTINE GAUSSIAN_OVERLAP

   !-----------------------------------------------------------------------
   SUBROUTINE GAUSSIAN_OVERLAP_DERIVATIVES ( SX, SY, SZ, PISX, PISY, PISZ )
   !-----------------------------------------------------------------------

   ! . Argument declarations.
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT) ::   SX,   SY,   SZ
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: PISX, PISY, PISZ

   ! . OEI_XYZDER.DAT.
   INTEGER            :: NI, NJ
   REAL ( KIND = DP ) :: T, XI, XINT, XJ, X0, YI, YINT, YJ, Y0, ZI, ZINT, ZJ, Z0

   ! . Local scalars.
   INTEGER :: IBEAD, NNOTRI

   ! . Calculate the overlap.
   CALL GET_OVERLAP_DERIVATIVES ( SX, SY, SZ, 1, NSHELL - NPISHELL, 0 )

   ! . There are PI atoms.
   IF ( NPIBEADS > 1 ) THEN

      ! . Set the element number counter.
      NNOTRI = SIZE ( SX )

      ! . Loop over the polymers.
      DO IBEAD = 1,NPIBEADS

         ! . Make sure KATOM is up to date.
         KATOM(NSHELL-NPISHELL+1:NSHELL) = PIKATOM(1:NPISHELL,IBEAD)

         ! . Get the overlap terms.
         CALL GET_OVERLAP_DERIVATIVES ( PISX(:,IBEAD), PISY(:,IBEAD), PISZ(:,IBEAD), NSHELL - NPISHELL + 1, NSHELL, NNOTRI )
      END DO

   END IF

   !===============================================================================
   CONTAINS
   !===============================================================================

      !-----------------------------------------------------------------------
      SUBROUTINE GET_OVERLAP_DERIVATIVES ( SX, SY, SZ, NSTART, NSTOP, INDINC )
      !-----------------------------------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(IN) :: INDINC, NSTART, NSTOP

      ! . Array arguments.
      REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: SX, SY, SZ

      ! . Local variables.
      INTEGER            :: I, IG, II, IJ, IN, ITOL, I1, I2, J, JG, JJ, JN, J1, J2, LIT, LOCI, LJT, LOCJ, &
                            MAXI, MAXJ, MINI, MINJ, N, NN, NX, NY, NZ, N0
      LOGICAL            :: IANDJ
      REAL ( KIND = DP ) :: AA, AA1, AI, AJ, ARRI, AX, AXI, AY, AYI, AZ, AZI, CDI, CDJ, CFI, CPI, CPJ, CSI, CSJ, &
                            DUM, DUM1, DUM2, FAC, RR, TOL

      ! . Local arrays.
      INTEGER,            DIMENSION(1:78) :: IJG, IJX, IJY, IJZ
      LOGICAL,            DIMENSION(1:20) :: ISKIP
      REAL ( KIND = DP ), DIMENSION(1:3)  :: RI, RJ
      REAL ( KIND = DP ), DIMENSION(1:78) :: DIJ, S
      REAL ( KIND = DP ), DIMENSION(1:36) :: XIN, YIN, ZIN

      ! . Local data declarations.
      INTEGER, DIMENSION(1:40) :: INDX
      INTEGER, DIMENSION(1:20) :: IX, IY, IZ
      INTEGER, DIMENSION(1:10) :: JX, JY, JZ
      DATA INDX / 1, 2, 3, 4, 5, 6, 7, 8, 9,10,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0, &
                  -0, 1, 2, 3,-0,-0,-0,-0,-0,-0,4, 5, 6, 7, 8, 9,10,11,12,13/
      DATA IX   / 1, 4, 1, 1, 7, 1, 1, 4, 4, 1,10, 1, 1, 7, 7, 4, 1, 4, 1, 4/
      DATA IY   / 1, 1, 4, 1, 1, 7, 1, 4, 1, 4,1,10, 1, 4, 1, 7, 7, 1, 4, 4/
      DATA IZ   / 1, 1, 1, 4, 1, 1, 7, 1, 4, 4,1, 1,10, 1, 4, 1, 4, 7, 7, 4/
      DATA JX   / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0/
      DATA JY   / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1/
      DATA JZ   / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1/

      ! . Initialization.
      SX = 0.0_DP ; SY = 0.0_DP ; SZ = 0.0_DP

      ITOL = 15
      TOL  = RLN10 * REAL ( ITOL, DP )

!     ----- ISHELL
      DO 9000 II=NSTART,NSTOP

      ! . Get information about the shell.
      RI = MOPAC_DATA_ATOM_POSITION ( KATOM(II), ATMCRD ) * ANGSTROMS_TO_BOHRS
      XI = RI(1)
      YI = RI(2)
      ZI = RI(3)

      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)+1
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
      ISKIP(1:20)=.TRUE.
      DO 20 I=MINI,MAXI
      GOTO (11,13,20,20,15,20,20,20,20,20),I
   11 ISKIP(2:4)=.FALSE.
      GOTO 20
   13 ISKIP(5:10)=.FALSE.
      ISKIP(1)=.FALSE.
      GOTO 20
   15 ISKIP(2:4)=.FALSE.
      ISKIP(11:20)=.FALSE.
   20 CONTINUE

!     ----- JSHELL
      DO 8000 JJ=1,II

      ! . Get information about the shell.
      RJ = MOPAC_DATA_ATOM_POSITION ( KATOM(JJ), ATMCRD ) * ANGSTROMS_TO_BOHRS
      XJ = RJ(1)
      YJ = RJ(2)
      ZJ = RJ(3)

      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IANDJ=II.EQ.JJ
!     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
      N0=0
      IF(LIT.EQ.4) N0=20
      IJ=0
      DO 60 I=1,20
      IF(ISKIP(I)) GOTO 60
      IN=INDX(I+N0)
      NX=IX(I)
      NY=IY(I)
      NZ=IZ(I)
      DO 50 J=MINJ,MAXJ
      IJ=IJ+1
      IJX(IJ)=NX+JX(J)
      IJY(IJ)=NY+JY(J)
      IJZ(IJ)=NZ+JZ(J)
      IJG(IJ)=IN+13*(J-MINJ)
   50 CONTINUE
   60 CONTINUE
      DO I=1,IJ
         N=IJG(I)
         S(N)=0.0_DP
      END DO
!     ----- I PRIMITIVE
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      DUM=AI+AI
      CSI=CP(IG)
      CPI=CS(IG)*DUM
      IF(LIT.EQ.4) CPI=CD(IG)
      CDI=CP(IG)*DUM
      CFI=CD(IG)*DUM
!     ----- J PRIMTIVE
      DO 6000 JG=J1,J2
      AJ=EX(JG)
      AA =AI+AJ
      AA1=1.0_DP/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GOTO 6000
      FAC=DEXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
!     ----- DENSITY FACTOR
      NN=0
      DO 180 I=1,20
      IF(ISKIP(I)) GOTO 180
      GOTO ( 70, 80,110,110, 90,110,110,110,110,110, 100,110,110,110,110,110,110,110,110,110),I
   70 DUM1=CSI*FAC
      GOTO 110
   80 DUM1=CPI*FAC
      GOTO 110
   90 DUM1=CDI*FAC
      GOTO 110
  100 DUM1=CFI*FAC
  110 CONTINUE
      DO J=MINJ,MAXJ
         GOTO (120,130,160,160,140,160,160,150,160,160),J
  120    DUM2=DUM1*CSJ
         GOTO 160
  130    DUM2=DUM1*CPJ
         GOTO 160
  140    DUM2=DUM1*CDJ
         GOTO 160
  150    DUM2=DUM2*SQRT3
  160    NN=NN+1
         DIJ(NN)=DUM2
      END DO
  180 CONTINUE
!     ----- OVERLAP
      T=SQRT(AA1)
      X0=AX
      Y0=AY
      Z0=AZ
      IN=-3
      DO I=1,LIT
         IN=IN+3
         NI=I
         DO J=1,LJT
            JN=IN+J
            NJ=J
            CALL DERXYZ
            XIN(JN)=XINT*T
            YIN(JN)=YINT*T
            ZIN(JN)=ZINT*T
         END DO
      END DO
      DO 350 I=1,IJ
      N=IJG(I)
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      S(N)=S(N)+DIJ(I)*XIN(NX)*YIN(NY)*ZIN(NZ)
  350 CONTINUE
!     ----- END OF PRIMITIVE LOOPS -----
 6000 CONTINUE
 7000 CONTINUE
!     ----- FORM INTEGRALS OVER DERIVATIVES -----
      NN=0
      N=1
      DO 7301 J=MINJ,MAXJ
      IF(MINI.GT.1) GOTO 7100
      NN=NN+1
      XIN(NN)= S(N+ 1)
      YIN(NN)= S(N+ 2)
      ZIN(NN)= S(N+ 3)
      IF(MAXI.EQ.1) GOTO 7300
 7100 IF(MINI.GT.2) GOTO 7200
      NN=NN+1
      XIN(NN)=(S(N+ 4)-S(N   ))
      YIN(NN)= S(N+ 7)
      ZIN(NN)= S(N+ 8)
      NN=NN+1
      XIN(NN)= S(N+ 7)
      YIN(NN)=(S(N+ 5)-S(N   ))
      ZIN(NN)= S(N+ 9)
      NN=NN+1
      XIN(NN)= S(N+ 8)
      YIN(NN)= S(N+ 9)
      ZIN(NN)=(S(N+ 6)-S(N   ))
      IF(MAXI.EQ.4) GOTO 7300
 7200 CONTINUE
      NN=NN+1
      XIN(NN)=(S(N+ 3)-S(N   )-S(N   ))
      YIN(NN)= S(N+ 6)
      ZIN(NN)= S(N+ 7)
      NN=NN+1
      XIN(NN)= S(N+ 8)
      YIN(NN)=(S(N+ 4)-S(N+ 1)-S(N+ 1))
      ZIN(NN)= S(N+ 9)
      NN=NN+1
      XIN(NN)= S(N+10)
      YIN(NN)= S(N+11)
      ZIN(NN)=(S(N+ 5)-S(N+ 2)-S(N+ 2))
      NN=NN+1
      DUM=SQRT3
      XIN(NN)=DUM*(S(N+ 6)-S(N+ 1))
      YIN(NN)=DUM*(S(N+ 8)-S(N   ))
      ZIN(NN)=DUM* S(N+12)
      NN=NN+1
      XIN(NN)=DUM*(S(N+ 7)-S(N+ 2))
      YIN(NN)=DUM* S(N+12)
      ZIN(NN)=DUM*(S(N+10)-S(N   ))
      NN=NN+1
      XIN(NN)=DUM* S(N+12)
      YIN(NN)=DUM*(S(N+ 9)-S(N+ 2))
      ZIN(NN)=DUM*(S(N+11)-S(N+ 1))
 7300 N=N+13
 7301 CONTINUE

!     ----- CALCULATE DERIVATIVES OF OVERLAP MATRIX -----
      N=0
      DO 7500 J=MINJ,MAXJ
      JN=LOCJ+J
      DO 7400 I=MINI,MAXI
      N=N+1
      IN=LOCI+I
      IF(JN.GT.IN) GOTO 7400
      NN=BFINDEX(IN)+JN - INDINC
      SX(NN)=XIN(N)
      SY(NN)=YIN(N)
      SZ(NN)=ZIN(N)
 7400 CONTINUE
 7500 CONTINUE
 8000 CONTINUE
 9000 CONTINUE
!     ----- END OF SHELL LOOPS -----

      END SUBROUTINE GET_OVERLAP_DERIVATIVES

      !----------------
      SUBROUTINE DERXYZ
      !----------------

      ! . Local variables.
      INTEGER            ::  I, IMIN, IMAX, NPTS
      REAL ( KIND = DP ) :: AX, AY, AZ, BX, BY, BZ, DUM, PTX, PTY, PTZ, PX, PY, PZ

      ! . Local data declarations.
      INTEGER, DIMENSION(1:6) :: MAX, MIN
      DATA MAX / 1, 3, 6, 10, 15, 21/
      DATA MIN / 1, 2, 4,  7, 11, 16/

      XINT=0.0_DP
      YINT=0.0_DP
      ZINT=0.0_DP
      NPTS=(NI+NJ-2)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 13 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)*T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GOTO (5,4,3,2,1),NI
    1 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    2 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    3 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    4 PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
    5 GOTO (12,11,10,9,8,7,6),NJ
    6 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    7 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    8 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
    9 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   10 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   11 PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
   12 CONTINUE
      XINT=XINT+PX
      YINT=YINT+PY
      ZINT=ZINT+PZ
   13 CONTINUE

      END SUBROUTINE DERXYZ

   END SUBROUTINE GAUSSIAN_OVERLAP_DERIVATIVES

   !---------------------------------------------------------------
   SUBROUTINE GAUSSIAN_SETUP ( NHEAVY, NLIGHT, NPIHEAVY, NPILIGHT )
   !---------------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: NHEAVY, NLIGHT, NPIHEAVY, NPILIGHT

   ! . Initialize the basis set data structure.
   CALL GAUSSIAN_INITIALIZE

   ! . Get the basis set.
   CALL GET_BASIS_SET   

   ! . Get the basis function labels.
   CALL GET_BASIS_FUNCTION_LABELS

   ! . Normalize the basis functions.
   CALL NORMALIZE_BASIS_FUNCTIONS

   !============================================================================
   CONTAINS
   !============================================================================

      !-----------------------
      SUBROUTINE GET_BASIS_SET
      !-----------------------

      ! . Local scalars.
      INTEGER            :: I, IATOM, ICOPY, IGAUSS, IQM, LOC, NGAUSS, NI
      REAL ( KIND = DP ) :: EE, FACP, FACS, SCALE, ZETAP, ZETAS

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:6) :: CCP, CCS, EXP, EXS

      ! . Calculate the number of primitives and shells.
      NPISHELL = ( 2 * NPIHEAVY ) + NPILIGHT
      NSHELL   = ( 2 * NHEAVY   ) + NLIGHT
      NGAUSS   = 6 * NSHELL

      ! . Allocate space for the basis function arrays.
      ALLOCATE ( KSTART(1:NSHELL), KATOM(1:NSHELL), KTYPE(1:NSHELL), KNG(1:NSHELL), &
                 KLOC(1:NSHELL),   KMIN(1:NSHELL),  KMAX(1:NSHELL), PIKATOM(1:NPISHELL,1:NPIBEADS) )
      ALLOCATE ( EX(1:NGAUSS), CS(1:NGAUSS), CP(1:NGAUSS), CD(1:NGAUSS), CF(1:NGAUSS), CG(1:NGAUSS) )

      ! . Initialize the basis set coefficient and exponent arrays.
      CD = 0.0_DP
      CF = 0.0_DP
      CG = 0.0_DP
      CP = 0.0_DP
      CS = 0.0_DP
      EX = 0.0_DP

      ! . Initialize some counters.
      LOC    = 0
      NGAUSS = 0
      NSHELL = 0

      ! . Loop over the QM atoms.
      DO IQM = 1,NATOMSQM

         ! . Get some information about the atom.
	 IATOM = QMATOM(IQM)%ATOM
         ICOPY = QMATOM(IQM)%COPY
         NI    = QMATOM(IQM)%NUMBER

         ! . Skip PI atoms not in the first copy.
         IF ( ICOPY > 1 ) CYCLE

         ! . Check for basis functions.
         IF ( NATORB(NI) > 0 ) THEN

         ! . Get the exponent scaling factors.
         ZETAS = ZS(NI)
         ZETAP = ZP(NI)

         ! . Fill the coefficient and exponent arrays.
         SELECT CASE ( NI )
         ! . 1s orbitals.
         CASE ( 1:2 )
            ! . 1s expansion.
            EXS(1) =  2.310303149E+01_DP
            CCS(1) =  9.163596280E-03_DP
            EXS(2) =  4.235915534E+00_DP
            CCS(2) =  4.936149294E-02_DP
            EXS(3) =  1.185056519E+00_DP
            CCS(3) =  1.685383049E-01_DP
            EXS(4) =  4.070988982E-01_DP
            CCS(4) =  3.705627997E-01_DP
            EXS(5) =  1.580884151E-01_DP
            CCS(5) =  4.164915298E-01_DP
            EXS(6) =  6.510953954E-02_DP
            CCS(6) =  1.303340841E-01_DP
         ! . 2s and 2p orbitals.
         CASE ( 3:10 )
            ! . 2s expansion.
            EXS(1) =  2.768496241E+01
            CCS(1) = -4.151277819E-03_DP
            EXS(2) =  5.077140627E+00_DP
            CCS(2) = -2.067024148E-02_DP
            EXS(3) =  1.426786050E+00_DP
            CCS(3) = -5.150303337E-02_DP
            EXS(4) =  2.040335729E-01_DP
            CCS(4) =  3.346271174E-01_DP
            EXS(5) =  9.260298399E-02_DP
            CCS(5) =  5.621061301E-01_DP
            EXS(6) =  4.416183978E-02_DP
            CCS(6) =  1.712994697E-01_DP
            ! . 2p expansion.
            EXP(1) =  5.868285913E+00_DP
            CCP(1) =  7.924233646E-03_DP
            EXP(2) =  1.530329631E+00_DP
            CCP(2) =  5.144104825E-02_DP
            EXP(3) =  5.475665231E-01_DP
            CCP(3) =  1.898400060E-01_DP
            EXP(4) =  2.288932733E-01_DP
            CCP(4) =  4.049863191E-01_DP
            EXP(5) =  1.046655969E-01_DP
            CCP(5) =  4.012362861E-01_DP
            EXP(6) =  4.948220127E-02_DP
            CCP(6) =  1.051855189E-01_DP
         ! . 3s and 3p orbitals.
         CASE ( 11:18 )
            ! . 3s expansion.
            EXS(1) =  3.273031938E+00_DP
            CCS(1) = -6.775596947E-03_DP
            EXS(2) =  9.200611311E-01_DP
            CCS(2) = -5.639325779E-02_DP
            EXS(3) =  3.593349765E-01_DP
            CCS(3) = -1.587856086E-01_DP
            EXS(4) =  8.636686991E-02_DP
            CCS(4) =  5.534527651E-01_DP
            EXS(5) =  4.797373812E-02_DP
            CCS(5) =  5.015351020E-01_DP
            EXS(6) =  2.724741144E-02_DP
            CCS(6) =  7.223633674E-02_DP
            ! . 3p expansion.
            EXP(1) =  5.077973607E+00_DP
            CCP(1) = -3.329929840E-03_DP
            EXP(2) =  1.340786940E+00_DP
            CCP(2) = -1.419488340E-02_DP
            EXP(3) =  2.248434849E-01_DP
            CCP(3) =  1.639395770E-01_DP
            EXP(4) =  1.131741848E-01_DP
            CCP(4) =  4.485358256E-01_DP
            EXP(5) =  6.076408893E-02_DP
            CCP(5) =  3.908813050E-01_DP
            EXP(6) =  3.315424265E-02_DP
            CCP(6) =  7.411456232E-02_DP
         ! . 4s and 4p orbitals.
         CASE ( 19:36 )
            ! . 4s expansion.
            EXS(1) =  3.232838646E+00_DP
            CCS(1) =  1.374817488E-03_DP
            EXS(2) =  3.605788802E-01_DP
            CCS(2) = -8.666390043E-02_DP
            EXS(3) =  1.717905487E-01_DP
            CCS(3) = -3.130627309E-01_DP
            EXS(4) =  5.277666487E-02_DP
            CCS(4) =  7.812787397E-01_DP
            EXS(5) =  3.163400284E-02_DP
            CCS(5) =  4.389247988E-01_DP
            EXS(6) =  1.874093091E-02_DP
            CCS(6) =  2.487178756E-02_DP
            ! . 4p expansion.
            EXP(1) =  2.389722618E+00_DP
            CCP(1) = -1.665913575E-03_DP
            EXP(2) =  7.960947826E-01_DP
            CCP(2) = -1.657464971E-02_DP
            EXP(3) =  3.415541380E-01_DP
            CCP(3) = -5.958513378E-02_DP
            EXP(4) =  8.847434525E-02_DP
            CCP(4) =  4.053115554E-01_DP
            EXP(5) =  4.958248334E-02_DP
            CCP(5) =  5.433958189E-01_DP
            EXP(6) =  2.816929784E-02_DP
            CCP(6) =  1.204970491E-01_DP
         ! . 5s and 5p orbitals.
         CASE ( 37:54 )
            ! . 5s expansion.
            EXS(1) =  1.410128298E+00_DP
            CCS(1) =  2.695439582E-03_DP
            EXS(2) =  5.077878915E-01_DP
            CCS(2) =  1.850157487E-02_DP
            EXS(3) =  1.847926858E-01_DP
            CCS(3) = -9.588628125E-02_DP
            EXS(4) =  1.061070594E-01_DP
            CCS(4) = -5.200673560E-01_DP
            EXS(5) =  3.669584901E-02_DP
            CCS(5) =  1.087619490E+00_DP
            EXS(6) =  2.213558430E-02_DP
            CCS(6) =  3.103964343E-01_DP
            ! . 5p expansion.
            EXP(1) =  3.778623374E+00_DP
            CCP(1) =  1.163246387E-04_DP
            EXP(2) =  3.499121109E-01_DP
            CCP(2) = -2.920771322E-02_DP
            EXP(3) =  1.683175469E-01_DP
            CCP(3) = -1.381051233E-01_DP
            EXP(4) =  5.404070736E-02_DP
            CCP(4) =  5.706134877E-01_DP
            EXP(5) =  3.328911801E-02_DP
            CCP(5) =  4.768808140E-01_DP
            EXP(6) =  2.063815019E-02_DP
            CCP(6) =  6.021665516E-02_DP
         ! . Higher shells.
         CASE DEFAULT
            CALL PRINT_ERROR ( "GAUSSIAN_SETUP", "Gaussian expansions unavailable for atomic number > 54.", IATOM )
         END SELECT

         ! . Fill the basis function data structure with S-function data.
         ! . Increment the shell number.
         NSHELL = NSHELL + 1
         ! . Fill the shell arrays.
         KATOM  (NSHELL) = IQM
         KLOC   (NSHELL) = LOC + 1
         KMAX   (NSHELL) = 1
         KMIN   (NSHELL) = 1
         KNG    (NSHELL) = 6
         KSTART (NSHELL) = NGAUSS + 1
         KTYPE  (NSHELL) = 1
         ! . Increment the number of primitives.
         NGAUSS = NGAUSS + 6
         ! . Fill the coefficient and exponent arrays.
         IGAUSS = KSTART(NSHELL) - 1
         SCALE  = ZETAS * ZETAS
         CS(IGAUSS+1:IGAUSS+6) = CCS(1:6)
         EX(IGAUSS+1:IGAUSS+6) = SCALE * EXS(1:6)
         ! . Increment the LOC counter.
         LOC = LOC + 1
         ! . Unnormalize the s function primitives.
         DO I = 1,6
            EE   = EX(IGAUSS+I) + EX(IGAUSS+I)
            FACS = PI32 / ( EE * SQRT ( EE ) )
            CS(IGAUSS+I) = CS(IGAUSS+I) / SQRT ( FACS )
         END DO

         ! . Fill the basis function data structure with P-function data.
         IF ( NATORB(NI) > 1 ) THEN
            ! . Increment the shell number.
            NSHELL = NSHELL + 1
            ! . Fill the shell arrays.
            KATOM  (NSHELL) = IQM
            KLOC   (NSHELL) = LOC + 1
            KMAX   (NSHELL) = 4
            KMIN   (NSHELL) = 2
            KNG    (NSHELL) = 6
            KSTART (NSHELL) = NGAUSS + 1
            KTYPE  (NSHELL) = 2
            ! . Increment the number of primitives.
            NGAUSS = NGAUSS + 6
            ! . Fill the coefficient and exponent arrays.
            IGAUSS = KSTART(NSHELL) - 1
            SCALE  = ZETAP * ZETAP
            CP(IGAUSS+1:IGAUSS+6) = CCP(1:6)
            EX(IGAUSS+1:IGAUSS+6) = SCALE * EXP(1:6)
            ! . Increment the LOC counter.
            LOC = LOC + 3
            ! . Unnormalise the p function primitives.
            DO I = 1,6
               EE   = EX(IGAUSS+I) + EX(IGAUSS+I)
               FACS = PI32 / ( EE * SQRT ( EE ) )
               FACP = 0.5_DP * FACS / EE
               CP(IGAUSS+I) = CP(IGAUSS+I) / SQRT ( FACP )
            END DO
         END IF
         END IF
      END DO

      ! . Fill the PIKATOM array for PI atoms.
      IF ( NPIBEADS > 1 ) THEN

         ! . Fill the array for the first polymer.
         PIKATOM(1:NPISHELL,1) = KATOM(NSHELL-NPISHELL+1:NSHELL)

         ! . Fill the remaining elements using the first polymer.
         DO I = 1,NPISHELL
            DO IATOM = 1,NPIATOMS
               IF ( PIKATOM(I,1) == ATMQMI(PILIST(IATOM,1)) ) THEN
	          DO ICOPY = 2,NPIBEADS
                     PIKATOM(I,ICOPY) = ATMQMI(PILIST(IATOM,ICOPY))
		  END DO
                  EXIT
               END IF
            END DO
         END DO

      END IF

      END SUBROUTINE GET_BASIS_SET

      !-----------------------------------
      SUBROUTINE GET_BASIS_FUNCTION_LABELS
      !-----------------------------------

      ! . Local scalars.
      CHARACTER ( LEN = 2 ) :: BFNAM
      INTEGER               :: I, IATOM, ISHELL, NFUNCT

      ! . Allocate the basis function label array.
      ALLOCATE ( BFLAB(1:NBASIS) )

      ! . Loop over the shells.
      NFUNCT = 0
      DO ISHELL = 1,NSHELL
         ! . Get the atom for the shell.
         IATOM = QMATOM(KATOM(ISHELL))%ATOM
         ! . Loop over the shell types.
         DO I = KMIN(ISHELL),KMAX(ISHELL)
            NFUNCT = NFUNCT + 1
            ! . Find the function/shell name.
            SELECT CASE ( I )
            CASE (     1 ) ; BFNAM = 'S '
            CASE (     2 ) ; BFNAM = 'X '
            CASE (     3 ) ; BFNAM = 'Y '
            CASE (     4 ) ; BFNAM = 'Z '
            CASE (     5 ) ; BFNAM = 'XX'
            CASE (     6 ) ; BFNAM = 'YY'
            CASE (     7 ) ; BFNAM = 'ZZ'
            CASE (     8 ) ; BFNAM = 'XY'
            CASE (     9 ) ; BFNAM = 'XZ'
            CASE (    10 ) ; BFNAM = 'YZ'
            CASE ( 11:20 ) ; BFNAM = 'F '
            CASE ( 21:35 ) ; BFNAM = 'G '
            END SELECT
            ! . Initialise the label.
            BFLAB(NFUNCT) = "            "
            ! . Encode the number of the atom for the label.
            CALL ENCODE_INTEGER ( IATOM, BFLAB(NFUNCT)(1:6) )
            ! . Build the remainder of the label.
            BFLAB(NFUNCT)(8:9)   = SYMBOL ( ATMNUM(IATOM) )
            BFLAB(NFUNCT)(11:12) = BFNAM
         END DO
      END DO

      ! . Check the number of basis functions.
      IF ( NFUNCT /= NBASIS ) CALL PRINT_ERROR ( "GAUSSIAN_SETUP", "Basis function number mismatch." )

      END SUBROUTINE GET_BASIS_FUNCTION_LABELS

      !-----------------------------------
      SUBROUTINE NORMALIZE_BASIS_FUNCTIONS
      !-----------------------------------

      ! . Local scalars.
      INTEGER            :: I, IG, II, I1, I2, JG, LIT, LOCI, MAXI, MINI
      REAL ( KIND = DP ) :: DNRM, DNRMI, DUM, DUMD, DUMF, DUMG, DUMP, DUMS, EE, &
                            FNRM, FNRMI, GNRM, GNRMI, PNRM, PNRMI, SNRM, SNRMI

      ! . Loop over the shells.
      DO II = 1,NSHELL
         I1 = KSTART(II)
         I2 = I1 + KNG(II) - 1
         LIT  = KTYPE(II)
         MINI = KMIN(II)
         MAXI = KMAX(II)
         LOCI = KLOC(II) - MINI
         ! . Initialize the normalization factors.
         SNRM = 0.0_DP
         PNRM = 0.0_DP
         DNRM = 0.0_DP
         FNRM = 0.0_DP
         GNRM = 0.0_DP
         DO I = MINI,MAXI
            SELECT CASE ( I )
            CASE (  1 ) ; SNRM = 1.0_DP
            CASE (  2 ) ; PNRM = 1.0_DP
            CASE (  5 ) ; DNRM = 1.0_DP
            CASE ( 11 ) ; FNRM = 1.0_DP
            CASE ( 21 ) ; GNRM = 1.0_DP
            END SELECT
         END DO
         ! . Calculate the normalization factors
         SNRMI = 0.0_DP
         PNRMI = 0.0_DP
         DNRMI = 0.0_DP
         FNRMI = 0.0_DP
         GNRMI = 0.0_DP
         DO IG = I1,I2
            DO JG = I1,IG
               EE   = EX(IG) + EX(JG)
               DUM  = EE * SQRT (EE)
               DUMS =             CS(IG) * CS(JG) / DUM
               DUMP = 0.5_DP    * CP(IG) * CP(JG) / ( EE      * DUM )
               DUMD = 0.75_DP   * CD(IG) * CD(JG) / ( EE ** 2 * DUM )
               DUMF = 1.875_DP  * CF(IG) * CF(JG) / ( EE ** 3 * DUM )
               DUMG = 6.5625_DP * CG(IG) * CG(JG) / ( EE ** 4 * DUM )
               IF ( JG /= IG ) THEN
                  DUMS = DUMS + DUMS
                  DUMP = DUMP + DUMP
                  DUMD = DUMD + DUMD
                  DUMF = DUMF + DUMF
                  DUMG = DUMG + DUMG
               END IF
               SNRMI = SNRMI + DUMS
               PNRMI = PNRMI + DUMP
               DNRMI = DNRMI + DUMD
               FNRMI = FNRMI + DUMF
               GNRMI = GNRMI + DUMG
            END DO
         END DO
         IF ( SNRMI > 1.0E-10_DP ) SNRM = 1.0_DP / SQRT ( SNRMI * PI32 )
         IF ( PNRMI > 1.0E-10_DP ) PNRM = 1.0_DP / SQRT ( PNRMI * PI32 )
         IF ( DNRMI > 1.0E-10_DP ) DNRM = 1.0_DP / SQRT ( DNRMI * PI32 )
         IF ( FNRMI > 1.0E-10_DP ) FNRM = 1.0_DP / SQRT ( FNRMI * PI32 )
         IF ( GNRMI > 1.0E-10_DP ) GNRM = 1.0_DP / SQRT ( GNRMI * PI32 )

         ! . Normalise the basis function coefficients.
         DO IG = I1,I2
            CS(IG) = CS(IG) * SNRM
            CP(IG) = CP(IG) * PNRM
            CD(IG) = CD(IG) * DNRM
            CF(IG) = CF(IG) * FNRM
            CG(IG) = CG(IG) * GNRM
         END DO
      END DO

      END SUBROUTINE NORMALIZE_BASIS_FUNCTIONS

   END SUBROUTINE GAUSSIAN_SETUP

END MODULE GAUSSIAN_BASIS

