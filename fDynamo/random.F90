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
!                              Random Number Module
!===============================================================================
!
! . Functions:
!
!   RANDOM                          Get a random number in the range [0,1].
!   RANDOM_GAUSS                    Get a Gaussianly distributed random number.
!   RANDOM_VECTOR                   Get a vector of random numbers.
!
! . Subroutines:
!
!   RANDOM_INITIALIZE               Initialize the random data structure.
!
! . Notes:
!
!   It is advisable that the random number seeds always be initialized explicitly
!   using RANDOM_INITIALIZE before the random number generator is used. The
!   routines should work without this but the random number seeds will always
!   start with the same values.
!
!   The sources for RANDOM and RANDOM_GAUSS are taken from the modules DUNI
!   and DNOR, respectively, from the package NMS:
!
!                KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
!           MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U.
!
!                COMPUTES DOUBLE PRECISION UNIFORM NUMBERS ON [0,1).
!                FROM THE BOOK, "NUMERICAL METHODS AND SOFTWARE" BY
!                D. KAHANER, C. MOLER, S. NASH, PRENTICE HALL, 1988.
!
!===============================================================================
MODULE RANDOM_NUMBERS

! . Module declarations.
USE DEFINITIONS,    ONLY : DP
USE LINEAR_ALGEBRA, ONLY : NORM, PROJECT_OUT, SCHMIDT_ORTHOGONALIZE 
USE PRINTING,       ONLY : PRINT_LINE, PRINT_PARAGRAPH

IMPLICIT NONE
PRIVATE
PUBLIC :: RANDOM, RANDOM_GAUSS, RANDOM_INITIALIZE, RANDOM_VECTOR
#ifndef PGPC
SAVE
#endif

! . Module parameters.
INTEGER,            PARAMETER :: K = 47 ! . K is the number of bits in the mantissa of the floating point word.
REAL ( KIND = DP ), PARAMETER :: CSAVE = 0.9162596898123E+13_DP   / 0.140737488355328E+15_DP, &
                                 CD    = 0.76543212345678E+14_DP  / 0.140737488355328E+15_DP, &
                                 CM    = 0.140737488355213E+15_DP / 0.140737488355328E+15_DP, &
				 AA    = 0.123758602991705622657860318807E+02_DP, &
				 B     = 0.487899177760378968003825536710E+00_DP, &
				 C     = 0.126770580788654778410032042685E+02_DP, &
				 C1    = 0.9689279_DP, &
				 C2    = 1.301198_DP,  &
				 PC    = 0.195830333955546914251231662871E-01_DP, &
				 XN    = 0.277699426966287548981739308903E+01_DP
REAL ( KIND = DP ), DIMENSION(1:65), PARAMETER :: V = (/ &
                        0.340945028703977890881422029377E+00_DP, 0.457314591866933536170242641762E+00_DP, &
                        0.539779281611666939596931167481E+00_DP, 0.606242679653048904661174385559E+00_DP, &
                        0.663169062764524854743428299352E+00_DP, 0.713697459056025891183276150202E+00_DP, &
                        0.759612474933920604605610034675E+00_DP, 0.802035600355531312751497342081E+00_DP, &
                        0.841722667978955421276418428136E+00_DP, 0.879210223208313639290346470191E+00_DP, &
                        0.914894804386750587541168254518E+00_DP, 0.949079113753090250631877133376E+00_DP, &
                        0.982000481239888200218207508382E+00_DP, 0.101384923802994173461911276018E+01_DP, &
                        0.104478103674017351825847605485E+01_DP, 0.107492538202855351339149779813E+01_DP, &
                        0.110439170226812581109973656162E+01_DP, 0.113327377624394079212251428682E+01_DP, &
                        0.116165303013393172842858957666E+01_DP, 0.118960104083873798956793871425E+01_DP, &
                        0.121718147070087121538258873613E+01_DP, 0.124445158789824683238161436879E+01_DP, &
                        0.127146348057211969375402099579E+01_DP, 0.129826504188319751190626947962E+01_DP, &
                        0.132490078218086109537654808436E+01_DP, 0.135141250993337129690631764473E+01_DP, &
                        0.137783991287001181384096757263E+01_DP, 0.140422106355997540689484486002E+01_DP, &
                        0.143059286850269131403410180874E+01_DP, 0.145699147613767157824869156623E+01_DP, &
                        0.148345265660321931119703498108E+01_DP, 0.151001216431851991531882612256E+01_DP, &
                        0.153670609335952099134607533122E+01_DP, 0.156357123503769104042967185962E+01_DP, &
                        0.159064544701425352365935513885E+01_DP, 0.161796804367444698360816323707E+01_DP, &
                        0.164558021836908161542865488149E+01_DP, 0.167352550956703867146009214486E+01_DP, &
                        0.170185032506274055254533570699E+01_DP, 0.173060454131778319060975251429E+01_DP, &
                        0.175984219903830120010946138955E+01_DP, 0.178962232156657450014351836107E+01_DP, &
                        0.182000989013069176863519209140E+01_DP, 0.185107702023027589942628767312E+01_DP, &
                        0.188290439759287281399927405628E+01_DP, 0.191558305194303202395065401364E+01_DP, &
                        0.194921657491636060191700129156E+01_DP, 0.198392392890568577258506733664E+01_DP, &
                        0.201984305290623555306662745946E+01_DP, 0.205713555999009616804474181513E+01_DP, &
                        0.209599295624939161758467989658E+01_DP, 0.213664502254438986524966622832E+01_DP, &
                        0.217937134039813565892460111431E+01_DP, 0.222451750721601784110056845259E+01_DP, &
                        0.227251855485014779874266158018E+01_DP, 0.232393382009430256940425938218E+01_DP, &
                        0.237950077408282829688673722776E+01_DP, 0.244022179797994340264326423618E+01_DP, &
                        0.250751170186531701106382130475E+01_DP, 0.258346583522542956831304962942E+01_DP, &
                        0.267139159032083601869533973173E+01_DP, 0.277699426966286466722522163764E+01_DP, &
                        0.277699426966286466722522163764E+01_DP, 0.277699426966286466722522163764E+01_DP, &
                        0.277699426966286466722522163764E+01_DP /)

! . Module scalars.
INTEGER            :: IG = 17, IR = 17, JG = 5, JR = 5
REAL ( KIND = DP ) :: CR = CSAVE

! . Module arrays.
REAL ( KIND = DP ), DIMENSION(1:17) :: UG = (/ 0.471960981577884755837789724978E+00_DP, &
					       0.930323453205669578433639632431E+00_DP, &
					       0.110161790933730836587127944899E+00_DP, &
					       0.571501996273139518362638757010E-01_DP, &
					       0.402467554779738266237538503137E+00_DP, &
					       0.451181953427459489458279456915E+00_DP, &
					       0.296076152342721102174129954053E+00_DP, &
					       0.128202189325888116466879622359E-01_DP, &
					       0.314274693850973603980853259266E+00_DP, &
					       0.335521366752294932468163594171E-02_DP, &
					       0.488685045200439371607850367840E+00_DP, &
					       0.195470426865656758693860613516E+00_DP, &
					       0.864162706791773556901599326053E+00_DP, &
					       0.335505955815259203596381170316E+00_DP, &
					       0.377190200199058085469526470541E+00_DP, &
					       0.400780392114818314671676525916E+00_DP, &
					       0.374224214182207466262750307281E+00_DP /), &
                                       UR = (/ 0.471960981577884755837789724978E+00_DP, &
					       0.930323453205669578433639632431E+00_DP, &
				               0.110161790933730836587127944899E+00_DP, &
				   	       0.571501996273139518362638757010E-01_DP, &
      					       0.402467554779738266237538503137E+00_DP, &
					       0.451181953427459489458279456915E+00_DP, &
					       0.296076152342721102174129954053E+00_DP, &
					       0.128202189325888116466879622359E-01_DP, &
					       0.314274693850973603980853259266E+00_DP, &
					       0.335521366752294932468163594171E-02_DP, &
					       0.488685045200439371607850367840E+00_DP, &
					       0.195470426865656758693860613516E+00_DP, &
					       0.864162706791773556901599326053E+00_DP, &
					       0.335505955815259203596381170316E+00_DP, &
					       0.377190200199058085469526470541E+00_DP, &
					       0.400780392114818314671676525916E+00_DP, &
					       0.374224214182207466262750307281E+00_DP /)

!==============================================================================
CONTAINS
!==============================================================================

   !------------------
   FUNCTION RANDOM ( )
   !------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: RANDOM

   ! . Local scalars.
   REAL ( KIND = DP ) :: DUNI

   ! . Basic generator is Fibonacci.
   DUNI = UR(IR) - UR(JR)
   IF ( DUNI < 0.0_DP ) DUNI = DUNI + 1.0_DP
   UR(IR) = DUNI
   IR = IR - 1
   IF ( IR == 0 ) IR = 17
   JR = JR - 1
   IF ( JR == 0 ) JR = 17

   ! . Second generator is congruential.
   CR = CR - CD
   IF ( CR < 0.0_DP ) CR = CR + CM

   ! . Combination generator.
   DUNI = DUNI - CR
   IF ( DUNI < 0.0_DP ) DUNI = DUNI + 1.0_DP

   ! . Set the seed.
   RANDOM = DUNI

   END FUNCTION RANDOM

   !---------------------------------
   FUNCTION RANDOM_GAUSS ( MEAN, SD )
   !---------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ) :: RANDOM_GAUSS

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: MEAN, SD

   ! . Local scalars.
   INTEGER            :: J
   REAL ( KIND = DP ) :: DRNOR, S, UN, VNI, X, Y

   ! . Basic generator is fibonacci.
   UN = UG(IG) - UG(JG)
   IF ( UN < 0.0_DP ) UN = UN + 1.0_DP
   UG(IG) = UN

   ! . UG(IG) and un are uniform on [0,1], vni is uniform on [-1,1].
   VNI = UN + UN - 1.0_DP
   IG  = IG - 1
   IF ( IG == 0 ) IG = 17
   JG = JG - 1
   IF ( JG == 0 ) JG = 17

   ! . INT(UN(IG)*128) in range [0,127],  J is in range [1,64].
   J = MOD ( INT ( UG(IG) * 128 ), 64 ) + 1

   ! . Pick sign as VNI is positive or negative.
   DRNOR = VNI * V(J+1)
   IF ( ABS ( DRNOR ) <= V(J) ) GO TO 999

   ! . Slow part: AA is A*F(0).
   X = ( ABS ( DRNOR ) - V(J) ) / ( V(J+1) - V(J) )

   ! . Y is uniform on [0,1].
   Y = UG(IG) - UG(JG)
   IF ( Y < 0.0_DP ) Y = Y + 1.0_DP 
   UG(IG) = Y
   IG = IG - 1
   IF ( IG == 0 ) IG = 17
   JG = JG - 1
   IF ( JG == 0 ) JG = 17

   S = X + Y
   IF ( S  > C2 ) GO TO 11
   IF ( S <= C1 ) GO TO 999
   IF ( Y > C - AA * EXP ( -0.5_DP * ( B - B * X )**2 ) ) GO TO 11 
   IF ( EXP ( -0.5_DP * V(J+1)**2 ) + Y * PC / V(J+1) <= EXP ( -0.5_DP * DRNOR**2 ) ) GO TO 999

   ! . Tail part:  .36010157... is 1.0_DP / XN.
   ! . Y is uniform on [0,1].
   22 Y = UG(IG) - UG(JG)
   IF ( Y <= 0.0_DP )  Y = Y + 1.0_DP
   UG(IG) = Y
   IG = IG - 1
   IF ( IG == 0 ) IG = 17
   JG = JG - 1
   IF ( JG == 0 ) JG = 17

   X = 0.360101571301190680192994239651_DP * LOG ( Y )

   ! . Y is uniform on [0,1].
   Y = UG(IG) - UG(JG)
   IF ( Y <= 0.0_DP )  Y = Y + 1.0_DP
   UG(IG) = Y
   IG = IG - 1
   IF ( IG == 0 ) IG = 17
   JG = JG - 1
   IF ( JG == 0 ) JG = 17
   IF ( -2.0_DP * LOG ( Y ) <= X**2 ) GO TO 22
   DRNOR = SIGN ( XN - X, DRNOR )
   GO TO 999

   11 DRNOR = SIGN ( B - B * X, DRNOR )

   ! . Calculate the deviate with the appropriate mean and standard deviation.
   999 CONTINUE
   RANDOM_GAUSS = SD * DRNOR + MEAN

   END FUNCTION RANDOM_GAUSS

   !----------------------------------------------------
   SUBROUTINE RANDOM_INITIALIZE ( SEED0, GAUSSIAN_SEED )
   !----------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN)           :: SEED0
   INTEGER, INTENT(IN), OPTIONAL :: GAUSSIAN_SEED

   ! . Local scalars.
   INTEGER            :: I1, II, J1, JJ, K1, L1, M1, SEEDG
   REAL ( KIND = DP ) :: S, T

   ! . Set up the calculation for uniformly-distributed numbers.
   IR = 17
   JR =  5
   CR = CSAVE

   IF ( SEED0 /= 0 ) THEN

      ! . Calculate four smallish positive integers.
      I1 = MOD(ABS(SEED0),177)+1
      J1 = MOD(ABS(SEED0),167)+1
      K1 = MOD(ABS(SEED0),157)+1
      L1 = MOD(ABS(SEED0),147)+1

      ! . Generate a random bit pattern based upon the input seed.
      DO II = 1,17
         S = 0.0_DP 
         T = 0.5_DP
         DO JJ = 1,K
            M1 = MOD(MOD(I1*J1,179)*K1,179) 
            I1 = J1
            J1 = K1
            K1 = M1
            L1 = MOD(53*L1+1,169) 
            IF ( MOD(L1*M1,64) >= 32 ) S = S + T
            T = 0.5_DP * T
         END DO
         UR(II) = S
      END DO
   END IF

   ! . Get the initial Gaussian seed.
   IF ( PRESENT ( GAUSSIAN_SEED ) ) THEN
      SEEDG = GAUSSIAN_SEED
   ELSE
      SEEDG = SEED0
   END IF

   ! . Set up the calculation for normally-distributed numbers.
   IG = 17
   JG =  5

   IF ( SEEDG /= 0 ) THEN

      ! . Generate a random bit pattern in array based upon the input seed.
      I1 = MOD ( ABS ( SEEDG ), 32707 )
      J1 = 1111
      K1 = 1947

      DO II = 1,17
         S = 0.0_DP 
         T = 0.5_DP 
         DO JJ = 1,95
            L1 = K1-I1
            IF ( L1 < 0 ) THEN
               L1 = L1 + 32707
               S  = S + T
	    END IF
            I1 = J1
            J1 = K1
            K1 = L1
            T = 0.5_DP * T
	 END DO
         UG(II) = S
      END DO
   END IF

   ! . Write out some information.
   WRITE ( PRINT_LINE, "(A,I10,A,I10,A)" ) "Random number seeds initialized to ", ABS ( SEED0 ), " and ", ABS ( SEEDG ), "."
   CALL PRINT_PARAGRAPH

   END SUBROUTINE RANDOM_INITIALIZE

   !-----------------------------------------------------
   FUNCTION RANDOM_VECTOR ( N, NORMALIZE, ORTHOGONAL_TO )
   !-----------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN)           :: N
   LOGICAL, INTENT(IN), OPTIONAL :: NORMALIZE

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN), OPTIONAL :: ORTHOGONAL_TO

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:N) :: RANDOM_VECTOR

   ! . Local parameters.
   REAL ( KIND = DP ), PARAMETER :: TOL = 1.0E-10_DP

   ! . Local scalars.
   INTEGER            :: I, NVEC
   LOGICAL            :: QNORMALIZE, QPROJECT
   REAL ( KIND = DP ) :: FACT

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: VECTORS

   ! . Get the normalization flag.
   IF ( PRESENT ( NORMALIZE ) ) THEN
      QNORMALIZE = NORMALIZE
   ELSE
      QNORMALIZE = .FALSE.
   END IF

   ! . Initialize the projection flag.
   QPROJECT = .FALSE.

   ! . See if a set of vectors has been input.
   IF ( PRESENT ( ORTHOGONAL_TO ) ) THEN

      ! . Check the vector dimensions.
      IF ( ( SIZE ( ORTHOGONAL_TO, 1 ) == N ) .AND. ( SIZE ( ORTHOGONAL_TO, 2 ) > 0 ) ) THEN

         ! . Store the vectors.
         ALLOCATE ( VECTORS(1:N,1:SIZE ( ORTHOGONAL_TO, 2 )) ) ; VECTORS = ORTHOGONAL_TO

         ! . Schmidt orthogonalize the vectors.
         CALL SCHMIDT_ORTHOGONALIZE ( VECTORS, NVEC )

         ! . NVEC already spans the full space.
         IF ( NVEC >= N ) THEN
            RANDOM_VECTOR = 0.0_DP ; RETURN
         ! . There are some independent vectors.
         ELSE IF ( NVEC > 0 ) THEN
            QPROJECT = .TRUE.
         END IF

      END IF
   END IF

   ! . Top of the loop over vector generation.
   10 CONTINUE

   ! . Return the random numbers.
   DO I = 1,N
      RANDOM_VECTOR(I) = RANDOM ( )
   END DO

   ! . Project out any vectors.
   IF ( QPROJECT ) CALL PROJECT_OUT ( RANDOM_VECTOR, VECTORS(1:N,1:NVEC) )

   ! . Check the size of the vector.
   IF ( QNORMALIZE .OR. QPROJECT ) THEN

      ! . Get the norm of the vector.
      FACT = NORM ( RANDOM_VECTOR )

      ! . Check the size of the vector.
      IF ( FACT == 0.0_DP ) GO TO 10

   END IF

   ! . Normalize the vector.
   IF ( QNORMALIZE ) RANDOM_VECTOR = RANDOM_VECTOR / FACT

   ! . Release any temporary arrays.
   IF ( QPROJECT ) DEALLOCATE ( VECTORS )

   END FUNCTION RANDOM_VECTOR

END MODULE RANDOM_NUMBERS
