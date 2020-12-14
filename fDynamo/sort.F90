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
!                               The Sorting Module
!===============================================================================
!
! . Subroutines:
!
!   SORT_INTEGER                      Order one or two integer arrays.
!   SORT_REAL_RANK                    Rank a real array (ascending order).
!
! . Notes:
!
!   SORT_INTEGER makes use of a modification of the public domain subroutine
!   ISORT written by R. E. Jones, D. K. Kahaner and J. A. Wisniewski. See:
!
!   Singleton,R.C., Algorithm 347, "An Efficient Algorithm For Sorting With
!   Minimal Storage", CACM, 12(3), 1969, 185-7.
!
!   SORT_REAL_RANK is a modification of the public domain subroutine SORT2
!   from NAPACK (obtained from NETLIB).
!
!===============================================================================
MODULE SORT

USE DEFINITIONS, ONLY : DP

IMPLICIT NONE
PRIVATE
PUBLIC :: SORT_INTEGER, SORT_REAL_RANK

!===============================================================================
CONTAINS
!===============================================================================

   !------------------------------------------
   SUBROUTINE SORT_INTEGER ( X, Y, ASCENDING )
   !------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: ASCENDING

   ! . Array arguments.
   INTEGER, DIMENSION(:),         INTENT(INOUT)           :: X
   INTEGER, DIMENSION(1:SIZE(X)), INTENT(INOUT), OPTIONAL :: Y

   ! . Local parameters.
   INTEGER, PARAMETER :: NDIM = 21

   ! . Local scalars.
   INTEGER            :: I, IJ, J, K, L, M, NN, T, TT, TTY, TY
   LOGICAL            :: QUP
   REAL ( KIND = DP ) :: R

   ! . Local arrays.
   INTEGER, DIMENSION(1:NDIM) :: IL, IU

   ! . Get the size of the array.
   NN = SIZE ( X )
   IF ( NN <= 1 ) RETURN

   ! . Check NN versus NDIM (do nothing if the array is too big).
   IF ( NN > ( 2**(NDIM+1) - 1 ) ) RETURN

   ! . Check to see if the array is to be ordered in ascending order.
   QUP = .TRUE.
   IF ( PRESENT ( ASCENDING ) ) QUP = ASCENDING

   ! . Alter X to get descending order.
   IF ( .NOT. QUP ) X = -X

   ! . Sort X only.
   IF ( .NOT. PRESENT ( Y ) ) THEN

      M = 1
      I = 1
      J = NN
      R = 0.375_DP

   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437_DP) THEN
         R = R+3.90625E-2_DP
      ELSE
         R = R-0.21875_DP
      ENDIF

   30 K = I

      ! . Select a central element of the array and save it in location T

      IJ = I + INT((J-I)*R)
      T = X(IJ)

      ! . If first element of array is greater than T, interchange with T

      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J

      ! . If last element of array is less than than T, interchange with T

      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)

         ! . If first element of array is greater than T, interchange with T

         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF

      ! . Find an element in the second half of the array which is smaller than T

   40 L = L-1
      IF (X(L) .GT. T) GO TO 40

      !  . Find an element in the first half of the array which is greater than T

   50 K = K+1
      IF (X(K) .LT. T) GO TO 50

      ! . Interchange these elements

      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF

      ! . Save upper and lower subscripts of the array yet to be sorted

      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70

      ! . Begin again on another portion of the unsorted array

   60 M = M-1
      IF (M .EQ. 0) GO TO 300
      I = IL(M)
      J = IU(M)

   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1

   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I

   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80

   ! . Sort X and Y.
   ELSE

      M = 1
      I = 1
      J = NN
      R = 0.375_DP

  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437_DP) THEN
         R = R+3.90625E-2_DP
      ELSE
         R = R-0.21875_DP
      ENDIF

  120 K = I

      ! . Select a central element of the array and save it in location T

      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)

      ! . If first element of array is greater than T, interchange with T

      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
      ENDIF
      L = J

      ! . If last element of array is less than T, interchange with T

      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)

      ! .    If first element of array is greater than T, interchange with T

         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
         ENDIF
      ENDIF

      ! . Find an element in the second half of the array which is smaller than T

  130 L = L-1
      IF (X(L) .GT. T) GO TO 130

      ! . Find an element in the first half of the array which is greater than T

  140 K = K+1
      IF (X(K) .LT. T) GO TO 140

      ! . Interchange these elements

      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
         GO TO 130
      ENDIF

      ! . Save upper and lower subscripts of the array yet to be sorted

      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160

      ! . Begin again on another portion of the unsorted array

  150 M = M-1
      IF (M .EQ. 0) GO TO 300
      I = IL(M)
      J = IU(M)

  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1

  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I

  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      GO TO 170

   END IF

   ! . End of the subroutine.
  300 CONTINUE

   ! . Reset X if necessary.
   IF ( .NOT. QUP ) X = -X

   END SUBROUTINE SORT_INTEGER

   !---------------------------------
   SUBROUTINE SORT_REAL_RANK ( X, Y )
   !---------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:),         INTENT(IN)  :: X
   INTEGER,            DIMENSION(1:SIZE(X)), INTENT(OUT) :: Y

   ! . Local scalars.
   INTEGER            :: I, J, K, L, M, N, P, Q
   REAL ( KIND = DP ) :: S, T

   ! . Local arrays.
   INTEGER, DIMENSION(1:SIZE(X)) :: W

   ! . Set N.
   N = SIZE ( X )

   ! . Do the sorting.
      I = 1
10    K = I
20    J = I
      Y(I) = I
      I = I + 1
      IF ( J .EQ. N ) GOTO 30
      IF ( X(I) .GE. X(J) ) GOTO 20
      W(K) = I
      GOTO 10
30    IF ( K .EQ. 1 ) RETURN
      W(K) = N + 1
40    M = 1
      L = 1
50    I = L
      IF ( I .GT. N ) GOTO 120
      P = Y(I)
      S = X(P)
      J = W(I)
      K = J
      IF ( J .GT. N ) GOTO 100
      Q = Y(J)
      T = X(Q)
      L = W(J)
      Y(I) = L
60    IF ( S .GT. T ) GOTO 70
      W(M) = P
      M = M + 1
      I = I + 1
      IF ( I .EQ. K ) GOTO 80
      P = Y(I)
      S = X(P)
      GOTO 60
70    W(M)= Q
      M = M + 1
      J = J + 1
      IF ( J .EQ. L ) GOTO 110
      Q = Y(J)
      T = X(Q)
      GOTO 60
80    W(M) = Q
      K = M + L - J
      I = J - M
90    M = M + 1
      IF ( M .EQ. K ) GOTO 50
      W(M) = Y(M+I)
      GOTO 90
100   Y(I) = J
      L = J
110   W(M) = P
      K = M + K - I
      I = I - M
      GOTO 90
120   I = 1
130   K = I
      J = Y(I)
140   Y(I) = W(I)
      I = I + 1
      IF ( I .LT. J ) GOTO 140
      W(K) = I
      IF ( I .LE. N ) GOTO 130
      IF ( K .GT. 1 ) GOTO 40

   END SUBROUTINE SORT_REAL_RANK

END MODULE SORT
