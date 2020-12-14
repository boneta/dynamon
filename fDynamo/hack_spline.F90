MODULE SPLINE

#ifdef	HACK_SPLN

! .Module declarations.
USE CONSTANTS,        ONLY : TO_DEGREES, PI
USE DEFINITIONS,      ONLY : DP
USE FILES,            ONLY : NEXT_UNIT
USE LINEAR_ALGEBRA,   ONLY : CROSS_PRODUCT
USE ATOMS,            ONLY : ATMCRD, NATOMS, NFREE
USE PARSING,          ONLY : POP_UNIT, PUSH_UNIT, GET_LINE, GET_INTEGER, GET_REAL
USE PRINTING,         ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, &
                             PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, PRINT_SUMMARY_STOP,  &
                             PRINT_TABLE_ELEMENT, PRINT_TABLE_OPTIONS, PRINT_TABLE_START,     &
			                 PRINT_TABLE_STOP

IMPLICIT NONE
PRIVATE
PUBLIC :: SPLINE_DEFINE, SPLINE_INITIALIZE, SPLINE_POINT_DEFINE, ENERGY_SPLINE
SAVE

! . The point type definition.
TYPE SPOINT_TYPE
   INTEGER                                   :: NATOMS
   INTEGER,            DIMENSION(:), POINTER :: INDICES
END TYPE SPOINT_TYPE

! . The linked point type definition.
TYPE LINKED_SPOINT_TYPE
   TYPE(SPOINT_TYPE),        POINTER :: POINT
   TYPE(LINKED_SPOINT_TYPE), POINTER :: NEXT_POINT
END TYPE LINKED_SPOINT_TYPE

! . The spline type definition.
TYPE SPLINE_TYPE
   CHARACTER ( LEN = 32 )                    :: TYPE
   INTEGER                                   :: NPOINTS
   TYPE(SPOINT_TYPE),  DIMENSION(:), POINTER :: POINTS
   TYPE(SPLINE_TYPE),                POINTER :: NEXT_POINT
   INTEGER                                   :: NX, NY
   REAL ( KIND = DP ), DIMENSION(:), POINTER :: X, Y, Z, Z21, Z22
   REAL ( KIND = DP ), DIMENSION(:), POINTER :: WEIGHTS
END TYPE SPLINE_TYPE

! . Point data.
INTEGER :: NSPOINTS = 0

! . Type data.
TYPE(LINKED_SPOINT_TYPE), POINTER :: SPOINT_FIRST, SPOINT_LAST

! . Spline data.
INTEGER :: NSPLINE = 0

! . Type data.
TYPE(SPLINE_TYPE), POINTER :: SPLINE_FIRST, SPLINE_LAST


!===============================================================================
CONTAINS
!===============================================================================
   !------------------------------------------------------------------------
    recursive subroutine quick_sort( n, lo0, hi0, idx, var )
        implicit none
        integer                          , intent(in) :: n, lo0, hi0
        integer, dimension(1:n)       , intent(inout) :: idx
        real( kind=8 ), dimension(1:n)   , intent(in) :: var
    
        integer                        :: lo, hi, swp
        real( kind=8 )                 :: mid
    
        lo = lo0
        hi = hi0
        if( hi0 > lo0 ) then
            mid = var( idx( ( lo0 + hi0 ) / 2 ) )
            do while( lo <= hi )
                do while( lo < hi0 .and. var( idx( lo ) ) < mid )
                    lo = lo + 1
                end do
                do while( hi > lo0 .and. var( idx( hi ) ) > mid )
                    hi = hi - 1
                end do
                if( lo <= hi ) then
                    swp = idx( lo )
                    idx( lo ) = idx( hi )
                    idx( hi ) = swp
                    lo = lo + 1
                    hi = hi - 1
                end if
            end do
            if( lo0 < hi ) call quick_sort( n, lo0, hi, idx, var )
            if( lo < hi0 ) call quick_sort( n, lo, hi0, idx, var )
        end if
    end subroutine
   !------------------------------------------------------------------------
   ! - Taken from numerical recipes: f3-3.pdf
   SUBROUTINE INIT_SPLINE1D( N, X, Z, Z2 )
       IMPLICIT NONE
       INTEGER, INTENT(IN)                            :: N
       REAL( KIND=DP ), DIMENSION(1:N), INTENT(IN)    :: X, Z
       REAL( KIND=DP ), DIMENSION(1:N), INTENT(INOUT) :: Z2

       INTEGER                         :: I
       REAL( KIND=DP )                 :: P, SIG
       REAL( KIND=DP ), DIMENSION(1:N) :: U

       Z2 = .0_DP
       U  = .0_DP
       DO I = 2, N - 1
           SIG   = ( X(I) - X(I-1) )/( X(I+1) - X(I-1) )
           P     = SIG * Z2(I-1) + 2._DP
           Z2(I) = ( SIG - 1._DP ) / P
           U(I)  = ( 6._DP * ( ( Z(I+1) - Z(I) ) / ( X(I+1) - X(I) ) - &
                   ( Z(I) - Z(I-1) ) / ( X(I) - X(I-1) ) ) / ( X(I+1) - &
                     X(I-1) ) - SIG * U(I-1) ) / p
       END DO
       DO I = N - 1, 1, -1
           Z2(I) = Z2(I) * Z2(I+1) + U(I)
       END DO
   END SUBROUTINE INIT_SPLINE1D

   SUBROUTINE CALC_SPLINE1D( N, X, Z, Z2, X0, Z0, DZ )
       IMPLICIT NONE
       INTEGER, INTENT(IN)                         :: N
       REAL( KIND=DP ), DIMENSION(1:N), INTENT(IN) :: X, Z, Z2
       REAL( KIND=DP ), INTENT(IN)                 :: X0
       REAL( KIND=DP ), INTENT(OUT)                :: Z0, DZ

       INTEGER         :: K, KHI, KLO
       REAL( KIND=DP ) :: A, B, H, RX

       IF( X(1) > X(N) ) CALL PRINT_ERROR ( "CALC_SPLINE", "Independent variable is not in ascending order." )

       RX = X0

       IF( X0 < X(1) ) THEN
           RX = X(1)
!           Z0 = Z(1)
!           DZ = .0_DP
!           RETURN
            
       END IF

       IF( X0 > X(N) ) THEN
           RX = X(N)
!           Z0 = Z(N)
!           DZ = .0_DP
!           RETURN
       END IF

       KLO = 1
       KHI = N
!       RX = X0
 10    IF( KHI - KLO > 1 ) THEN
           K = ( KHI + KLO ) / 2
           IF( X(K) > RX ) THEN
               KHI = K
           ELSE
               KLO = K
           END IF
           GOTO 10
       END IF

       H  = X(KHI) - X(KLO)
       A  = ( X(KHI) - RX ) / H
       B  = ( RX - X(KLO) ) / H
       Z0 = A * Z(KLO) + B * Z(KHI) + &
            ( ( A ** 3 - A ) * Z2(KLO) + ( B ** 3 - B ) * Z2(KHI) ) * ( H ** 2 ) / 6._DP
       DZ = ( Z(KHI) - Z(KLO) ) / ( X(KHI) - X(KLO) ) + ( X(KHI) - X(KLO) ) * &
            ( ( 3._DP * B * B - 1._DP ) * Z2(KHI) - ( 3._DP *  A * A - 1._DP ) * Z2(KLO) ) / 6._DP
   END SUBROUTINE CALC_SPLINE1D


   SUBROUTINE INIT_SPLINE2D1( NX, NY, Y, Z, Z21 )
       IMPLICIT NONE
       INTEGER, INTENT(IN)                              :: NY, NX
       REAL( KIND=DP ), DIMENSION(1:NY), INTENT(IN)     :: Y
       REAL( KIND=DP ), DIMENSION(1:NY*NX), INTENT(IN)  :: Z
       REAL( KIND=DP ), DIMENSION(1:NY*NX), INTENT(OUT) :: Z21

       INTEGER                          :: J, K
       REAL( KIND=DP ), DIMENSION(1:NY) :: TMP, TMP2

       DO J = 1, NX
           DO K = 1, NY
               TMP(K) = Z( NY * ( J -1 ) + K )
           END DO
           CALL INIT_SPLINE1D( NY, Y, TMP, TMP2 )
           DO K = 1, NY
               Z21( NY * ( J - 1 ) + K ) = TMP2(K)
           END DO
       END DO
   END SUBROUTINE INIT_SPLINE2D1

   SUBROUTINE INIT_SPLINE2D2( NX, X, NY, Z, Z22 )
       IMPLICIT NONE
       INTEGER, INTENT(IN)                              :: NY, NX
       REAL( KIND=DP ), DIMENSION(1:NX), INTENT(IN)     :: X
       REAL( KIND=DP ), DIMENSION(1:NY*NX), INTENT(IN)  :: Z
       REAL( KIND=DP ), DIMENSION(1:NY*NX), INTENT(OUT) :: Z22

       INTEGER                          :: J, K
       REAL( KIND=DP ), DIMENSION(1:NX) :: TMP, TMP2

       DO K = 1, NY
           DO J = 1, NX
               TMP(J) = Z( NY * ( J -1 ) + K )
           END DO
           CALL INIT_SPLINE1D( NX, X, TMP, TMP2 )
           DO J = 1, NX
               Z22( NY * ( J - 1 ) + K ) = TMP2(J)
           END DO
       END DO
   END SUBROUTINE INIT_SPLINE2D2

   SUBROUTINE CALC_SPLINE2D1( NX, X, NY, Y, Z, Z2, X0, Y0, Z0, DZ )
       IMPLICIT NONE
       INTEGER, INTENT(IN)                             :: NX, NY
       REAL( KIND=DP ), DIMENSION(1:NX), INTENT(IN)    :: X
       REAL( KIND=DP ), DIMENSION(1:NY), INTENT(IN)    :: Y
       REAL( KIND=DP ), DIMENSION(1:NX*NY), INTENT(IN) :: Z, Z2
       REAL( KIND=DP ), INTENT(IN)                     :: X0, Y0
       REAL( KIND=DP ), INTENT(OUT)                    :: Z0, DZ

       INTEGER                          :: J, K
       REAL( KIND=DP )                  :: TMP
       REAL( KIND=DP ), DIMENSION(1:NX) :: TMP2, TMP3
       REAL( KIND=DP ), DIMENSION(1:NY) :: TMP4, TMP5

       DO J = 1, NX
           DO K = 1, NY
               TMP4(K) =  Z( NY * ( J - 1 ) + K )
               TMP5(K) = Z2( NY * ( J - 1 ) + K )
           END DO
           CALL CALC_SPLINE1D( NY, Y, TMP4, TMP5, Y0, TMP2(J), TMP )
       END DO
       CALL INIT_SPLINE1D( NX, X, TMP2, TMP3 )
       CALL CALC_SPLINE1D( NX, X, TMP2, TMP3, X0, Z0, DZ )
   END SUBROUTINE CALC_SPLINE2D1

   SUBROUTINE CALC_SPLINE2D2( NX, X, NY, Y, Z, Z2, X0, Y0, Z0, DZ )
       IMPLICIT NONE
       INTEGER, INTENT(IN)                             :: NX, NY
       REAL( KIND=DP ), DIMENSION(1:NX), INTENT(IN)    :: X
       REAL( KIND=DP ), DIMENSION(1:NY), INTENT(IN)    :: Y
       REAL( KIND=DP ), DIMENSION(1:NX*NY), INTENT(IN) :: Z, Z2
       REAL( KIND=DP ), INTENT(IN)                     :: X0, Y0
       REAL( KIND=DP ), INTENT(OUT)                    :: Z0, DZ

       INTEGER                          :: J, K
       REAL( KIND=DP )                  :: TMP
       REAL( KIND=DP ), DIMENSION(1:NX) :: TMP2, TMP3
       REAL( KIND=DP ), DIMENSION(1:NY) :: TMP4, TMP5

       DO K = 1, NY
           DO J = 1, NX
               TMP2(J) =  Z( NY * ( J - 1 ) + K )
               TMP3(J) = Z2( NY * ( J - 1 ) + K )
           END DO
           CALL CALC_SPLINE1D( NX, X, TMP2, TMP3, X0, TMP4(K), TMP )
       END DO
       CALL INIT_SPLINE1D( NY, Y, TMP4, TMP5 )
       CALL CALC_SPLINE1D( NY, Y, TMP4, TMP5, Y0, Z0, DZ )
   END SUBROUTINE CALC_SPLINE2D2
   !
   !------------------------------------------------------------------------

   !----------------------------------------------------
   SUBROUTINE ENERGY_SPLINE ( ECORR, GRADIENT, HESSIAN )
   !----------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(INOUT) :: ECORR

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS),              INTENT(INOUT), OPTIONAL :: GRADIENT
   REAL ( KIND = DP ), DIMENSION(1:(3*NFREE*(3*NFREE+1))/2), INTENT(INOUT), OPTIONAL :: HESSIAN

   ! . Initialization.
   ECORR = 0.0_DP

   ! . Check to see if any energy terms need to be calculated.
   IF ( NSPLINE <= 0 ) RETURN

   ! . Check the consistency of the derivative options.
   IF ( PRESENT ( HESSIAN ) ) CALL PRINT_ERROR ( "ENERGY_SPLINE", "Second derivatives unavailable for splines." )

   ! . Call the various energy subroutines.
   CALL ENERGY_SPLINEX ( ECORR, GRADIENT )

   END SUBROUTINE ENERGY_SPLINE


   !------------------------------------------------------------
   SUBROUTINE SPLINE_DEFINE ( TYPE, XY_FILE, XYZ_FILE, WEIGHTS )
   !------------------------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: TYPE

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: XY_FILE, XYZ_FILE
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: WEIGHTS

   ! . Local scalars.
   INTEGER                          :: I, IOSTAT, NCOEF, IOUNIT, MDF
   LOGICAL                          :: QOK
   TYPE(SPLINE_TYPE),   POINTER     :: CLOCAL
   TYPE(LINKED_SPOINT_TYPE), POINTER :: PCURRENT, PNEXT

   INTEGER, DIMENSION(:), ALLOCATABLE           :: IDX
   REAL( KIND = DP ), DIMENSION(:), ALLOCATABLE :: SWP


   ! . Check to see if there are any atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Initialization.
   NCOEF = 0
   QOK   = .TRUE.
   MDF   = 0

   ! . Branch on the type of constraint.
   SELECT CASE ( TYPE )
   CASE ( "ANGLE"             ) ; QOK = ( NSPOINTS == 3 )
   CASE ( "DIHEDRAL"          ) ; QOK = ( NSPOINTS == 4 )
   CASE ( "DISTANCE"          ) ; QOK = ( NSPOINTS == 2 )
   CASE ( "MULTIPLE_DISTANCE" ) ; NCOEF = NSPOINTS / 2
       IF( PRESENT( XY_FILE ) ) THEN
           MDF = 1
           IF( .NOT. ( NSPOINTS > 0 .AND. 2 * NCOEF == NSPOINTS ) ) &
               CALL PRINT_ERROR ( "SPLINE_DEFINE", "The spline has an invalid number of points (1D)." )
       ELSE
           NCOEF = 1
           MDF   = 2
           IF( NSPOINTS /= 4 ) &
               CALL PRINT_ERROR ( "SPLINE_DEFINE", "The spline has an invalid number of points (2D)." )
       END IF
   CASE ( "MULTIPLE_ANGLE"    ) ; QOK = ( NSPOINTS == 6 )
! ---------------------------------------------------------------------------------------------------------
! - Esto corresponde a un par de coordenadas: antisimétrica (A -- B -- C) [X] vs distancia (D -- E) [Y]
   CASE ( "HOMEBREW" )
       NCOEF = 1
       QOK = .TRUE.
       NSPOINTS = 5
       CALL PRINT_PARAGRAPH( TEXT = "[HOMEBREW]: ANTISYMMETRIC DIST. (WITH COMMON NODE) + DISTANCE: 5 ATOMS" )
! ---------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------
! - Esto corresponde a un par de coordenadas antisimétrica (A -- B -- C) [X] vs (D -- E -- F) [Y]
   CASE ( "BURRITO_LOCO" )
       NCOEF = 1
       QOK = .TRUE.
       NSPOINTS = 6
       CALL PRINT_PARAGRAPH( TEXT = "[BURRITO_LOCO]: BOTH ANTISYMMETRIC DIST. (WITH COMMON NODE): 6 ATOMS" )
! ---------------------------------------------------------------------------------------------------------
   CASE DEFAULT ; CALL PRINT_ERROR ( "SPLINE_DEFINE", "Unknown constraint type." )
   END SELECT

   ! . Check that the number of constraint points is valid.
   IF ( .NOT. QOK ) CALL PRINT_ERROR ( "SPLINE_DEFINE", "The constraint has an invalid number of points." )

   ! . Allocate the new constraint.
   ALLOCATE ( CLOCAL )

   ! . Fill the local type variables.
   CLOCAL%TYPE    = TYPE
   IF( MDF > 0 ) THEN
       IF( MDF == 1 ) THEN
           CLOCAL%TYPE = "MULTIPLE_DISTANCE_1D"
       ELSE
           CLOCAL%TYPE = "MULTIPLE_DISTANCE_2D"
       END IF
   END IF
   CLOCAL%NPOINTS = NSPOINTS

   ALLOCATE ( CLOCAL%WEIGHTS(1:NCOEF) )
   IF( PRESENT( WEIGHTS ) ) THEN
       IF( SIZE( WEIGHTS ) /= NCOEF ) THEN
           CALL PRINT_ERROR ( "SPLINE_DEFINE", "The WEIGHTS array has the wrong size." )
       ELSE
           CLOCAL%WEIGHTS = WEIGHTS
       END IF
   ELSE
       CLOCAL%WEIGHTS = 1._DP
   END IF

   ! . Allocate the points in the constraint.
   ALLOCATE ( CLOCAL%POINTS(1:NSPOINTS) )

   ! . Loop over the points.
   DO I = 1, NSPOINTS

      ! . Get the next constraint in the list.
      IF ( I == 1 ) THEN
         PCURRENT => SPOINT_FIRST
         PNEXT    => SPOINT_FIRST%NEXT_POINT
      ELSE
         PCURRENT => PNEXT
         PNEXT    => PNEXT%NEXT_POINT
      END IF

      ! . Set the number of atoms in the point.
      CLOCAL%POINTS(I)%NATOMS = PCURRENT%POINT%NATOMS

      ! . Set the point arrays.
      CLOCAL%POINTS(I)%INDICES => PCURRENT%POINT%INDICES

      ! . Nullify the spline point arrays
      NULLIFY ( PCURRENT%POINT%INDICES )

      ! . Deallocate the point.
      DEALLOCATE ( PCURRENT%POINT )

      ! . Nullify the pointer in PCURRENT.
      NULLIFY ( PCURRENT%NEXT_POINT )

      ! . Deallocate the constraint.
      DEALLOCATE ( PCURRENT )

  END DO

   ! . Reset the number of constraint points.
   NSPOINTS = 0

   ! . Nullify the first and last constraint point pointers.
   NULLIFY ( SPLINE_FIRST, SPLINE_LAST )

   ! . Nullify the pointer in CLOCAL.
   NULLIFY ( CLOCAL%NEXT_POINT )

   ! . Increment the number of constraints.
   NSPLINE = NSPLINE + 1

   ! . This is the first constraint on the list.
   IF ( NSPLINE == 1 ) THEN

      ! . Define the first constraint pointer.
      SPLINE_FIRST => CLOCAL

   ! . There are other constraints in the list.
   ELSE

      ! . Set NEXT_CONSTRAINT of the old last constraint.
      SPLINE_LAST%NEXT_POINT => CLOCAL

   END IF

   ! . Update the last constraint pointer.
   SPLINE_LAST => CLOCAL

   ! . Read the spline table/s
   CLOCAL%NX = 0
   CLOCAL%NY = 0
   IF ( PRESENT ( XY_FILE ) ) THEN
       IOUNIT = NEXT_UNIT()
       OPEN( UNIT = IOUNIT, FILE = XY_FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )
       IF( IOSTAT /= 0 ) CALL PRINT_ERROR( "SPLINE_DEFINE", "I/O Error on correction spline file.", IOSTAT )
       CALL PUSH_UNIT( IOUNIT )
       CALL GET_LINE
       ! . Monodimensional SPLINE
       CALL GET_INTEGER( CLOCAL%NX )
       ALLOCATE( CLOCAL%X(1:CLOCAL%NX), CLOCAL%Z(1:CLOCAL%NX), CLOCAL%Z21(1:CLOCAL%NX) )
       ALLOCATE( IDX(1:CLOCAL%NX), SWP(1:CLOCAL%NX) )
       DO I = 1, CLOCAL%NX
           IDX(I) = I
           CALL GET_LINE
           CALL GET_REAL( CLOCAL%X(I) )
           CALL GET_REAL( CLOCAL%Z(I) )
       END DO
       CLOSE( IOUNIT )
       CALL POP_UNIT
       ! . Sort independent variables
       CALL QUICK_SORT( CLOCAL%NX, 1, CLOCAL%NX, IDX, CLOCAL%X )
       DO I = 1, CLOCAL%NX
           SWP(I) = CLOCAL%X(IDX(I))
       END DO
       CLOCAL%X = SWP
       DO I = 1, CLOCAL%NX
           SWP(I) = CLOCAL%Z(IDX(I))
       END DO
       CLOCAL%Z = SWP
       DEALLOCATE( IDX, SWP )
       ! . Initialize the monodimensional spline tables...
       CALL INIT_SPLINE1D( CLOCAL%NX, CLOCAL%X, CLOCAL%Z, CLOCAL%Z21 )
   END IF

   ! - No se permiten mallas no uniformes en sampleo (NX != NY, pero valores fijos)
   !      al generar los splines intermedios se encontrarian inconsistencias (pocos puntos...)
   !
   ! Haz buen uso de spline.py para pasar los datos de Excel (1:nx*ny,1:nx*ny,1:nx*ny) a dynamo
   IF ( PRESENT ( XYZ_FILE ) ) THEN
       IOUNIT = NEXT_UNIT()
       OPEN( UNIT = IOUNIT, FILE = XYZ_FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )
       IF( IOSTAT /= 0 ) CALL PRINT_ERROR( "SPLINE_DEFINE", "I/O Error on correction spline file.", IOSTAT )

       ! . Get X
       CALL PUSH_UNIT( IOUNIT )
       CALL GET_LINE
       CALL GET_INTEGER( CLOCAL%NX )
       ALLOCATE( CLOCAL%X(1:CLOCAL%NX) )
       DO I = 1, CLOCAL%NX
           CALL GET_LINE
           CALL GET_REAL( CLOCAL%X(I) )
       END DO

       ! . Get Y
       CALL GET_LINE
       CALL GET_INTEGER( CLOCAL%NY )
       ALLOCATE( CLOCAL%Y(1:CLOCAL%NY) )
       DO I = 1, CLOCAL%NY
           CALL GET_LINE
           CALL GET_REAL( CLOCAL%Y(I) )
       END DO

       ! . Get Z (Y varia para valores fijos dados de X, Ojo! que hay una linea en blanco)
       CALL GET_LINE
       CALL GET_INTEGER( I )
if(i/=clocal%nx*clocal%ny) call print_error('spline_define','No coinciden las dimensiones')
       ALLOCATE( CLOCAL%Z(1:CLOCAL%NX * CLOCAL%NY), &
                 CLOCAL%Z21(1:CLOCAL%NX * CLOCAL%NY), CLOCAL%Z22(1:CLOCAL%NX * CLOCAL%NY) )
       DO I = 1, CLOCAL%NX * CLOCAL%NY
           CALL GET_LINE
           CALL GET_REAL( CLOCAL%Z(I) )
       END DO
       CLOSE( IOUNIT )
       CALL POP_UNIT
       ! . Initialize the bidimensional spline tables...
       CALL INIT_SPLINE2D1( CLOCAL%NX, CLOCAL%NY, CLOCAL%Y, CLOCAL%Z, CLOCAL%Z21 )
       CALL INIT_SPLINE2D2( CLOCAL%NX, CLOCAL%X, CLOCAL%NY, CLOCAL%Z, CLOCAL%Z22 )
   END IF

   IF( CLOCAL%NX == 0 ) CALL PRINT_ERROR( "SPLINE_DEFINE", "No external correction spline file found!!!" )

   WRITE( PRINT_LINE, '(A,A,A,I5,A,I5)') 'Spline_Define: ', TRIM( CLOCAL%TYPE ), '  Nx:', CLOCAL%NX, '  Ny:', CLOCAL%NY
   CALL PRINT_PARAGRAPH

   END SUBROUTINE SPLINE_DEFINE

   !---------------------------
   SUBROUTINE SPLINE_INITIALIZE
   !---------------------------

   INTEGER                    :: I, P
   TYPE(SPLINE_TYPE), POINTER :: CURRENT, NEXT

   IF ( NSPLINE <= 0 ) RETURN

   DO I = 1,NSPLINE

      ! . Get the next constraint in the list.
      IF ( I == 1 ) THEN
         CURRENT => SPLINE_FIRST
         NEXT    => SPLINE_FIRST%NEXT_POINT
      ELSE
         CURRENT => NEXT
         NEXT    => NEXT%NEXT_POINT
      END IF

      DO P = 1,CURRENT%NPOINTS
          DEALLOCATE ( CURRENT%POINTS(P)%INDICES )
      END DO


      DEALLOCATE ( CURRENT%POINTS  )
      DEALLOCATE ( CURRENT%WEIGHTS )

      DEALLOCATE( CURRENT%X, CURRENT%Z, CURRENT%Z21 )
      IF( CURRENT%NY > 0 ) THEN
          DEALLOCATE( CURRENT%Y, CURRENT%Z22 )
      END IF

      NULLIFY ( CURRENT%NEXT_POINT )

      DEALLOCATE ( CURRENT )

  END DO

   NSPLINE = 0

   NULLIFY ( SPLINE_FIRST, SPLINE_LAST )

   CALL PRINT_PARAGRAPH( TEXT = "All splines have been removed." )

   END SUBROUTINE SPLINE_INITIALIZE

   !------------------------------------------------------------
   SUBROUTINE ENERGY_SPLINEX ( ECORR, GRADIENT )
   !------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(OUT)   :: ECORR

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: GRADIENT

   ! . Local scalars.
   INTEGER                    :: I, IC, N, P
   LOGICAL                    :: QGRADIENT
   TYPE(SPLINE_TYPE), POINTER :: CURRENT

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: CENTERS, DEDC

   ! . Initialize the constraint energy.
   ECORR = 0.0_DP

   ! . Check whether there are constraints.
   IF ( NSPLINE <= 0 ) RETURN

   ! . Set the gradient flag.
   QGRADIENT = PRESENT ( GRADIENT )

   ! . Loop over the constraints.
   DO IC = 1,NSPLINE

      ! . Get the next constraint on the list.
      IF ( IC == 1 ) THEN
         CURRENT => SPLINE_FIRST
      ELSE
         CURRENT => CURRENT%NEXT_POINT
      END IF

      ! . Allocate space for the points and their derivatives.
      ALLOCATE ( CENTERS(1:3,1:CURRENT%NPOINTS) )
      IF ( QGRADIENT ) ALLOCATE ( DEDC(1:3,1:CURRENT%NPOINTS) )

      ! . Calculate the positions of the points.
      CENTERS = 0.0_DP
      DO P = 1,CURRENT%NPOINTS
         DO I = 1,CURRENT%POINTS(P)%NATOMS
            N              = CURRENT%POINTS(P)%INDICES(I)
            CENTERS(1:3,P) = CENTERS(1:3,P) + ATMCRD(1:3,N)
         END DO
      END DO

      ! . Branch on the type of constraint.
      SELECT CASE ( CURRENT%TYPE )
      CASE ( "ANGLE"                ) ; CALL ENERGY_ANGLE
      CASE ( "DIHEDRAL"             ) ; CALL ENERGY_DIHEDRAL
      CASE ( "DISTANCE"             ) ; CALL ENERGY_DISTANCE
      CASE ( "MULTIPLE_DISTANCE_1D" ) ; CALL ENERGY_MULTIPLE1D
      CASE ( "MULTIPLE_DISTANCE_2D" ) ; CALL ENERGY_MULTIPLE2D
      CASE ( "MULTIPLE_ANGLE"       ) ; CALL ENERGY_MULTIPLE_ANGLE
! --------------------------------------------------------------------------------------
      CASE ( "HOMEBREW" ) ; CALL HOMEBREW
      CASE ( "BURRITO_LOCO" ) ; CALL BURRITO_LOCO
! --------------------------------------------------------------------------------------
      END SELECT

      ! . Deallocate the center array.
      DEALLOCATE ( CENTERS )

      ! . Convert the gradients with respect to points to those with respect to atoms.
      IF ( QGRADIENT ) THEN

         ! . Loop over the points.
	 DO P = 1,CURRENT%NPOINTS
            DO I = 1,CURRENT%POINTS(P)%NATOMS
               N               = CURRENT%POINTS(P)%INDICES(I)
               GRADIENT(1:3,N) = GRADIENT(1:3,N) + DEDC(1:3,P)
	    END DO
	 END DO

         ! . Deallocate the derivative array.
	 DEALLOCATE ( DEDC )

      END IF
   END DO

   !============================================================================
   CONTAINS
   !============================================================================


      !----------------------
      SUBROUTINE ENERGY_ANGLE
      !----------------------

      ! . Parameter declarations.
      REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP - 1.0E-6_DP

      ! . Local scalars.
      REAL ( KIND = DP ) :: DF, DOTFAC, RIJ, RKJ, E0, VALUE

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ, DRKJ, DTI, DTJ, DTK

      ! . Get the interatomic vectors.
      DRIJ = CENTERS(1:3,1) - CENTERS(1:3,2)
      DRKJ = CENTERS(1:3,3) - CENTERS(1:3,2)

      ! . Calculate their sizes.
      RIJ = SQRT ( DOT_PRODUCT ( DRIJ, DRIJ ) )
      RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

      ! . Normalize the vectors.
      DRIJ = DRIJ / RIJ
      DRKJ = DRKJ / RKJ

      ! . Calculate the dot product between the two vectors.
      DOTFAC = DOT_PRODUCT ( DRIJ, DRKJ )

      ! . Ensure DOTFAC is between -1 and 1.
      DOTFAC = SIGN ( MIN ( ABS ( DOTFAC ), DOT_LIMIT ), DOTFAC )

      ! . Get the angle.
      VALUE = ACOS ( DOTFAC ) * TO_DEGREES

      CALL CALC_SPLINE1D( CURRENT%NX, CURRENT%X, CURRENT%Z, CURRENT%Z21, VALUE, E0, DF )

      ! . Calculate the contribution to the energy.
      ECORR = ECORR + E0

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) RETURN

      ! . Calculate DF for the first derivative calculation.
      DF = - TO_DEGREES * DF / SQRT ( 1.0_DP - DOTFAC * DOTFAC )

      ! . Calculate the components of the derivatives.
      DTI  = ( DRKJ - DOTFAC * DRIJ ) / RIJ
      DTK  = ( DRIJ - DOTFAC * DRKJ ) / RKJ
      DTJ  = - ( DTI + DTK )

      ! . Calculate the gradients.
      DEDC(1:3,1) = DF * DTI
      DEDC(1:3,2) = DF * DTJ
      DEDC(1:3,3) = DF * DTK

      END SUBROUTINE ENERGY_ANGLE

      !-------------------------
      SUBROUTINE ENERGY_DIHEDRAL
      !-------------------------

      ! . Parameter declarations.
      REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP !- 1.0E-6_DP

      ! . Local scalars.
      REAL ( KIND = DP ) :: DF, DOTFAC, FACTIJ, FACTKL, MSIZE, NSIZE, RKJ, SGNFAC, E0, VALUE

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DRIJ, DRKJ, DRKL, DTI, DTJ, DTK, DTL, M, N

      ! . Get the interatomic vectors.
      DRIJ = CENTERS(1:3,1) - CENTERS(1:3,2)
      DRKJ = CENTERS(1:3,3) - CENTERS(1:3,2)
      DRKL = CENTERS(1:3,3) - CENTERS(1:3,4)

      ! . Get M and N.
      M = CROSS_PRODUCT ( DRIJ, DRKJ )
      N = CROSS_PRODUCT ( DRKJ, DRKL )

      ! . Get the size of M and N.
      MSIZE = SQRT ( DOT_PRODUCT ( M, M ) )
      NSIZE = SQRT ( DOT_PRODUCT ( N, N ) )

      ! . Normalize M and N.
      M = M / MSIZE
      N = N / NSIZE

      ! . Ensure DOTFAC is between -1 and 1.
      DOTFAC = DOT_PRODUCT ( M, N )
      DOTFAC = SIGN ( MIN ( ABS ( DOTFAC ), DOT_LIMIT ), DOTFAC )

      ! . Get the sign of the angle.
      SGNFAC = 1.0_DP
      IF ( DOT_PRODUCT ( DRIJ, N ) < 0.0_DP ) SGNFAC = -1.0_DP

      ! . Determine the dihedral.
      VALUE = SGNFAC * TO_DEGREES * ACOS ( DOTFAC )

      CALL CALC_SPLINE1D( CURRENT%NX, CURRENT%X, CURRENT%Z, CURRENT%Z21, VALUE, E0, DF )

      ! . Calculate the contribution to the energy.
      ECORR = ECORR + E0

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) RETURN

      ! . Scale DF.
      DF = TO_DEGREES * DF

      ! . Calculate RKJ.
      RKJ = SQRT ( DOT_PRODUCT ( DRKJ, DRKJ ) )

      ! . Calculate DTI and DTL.
      DTI =   DF * RKJ * M / MSIZE
      DTL = - DF * RKJ * N / NSIZE

      ! . Calculate some additional factors.
      FACTIJ = DOT_PRODUCT ( DRIJ, DRKJ ) / ( RKJ * RKJ )
      FACTKL = DOT_PRODUCT ( DRKL, DRKJ ) / ( RKJ * RKJ )

      ! . Calculate DTJ and DTK.
      DTJ = DTI * ( FACTIJ - 1.0_DP ) - FACTKL * DTL
      DTK = DTL * ( FACTKL - 1.0_DP ) - FACTIJ * DTI

      ! . Calculate the gradients.
      DEDC(1:3,1) = DTI
      DEDC(1:3,2) = DTJ
      DEDC(1:3,3) = DTK
      DEDC(1:3,4) = DTL

      END SUBROUTINE ENERGY_DIHEDRAL

      !-------------------------
      SUBROUTINE ENERGY_DISTANCE
      !-------------------------

      ! . Local scalars.
      REAL ( KIND = DP ) :: DF, E0, VALUE

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3) :: DR

      ! . Calculate the distance between the centers.
      DR    = CENTERS(1:3,1) - CENTERS(1:3,2)
      VALUE = SQRT ( DOT_PRODUCT ( DR, DR ) )

      CALL CALC_SPLINE1D( CURRENT%NX, CURRENT%X, CURRENT%Z, CURRENT%Z21, VALUE, E0, DF )

      ! . Calculate the contribution to the energy.
      ECORR = ECORR + E0

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) RETURN

      ! . Calculate the gradients.
      DF = DF / VALUE
      DEDC(1:3,1) =   DF * DR
      DEDC(1:3,2) = - DF * DR

      END SUBROUTINE ENERGY_DISTANCE

      !---------------------------
      SUBROUTINE ENERGY_MULTIPLE1D
      !---------------------------

      ! . Local scalars.
      INTEGER            :: I, II, NDIST
      REAL ( KIND = DP ) :: DF, E0, VALUE

      ! . Local arrays.
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: R
      REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: DR

      ! . Get the number of distances involved in the constraint.
      NDIST = SIZE ( CURRENT%WEIGHTS )

      ! . Allocate space for the distances and intercenter vectors.
      ALLOCATE ( DR(1:3,NDIST), R(1:NDIST) )

      ! . Calculate the distances between the centers.
      DO I = 1,NDIST
         II        = 2 * ( I - 1 )
         DR(1:3,I) = CENTERS(1:3,II+1) - CENTERS(1:3,II+2)
         R(I)      = SQRT ( DOT_PRODUCT ( DR(1:3,I), DR(1:3,I) ) )
      END DO

      ! . Calculate the value of the constraint.
      VALUE = DOT_PRODUCT ( CURRENT%WEIGHTS, R )


      CALL CALC_SPLINE1D( CURRENT%NX, CURRENT%X, CURRENT%Z, CURRENT%Z21, VALUE, E0, DF )

      ! . Calculate the contribution to the energy.
      ECORR = ECORR + E0

!write(*,'(/3f20.10)' ) (dr(1:3,i),i=1,ndist)
!write(*,'(/f20.10)'  ) (r(i),i=1,ndist)
!write(*,'(/f20.10)'  ) value
!write(*,'(/2f20.10)' ) e0, df

      ! . Check for a gradient calculation.
      IF ( .NOT. QGRADIENT ) RETURN

      DO I = 1,NDIST
            II             = 2 * ( I - 1 )
            DEDC(1:3,II+1) =   CURRENT%WEIGHTS(I) * DF * DR(1:3,I) / R(I)
            DEDC(1:3,II+2) = - CURRENT%WEIGHTS(I) * DF * DR(1:3,I) / R(I)
      END DO

      DEALLOCATE ( DR, R )

      END SUBROUTINE ENERGY_MULTIPLE1D

      !---------------------------
      SUBROUTINE ENERGY_MULTIPLE2D
      !---------------------------

      ! . Local scalars.
      INTEGER            :: I,  II
      REAL ( KIND = DP ), DIMENSION(1:2) :: DF, E0, VALUE

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:2,1:3) :: DR

      ! . Only Bidimensional stuff allowed
      DR(1,1:3) = CENTERS(1:3,1) - CENTERS(1:3,2)
      VALUE(1)  = SQRT ( DOT_PRODUCT ( DR(1,1:3), DR(1,1:3) ) )
      DR(2,1:3) = CENTERS(1:3,3) - CENTERS(1:3,4)
      VALUE(2)  = SQRT ( DOT_PRODUCT ( DR(2,1:3), DR(2,1:3) ) )

      DF = .0_DP
      CALL CALC_SPLINE2D1( CURRENT%NX, CURRENT%X, CURRENT%NY, CURRENT%Y, &
                           CURRENT%Z, CURRENT%Z21, VALUE(1), VALUE(2), E0(1), DF(1) )
      CALL CALC_SPLINE2D2( CURRENT%NX, CURRENT%X, CURRENT%NY, CURRENT%Y, &
                           CURRENT%Z, CURRENT%Z22, VALUE(1), VALUE(2), E0(2), DF(2) )
      ECORR = ECORR + E0(1) * .5_DP + E0(2) * .5_DP

      ! . If one of the gradient vectors is zero (out of range), ignore the other ...
      IF( DF(1) == .0_DP ) THEN
          DF(2) = .0_DP
      ELSE 
          IF( DF(2) == .0_DP ) DF(1) = .0_DP
      END IF

write(*,'(a,6f12.6)') 'S2D', value(1), value(2), e0(1), e0(2), df(1), df(2)

      IF( .NOT. QGRADIENT ) RETURN
      DF(1)       =   DF(1) / VALUE(1)
      DEDC(1:3,1) =   DF(1) * DR(1,1:3)
      DEDC(1:3,2) = - DF(1) * DR(1,1:3)
      DF(2)       =   DF(2) / VALUE(2)
      DEDC(1:3,3) =   DF(2) * DR(2,1:3)
      DEDC(1:3,4) = - DF(2) * DR(2,1:3)

      END SUBROUTINE ENERGY_MULTIPLE2D

! --------------------------------------------------------------------------------------
      SUBROUTINE HOMEBREW

      ! . Local scalars.
      INTEGER            :: I,  II
      REAL ( KIND = DP ), DIMENSION(1:2) :: DF, E0, VALUE, DR2

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:3,1:3) :: DR

      DR(1,1:3) = CENTERS(1:3,1) - CENTERS(1:3,2)
      DR2(1)    = SQRT ( DOT_PRODUCT ( DR(1,1:3), DR(1,1:3) ) )

      DR(2,1:3) = CENTERS(1:3,2) - CENTERS(1:3,3)
      DR2(2)    = SQRT ( DOT_PRODUCT ( DR(2,1:3), DR(2,1:3) ) )

      VALUE(1)  = DR2(1) - DR2(2)

      DR(3,1:3) = CENTERS(1:3,4) - CENTERS(1:3,5)
      VALUE(2)  = SQRT ( DOT_PRODUCT ( DR(3,1:3), DR(3,1:3) ) )

      DF = .0_DP
      CALL CALC_SPLINE2D1( CURRENT%NX, CURRENT%X, CURRENT%NY, CURRENT%Y, &
                           CURRENT%Z, CURRENT%Z21, VALUE(1), VALUE(2), E0(1), DF(1) )
      CALL CALC_SPLINE2D2( CURRENT%NX, CURRENT%X, CURRENT%NY, CURRENT%Y, &
                           CURRENT%Z, CURRENT%Z22, VALUE(1), VALUE(2), E0(2), DF(2) )
      ECORR = ECORR + E0(1) * .5_DP + E0(2) * .5_DP

      ! . If one of the gradient vectors is zero (out of range), ignore the other ...
      IF( DF(1) == .0_DP ) THEN
          DF(2) = .0_DP
      ELSE 
          IF( DF(2) == .0_DP ) DF(1) = .0_DP
      END IF

write(*,'(a,6f12.6)') 'S2D', value(1), value(2), e0(1), e0(2), df(1), df(2)

      IF( .NOT. QGRADIENT ) RETURN
      DEDC(1:3,1) =   DF(1) / DR2(1) * DR(1,1:3)
      DEDC(1:3,3) =   DF(1) / DR2(2) * DR(2,1:3)
      DEDC(1:3,2) = - DF(1) / DR2(1) * DR(1,1:3) - DF(1) / DR2(2) * DR(2,1:3)

      DF(2)       =   DF(2) / VALUE(2)
      DEDC(1:3,4) =   DF(2) * DR(3,1:3)
      DEDC(1:3,5) = - DF(2) * DR(3,1:3)

      END SUBROUTINE HOMEBREW
! --------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------
      SUBROUTINE BURRITO_LOCO

      ! . Local scalars.
      INTEGER            :: I,  II
      REAL ( KIND = DP ), DIMENSION(1:2) :: DF, E0, VALUE
      REAL ( KIND = DP ), DIMENSION(1:4) :: DR2

      ! . Local arrays.
      REAL ( KIND = DP ), DIMENSION(1:4,1:3) :: DR

      DR(1,1:3) = CENTERS(1:3,1) - CENTERS(1:3,2)
      DR2(1)    = SQRT ( DOT_PRODUCT ( DR(1,1:3), DR(1,1:3) ) )

      DR(2,1:3) = CENTERS(1:3,2) - CENTERS(1:3,3)
      DR2(2)    = SQRT ( DOT_PRODUCT ( DR(2,1:3), DR(2,1:3) ) )

      VALUE(1)  = DR2(1) - DR2(2)

      DR(3,1:3) = CENTERS(1:3,4) - CENTERS(1:3,5)
      DR2(3)    = SQRT ( DOT_PRODUCT ( DR(3,1:3), DR(3,1:3) ) )

      DR(4,1:3) = CENTERS(1:3,5) - CENTERS(1:3,6)
      DR2(4)    = SQRT ( DOT_PRODUCT ( DR(4,1:3), DR(4,1:3) ) )

      VALUE(2)  = DR2(3) - DR2(4)

      DF = .0_DP
      CALL CALC_SPLINE2D1( CURRENT%NX, CURRENT%X, CURRENT%NY, CURRENT%Y, &
                           CURRENT%Z, CURRENT%Z21, VALUE(1), VALUE(2), E0(1), DF(1) )
      CALL CALC_SPLINE2D2( CURRENT%NX, CURRENT%X, CURRENT%NY, CURRENT%Y, &
                           CURRENT%Z, CURRENT%Z22, VALUE(1), VALUE(2), E0(2), DF(2) )
      ECORR = ECORR + E0(1) * .5_DP + E0(2) * .5_DP

      ! . If one of the gradient vectors is zero (out of range), ignore the other ...
      IF( DF(1) == .0_DP ) THEN
          DF(2) = .0_DP
      ELSE 
          IF( DF(2) == .0_DP ) DF(1) = .0_DP
      END IF

write(*,'(a,6f12.6)') 'S2D', value(1), value(2), e0(1), e0(2), df(1), df(2)

      IF( .NOT. QGRADIENT ) RETURN
      DEDC(1:3,1) =   DF(1) / DR2(1) * DR(1,1:3)
      DEDC(1:3,3) =   DF(1) / DR2(2) * DR(2,1:3)
      DEDC(1:3,2) = - DF(1) / DR2(1) * DR(1,1:3) - DF(1) / DR2(2) * DR(2,1:3)

      DEDC(1:3,4) =   DF(2) / DR2(3) * DR(3,1:3)
      DEDC(1:3,6) =   DF(2) / DR2(4) * DR(4,1:3)
      DEDC(1:3,5) = - DF(2) / DR2(3) * DR(3,1:3) - DF(2) / DR2(4) * DR(4,1:3)

      END SUBROUTINE BURRITO_LOCO
! --------------------------------------------------------------------------------------

      !-------------------------------
      SUBROUTINE ENERGY_MULTIPLE_ANGLE
      !-------------------------------

      REAL ( KIND = DP ), PARAMETER :: DOT_LIMIT = 1.0_DP - 1.0E-6_DP
      REAL ( KIND = DP ), DIMENSION(1:2) :: DF, DOTFAC, RIJ, RKJ, E0, VALUE
      REAL ( KIND = DP ), DIMENSION(1:2,1:3) :: DRIJ, DRKJ, DTI, DTJ, DTK

      ! . Just Bidimensional stuff
      DRIJ(1,1:3) = CENTERS(1:3,1) - CENTERS(1:3,2)
      DRKJ(1,1:3) = CENTERS(1:3,3) - CENTERS(1:3,2)
      RIJ(1) = SQRT ( DOT_PRODUCT ( DRIJ(1,1:3), DRIJ(1,1:3) ) )
      RKJ(1) = SQRT ( DOT_PRODUCT ( DRKJ(1,1:3), DRKJ(1,1:3) ) )
      DRIJ(1,1:3) = DRIJ(1,1:3) / RIJ(1)
      DRKJ(1,1:3) = DRKJ(1,1:3) / RKJ(1)
      DOTFAC(1) = DOT_PRODUCT ( DRIJ(1,1:3), DRKJ(1,1:3) )
      DOTFAC(1) = SIGN ( MIN ( ABS ( DOTFAC(1) ), DOT_LIMIT ), DOTFAC(1) )
      VALUE(1) = ACOS ( DOTFAC(1) ) * TO_DEGREES

      DRIJ(2,1:3) = CENTERS(1:3,4) - CENTERS(1:3,5)
      DRKJ(2,1:3) = CENTERS(1:3,6) - CENTERS(1:3,5)
      RIJ(2) = SQRT ( DOT_PRODUCT ( DRIJ(2,1:3), DRIJ(2,1:3) ) )
      RKJ(2) = SQRT ( DOT_PRODUCT ( DRKJ(2,1:3), DRKJ(2,1:3) ) )
      DRIJ(2,1:3) = DRIJ(2,1:3) / RIJ(2)
      DRKJ(2,1:3) = DRKJ(2,1:3) / RKJ(2)
      DOTFAC(2) = DOT_PRODUCT ( DRIJ(2,1:3), DRKJ(2,1:3) )
      DOTFAC(2) = SIGN ( MIN ( ABS ( DOTFAC(2) ), DOT_LIMIT ), DOTFAC(2) )
      VALUE(2) = ACOS ( DOTFAC(2) ) * TO_DEGREES

      CALL CALC_SPLINE2D1( CURRENT%NX, CURRENT%X, CURRENT%NY, CURRENT%Y, &
                           CURRENT%Z, CURRENT%Z21, VALUE(1), VALUE(2), E0(1), DF(1) )
      CALL CALC_SPLINE2D2( CURRENT%NX, CURRENT%X, CURRENT%NY, CURRENT%Y, &
                           CURRENT%Z, CURRENT%Z22, VALUE(1), VALUE(2), E0(2), DF(2) )
      ECORR = ECORR + E0(1) * .5_DP + E0(2) * .5_DP

      IF( .NOT. QGRADIENT ) RETURN
      DF(1)       = - TO_DEGREES * DF(1) / SQRT ( 1.0_DP - DOTFAC(1) * DOTFAC(1) )
      DTI(1,1:3)  = ( DRKJ(1,1:3) - DOTFAC(1) * DRIJ(1,1:3) ) / RIJ(1)
      DTK(1,1:3)  = ( DRIJ(1,1:3) - DOTFAC(1) * DRKJ(1,1:3) ) / RKJ(1)
      DTJ(1,1:3)  = - ( DTI(1,1:3) + DTK(1,1:3) )
      DEDC(1:3,1) = DF(1) * DTI(1,1:3)
      DEDC(1:3,2) = DF(1) * DTJ(1,1:3)
      DEDC(1:3,3) = DF(1) * DTK(1,1:3)

      DF(2)       = - TO_DEGREES * DF(2) / SQRT ( 1.0_DP - DOTFAC(2) * DOTFAC(2) )
      DTI(2,1:3)  = ( DRKJ(2,1:3) - DOTFAC(2) * DRIJ(2,1:3) ) / RIJ(2)
      DTK(2,1:3)  = ( DRIJ(2,1:3) - DOTFAC(2) * DRKJ(2,1:3) ) / RKJ(2)
      DTJ(2,1:3)  = - ( DTI(2,1:3) + DTK(2,1:3) )
      DEDC(1:3,4) = DF(2) * DTI(2,1:3)
      DEDC(1:3,5) = DF(2) * DTJ(2,1:3)
      DEDC(1:3,6) = DF(2) * DTK(2,1:3)

      END SUBROUTINE ENERGY_MULTIPLE_ANGLE

   END SUBROUTINE ENERGY_SPLINEX

   !--------------------------------------------------
   SUBROUTINE SPLINE_POINT_DEFINE ( SELECTION, PRINT )
   !--------------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Array arguments.
   LOGICAL,            DIMENSION(:), INTENT(IN)           :: SELECTION

   ! . Local scalars.
   INTEGER                          :: I, N, NSELECTED
   LOGICAL                          :: QPRINT
   TYPE(LINKED_SPOINT_TYPE), POINTER :: PLOCAL

   ! . Check to see if there are any atoms.
   IF ( NATOMS <= 0 ) RETURN

   ! . Check the selection array size.
   IF ( SIZE ( SELECTION ) /= NATOMS ) CALL PRINT_ERROR ( "SPLINE_POINT_DEFINE", "Invalid selection array size." )

   ! . Get the number of atoms in each selection.
   NSELECTED = COUNT ( SELECTION )

   ! . Check the number of atoms.
   IF ( NSELECTED <= 0 ) CALL PRINT_ERROR ( "SPLINE_POINT_DEFINE", "No atoms have been selected." )

   ! . Allocate the new constraint point.
   ALLOCATE ( PLOCAL )

   ! . Allocate the point type.
   ALLOCATE ( PLOCAL%POINT )

   ! . Set the number of atoms in the point.
   PLOCAL%POINT%NATOMS = NSELECTED

   ! . Allocate the constraint point arrays.
   ALLOCATE ( PLOCAL%POINT%INDICES(1:NSELECTED) )

   ! . Fill the point index array.
   N = 0
   DO I = 1,NATOMS
      IF ( SELECTION(I) ) THEN
         N = N + 1
         PLOCAL%POINT%INDICES(N) = I
      END IF
   END DO

   ! . Nullify the pointer in WLOCAL.
   NULLIFY ( PLOCAL%NEXT_POINT )

   ! . Increment the number of constraint points.
   NSPOINTS = NSPOINTS + 1

   ! . This is the first point on the list.
   IF ( NSPOINTS == 1 ) THEN

      ! . Define the first constraint point pointer.
      SPOINT_FIRST => PLOCAL

   ! . There are other constraint points on the list.
   ELSE

      ! . Set NEXT_POINT of the old last constraint point.
      SPOINT_LAST%NEXT_POINT => PLOCAL

   END IF

   ! . Update the last constraint point pointer.
   SPOINT_LAST => PLOCAL

   ! . Get the print flag.
   IF ( PRESENT ( PRINT ) ) THEN
      QPRINT = PRINT
   ELSE
      QPRINT = .TRUE.
   END IF

   ! . Write out some information about the point.
   IF ( QPRINT ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "0000FF" )
      CALL PRINT_SUMMARY_START ( "Spline Data" )
      WRITE( PRINT_LINE, "(I6)" ) NSPOINTS;  CALL PRINT_SUMMARY_ELEMENT( "Current Point Number" )
      WRITE( PRINT_LINE, "(I6)" ) NSELECTED; CALL PRINT_SUMMARY_ELEMENT( "Number of Atoms in Point" )
      CALL PRINT_SUMMARY_STOP
   END IF

   END SUBROUTINE SPLINE_POINT_DEFINE

#endif

END MODULE SPLINE
