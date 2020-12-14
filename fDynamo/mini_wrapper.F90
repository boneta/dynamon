MODULE MINI_WRAPPER

USE DEFINITIONS,           ONLY : DP
USE ATOMS,                 ONLY : ATMCRD, ATMFIX, NATOMS, NFIXED, NFREE
USE POTENTIAL_ENERGY,      ONLY : ATMDER, GRADIENT, ENERGY, ETOTAL
USE IO_UNITS,              ONLY : OUTPUT
USE FILES

IMPLICIT NONE
PRIVATE
PUBLIC  :: OPTIMIZE_LBFGSB, LBFGSB_MINIMIZE, OPTIMIZE_CGP, CGP_MINIMIZE


CONTAINS

subroutine start_trj( idcd, name, iframes )
	integer, intent(in)          :: idcd
	character(len=*), intent(in) :: name
	integer, intent(in)          :: iframes
	integer, dimension(1:20)     :: icntrl
	open( unit=idcd, file=trim( name ), action="write", form="unformatted" )
	icntrl = (/ iframes, 1, 1, iframes, 0, 0, 0, 3*natoms, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
	write( idcd ) "CORD", icntrl(1:20)
	write( idcd ) 1, repeat( " ", 80 )
	write( idcd ) natoms
	call flush( idcd )
end subroutine
subroutine write_trj( idcd )
	integer, intent(in) :: idcd
	integer             :: i
	write( idcd ) (real(atmcrd(1,i),4),i=1,natoms)
	write( idcd ) (real(atmcrd(2,i),4),i=1,natoms)
	write( idcd ) (real(atmcrd(3,i),4),i=1,natoms)
	call flush( idcd )
end subroutine
subroutine close_trj( idcd )
	integer, intent(in) :: idcd
	call flush( idcd )
	close( idcd )
end subroutine

   SUBROUTINE OPTIMIZE_LBFGSB( STEP_NUMBER, GRADIENT_TOLERANCE, PRINT_FREQUENCY, BRACKET )
   IMPLICIT NONE
   INTEGER,                           INTENT(IN)   :: STEP_NUMBER, PRINT_FREQUENCY
   REAL ( KIND=DP ),                  INTENT(IN)   :: GRADIENT_TOLERANCE
   REAL( KIND=DP ),         INTENT(IN), OPTIONAL   :: BRACKET

   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: X
   IF ( NFREE < 0 ) RETURN
   ALLOCATE ( X(1:3*NFREE) )
   CALL DAT_TO_X( ATMCRD, X )
   CALL LBFGSB_MINIMIZE( EGCALC, X, STEP_NUMBER, GRADIENT_TOLERANCE, PRINT_FREQUENCY, BRACKET )
   CALL X_TO_DAT( X, ATMCRD )
   END SUBROUTINE OPTIMIZE_LBFGSB

   SUBROUTINE LBFGSB_MINIMIZE( FCALC, X, STEP_NUMBER, GRADIENT_TOLERANCE, PRINT_FREQUENCY, BRACKET )
   IMPLICIT NONE
   INTERFACE
      SUBROUTINE FCALC ( X, F, G )
         USE DEFINITIONS, ONLY : DP
         REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)  :: X
         REAL ( KIND = DP ),               INTENT(OUT) :: F
         REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: G
      END SUBROUTINE FCALC
   END INTERFACE
 
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: X
   INTEGER,                           INTENT(IN)   :: STEP_NUMBER, PRINT_FREQUENCY
   REAL ( KIND=DP ),                  INTENT(IN)   :: GRADIENT_TOLERANCE
   REAL( KIND=DP ),         INTENT(IN), OPTIONAL   :: BRACKET
 
   CHARACTER ( LEN=60 )                        :: TASK, CSAVE
   LOGICAL, DIMENSION(1:4)                     :: LSAVE
   INTEGER                                     :: NMAX, MMAX, ISTEP
   INTEGER, DIMENSION(1:44)                    :: ISAVE
   INTEGER, ALLOCATABLE, DIMENSION(:)          :: NBD, IWA
   REAL ( KIND=DP )                            :: FACTR, CGRMS, CENER, GRMS_BEST, ENER_BEST, SQRTn
   REAL ( KIND=DP ), DIMENSION(1:29)           :: DSAVE
   REAL ( KIND=DP ), ALLOCATABLE, DIMENSION(:) :: X_BEST, L, U, G, WA
 
   NMAX = SIZE( X )
   ALLOCATE( NBD(1:NMAX), IWA(1:3*NMAX) )
   ALLOCATE( X_BEST(1:NMAX), L(1:NMAX), U(1:NMAX), G(1:NMAX) )
   SQRTn = SQRT( REAL( NMAX, KIND=DP ) )
 
   MMAX = 5
   MMAX = 9
   ALLOCATE( WA(1:2*MMAX*NMAX+4*NMAX+12*MMAX*MMAX+12*MMAX) )
 
   CALL FCALC( X, CENER, G )
   X_BEST    = X
   ENER_BEST = CENER
   FACTR     = 1.E+3_DP
!   FACTR     = 1.E+1_DP
   ISTEP     =  1
   IF( PRESENT( BRACKET ) ) THEN
       NBD =  2
       L   = X - BRACKET
       U   = X + BRACKET
       WRITE( OUTPUT, "(/10('-'),A,F8.3)" ) " Optimization (L-BFGS-B) [Constrained]: ", BRACKET
   ELSE
       NBD =  0
       L   = .0_DP
       U   = .0_DP
       WRITE( OUTPUT, "(/10('-'),A)" ) " Optimization (L-BFGS-B) [UNConstrained]"
   END IF
   CGRMS     = SQRT( SUM( G * G ) ) / SQRTn
   GRMS_BEST = CGRMS
   TASK      = "START"
 
   WRITE( OUTPUT, "(A10,2A20)" ) "Step", "Function", "RMS Gradient"
   WRITE( OUTPUT, "(50('-'))" )
   WRITE( OUTPUT, "(I10,2F20.8)" ) 0, CENER, CGRMS
 
   LBFGSB_LOOP: DO WHILE( ISTEP <= STEP_NUMBER       .AND. &
                           CGRMS > GRADIENT_TOLERANCE .AND. &
                           TASK(1:4) /= "CONV"        .AND. &
                           TASK(1:4) /= "STOP"        .AND. &
                           TASK(1:5) /= "ERROR"       .AND. &
                           TASK(1:4) /= "ABNO"              )
 
       CALL SETULB( NMAX, MMAX, X, L, U, NBD, CENER, G, FACTR, GRADIENT_TOLERANCE, &
                     WA, IWA, TASK, -1, CSAVE, LSAVE, ISAVE, DSAVE )
 
       IF( TASK(1:2) == "FG" ) THEN
           CALL FCALC( X, CENER, G )
           CGRMS = SQRT( SUM( G * G ) ) / SQRTn
!           IF( CGRMS < GRMS_BEST ) THEN
           IF( CENER < ENER_BEST ) THEN
               X_BEST    = X
               GRMS_BEST = CGRMS
               ENER_BEST = CENER
           END IF
           IF( MOD( ISTEP, PRINT_FREQUENCY ) == 0 ) WRITE( OUTPUT, "(I10,2F20.8)" ) ISTEP, CENER, CGRMS
           ISTEP = ISTEP + 1
       END IF
 
   END DO LBFGSB_LOOP
 
   IF ( MOD( ISTEP, PRINT_FREQUENCY ) /= 0 ) THEN
       WRITE( OUTPUT, "(I10,2F20.8)" ) ISTEP, CENER, CGRMS
   END IF
 
   WRITE( OUTPUT, "(50('-'))" )
  
   IF( TASK(1:4) == "ABNO"  .OR. &
       TASK(1:5) == "ERROR" .OR. &
       TASK(1:4) == "STOP"       ) THEN
       WRITE( OUTPUT, "(/A)" ) TASK
       WRITE( OUTPUT, "(A)" ) "LBFGSB_MINIMIZE: Restoring best function movement..."
       WRITE( OUTPUT, "(10X,2F20.8/)" ) ENER_BEST, GRMS_BEST
       X = X_BEST
   END IF
 
   DEALLOCATE( NBD, IWA )
   DEALLOCATE( X_BEST, L, U, G )
   DEALLOCATE( WA )
 
   END SUBROUTINE LBFGSB_MINIMIZE


   SUBROUTINE OPTIMIZE_CGP( STEP_NUMBER, GRADIENT_TOLERANCE, PRINT_FREQUENCY )
   IMPLICIT NONE
   INTEGER,                           INTENT(IN)   :: STEP_NUMBER, PRINT_FREQUENCY
   REAL ( KIND=DP ),                  INTENT(IN)   :: GRADIENT_TOLERANCE

   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: X
   IF ( NFREE < 0 ) RETURN
   ALLOCATE ( X(1:3*NFREE) )
   CALL DAT_TO_X( ATMCRD, X )
   CALL CGP_MINIMIZE( EGCALC, X, STEP_NUMBER, GRADIENT_TOLERANCE, PRINT_FREQUENCY )
   CALL X_TO_DAT( X, ATMCRD )
   END SUBROUTINE OPTIMIZE_CGP


   SUBROUTINE CGP_MINIMIZE( FCALC, X, STEP_NUMBER, GRADIENT_TOLERANCE, PRINT_FREQUENCY )
   IMPLICIT NONE
   INTERFACE
      SUBROUTINE FCALC ( X, F, G )
         USE DEFINITIONS, ONLY : DP
         REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)  :: X
         REAL ( KIND = DP ),               INTENT(OUT) :: F
         REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: G
      END SUBROUTINE FCALC
   END INTERFACE
 
   REAL ( KIND = DP ), DIMENSION(:), INTENT(INOUT) :: X
   INTEGER,                           INTENT(IN)   :: STEP_NUMBER, PRINT_FREQUENCY
   REAL ( KIND=DP ),                  INTENT(IN)   :: GRADIENT_TOLERANCE

   INTEGER          :: NX, ISTEP, IFLAG, IREST, IMETH
   REAL ( KIND=DP ) :: GRMS, FUNC, SQRTn, FB
   REAL ( KIND = DP ), DIMENSION(:), ALLOCATABLE :: XB, G, D, O, W

   NX    = SIZE( X )
   SQRTn = SQRT( REAL( NX, KIND=DP ) )
   ALLOCATE( G(1:NX), D(1:NX), O(1:NX), W(1:NX) )
   CALL FCALC( X, FUNC, G )
   XB    = X
   FB    = FUNC
   GRMS  = SQRT( SUM( G * G ) ) / SQRTn
   WRITE( OUTPUT, "(/10('-'),A)" ) " Optimization (CG+)"
   WRITE( OUTPUT, "(A10,2A20)" ) "Step", "Function", "RMS Gradient"
   WRITE( OUTPUT, "(50('-'))" )
   WRITE( OUTPUT, "(I10,2F20.8)" ) 0, FUNC, GRMS
   ISTEP = 0
   IFLAG = 0
!     IREST  =  0 (NO RESTARTS); 1 (RESTART EVERY N STEPS)
   IREST = 1
!     METHOD =  1 : FLETCHER-REEVES 
!               2 : POLAK-RIBIERE
!               3 : POSITIVE POLAK-RIBIERE ( BETA=MAX{BETA,0} )
   IMETH = 3
   CGP_LOOP: DO WHILE( ISTEP <= STEP_NUMBER       .AND. &
                       GRMS > GRADIENT_TOLERANCE )
      CALL CGFAM( NX, X, FUNC, G, D, O, GRADIENT_TOLERANCE, W, IFLAG, IREST, IMETH )
      IF( IFLAG == -3 ) THEN
          WRITE( OUTPUT, "(/A/)" ) "CGP_MINIMIZE: Improper input parameters... restoring best try."
          X = XB
          ISTEP = STEP_NUMBER + 1
      ELSE IF( IFLAG == -2 ) THEN
          WRITE( OUTPUT, "(/A/)" ) "CGP_MINIMIZE: Descent was not obtained... restoring best try."
          X = XB
          ISTEP = STEP_NUMBER + 1
      ELSE IF( IFLAG == -2 ) THEN
          WRITE( OUTPUT, "(/A/)" ) "CGP_MINIMIZE: Line search failure... restoring best try."
          X = XB
          ISTEP = STEP_NUMBER + 1
      ELSE
          DO WHILE( IFLAG == 2 )
              CALL CGFAM( NX, X, FUNC, G, D, O, GRADIENT_TOLERANCE, W, IFLAG, IREST, IMETH )
          END DO
          CALL FCALC( X, FUNC, G )
          GRMS  = SQRT( SUM( G * G ) ) / SQRTn
          IF( FUNC < FB ) THEN
              XB = X
              FB = FUNC
          END IF
      END IF
      IF( MOD( ISTEP, PRINT_FREQUENCY ) == 0 ) WRITE( OUTPUT, "(I10,2F20.8)" ) ISTEP, FUNC, GRMS
      ISTEP = ISTEP + 1
   END DO CGP_LOOP
   DEALLOCATE( XB, G, D, O, W )
   END SUBROUTINE CGP_MINIMIZE


! - common subroutines...

   SUBROUTINE X_TO_DAT( X, ATMDAT )
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(IN)  :: X
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT) :: ATMDAT
   INTEGER :: IATOM, II
   II = -3
   DO IATOM = 1,NATOMS
      IF ( .NOT. ATMFIX(IATOM) ) THEN
         II = II + 3
         ATMDAT(1:3,IATOM) = X(II+1:II+3)
      END IF
   END DO
   END SUBROUTINE X_TO_DAT


   SUBROUTINE DAT_TO_X( ATMDAT, X )
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(IN)  :: ATMDAT
   REAL ( KIND = DP ), DIMENSION(:),   INTENT(OUT) :: X
   INTEGER :: IATOM, II
   II = -3
   DO IATOM = 1,NATOMS
      IF ( .NOT. ATMFIX(IATOM) ) THEN
         II = II + 3
         X(II+1:II+3) = ATMDAT(1:3,IATOM)
      END IF
   END DO
   END SUBROUTINE DAT_TO_X


   SUBROUTINE EGCALC ( X, E, G )
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN)  :: X
   REAL ( KIND = DP ),               INTENT(OUT) :: E
   REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: G

   CALL X_TO_DAT( X, ATMCRD )
   CALL GRADIENT ( PRINT = .FALSE. )
   E = ETOTAL
   CALL DAT_TO_X( ATMDER, G )
   END SUBROUTINE EGCALC


END MODULE MINI_WRAPPER
