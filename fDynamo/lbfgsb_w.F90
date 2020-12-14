MODULE LBFGSB_WRAPPER

USE DEFINITIONS,           ONLY : DP
USE ATOMS,                 ONLY : ATMCRD, ATMFIX, NATOMS, NFIXED, NFREE
USE POTENTIAL_ENERGY,      ONLY : ATMDER, GRADIENT, ENERGY, ETOTAL
!USE PRINTING,              ONLY : PRINT_LINE, PRINT_PARAGRAPH, PRINT_TABLE_ELEMENT, &
!			                      PRINT_TABLE_OPTIONS, PRINT_TABLE_START, PRINT_TABLE_STOP
USE IO_UNITS,              ONLY : OUTPUT
USE FILES

IMPLICIT NONE
PRIVATE
PUBLIC  :: OPTIMIZE_LBFGSB


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

  SUBROUTINE OPTIMIZE_LBFGSB( STEP_NUMBER, GRADIENT_TOLERANCE, PRINT_FREQUENCY, QMMAX, DCD, BRACKET )

  IMPLICIT NONE

  INTEGER,                           INTENT(IN) :: STEP_NUMBER, PRINT_FREQUENCY
  REAL ( KIND=DP ),                  INTENT(IN) :: GRADIENT_TOLERANCE
  INTEGER,                 INTENT(IN), OPTIONAL :: QMMAX
  CHARACTER( LEN=* ),      OPTIONAL, INTENT(IN) :: DCD
  REAL( KIND=DP ),         INTENT(IN), OPTIONAL :: BRACKET

  CHARACTER ( LEN=60 )                          :: TASK, CSAVE
  LOGICAL, DIMENSION(1:4)                       :: LSAVE
  INTEGER                                       :: NMAX, MMAX, CI, CJ, ISTEP, IDCD
  INTEGER, DIMENSION(1:44)                      :: ISAVE
  INTEGER, ALLOCATABLE, DIMENSION(:)            :: NBD, IWA
  REAL ( KIND=DP )                              :: FACTR, CGRMS, CENER, GRMS_BEST, ENER_BEST, SQRTn
  REAL ( KIND=DP ), DIMENSION(1:29)             :: DSAVE
  REAL ( KIND=DP ), ALLOCATABLE, DIMENSION(:)   :: X, X_BEST, L, U, G, WA
  LOGICAL                                       :: QDCD

  IF ( NFREE == 0 ) RETURN

  NMAX = NFREE * 3
  ALLOCATE( NBD(1:NMAX), IWA(1:3*NMAX) )
  ALLOCATE( X(1:NMAX), X_BEST(1:NMAX), L(1:NMAX), U(1:NMAX), G(1:NMAX) )
  SQRTn = SQRT( REAL( NMAX, KIND=DP ) )

  MMAX = 5
  IF( PRESENT( QMMAX ) ) MMAX = QMMAX
  ALLOCATE( WA(1:2*MMAX*NMAX+4*NMAX+12*MMAX*MMAX+12*MMAX) )

  CALL GRADIENT( .FALSE. )

     CJ = 1
  DO CI = 1, NATOMS
     IF( .NOT. ATMFIX(CI) )  THEN
         X(CJ  ) = ATMCRD(1,CI)
         X(CJ+1) = ATMCRD(2,CI)
         X(CJ+2) = ATMCRD(3,CI)
         G(CJ  ) = ATMDER(1,CI)
         G(CJ+1) = ATMDER(2,CI)
         G(CJ+2) = ATMDER(3,CI)
         CJ = CJ + 3
     END IF
  END DO

  X_BEST    = X

  CENER     = ETOTAL
  ENER_BEST = CENER
  FACTR     = 1.E+3_DP
  FACTR     = 1.E+1_DP
  ISTEP     =  1
  IF( PRESENT( BRACKET ) ) THEN
      NBD       =  2
      L         = X - BRACKET
      U         = X + BRACKET
      WRITE( 6, "(/10('-'),A,F8.3)" ) " Optimization (L-BFGS-B) [Constrained]: ", BRACKET
!      WRITE( PRINT_LINE, "(10('-'),A,F8.3)" ) " Optimization (L-BFGS-B) [Constrained]: ", BRACKET
!      CALL PRINT_PARAGRAPH( BOLD = .TRUE. )
  ELSE
      NBD       =  0
      L         = .0_DP
      U         = .0_DP
      WRITE( 6, "(/10('-'),A)" ) " Optimization (L-BFGS-B) [UNConstrained]"
!      WRITE( PRINT_LINE, "(10('-'),A)" ) " Optimization (L-BFGS-B) [UNConstrained]"
!      CALL PRINT_PARAGRAPH( BOLD = .TRUE. )
  END IF
  CGRMS     = SQRT( SUM( G * G ) ) / SQRTn
  GRMS_BEST = CGRMS

  TASK      = "START"

!  CALL PRINT_TABLE_OPTIONS ( COLUMNS = 3, HEADER_COLOR = "#55FF55", PAGEWIDTH = 50, VARIABLEWIDTHS = (/ 10, 20, 20 /) )
!  CALL PRINT_TABLE_START
!  CALL PRINT_TABLE_ELEMENT ( TEXT = "Step",         HEADER = .TRUE. )
!  CALL PRINT_TABLE_ELEMENT ( TEXT = "Function",     HEADER = .TRUE. )
!  CALL PRINT_TABLE_ELEMENT ( TEXT = "RMS Gradient", HEADER = .TRUE. )
!  WRITE ( PRINT_LINE, "(I10)"   ) 0     ; CALL PRINT_TABLE_ELEMENT
!  WRITE ( PRINT_LINE, "(F20.8)" ) CENER ; CALL PRINT_TABLE_ELEMENT
!  WRITE ( PRINT_LINE, "(F20.8)" ) CGRMS ; CALL PRINT_TABLE_ELEMENT
  WRITE( 6, "(A10,2A20)" ) "Step", "Function", "RMS Gradient"
  WRITE( 6, "(50('-'))" )
  WRITE( 6, "(I10,2F20.8)" ) 0, CENER, CGRMS

  QDCD = .FALSE.
  IF( PRESENT( DCD ) ) THEN
      QDCD = .TRUE.
      IDCD = NEXT_UNIT()
      call start_trj( IDCD, TRIM( DCD ), STEP_NUMBER )
  END IF

  LBFGSB_LOOP: DO WHILE( ISTEP <= STEP_NUMBER       .AND. &
                          CGRMS > GRADIENT_TOLERANCE .AND. &
                          TASK(1:4) /= "CONV"        .AND. &
                          TASK(1:4) /= "STOP"        .AND. &
                          TASK(1:5) /= "ERROR"       .AND. &
                          TASK(1:4) /= "ABNO"              )

      CALL SETULB( NMAX, MMAX, X, L, U, NBD, CENER, G, FACTR, GRADIENT_TOLERANCE, &
                    WA, IWA, TASK, -1, CSAVE, LSAVE, ISAVE, DSAVE )

      IF( TASK(1:2) == "FG" ) THEN
             CJ = 1
          DO CI = 1, NATOMS
             IF( .NOT. ATMFIX(CI) )  THEN
                 ATMCRD(1,CI) = X(CJ  )
                 ATMCRD(2,CI) = X(CJ+1)
                 ATMCRD(3,CI) = X(CJ+2)
                 CJ = CJ + 3
             END IF
          END DO
          CALL GRADIENT( .FALSE. )
             CJ = 1
          DO CI = 1, NATOMS
             IF( .NOT. ATMFIX(CI) )  THEN
                 G(CJ  ) = ATMDER(1,CI)
                 G(CJ+1) = ATMDER(2,CI)
                 G(CJ+2) = ATMDER(3,CI)
                 CJ = CJ + 3
             END IF
          END DO
          CENER = ETOTAL
          CGRMS = SQRT( SUM( G * G ) ) / SQRTn
!          IF( CGRMS < GRMS_BEST ) THEN
          IF( CENER < ENER_BEST ) THEN
              X_BEST    = X
              GRMS_BEST = CGRMS
              ENER_BEST = CENER
          END IF
          IF( MOD( ISTEP, PRINT_FREQUENCY ) == 0 ) THEN
!              WRITE ( PRINT_LINE, "(I10)"   ) ISTEP ; CALL PRINT_TABLE_ELEMENT
!              WRITE ( PRINT_LINE, "(F20.8)" ) CENER ; CALL PRINT_TABLE_ELEMENT
!              WRITE ( PRINT_LINE, "(F20.8)" ) CGRMS ; CALL PRINT_TABLE_ELEMENT
              WRITE( 6, "(I10,2F20.8)" ) ISTEP, CENER, CGRMS
              IF( QDCD ) CALL write_trj( IDCD )
          END IF
          ISTEP = ISTEP + 1
      END IF

  END DO LBFGSB_LOOP

  IF( QDCD ) CALL close_trj( IDCD )

  IF ( MOD( ISTEP, PRINT_FREQUENCY ) /= 0 ) THEN
!      WRITE ( PRINT_LINE, "(I10)"   ) ISTEP ; CALL PRINT_TABLE_ELEMENT
!      WRITE ( PRINT_LINE, "(F20.8)" ) CENER ; CALL PRINT_TABLE_ELEMENT
!      WRITE ( PRINT_LINE, "(F20.8)" ) CGRMS ; CALL PRINT_TABLE_ELEMENT
      WRITE( 6, "(I10,2F20.8)" ) ISTEP, CENER, CGRMS
  END IF

!  CALL PRINT_TABLE_STOP
  WRITE( 6, "(50('-'))" )
 
  IF( TASK(1:4) == "ABNO"  .OR. &
       TASK(1:5) == "ERROR" .OR. &
       TASK(1:4) == "STOP"       ) THEN
!      CALL PRINT_PARAGRAPH( TEXT = TASK )
!      CALL PRINT_PARAGRAPH( TEXT = "MINIMIZE_LBFGSB: Restoring best movement.", BOLD = .TRUE. )
       WRITE( 6, "(/A)" ) TASK
       WRITE( 6, "(A/)" ) "MINIMIZE_LBFGSB: Restoring best function movement."
         CJ = 1
      DO CI = 1, NATOMS
         IF( .NOT. ATMFIX(CI) )  THEN
             ATMCRD(1,CI) = X_BEST(CJ  )
             ATMCRD(2,CI) = X_BEST(CJ+1)
             ATMCRD(3,CI) = X_BEST(CJ+2)
             CJ = CJ + 3
         END IF
      END DO
      ETOTAL = ENER_BEST
  ELSE
         CJ = 1
      DO CI = 1, NATOMS
         IF( .NOT. ATMFIX(CI) )  THEN
             ATMCRD(1,CI) = X(CJ  )
             ATMCRD(2,CI) = X(CJ+1)
             ATMCRD(3,CI) = X(CJ+2)
             CJ = CJ + 3
         END IF
      END DO
      ETOTAL = CENER
  END IF

  DEALLOCATE( NBD, IWA )
  DEALLOCATE( X, X_BEST, L, U, G )
  DEALLOCATE( WA )

  END SUBROUTINE OPTIMIZE_LBFGSB

END MODULE LBFGSB_WRAPPER
