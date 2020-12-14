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
!                        The Normal Mode Analysis Module
!===============================================================================
!
! . Subroutines:
!
!   NORMAL_MODE_FREQUENCIES         Calculate the harmonic frequencies.
!   NORMAL_MODE_INITIALIZE          Initialize the normal mode data structures.
!   NORMAL_MODE_INTENSITIES         Calculate the infrared intensities.
!   NORMAL_MODE_PRINT               Print information about the modes.
!   NORMAL_MODE_TRAJECTORY          Generate a trajectory for a single mode.
!
!===============================================================================
MODULE NORMAL_MODE

! . Module declarations.
USE CONSTANTS,             ONLY : KBOLTZ, NAVOGADRO, PI, TO_HZ, TO_KM_PER_MOLE, TO_WAVENUMBERS
USE DEFINITIONS,           ONLY : DP
USE ELEMENTS,              ONLY : SYMBOL
USE IO_UNITS,              ONLY : OUTPUT
USE PRINTING,              ONLY : PRINT_ERROR, PRINT_LINE, PRINT_NOFORMAT_START, PRINT_NOFORMAT_STOP, PRINT_PARAGRAPH

USE DIAGONALIZATION,       ONLY : SYMMETRIC_UPPER

USE ATOMS,                 ONLY : ATMFIX, ATMMAS, ATMNUM, NATOMS, NFIXED, NFREE
USE DCD_IO
USE MULTIPOLES,            ONLY : DIPOLE_DERIVATIVES
USE NORMAL_MODE_UTILITIES
USE POTENTIAL_ENERGY,      ONLY : ATMHES, CALCULATE_HESSIAN => HESSIAN
USE FILES
USE STRING,                ONLY : ENCODE_INTEGER

IMPLICIT NONE
PRIVATE
PUBLIC :: FREQUENCIES, INTENSITIES, NMODES, MODES, NORMAL_MODE_FREQUENCIES, NORMAL_MODE_INTENSITIES, &
                                                   NORMAL_MODE_PRINT,       NORMAL_MODE_TRAJECTORY, &
                                                   NORMAL_MODE_VIEW
#ifndef PGPC
SAVE
#endif

! . Module parameters.
REAL ( KIND = DP ), PARAMETER :: LOW_FREQUENCY = 0.1_DP

! . Module scalars.
INTEGER :: NMODES = 0

! . Module arrays.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: FREQUENCIES, INTENSITIES
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: MODES

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------------------------------
   SUBROUTINE NORMAL_MODE_FREQUENCIES ( HESSIAN, MODIFY, PRINT )
   !------------------------------------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: MODIFY
   LOGICAL,               INTENT(IN), OPTIONAL :: PRINT

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN), OPTIONAL :: HESSIAN

   ! . Local scalars.
   INTEGER :: I, IATOM, IFREE, II, N
   LOGICAL :: QMASS, QPRINT

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:) :: EIGVAL, H, WEIGHTS

   ! . Initialize the data structure.
   CALL NORMAL_MODE_INITIALIZE

   ! . Check the number of free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Get the dimension of the problem.
   NMODES = 3 * NFREE
   N      = ( NMODES * ( NMODES + 1 ) ) / 2

   ! . Allocate space for the arrays.
   ALLOCATE ( EIGVAL(1:NMODES), FREQUENCIES(1:NMODES), H(1:N), MODES(1:NMODES,1:NMODES), WEIGHTS(1:NMODES) )

   ! . A HESSIAN has been supplied.
   IF ( PRESENT ( HESSIAN ) ) THEN

      ! . Check the size of HESSIAN.
      IF ( SIZE ( HESSIAN ) /= N ) CALL PRINT_ERROR ( "NORMAL_MODE_FREQUENCIES", "Existing HESSIAN invalid." )

      ! . Save the HESSIAN.
      H = HESSIAN

   ! . A HESSIAN needs to be calculated.
   ELSE

      ! . Calculate and save the HESSIAN.
      CALL CALCULATE_HESSIAN ; H = ATMHES

   END IF

   ! . Mass-weighting is always true.
   QMASS = .TRUE.

   ! . Fill the weights array.
   IFREE = 0
   DO IATOM = 1,NATOMS
      IF ( .NOT. ATMFIX(IATOM) ) THEN
         IFREE = IFREE + 1
         II    = 3 * ( IFREE - 1 )
         WEIGHTS(II+1:II+3) = 1.0_DP / SQRT ( ATMMAS(IATOM) )
      END IF
   END DO

   ! . Scale the matrix elements.
   CALL MASS_WEIGHT_HESSIAN ( H, WEIGHTS )

   ! . Modify the second derivative matrix.
   IF ( PRESENT ( MODIFY ) ) THEN
      SELECT CASE ( MODIFY )
      CASE ( "PROJECT" ) ; CALL PROJECT_ROTATION_TRANSLATION ( H, QMASS )
      CASE ( "RAISE"   ) ; CALL   RAISE_ROTATION_TRANSLATION ( H, QMASS )
      CASE DEFAULT ; CALL PRINT_ERROR ( "NORMAL_MODE_FREQUENCIES", "Invalid MODIFY option." )
      END SELECT
   END IF

   ! . Diagonalize the matrix.
   CALL SYMMETRIC_UPPER ( H, EIGVAL, MODES )

   ! . Convert the eigenvalues to frequencies.
   FREQUENCIES = TO_WAVENUMBERS * SIGN ( SQRT ( ABS ( EIGVAL ) ), EIGVAL )

   ! . Un mass-weight the eigenvectors.
   DO I = 1,NMODES
      MODES(1:NMODES,I) = MODES(1:NMODES,I) * WEIGHTS(1:NMODES)
   END DO

   ! . Deallocate all temporary arrays.
   DEALLOCATE ( EIGVAL, H, WEIGHTS )

   ! . Check the print option.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Print out the frequencies if required.
   IF ( QPRINT ) THEN
      CALL PRINT_NOFORMAT_START
      WRITE ( OUTPUT, "(/24('-'),A,24('-'))" ) " Harmonic Frequencies (cm^(-1)) "
      WRITE ( OUTPUT, "(8F10.3)" ) FREQUENCIES(1:NMODES)
      WRITE ( OUTPUT, "(80('-'))" )
      CALL PRINT_NOFORMAT_STOP
   END IF

   END SUBROUTINE NORMAL_MODE_FREQUENCIES

   !--------------------------------
   SUBROUTINE NORMAL_MODE_INITIALIZE
   !--------------------------------

   ! . Deallocate the arrays.
   IF ( ALLOCATED ( FREQUENCIES ) ) DEALLOCATE ( FREQUENCIES )
   IF ( ALLOCATED ( INTENSITIES ) ) DEALLOCATE ( INTENSITIES )
   IF ( ALLOCATED ( MODES       ) ) DEALLOCATE ( MODES       )

   ! . Reset the number of modes counter.
   NMODES = 0

   END SUBROUTINE NORMAL_MODE_INITIALIZE

   !-------------------------------------------
   SUBROUTINE NORMAL_MODE_INTENSITIES ( PRINT )
   !-------------------------------------------

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT

   ! . Local scalars.
   INTEGER            :: I, IMODE
   LOGICAL            :: QPRINT
   REAL ( KIND = DP ) :: SUM

   ! . Local allocatable arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: DMUDR

   ! . Initialize the data structure.
   IF ( ALLOCATED ( INTENSITIES ) ) DEALLOCATE ( INTENSITIES )

   ! . Check the number of free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Check NMODES.
   IF ( ( NMODES /= 3 * NFREE ) .OR. ( NMODES == 0 ) ) THEN
      CALL PRINT_ERROR ( "NORMAL_MODE_INTENSITIES", "There are no modes available." )
   END IF

   ! . Allocate and initialize the intensity array.
   ALLOCATE ( INTENSITIES(1:NMODES) )
   INTENSITIES = 0.0_DP

   ! . Allocate and calculate the dipole derivatives array.
   ALLOCATE ( DMUDR(1:3*NATOMS,1:3) )
   CALL DIPOLE_DERIVATIVES ( DMUDR )

   ! . Contract the dipole derivative arrays.
   IF ( NFIXED > 0 ) THEN
      DO I = 1,3
         DMUDR(1:NMODES,I) = PACK ( DMUDR(1:3*NATOMS,I), MASK = ( DMUDR(1:3*NATOMS,I) /= 0.0_DP ) )
      END DO
   END IF

   ! . Loop over the modes.
   DO IMODE = 1,NMODES

      ! . Check the value of the frequency.
      IF ( ABS ( FREQUENCIES(IMODE) ) < LOW_FREQUENCY ) CYCLE

      ! . Loop over the components of the dipole derivative array.
      SUM = 0.0_DP
      DO I = 1,3
         SUM = SUM + DOT_PRODUCT ( MODES(1:NMODES,IMODE), DMUDR(1:NMODES,I) )**2
      END DO

      ! . Calculate the intensity.
      INTENSITIES(IMODE) = TO_KM_PER_MOLE * SUM

   END DO

   ! . Deallocate the temporary array.
   DEALLOCATE ( DMUDR )

   ! . Check the print option.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Print out the frequencies if required.
   IF ( QPRINT ) THEN
      CALL PRINT_NOFORMAT_START
      WRITE ( OUTPUT, "(/23('-'),A,23('-'))" ) " Infrared Intensities (km mol^-1) "
      WRITE ( OUTPUT, "(8F10.3)" ) INTENSITIES(1:NMODES)
      WRITE ( OUTPUT, "(80('-'))" )
      CALL PRINT_NOFORMAT_STOP
   END IF

   END SUBROUTINE NORMAL_MODE_INTENSITIES

   !------------------------------------------------------
   SUBROUTINE NORMAL_MODE_PRINT ( START, STOP, SELECTION )
   !------------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN), OPTIONAL :: START, STOP

   ! . Array arguments.
   LOGICAL, DIMENSION(1:NATOMS), INTENT(IN), OPTIONAL :: SELECTION

   ! . Local parameters.
   INTEGER, PARAMETER :: NCOL = 6

   ! . Local scalars.
   INTEGER :: BEG, END, I, IATOM, IFREE, II, N, NOUT
   LOGICAL :: QINTENS

   ! . Check NMODES.
   IF ( NMODES <= 0 ) RETURN

   ! . Check the starting and stopping positions.
   BEG = 1
   END = NMODES   
   IF ( PRESENT ( START ) ) BEG = START
   IF ( PRESENT ( STOP  ) ) END = STOP
   NOUT = END - BEG + 1

   ! . Check the number of modes to print.
   IF ( NOUT <= 0 ) RETURN

   ! . Check for the presence of intensities.
   QINTENS = ALLOCATED ( INTENSITIES )

   ! . Print out the header.
   CALL PRINT_NOFORMAT_START
   WRITE ( OUTPUT, "(/27('-'),A,27('-'))" ) " Normal Mode Eigenvectors "

   ! . Initialization.
   N   = BEG - 1
   BEG = BEG - NCOL

   ! . Top of the loop over printing.
   10 CONTINUE
   BEG = BEG + NCOL
   N   = MIN ( ( N + NCOL ), END )

   ! . Write out the first few rows.
   WRITE ( OUTPUT, "(/A,13X,6(6X,I4))" ) "Mode", ( I, I = BEG,N )
   WRITE ( OUTPUT,    "(A,6F11.3)" ) "Freq. (cm^-1) ", ( FREQUENCIES(I),         I = BEG,N )
   WRITE ( OUTPUT, "(6X,A,6F11.3)" )       "(ps^-1) ", ( FREQUENCIES(I) * TO_HZ, I = BEG,N )
   IF ( QINTENS ) THEN
      WRITE ( OUTPUT, "(A,5X,6F11.3)" ) "Intensity", ( INTENSITIES(I), I = BEG,N )
      WRITE ( OUTPUT, "(A)" ) "(km mol^-1)"
   END IF

   ! . Print out the modes requested.
   IF ( PRESENT ( SELECTION ) ) THEN
      IFREE = 0
      DO IATOM = 1,NATOMS
         IF ( .NOT. ATMFIX(IATOM) ) THEN
            IFREE = IFREE + 1
            IF ( SELECTION(IATOM) ) THEN
               II = 3 * ( IFREE - 1 )
               WRITE ( OUTPUT, "(I4,2X,A2,2X,A,3X,6F11.5)" ) IATOM, SYMBOL(ATMNUM(IATOM)), "x", ( MODES(II+1,I), I = BEG,N )
               WRITE ( OUTPUT,         "(10X,A,3X,6F11.5)" )                               "y", ( MODES(II+2,I), I = BEG,N )
               WRITE ( OUTPUT,         "(10X,A,3X,6F11.5)" )                               "z", ( MODES(II+3,I), I = BEG,N )
            END IF
         END IF
      END DO
   ELSE
      IFREE = 0
      DO IATOM = 1,NATOMS
         IF ( .NOT. ATMFIX(IATOM) ) THEN
            IFREE = IFREE + 1
            II    = 3 * ( IFREE - 1 )
            WRITE ( OUTPUT, "(I4,2X,A2,2X,A,3X,6F11.5)" ) IATOM, SYMBOL(ATMNUM(IATOM)), "x", ( MODES(II+1,I), I = BEG,N )
            WRITE ( OUTPUT,         "(10X,A,3X,6F11.5)" )                               "y", ( MODES(II+2,I), I = BEG,N )
            WRITE ( OUTPUT,         "(10X,A,3X,6F11.5)" )                               "z", ( MODES(II+3,I), I = BEG,N )
         END IF
      END DO
   END IF

   ! . Check to see if there are more modes to print.
   IF ( N < NOUT ) GO TO 10

   ! . Print out the terminator.
   WRITE ( OUTPUT, "(80('-'))" )
   CALL PRINT_NOFORMAT_STOP

   END SUBROUTINE NORMAL_MODE_PRINT

   !--------------------------------------------------------------------
   SUBROUTINE NORMAL_MODE_TRAJECTORY ( FILE, MODE, NCYCLES, NFRAMES, T )
   !--------------------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE
   INTEGER,               INTENT(IN) :: MODE, NCYCLES, NFRAMES
   REAL ( KIND = DP ),    INTENT(IN) :: T

   ! . Local scalars.
   INTEGER            :: IATOM, ICYCLE, IFRAME, IFREE, II, NTOTAL
   REAL ( KIND = DP ) :: AMPLITUDE, FACTOR, OMEGA

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: FULL_MODE

   ! . Local types.
   TYPE(DCD_TYPE) :: DCD

   ! . Calculate the number of frames.
   NTOTAL = NCYCLES * NFRAMES

   ! . Check the modes and the number of frames.
   IF ( ( NMODES <= 0 ) .OR. ( NMODES /= 3 * NFREE ) .OR. ( MODE <= 0 ) .OR. ( MODE > NMODES ) .OR. ( NTOTAL <= 0 ) ) RETURN

   ! . Check the frequency of the mode.
   OMEGA = ABS ( FREQUENCIES(MODE) )
   IF ( OMEGA < LOW_FREQUENCY ) RETURN

   ! . Initialize and activate the trajectory.
   CALL DCD_INITIALIZE ( DCD )
   CALL DCD_ACTIVATE_WRITE ( FILE, DCD, "CORD", NATOMS, NFIXED, NTOTAL, QFIX = ATMFIX )

   ! . Calculate the amplitude (in Angstroms).
   AMPLITUDE = SQRT ( 2.0_DP * 1.0E-3_DP * KBOLTZ * NAVOGADRO * T ) * ( TO_WAVENUMBERS / OMEGA )

   ! . Fill FULL_MODE.
   FULL_MODE = 0.0_DP
   IFREE     = 0
   DO IATOM = 1,NATOMS
      IF ( .NOT. ATMFIX(IATOM) ) THEN
         IFREE = IFREE + 1
         II    = 3 * ( IFREE - 1 )
         FULL_MODE(1:3,IATOM) = MODES(II+1:II+3,MODE)
      END IF
   END DO

   ! . Loop over the cycles.
   DO ICYCLE = 1,NCYCLES

      ! . Loop over the frames.
      DO IFRAME = 1,NFRAMES

         ! . Calculate the displacement prefactor (use sine instead of cosine).
         FACTOR = AMPLITUDE * SIN ( 2.0_DP * PI * REAL ( ( IFRAME - 1 ), DP ) / REAL ( NFRAMES, DP ) )

         ! . Write out the coordinates.
         CALL DCD_WRITE ( DCD, ( ATMCRD + FACTOR * FULL_MODE ) )

      END DO

   END DO

   ! . Deactivate the trajectory.
   CALL DCD_DEACTIVATE ( DCD )

   ! . Write out some information.
   WRITE ( PRINT_LINE, "(A,I4,A)" ) "Trajectory written for mode ", MODE, " to "//FILE(1:LEN_TRIM(FILE))
   CALL PRINT_PARAGRAPH

   END SUBROUTINE NORMAL_MODE_TRAJECTORY


   ! ------------------------------------------------
   SUBROUTINE NORMAL_MODE_VIEW ( MODE, TNEED, AFACT )
   ! ------------------------------------------------
   IMPLICIT NONE
 
   INTEGER,                         INTENT(IN) :: MODE
   REAL ( KIND = DP ),    OPTIONAL, INTENT(IN) :: TNEED
   REAL ( KIND = DP ),    OPTIONAL, INTENT(IN) :: AFACT
 
   INTEGER                                     :: IATOM, ICYCLE, IFRAME, IFREE, II, FD
   REAL ( KIND = DP )                          :: AMPLITUDE, FACTOR, OMEGA, TEMP, RF
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS) :: FULL_MODE
   CHARACTER( LEN=256 )                        :: FN
 
   INTEGER,            PARAMETER :: NCYCLES = 10
   INTEGER,            PARAMETER :: NFRAMES = 10
   REAL ( KIND = DP ), PARAMETER :: LOW_FREQUENCY = 0.1_DP
 
   RF = 1._DP
   IF( PRESENT( AFACT ) ) RF = AFACT
   TEMP = 298.15_DP
   IF( PRESENT( TNEED ) ) TEMP = TNEED
   OMEGA  = ABS ( FREQUENCIES(MODE) )
   IF ( OMEGA < LOW_FREQUENCY ) RETURN
   AMPLITUDE = SQRT ( 2.0_DP * 1.0E-3_DP * KBOLTZ * NAVOGADRO * TEMP ) * &
       ( TO_WAVENUMBERS / OMEGA )
   FULL_MODE = 0.0_DP
   IFREE     = 0
   DO IATOM = 1, NATOMS
      IF ( .NOT. ATMFIX(IATOM) ) THEN
         IFREE = IFREE + 1
         II    = 3 * ( IFREE - 1 )
         FULL_MODE(1:3,IATOM) = MODES(II+1:II+3,MODE)
      END IF
   END DO
   FD = NEXT_UNIT()
   CALL ENCODE_INTEGER( MODE, FN, "(I10)" )
   FN = "nmode." // TRIM( FN )
   OPEN( UNIT = FD, FILE = TRIM( FN ), ACTION = "WRITE", FORM = "FORMATTED", STATUS = "UNKNOWN" )
   DO ICYCLE = 1, NCYCLES
      DO IFRAME = 1, NFRAMES
         WRITE ( FD, "(I10)" ) NFREE
         WRITE ( FD, "(F10.3,A)" ) FREQUENCIES(MODE), " cm^-1"
         FACTOR = RF * AMPLITUDE * &
             SIN ( 2.0_DP * PI * REAL ( ( IFRAME - 1 ), DP ) / REAL ( NFRAMES, DP ) )
         DO IATOM = 1, NATOMS
             IF ( .NOT. ATMFIX(IATOM) ) THEN
                 WRITE ( FD, "(A,8X,3F20.12)" ) &
                     SYMBOL(ATMNUM(IATOM)),          &
                     ATMCRD(1:3,IATOM) + FACTOR * FULL_MODE(1:3,IATOM)
             END IF
         END DO
      END DO
   END DO
   CLOSE( FD )
   END SUBROUTINE NORMAL_MODE_VIEW

END MODULE NORMAL_MODE
