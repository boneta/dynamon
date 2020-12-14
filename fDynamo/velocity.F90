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
!                              The Velocity Module
!===============================================================================
!
! . Scalars:
!
!   NDEGF                           The current number of degrees of freedom.
!   NPROJECT                        The number of translational projection
!                                   vectors.
!
!   EKE                             The kinetic energy.
!   TEMPERATURE                     The temperature.
!
! . Arrays:
!
!   ATMVEL                          The current atom velocities.
!   PROJECT                         The projection vectors.
!
!   TOTANG                          The total angular momentum.
!   TOTMOM                          The total momentum.
!
! . Functions:
!
!   VELOCITY_ANGULAR_MOMENTUM       Calculate the total system angular momentum.
!   VELOCITY_MOMENTUM               Calculate the total system momentum.
!
! . Subroutines:
!
!   VELOCITY_ASSIGN                 Assign the atom velocities.
!   VELOCITY_DEGREES_OF_FREEDOM     Set the number of degrees of freedom for the
!                                   system and calculate the projection vectors.
!   VELOCITY_INITIALIZE             Initialize the velocity data structure.
!   VELOCITY_PROJECT                Project out the translational vectors from
!                                   an atom data set.
!   VELOCITY_READ                   Read in the velocities.
!   VELOCITY_SCALE                  Scale velocities to a specific temperature.
!   VELOCITY_TEMPERATURE            Calculate the KE and the temperature.
!   VELOCITY_WRITE                  Write out the velocities.
!
!===============================================================================
MODULE VELOCITY

! . Data structure declarations.
USE CONSTANTS,      ONLY : AMUA2PS2_TO_K, AMUA2PS2_TO_KJMOL, AMU_TO_KG, KBOLTZ, MS_TO_APS
USE DEFINITIONS,    ONLY : DP
USE PRINTING,       ONLY : PRINT_ERROR, PRINT_LINE, PRINT_PARAGRAPH, PRINT_SUMMARY_ELEMENT, &
                           PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START, PRINT_SUMMARY_STOP
USE RANDOM_NUMBERS, ONLY : RANDOM_GAUSS

USE ATOMS,          ONLY : ATMCRD, ATMFIX, ATMMAS, NATOMS, NFIXED, NFREE
USE CONSTRAINT,     ONLY : QTETHER
USE LINEAR_ALGEBRA, ONLY : CROSS_PRODUCT, SCHMIDT_ORTHOGONALIZE
USE XYZ_IO,         ONLY : XYZ_READ, XYZ_WRITE

IMPLICIT NONE
PRIVATE
PUBLIC :: ATMVEL, EKE, TEMPERATURE, &
          VELOCITY_ASSIGN, VELOCITY_PROJECT, VELOCITY_READ, VELOCITY_SCALE, VELOCITY_TEMPERATURE, VELOCITY_WRITE
#ifndef PGPC
SAVE
#endif

! . Module scalars.
INTEGER            :: NDEGF = 0,      NPROJECT    = 0
REAL ( KIND = DP ) :: EKE   = 0.0_DP, TEMPERATURE = 0.0_DP

! . Module arrays.
REAL ( KIND = DP ),              DIMENSION(1:3)   :: TOTANG = 0.0_DP, TOTMOM = 0.0_DP
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:)   :: ATMVEL
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:,:) :: PROJECT

!==============================================================================
CONTAINS
!==============================================================================

   !-------------------------------------
   FUNCTION VELOCITY_ANGULAR_MOMENTUM ( )
   !-------------------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3) :: VELOCITY_ANGULAR_MOMENTUM

   ! . Local scalars.
   INTEGER :: I

   ! . Initialize the function value.
   VELOCITY_ANGULAR_MOMENTUM = 0.0_DP

   ! . Loop over the atoms.
   DO I = 1,NATOMS
      IF ( .NOT. ATMFIX(I) ) THEN
         VELOCITY_ANGULAR_MOMENTUM = VELOCITY_ANGULAR_MOMENTUM + ATMMAS(I) * CROSS_PRODUCT ( ATMCRD(1:3,I), ATMVEL(1:3,I) )
      END IF
   END DO

   END FUNCTION VELOCITY_ANGULAR_MOMENTUM

   !-----------------------------
   FUNCTION VELOCITY_MOMENTUM ( )
   !-----------------------------

   ! . Function declarations.
   REAL ( KIND = DP ), DIMENSION(1:3) :: VELOCITY_MOMENTUM

   ! . Local scalars.
   INTEGER :: I

   ! . Initialize the function value.
   VELOCITY_MOMENTUM = 0.0_DP

   ! . Loop over the atoms.
   DO I = 1,NATOMS
      IF ( .NOT. ATMFIX(I) ) VELOCITY_MOMENTUM = VELOCITY_MOMENTUM + ATMMAS(I) * ATMVEL(1:3,I)
   END DO

   END FUNCTION VELOCITY_MOMENTUM

   !----------------------------------------------------------------
   SUBROUTINE VELOCITY_ASSIGN ( TNEEDED, REMOVE_TRANSLATION, PRINT )
   !----------------------------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: TNEEDED

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT, REMOVE_TRANSLATION

   ! . Local scalars.
   INTEGER            :: I, J
   LOGICAL            :: QPRINT
   REAL ( KIND = DP ) :: KT, SD

   ! . Initialize the velocity data structure.
   CALL VELOCITY_INITIALIZE

   ! . Check the number of free atoms in the system.
   IF ( NFREE <= 0 ) RETURN

   ! . Determine the number of degrees of freedom in the system.
   CALL VELOCITY_DEGREES_OF_FREEDOM ( REMOVE_TRANSLATION, PRINT )

   ! . Allocate the velocity array.
   ALLOCATE ( ATMVEL(1:3,1:NATOMS) ) ; ATMVEL = 0.0_DP

   ! . Assign the velocities for a Maxwell-Boltzmann distribution at TNEEDED Kelvin.
   KT = KBOLTZ * TNEEDED
   DO I = 1,NATOMS
      IF ( .NOT. ATMFIX(I) ) THEN
         SD = SQRT ( KT / ( AMU_TO_KG * ATMMAS(I) ) )
         DO J = 1,3
            ATMVEL(J,I) = RANDOM_GAUSS ( 0.0_DP, SD )
         END DO
      END IF
   END DO

   ! . Convert the velocities from m s^-1 to A ps^-1.
   ATMVEL = ATMVEL * MS_TO_APS

   ! . Remove the translational contributions from the velocities.
   CALL VELOCITY_PROJECT ( ATMVEL )

   ! . Calculate the actual temperature.
   CALL VELOCITY_TEMPERATURE

   ! . Scale all the velocities uniformly to get the desired temperature.
   CALL VELOCITY_SCALE ( TNEEDED )

   ! . Get the print flag.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Write out some information.
   IF ( QPRINT ) THEN
      WRITE ( PRINT_LINE, "(A,F10.3,A)" ) "Velocities assigned to the system at a temperature of ", TNEEDED, " K."
      CALL PRINT_PARAGRAPH
   END IF

   END SUBROUTINE VELOCITY_ASSIGN

   !-------------------------------------------------------------------
   SUBROUTINE VELOCITY_DEGREES_OF_FREEDOM ( REMOVE_TRANSLATION, PRINT )
   !-------------------------------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: PRINT, REMOVE_TRANSLATION

   ! . Local scalars.
   INTEGER :: I, II, NVECTORS
   LOGICAL :: QPRINT, QREMOVET

   ! . Local arrays.
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: VECTORS

   ! . Initialize the number of degrees of freedom.
   NDEGF = 3 * NFREE

   ! . Determine the number of translational vectors.
   ! . There are fixed atoms or harmonic position constraints.
   IF ( ( NFIXED > 0 ) .OR. QTETHER ) THEN

      ! . There are no such modes.
      QREMOVET = .FALSE.

   ! . There are no fixed atoms.
   ELSE

      ! . Remove the translations if there is more than one atom (it does not matter if there are PBCs).
      QREMOVET = ( NATOMS > 1 )

   END IF

   ! . Check the passed remove translation option.
   IF ( PRESENT ( REMOVE_TRANSLATION ) ) THEN
      IF ( REMOVE_TRANSLATION ) THEN
         IF ( .NOT. QREMOVET ) THEN
	    CALL PRINT_ERROR ( "VELOCITY_DEGREES_OF_FREEDOM", "It is not possible to remove translations in this case." )
	 END IF
      ELSE
         QREMOVET = .FALSE.
      END IF
   END IF

   ! . Determine the translational modes.
   IF ( QREMOVET ) THEN

      ! . Initialize the number of projection vectors.
      NVECTORS = 3

      ! . Allocate and initialize the VECTORS array.
      ALLOCATE ( VECTORS(1:3*NATOMS,1:3) ) ; VECTORS = 0.0_DP

      ! . Loop over the atoms.
      DO I = 1,NATOMS

         ! . Calculate the array index.
         II = 3 * ( I - 1 )

         ! . The translations.
         VECTORS(II+1,1) = ATMMAS(I)
         VECTORS(II+2,2) = ATMMAS(I)
         VECTORS(II+3,3) = ATMMAS(I)

      END DO

      ! . Orthogonalize the vectors.
      CALL SCHMIDT_ORTHOGONALIZE ( VECTORS(1:3*NATOMS,1:NVECTORS), NPROJECT )

      ! . There are independent vectors.
      IF ( NPROJECT > 0 ) THEN

         ! . Allocate PROJECT.
         IF ( ALLOCATED ( PROJECT ) ) DEALLOCATE ( PROJECT )
	 ALLOCATE ( PROJECT(1:3,1:NATOMS,1:NPROJECT) )

         ! . Copy VECTORS to NPROJECT.
	 PROJECT = RESHAPE ( VECTORS(1:3*NATOMS,1:NPROJECT), (/ 3, NATOMS, NPROJECT /) )

      END IF

      ! . Deallocate temporary space.
      DEALLOCATE ( VECTORS )

   ! . There are no projection vectors.
   ELSE
      NPROJECT = 0
   END IF

   ! . Reset the number of degrees of freedom.
   NDEGF = NDEGF - NPROJECT

   ! . Check the number of degrees of freedom.
   IF ( ( NDEGF <= 0 ) .OR. ( NDEGF > 3 * NFREE ) ) THEN
      CALL PRINT_ERROR ( "VELOCITY_DEGREES_OF_FREEDOM", "The number of degrees of freedom is invalid." )
   END IF

   ! . Get the print flag.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Do some printing.
   IF ( QPRINT ) THEN
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#005500", VARIABLEWIDTH = 8 )
      CALL PRINT_SUMMARY_START ( "Degrees of Freedom Data" )
      WRITE ( PRINT_LINE, "(I8)" ) NDEGF    ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Degrees of Freedom" )
      WRITE ( PRINT_LINE, "(I8)" ) NPROJECT ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Projection Vectors" )
      CALL PRINT_SUMMARY_STOP
   END IF

   END SUBROUTINE VELOCITY_DEGREES_OF_FREEDOM

   !-----------------------------
   SUBROUTINE VELOCITY_INITIALIZE
   !-----------------------------

   ! . Initialize some variables.
   NDEGF       = 0
   NPROJECT    = 0
   EKE         = 0.0_DP
   TEMPERATURE = 0.0_DP
   TOTANG      = 0.0_DP
   TOTMOM      = 0.0_DP

   ! . Deallocate the velocity array.
   IF ( ALLOCATED ( ATMVEL  ) ) DEALLOCATE ( ATMVEL  )
   IF ( ALLOCATED ( PROJECT ) ) DEALLOCATE ( PROJECT )

   END SUBROUTINE VELOCITY_INITIALIZE

   !-------------------------------------
   SUBROUTINE VELOCITY_PROJECT ( ATMDAT )
   !-------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(1:3,1:NATOMS), INTENT(INOUT) :: ATMDAT

   ! . Local scalars.
   INTEGER            :: I
   REAL ( KIND = DP ) :: DOTFAC

   ! . There are projection vectors.
   IF ( NPROJECT > 0 ) THEN

      ! . Project out PROJECT from ATMDAT.
      DO I = 1,NPROJECT
         DOTFAC = SUM ( ATMDAT * PROJECT(:,:,I) )
         ATMDAT = ATMDAT - DOTFAC * PROJECT(:,:,I)
      END DO

   END IF

   END SUBROUTINE VELOCITY_PROJECT

   !----------------------------------------------------
   SUBROUTINE VELOCITY_READ ( FILE, REMOVE_TRANSLATION )
   !----------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Optional scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: REMOVE_TRANSLATION

   ! . Reinitialize the velocity data structure.
   CALL VELOCITY_INITIALIZE

   ! . Check the number of free atoms.
   IF ( NFREE <= 0 ) RETURN

   ! . Determine the number of degrees of freedom in the system.
   CALL VELOCITY_DEGREES_OF_FREEDOM ( REMOVE_TRANSLATION )

   ! . Allocate the velocity array.
   ALLOCATE ( ATMVEL(1:3,1:NATOMS) ) ; ATMVEL = 0.0_DP

   ! . Read in the velocities.
   CALL XYZ_READ ( FILE, ATMVEL )

   ! . Make sure they satisfy the constraints.
   CALL VELOCITY_PROJECT ( ATMVEL )

   ! . Calculate the temperature.
   CALL VELOCITY_TEMPERATURE

   END SUBROUTINE VELOCITY_READ

   !------------------------------------
   SUBROUTINE VELOCITY_SCALE ( TNEEDED )
   !------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: TNEEDED

   ! . Local scalars.
   REAL ( KIND = DP ) :: SCALE

   ! . Calculate the scale factor.
   SCALE = TNEEDED / TEMPERATURE

   ! . The scaling factor is invalid.
   IF ( SCALE <= 0.0_DP ) THEN

      ! . Set the velocities to zero.
      ATMVEL = 0.0_DP

      ! . Reset the kinetic energy, temperature, momentum and angular momentum.
      EKE         = 0.0_DP
      TEMPERATURE = 0.0_DP
      TOTANG      = 0.0_DP
      TOTMOM      = 0.0_DP

   ! . The scaling factor is OK.
   ELSE

      ! . Scale all the velocities uniformly to get the desired temperature.
      ATMVEL = SQRT ( SCALE ) * ATMVEL

      ! . Reset the actual kinetic energy, temperature, momentum and angular momentum.
      EKE         = SCALE * EKE
      TEMPERATURE = TNEEDED
      TOTANG      = SQRT ( SCALE ) * TOTANG
      TOTMOM      = SQRT ( SCALE ) * TOTMOM

   END IF

   END SUBROUTINE VELOCITY_SCALE

   !------------------------------
   SUBROUTINE VELOCITY_TEMPERATURE
   !------------------------------

   ! . Local scalars.
   INTEGER            :: I
   REAL ( KIND = DP ) :: KE

   ! . Initialization.
   EKE         = 0.0_DP
   TEMPERATURE = 0.0_DP
   TOTANG      = 0.0_DP
   TOTMOM      = 0.0_DP

   ! . Check that the number of atoms and degrees of freedom are OK and that there are velocities.
   IF ( ( NATOMS <= 0 ) .OR. ( NDEGF <= 0 ) .OR. ( .NOT. ALLOCATED ( ATMVEL ) ) ) RETURN

   ! . Loop over the atoms.
   KE = 0.0_DP
   DO I = 1,NATOMS
      IF ( .NOT. ATMFIX(I) ) KE = KE + ATMMAS(I) * DOT_PRODUCT ( ATMVEL(1:3,I), ATMVEL(1:3,I) )
   END DO

   ! . Calculate the kinetic energy and temperature.
   EKE         = 0.5_DP * AMUA2PS2_TO_KJMOL * KE
   TEMPERATURE =          AMUA2PS2_TO_K     * KE / REAL ( NDEGF, DP )

   ! . Calculate the momentum and angular momentum.
   ! . These are not really needed here and they produce a problem with the PGF90 compiler in any case.
!   TOTMOM = VELOCITY_MOMENTUM ( )
!   TOTANG = VELOCITY_ANGULAR_MOMENTUM ( )

   END SUBROUTINE VELOCITY_TEMPERATURE

   !---------------------------------
   SUBROUTINE VELOCITY_WRITE ( FILE )
   !---------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Check the coordinate and velocity arrays.
   IF ( .NOT. ALLOCATED ( ATMVEL ) .OR. ( SIZE ( ATMCRD ) /= SIZE ( ATMVEL ) ) ) RETURN

   ! . Write out the velocities.
   CALL XYZ_WRITE ( FILE, ATMVEL )

   END SUBROUTINE VELOCITY_WRITE

END MODULE VELOCITY
