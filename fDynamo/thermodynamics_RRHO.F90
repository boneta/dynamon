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
!                           The Thermodynamics Module
!===============================================================================
!
! . Subroutines:
!
!   THERMODYNAMICS                Calculate some thermodynamic properties.
!
! . Notes:
!
!   Only free atoms are counted in the thermodynamic analysis.
!
!===============================================================================
MODULE THERMODYNAMICS_RRHO

! . Module declarations.
USE CONSTANTS
USE DEFINITIONS, ONLY : DP
USE PRINTING,    ONLY : PRINT_LINE, PRINT_SUMMARY_ELEMENT, PRINT_SUMMARY_OPTIONS, PRINT_SUMMARY_START,   &
                        PRINT_SUMMARY_STOP, PRINT_TABLE_ELEMENT, PRINT_TABLE_OPTIONS, PRINT_TABLE_START, &
			PRINT_TABLE_STOP

USE ATOMS,       ONLY : ATMCRD, ATMFIX, ATMMAS, NATOMS, NFIXED, NFREE
USE NORMAL_MODE, ONLY : FREQUENCIES, NMODES
USE SYMMETRY,    ONLY : QBOX
USE TRANSFORMATION


IMPLICIT NONE
PRIVATE
PUBLIC :: THERMODYNAMICS

! . Type definitions.
TYPE TDICS
   REAL ( KIND = DP ) :: LNZ, A, CP, CV, G, H, S, U
END TYPE TDICS

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------------------------------------------------------
   SUBROUTINE THERMODYNAMICS ( TEMPERATURES, RESULTS, SYMMETRY_NUMBER, PRESSURE, PRINT )
   !------------------------------------------------------------------------------------

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(IN) :: TEMPERATURES

   ! . Optional scalar arguments.
   INTEGER,            INTENT(IN), OPTIONAL :: SYMMETRY_NUMBER
   LOGICAL,            INTENT(IN), OPTIONAL :: PRINT
   REAL ( KIND = DP ), INTENT(IN), OPTIONAL :: PRESSURE

   ! . Optional array arguments.
   REAL ( KIND = DP ), DIMENSION(:,:), INTENT(OUT), OPTIONAL :: RESULTS

   ! . Local parameters.
   REAL ( KIND = DP ) :: UNDERFLOW = 75.0_DP

   ! . Local scalars.
   INTEGER            :: I, ITEMP, ITR, MODE, NIMAG, NTEMP, NTRROT
   LOGICAL            :: QLINEAR, QPRINT, QSAVE
   REAL ( KIND = DP ) :: EXPM, EXPP, FACTOR, HVKT, P, R, RT, SMALL, SYMNUM, T, TOTMAS, VOLUME, ZPE

   ! . Local TDics scalars.
   REAL ( KIND = DP ) :: A, CP, CV, G, H, LNZ, S, U, Z

   ! . Local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3)      :: MOMENTS
   REAL ( KIND = DP ), DIMENSION(1:NATOMS) :: MASSES
   REAL ( KIND = DP ), DIMENSION(1:NMODES) :: OMEGA

   ! . Local types.
   TYPE(TDICS) :: ELECTRONIC, ROTATION, TRANSLATION, VIBRATION

   !---------------------------------------------------------------------------
   ! . Initialization.
   !---------------------------------------------------------------------------
   ! . Check the number of free atoms and the number of modes.
   IF ( ( NFREE <=0 ) .OR. ( NMODES <= 0 ) .OR. ( NMODES /= 3 * NFREE ) ) RETURN

   ! . Get the number of temperatures.
   NTEMP = SIZE ( TEMPERATURES )
   IF ( NTEMP <= 0 ) RETURN

   ! . Check the RESULTS argument.
   IF ( PRESENT ( RESULTS ) ) THEN
      IF ( ( SIZE ( RESULTS, 1 ) /= 7 ) .OR. ( SIZE ( RESULTS, 2 ) /= NTEMP ) ) RETURN
      QSAVE = .TRUE.
   ELSE
      QSAVE = .FALSE.
   END IF

   ! . Get the symmetry number.
   SYMNUM = 1.0_DP
   IF ( PRESENT ( SYMMETRY_NUMBER ) ) SYMNUM = REAL ( SYMMETRY_NUMBER, DP )

   ! . Calculate the pressure.
   P = ATM_TO_PASCALS
   IF ( PRESENT ( PRESSURE ) ) P = P * PRESSURE

   ! . Calculate R (kJ mol^-1 K^-1).
   R = 1.0E-3_DP * KBOLTZ * NAVOGADRO

   ! . Get the masses array.
   MASSES = ATMMAS
   WHERE ( ATMFIX ) MASSES = 0.0_DP

   ! . Calculate the total mass of the free atoms in the system in kg.
   TOTMAS = AMU_TO_KG * SUM ( MASSES )

   ! . Do some processing for the rotational and vibrational contributions.
   IF ( NATOMS > 1 ) THEN

      ! . Get the moments of inertia of the system (in amu Angstroms**2).
      CALL MOMENTS_OF_INERTIA_AT_CENTER ( ATMCRD, MOMENTS, WEIGHTS = MASSES )

      ! . Convert the moments to SI units.
      MOMENTS = AMU_TO_KG * 1.0E-20_DP * MOMENTS

      ! . Determine if the molecule is linear.
      QLINEAR = COUNT ( MOMENTS == 0.0_DP ) >= 2

      ! . Save the frequencies.
      OMEGA = FREQUENCIES

      ! . Determine the number of rotational and translational modes.
      ! . All atoms are free.
      IF ( NFIXED == 0 ) THEN

         ! . Single atom.
	 IF ( NATOMS == 1 ) THEN
	    NTRROT = 3
	 ! . Multiple atoms.
	 ELSE

            ! . Linear systems.
	    IF ( QLINEAR ) THEN
	       NTRROT = 5
	    ! . Periodic systems.
	    ELSE IF ( QBOX ) THEN
	       NTRROT = 3
	    ! . General case.
	    ELSE
               NTRROT = 6
	    END IF

	 END IF

      ! . There are fixed atoms (this result will not be strictly true for small numbers of fixed atoms).
      ELSE
         NTRROT = 0
      END IF

      ! . Zero out the frequencies of the rotational and translational modes.
      DO ITR = 1,NTRROT
         SMALL = HUGE ( 0.0_DP )
         DO I = 1,NMODES
            IF ( ( ABS ( OMEGA(I) ) < SMALL ) .AND. ( OMEGA(I) /= 0.0_DP ) ) THEN
               SMALL = ABS ( OMEGA(I) )
               MODE  = I
            END IF
         END DO
         OMEGA(MODE) = 0.0_DP
      END DO

      ! . Count the number of imaginary modes.
      NIMAG = COUNT ( OMEGA < 0.0_DP )

      ! . Convert the frequencies to J.
      OMEGA = 1.0E+12_DP * PLANCK * TO_HZ * OMEGA

      ! . Calculate the zero-point energy in kJ mol^-1.
      ZPE = 0.5_DP * 1.0E-3_DP * NAVOGADRO * SUM ( OMEGA(NIMAG+NTRROT+1:NMODES) )

   ! . The system is monatomic.
   ELSE
      ZPE = 0.0_DP
   END IF

   ! . Set the QPRINT flag.
   QPRINT = .TRUE.
   IF ( PRESENT ( PRINT ) ) QPRINT = PRINT

   ! . Do some printing.
   IF ( QPRINT ) THEN

      ! . Print out some basic data.
      CALL PRINT_SUMMARY_OPTIONS ( HEADER_COLOR = "#AA6666", VARIABLEWIDTH = 16 )
      CALL PRINT_SUMMARY_START ( "Ideal Gas/Rigid Rotor/Harmonic Oscillator Thermodynamics" )
      WRITE ( PRINT_LINE, "(I16)"   ) NTEMP      ; CALL PRINT_SUMMARY_ELEMENT ( "Number of Temps."   )
      WRITE ( PRINT_LINE, "(F16.0)" ) SYMNUM     ; CALL PRINT_SUMMARY_ELEMENT ( "Symmetry Number"    )
      WRITE ( PRINT_LINE, "(G16.6)" ) P          ; CALL PRINT_SUMMARY_ELEMENT ( "Pressure (Pa)"      )
      WRITE ( PRINT_LINE, "(G16.6)" ) TOTMAS     ; CALL PRINT_SUMMARY_ELEMENT ( "Total Mass (kg)"    )
      WRITE ( PRINT_LINE, "(G16.6)" ) MOMENTS(1) ; CALL PRINT_SUMMARY_ELEMENT ( "Moment A (kg m^-2)" )
      WRITE ( PRINT_LINE, "(G16.6)" ) MOMENTS(2) ; CALL PRINT_SUMMARY_ELEMENT ( "Moment B"           )
      WRITE ( PRINT_LINE, "(G16.6)" ) MOMENTS(3) ; CALL PRINT_SUMMARY_ELEMENT ( "Moment C"           )
      WRITE ( PRINT_LINE, "(G16.6)" ) ZPE        ; CALL PRINT_SUMMARY_ELEMENT ( "ZPE (kJ mol^-1)"    )
      WRITE ( PRINT_LINE, "(I16)"   ) NIMAG      ; CALL PRINT_SUMMARY_ELEMENT ( "Imaginary Modes"    )
      WRITE ( PRINT_LINE, "(I16)"   ) NTRROT     ; CALL PRINT_SUMMARY_ELEMENT ( "Rot./Tr. Modes"     )
      CALL PRINT_SUMMARY_STOP

      ! . Print the table header.
      CALL PRINT_TABLE_OPTIONS ( COLUMNS = 10, HEADER_COLOR = "#AA6666", PAGEWIDTH = 120 )
      CALL PRINT_TABLE_START
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Temperature Results (K, m^3, kJ mol^-1)", COLSPAN = 10, HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Temperature", HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Volume",      HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Ln z",        HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "S",           HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "U",           HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "H",           HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "A",           HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "G",           HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Cv",          HEADER = .TRUE. )
      CALL PRINT_TABLE_ELEMENT ( TEXT = "Cp",          HEADER = .TRUE. )
   END IF


   !---------------------------------------------------------------------------
   ! . Loop over the temperatures.
   !---------------------------------------------------------------------------
   DO ITEMP = 1,NTEMP

      ! . Get the temperature.
      T = TEMPERATURES(ITEMP)

      ! . Calculate RT (kJ mol^-1).
      RT = R * T

      ! . Calculate the volume.
      VOLUME = KBOLTZ * T / P

      !------------------------------------------------------------------------
      ! . Electronic contributions.
      !------------------------------------------------------------------------
      ! . Assume there is one state with a multiplicity of one.
      ELECTRONIC = TDICS( 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP )

      !------------------------------------------------------------------------
      ! . Rotational contributions.
      !------------------------------------------------------------------------
      ! . Monatomic system.
      IF ( NFREE == 1 ) THEN

         ! . There are no rotational states.
         ROTATION = TDICS( 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP )

      ! . Polyatomic system.
      ELSE

         ! . Calculate some factors.
         FACTOR = ( 8.0_DP * PI * PI * KBOLTZ * T ) / ( PLANCK**2 )

         ! . Linear molecule.
         IF ( QLINEAR ) THEN
            Z            = FACTOR * MOMENTS(3) / SYMNUM
            ROTATION%LNZ = LOG ( Z )
            ROTATION%CV  = R
            ROTATION%S   = R * ( ROTATION%LNZ + 1.0_DP )
            ROTATION%U   = RT
         ! . Non-linear molecule.
         ELSE
            Z            = SQRT ( PI * FACTOR**3 * PRODUCT ( MOMENTS ) ) / SYMNUM
            ROTATION%LNZ = LOG ( Z )
            ROTATION%CV  = 1.5_DP * R
            ROTATION%S   = R * ( ROTATION%LNZ + 1.5_DP )
            ROTATION%U   = 1.5_DP * RT
         END IF

         ! . Calculate the common terms.
         ROTATION%A   = - RT * ROTATION%LNZ
         ROTATION%CP  = ROTATION%CV
         ROTATION%H   = ROTATION%U
         ROTATION%G   = ROTATION%H - T * ROTATION%S

      END IF

      !------------------------------------------------------------------------
      ! . Translational contributions.
      !------------------------------------------------------------------------
      ! . Calculate the terms.
      Z               = ( SQRT ( 2.0_DP * PI * TOTMAS * KBOLTZ * T ) / PLANCK )**3 * VOLUME
      TRANSLATION%LNZ = LOG ( Z )
      TRANSLATION%A   = - RT * TRANSLATION%LNZ
      TRANSLATION%CP  = 2.5_DP * R
      TRANSLATION%CV  = 1.5_DP * R
      TRANSLATION%H   = 2.5_DP * RT
      TRANSLATION%S   = R * ( TRANSLATION%LNZ + 1.5_DP )
      TRANSLATION%U   = 1.5_DP * RT
      TRANSLATION%G   = TRANSLATION%H - T * TRANSLATION%S

      !------------------------------------------------------------------------
      ! . Vibrational contributions.
      !------------------------------------------------------------------------
      ! . Monatomic system with no environment.
      IF ( NATOMS == 1 ) THEN

         ! . There are no vibrational states.
         VIBRATION = TDICS ( 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP )

      ! . Polyatomic system.
      ELSE

         ! . Loop over the frequencies.
         CV  = 0.0_DP
         LNZ = 0.0_DP
         U   = 0.0_DP
         DO I = (NIMAG+NTRROT+1),NMODES 
            HVKT = OMEGA(I) / ( KBOLTZ * T )
            LNZ  = LNZ - 0.5_DP * HVKT
            IF ( HVKT <= UNDERFLOW ) THEN
               EXPM = EXP ( - HVKT )
               EXPP = 1.0_DP - EXPM
               LNZ  = LNZ - LOG ( EXPP )
               EXPP = EXPM / EXPP
               CV   =  CV + R * HVKT * HVKT * EXPP * ( 1.0_DP + EXPP )
               U    =  U + RT * HVKT * EXPP
            END IF
         END DO 
    
         ! . Assign the terms.
         VIBRATION%LNZ = LNZ
         VIBRATION%A   = - RT * VIBRATION%LNZ
         VIBRATION%CP  = CV
         VIBRATION%CV  = CV
         VIBRATION%U   = U + ZPE
         VIBRATION%S   = R * VIBRATION%LNZ + VIBRATION%U / T
         VIBRATION%G   = VIBRATION%U - T * VIBRATION%S
         VIBRATION%H   = VIBRATION%U

      END IF

      !------------------------------------------------------------------------
      ! . Find the totals.
      !------------------------------------------------------------------------
      ! . Find the local totals.
      LNZ = ELECTRONIC%LNZ + ROTATION%LNZ + TRANSLATION%LNZ + VIBRATION%LNZ
      A   = ELECTRONIC%A   + ROTATION%A   + TRANSLATION%A   + VIBRATION%A
      CP  = ELECTRONIC%CP  + ROTATION%CP  + TRANSLATION%CP  + VIBRATION%CP
      CV  = ELECTRONIC%CV  + ROTATION%CV  + TRANSLATION%CV  + VIBRATION%CV
      G   = ELECTRONIC%G   + ROTATION%G   + TRANSLATION%G   + VIBRATION%G
      H   = ELECTRONIC%H   + ROTATION%H   + TRANSLATION%H   + VIBRATION%H
      S   = ELECTRONIC%S   + ROTATION%S   + TRANSLATION%S   + VIBRATION%S
      U   = ELECTRONIC%U   + ROTATION%U   + TRANSLATION%U   + VIBRATION%U

      ! . Add in some extra terms pertaining to the ideal gas.
      FACTOR = LOG ( NAVOGADRO ) - 1.0_DP
      A = A + RT * FACTOR
      G = G - RT * FACTOR
      S = S - R  * FACTOR


      !------------------------------------------------------------------------
      ! . Print out the results.
      !------------------------------------------------------------------------
      ! . Do the printing.
      IF ( QPRINT ) THEN
         WRITE ( PRINT_LINE, "(G12.4)" ) T                  ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) NAVOGADRO * VOLUME ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) LNZ                ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) S                  ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) U                  ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) H                  ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) A                  ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) G                  ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) CV                 ; CALL PRINT_TABLE_ELEMENT
         WRITE ( PRINT_LINE, "(G12.4)" ) CP                 ; CALL PRINT_TABLE_ELEMENT
      END IF

      ! . Save the results if required.
      IF ( QSAVE ) THEN
         RESULTS(1,ITEMP) = A
         RESULTS(2,ITEMP) = CP
         RESULTS(3,ITEMP) = CV
         RESULTS(4,ITEMP) = G
         RESULTS(5,ITEMP) = H
         RESULTS(6,ITEMP) = S
         RESULTS(7,ITEMP) = U
      END IF

   END DO

   ! . Write out the terminator.
   IF ( QPRINT ) CALL PRINT_TABLE_STOP   

   END SUBROUTINE THERMODYNAMICS

END MODULE THERMODYNAMICS_RRHO
