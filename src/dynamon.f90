!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                   _   _           !!
!!                                                  |  \| |          !!
!!               _                                  | . ` |          !!
!!            __| |_   _ _ __   __ _ _ __ ___   ___ |_|\ _|          !!
!!           / _` | | | | '_ \ / _` | '_ ` _ \ / _ \                 !!
!!          | (_| | |_| | | | | (_| | | | | | | (_) |                !!
!!           \__,_|\__, |_| |_|\__,_|_| |_| |_|\___/                 !!
!!                 |___/                                             !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                      A general-purpose script                       !
!                      for common calculations                        !
!                            with fDynamo                             !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 2020, Sergio Boneta
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see https://www.gnu.org/licenses/
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program DYNAMON

  use dynamo

  ! dynamon depedencies
  use initialization
  use common
  use calculation_modes

  CALL dynamon_header
  CALL dynamo_header

  ! read .dynn options file
  CALL read_options
  skip_abinitio = .not. gauss_flg

  ! read .bin and .crd
  CALL mm_system_read(trim(sys_bin))
  CALL coordinates_read(trim(coord))

  ! QM atoms and setup
  allocate(qm_sele(natoms), nofix_sele(natoms))
  CALL qm_selection(qm_sele)
  CALL qm_initialize(gauss_flg)

  ! calculation mode selection
  select case (mode)
    case ('SP','CORR')
      CALL print_mode('SINGLE POINT CALCULATION')
      CALL dynamon_sp
    case ('MINI','SCAN','PES')
      CALL print_mode('MINIMIZATION CALCULATION')
      CALL dynamon_minimization
    case ('LOCATE')
      CALL print_mode('LOCATE CALCULATION')
      CALL dynamon_locate
    case ('IRC')
      CALL print_mode('IRC CALCULATION')
      CALL dynamon_irc
    case ('MD','PMF')
      CALL print_mode('MOLECULAR DYNAMICS CALCULATION')
      CALL dynamon_md
    case ('INTERACTION')
      CALL print_mode('INTERACTION CALCULATION')
      CALL dynamon_interaction
    case ('KIE')
      CALL print_mode('KIE CALCULATION')
      CALL dynamon_kie
    case DEFAULT
      write(*,fmt='(A)') 'ERROR: UNKOWN MODE'
      STOP
  end select

  deallocate(qm_sele, nofix_sele)
  if (gauss_flg) CALL abinitio_exit

  CALL dynamo_footer
  CALL dynamon_footer

end program

include 'nofix_qm.f90'
include 'with_gaussian.f90'
