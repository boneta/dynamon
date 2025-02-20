!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                          BAKER ALGORITHM                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutines
!  -----------
!   SETUP_OPTIMIZATION            ?
!   CLEANUP_OPTIMIZATION          ?
!   EGH_CALC                      ?
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module panadero

  use initialization,        only : hessfile

  use definitions,           only : dp
  use atoms,                 only : natoms, nfree, nfixed, atmcrd, atmfix, atoms_fix, atmnam
  use mm_terms,              only : atmchg
  use dcd_io,                only : dcd_type, dcd_initialize, dcd_activate_write, dcd_write
  use potential_energy,      only : energy, gradient, hessian, etotal, atmder, atmhes, grms
  use symmetry,              only : boxl
  use normal_mode_utilities, only : raise_rotation_translation
  use lbfgsb_wrapper,        only : optimize_lbfgsb
  use optimize_coordinates,  only : optimize_conjugate_gradient

  implicit none
  private
  save

  type( dcd_type )                   :: dcd_file
  logical, dimension(:), allocatable :: core, envi
  integer                            :: mm_step_number, mm_print_frequency
  real( kind=dp )                    :: mm_gradient_tolerance

  public :: setup_optimization, egh_calc, cleanup_optimization

  contains

    !  SETUP_OPTIMIZATION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine setup_optimization( core_atoms, envi_atoms, &
                                   step_number, gradient_tolerance, &
                                   print_frequency, trajectory )

      logical, dimension(1:natoms), intent( in ) :: core_atoms, envi_atoms
      integer, intent( in )                      :: step_number, print_frequency
      real( kind=dp ), intent( in )              :: gradient_tolerance
      character( len=* ), intent( in )           :: trajectory

      allocate( core(1:natoms), envi(1:natoms) )
      core                  = core_atoms
      envi                  = envi_atoms
      mm_step_number        = step_number
      mm_print_frequency    = print_frequency
      mm_gradient_tolerance = gradient_tolerance
      write( 6, "(80('!'))" )
      write( 6, "(a20,i12)" ) " Core atoms:        ", count( core )
      write( 6, "(a20,i12)" ) " Environment atoms: ", count( envi )
      write( 6, "(a20,i12)" ) " Step Number:       ", mm_step_number
      write( 6, "(a20,i12)" ) " Print Frequency:   ", print_frequency
      write( 6, "(a20,f12.6)" ) " Gradient Tolerance:", mm_gradient_tolerance
      CALL atoms_fix( .not. ( core .or. envi ) )
      CALL dcd_initialize( dcd_file )
      CALL dcd_activate_write( &
        file    = trim( trajectory ), &
        dcd     = dcd_file, &
        type    = "CORD", &
        natoms  = natoms, &
        nfixed  = nfixed, &
        nframes = 10000, &
        qfix    = atmfix )
    end subroutine setup_optimization

    !  CLEANUP_OPTIMIZATION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cleanup_optimization
      deallocate( core, envi )
      CALL flush( dcd_file%unit )
      close( dcd_file%unit )
      write( 6, "(80('!'))" )
    end subroutine cleanup_optimization

    !  EGH_CALC  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine egh_calc( x, e, g, h )
      real( kind=dp ), dimension(:), intent(in)  :: x
      real( kind=dp ), intent(out)               :: e
      real( kind=dp ), dimension(:), intent(out) :: g, h

      integer :: i, j
      real( kind=dp ) :: t

      j = -3
      do i = 1, natoms
        if( core(i) ) then
          j = j + 3
          atmcrd(1:3,i) = x(j+1:j+3)
        end if
      end do

  ! - technically we shoud ask for the charges on the new QM geometry
  !
  !    skip_abinitio = .false.
  !    CALL atoms_fix( .not. acs )
  !    CALL energy
  !
  ! - or we could provide the global GRMS for the system...
  !
      CALL atoms_fix( .not. ( envi .or. core ) )
      CALL gradient
      t = dot_product( atmder(1:3,1), atmder(1:3,1) )
      do i = 2, natoms
        t = max( t, dot_product( atmder(1:3,i), atmder(1:3,i) ) )
      end do
      write( *, "(a,f20.10)" ) ">> Cur GRMS: ", grms
      write( *, "(a,f20.10)" ) ">> Max GRMS: ", dsqrt( t / 3._dp )

      CALL atoms_fix( .not. envi .or. core )
      CALL optimize_lbfgsb( &
        print_frequency    = mm_print_frequency, &
        gradient_tolerance = mm_gradient_tolerance, &
        step_number        = mm_step_number )
      CALL dcd_write( dcd_file, atmcrd, boxl )

      CALL atoms_fix( .not. core )
      CALL hessian( print = .false., file = hessfile )
      e = etotal
      j = -3
      do i = 1, natoms
        if( core(i) ) then
          j = j + 3
          g(j+1:j+3) = atmder(1:3,i)
        end if
      end do
      h = atmhes
      CALL raise_rotation_translation( h, .false. )
    end subroutine egh_calc

end module panadero
