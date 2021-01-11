!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                      INTERFACE TO GAUSSIAN 09                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutines
!  -----------
!   MY_SELE                       Returns NOFIX atoms (required by hack_chrg.F90)
!   CHRG_SERVER                   Call Gaussian 09 to perform a calculation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine my_sele(sele)

    use initialization, only : nofix_sele

    logical, intent(inout)  :: sele( SIZE(nofix_sele,1) )

    sele(:) = nofix_sele(:)

end subroutine

subroutine chrg_server(code, natm, naqm, atmn, cord, ener, qfit, grad, hess)

  use initialization

  implicit none

  integer, intent( in )                         :: code, natm, naqm
  real*8, dimension(1:natm), intent( in )       :: atmn
  real*8, dimension(1:3*natm), intent( in )     :: cord
  real*8, intent( inout )                       :: ener
  real*8, dimension(1:naqm), intent( inout )    :: qfit
  real*8, dimension(1:*), intent( inout )       :: grad
  real*8, dimension(1:*), intent( inout )       :: hess

  integer                                       :: i, j, k, l, n3, nh
  real*8                                        :: smm
  character( len=256 )                          :: str

  real*8, parameter :: bohr = 0.529177249d0
  character( len=2 ), dimension(1:109), parameter :: smb = (/ &
      "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", "Na", "Mg", &
      "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", &
      "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
      "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", &
      "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", &
      "Tm", "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", &
      "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", &
      "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt" /)


  open( unit=999, file = "calc.com", action = "write", form = "formatted" )
  write( 999, "(a/a/a/a)" ) "%Chk=calc.chk", &
    "%NProcShared=" // trim(cores), &
    "%Mem=" // trim(memory), &
    "#p " // trim(dft_func) // "/" // trim(dft_basis) // " charge guess=(read,tcheck) scf=(direct,conver=6) nosymm"
  if( code == 1 ) then
      write( 999, "(a/)" ) "force pop=chelpg fchk"
  else if( code == 2 ) then
      write( 999, "(a/)" ) "freq=noraman pop=chelpg fchk"
  else
      write( 999, "(a/)" ) "pop=chelpg fchk"
  end if
  write( 999, "(a//2i4)" ) "- slave -", qm_charge, qm_multi
  do i = 1, naqm
    j = ( i - 1 ) * 3
    write( 999, "(a8,3f20.10)" ) smb(abs(int(atmn(i)))), cord(j+1), cord(j+2), cord(j+3)
  end do
  write( 999, "(a)" ) ""
  smm = 0.0d0
  if( naqm == natm ) then
    write( 999, "(f28.10,2f20.10,f8.3)" ) 999.0d0, 999.0d0, 999.0d0, 0.0d0
  else
    do i = naqm + 1, natm
      j = ( i - 1 ) * 3
      write( 999, "(f28.10,2f20.10,f8.3)" ) cord(j+1), cord(j+2), cord(j+3), atmn(i)
      do k = i + 1, natm
        l = ( k - 1 ) * 3
        smm = smm + atmn(i) * atmn(k) / dsqrt( sum( ( cord(j+1:j+3) - cord(l+1:l+3) ) ** 2 ) ) * bohr
      end do
    end do
  end if
  write( 999, "(//)" )
  close( 999 )

  ! write(*,"(a,f20.10)") "self_charges:", smm

  CALL system( "g09 < calc.com > calc.log" )

  ener = 0.0d0
  qfit = 0.0d0
  open( unit = 999, file = "Test.FChk", action = "read", form = "formatted" )
  read( 999, "(a)", end = 999 ) str
  do while( .true. )
    if( str(1:12) == "Total Energy" ) read( str(50:71), "(f22.15)", end = 999 ) ener

    if( str(1:11) == "ESP Charges" ) then
      k = 1
      do i = 1, int( naqm / 5 )
        read( 999, "(a)", end = 999 ) str
        do j = 1, 5
          read( str((j-1)*16+1:j*16), "(f16.8)" ) qfit(k)
          k = k + 1
        end do
      end do
      i = mod( naqm, 5 )
      if( i > 0 ) then
        read( 999, "(a)", end = 999 ) str
        do j = 1, i
          read( str((j-1)*16+1:j*16), "(f16.8)" ) qfit(k)
          k = k + 1
        end do
      end if
    end if

    if( str(1:18) == "Cartesian Gradient" .and. ( code == 1 .or. code == 2 ) ) then
      n3 = 3 * naqm
      grad(1:n3) = 0.0d0
      k = 1
      do i = 1, int( n3 / 5 )
        read( 999, "(a)", end = 999 ) str
        do j = 1, 5
          read( str((j-1)*16+1:j*16), "(f16.8)" ) grad(k)
          k = k + 1
        end do
      end do
      i = mod( n3, 5 )
      if( i > 0 ) then
        read( 999, "(a)", end = 999 ) str
        do j = 1, i
          read( str((j-1)*16+1:j*16), "(f16.8)" ) grad(k)
          k = k + 1
        end do
      end if
    end if

    if( str(1:25) == "Cartesian Force Constants" .and. code == 2 ) then
      nh = 3 * naqm * ( 3 * naqm + 1 ) / 2
      hess(1:nh) = 0.0d0
      k = 1
      do i = 1, int( nh / 5 )
        read( 999, "(a)", end = 999 ) str
        do j = 1, 5
          read( str((j-1)*16+1:j*16), "(f16.8)" ) hess(k)
          k = k + 1
        end do
      end do
      i = mod( nh, 5 )
      if( i > 0 ) then
        read( 999, "(a)", end = 999 ) str
        do j = 1, i
          read( str((j-1)*16+1:j*16), "(f16.8)" ) hess(k)
          k = k + 1
        end do
      end if
    end if

    read( 999, "(a)", end = 999 ) str
  end do
  999  continue
  close( 999 )

  ener = ener - smm

  CALL system( "/bin/rm -f Test.FChk" )

end subroutine
