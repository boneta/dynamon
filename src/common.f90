!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                COMMON GENERAL-PURPOSE SUBROUTINES                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutines
!  -----------
!   PRINT_MODE                    Print the calculation mode in a formatted fashion
!   QM_INITIALIZE                 Set-up the QM calculation type, fix atoms and initialize energy
!   GET_ATOM_NUMBERS              Obtain the atom numbers for constraints references
!   DEFINE_CONSTRAINTS            Set-up constraints
!   OUT_DIST_ENERGY               Write distances corresponding to constraint and energy
!   CALCULATE_GRADIENT            Calculate stable gradient
!   CALCULATE_MINIMIZATION        Perform CG and L-BFGS-B structural optimizations
!
!  Functions
!  ---------
!   DISTANCE_CRD                  Distance for the 'n' constraint from coordinates
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module COMMON


    use dynamo
    use initialization

    contains

    !  PRINT_MODE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine print_mode(mode)

        !--------------------------------------------------------------
        ! Print the calculation mode
        !--------------------------------------------------------------

        character(*), intent(in)           :: mode

        write(*,*)
        write(*,fmt='(A80)') REPEAT('>',80)
        write(*,fmt='(A3,7X,A)') '>>>', trim(mode)
        write(*,fmt='(A80)') REPEAT('>',80)
        write(*,*)

    end subroutine

    !  QM_INITIALIZE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine qm_initialize(abinitio_flg)

        !--------------------------------------------------------------
        ! Set-up QM calculation, fix atoms and initialize energy
        !--------------------------------------------------------------

        logical, intent(in)                :: abinitio_flg

        if (abinitio_flg) then
            CALL abinitio_init
            CALL abinitio_setup(qm_sele)
        else
            CALL mopac_setup( &
                method       = semiemp, &
                charge       = qm_charge, &
                multiplicity = qm_multi, &
                selection    = qm_sele )
            CALL mopac_scf_options(iterations=500000)
        end if

        CALL my_sele( nofix_sele )
        CALL atoms_fix( .not. nofix_sele )

        CALL energy_initialize
        CALL energy_non_bonding_options ( &
            list_cutoff   = 18.0_dp, &
            outer_cutoff  = 16.0_dp, &
            inner_cutoff  = 14.5_dp, &
            minimum_image = pbc )

    end subroutine

    !  GET_ATOM_NUMBERS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_atom_numbers()

        !--------------------------------------------------------------
        ! Obtain the atom numbers for constraints references
        !--------------------------------------------------------------

        integer                            :: i
        integer                            :: nres_tmp

        allocate(a_anum(a_natoms))

        do i=1, a_natoms
            read(a_atoms(i,2),*) nres_tmp
            a_anum(i) = atom_number(subsystem      = a_atoms(i,1), &
                                    residue_number = nres_tmp, &
                                    atom_name      = a_atoms(i,3) )
        end do

    end subroutine

    !  DISTANCE_CRD  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function distance_crd(n)

        !--------------------------------------------------------------
        ! Read distance for the 'n' constraint from coordinates
        !--------------------------------------------------------------

        implicit none

        integer, intent(in)                :: n
        real(8)                            :: distance_crd

        select case(c_type(n))
            case (1)
                distance_crd = geometry_distance(atmcrd, a_anum(c_atoms(n,1)), a_anum(c_atoms(n,2)))
            case (2)
                distance_crd = geometry_distance(atmcrd, a_anum(c_atoms(n,1)), a_anum(c_atoms(n,2))) &
                 + c_symm(n) * geometry_distance(atmcrd, a_anum(c_atoms(n,3)), a_anum(c_atoms(n,4)))
        end select

    end function

    !  DEFINE_CONSTRAINTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine define_constraints(print_file)

        !--------------------------------------------------------------
        ! Set-up constraints
        !--------------------------------------------------------------

        implicit none

        logical, intent(in)                :: print_file               ! Use file to print constraint value (dat_)

        integer                            :: i, j
        character(len=20)                  :: str_tmp

        if (.not. constr_flg) return
        if (.not. ALLOCATED(a_anum)) CALL get_atom_numbers

        CALL constraint_initialize

        do i=1, c_nconstr

            ! build file name (dat_)
            if (len_trim(c_file(i))==0) then
                write(str_tmp,*) i
                c_file(i) = "dat_" // trim(adjustl(str_tmp))
                do j=1, c_nconstr
                    write(str_tmp,*) c_indx(j)
                    c_file(i) = trim(c_file(i)) // "." // trim(str_tmp)
                end do
            end if

            ! define points
            do j=1, 4
                if (c_atoms(i,j) /= 0) CALL constraint_point_define( atom_selection(atom_number=(/ a_anum( c_atoms(i,j) ) /)) )
            end do

            ! use distance from .crd if requested
            if (c_dcrd(i)) then
                c_dist(i) = distance_crd(i)
            ! calculate distance if not from input
            else if (.not. c_dist_flg(i)) then
                c_dist(i) = c_dini(i) + c_step(i) * c_indx(i)
            end if

            ! define constraint
            select case (c_type(i))
                case (1)
                    if (print_file) then
                        CALL constraint_define(type = 'DISTANCE', &
                                               fc   = c_forc(i), &
                                               eq   = c_dist(i), &
                                               file = c_file(i))
                    else
                        CALL constraint_define(type = 'DISTANCE', &
                                               fc   = c_forc(i), &
                                               eq   = c_dist(i))
                    end if
                case (2)
                    if (print_file) then
                        CALL constraint_define(type    = 'MULTIPLE_DISTANCE', &
                                               fc      = c_forc(i), &
                                               eq      = c_dist(i), &
                                               weights = cof_sym(i,:), &
                                               file    = c_file(i))
                    else
                        CALL constraint_define(type    = 'MULTIPLE_DISTANCE', &
                                               fc      = c_forc(i), &
                                               eq      = c_dist(i), &
                                               weights = cof_sym(i,:))
                    end if
            end select

        end do

    end subroutine

    !  OUT_DIST_ENERGY  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine out_dist_energy(out_file)

        !--------------------------------------------------------------
        ! Write distances corresponding to constraint and energy
        !--------------------------------------------------------------

        character(len=256), intent(in)     :: out_file

        integer                            :: i
        logical                            :: f_exist

        if (.not. constr_flg) return
        if (.not. ALLOCATED(a_anum)) CALL get_atom_numbers

        ! check file existence to append or create new
        INQUIRE(file=out_file, exist=f_exist)
        if (f_exist) then
          open(200, file=out_file, form='formatted', status='old', position='append')
        else
          open(200, file=out_file, form='formatted', status='new')
        endif

        ! current distances
        do i=1, c_nconstr
          write(200, fmt='(f12.4,2X)', advance='no') distance_crd(i)
        end do

        ! total energy
        write(200, fmt='(f20.10,2X, f20.10,2X)', advance='no') etotal, eqm

        ! index number
        do i=1, c_nconstr
          write(200, fmt='(I5,2X)', advance='no') c_indx(i)
        end do

        ! reference distances
        do i=1, c_nconstr
          write(200, fmt='(f12.4,2X)', advance='no') c_dist(i)
        end do

        close(200)

    end subroutine

    !  CALCULATE_GRADIENT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculate_gradient()

        !--------------------------------------------------------------
        ! Calculate stable gradient
        !--------------------------------------------------------------

        integer                            :: i

        CALL gradient
        do i = 1, 1000
            if( grms < 40._dp ) cycle
            atmcrd = atmcrd - atmder / DSQRT( SUM( atmder **2 ) ) * 0.01_dp
            CALL gradient
        end do

    end subroutine

    !  CALCULATE_MINIMIZATION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculate_minimization()

        !--------------------------------------------------------------
        ! Perform CG and L-BFGS-B optimizations
        !--------------------------------------------------------------

        ! Minimize Conjugate-Gradient method
        CALL optimize_conjugate_gradient ( &
            step_number        = cg_steps, &
            print_frequency    = 2, &
            gradient_tolerance = cg_tolerance )

        ! Minimize L-BFGS-B
        CALL optimize_lbfgsb( &
            step_number        = lbfgsb_steps, &
            print_frequency    = 1, &
            gradient_tolerance = lbfgsb_tolerance )

    end subroutine

end module
