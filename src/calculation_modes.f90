!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                         CALCULATION MODES                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutines
!  -----------
!   PRINT_MODE                    Print the calculation mode in a formatted fashion
!   BUILD_BIN                     Build a binary file from sequence and FF
!   DYNAMON_SP                    Single-point energy calculation
!   DYNAMON_MINIMIZATION          Structure optimization
!   DYNAMON_LOCATE                Localize and characterize minima or saddle point
!   DYNAMON_IRC                   Internal Reaction Coordinate from TS
!   DYNAMON_MD                    BO Molecular Dynamics simulation
!   DYNAMON_INTERACTION           Electrostatic and Lennard-Jones interactions
!   DYNAMON_KIE                   Kinetic Isotope Effect
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module CALCULATION_MODES

    use dynamo
    use initialization
    use common

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

    !  BUILD_BIN  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine build_bin(exit_after)

        logical, intent(in)                :: exit_after

        character(len=256)                 :: ffname
        integer                            :: dot_position

        ! get FF name without suffix
        ffname = fffile
        dot_position = SCAN(trim(ffname), '.', back=.true.)
        if (dot_position > 0) ffname = ffname(1:dot_position-1)

        ! name of formatted FF file
        ffname = trim(ffname) // '_processed.ff'

        CALL mm_file_process(trim(ffname), trim(fffile))
        CALL mm_system_construct(trim(ffname), trim(seqfile))
        CALL mm_system_write(trim(binfile))

        if (Len_Trim(coord) > 0) CALL coordinates_read(trim(coord))

        if (exit_after) then
            CALL dynamo_footer
            CALL dynamon_footer
            STOP
        end if

    end subroutine

    !  DYNAMON_SP  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dynamon_sp()

        if (len_trim(name)==0) name = trim(coord_name) // "-sp"

        CALL energy

        if (mode == 'CORR') CALL out_dist_energy(trim(name)//".out")

    end subroutine

    !  DYNAMON_MINIMIZATION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dynamon_minimization()

        if (len_trim(name)==0) name = trim(coord_name) // "-mini"

        CALL define_constraints(print_file=.false.)
        CALL calculate_gradient
        CALL calculate_minimization

        if (mode == 'SCAN' .or. mode == 'PES') CALL out_dist_energy(trim(name)//".out")
        CALL coordinates_write(trim(name)//".crd")
        CALL pdb_write(trim(name)//".pdb")

    end subroutine

    !  DYNAMON_LOCATE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dynamon_locate()

        use panadero

        integer                            :: i, j
        real(8), allocatable               :: x(:)

        if (len_trim(name)==0) name = trim(coord_name) // "-loc"

        use_hessian_numerical = .false.
        use_hessian_recalc = 1000
        use_hessian_method = "BOFILL"

        nofix_sele = nofix_sele .and. .not. qm_sele
        CALL atoms_fix( .not. ( nofix_sele .or. qm_sele ) )
        CALL energy

        CALL setup_optimization( &
            core_atoms         = qm_sele, &
            envi_atoms         = nofix_sele, &
            step_number        = 1000, &
            gradient_tolerance = 0.5_dp, &
            print_frequency    = 1, &
            trajectory         = "snapshot.dcd" )

        allocate( x(1:3*count(qm_sele)) )
        j = -3
        do i = 1, natoms
            if( qm_sele(i) ) then
            j = j + 3
            x(j+1:j+3) = atmcrd(1:3,i)
            end if
        end do

        CALL baker_search( &
            fcalc              = egh_calc, &
            x                  = x, &
            print_frequency    = 1, &
            step_number        = loc_steps, &
            maximum_step       = 0.1_dp, &
            locate_saddle      = ts_search, &
            gradient_tolerance = loc_tolerance )

        j = -3
        do i = 1, natoms
            if( qm_sele(i) ) then
                j = j + 3
                atmcrd(1:3,i) = x(j+1:j+3)
            end if
        end do

        CALL cleanup_optimization

        CALL coordinates_write(trim(name)//".crd")
        CALL pdb_write(trim(name)//".pdb")

        CALL atoms_fix( .not. ( nofix_sele .or. qm_sele ) )
        CALL gradient
        j = next_unit()
        open( unit = j, file = "forces.dump", action = "write", form = "formatted" )
        do i = 1, natoms
            write( j, "(4f20.10)") atmder(1:3,i), sqrt(sum(atmder(1:3,i)*atmder(1:3,i))/3._dp)
        end do
        close( j )

        CALL atoms_fix( .not. qm_sele )
        use_hessian_recalc = 1
        CALL hessian
        open( unit = j, file = "hessian.dump", action = "write", form = "unformatted" )
        write( j ) atmhes
        close( j )
        CALL project_rotation_translation( atmhes, .true. )
        CALL normal_mode_frequencies( atmhes )
        do i = 1, 10
            CALL normal_mode_view( i, afact = 4._dp )
        end do

    end subroutine

    !  DYNAMON_IRC  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dynamon_irc()

        implicit none

        integer                            :: i, j, k
        integer                            :: fd
        integer                            :: it1, it2
        integer                            :: it1_max, it2_max
        logical                            :: xlg
        integer                            :: acs_n, acs_n3, acs_nh
        integer, allocatable               :: acs_i(:)
        real( kind=dp )                    :: max_rms, cur_rms, tmp, g_nrm
        real( kind=dp )                    :: lamb, l1, l2
        real( kind=dp ), allocatable       :: evec(:,:)
        real( kind=dp ), allocatable       :: eval(:), x_cur(:), g_cur(:), vec(:), mw(:), x_ref(:), x_1st(:)
        type( dcd_type )                   :: dcd
        real( kind=dp ), parameter         :: large = 1000000._dp
        real( kind=dp ), parameter         :: stp = 150._dp
        real( kind=dp ), parameter         :: l_tol = 0.00000001_dp

        it1_max = irc_steps
        it2_max = irc_steps + 1000

        if (len_trim(name)==0) then
            if (irc_dir == -1) then
                name = trim(coord_name) // "-irc-back"
            else if (irc_dir == 1) then
                name = trim(coord_name) // "-irc-for"
            else
                CALL dynn_log(1, 'Unknown IRC direction')
            end if
        end if

        ! file to dump energies
        open(unit=900, file=trim(name)//".dat", form='formatted')

        CALL dcd_initialize( dcd )
        CALL dcd_activate_write( &
            file    = trim(name)//".dcd", &
            dcd     = dcd, &
            type    = "CORD", &
            natoms  = natoms, &
            nfixed  = nfixed, &
            nframes = it1_max, &
            qfix    = atmfix )

        acs_n = count( qm_sele )
        acs_n3 = 3 * acs_n
        acs_nh = acs_n3 * ( acs_n3 + 1 ) / 2
        allocate( acs_i(1:acs_n), mw(1:acs_n3) )

        j = 0
        do i = 1, natoms
            if( qm_sele(i) ) then
                j = j + 1
                acs_i(j) = i
                k = 3 * ( j - 1 )
                mw(k+1:k+3) = 1._dp / sqrt( atmmas(i) )
            end if
        end do

        CALL atoms_fix( .not. nofix_sele )
        CALL energy
        tmp = etotal

        use_hessian_numerical = .false.
        use_hessian_recalc = 1000
        use_hessian_method = "BOFILL"

        allocate( eval(1:acs_n3), evec(1:acs_n3,1:acs_n3) )
        allocate( x_cur(1:acs_n3), g_cur(1:acs_n3), vec(1:acs_n3), x_ref(1:acs_n3), x_1st(1:acs_n3) )

        CALL atoms_fix( .not. qm_sele )
        CALL hessian
        k = 0
        do i = 1, acs_n3
            do j = 1, i
                k = k + 1
                atmhes(k) = atmhes(k) * mw(i) * mw(j)
            end do
        end do
        do i = 1, acs_n
            j = 3 * ( i - 1 )
            x_cur(j+1:j+3) = atmcrd(1:3,acs_i(i)) / mw(j+1)
        end do
        CALL project_hrt( acs_n3, x_cur, mw, atmhes )
        CALL symmetric_upper( atmhes, eval, evec )
        write(*,"(8f10.2)") eval
        write(*,"(a20,2f20.10)") ">>>  RP:", 0.0_dp, tmp

        CALL atoms_fix( .not. nofix_sele )
        CALL gradient
        max_rms = .0_dp
        cur_rms = .0_dp
        do i = 1, natoms
            tmp = sum( atmder(1:3,i) ** 2 )
            cur_rms = cur_rms + tmp
            if( tmp > max_rms ) max_rms = tmp
        end do
        cur_rms = sqrt( cur_rms / real( 3 * nfree, kind = dp ) )
        max_rms = sqrt( max_rms / 3._dp )
        write(*,"(a20,2f20.10)") ">>> RMS:", cur_rms, max_rms
        CALL dcd_write( dcd, atmcrd, boxl )
        do i = 1, acs_n
            j = 3 * ( i - 1 )
            x_ref(j+1:j+3) = atmcrd(1:3,acs_i(i)) / mw(j+1)
        end do

        x_cur(1:acs_n3) = evec(1:acs_n3,1) / sqrt( sum( evec(1:acs_n3,1) * evec(1:acs_n3,1) ) ) * irc_dsp * irc_dir
                x_1st(1:acs_n3) = x_cur(1:acs_n3)
        do i = 1, acs_n
            j = 3 * ( i - 1 )
            atmcrd(1:3,acs_i(i)) = atmcrd(1:3,acs_i(i)) + x_cur(j+1:j+3) * mw(j+1)
        end do

        CALL atoms_fix( .not. nofix_sele .or. qm_sele )
        CALL optimize_conjugate_gradient( &
            step_size          = 0.01_dp, &
            print_frequency    = 10, &
            gradient_tolerance = 1.0_dp, &
            step_number        = 1000 )

        CALL atoms_fix( .not. nofix_sele )
        CALL gradient
        max_rms = .0_dp
        cur_rms = .0_dp
        do i = 1, natoms
            tmp = sum( atmder(1:3,i) ** 2 )
            cur_rms = cur_rms + tmp
            if( tmp > max_rms ) max_rms = tmp
        end do
        cur_rms = sqrt( cur_rms / real( 3 * nfree, kind = dp ) )
        max_rms = sqrt( max_rms / 3._dp )
        write(*,"(a20,2f20.10)") ">>> RMS:", cur_rms, max_rms
        CALL dcd_write( dcd, atmcrd, boxl )
        write(*,"(a20,2f20.10)") ">>>  RP:", sqrt( sum( x_cur ** 2 ) ), etotal

        write(900, fmt='(A6,3X,F20.10)') 0, etotal
        c_indx = 0
        CALL out_dist_energy(trim(name)//".out")

        it1 = 1
        do while( ( it1 < it1_max ) .and. &
            ( ( cur_rms > 1.49_dp .or. max_rms > 2.23_dp ) .or. it1 < 10 ) )

            do i = 1, acs_n
                j = 3 * ( i - 1 )
                g_cur(j+1:j+3) = atmder(1:3,acs_i(i)) * mw(j+1)
                x_cur(j+1:j+3) = atmcrd(1:3,acs_i(i)) / mw(j+1)
            end do
            CALL project_grt( acs_n3, x_cur, mw, g_cur )
            CALL atoms_fix( .not. qm_sele )
            CALL hessian( .false. )
            k = 0
            do i = 1, acs_n3
                do j = 1, i
                    k = k + 1
                    atmhes(k) = atmhes(k) * mw(i) * mw(j)
                end do
            end do
            CALL project_hrt( acs_n3, x_cur, mw, atmhes )
            CALL symmetric_upper( atmhes, eval, evec )
            g_nrm = sqrt( dot_product( g_cur, g_cur ) )
            g_cur = g_cur / g_nrm
            eval = eval / g_nrm
            write(*,"(8f10.2)") eval
            do i = 1, acs_n3
                vec(i) = sum( evec(1:acs_n3,i) * g_cur(1:acs_n3) )
            end do
            lamb = .0_dp
            if( eval(1) < 0.0_dp ) then
                lamb = eval(1) - stp
                l1 = eval(1)
                l2 = - large
            end if

            it2 = 1
            xlg = .true.
            do while( it2 < it2_max .and. xlg )
                tmp = .0_dp
                do i = 1, acs_n3
                    tmp = tmp + vec(i) * vec(i) / ( lamb - eval(i) )
                end do
                xlg = dabs( lamb - tmp ) > l_tol
                if( xlg ) then
                    if( eval(1) > 0.0_dp ) then
                        lamb = tmp
                    else
                        if( tmp < lamb ) l1 = lamb
                        if( tmp > lamb ) l2 = lamb
                        if( l2 > - large ) then
                            lamb = 0.5_dp * ( l1 + l2 )
                        else if( l2 == - large ) then
                            lamb = lamb - stp
                        end if
                    end if
                end if
                it2 = it2 + 1
            end do
            write(*,"(2i6,g16.8)") it1, it2, lamb
            x_cur = matmul( evec, vec / ( lamb - eval ) )
            tmp = sqrt( dot_product( x_cur, x_cur ) )
            if( tmp > irc_dsp ) x_cur = irc_dsp / tmp * x_cur
            if( sum( x_cur * x_1st ) < .0_dp ) x_cur = -x_cur
            do i = 1, acs_n
                j = 3 * ( i - 1 )
                atmcrd(1:3,acs_i(i)) = atmcrd(1:3,acs_i(i)) + x_cur(j+1:j+3) * mw(j+1)
                x_cur(j+1:j+3) = atmcrd(1:3,acs_i(i)) / mw(j+1)
            end do

            CALL atoms_fix( .not. nofix_sele .or. qm_sele )
            CALL optimize_conjugate_gradient( &
                step_size          = 0.01_dp, &
                print_frequency    = 10, &
                gradient_tolerance = 1.0_dp, &
                step_number        = 1000 )

            CALL atoms_fix( .not. nofix_sele )
            CALL gradient
            max_rms = .0_dp
            cur_rms = .0_dp
            do i = 1, natoms
                tmp = sum( atmder(1:3,i) ** 2 )
                cur_rms = cur_rms + tmp
                if( tmp > max_rms ) max_rms = tmp
            end do
            cur_rms = sqrt( cur_rms / real( 3 * nfree, kind = dp ) )
            max_rms = sqrt( max_rms / 3._dp )
            write(*,"(a20,2f20.10)") ">>> RMS:", cur_rms, max_rms
            write(*,"(a20,2f20.10)") ">>>  RP:", sqrt( sum( ( x_cur - x_ref ) ** 2 ) ), etotal
            CALL dcd_write( dcd, atmcrd, boxl )

            write(900, fmt='(A6,3X,F20.10)') it1, etotal
            c_indx = irc_dir * it1
            CALL out_dist_energy(trim(name)//".out")

            it1 = it1 + 1
        end do

        CALL flush( dcd%unit )
        close( dcd%unit )

        CALL coordinates_write( trim(name)//".crd" )
        CALL pdb_write( trim(name)//".pdb" )

        deallocate( eval, evec )
        deallocate( x_cur, g_cur, vec, mw, x_ref )

        contains                                                       ! MEP_PROJECT

        subroutine project_grt( n, x, w, g )
            implicit none
            integer, intent(in) :: n
            real*8, dimension(1:n), intent(in) :: x, w
            real*8, dimension(1:n), intent(inout) :: g

            real*8, dimension(:,:), allocatable :: gc
            real*8, dimension(1:3) :: mc
            real*8 :: t
            integer :: i, j, k

            allocate( gc(1:n,1:6) )
            mc = 0.0d0
            t = 0.0d0
            do i = 1, n / 3
                j = 3 * ( i - 1 )
                t = t + ( 1.d0 / w(j+1) ) ** 2
                mc(1:3) = mc(1:3) + x(j+1:j+3) / w(j+1)
            end do
            mc = mc / t
            gc = 0.0d0
            do i = 1, n / 3
                j = 3 * ( i - 1 )
                gc(j+1,1) = 1.d0 / w(j+1)
                gc(j+2,2) = 1.d0 / w(j+2)
                gc(j+3,3) = 1.d0 / w(j+3)
                gc(j+2,4) = -( x(j+3) - mc(3) / w(j+1) )
                gc(j+3,4) =  ( x(j+2) - mc(2) / w(j+1) )
                gc(j+1,5) =  ( x(j+3) - mc(3) / w(j+1) )
                gc(j+3,5) = -( x(j+1) - mc(1) / w(j+1) )
                gc(j+1,6) = -( x(j+2) - mc(2) / w(j+1) )
                gc(j+2,6) =  ( x(j+1) - mc(1) / w(j+1) )
            end do
            do i = 1, 6
                do j = 1, i - 1
                    t = dot_product( gc(1:n,i), gc(1:n,j) )
                    gc(1:n,i) = gc(1:n,i) - t * gc(1:n,j)
                end do
                t = sqrt( dot_product( gc(1:n,i), gc(1:n,i) ) )
                gc(1:n,i) = gc(1:n,i) / t
            end do
            ! G' = G - Tx * G � Tx - ... - Rx * G � Rx - ...
            do i = 1, 6
                g(1:n) = g(1:n) - dot_product( gc(1:n,i), g(1:n) ) * gc(1:n,i)
            end do
            deallocate( gc )
        end subroutine

        subroutine project_hrt( n, x, w, h )
            implicit none
            integer, intent(in) :: n
            real*8, dimension(1:n), intent(in) :: x, w
            real*8, dimension(1:n*(n+1)/2), intent(inout) :: h

            real*8, dimension(:,:), allocatable :: gc, hv
            real*8, dimension(1:6,1:6) :: vhv
            real*8, dimension(1:3) :: mc
            real*8 :: t
            integer :: i, ii, j, jj, k

            allocate( gc(1:n,1:6), hv(1:n,1:6) )
            mc = 0.0d0
            t = 0.0d0
            do i = 1, n / 3
                j = 3 * ( i - 1 )
                t = t + ( 1.d0 / w(j+1) ) ** 2
                mc(1:3) = mc(1:3) + x(j+1:j+3) / w(j+1)
            end do
            mc = mc / t
            gc = 0.0d0
            do i = 1, n / 3
                j = 3 * ( i - 1 )
                gc(j+1,1) = 1.d0 / w(j+1)
                gc(j+2,2) = 1.d0 / w(j+1)
                gc(j+3,3) = 1.d0 / w(j+1)
                gc(j+2,4) = -( x(j+3) - mc(3) / w(j+1) )
                gc(j+3,4) =  ( x(j+2) - mc(2) / w(j+1) )
                gc(j+1,5) =  ( x(j+3) - mc(3) / w(j+1) )
                gc(j+3,5) = -( x(j+1) - mc(1) / w(j+1) )
                gc(j+1,6) = -( x(j+2) - mc(2) / w(j+1) )
                gc(j+2,6) =  ( x(j+1) - mc(1) / w(j+1) )
            end do
            do i = 1, 6
                do j = 1, i - 1
                    t = dot_product( gc(1:n,i), gc(1:n,j) )
                    gc(1:n,i) = gc(1:n,i) - t * gc(1:n,j)
                end do
                t = sqrt( dot_product( gc(1:n,i), gc(1:n,i) ) )
                gc(1:n,i) = gc(1:n,i) / t
            end do
            !  P = I - Tx � Tx - ... - Rx � Rx - ...
            ! H' = P � H � P
            do ii = 1, 6
                do i = 1, n
                    k = i * ( i - 1 ) / 2
                    t = 0.0d0
                    do j = 1, i
                        k = k + 1
                        t = t + h(k) * gc(j,ii)
                    end do
                    do j = i + 1, n
                        k = k + j - 1
                        t = t + h(k) * gc(j,ii)
                    end do
                    hv(i,ii) = t
                end do
                do jj = 1, ii
                    vhv(jj,ii) = dot_product( gc(1:n,ii), hv(1:n,jj) )
                    vhv(ii,jj) = vhv(jj,ii)
                end do
            end do
            k = 0
            do i = 1, n
                do j = 1, i
                    t = 0.0d0
                    do ii = 1, 6
                        t = t - hv(i,ii) * gc(j,ii) - hv(j,ii) * gc(i,ii)
                        do jj = 1, 6
                            t = t + gc(j,jj) * vhv(jj,ii) * gc(i,ii)
                        end do
                    end do
                    k = k + 1
                    h(k) = h(k) + t
                end do
            end do
            deallocate( gc, hv )
        end subroutine

    end subroutine

    !  DYNAMON_MD  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dynamon_md()

        integer                            :: my_random
        real(8)                            :: fortran_random

        if (len_trim(name)==0) name = trim(coord_name) // "-md"

        CALL define_constraints(print_file=.true.)
        CALL calculate_gradient

        ! read velocities or random
        if (len_trim(velocities) /= 0) then
            CALL velocity_read( trim(velocities) )
        else
            CALL RANDOM_NUMBER(fortran_random)
            CALL random_initialize( INT(fortran_random) + my_random() )
            CALL velocity_assign( temp, .false. )
        end if

        ! equilibration
        if (equilibration > 0) then
            CALL dynamics_options( &
                time_step       = md_step, &
                print_frequency = 100, &
                steps           = equilibration )
            CALL langevin_verlet_dynamics( temp, 100._dp )
        end if

        ! production
        CALL dynamics_options( &
            time_step       = md_step, &
            print_frequency = 100, &
            save_frequency  = dcd_freq, &
            coordinate_file = trim(name)//".dcd", &
            steps           = production )

        if (constr_flg) CALL constraint_writing_start
        CALL langevin_verlet_dynamics( temp, 100._dp )
        if (constr_flg) CALL constraint_writing_stop

        CALL coordinates_write(trim(name)//".crd")
        CALL velocity_write(trim(name)//".vel")

    end subroutine

    !  DYNAMON_INTERACTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dynamon_interaction()

        use utils, only                    : append_1d

        implicit none

        integer                            :: i, j
        character(len=128)                 :: fmtString
        character(len=128)                 :: dens_file
        type(dcd_type)                     :: trj
        integer                            :: inter
        integer                            :: ninter
        logical                            :: wbox_atoms(1:natoms), ions_atoms(1:natoms)
        integer, allocatable               :: qm_res(:), wbox_res(:), ions_res(:)
        real(kind=dp), allocatable         :: my_eqm(:), my_lj(:)
        real(kind=dp)                      :: my_atmeps(natoms), &     ! atom LJ epsilon
                                              my_atmchg(natoms)        ! atom charges

        ! add prefix to output file if name set
        if (len_trim(name)>0) name = trim(name) // "-"

        ! output files
        dens_file = trim(name)//'dens.mat'
        open( 500, file=trim(name)//'intmatrix_eqm.dat' )
        open( 600, file=trim(name)//'intmatrix_lj.dat'  )
        open( 700, file=trim(name)//'intmatrix_tot.dat' )
        open( 800, file=trim(name)//'eqm_trj.dat' )
        open( 900, file=trim(name)//'seq.resid' )

        ! get global residue numbers that contain QM atoms
        CALL res_in_atom_array(qm_sele, qm_res)
        write(*,*) "@@ RESIDUES TREATED QM: ", qm_res

        ! get residue number that match subsystem name for wbox and ions
        wbox_atoms = ATOM_SELECTION( subsystem = wbox_sname )
        CALL res_in_atom_array(wbox_atoms, wbox_res)
        ions_atoms = ATOM_SELECTION( subsystem = ions_sname )
        CALL res_in_atom_array(ions_atoms, ions_res)

        ! count number of interactions
        ! #res_total - #qm_res - #wbox_res - #ions_res + 2
        ninter = nresid - SIZE(qm_res) - SIZE(wbox_res) - SIZE(ions_res) + 2
        write(*,*) "@@ NUMBER OF INTERACTIONS: ", ninter

        ! write sequence of residues to calculate interactions
        do i = 1, ninter
          if( ANY(qm_res==i) .or. ANY(wbox_res==i) .or. ANY(ions_res==i)  ) CYCLE
          write(900,'(A,4X,I6)') trim(resnam(i)), i
        end do
        write(900,'(A,X,I6)') "WATERS", SIZE(wbox_res)
        write(900,'(A,X,I6)') "IONS  ", SIZE(ions_res)
        flush(900)

        allocate(my_eqm(ninter), my_lj(ninter))
        write(fmtString,*) ninter
        fmtString = '('//trim(adjustl(fmtString))//'f20.10)'

        CALL atoms_fix( .not. qm_sele )

        atmeps14  = .0_dp
        my_atmeps = atmeps
        atmchg14  = .0_dp
        my_atmchg = atmchg

        CALL dcd_initialize( trj )
        CALL dcd_activate_read( trim(int_dcd), trj )
        CALL dcd_read( trj, atmcrd, boxl )

        write(*,*) "@@ NUMBER OF FRAMES: ", trj%nframes
        do i = 2, trj%nframes
            ! read next trajectory frame
            CALL dcd_read( trj, atmcrd, boxl )
            ! skip calculation of frame if not multiple of stride
            if (MOD(i,dcd_stride) /= 0) CYCLE

            write(*,*) "@@ FRAME: ", i

            atmeps = my_atmeps
            atmchg = my_atmchg
            CALL mopac_scf_initialize
            CALL energy
            CALL density_write(trim(dens_file))
            CALL mopac_scf_options( iterations = -1 )
            atmeps = .0_dp
            atmchg = .0_dp
            CALL density_read(trim(dens_file))
            CALL energy
            write( 800, '(f20.10)' ) eqm
            CALL flush( 800 )

            inter = 1
            do j = 1, nresid
                if( ANY(qm_res==j) ) CYCLE    ! skip QM residues
                if( ANY(wbox_res==j) ) CYCLE  ! skip WBulk residues
                if( ANY(ions_res==j) ) CYCLE  ! skip ions residues
                write(*,*) "## INT #  FR: ", i, " INT: ", inter, " RES: ", j
                ! == Calculate interaction ==
                atmeps = .0_dp
                atmchg = .0_dp
                ! Recover QM(Eps) vector
                ! FIXME: To be checked. Now done only for QM-atoms. Better with whole qm-residue??
                atmeps = MERGE(my_atmeps, atmeps, qm_sele)
                ! Recover MM(Chrg/Eps) vector for the selected residue
                atmchg( resind(j)+1 : resind(j+1) ) = my_atmchg( resind(j)+1 : resind(j+1) )
                atmeps( resind(j)+1 : resind(j+1) ) = my_atmeps( resind(j)+1 : resind(j+1) )
                ! Load Molecular Orbital and calc energy (without SCF)
                CALL density_read(trim(dens_file))
                CALL energy
                my_eqm(inter) = eqm
                my_lj(inter)  = qmlj
                inter = inter + 1
            end do

            ! WBulk residues
            write(*,*) "## INT #  FR: ", i, " INT: ", inter, " RES: WATERS"
            atmeps = .0_dp
            atmchg = .0_dp
            atmeps = MERGE(my_atmeps, atmeps, qm_sele)
            atmchg = MERGE(my_atmchg, atmchg, wbox_atoms)
            atmeps = MERGE(my_atmeps, atmeps, wbox_atoms)
            ! -- en princpio no necesario por que al definirlo QM las cargas se hacen cero... (??)
            CALL density_read(trim(dens_file))
            CALL energy
            my_eqm(inter) = eqm
            my_lj(inter)  = qmlj
            inter = inter + 1

            ! ions residues
            write(*,*) "## INT #  FR: ", i, " INT: ", inter, " RES: IONS"
            atmeps = .0_dp
            atmchg = .0_dp
            ! atmeps = MERGE(my_atmeps, atmeps, qm_sele)
            atmeps = MERGE(my_atmeps, .0_dp, qm_sele)
            atmchg = MERGE(my_atmchg, atmchg, ions_atoms)
            atmeps = MERGE(my_atmeps, atmeps, ions_atoms)
            CALL density_read(trim(dens_file))
            CALL energy
            my_eqm(inter) = eqm
            my_lj(inter)  = qmlj

            write( 500, fmtString ) my_eqm(1:ninter)
            write( 600, fmtString ) my_lj(1:ninter)
            write( 700, fmtString ) my_eqm(1:ninter)+my_lj(1:ninter)
            CALL flush( 500 )
            CALL flush( 600 )
            CALL flush( 700 )

        end do

        close(500)
        close(600)
        close(700)
        close(800)
        close(900)

        CALL dcd_deactivate( trj )

        contains

        subroutine res_in_atom_array(atom_array, res_array)

            !----------------------------------------------------------
            ! Return an array of global residue numbers that contain
            ! any atom that is True in an atom array
            !----------------------------------------------------------

            implicit none
            logical, intent(in)                  :: atom_array(natoms)
            integer, allocatable, intent(out)    :: res_array(:)
            integer                              :: i, j

            if (.not. allocated(res_array)) allocate(res_array(0))
            do i = 1, nresid
                do j = resind(i)+1, resind(i+1)
                    if (atom_array(j)) then
                        CALL append_1d(res_array, i)
                        EXIT
                    end if
                end do
            end do

        end subroutine

    end subroutine

    !  DYNAMON_KIE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dynamon_kie()

        implicit none

        integer                            :: n
        integer                            :: nres_tmp
        integer                            :: kie_anum
        real( kind=dp )                    :: pres = 1.D0
        real( kind=dp )                    :: gxh, gxd

        ! obtain atom number
        read(kie_atom(2),*) nres_tmp
        kie_anum = atom_number(subsystem      = kie_atom(1), &
                               residue_number = nres_tmp, &
                               atom_name      = kie_atom(3) )

        n = count( qm_sele )
        allocate( atmhes(1:3*n*(3*n+1)/2) )
        CALL atoms_fix( .not. qm_sele )

        ! standard energy
        open(unit=999, file=trim( kie_hess ), action="read", form="unformatted")
        read( 999 ) atmhes
        close( 999 )
        CALL gibbs_energy( pres, temp, 1.0d0, kie_skip, gxh )

        ! change masses
        atmmas(kie_anum) = kie_mass

        ! isotopic reactants/ts energy
        CALL coordinates_read( trim( coord ) )
        open(unit=999, file=trim( kie_hess ), action="read", form="unformatted")
        read( 999 ) atmhes
        close( 999 )
        CALL gibbs_energy( pres, temp, 1.0d0, kie_skip, gxd )

        deallocate( atmhes )

    end subroutine

end module
