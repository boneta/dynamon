!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            OPTIONS, PARAMETERS AND INITIAL SUBROUTINES            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutines
!  -----------
!   DYNAMON_HEADER                Print the DYNAMON greeting and starts time
!   DYNAMON_FOOTER                Print elapsed time and normal ending
!   READ_OPTIONS                  Read options from an file in argument
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module INITIALIZATION

  implicit none

  !  DYNAMON VARIABLES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character(len=128), parameter      :: dynamon_version = '0.2.2'
  character(len=512)                 :: dynamon_path                  ! Installation path read from $DYNAMON env variable
  character(len=128)                 :: binaries_path = '/user/binaries/'  ! Relative location from dynamon_path to the binary files (.bin)
  integer                            :: t_ini, t_end, clock_rate      ! Elapsed time measurement

  !  OPTIONS & DEFAULTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character(len=128)                 :: MODE = ''                     ! Calculation mode

  character(len=256)                 :: name = ''                     ! Calculation name/title
  character(len=256)                 :: sys_bin = ''                  ! System binary file
  character(len=256)                 :: coord = ''                    ! Coordinate file (.crd)
  character(len=256)                 :: coord_name = ''               ! Coordinate file name without suffix

  character(len=128)                 :: cores = '1'                   ! Total number of cores (for Gaussian)
  character(len=128)                 :: memory = '3000MB'             ! Total RAM memory (for Gaussian)

  integer                            :: qm_charge = 0                 ! QM-region charge
  integer                            :: qm_multi = 1                  ! QM-region multiplicity

  character(len=50)                  :: semiemp = 'AM1'               ! Semi-empirical method
  logical                            :: gauss_flg = .false.           ! Call GAUSSIAN
  character(len=50)                  :: dft_func = 'M062X'            ! DFT functional
  character(len=50)                  :: dft_basis = '6-31+G(d,p)'     ! DFT basis set

  integer                            :: cg_steps = 10000              ! CG max number of steps
  real(8)                            :: cg_tolerance = 0.2            ! CG convergence
  integer                            :: lbfgsb_steps = 1000           ! L-BFGSB max number of steps
  real(8)                            :: lbfgsb_tolerance = 0.1        ! L-BFGSB convergence

  real(8)                            :: temp = 298.                   ! Temperature [k]
  real(8)                            :: md_step = 0.001               ! Time-step [ps]
  integer                            :: equilibration = 0             ! Equilibration [#steps]
  integer                            :: production = 1000             ! Production [#steps]
  integer                            :: dcd_freq = 100                ! Frequency to save trajectory
  character(len=256)                 :: velocities = ''               ! Velocities file to read

  logical                            :: pbc = .false.                 ! Periodic Boundary Conditions (minimum_image)

  integer                            :: loc_steps = 100               ! Number of steps for baker location
  real(8)                            :: loc_tolerance = 1.0           ! Location convergence
  logical                            :: ts_search = .false.           ! Search for a saddle point

  integer                            :: irc_dir = 0                   ! IRC direction {-1,1}
  integer                            :: irc_steps = 400               ! Maximum number of steps for IRC
  real(8)                            :: irc_dsp = 0.01                ! Displacement on every step of IRC

  integer                            :: int_nint = 0                  ! Number of residues to compute interactions with (aa + 1 H2O + ions)
  integer                            :: int_nres = 0                  ! Number of protein residues + ligands
  integer                            :: int_wbox(2) = 0               ! First and last residue number for Water Box atoms
  integer                            :: int_ions(2) = 0               ! First and last residue number for Ions atoms

  integer                            :: kie_atomn = 0                 ! Atom number to calculate KIE
  integer                            :: kie_skip = 0                  ! Number of frequencies to skip
  real(8)                            :: kie_mass = 2.01410177812D0    ! New mass to substitute in the atom
  character(len=256)                 :: kie_hess = 'update.dump'      ! Hessian file

  integer                            :: a_natoms = 0                  ! Number of input atoms
  character(10), allocatable         :: a_atoms(:,:)                  ! Input atoms (subsystem, residue_number, atom_name)
  integer, allocatable               :: a_anum(:)                     ! Atom numbers (from atom_number)

  logical                            :: constr_flg = .false.          ! Apply constraints
  integer                            :: c_nconstr = 0                 ! Number of contraints
  integer, allocatable               :: c_indx(:)                     ! Constraints indexes (i,j)
  integer, allocatable               :: c_type(:)                     ! Constraints type {1:'d', 2:'m'}
  integer, allocatable               :: c_symm(:)                     ! Constraints symmetry {-1,0,1}
  real(8), allocatable               :: c_forc(:)                     ! Constraints force [J·A^-n]
  integer, allocatable               :: c_atoms(:,:)                  ! Constrained atoms indexes
  logical, allocatable               :: c_dcrd(:)                     ! Use distance from structure
  real(8), allocatable               :: c_dini(:)                     ! Initial distance [A]
  real(8), allocatable               :: c_dend(:)                     ! Ending distance [A]
  real(8), allocatable               :: c_dist(:)                     ! Distance [A]
  logical, allocatable               :: c_dist_flg(:)                 ! Distance read from input
  real(8), allocatable               :: c_step(:)                     ! Scan step distance [A]
  character(256), allocatable        :: c_file(:)                     ! File naming (dat_)

  !  DYNAMO RECURRING VARIABLES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, allocatable               :: qm_sele(:)                    ! QM region atoms (aka 'acs')
  logical, allocatable               :: nofix_sele(:)                 ! No-fixed atoms (aka 'flg')

  real(8), parameter                 :: cof_sym(-1:1,2) = RESHAPE([1.D0, 0.D0, 1.D0, -1.D0, 0.D0, 1.D0],[3,2])

  contains

  !  DYNAMON_HEADER  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dynamon_header()

    write(*,fmt='(A80)') REPEAT('-',80)
    write(*,fmt='(30X,A22)')           '                    _ '
    write(*,fmt='(30X,A22)')           ' _|   _  _  _ _  _ | |'
    write(*,fmt='(30X,A22)')           '(_|\/| |(_|| | |(_)   '
    write(*,fmt='(30X,A22,20X,A1,A8)') '   /                  ','v',dynamon_version

    ! start time counting
    CALL SYSTEM_CLOCK(t_ini, clock_rate)

  end subroutine

  !  DYNAMON_FOOTER  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dynamon_footer()

    implicit none
    real(8)                            :: delta_t
    integer                            :: sec, min, hour, days
    real(8)                            :: msec

    ! stop time counting
    CALL SYSTEM_CLOCK(t_end)
    delta_t = real(t_end-t_ini)/clock_rate

    ! calculate times
    msec = MOD(delta_t, 1.)
    sec  = MOD(delta_t, 60.)
    min  = MOD(delta_t, 3600.) / 60
    hour = MOD(delta_t, 86400.) / 3600
    days = delta_t / 86400

    ! elapsed time and normal ending marker
    write(*,fmt='(/,A13,9X,I2,A,I0.2,A,I0.2,A,I0.2,F0.3)') 'Elapsed time:', days, '-', hour, ':', min, ':', sec, msec
    write(*,fmt='(A80)') REPEAT('<',80)

  end subroutine

  !  READ_OPTIONS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_options()

    !------------------------------------------------------------------
    ! Read options from a .dynn file pased as first argument
    !------------------------------------------------------------------

    implicit none

    integer                            :: i
    character(len=512)                 :: input_file
    character(len=20)                  :: option
    character(len=512)                 :: arg
    integer                            :: argn
    logical                            :: f_exist
    integer                            :: io_stat
    integer                            :: dot_position

    character(10), allocatable         :: a_atoms_tmp(:,:)
    integer, allocatable               :: c_indx_tmp(:)
    integer, allocatable               :: c_type_tmp(:)
    integer, allocatable               :: c_symm_tmp(:)
    real(8), allocatable               :: c_forc_tmp(:)
    integer, allocatable               :: c_atoms_tmp(:,:)
    logical, allocatable               :: c_dcrd_tmp(:)
    real(8), allocatable               :: c_dini_tmp(:)
    real(8), allocatable               :: c_dend_tmp(:)
    real(8), allocatable               :: c_dist_tmp(:)
    logical, allocatable               :: c_dist_flg_tmp(:)
    real(8), allocatable               :: c_step_tmp(:)
    character(256), allocatable        :: c_file_tmp(:)

    ! check no argument
    if (COMMAND_ARGUMENT_COUNT() == 0) then
      write(*,fmt='(A)') 'ERROR: No arguments found'
      STOP
    end if

    ! first argument as input file name
    CALL GETARG(1,input_file)
    ! check input_file not found
    INQUIRE(file=trim(input_file),exist=f_exist)
    if (f_exist) then
      write(*,fmt='(/,A,A)') "Options read from ", trim(input_file)
      open(unit=100, file=input_file, status='old', action='read', form='formatted')
      ! allocate arrays
      allocate( a_atoms(a_natoms,3), &
                c_indx(c_nconstr), &
                c_type(c_nconstr), &
                c_symm(c_nconstr), &
                c_forc(c_nconstr), &
                c_dcrd(c_nconstr), &
                c_dini(c_nconstr), &
                c_dend(c_nconstr), &
                c_dist(c_nconstr), &
                c_dist_flg(c_nconstr), &
                c_step(c_nconstr), &
                c_file(c_nconstr), &
                c_atoms(c_nconstr,4) )
      ! read line by line and asign general variables
      do
        read(100, *, iostat=io_stat) option
        if (io_stat/=0) exit
        if (option(1:1) == '!' .or. option(1:1) == '#') cycle
        backspace(100)
        select case (trim(option))   !FIXME: Make it case insensitive
          case ('MODE')
            read(100,*,iostat=io_stat) option, mode
          case ('NAME')
            read(100,*,iostat=io_stat) option, name
          case ('BIN')
            read(100,*,iostat=io_stat) option, sys_bin
          case ('COORD')
            read(100,*,iostat=io_stat) option, coord
          case ('CORES')
            read(100,*,iostat=io_stat) option, cores
          case ('MEMORY')
            read(100,*,iostat=io_stat) option, memory
          case ('CHARGE')
            read(100,*,iostat=io_stat) option, qm_charge
          case ('MULTI')
            read(100,*,iostat=io_stat) option, qm_multi
          case ('SEMIEMP')
            read(100,*,iostat=io_stat) option, semiemp
          case ('GAUSS')
            read(100,*,iostat=io_stat) option, gauss_flg
          case ('FUNC')
            read(100,*,iostat=io_stat) option, dft_func
          case ('BASIS')
            read(100,*,iostat=io_stat) option, dft_basis
          case ('CG_STEPS')
            read(100,*,iostat=io_stat) option, cg_steps
          case ('CG_TOLERANCE')
            read(100,*,iostat=io_stat) option, cg_tolerance
          case ('LBFGSB_STEPS')
            read(100,*,iostat=io_stat) option, lbfgsb_steps
          case ('LBFGSB_TOLERANCE')
            read(100,*,iostat=io_stat) option, lbfgsb_tolerance
          case ('TEMP')
            read(100,*,iostat=io_stat) option, temp
          case ('MD_STEP')
            read(100,*,iostat=io_stat) option, md_step
          case ('EQUI')
            read(100,*,iostat=io_stat) option, equilibration
          case ('PROD')
            read(100,*,iostat=io_stat) option, production
          case ('DCD_FREQ')
            read(100,*,iostat=io_stat) option, dcd_freq
          case ('VEL')
            read(100,*,iostat=io_stat) option, velocities
          case ('PBC')
            read(100,*,iostat=io_stat) option, pbc
          case ('LOC_STEPS')
            read(100,*,iostat=io_stat) option, loc_steps
          case ('LOC_TOLERANCE')
            read(100,*,iostat=io_stat) option, loc_tolerance
          case ('TS')
            read(100,*,iostat=io_stat) option, ts_search
          case ('IRC_DIR')
            read(100,*,iostat=io_stat) option, irc_dir
          case ('IRC_STEPS')
            read(100,*,iostat=io_stat) option, irc_steps
          case ('IRC_DSP')
            read(100,*,iostat=io_stat) option, irc_dsp
          case ('INT_NINT')
            read(100,*,iostat=io_stat) option, int_nint
          case ('INT_NRES')
            read(100,*,iostat=io_stat) option, int_nres
          case ('INT_WBOX')
            read(100,*,iostat=io_stat) option, int_wbox(:)
          case ('INT_IONS')
            read(100,*,iostat=io_stat) option, int_ions(:)
          case ('KIE_ATOMN')
            read(100,*,iostat=io_stat) option, kie_atomn
          case ('KIE_SKIP')
            read(100,*,iostat=io_stat) option, kie_skip
          case ('KIE_MASS')
            read(100,*,iostat=io_stat) option, kie_mass
          case ('KIE_HESS')
            read(100,*,iostat=io_stat) option, kie_hess
          case ('ATOM','A')
            a_natoms = a_natoms + 1
            allocate(a_atoms_tmp(a_natoms,3))
            a_atoms_tmp(:,:) = a_atoms(:,:)
            CALL move_alloc(a_atoms_tmp,a_atoms)
            read(100,*,iostat=io_stat) option, a_atoms(a_natoms,:)
          case ('CONSTR','C')
            constr_flg = .true.
            c_nconstr = c_nconstr + 1
            ! increment arrays
            allocate( c_indx_tmp(c_nconstr), &
                      c_type_tmp(c_nconstr), &
                      c_symm_tmp(c_nconstr), &
                      c_forc_tmp(c_nconstr), &
                      c_dcrd_tmp(c_nconstr), &
                      c_dini_tmp(c_nconstr), &
                      c_dend_tmp(c_nconstr), &
                      c_dist_tmp(c_nconstr), &
                      c_dist_flg_tmp(c_nconstr), &
                      c_step_tmp(c_nconstr), &
                      c_file_tmp(c_nconstr), &
                      c_atoms_tmp(c_nconstr,4) )
            c_indx_tmp(:) = 0
            c_type_tmp(:) = 0
            c_symm_tmp(:) = 0
            c_forc_tmp(:) = 0
            c_dcrd_tmp(:) = .false.
            c_dini_tmp(:) = 0
            c_dend_tmp(:) = 0
            c_dist_tmp(:) = 0
            c_dist_flg_tmp(:) = .false.
            c_step_tmp(:) = 0
            c_file_tmp(:) = ''
            c_atoms_tmp(:,:) = 0
            c_indx_tmp(:) = c_indx(:)
            c_type_tmp(:) = c_type(:)
            c_symm_tmp(:) = c_symm(:)
            c_forc_tmp(:) = c_forc(:)
            c_dcrd_tmp(:) = c_dcrd(:)
            c_dini_tmp(:) = c_dini(:)
            c_dend_tmp(:) = c_dend(:)
            c_dist_tmp(:) = c_dist(:)
            c_dist_flg_tmp(:) = c_dist_flg(:)
            c_step_tmp(:) = c_step(:)
            c_file_tmp(:) = c_file(:)
            c_atoms_tmp(:,:) = c_atoms(:,:)
            CALL move_alloc(c_indx_tmp,c_indx)
            CALL move_alloc(c_type_tmp,c_type)
            CALL move_alloc(c_symm_tmp,c_symm)
            CALL move_alloc(c_forc_tmp,c_forc)
            CALL move_alloc(c_dcrd_tmp,c_dcrd)
            CALL move_alloc(c_dini_tmp,c_dini)
            CALL move_alloc(c_dend_tmp,c_dend)
            CALL move_alloc(c_dist_tmp,c_dist)
            CALL move_alloc(c_dist_flg_tmp,c_dist_flg)
            CALL move_alloc(c_step_tmp,c_step)
            CALL move_alloc(c_file_tmp,c_file)
            CALL move_alloc(c_atoms_tmp,c_atoms)
            ! type of constraint
            read(100,*,iostat=io_stat) option, arg
            select case (trim(arg))
              case ('D','d')
                c_type(c_nconstr) = 1
              case ('M','m')
                c_type(c_nconstr) = 2
            end select
            ! read all constraint options
            do
              read(100, *, iostat=io_stat) option
              if (io_stat/=0) exit
              if (trim(option)=='CONSTR' .or. trim(option)=='C') exit
              if (option(1:1) == '!' .or. option(1:1) == '#') cycle
              backspace(100)
              select case (trim(option))
                case ('SYMM')
                  read(100,*,iostat=io_stat) option, c_symm(c_nconstr)
                case ('N')
                  read(100,*,iostat=io_stat) option, c_indx(c_nconstr)
                case ('ATOMS')
                  read(100,*,iostat=io_stat) option, c_atoms(c_nconstr,1:c_type(c_nconstr)*2)
                case ('FORCE')
                  read(100,*,iostat=io_stat) option, c_forc(c_nconstr)
                case ('DCRD')
                  read(100,*,iostat=io_stat) option, c_dcrd(c_nconstr)
                case ('DINIT')
                  read(100,*,iostat=io_stat) option, c_dini(c_nconstr)
                case ('DEND')
                  read(100,*,iostat=io_stat) option, c_dend(c_nconstr)
                case ('DIST','D')
                  c_dist_flg(c_nconstr) = .true.
                  read(100,*,iostat=io_stat) option, c_dist(c_nconstr)
                case ('STEP')
                  read(100,*,iostat=io_stat) option, c_step(c_nconstr)
                case ('DFILE')
                  read(100,*,iostat=io_stat) option, c_file(c_nconstr)
                case DEFAULT
                  read(100,*) arg
                  write(*,fmt='(A,A)') "Omitting unknown option: ", arg
                end select
              end do
          case DEFAULT
            read(100,*) arg
            write(*,fmt='(A,A)') "Omitting unknown option: ", arg
        end select
      end do
      close(unit=100)
      argn = 2
    else
      argn = 1
    end if

    ! read command line arguments
    write(*,*)
    do while (argn <= COMMAND_ARGUMENT_COUNT())
      CALL GETARG(argn, option)
      CALL GETARG(argn+1, arg)
      write(*,fmt='(A,5X,A8,2X,A)') "Argument:", ADJUSTL(option(3:)), trim(arg)
      select case (trim(option))
        case ('--MODE')
          read(arg,'(A)',iostat=io_stat) mode
        case ('--NAME')
          read(arg,'(A)',iostat=io_stat) name
        case ('--BIN')
          read(arg,'(A)',iostat=io_stat) sys_bin
        case ('--COORD')
          read(arg,'(A)',iostat=io_stat) coord
        case ('--CORES')
          read(arg,*,iostat=io_stat) cores
        case ('--MEMORY')
          read(arg,*,iostat=io_stat) memory
        case ('--CHARGE')
          read(arg,*,iostat=io_stat) qm_charge
        case ('--MULTI')
          read(arg,*,iostat=io_stat) qm_multi
        case ('--SEMIEMP')
          read(arg,'(A)',iostat=io_stat) semiemp
        case ('--GAUSS')
          read(arg,*,iostat=io_stat) gauss_flg
        case ('--FUNC')
          read(arg,'(A)',iostat=io_stat) dft_func
        case ('--BASIS')
          read(arg,'(A)',iostat=io_stat) dft_basis
        case ('--TEMP')
          read(arg,*,iostat=io_stat) temp
        case ('--EQUI')
          read(arg,*,iostat=io_stat) equilibration
        case ('--PROD')
          read(arg,*,iostat=io_stat) production
        case ('--VEL')
          read(arg,'(A)',iostat=io_stat) velocities
        case ('--TS')
          read(arg,*,iostat=io_stat) ts_search
        case ('--IRC_DIR')
          read(arg,*,iostat=io_stat) irc_dir
        case ('--N')
          do i=1, c_nconstr
            CALL GETARG(argn+i, arg)
            read(arg,*,iostat=io_stat) c_indx(i)
          end do
          argn = argn + c_nconstr - 1
        case ('--DIST')
          c_dist_flg = .true.
          do i=1, c_nconstr
            CALL GETARG(argn+i, arg)
            read(arg,*,iostat=io_stat) c_dist(i)
          end do
          argn = argn + c_nconstr - 1
      end select
      argn = argn + 2
    end do

    ! get DYNAMON installation path from environment variable
    CALL GET_ENVIRONMENT_VARIABLE('DYNAMON', dynamon_path)

    ! if binary not found locally, check in binaries folder
    if (Len_Trim(sys_bin) > 0) then
      INQUIRE(file=trim(sys_bin),exist=f_exist)
      if (.not. f_exist .and. Len_Trim(dynamon_path) > 0) then
        sys_bin = trim(dynamon_path)//trim(binaries_path)//trim(sys_bin)
      end if
    end if

    ! check neccesary inputs
    if (Len_Trim(mode) == 0 .or. Len_Trim(sys_bin) == 0 .or. Len_Trim(coord) == 0) then
      write(*,fmt='(A)') 'ERROR: Missing mandatory options (MODE/BIN/COORD)'
      STOP
    end if

    ! get coord without extension
    dot_position = SCAN(trim(coord),'.', back= .true.)
    if (dot_position > 0) then
      coord_name = coord(1:dot_position-1)
    else
      coord_name = coord
    end if

  end subroutine

end module
