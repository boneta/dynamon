!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                        SELECTION PROCEDURE                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutines
!  -----------
!   BUILD_NOFIX_QM               Create selection of QM atoms and NOFIX from .dynn file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module NOFIX_QM

    use initialization
    use atoms, only             : NATOMS
    use atom_manipulation, only : ATOM_SELECTION

    contains

    !  BUILD_NOFIX_QM  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine build_nofix_qm(filename)

        !--------------------------------------------------------------
        ! Create selection of QM atoms and NOFIX from .dynn file
        !--------------------------------------------------------------

        character(len=*), intent(in)       :: filename

        character(len=32)                  :: option
        character(len=32)                  :: sele_mode
        integer                            :: io_stat

        character(len=64)                  :: subsystem(1)
        integer                            :: residue(1)
        character(len=64)                  :: atom(1)

        logical, allocatable               :: sele(:)
        logical, allocatable               :: subsystems(:)
        logical, allocatable               :: residues(:)
        logical, allocatable               :: atoms(:)

        CALL read_file(filename, .false.)

        ! initializate qm selection and nofix selection
        if (.not. allocated(qm_sele)) allocate(qm_sele(1:NATOMS))
        if (.not. allocated(nofix_sele)) allocate(nofix_sele(1:NATOMS))

        allocate( sele(1:NATOMS), &
                  subsystems(1:NATOMS), &
                  residues(1:NATOMS), &
                  atoms(1:NATOMS) )

        ! read selection from .dynn file, ignore not selection statements
        write(*,fmt='(/,A,A)') "QM and NOFIX read from ", trim(filename)
        open(unit=200, file=filename, status='old', action='read', form='formatted')
        do
            read(200, *, iostat=io_stat) option
            if (io_stat/=0) EXIT
            if (option(1:1) == '!' .or. option(1:1) == '#') CYCLE
            select case (trim(option))
                case ('QM', 'NOFIX')
                    sele_mode = trim(option)
                    sele = .false.
                    ! read selection block until same flag
                    do
                        read(200, *, iostat=io_stat) option
                        if (io_stat/=0 .or. option==sele_mode) EXIT
                        if (option(1:1) == '!' .or. option(1:1) == '#') CYCLE
                        backspace(200)
                        select case(option)
                            case('S')
                                read(200,*) option, subsystem(1)
                                subsystems = ATOM_SELECTION( subsystem=subsystem )
                                ! select whole subsystem, later removed if residues specified
                                sele = sele .or. subsystems
                                residues = .false.
                            case('R')
                                read(200,*) option, residue(1)
                                residues = ATOM_SELECTION( subsystem=subsystem, &
                                                           residue_number=residue )
                                ! update subsystems selection and use it as
                                ! negative -> remove from subsystems the residue
                                subsystems = subsystems .and. .not. residues
                                ! update total selection only including only
                                ! the falses in subsystems and the new residue
                                sele = sele .and. .not. subsystems .or. residues
                                atoms = .false.
                            case('A')
                                read(200,*) option, atom(1)
                                atoms = ATOM_SELECTION( subsystem=subsystem, &
                                                        residue_number=residue, &
                                                        atom_name=atom )
                                ! same procedure as before but with atoms
                                residues = residues .and. .not. atoms
                                sele = sele .and. .not. residues .or. atoms
                            case DEFAULT
                                read(100,*) option
                                write(*,fmt='(A,A)') "Omitting unknown option: ", option
                        end select
                    end do
                    ! assign built array to corresponding selection
                    select case(sele_mode)
                        case ('QM')
                            qm_sele(:) = sele(:)
                        case ('NOFIX')
                            nofix_sele(:) = sele(:)
                    end select
            end select
        end do
        close(unit=200)

    end subroutine

end module
