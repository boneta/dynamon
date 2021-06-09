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
        character(len=512)                 :: sele_name
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
            if (option=='SELECTION') then
                backspace(200)
                read(200, *) option, sele_name
                sele = .false.
                ! read selection block until 'SELECTION'
                do
                    read(200, *, iostat=io_stat) option
                    if (io_stat/=0 .or. option=='SELECTION') EXIT
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
                select case(sele_name)
                    case ('QM')
                        qm_sele(:) = sele(:)
                    case ('NOFIX')
                        nofix_sele(:) = sele(:)
                    case DEFAULT
                        write(*,fmt='(A,A)') "Omitting selection: ", trim(sele_name)
                end select
            end if
        end do
        close(unit=200)

        ! warning message if not QM/NOFIX read
        if (.NOT. ANY(qm_sele)) CALL dynn_log(2, "Missing crucial selection", "QM")
        if (.NOT. ANY(nofix_sele)) CALL dynn_log(2, "Missing crucial selection", "NOFIX")

    end subroutine

end module
