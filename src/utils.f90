!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                        UTILS & MISCELLANY                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Subroutines
!  -----------
!   APPEND_1D                    Increment by one a 1D array and append
!                                a value to the end if present
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module UTILS

    implicit none

    interface append_1d
        !--------------------------------------------------------------
        ! Increment 1D array by one and append value if present
        !--------------------------------------------------------------
        module procedure append_1d_integer
        module procedure append_1d_real4
        module procedure append_1d_real8
        module procedure append_1d_logical
        module procedure append_1d_character
    end interface

    private
    public         :: append_1d

    contains

    subroutine append_1d_integer(array, val)

        integer, allocatable, intent(inout) :: array(:)
        integer, intent(in), optional       :: val

        integer                             :: n
        integer, allocatable                :: array_tmp(:)

        n = SIZE(array, 1) + 1
        allocate(array_tmp(n))
        array_tmp = 0
        array_tmp(:) = array(:)
        CALL MOVE_ALLOC(array_tmp, array)
        if (PRESENT(val)) array(n) = val

    end subroutine

    subroutine append_1d_real4(array, val)

        real(4), allocatable, intent(inout) :: array(:)
        real(4), intent(in), optional       :: val

        integer                             :: n
        real(4), allocatable                :: array_tmp(:)

        n = SIZE(array, 1) + 1
        allocate(array_tmp(n))
        array_tmp = 0.E0
        array_tmp(:) = array(:)
        CALL MOVE_ALLOC(array_tmp, array)
        if (PRESENT(val)) array(n) = val

    end subroutine

    subroutine append_1d_real8(array, val)

        real(8), allocatable, intent(inout) :: array(:)
        real(8), intent(in), optional       :: val

        integer                             :: n
        real(8), allocatable                :: array_tmp(:)

        n = SIZE(array, 1) + 1
        allocate(array_tmp(n))
        array_tmp = 0.D0
        array_tmp(:) = array(:)
        CALL MOVE_ALLOC(array_tmp, array)
        if (PRESENT(val)) array(n) = val

    end subroutine

    subroutine append_1d_logical(array, val)

        logical, allocatable, intent(inout) :: array(:)
        logical, intent(in), optional       :: val

        integer                             :: n
        logical, allocatable                :: array_tmp(:)

        n = SIZE(array, 1) + 1
        allocate(array_tmp(n))
        array_tmp = .false.
        array_tmp(:) = array(:)
        CALL MOVE_ALLOC(array_tmp, array)
        if (PRESENT(val)) array(n) = val

    end subroutine

    subroutine append_1d_character(array, val)

        character(len=*), allocatable, intent(inout) :: array(:)
        character(len=*), intent(in), optional       :: val

        integer                                      :: n
        character(len=LEN(array)), allocatable       :: array_tmp(:)

        n = SIZE(array, 1) + 1
        allocate(array_tmp(n))
        array_tmp = ''
        array_tmp(:) = array(:)
        CALL MOVE_ALLOC(array_tmp, array)
        if (PRESENT(val)) array(n) = val

    end subroutine

end module
