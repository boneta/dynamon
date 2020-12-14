subroutine flush( unit )
        include "fiosetup_.h"
        implicit none
        integer, intent(in) :: unit
        integer :: ires, fiosetup_
        ires = fiosetup_( unit, IO_CMD_FLUSH_AFTER_WRITE, IO_ARG_FLUSH_YES )
end subroutine flush
