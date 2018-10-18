!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
module spllt_error_mod
  use spllt_data_mod
  implicit none

contains
  
  subroutine spllt_print_err(iflag, context, st)
    ! use hsl_ma87_double
    use spllt_data_mod
    implicit none

    integer, intent(in) :: iflag
    character (len=*), optional, intent(in) :: context
    integer, intent(in), optional :: st

    select case(iflag)
    case(spllt_error_allocation)
       write(*,*) 'allocation error'
    case default
       write(*,*) 'unknown error'
    end select

    return
  end subroutine spllt_print_err
  
end module spllt_error_mod
