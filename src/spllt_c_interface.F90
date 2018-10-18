!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
module spllt_c_interface
  use iso_c_binding

contains

  function is_blk_diag(blk_c) bind(C)
    use iso_c_binding
    use spllt_data_mod, only: spllt_block
    implicit none

    type(c_ptr), value, target :: blk_c
    integer(c_int)             :: is_blk_diag

    type(spllt_block), pointer :: blk ! blocks

    call c_f_pointer(blk_c, blk)

    is_blk_diag = 0

    if (blk%dblk .eq. blk%id) is_blk_diag = 1

  end function is_blk_diag
  
  function is_diag(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_data_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value             :: id
    integer(c_int)             :: is_diag

    type(spllt_block), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    is_diag = 0

    if (bc(id)%dblk .eq. id) is_diag = 1

  end function is_diag
  
end module spllt_c_interface
