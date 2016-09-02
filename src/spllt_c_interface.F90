module spllt_c_interface
  use iso_c_binding

contains

  function is_blk_diag(blk_c) bind(C)
    use iso_c_binding
    use hsl_ma87_double, only: block_type
    implicit none

    type(c_ptr), value, target :: blk_c
    integer(c_int)             :: is_blk_diag

    type(block_type), pointer :: blk ! blocks

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

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    is_diag = 0

    if (bc(id)%blk%dblk .eq. id) is_diag = 1

  end function is_diag
  
end module spllt_c_interface