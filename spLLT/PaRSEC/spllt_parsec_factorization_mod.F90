module spllt_parsec_factorization_mod

  interface
     function spllt_parsec_factorize(num_nodes) bind(C)
       use iso_c_binding
       use dague_f08_interfaces
       integer(c_int), value   :: num_nodes
       type(dague_handle_t)    :: spllt_parsec_factorize
     end function spllt_parsec_factorize
  end interface

contains

  function get_blk_ptr(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), target :: bcs_c
    integer :: id
    type(c_ptr) :: get_blk_ptr

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))
    
    get_blk_ptr = c_loc(bc(id)%c(1))
    
  end function get_blk_ptr

  function get_blk_sze(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), target :: bcs_c
    integer :: id
    integer(c_size_t) :: get_blk_sze

    type(spllt_bc_type), pointer :: bc(:) ! blocks
    integer m, n

    call c_f_pointer(bcs_c, bc,(/id/))
    
    m = bc(id)%blk%blkm
    n = bc(id)%blk%blkn

    get_blk_sze = m*n * wp
    
  end function get_blk_sze

  function get_lcol_ptr(lfact_c, bcol) bind(C)
    use iso_c_binding
    use hsl_ma87_double
    implicit none

    type(c_ptr), target :: lfact_c
    integer :: bcol
    type(c_ptr) :: get_lcol_ptr

    type(lfactor), dimension(:), pointer :: lfact    

    call c_f_pointer(lfact_c, lfact,(/bcol/))
    
    get_lcol_ptr = c_loc(lfact(bcol)%lcol(1))

  end function get_lcol_ptr

  function get_lcol_sze(lfact_c, bcol) bind(C)
    use iso_c_binding
    use spllt_mod
    use hsl_ma87_double
    implicit none

    type(c_ptr), target :: lfact_c
    integer :: bcol
    integer(c_size_t) :: get_lcol_sze

    type(lfactor), dimension(:), pointer :: lfact    

    call c_f_pointer(lfact_c, lfact,(/bcol/))
    
    get_lcol_sze = size(lfact(bcol)%lcol) * wp

  end function get_lcol_sze

end module spllt_parsec_factorization_mod
