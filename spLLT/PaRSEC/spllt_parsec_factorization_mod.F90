module spllt_parsec_factorization_mod

  interface
     function spllt_parsec_factorize(nodes, num_nodes, bcs, nbc) bind(C)
       use iso_c_binding
       use dague_f08_interfaces
       type(c_ptr), value      :: nodes, bcs
       integer(c_int), value   :: num_nodes, nbc
       type(dague_handle_t)    :: spllt_parsec_factorize
     end function spllt_parsec_factorize

     ! function spllt_parsec_factorize(num_nodes) bind(C)
     !   use iso_c_binding
     !   use dague_f08_interfaces
     !   integer(c_int), value   :: num_nodes
     !   type(dague_handle_t)    :: spllt_parsec_factorize
     ! end function spllt_parsec_factorize
  end interface

contains

  function get_blk_sa(nodes_c, nnodes, bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    use hsl_ma87_double
    implicit none

    type(c_ptr), value, target :: nodes_c, bcs_c
    integer(c_int)             :: nnodes
    integer(long), value       :: id
    integer(long)              :: get_blk_sa

    type(spllt_bc_type), pointer :: bc(:) ! blocks
    type(node_type), pointer :: node, nodes(:)
    
    call c_f_pointer(nodes_c, nodes, (/nnodes/))    
    call c_f_pointer(bcs_c, bc,(/id/))    

    node => nodes(bc(id)%blk%node)
    get_blk_sa = node%blk_sa 
    
  end function get_blk_sa

  function get_nc(nodes_c, nnodes, bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    use hsl_ma87_double
    implicit none

    type(c_ptr), value, target :: nodes_c, bcs_c
    integer(c_int)             :: nnodes
    integer(long), value       :: id
    integer(long)              :: get_nc

    type(spllt_bc_type), pointer :: bc(:) ! blocks
    type(node_type), pointer :: node, nodes(:)
    integer :: sa, en, numcol, s_nb
    
    call c_f_pointer(nodes_c, nodes, (/nnodes/))    
    call c_f_pointer(bcs_c, bc,(/id/))    

    node => nodes(bc(id)%blk%node)
    sa = node%sa
    en = node%en
    numcol = en - sa + 1
    s_nb = node%nb
    get_nc = (numcol-1) / s_nb + 1
    
  end function get_nc

  function get_dest_blk_id(bcs_c, nbc, id_jk, id_ik) bind(C)
    use iso_c_binding
    use spllt_mod
    use hsl_ma87_double
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(c_int), value      :: nbc
    integer(long), value       :: id_jk, id_ik
    integer(long)              :: get_dest_blk_id

    type(spllt_bc_type), pointer :: bc(:) ! blocks
    
    call c_f_pointer(bcs_c, bc,(/nbc/))    
    
    get_dest_blk_id = get_dest_block(bc(id_jk)%blk, bc(id_ik)%blk)

  end function get_dest_blk_id

  function get_next_dblk(bcs_c, diag) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: diag
    integer(c_int)             :: get_prev_dblk

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/diag/))    
    
    get_prev_dblk = bc(diag-1)%blk%dblk 
      
  end function get_next_dblk

  function get_prev_dblk(bcs_c, diag) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: diag
    integer(c_int)             :: get_prev_dblk

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/diag/))    
    
    get_prev_dblk = 0

    if (diag .gt. 1) then
       get_prev_dblk = bc(diag-1)%blk%dblk 
    end if
       
  end function get_prev_dblk

  function get_last_blk(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int)             :: get_last_blk

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_last_blk = bc(id)%blk%last_blk 
    
  end function get_last_blk

  function get_bcol(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int)             :: get_bcol

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))

    get_bcol = bc(id)%blk%bcol 
    
  end function get_bcol

  function get_num_subdiag(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int)             :: get_num_subdiag

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_num_subdiag = 0
    ! check that blk id is diag
    if (bc(id)%blk%dblk .eq. id) then
       get_num_subdiag = bc(id)%blk%last_blk - id 
    end if

  end function get_num_subdiag

  function is_diag(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer, value             :: id
    integer(c_int)             :: is_diag

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    
    
    is_diag = 0

    if (bc(id)%blk%dblk .eq. id) is_diag = 1

  end function is_diag

  function get_blk_m(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(c_int), value      :: id
    integer(c_int)             :: get_blk_m

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))

    get_blk_m = bc(id)%blk%blkm

  end function get_blk_m

  function get_blk_n(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(c_int), value      :: id
    integer(c_int)             :: get_blk_n

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))

    get_blk_n = bc(id)%blk%blkn

  end function get_blk_n

  function get_blk_ptr(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(c_int), value      :: id
    type(c_ptr)                :: get_blk_ptr

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    ! write(*,*)'[get_blk_ptr] id: ', id
    get_blk_ptr = c_null_ptr

    call c_f_pointer(bcs_c, bc,(/id/))
    
    get_blk_ptr = c_loc(bc(id)%c(1))

    ! write(*,*)'[get_blk_ptr] id:', id, ', ptr:', get_blk_ptr
    
  end function get_blk_ptr

  function get_blk_sze(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer, value             :: id
    integer(c_size_t)          :: get_blk_sze

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
