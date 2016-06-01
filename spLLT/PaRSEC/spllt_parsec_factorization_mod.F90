module spllt_parsec_factorization_mod

  interface
     function spllt_parsec_factorize(nodes, num_nodes, bcs, nbc, diags, ndiag, &
          & min_width_blas) bind(C)
       use iso_c_binding
       use dague_f08_interfaces
       type(c_ptr), value      :: nodes, bcs, diags
       integer(c_int), value   :: num_nodes, nbc, ndiag
       integer(c_int), value   :: min_width_blas
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

  function get_node(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    use hsl_ma87_double
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int)              :: get_node

    type(spllt_bc_type), pointer :: bc(:) ! blocks
    
    call c_f_pointer(bcs_c, bc,(/id/))    

    get_node = bc(id)%blk%node 
    
  end function get_node

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

  ! get number of internode updates to be performed on block id 
  function get_dep_out_count(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int)             :: get_dep_out_count

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_out_count = 0
    
    if (associated(bc(id)%dep_out)) then
       get_dep_out_count = size(bc(id)%dep_out)
    end if

  end function get_dep_out_count

  function get_dep_out_p(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(c_int)             :: get_dep_out_p

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_out_p = 0

    if (associated(bc(id)%dep_out)) then
       get_dep_out_p = bc(id)%dep_out(i)%p
    end if

  end function get_dep_out_p

  ! get number of internode updates to be performed on block id 
  function get_dep_out_id_ij(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(long)              :: get_dep_out_id_ij

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_out_id_ij = 0

    if (associated(bc(id)%dep_out)) then
       get_dep_out_id_ij = bc(id)%dep_out(i)%id_ij
    end if

  end function get_dep_out_id_ij

  function get_dep_out_flow(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer                    :: get_dep_out_flow

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_out_flow = 0

    if (associated(bc(id)%dep_out)) then
       get_dep_out_flow = bc(id)%dep_out(i)%flow
    end if

  end function get_dep_out_flow

  function get_dep_in_csrc(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(c_int)             :: get_dep_in_csrc

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_in_csrc = 0
    ! write(*,*) 'i: ', i
    if (associated(bc(id)%dep_in)) then
        get_dep_in_csrc = bc(id)%dep_in(i)%csrc
    end if

  end function get_dep_in_csrc

  function get_dep_in_rsrc(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(c_int)             :: get_dep_in_rsrc

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_in_rsrc = 0
    ! write(*,*) 'i: ', i
    if (associated(bc(id)%dep_in)) then
        get_dep_in_rsrc = bc(id)%dep_in(i)%rsrc
    end if

  end function get_dep_in_rsrc

  function get_dep_in_sync(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(c_int)             :: get_dep_in_sync

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_in_sync = 0
    ! write(*,*) 'i: ', i
    if (associated(bc(id)%dep_in)) then
        if (bc(id)%dep_in(i)%sync) get_dep_in_sync = 1
    end if

  end function get_dep_in_sync
  
  function get_dep_in_p1(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(c_int)             :: get_dep_in_p1

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_in_p1 = 0
    ! write(*,*) 'i: ', i
    if (associated(bc(id)%dep_in)) then
       get_dep_in_p1 = bc(id)%dep_in(i)%p1
    end if

  end function get_dep_in_p1

  function get_dep_in_id_jk(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(long)              :: get_dep_in_id_jk

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_in_id_jk = 0

    if (associated(bc(id)%dep_in)) then
       get_dep_in_id_jk = bc(id)%dep_in(i)%id_jk
    end if

  end function get_dep_in_id_jk

  function get_dep_in_p2(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(c_int)             :: get_dep_in_p2

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_in_p2 = 0

    if (associated(bc(id)%dep_in)) then
       get_dep_in_p2 = bc(id)%dep_in(i)%p2
    end if

  end function get_dep_in_p2

  function get_dep_in_id_ik(bcs_c, id, i) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int), value      :: i
    integer(long)              :: get_dep_in_id_ik

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_dep_in_id_ik = 0

    if (associated(bc(id)%dep_in)) then
       get_dep_in_id_ik = bc(id)%dep_in(i)%id_ik
    end if

  end function get_dep_in_id_ik

  ! get number of internode updates to be performed on block id 
  function get_upd_count(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: id
    integer(c_int)             :: get_upd_count

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))    

    get_upd_count = 0
    
    if (associated(bc(id)%dep_in)) then
       get_upd_count = size(bc(id)%dep_in)
    end if

  end function get_upd_count

  function get_next_dblk(bcs_c, diag) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value       :: diag
    integer(c_int)             :: get_next_dblk

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/diag/))    
    
    get_next_dblk = bc(diag)%blk%last_blk + 1 
      
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
    integer(long), value             :: id
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
    integer(long), value      :: id
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
    integer(long), value      :: id
    integer(c_int)             :: get_blk_n

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    call c_f_pointer(bcs_c, bc,(/id/))

    get_blk_n = bc(id)%blk%blkn

  end function get_blk_n

  function get_blk(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value      :: id
    type(c_ptr)                :: get_blk

    type(spllt_bc_type), pointer :: bc(:) ! blocks

    ! write(*,*)'[get_blk_ptr] id: ', id
    get_blk = c_null_ptr

    call c_f_pointer(bcs_c, bc,(/id/))
    
    get_blk = c_loc(bc(id)%blk)

    ! write(*,*)'[get_blk_ptr] id:', id, ', ptr:', get_blk_ptr
    
  end function get_blk

  function get_blk_ptr(bcs_c, id) bind(C)
    use iso_c_binding
    use spllt_mod
    implicit none

    type(c_ptr), value, target :: bcs_c
    integer(long), value      :: id
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
    integer(long), value             :: id
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
    integer(c_int) :: bcol
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
    integer(c_int) :: bcol
    integer(c_size_t) :: get_lcol_sze

    type(lfactor), dimension(:), pointer :: lfact    

    call c_f_pointer(lfact_c, lfact,(/bcol/))
    
    get_lcol_sze = size(lfact(bcol)%lcol) * wp

  end function get_lcol_sze

end module spllt_parsec_factorization_mod
