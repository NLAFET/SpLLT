module spllt_ptg_mod

contains

  subroutine spllt_ptg_factorize(adata, val, keep, cntl, fdata, info)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod, only : spllt_activate_node, spllt_init_node
!    use spllt_factorization_mod
#if defined(SPLLT_USE_PARSEC)
    use iso_c_binding
    use dague_f08_interfaces
    use spllt_parsec_factorization_mod
#endif
    implicit none

    type(spllt_adata_type), intent(in) :: adata
    real(wp), intent(in) :: val(:) ! matrix values
    type(MA87_keep), target, intent(inout) :: keep 
    type(spllt_cntl)      :: cntl
    type(spllt_data_type), target, intent(inout) :: fdata
    type(MA87_info), intent(out) :: info 

    integer :: n ! matrix order
    integer :: snode, num_nodes
    integer, dimension(:), allocatable ::  map ! allocated to have size n.

    integer :: st ! stat parameter

#if defined(SPLLT_USE_PARSEC)
    ! PaRSEC 
    type(dague_handle_t)            :: fac_hdl
    type(c_ptr) :: bc_c, nodes_c, diags_c
    integer(c_int) :: nbc
#endif

    write(*,'("[spllt_ptg_factorize]")')

    num_nodes = adata%nnodes
    n = adata%n

    write(*,'("[spllt_ptg_factorize] num_nodes: ", i8)')num_nodes
    write(*,'("[spllt_ptg_factorize]         n: ", i8)')n

    deallocate (keep%lfact,stat=st)
    allocate (keep%lfact(keep%nbcol),stat=st)
    if(st.ne.0) go to 10

    allocate(map(n),stat=st)
    if(st.ne.0) go to 10

    do snode = 1, num_nodes ! loop over nodes
       ! activate node: allocate factors
       call spllt_activate_node(snode, keep, fdata)
       ! init node 
       ! TODO parallelize node init
       call spllt_init_node(snode, n, val, map, keep)
    end do

#if defined(SPLLT_USE_PARSEC)

    nodes_c = c_loc(keep%nodes(1))

    bc_c = c_loc(fdata%bc(1))
    nbc = size(fdata%bc,1)
    diags_c = c_loc(fdata%diags(1))

    fac_hdl = spllt_parsec_factorize(nodes_c, num_nodes,bc_c, nbc, diags_c, keep%nbcol, &
         & cntl%min_width_blas)

    ! add factorization DAG to PaRSEC
    call dague_enqueue(ctx, fac_hdl)

    ! start factorization
    call dague_context_start(ctx)
#endif

10  if(st.ne.0) then
       info%flag = spllt_error_allocation
       info%stat = st
       call spllt_print_err(info%flag, context='spllt_ptg_factorize',st=st)
       return
    endif

    return
  end subroutine spllt_ptg_factorize

end module spllt_ptg_mod
