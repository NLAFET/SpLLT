module spllt_ptg_mod

contains

  subroutine spllt_ptg_factorize(adata, val, keep, cntl, fdata, info)
    use spllt_mod
    use hsl_ma87_double
    use spllt_error_mod
    use spllt_kernels_mod, only : spllt_activate_node, spllt_init_node, spllt_init_blk
!    use spllt_factorization_mod
#if defined(SPLLT_USE_PARSEC)
    use iso_c_binding
    use dague_f08_interfaces
    use spllt_parsec_mod
    use spllt_parsec_factorization_mod
#endif
    implicit none

    type(spllt_adata_type), intent(in) :: adata
    real(wp), target, intent(in) :: val(:) ! matrix values
    type(MA87_keep), target, intent(inout) :: keep 
    type(spllt_cntl)      :: cntl
    type(spllt_data_type), target, intent(inout) :: fdata
    type(MA87_info), intent(out) :: info 

    integer :: n ! matrix order
    integer :: snode, num_nodes

    integer :: st ! stat parameter

#if defined(SPLLT_USE_PARSEC)
    ! PaRSEC 
    type(dague_handle_t)            :: fac_hdl
    type(c_ptr)                     :: bc_c, nodes_c, diags_c, keep_c, val_c
    integer(c_int)                  :: nbc, nval
    integer                         :: start_setup_t, stop_setup_t, rate_setup_t
#endif

    integer(long) :: id

    write(*,'("[spllt_ptg_factorize]")')

    num_nodes = adata%nnodes
    n = adata%n

    write(*,'("[spllt_ptg_factorize] num_nodes: ", i8)')num_nodes
    write(*,'("[spllt_ptg_factorize]         n: ", i8)')n

    deallocate (keep%lfact,stat=st)
    allocate (keep%lfact(keep%nbcol),stat=st)
    if(st.ne.0) go to 10

    do snode = 1, num_nodes ! loop over nodes
       ! activate node: allocate factors
       call spllt_activate_node(snode, keep, fdata)
       ! init node 
       ! TODO parallelize node init
       ! call spllt_init_node(snode, val, keep)
    end do

    ! do id = 1, keep%final_blk
    !    call spllt_init_blk(id, val, keep)
    ! end do

#if defined(SPLLT_USE_PARSEC)

    nodes_c = c_loc(keep%nodes(1))

    bc_c = c_loc(fdata%bc(1))
    nbc = size(fdata%bc,1)
    diags_c = c_loc(fdata%diags(1))
    keep_c = c_loc(keep)
    val_c = c_loc(val(1))
    nval = size(val, 1)

    call system_clock(start_setup_t, rate_setup_t)

    ! initialize block data descriptor

    fdata%ddesc = spllt_alloc_blk_desc()

    call spllt_parsec_blk_data_init(fdata%ddesc, bc_c, nbc, nds, rank) 
    
    fac_hdl = spllt_parsec_factorize(fdata%ddesc, nodes_c, num_nodes,bc_c, nbc, diags_c, keep%nbcol, &
         & cntl%min_width_blas, keep%maxmn, val_c, nval, keep_c)
    
    ! add factorization DAG to PaRSEC
    call dague_enqueue(ctx, fac_hdl)
    call system_clock(stop_setup_t)
    write(*,'("[>] [spllt_ptg_factorize]   setup time: ", es10.3, " s")') (stop_setup_t - start_setup_t)/real(rate_setup_t)

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
