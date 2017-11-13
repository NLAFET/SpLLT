module spllt_ptg_mod

contains

  subroutine spllt_ptg_factorize(adata, val, cntl, fdata, info)
    use spllt_mod
    ! use hsl_ma87_double
    use spllt_error_mod
    use spllt_kernels_mod, only : spllt_activate_node, spllt_init_node, spllt_init_blk
!    use spllt_factorization_mod
#if defined(SPLLT_USE_PARSEC)
    use iso_c_binding
    use parsec_f08_interfaces
    use spllt_parsec_mod
    use spllt_parsec_factorization_mod
#endif
    implicit none

    type(spllt_adata_type), intent(in) :: adata
    real(wp), target, intent(in) :: val(:) ! matrix values
    ! type(spllt_keep), target, intent(inout) :: keep 
    type(spllt_cntl)      :: cntl
    type(spllt_fdata_type), target, intent(inout) :: fdata
    type(spllt_info), intent(out) :: info 

    integer :: n ! matrix order
    integer :: snode, num_nodes

    integer :: st ! stat parameter

#if defined(SPLLT_USE_PARSEC)
    ! PaRSEC 
    type(parsec_taskpool_t)            :: fac_tp
    type(c_ptr)                     :: bc_c, nodes_c, diags_c, fdata_c, val_c
    integer(c_int)                  :: nbc, nval
    integer                         :: start_setup_t, stop_setup_t, rate_setup_t
#endif

    integer(long) :: id

    ! write(*,'("[spllt_ptg_factorize]")')

    num_nodes = adata%nnodes
    n = adata%n

    ! write(*,'("[spllt_ptg_factorize] num_nodes: ", i8)')num_nodes
    ! write(*,'("[spllt_ptg_factorize]         n: ", i8)')n

    deallocate (fdata%lfact,stat=st)
    allocate (fdata%lfact(fdata%nbcol),stat=st)
    if(st.ne.0) go to 10

    do snode = 1, num_nodes ! loop over nodes
       ! activate node: allocate factors
       call spllt_activate_node(snode, fdata, adata)
       ! init node 
       ! TODO parallelize node init
       ! call spllt_init_node(snode, val, fdata)
    end do

    ! do id = 1, fdata%final_blk
    !    call spllt_init_blk(id, val, fdata)
    ! end do

#if defined(SPLLT_USE_PARSEC)

    nodes_c = c_loc(fdata%nodes(1))

    bc_c = c_loc(fdata%bc(1))
    nbc = size(fdata%bc,1)
    diags_c = c_loc(fdata%diags(1))
    fdata_c = c_loc(fdata)
    val_c = c_loc(val(1))
    nval = size(val, 1)

    call system_clock(start_setup_t, rate_setup_t)

    ! initialize block data descriptor

    fdata%ddesc = spllt_alloc_blk_desc()

#if defined(SPLLT_USE_MPI)
    call spllt_parsec_blk_data_init(fdata%ddesc, bc_c, nbc, nds, rank)
#else
    call data_init(fdata%ddesc, bc_c, nbc, nds, rank)
#endif    

    fac_tp = spllt_parsec_factorize(fdata%ddesc, nodes_c, num_nodes,bc_c, nbc, diags_c, fdata%nbcol, &
         & cntl%min_width_blas, fdata%maxmn, val_c, nval, fdata_c)
    
    ! add factorization DAG to PaRSEC
    call parsec_enqueue(ctx, fac_tp)
    call system_clock(stop_setup_t)
    ! write(*,'("[>] [spllt_ptg_factorize]   setup time: ", es10.3, " s")') (stop_setup_t - start_setup_t)/real(rate_setup_t)

    ! start factorization
    call parsec_context_start(ctx)
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
