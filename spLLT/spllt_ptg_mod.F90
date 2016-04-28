module spllt_ptg_mod

contains

  subroutine spllt_ptg_factorize(adata, keep, cntl, fdata, info)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod, only : spllt_activate_node
!    use spllt_factorization_mod
#if defined(SPLLT_USE_PARSEC)
    use dague_f08_interfaces
    use spllt_parsec_factorization_mod
#endif
    implicit none

    type(spllt_adata_type), intent(in) :: adata
    type(MA87_keep), target, intent(inout) :: keep 
    type(spllt_cntl)      :: cntl
    type(spllt_data_type), intent(inout) :: fdata
    type(MA87_info), intent(out) :: info 

    integer :: n ! matrix order
    integer :: snode, num_nodes
    integer, dimension(:), allocatable ::  map ! allocated to have size n.

    integer :: st ! stat parameter

#if defined(SPLLT_USE_PARSEC)
    ! PaRSEC 
    type(dague_handle_t)            :: fac_hdl
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
       ! activate node: allocate factors, register handles
       call spllt_activate_node(snode, keep, fdata)
    end do

#if defined(SPLLT_USE_PARSEC)
    fac_hdl = spllt_parsec_factorize(num_nodes)

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
