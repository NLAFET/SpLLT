module spllt_ptg_mod

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Perform the task-based Cholesky factoization using a PTG
  !> model. Note that this call is asynchronous.
  !>
  !> @param akeep Symbolic factorization data.
  !> @param fkeep Factorization data.
  subroutine spllt_ptg_factorize(akeep, fkeep, options, val, info)
    use spllt_error_mod
    use spllt_kernels_mod, only : spllt_activate_node, &
         spllt_init_node, spllt_init_blk
#if defined(SPLLT_USE_PARSEC)
    use iso_c_binding
    use parsec_f08_interfaces
    use spllt_parsec_mod
    use spllt_parsec_factorization_mod
#endif
    implicit none

    type(spllt_akeep), intent(in)            :: akeep
    type(spllt_fkeep), target, intent(inout) :: fkeep
    type(spllt_options), intent(in)          :: options ! User-supplied options
    real(wp), target, intent(in)             :: val(:) ! Matrix values
    type(spllt_inform), intent(out)          :: info 

    integer :: n ! matrix order
    integer :: snode
    integer :: nnodes ! Number of nodes in the assembly tree
    integer :: st ! stat parameter

#if defined(SPLLT_USE_PARSEC)
    ! Parsec specific variables 
    type(parsec_taskpool_t)         :: fac_tp
    type(c_ptr)                     :: bc_c, nodes_c, diags_c, fkeep_c, val_c
    integer(c_int)                  :: nbc, nval
    integer                         :: start_setup_t, stop_setup_t, rate_setup_t
#endif

    integer(long) :: id

    nnodes = akeep%nnodes
    n = akeep%n

    deallocate (fkeep%lfact,stat=st)
    allocate (fkeep%lfact(fkeep%nbcol),stat=st)
    if(st.ne.0) go to 10

    ! Allocate factors
    do snode = 1, nnodes ! loop over nodes
       call spllt_activate_node(snode, fkeep)
    end do

#if defined(SPLLT_USE_PARSEC)

    nodes_c = c_loc(fkeep%nodes(1))

    bc_c = c_loc(fkeep%bc(1))
    nbc = size(fkeep%bc,1)
    diags_c = c_loc(fkeep%diags(1))
    fkeep_c = c_loc(fkeep)
    val_c = c_loc(val(1))
    nval = size(val, 1)

    call system_clock(start_setup_t, rate_setup_t)

    ! initialize block data descriptor

    fkeep%ddesc = spllt_alloc_blk_desc()

#if defined(SPLLT_USE_MPI)
    call spllt_parsec_blk_data_init(fkeep%ddesc, bc_c, nbc, nds, rank)
#else
    call data_init(fkeep%ddesc, bc_c, nbc, nds, rank)
#endif    

    ! Create Parsec taskpool for factorization.
    fac_tp = spllt_parsec_factorize(&
         fkeep%ddesc, nodes_c, nnodes, bc_c, nbc, diags_c, fkeep%nbcol, &
         & options%min_width_blas, fkeep%maxmn, val_c, nval, fkeep_c)
    
    ! Add factorization taskpool to Parsec runtime.
    call parsec_enqueue(ctx, fac_tp)
    call system_clock(stop_setup_t)
    ! write(*,'("[>] [spllt_ptg_factorize]   setup time: ", es10.3, " s")') (stop_setup_t - start_setup_t)/real(rate_setup_t)

    ! Start factorization
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
