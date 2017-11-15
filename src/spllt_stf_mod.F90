module spllt_stf_mod

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Perform the task-based Cholesky factoization using a STF
  !> model. Note that this call is asynchronous.
  !>
  !> @param adata Symbolic factorization data.
  !> @param fdata Factorization data.
  !> @param options User-supplied options.
  !> @param val Matrix values
  !> @param info 
  subroutine spllt_stf_factorize(adata, fdata, options, val, info)
    use spllt_data_mod
    use spllt_error_mod
    use spllt_factorization_task_mod
#if defined(SPLLT_USE_STARPU) 
    use iso_c_binding
    use starpu_f_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#endif
    use spllt_kernels_mod
    use spllt_factorization_mod
    implicit none

    type(spllt_adata), target, intent(in)    :: adata ! Symbolic factorization data.
    type(spllt_fdata), target, intent(inout) :: fdata  ! Factorization data.
    type(spllt_options), intent(in)          :: options ! User-supplied options
    real(wp), intent(in)                     :: val(:) ! matrix values
    type(spllt_inform), intent(out)          :: info

    ! used in copying entries of user's matrix a into factor storage 
    ! (keep%fact).
    integer, dimension(:), pointer ::  map
    ! shortcuts
    type(spllt_node), pointer     :: node ! node in the atree    
    type(spllt_block), pointer :: bc_kk, bc_ik, bc_jk, bc_ij
    type(spllt_block), dimension(:), pointer :: blocks ! block info. 
    type(spllt_node), dimension(:), pointer :: nodes 

    ! local scalars
    integer :: n ! Order of A
    integer(long) :: blk, blk1, blk2 ! block identity
    integer(long) :: dblk ! diagonal block within block column
    integer :: en ! holds keep%nodes(snode)%en
    integer :: sa ! holds keep%nodes(snode)%sa
    ! integer :: s_nb ! set to block size of snode (keep%nodes(snode)%nb)
    integer :: nc, nr ! number of block column/row
    integer :: snode, num_nodes
    integer :: st = 0 ! stat parameter
    integer :: numrow, numcol
    integer :: ii, jj, kk
    integer :: prio

#if defined(SPLLT_USE_OMP)
    integer :: nt
#endif
    ! timing
    integer :: start_t, stop_t, rate_t
    integer :: stf_start_t, stf_stop_t, stf_rate_t
    integer :: start_nosub_t, rate_nosub_t
    integer :: start_setup_t, stop_setup_t, rate_setup_t
    integer :: subtree_start_t, subtree_stop_t, subtree_rate_t

    ! subtree factor variable
    type(spllt_block) :: buffer

    ! start measuring setup time
    call system_clock(start_setup_t, rate_setup_t)

    n = adata%n
    
    ! shortcut
    blocks => fdata%bc
    nodes  => fdata%nodes
    num_nodes = fdata%info%num_nodes

    ! init facto, allocate map array, factor blocks
    call spllt_factorization_init(fdata, map)    

    if (options%prune_tree) allocate(buffer%c(1))

    ! start measuring time for submitting tasks
    call system_clock(stf_start_t, stf_rate_t)

    do snode = 1, num_nodes ! loop over nodes
#if defined(SPLLT_USE_STARPU)
       ! TODO put in activate routine
       call starpu_f_void_data_register(fdata%nodes(snode)%hdl)
#if defined(SPLLT_USE_NESTED_STF)
       call starpu_f_void_data_register(fdata%nodes(snode)%hdl2)
#endif
#endif
       ! activate node: allocate factors, register handles
       call spllt_activate_node(snode, fdata, adata)
    end do

    call system_clock(stop_setup_t)
    write(*,'("[>] [spllt_stf_factorize]   setup and activate nodes time: ", es10.3, " s")') &
         & (stop_setup_t - start_setup_t)/real(rate_setup_t)

#if defined(SPLLT_USE_STARPU) && defined(SPLLT_STARPU_NOSUB)
    call starpu_f_pause()
#endif

    ! #if defined(SPLLT_USE_STARPU)
    !     call starpu_f_task_wait_for_all()
    ! #endif    
    ! call system_clock(start_cpya2l_t, rate_cpya2l_t)
    ! call copy_a_to_l(n,num_nodes,val,map,keep)
    ! write(*,*)"num_nodes: ", num_nodes
    do snode = 1, num_nodes ! loop over nodes
       prio = 5 ! max priority 
       ! prio = huge(1)
       ! init node
       if (adata%small(snode) .ne. 0) cycle

       call spllt_init_node_task(fdata, fdata%nodes(snode), val, prio)
    end do
    ! call system_clock(stop_cpya2l_t)

#if defined(SPLLT_USE_OMP)
    ! @todo remove synchronization between init and facto phases.
    !$omp taskwait
#endif

    ! factorize nodes
    do snode = 1, num_nodes

       if (adata%small(snode) .lt. 0) cycle
       if (adata%small(snode) .eq. 1) then
        
          ! Factorize the subtree rooted at node snode. All
          ! contributions to ancestors above the root node are
          ! accumulated into a buffer that is scatered after the
          ! subtree factorization
          call spllt_subtree_factorize_apply(fdata, options, snode, val, map, buffer)

       else
          ! if (adata%small())
#if defined(SPLLT_USE_NESTED_STF)
          ! prio = 5 ! max priority 
          prio = -5 ! min priority 
          ! prio = huge(1)

          call spllt_factorize_apply_node_task(fdata, options, fdata%nodes(snode), prio)
#else
          call spllt_factorize_apply_node(fdata, options, fdata%nodes(snode), map)
#endif

       end if

    end do

    ! @todo unregister handles when using nested stf 
#ifndef SPLLT_USE_NESTED_STF
    ! deinit factorization
    ! clean up data structure, unregister data handles
    call spllt_factorization_fini(fdata, map, adata)
#endif

    call system_clock(stf_stop_t)
    write(*,'("[>] [spllt_stf_factorize] task insert time: ", es10.3, " s")') (stf_stop_t - stf_start_t)/real(stf_rate_t)

#if defined(SPLLT_USE_STARPU) && defined(SPLLT_STARPU_NOSUB)
    call system_clock(start_nosub_t, rate_nosub_t)
    call starpu_f_resume()
    ! wait for task completion
    call starpu_f_task_wait_for_all()
#endif

    call system_clock(stop_t)
    ! write(*,'("[>] [spllt_stf_factorize] time: ", es10.3, " s")') (stop_t - start_t)/real(rate_t)
#if defined(SPLLT_USE_STARPU) && defined(SPLLT_STARPU_NOSUB)
    write(*,'("[>] [spllt_stf_factorize] nosub time: ", es10.3, " s")') (stop_t - start_nosub_t)/real(rate_nosub_t)
#endif

10  if(st.ne.0) then
       info%flag = spllt_error_allocation
       info%stat = st
       call spllt_print_err(info%flag, context='spllt_stf_factorize',st=st)
       return
    endif

    return
  end subroutine spllt_stf_factorize

  ! left looking variant of the STF algorithm
  !   subroutine spllt_stf_ll_factorize(n, ptr, row, val, order, keep, control, info, fdata, options)
  !     use spllt_data_mod
  !     use spllt_error_mod
  !     use hsl_ma87_double
  !     use spllt_factorization_task_mod
  ! #if defined(SPLLT_USE_STARPU) 
  !     use iso_c_binding
  !     use starpu_f_mod
  ! #elif defined(SPLLT_USE_OMP)
  !     !$ use omp_lib
  ! #endif
  !     use spllt_kernels_mod
  !     implicit none

  !     integer, intent(in) :: n ! order of A
  !     integer, intent(in) :: row(:) ! row indices of lower triangular part
  !     integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
  !     real(wp), intent(in) :: val(:) ! matrix values
  !     integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
  !     ! since the analyse phase)

  !     type(MA87_keep), target, intent(inout) :: keep 
  !     type(MA87_control), intent(in) :: control 
  !     type(MA87_info), intent(out) :: info 

  !     type(spllt_fdata), target :: fdata
  !     type(spllt_options)      :: options

  !     type(block_type), dimension(:), pointer :: blocks ! block info. 
  !     type(spllt_node), dimension(:), pointer :: nodes 
  !     integer :: snode, num_nodes
  !     type(spllt_node), pointer     :: node, dnode ! node in the atree    

  !     integer, dimension(:), allocatable ::  tmpmap

  !     integer :: st ! stat parameter
  !     integer :: least_desc, d_num, d_sa, d_en, d_numcol, d_numrow, d_nb, d_nc
  !     ! logical :: map_done
  !     integer :: cptr
  !     integer :: cptr2
  !     integer :: cb, jb
  !     integer(long) :: dblk, blk 
  !     integer :: jlast
  !     integer :: stf_start_t, stf_stop_t, stf_rate_t
  !     integer :: ii, i, ilast, k
  !     type(spllt_block), pointer :: bc, d_bc_kk
  !     integer(long) :: d_dblk 
  !     integer :: kk

  !     integer :: prio
  !     integer :: sa, en
  !     integer :: numcol, numrow
  !     integer :: nc, nr
  !     integer :: s_nb
  !     integer(long) :: blk1, blk2 
  !     type(spllt_block), pointer :: bc_kk, bc_ik, bc_jk, bc_ij
  !     integer :: jj

  ! #if defined(SPLLT_USE_OMP)
  !     integer :: nt
  ! #endif

  !     ! shortcut
  !     blocks => keep%blocks
  !     nodes  => keep%nodes
  !     num_nodes = keep%info%num_nodes

  !     ! allocate L factors
  !     deallocate (keep%lfact,stat=st)
  !     allocate (keep%lfact(keep%nbcol),stat=st)
  !     if(st.ne.0) go to 9999

  !     allocate(tmpmap(n),stat=st)
  !     if(st.ne.0) go to 9999

  !     ! init facto    

  ! #if defined(SPLLT_USE_STARPU)

  !     ! register workspace handle
  !     call starpu_f_vector_data_register(fdata%workspace%hdl, -1, c_null_ptr, &
  !          & int(keep%maxmn*keep%maxmn, kind=c_int), int(wp,kind=c_size_t))
  ! #elif defined(SPLLT_USE_OMP)

  !     nt = 1
  ! !$  nt = omp_get_num_threads()
  !     allocate(fdata%workspace(0:nt-1))

  !     do i=0,nt-1
  !        allocate(fdata%workspace(i)%c(keep%maxmn*keep%maxmn), stat=st)
  !     end do
  ! #else
  !     allocate(fdata%workspace%c(keep%maxmn*keep%maxmn), stat=st)
  !     if(st.ne.0) goto 9999
  ! #endif

  !     call system_clock(stf_start_t, stf_rate_t)

  !     do snode = 1, num_nodes ! loop over nodes
  ! #if defined(SPLLT_USE_STARPU)
  !        ! TODO put in activate routine
  !        call starpu_f_void_data_register(fdata%nodes(snode)%hdl)
  ! #endif
  !        ! activate node: allocate factors, register handles
  !        call spllt_activate_node(snode, keep, fdata)
  !     end do

  ! ! #if defined(SPLLT_USE_STARPU)
  ! !     call starpu_f_task_wait_for_all()
  ! ! #endif    
  !     ! write(*,*)"num_nodes: ", num_nodes
  !     do snode = 1, num_nodes ! loop over nodes
  !        ! init node
  !        call spllt_init_node_task(fdata, fdata%nodes(snode), val, keep, huge(1))
  !     end do
  ! ! #if defined(SPLLT_USE_STARPU)
  ! !     call starpu_f_task_wait_for_all()
  ! ! #endif

  ! #if defined(SPLLT_USE_OMP)
  ! !$omp taskwait
  ! #endif

  !     do snode = 1, num_nodes

  !        ! write(*,*)'snode: ', snode

  !        node => keep%nodes(snode)

  !        ! task priority
  !        prio = (num_nodes - snode + 1)*4

  !        ! initialize node data
  !        sa = node%sa
  !        en = node%en
  !        numcol = en - sa + 1
  !        numrow = size(node%index)

  !        ! s_nb is the size of the blocks
  !        s_nb = node%nb

  !        nc = (numcol-1) / s_nb + 1 
  !        nr = (numrow-1) / s_nb + 1 

  !        !  Loop over descendent of snode
  !        least_desc = node%least_desc  

  !        ! build a map between row and block index
  !        call spllt_build_rowmap(node, tmpmap) 

  !        do d_num=least_desc, snode-1

  !           ! write(*,*)'d_num: ', d_num

  !           dnode => keep%nodes(d_num)

  !           ! initialize node data
  !           d_sa = dnode%sa
  !           d_en = dnode%en
  !           d_numcol = d_en - d_sa + 1
  !           d_numrow = size(dnode%index)

  !           ! s_nb is the size of the blocks
  !           d_nb = dnode%nb
  !           d_nc = (d_numcol-1) / d_nb + 1 

  !           cptr = 1 + d_numcol

  !           do cptr = cptr, d_numrow
  !              if(dnode%index(cptr) .ge. node%sa) exit
  !           end do
  !           if(cptr .gt. d_numrow) cycle ! finished with dnode

  !           ! map_done = .false. ! We will only build a map when we need it

  !           ! Loop over affected block columns of anode
  !           bcols: do
  !              if(cptr .gt. d_numrow) exit
  !              if(dnode%index(cptr) .gt. node%en) exit

  !              ! compute local index of block column in anode and find the id of 
  !              ! its diagonal block
  !              cb = (dnode%index(cptr) - node%sa)/node%nb + 1
  !              dblk = node%blk_sa
  !              do jb = 2, cb
  !                 dblk = blocks(dblk)%last_blk + 1
  !              end do

  !              ! Find cptr2
  !              jlast = min(node%sa + cb*node%nb - 1, node%en)
  !              do cptr2 = cptr, d_numrow
  !                 if(dnode%index(cptr2) > jlast) exit
  !              end do
  !              cptr2 = cptr2 - 1
  !              ! write(*,*)'cptr: ', cptr, ', cptr2: ', cptr2
  !              ! if(.not.map_done) call spllt_build_rowmap(node, tmpmap) 
  !              ! Loop over the blocks of snode
  !              ii = tmpmap(dnode%index(cptr)) 
  !              ! ii = -1
  !              ilast = cptr ! Set start of current block

  !              do i = cptr, d_numrow
  !                 k = tmpmap(dnode%index(i))

  !                 if(k.ne.ii) then

  !                    blk = dblk + ii - cb
  !                    bc => fdata%bc(blk)

  !                    d_dblk = dnode%blk_sa
  !                    ! Loop over the block columns in node. 
  !                    do kk = 1, d_nc

  !                       d_bc_kk => fdata%bc(d_dblk)

  !                       call spllt_update_between_task( &
  !                            & fdata, &
  !                            & d_bc_kk, dnode, bc, fdata%nodes(snode), &
  !                            & cptr, cptr2, ilast, i-1, &
  !                            & fdata%row_list, fdata%col_list, fdata%workspace, &
  !                            & keep%lfact, keep%blocks, fdata%bc, &
  !                            & control, prio)

  !                          d_dblk = blocks(d_dblk)%last_blk + 1
  !                    end do

  !                    ii = k
  !                    ilast = i ! Update start of current block      
  !                 end if
  !              end do

  !              blk = dblk + ii - cb
  !              bc => fdata%bc(blk)

  !              d_dblk = dnode%blk_sa
  !              ! Loop over the block columns in node. 
  !              do kk = 1, d_nc

  !                 d_bc_kk => fdata%bc(d_dblk)

  !                 call spllt_update_between_task( &
  !                      & fdata, &
  !                      & d_bc_kk, dnode, bc, fdata%nodes(snode), &
  !                      & cptr, cptr2, ilast, i-1, &
  !                      & fdata%row_list, fdata%col_list, fdata%workspace, &
  !                      & keep%lfact, keep%blocks, fdata%bc, &
  !                      & control, prio)

  !                 d_dblk = blocks(d_dblk)%last_blk + 1
  !              end do

  !              ! Move cptr down, ready for next block column of anode
  !              cptr = cptr2 + 1
  !           end do bcols

  !        end do

  !        ! first block in node
  !        dblk = node%blk_sa

  !        ! Loop over the block columns in node. 
  !        do kk = 1, nc

  !           ! #if defined(SPLLT_USE_OMP)
  !           ! !$omp taskwait
  !           ! #endif

  !           ! A_kk          

  !           bc_kk => fdata%bc(dblk)
  !           call spllt_factorize_block_task(fdata, fdata%nodes(snode), bc_kk, keep%lfact, prio+3)

  !           ! #if defined(SPLLT_USE_OMP)
  !           ! !$omp taskwait
  !           ! #endif

  !           ! #if defined(SPLLT_USE_STARPU)
  !           ! call starpu_f_task_wait_for_all()
  !           ! #endif
  !           ! loop over the row blocks (that is, loop over blocks in block col)
  !           do ii = kk+1,nr
  !              ! do ii = nr,kk+1,-1             
  !              ! A_mk
  !              blk = dblk+ii-kk
  !              ! bc_ik => keep%blocks(blk)
  !              bc_ik => fdata%bc(blk)

  !              call spllt_solve_block_task(fdata, bc_kk, bc_ik, keep%lfact,prio+2)
  !           end do

  !           ! #if defined(SPLLT_USE_OMP)
  !           ! !$omp taskwait
  !           ! #endif

  !           ! #if defined(SPLLT_USE_STARPU)
  !           !          call starpu_f_task_wait_for_all()
  !           ! #endif

  !           do jj = kk+1,nc

  !              ! L_jk
  !              blk2 = dblk+jj-kk
  !              bc_jk => fdata%bc(blk2)

  !              do ii = jj,nr

  !                 ! L_ik
  !                 blk1 = dblk+ii-kk                
  !                 bc_ik => fdata%bc(blk1)

  !                 ! A_ij
  !                 ! blk = get_dest_block(keep%blocks(blk1), keep%blocks(blk2))
  !                 blk = get_dest_block(keep%blocks(blk2), keep%blocks(blk1))
  !                 bc_ij => fdata%bc(blk)
  !                 call spllt_update_block_task(fdata, bc_ik, bc_jk, bc_ij, keep%lfact, prio+1)

  !              end do
  !           end do

  !           ! move to next block column in snode
  !           dblk = blocks(dblk)%last_blk + 1
  !           ! numrow = numrow - s_nb
  !        end do



  !     end do

  ! #if defined(SPLLT_USE_STARPU)
  !     ! unregister workspace handle
  !     call starpu_f_data_unregister_submit(fdata%workspace%hdl)

  !     ! unregister data handles
  !     call spllt_deinit_task(keep, fdata)

  !     ! do snode = 1, num_nodes
  !     !    call starpu_f_void_unregister_submit(fdata%nodes(snode)%hdl)        
  !     ! end do
  ! #endif

  !     call system_clock(stf_stop_t)
  !     write(*,'("[>] [spllt_stf_factorize] task insert time: ", es10.3, " s")') (stf_stop_t - stf_start_t)/real(stf_rate_t)

  !     return

  ! 9999 if(st.ne.0) then
  !      info%flag = spllt_error_allocation
  !      info%stat = st
  !      call spllt_print_err(info%flag, context='spllt_stf_factorize',st=st)
  !      return
  !   endif

  !   end subroutine spllt_stf_ll_factorize

end module spllt_stf_mod
