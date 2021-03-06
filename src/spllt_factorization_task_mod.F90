!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
module spllt_factorization_task_mod
  use spllt_data_mod
  implicit none
  
contains

  !*************************************************
  !

  subroutine spllt_scatter_block_task(rptr, rptr2, cptr, cptr2,  buffer, root, &
       & a_rptr, a_cptr, dest, anode)
    use spllt_data_mod
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use iso_c_binding
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#endif
    implicit none

    integer, intent(in) :: rptr, rptr2, cptr, cptr2
    type(spllt_block), intent(in) :: buffer
    type(spllt_node), target :: root
    integer, intent(in) :: a_rptr, a_cptr
    type(spllt_block) :: dest
    type(spllt_node), target :: anode
    
    integer :: rm, rn
    integer :: b_sz
    integer :: m, n
    integer :: bsa, ben

#if defined(SPLLT_USE_OMP)
    type(spllt_node), pointer :: p_root
    type(spllt_node), pointer :: p_anode
    real(wp), pointer :: buffer_c(:), dest_c(:)
    integer :: blkn, blkm
#endif

#if defined(SPLLT_USE_STARPU)

    call spllt_starpu_insert_subtree_scatter_block_task_c(rptr, rptr2, cptr, cptr2, &
         buffer%hdl, c_loc(root), a_rptr, a_cptr, dest%hdl, c_loc(anode))

#elif defined(SPLLT_USE_OMP)

    ! print *, "spllt_scatter_block_task"

    rm = size(root%index)
    rn = root%en - root%sa + 1
    b_sz = rm-rn

    m = rptr2-rptr+1
    n = cptr2-cptr+1

    bsa = (rptr-rn-1)*b_sz+cptr-rn
    ben = (rptr2-rn-1)*b_sz+cptr2-rn

    p_root => root 

    buffer_c => buffer%c

    p_anode => anode 

    dest_c => dest%c
    blkm = dest%blkm
    blkn = dest%blkn

!$omp task firstprivate(m, n, rptr, rptr2, cptr, cptr2, b_sz, blkm, blkn) &
!$omp      & firstprivate(bsa, ben) &
!$omp      & firstprivate(a_rptr, a_cptr, p_root, p_anode) &
!$omp      & firstprivate(buffer_c, dest_c) &
!$omp      & depend(in:buffer_c(1)) &
!$omp      & depend(inout:dest_c(1))

    call spllt_scatter_block(m, n,& 
         p_root%index(rptr:rptr2),&
         p_root%index(cptr:cptr2),&
         buffer_c(bsa:ben), b_sz, &
         p_anode%index(a_rptr),   &
         p_anode%index(a_cptr),   &
         blkm,                    &
         blkn,                    &
         dest_c, blkn)

!$omp end task

#else
    rm = size(root%index)
    rn = root%en - root%sa + 1
    b_sz = rm-rn

    m = rptr2-rptr+1
    n = cptr2-cptr+1

    bsa = (rptr-rn-1)*b_sz+cptr-rn
    ben = (rptr2-rn-1)*b_sz+cptr2-rn
    
    call spllt_scatter_block(m, n, &
         root%index(rptr:rptr2), &
         root%index(cptr:cptr2), &
         buffer%c(bsa:ben), b_sz, &
         anode%index(a_rptr), &
         anode%index(a_cptr), &
         dest%c, &
         dest%blkn)
#endif

  end subroutine spllt_scatter_block_task

  !*************************************************
  !

  subroutine spllt_subtree_factorize_task(root, fkeep, val, buffer, cntl, map)
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#endif
    use iso_c_binding
    implicit none

    integer, intent(in) :: root
    type(spllt_fkeep), target, intent(inout)  :: fkeep
    real(wp), dimension(:), target, intent(in) :: val ! user's matrix values
    ! type(spllt_keep), target, intent(inout) :: keep
    type(spllt_block), target, intent(inout) :: buffer ! update_buffer workspace
    type(spllt_options), target, intent(in) :: cntl
    integer, pointer, intent(inout) :: map(:)

#if defined(SPLLT_USE_STARPU)
    type(c_ptr) :: val_c, fkeep_c, cntl_c 
#endif

#if defined(SPLLT_USE_OMP)
    real(wp), pointer :: p_val(:) => null()
    real(wp), dimension(:), pointer :: buffer_c => null()
    type(spllt_block), pointer :: p_workspace(:) => null()
    type(spllt_workspace_i), pointer :: p_rlst(:) => null(), p_clst(:) => null()
    real(wp), dimension(:), pointer :: work => null()
    integer, dimension(:), pointer :: rlst => null(), clst => null()
    integer :: th_id ! thread id
    type(spllt_fkeep), pointer :: p_fkeep => null()
    type(spllt_options), pointer :: p_cntl => null()
    type(spllt_workspace_i), dimension(:), pointer :: p_map => null()
#endif

#if defined(SPLLT_USE_STARPU)
    
    val_c  = c_loc(val(1)) 
    fkeep_c = c_loc(fkeep)
    cntl_c = c_loc(cntl)

    call spllt_insert_subtree_factorize_task_c(root, val_c, fkeep_c, buffer%hdl, &
         & cntl_c, fkeep%map%hdl, fkeep%row_list%hdl, fkeep%col_list%hdl, &
         & fkeep%workspace%hdl)

#elif defined(SPLLT_USE_OMP)

    p_val => val
    buffer_c => buffer%c

    p_workspace => fkeep%workspace

    p_rlst => fkeep%row_list
    p_clst => fkeep%col_list

    p_map => fkeep%map

    p_fkeep => fkeep
    p_cntl => cntl

    th_id = 0

    if(size(buffer_c).ne.0) then
!$omp task firstprivate(p_val, buffer_c, p_workspace, p_rlst, p_clst) & 
!$omp & firstprivate(root, p_fkeep, p_cntl, p_map) &
!$omp & private(work, rlst, clst) &
!$omp & private(th_id) &
!$omp & depend(out:buffer_c(1))

!$ th_id = omp_get_thread_num()

    work => p_workspace(th_id)%c

    rlst => p_rlst(th_id)%c
    clst => p_clst(th_id)%c

    call spllt_subtree_factorize(root, p_val, p_fkeep, buffer_c, p_cntl, p_map(th_id)%c, &
         & rlst, clst, work)
    
!$omp end task
    else
!$omp task firstprivate(p_val, buffer_c, p_workspace, p_rlst, p_clst) & 
!$omp & firstprivate(root, p_fkeep, p_cntl, p_map) &
!$omp & private(work, rlst, clst) &
!$omp & private(th_id)

!$ th_id = omp_get_thread_num()

    work => p_workspace(th_id)%c

    rlst => p_rlst(th_id)%c
    clst => p_clst(th_id)%c

    call spllt_subtree_factorize(root, p_val, p_fkeep, buffer_c, p_cntl, p_map(th_id)%c, &
         & rlst, clst, work)
    
!$omp end task
    end if
#else
    
    call spllt_subtree_factorize(root, val, fkeep, buffer%c, cntl, map, &
         & fkeep%row_list%c, fkeep%col_list%c , fkeep%workspace%c)

#endif

  end subroutine spllt_subtree_factorize_task

  !*************************************************
  !

  subroutine spllt_factor_subtree_task(root, fkeep, buffer)
    use spllt_data_mod
    use spllt_kernels_mod
    implicit none

    integer, intent(in) :: root
    type(spllt_fkeep), target, intent(inout) :: fkeep
    real(wp), dimension(:), allocatable :: buffer ! update_buffer workspace

    type(spllt_node), pointer :: node ! node in the atree    
    integer                   :: m, n
    integer(long)             :: blk

    node => fkeep%nodes(root)
    blk = node%blk_sa

    m = fkeep%bc(blk)%blkm
    n = fkeep%bc(blk)%blkn

    if (size(buffer).lt.(m-n)**2) then
       deallocate(buffer)
       allocate(buffer((m-n)**2))
    end if

    ! perform factorization of subtree rooted at snode
    call spllt_factor_subtree(root, fkeep%nodes, fkeep%bc, fkeep%lfact, buffer)
 
  end subroutine spllt_factor_subtree_task

#if defined(SPLLT_USE_STARPU)

  !*************************************************
  !
  ! Deinitialize factorization
  ! StarPU: unregister data handles (block handles) in StarPU
  subroutine spllt_data_unregister_task(akeep, fkeep)
    use spllt_data_mod
    use  starpu_f_mod
    implicit none

    type(spllt_akeep), intent(in) :: akeep
    type(spllt_fkeep), target, intent(inout) :: fkeep

    integer :: i
    integer :: snode, num_nodes
    type(spllt_node), pointer :: node ! node in the atree    
    integer(long) :: blk, dblk
    integer :: l_nb, sz, sa, en

    num_nodes = fkeep%info%num_nodes
    blk = 1

    do snode = 1, num_nodes

       ! if (akeep%small(snode) .ne. 0) cycle
       ! print *, "[data_unregister_task] node", snode, ", small: ", akeep%small(snode)

       node => fkeep%nodes(snode)

       l_nb = node%nb
       sz = (size(node%index) - 1) / l_nb + 1
       sa = node%sa
       en = node%en
       
       do i = sa, en, l_nb
          dblk = blk
          do blk = dblk, dblk+sz-1
             ! write(*,*) 'unregister blk: ', blk

             call starpu_f_data_unregister_submit(fkeep%bc(blk)%hdl)
          end do
          sz = sz - 1
       end do
    end do

  end subroutine spllt_data_unregister_task

  !*************************************************
  !

  subroutine spllt_data_unregister(fkeep)
    use spllt_data_mod
    use starpu_f_mod
    implicit none

    type(spllt_fkeep), target, intent(inout) :: fkeep

    integer :: i
    integer :: snode, num_nodes
    type(spllt_node), pointer :: node ! node in the atree    
    integer(long) :: blk, dblk
    integer :: l_nb, sz, sa, en

    num_nodes = fkeep%info%num_nodes
    blk = 1

    do snode = 1, num_nodes

       node => fkeep%nodes(snode)

       l_nb = node%nb
       sz = (size(node%index) - 1) / l_nb + 1
       sa = node%sa
       en = node%en

       do i = sa, en, l_nb
          dblk = blk
          do blk = dblk, dblk+sz-1
             ! write(*,*) 'unregister blk: ', blk
             call starpu_f_data_unregister(fkeep%bc(blk)%hdl)
          end do
          sz = sz - 1
       end do
    end do

  end subroutine spllt_data_unregister
#endif

  !*************************************************
  !
  ! Factorize block task 
  ! _potrf
  subroutine spllt_factorize_block_task(fkeep, node, bc, lfact, prio)
    use spllt_data_mod
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none
    
    type(spllt_fkeep), target, intent(in)  :: fkeep
    type(spllt_node), intent(in)          :: node
#if defined(SPLLT_USE_OMP)
    type(spllt_block), pointer, intent(inout) :: bc ! block to be factorized    
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
#else
    type(spllt_block), target, intent(inout) :: bc ! block to be factorized    
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    integer, optional :: prio

    integer :: m, n, bcol, sa
    integer(long) :: id
    ! type(block_type), pointer :: blk ! block to be factorized
    integer :: p
#if defined(SPLLT_USE_STARPU)
    type(c_ptr) :: node_hdl
#elif defined(SPLLT_USE_OMP)
    real(wp), dimension(:), pointer :: lcol
    real(wp), dimension(:), pointer :: bc_c
#if defined(SPLLT_OMP_TRACE)
    integer :: th_id
#endif
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)

    ! pass snode handle to tasks only in first block in snode
    node_hdl = c_null_ptr
    if (node%blk_sa .eq. bc%id) then
    node_hdl = node%hdl
    end if

    call spllt_starpu_insert_factorize_block_c(bc%hdl, node_hdl, p)

    ! call spllt_starpu_insert_factorize_block(bc, p)
#elif defined(SPLLT_USE_OMP)

    ! blk => bc%blk
    bcol = bc%bcol
    lcol => lfact(bcol)%lcol

    m    = bc%blkm
    n    = bc%blkn
    sa   = bc%sa
    id   = bc%id
    bc_c => fkeep%bc(id)%c

! !$omp taskwait

!$omp task firstprivate(m, n, sa) &
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(fac_blk_id) &
#endif
!$omp    & firstprivate(bc_c, bcol, lcol) &
! !$omp    & depend(inout:lcol(sa:sa+n*m-1))
! !$omp    & depend(inout:fkeep%bc(id)%c(:))
!$omp    & depend(inout:bc_c(1))
! !$omp    & priority(p)

! !$omp    & depend(inout:bc%c(1))
! !$omp    & depend(inout:bc_c)
! !$omp    & depend(inout:bc_c(1))

! !$omp critical

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(fac_blk_id, th_id)
#endif

    ! write(*,*)"bcol: ", bcol
    call spllt_factor_diag_block(m, n, &
         ! & lfact(bcol)%lcol(sa:sa+n*m-1) &
         & lcol(sa:sa+n*m-1) &
         &)

#if defined(SPLLT_OMP_TRACE)
    ! write(*,*)"fac_blk_id: ", fac_blk_id
    call trace_event_stop(fac_blk_id, th_id)
#endif

! !$omp end critical

!$omp end task

! !$omp taskwait

#else    

    ! blk => bc%blk

    m    = bc%blkm
    n    = bc%blkn
    id   = bc%id 
    bcol = bc%bcol

    sa   = bc%sa

    ! factorize_block
    ! call factor_diag_block(n, m, id, &
    !      & lfact(bcol)%lcol(sa:sa+n*m-1),   &
    !      & control, flag, detlog(0))

    call spllt_factor_diag_block(m, n, lfact(bcol)%lcol(sa:sa+n*m-1))
#endif
    
    return
  end subroutine spllt_factorize_block_task

  ! _trsm
  subroutine spllt_solve_block_task(fkeep, bc_kk, bc_ik, lfact, prio)
    use spllt_data_mod
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none

    type(spllt_fkeep), target, intent(inout)  :: fkeep
#if defined(SPLLT_USE_OMP)
    type(spllt_block), pointer, intent(inout) :: bc_kk, bc_ik ! block to be factorized    
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
    ! type(block_type), pointer :: blk_kk, blk_ik ! block to be factorized
#else
    type(spllt_block), intent(inout) :: bc_kk, bc_ik ! block to be factorized    
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    integer, optional :: prio 
    
    integer :: m, n, bcol, sa
    integer :: d_m, d_n, d_sa
    integer(long) :: id, d_id
    integer :: p
#if defined(SPLLT_USE_OMP)
    real(wp), dimension(:), pointer :: lcol
    real(wp), dimension(:), pointer :: bc_kk_c, bc_ik_c
    integer :: th_id
    integer :: snode
    integer(long) :: blk_sa
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)    
    call spllt_starpu_insert_solve_block_c(bc_kk%hdl, bc_ik%hdl, p)
    ! call spllt_starpu_insert_solve_block(bc_kk, bc_ik, 0)

#elif defined(SPLLT_USE_OMP)
    

    ! blk_kk => bc_kk%blk
    ! blk_ik => bc_ik%blk
    ! bcol is block column that blk and dblk belong to
    bcol = bc_kk%bcol    
    ! write(*,*)"bcol: ", bcol
    lcol => lfact(bcol)%lcol

    ! bc_kk
    d_m  = bc_kk%blkm
    d_n  = bc_kk%blkn
    d_sa = bc_kk%sa
    d_id = bc_kk%id

    ! bc_ik
    n  = bc_ik%blkn
    m  = bc_ik%blkm
    sa = bc_ik%sa
    id = bc_ik%id
    ! write(*,*)"solve id: ", id

    bc_kk_c => fkeep%bc(d_id)%c
    bc_ik_c => fkeep%bc(id)%c

    snode = bc_kk%node
    blk_sa = fkeep%nodes(snode)%blk_sa
    ! write(*,*)"blk_sa: ", blk_sa
! !$omp taskwait

!$omp task private(th_id) &
!$omp    & firstprivate(m, n, sa, d_m, d_n, d_sa) &
!$omp    & firstprivate(lcol, bcol) &
! !$omp    & firstprivate(blk_kk, blk_ik) &
!$omp    & firstprivate(bc_kk_c, bc_ik_c) &
!$omp    & firstprivate(d_id, id) &
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(slv_blk_id) &
#endif
!$omp    & depend(in:bc_kk_c(1)) &
!$omp    & depend(inout:bc_ik_c(1))

! !$omp    & depend(in:fkeep%bc(d_id)%c(1)) &
! !$omp    & depend(inout:fkeep%bc(id)%c)

! !$omp    & depend(in:bc_kk%c(1)) &
! !$omp    & depend(in:bc_kk_c) &
! !$omp    & depend(in:bc_kk_c(1:1)) &
! !$omp    & depend(inout:bc_ik_c(1:1))

! !$omp    & depend(inout:bc_ik%c(1))
! !$omp    & depend(in:lcol(d_sa:d_n*d_m)) &
! !$omp    & depend(inout:lcol(sa:n*m-1))

! !$omp critical

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(slv_blk_id, th_id)
#endif

    ! call dtrsm('Left', 'Upper', 'Transpose', 'Non-Unit', n, m, &
    !      & one, lcol(d_sa:d_sa+d_n*d_m), n, &
    !      & lcol(sa:sa+n*m-1), n)

    ! call dtrsm('Left', 'Upper', 'Transpose', 'Non-Unit', n, m, &
    !      & one, bc_kk_c(1:d_n*d_m), n, &
    !      & bc_ik_c(1:n*m), n)

    call spllt_solve_block(m, n, &
         & lcol(sa:sa+n*m-1), & 
         & lcol(d_sa:d_sa+d_n*d_m))

#if defined(SPLLT_OMP_TRACE)
    ! write(*,*)"fac_blk_id: ", fac_blk_id
    call trace_event_stop(slv_blk_id, th_id)
#endif

! !$omp end critical

!$omp end task

! !$omp taskwait

#else    

    ! bc_kk
    d_m  = bc_kk%blkm
    d_n  = bc_kk%blkn
    d_sa = bc_kk%sa
    d_id = bc_kk%id

    ! bc_ik
    n  = bc_ik%blkn
    m  = bc_ik%blkm
    sa = bc_ik%sa
    id = bc_ik%id
    
    ! bcol is block column that blk and dblk belong to
    bcol = bc_kk%bcol    

    ! solve_block task
    ! call solv_col_block(m, n, id, & 
    !      & lfact(bcol)%lcol(sa:sa+n*m-1), &
    !      & d_id, lfact(bcol)%lcol(d_sa:d_sa+d_n*d_m), &
    !      & control)

    call spllt_solve_block(m, n, &
         & lfact(bcol)%lcol(sa:sa+n*m-1), & 
         & lfact(bcol)%lcol(d_sa:d_sa+d_n*d_m))

#endif

    return
  end subroutine spllt_solve_block_task

  ! syrk/gemm (same node)
  ! A_ij <- A_ij - A_ik A_jk^T
  subroutine spllt_update_block_task(fkeep, bc_ik, bc_jk, bc_ij, lfact, prio)
    use spllt_data_mod
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none
    
    type(spllt_fkeep), target, intent(in)  :: fkeep
    ! type(block_type), intent(inout) :: bc_ik, bc_jk, bc_ij ! block to be updated    
#if defined(SPLLT_USE_OMP)
    type(spllt_block), pointer, intent(inout) :: bc_ik, bc_jk, bc_ij
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
#else
    type(spllt_block), intent(inout) :: bc_ik, bc_jk, bc_ij
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    integer, optional :: prio 

    integer :: n1, m1, sa1, n2, m2, sa2, n, m, sa
    integer :: bcol1, bcol, bcol2
    integer :: p
#if defined(SPLLT_USE_STARPU)
    integer :: d
#endif
#if defined(SPLLT_USE_OMP)
    logical :: is_diag
    real(wp), dimension(:), pointer :: lcol1, lcol2, lcol
    real(wp), dimension(:), pointer :: bc_ik_c, bc_jk_c, bc_ij_c
    ! type(block_type), pointer :: blk_ik, blk_jk, blk_ij ! block to be factorized
    integer :: th_id
    integer(long) :: id_ik, id_jk, id_ij
    integer :: snode
    integer(long) :: blk_sa
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)
    if (bc_ij%dblk .eq. bc_ij%id) then
       d = 1
    else
       d = 0
    end if

    call spllt_starpu_insert_update_block_c(bc_ik%hdl, bc_jk%hdl, bc_ij%hdl, d, p)
#elif defined(SPLLT_USE_OMP)

    ! blk_ik => bc_ik%blk
    ! blk_jk => bc_jk%blk
    ! blk_ij => bc_ij%blk

    id_ik = bc_ik%id
    id_jk = bc_jk%id
    id_ij = bc_ij%id

    bc_ik_c => fkeep%bc(id_ik)%c
    bc_jk_c => fkeep%bc(id_jk)%c
    bc_ij_c => fkeep%bc(id_ij)%c

    bcol1 = bc_ik%bcol
    lcol1 => lfact(bcol1)%lcol

    bcol2 = bc_jk%bcol
    lcol2 => lfact(bcol2)%lcol

    bcol = bc_ij%bcol
    lcol => lfact(bcol)%lcol

    snode = bc_ij%node
    blk_sa = fkeep%nodes(snode)%blk_sa
    
! !$omp taskwait

!$omp task private(n1, m1, sa1, n2, m2, sa2, &
!$omp    & n, m, sa, is_diag, th_id) &
!$omp    & firstprivate(bcol2, lcol2, bcol1, lcol1, bcol, lcol) &
!$omp    & firstprivate(bc_ik, bc_jk, bc_ij) &
!$omp    & firstprivate(bc_ik_c, bc_jk_c, bc_ij_c) &
! !$omp    & firstprivate(blk_ik, blk_jk, blk_ij) &
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(upd_blk_id) &
#endif
!$omp    & depend(in:bc_ik_c(1), bc_jk_c(1)) &
!$omp    & depend(inout: bc_ij_c(1))

! !$omp    & depend(in:fkeep%bc(id_ik)%c, fkeep%bc(id_jk)%c) &
! !$omp    & depend(inout: fkeep%bc(id_ij)%c)

! !$omp    & depend(in:bc_ik%c(1), bc_jk%c(1)) &
! !$omp    & depend(in:bc_ik_c(1:1), bc_jk_c(1:1)) &
! !$omp    & depend(in:bc_ik%c(1), bc_jk%c(1)) &
! !$omp    & depend(in:bc_ik_c(1:1), bc_jk_c(1:1)) &
! !$omp    & depend(inout: bc_ij_c(1:1))

! !$omp    & depend(in:lcol1(sa1:n1*m1-1)) &
! !$omp    & depend(in:lcol2(sa2:n1*m2-1)) &
! !$omp    & depend(inout: lcol(sa:n*m-1))

! !$omp critical

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(upd_blk_id, th_id)
#endif

    ! bc_ik
    n1  = bc_ik%blkn
    m1  = bc_ik%blkm
    sa1 = bc_ik%sa

    ! bc_jk
    n2  = bc_jk%blkn
    m2  = bc_jk%blkm
    sa2 = bc_jk%sa

    ! bc_ij
    n  = bc_ij%blkn
    m  = bc_ij%blkm
    sa = bc_ij%sa

    is_diag = bc_ij%dblk.eq.bc_ij%id 

    ! write(*,*)"bc_ik id: ", bc_ik%id, ", m: ", m1, ", n: ", n1
    ! write(*,*)"bc_jk id: ", bc_jk%id, ", m: ", m2, ", n: ", n2
    ! write(*,*)"bc_ij id: ", bc_ij%id, ", m: ", m, ", n: ", n

    ! write(*,*) "size lfact(bcol1)%lcol: ", size(lfact(bcol1)%lcol)
    ! write(*,*) "sa2+n1*m2-1: ", sa2+n1*m2-1

    call spllt_update_block(m, n, &
         & lcol(sa:sa+n*m-1), &
         & is_diag, n1, &
         & lcol1(sa2:sa2+n1*m2-1), &
         & lcol1(sa1:sa1+n1*m1-1))

#if defined(SPLLT_OMP_TRACE)
    ! write(*,*)"fac_blk_id: ", fac_blk_id
    call trace_event_stop(upd_blk_id, th_id)
#endif

! !$omp end critical

!$omp end task

! !$omp taskwait

#else

    ! bc_ik
    n1    = bc_ik%blkn
    m1    = bc_ik%blkm
    sa1   = bc_ik%sa
    bcol1 = bc_ik%bcol
    
    ! bc_jk
    n2    = bc_jk%blkn
    m2    = bc_jk%blkm
    sa2   = bc_jk%sa
    bcol2 = bc_jk%bcol
    
    ! bc_ij
    n  = bc_ij%blkn
    m  = bc_ij%blkm
    sa = bc_ij%sa
    
    bcol = bc_ij%bcol

    ! write(*,*)"bc_ik id: ", bc_ik%id, ", m: ", m1, ", n: ", n1
    ! write(*,*)"bc_jk id: ", bc_jk%id, ", m: ", m2, ", n: ", n2
    ! write(*,*)"bc_ij id: ", bc_ij%id, ", m: ", m, ", n: ", n
    
    ! write(*,*) "size lfact(bcol1)%lcol: ", size(lfact(bcol1)%lcol)
    ! write(*,*) "sa2+n1*m2-1: ", sa2+n1*m2-1

    call spllt_update_block(m, n, &
         & lfact(bcol)%lcol(sa:sa+n*m-1), &
         & bc_ij%dblk.eq.bc_ij%id, n1, &
         & lfact(bcol1)%lcol(sa2:sa2+n1*m2-1), &
         & lfact(bcol1)%lcol(sa1:sa1+n1*m1-1))

#endif

    return
  end subroutine spllt_update_block_task

! #if defined(SPLLT_USE_OMP)
!   subroutine spllt_omp_update_between_cpu_func(m, n, blk, dcol, dnode, n1, scol, snode, dest, csrc, rsrc, min_width_blas)
!     use spllt_data_mod
!     use spllt_kernels_mod
!     implicit none

!     integer, intent(in) :: m  ! number of rows in destination block
!     integer, intent(in) :: n  ! number of columns in destination block
!     ! integer(long), intent(in) :: blk ! identifier of destination block
!     type(block_type), intent(in) :: blk ! destination block
!     integer, intent(in) :: dcol ! index of block column that blk belongs to in dnode
!     type(spllt_node), intent(in) :: dnode ! Node to which blk belongs
!     integer :: n1 ! number of columns in source block column
!     ! integer(long), intent(in) :: src  ! identifier of block in source block col
!     integer, intent(in) :: scol ! index of block column that src belongs to in snode
!     type(spllt_node), intent(in) :: snode ! Node to which src belongs
!     real(wp), dimension(*), intent(inout) :: dest ! holds block in L
!     ! that is to be updated.
!     real(wp), dimension(*), intent(in) :: csrc ! holds csrc block
!     real(wp), dimension(*), intent(in) :: rsrc ! holds rsrc block
!     ! type(block_type), dimension(:), intent(inout) :: blocks
!     ! real(wp), dimension(:), allocatable :: buffer
!     integer, intent(in) :: min_width_blas      ! Minimum width of source block

!     real(wp), dimension(:), pointer :: buffer
!     integer, dimension(:), allocatable :: row_list ! reallocated to min size m
!     integer, dimension(:), allocatable :: col_list ! reallocated to min size n

!     allocate(buffer(m*n))
!     allocate(row_list(1), col_list(1))

!     call spllt_update_between(m, n, blk, dcol, dnode, &
!          & n1, scol, snode, &
!          & dest, &
!          & csrc, &
!          & rsrc, &
!          & row_list, col_list, buffer, &
!          & min_width_blas)

!     deallocate(buffer)
!     deallocate(row_list, col_list)

!     return

!   end subroutine spllt_omp_update_between_cpu_func
! #endif

  ! syrk/gemm (inter-node)
  subroutine spllt_update_between_task(fkeep, &
       & dbc, snode, a_bc, anode, &
       ! & csrc, rsrc, &
       & cptr, cptr2, rptr, rptr2, &
       & row_list, col_list, workspace, &
       & lfact, bcs, &
       & min_width_blas, prio)
    use spllt_data_mod
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
!$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none

    type(spllt_fkeep), target, intent(in)       :: fkeep
#if defined(SPLLT_USE_OMP)
    type(spllt_block), pointer, intent(inout)     :: a_bc ! dest block
    type(spllt_block), pointer, intent(inout)     :: dbc ! diag block in source node

    type(spllt_block), allocatable, target        :: workspace(:)
    type(spllt_block), target                     :: bcs(:) ! block info.
    type(spllt_workspace_i), allocatable, target    :: row_list(:), col_list(:)
#else
    type(spllt_block), target, intent(inout)              :: a_bc ! dest block
    type(spllt_block), intent(in)                 :: dbc ! diag block in source node

    type(spllt_block)                             :: workspace
    type(spllt_block)                             :: bcs(:) ! block info.
    type(spllt_workspace_i)                         :: row_list, col_list
#endif

#if defined(SPLLT_USE_STARPU)
    type(spllt_node), target                         :: snode ! src node
    type(spllt_node), target                   :: anode ! dest node
#elif defined(SPLLT_USE_OMP)
    type(spllt_node), pointer, intent(in)            :: snode ! src node
    type(spllt_node), pointer, intent(in)      :: anode ! dest node
#else
    type(spllt_node)                                 :: snode ! src node
    type(spllt_node)                           :: anode ! dest node
#endif
    integer                                         :: cptr, cptr2, rptr, rptr2 
    !    integer :: csrc(2), rsrc(2) ! used for update_between tasks to

    ! real(wp), dimension(:), allocatable :: buffer ! update_buffer workspace
    ! type(block_type), dimension(:)                  :: blocks ! block info. 

#if defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)
    type(lfactor), allocatable, target, intent(inout)  :: lfact(:)
#else
    type(lfactor), allocatable, intent(inout)          :: lfact(:)
#endif
    integer :: min_width_blas
    ! integer, intent(out) :: st ! TODO error managment
    integer, optional :: prio 

    ! integer :: info ! TODO error managment
    
    integer :: p
    integer :: m, n, n1, sa
    integer :: bcol, bcol1
!   integer(long) :: id, id1
    integer :: dcol, scol
    integer :: csrc, csrc2, rsrc, rsrc2
    integer :: s_nb

#if defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)
    integer(long) :: blk_sa, blk_en, nb_blk, dblk 
#endif  

#if defined(SPLLT_USE_STARPU)
    integer :: blkn, ljk_sa
    integer :: nhljk, nhlik
    integer(c_int) :: ljk_m, ljk_n, lik_m
    integer(long) :: blk
    type(c_ptr), dimension(:), allocatable :: lik_handles, ljk_handles
    type(c_ptr) :: ljk_hdl = c_null_ptr, lik_hdl = c_null_ptr
    type(c_ptr) :: snode_c, anode_c, a_bc_c
#endif

#if defined(SPLLT_USE_OMP)
    real(wp), dimension(:), pointer :: lcol1, lcol
    integer :: th_id
    integer :: csrc_sa, rsrc_sa, csrc_en, rsrc_en
    ! type(spllt_block), pointer :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    real(wp), dimension(:), pointer :: a_bc_c, dbc_c
    ! real(wp), dimension(:), pointer :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    ! type(spllt_block), pointer    :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    real(wp), dimension(:), pointer    :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    ! type(block_type), pointer :: blk_kk, a_blk
    type(spllt_block), pointer :: p_workspace(:) => null()
    type(spllt_workspace_i), pointer :: p_rlst(:) => null(), p_clst(:) => null()
    ! real(wp), dimension(:), allocatable :: work
    real(wp), dimension(:), pointer :: work
    integer, dimension(:), pointer :: rlst, clst ! row and column lists
    type(spllt_node), pointer :: p_snode ! src node
    type(spllt_node), pointer :: p_anode ! dest node
    integer(long) :: id_ik_sa, id_ik_en, id_jk_sa, id_jk_en, id_ij
    real(wp), dimension(:), pointer :: lcol_ik, lcol_jk
    type(lfactor), pointer :: p_lfact(:)
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

    ! print *, 'update_between_task'

    s_nb = snode%nb    ! block size in source node
    n1  = dbc%blkn ! width of column

    bcol1 = dbc%bcol
    scol = bcol1 - bcs(snode%blk_sa)%bcol + 1

    bcol = a_bc%bcol
    dcol = bcol - bcs(anode%blk_sa)%bcol + 1
    ! dcol = a_bc%blk%bcol ! blocks(anode%node%blk_sa)%bcol + 1

#if defined(SPLLT_USE_STARPU)
! #if 0
    ! write(*,*)"associated blk: ", associated(a_bc%blk)
    dblk = dbc%id

    ! ljk
    blk_sa = (cptr -1)/s_nb - (scol-1) + dblk
    blk_en = (cptr2-1)/s_nb - (scol-1) + dblk
    nb_blk = blk_en-blk_sa+1

    allocate(ljk_handles(nb_blk))
    ljk_m = 0 ! compute height of ljk factor
    nhljk=0
    do blk=blk_sa,blk_en
       nhljk=nhljk+1
       ljk_handles(nhljk) = bcs(blk)%hdl
       ljk_m = ljk_m + bcs(blk)%blkm
    end do
    ljk_n = dbc%blkn
    ! write(*,*) "ljk_hdl: ", ljk_hdl, ", ptr: ", c_loc(fkeep%bc(blk_sa)%c(1))
#if defined(SPLLT_USE_GPU)
    ! create temporary handle for gathering blocks in Ljk factor
    ! DEBUG
    call starpu_matrix_data_register(ljk_hdl, fkeep%bc(blk_sa)%mem_node, &
         & c_loc(fkeep%bc(blk_sa)%c(1)), ljk_m, ljk_m, ljk_n, &
         & int(wp,kind=c_size_t))
    ! call starpu_f_data_unregister_submit(ljk_hdl)
#endif
    ! return

    ! lik
    blk_sa = (rptr -1)/s_nb - (scol-1) + dblk
    blk_en = (rptr2-1)/s_nb - (scol-1) + dblk
    nb_blk = blk_en-blk_sa+1

    allocate(lik_handles(nb_blk))
    lik_m = 0 ! compute height of lik factor
    nhlik=0
    do blk=blk_sa,blk_en
       nhlik=nhlik+1
       lik_handles(nhlik) = bcs(blk)%hdl
       lik_m = lik_m + bcs(blk)%blkm
    end do

    ! write(*,*) "lik_m: ", lik_m

#if defined(SPLLT_USE_GPU)
    ! create temporary handle for gathering blocks in Ljk factor
    call starpu_matrix_data_register(lik_hdl, fkeep%bc(blk_sa)%mem_node, &
         & c_loc(fkeep%bc(blk_sa)%c(1)), lik_m, lik_m, ljk_n, &
         & int(wp,kind=c_size_t))
    ! call starpu_f_data_unregister_submit(lik_hdl)
#endif
    ! write(*,*) "ljk_hdl: ", ljk_hdl, "lik_hdl: ", lik_hdl
    ! write(*,*) "ljk_hdl: ", ljk_hdl, "lik_hdl: ", lik_hdl, ", ptr: ", c_loc(fkeep%bc(blk_sa)%c(1))
    ! return
    
    csrc  = 1 + (mod(cptr-1, s_nb))*n1
    csrc2 = (cptr2 - cptr + 1)*n1

    rsrc  = 1 + (mod(rptr-1, s_nb))*n1
    rsrc2 = (rptr2 - rptr + 1)*n1

    snode_c = c_loc(snode)
    anode_c = c_loc(anode)
    a_bc_c  = c_loc(a_bc)
    
    ! write(*,*)"dcol: ", dcol
#if defined(SPLLT_USE_GPU)

    call spllt_insert_unpartition_task_c(lik_hdl, lik_handles, nhlik, p)
    call spllt_insert_unpartition_task_c(ljk_hdl, ljk_handles, nhljk, p)

    call spllt_starpu_insert_update_between_c(&
         & lik_hdl, &
         & ljk_hdl, &
         & a_bc%hdl, &
         & snode_c, scol, &
         & anode_c, a_bc_c, dcol, &
         & csrc, csrc2, rsrc, rsrc2, &
         & min_width_blas, &
         & workspace%hdl, row_list%hdl, col_list%hdl, &
         & anode%hdl, &
         & p &
         &)
#else
    call spllt_starpu_insert_update_between_c(&
         & lik_handles, nhlik, &
         & ljk_handles, nhljk, &
         & a_bc%hdl, &
         & snode_c, scol, &
         & anode_c, a_bc_c, dcol, &
         & csrc, csrc2, rsrc, rsrc2, &
         & min_width_blas, &
         & workspace%hdl, row_list%hdl, col_list%hdl, &
         & anode%hdl, &
         & p &
         &)
#endif


#if defined(SPLLT_USE_GPU)
    ! unregister temporary handles when it is no longer needed
    call starpu_f_data_unregister_submit(ljk_hdl)
    call starpu_f_data_unregister_submit(lik_hdl)
#endif

    ! call test_insert_c(lik_handles, nhlik, ljk_handles, nhljk)

    ! call starpu_f_task_wait_for_all()

    deallocate(ljk_handles)
    deallocate(lik_handles)
    
#elif defined(SPLLT_USE_OMP)

    p_lfact => lfact
    ! a_blk => a_bc%blk

    ! ! lik
    ! blk_sa = (rptr -1)/s_nb - (scol-1) + dblk
    ! blk_en = (rptr2-1)/s_nb - (scol-1) + dblk
    ! nb_blk = blk_en-blk_sa+1
    
    ! write(*,*) "update_between_task"
    ! call update_between(m, n, id, anode, &
    !      & n1, id1, snode, &
    !      & lfact(bcol)%lcol(sa:sa+m*n-1), &
    !      & lfact(bcol1)%lcol(csrc(1):csrc(1)+csrc(2)-1), &
    !      & lfact(bcol1)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1), &
    !      & blocks, row_list, col_list, buffer, &
    !      & control, info, st)
    ! blk_kk => dbc%blk
    dblk = dbc%id

    ! ljk
    blk_sa = (cptr -1)/s_nb - (scol-1) + dblk
    blk_en = (cptr2-1)/s_nb - (scol-1) + dblk
    nb_blk = blk_en-blk_sa+1
    ! bc_jk_sa => bcs(blk_sa)%c
    ! bc_jk_en => bcs(blk_en)%c
    id_jk_sa = blk_sa 
    id_jk_en = blk_en
    ! write(*,*) "nb_blk:", nb_blk, "id_jk_sa: ", id_jk_sa, ", id_jk_en: ", id_jk_en
    ! if (nb_blk.gt.2) write(*,*) "nb_blk:", nb_blk, "id_jk_sa: ", id_jk_sa, ", id_jk_en: ", id_jk_en
    bc_jk_sa => fkeep%bc(id_jk_sa)%c ! bcs(blk_sa)%c
    bc_jk_en => fkeep%bc(id_jk_en)%c ! bcs(blk_en)%c
    ! write(*,*) "nb_blk: ", nb_blk
    ! write(*,*) "nb_blk:", nb_blk, "blk_sa: ", blk_sa, ", blk_en: ", blk_en
    csrc_sa = bcs(blk_sa)%sa
    csrc_en = bcs(blk_en)%sa

    ! csrc_sa = 1 + (cptr_sa-(scol-1)*s_nb-1)*n1
    ! csrc_le = (cptr_en - cptr_sa + 1)*n1

    ! lik
    blk_sa = (rptr -1)/s_nb - (scol-1) + dblk
    blk_en = (rptr2-1)/s_nb - (scol-1) + dblk
    nb_blk = blk_en-blk_sa+1
    ! bc_ik_sa => bcs(blk_sa)%c
    ! bc_ik_en => bcs(blk_en)%c
    id_ik_sa = blk_sa 
    id_ik_en = blk_en
    ! if (nb_blk.gt.2)  write(*,*) "nb_blk:", nb_blk, "id_ik_sa: ", id_ik_sa, ", id_ik_en: ", id_ik_en

    bc_ik_sa => fkeep%bc(id_ik_sa)%c ! bcs(blk_sa)%c
    bc_ik_en => fkeep%bc(id_ik_en)%c ! bcs(blk_en)%c

    rsrc_sa = bcs(blk_sa)%sa
    rsrc_en = bcs(blk_en)%sa

    id_ij = a_bc%id
    a_bc_c => fkeep%bc(id_ij)%c

    dbc_c => dbc%c

    bcol1 = dbc%bcol
    bcol  = a_bc%bcol

    lcol1 => lfact(bcol1)%lcol
    lcol  => lfact(bcol)%lcol

    p_workspace => workspace

    p_rlst => row_list
    p_clst => col_list

    ! min_width_blas = control%min_width_blas

    ! p_snode => snode
    ! p_anode => anode

    p_snode => fkeep%nodes(dbc%node)
    p_anode => fkeep%nodes(a_bc%node)

    s_nb = p_snode%nb    ! block size in source node
    n1  = dbc%blkn ! width of column

    csrc  = 1 + (cptr-(scol-1)*s_nb-1)*n1
    csrc2 = (cptr2 - cptr + 1)*n1

    rsrc  = 1 + (rptr-(scol-1)*s_nb-1)*n1
    rsrc2 = (rptr2 - rptr + 1)*n1

    lcol_ik => lcol1(rsrc:rsrc+rsrc2-1)
    lcol_jk => lcol1(csrc:csrc+csrc2-1)

    m  = a_bc%blkm
    n  = a_bc%blkn
    sa = a_bc%sa

    blk_sa = snode%blk_sa

! !$omp taskwait    

!$omp task firstprivate(m, n, a_bc, dcol, p_anode, n1, scol, p_snode, &
!$omp    & lcol, lcol1, min_width_blas, sa, csrc, csrc2, rsrc, rsrc2) &
!$omp    & private(rlst, clst, work) & 
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(upd_btw_id) &
#endif
!$omp    & firstprivate(bc_ik_sa, bc_ik_en, bc_jk_sa, bc_jk_en, a_bc_c) &
!$omp    & depend(in:bc_ik_sa(1), bc_ik_en(1)) &
!$omp    & depend(in:bc_jk_sa(1), bc_jk_en(1)) &
!$omp    & depend(inout: a_bc_c(1))

! !$omp    & depend(in:fkeep%bc(id_ik_sa)%c, fkeep%bc(id_ik_en)%c) &
! !$omp    & depend(in:fkeep%bc(id_jk_sa)%c, fkeep%bc(id_jk_en)%c) &
! !$omp    & depend(inout:fkeep%bc(id_ij)%c)

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"upd_btw_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(upd_btw_id, th_id)
#endif
     ! write(*,*)"min_width_blas: ", min_width_blas
    ! write(*,*)"thread id: ", th_id
    ! write(*,*)"m: ", m, ", n: ", n, ", n1: ", n1

    th_id = 0
    !$  th_id = omp_get_thread_num()

    ! call spllt_omp_update_between_cpu_func(m, n, a_bc, dcol, p_anode%node, &
    !      & n1, scol, p_snode, &
    !      & lcol(sa:sa+m*n-1), &
    !      & lcol1(csrc:csrc+csrc2-1), &
    !      & lcol1(rsrc:rsrc+rsrc2-1), &
    !      & min_width_blas)


    ! write(*,*)"th_id: ", th_id
    work => p_workspace(th_id)%c

    rlst => p_rlst(th_id)%c
    clst => p_clst(th_id)%c
    ! allocate(work(m*n))
    ! work(:) = 0
    ! allocate(rlst(1), clst(1))

    call spllt_update_between(m, n, a_bc, dcol, p_anode, &
         & n1, scol, p_snode, &
         & lcol(sa:sa+m*n-1), &
         & lcol1(csrc:csrc+csrc2-1), &
         & lcol1(rsrc:rsrc+rsrc2-1), &
         & rlst, clst, work, &
         & min_width_blas)

    ! write(*,*)"work(1): ", work(1) 
    ! deallocate(rlst, clst)
    ! deallocate(work)

#if defined(SPLLT_OMP_TRACE)
    ! write(*,*)"fac_blk_id: ", fac_blk_id
    call trace_event_stop(upd_btw_id, th_id)
#endif

!$omp end task

! !$omp taskwait    
    
#else
    
    m  = a_bc%blkm
    n  = a_bc%blkn
   !id = a_bc%id
   !sa = a_bc%sa
    
    ! write(*,*) "update_between_task"
    ! call update_between(m, n, id, anode, &
    !      & n1, id1, snode, &
    !      & lfact(bcol)%lcol(sa:sa+m*n-1), &
    !      & lfact(bcol1)%lcol(csrc(1):csrc(1)+csrc(2)-1), &
    !      & lfact(bcol1)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1), &
    !      & blocks, row_list, col_list, buffer, &
    !      & control, info, st)

    csrc  = 1 + (cptr-(scol-1)*s_nb-1)*n1
    csrc2 = (cptr2 - cptr + 1)*n1

    rsrc  = 1 + (rptr-(scol-1)*s_nb-1)*n1
    rsrc2 = (rptr2 - rptr + 1)*n1

    call spllt_update_between(m, n, a_bc, dcol, anode, &
         & n1, scol, snode, &
         & lfact(bcol)%lcol(sa:sa+m*n-1), &
         & lfact(bcol1)%lcol(csrc:csrc+csrc2-1), &
         & lfact(bcol1)%lcol(rsrc:rsrc+rsrc2-1), &
         & row_list%c, col_list%c, workspace%c, &
         & min_width_blas)

#endif

    return
  end subroutine spllt_update_between_task

  ! init node
  subroutine spllt_init_node_task(fkeep, node, val, prio)
    use spllt_data_mod
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
!$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none

    type(spllt_fkeep), target, intent(in)  :: fkeep    
    type(spllt_node), intent(in) :: node
    real(wp), dimension(:), target, intent(in) :: val ! user's matrix values

     ! so that, if variable (row) i is involved in node,
     ! map(i) is set to its local row index
    ! type(spllt_fkeep), target, intent(inout) :: fkeep ! on exit, matrix a copied
     ! into relevant part of keep%lfact
    integer, optional :: prio 

    integer :: p
#if defined(SPLLT_USE_STARPU)
    type(c_ptr) :: val_c, fkeep_c
    integer(c_int) :: nval
#endif

#if defined(SPLLT_USE_OMP)
    integer :: snum
    integer :: th_id
    ! type(MA87_keep), intent(inout) :: keep ! on exit, matrix a copied
    type(spllt_fkeep), pointer :: p_fkeep
    real(wp), pointer :: p_val(:)
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)
    nval = size(val,1)
    val_c  = c_loc(val(1)) 
    fkeep_c = c_loc(fkeep)

    call spllt_insert_init_node_task_c(node%hdl, &
         & node%num, val_c, nval, fkeep_c, p)
#elif defined(SPLLT_USE_OMP)

    snum = node%num
    p_fkeep => fkeep
    p_val => val

!$omp task firstprivate(snum) firstprivate(p_fkeep, p_val) &
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(ini_nde_id) &
#endif
!$omp    & private(th_id)

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(ini_nde_id, th_id)
#endif

    call spllt_init_node(snum, p_val, p_fkeep)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(ini_nde_id, th_id)
#endif

!$omp end task

#else
    call spllt_init_node(node%num, val, fkeep)
#endif
    
    return
  end subroutine spllt_init_node_task

end module spllt_factorization_task_mod
