module spllt_factorization_task_mod
  use spllt_data_mod
  implicit none
  
contains

  ! node factorization task
  subroutine spllt_factorize_node(snode, map, fdata, keep, control)
    use spllt_data_mod
    use hsl_ma87_double
    use spllt_kernels_mod
    implicit none

    type(spllt_node_type), target, intent(inout)        :: snode ! node to factorize (spllt)    
    integer, dimension(:), allocatable, intent(inout)   ::  map
    type(spllt_data_type), target, intent(inout)        :: fdata
    type(MA87_keep), target, intent(inout)              :: keep 
    type(MA87_control), intent(in) :: control 

    type(node_type), pointer :: node ! node to factorize (hsl_ma87)
    integer :: num_nodes ! number of nodes in etree
    integer :: snum ! node idx
    integer :: prio ! task priority
    ! shortcuts
    integer :: sa ! first column 
    integer :: en ! last column
    integer :: numcol ! number of columns in node 
    integer :: numrow ! number of rows in node 
    integer :: nc ! number of block columns
    integer :: nr ! number of block rows
    integer(long) :: dblk ! diagonal block
    integer :: s_nb ! block size
    integer :: ii, jj, kk ! indexes 
    type(spllt_bc_type), pointer :: bc_kk, bc_ik, bc_jk, bc_ij ! block pointers
    integer(long) :: blk, blk1, blk2 ! block id

    ! update between
    type(node_type), pointer :: anode ! ancestor node in the atree
    type(spllt_node_type), pointer :: a_node ! ancestor node in the atree
    ! locate source blocks
    integer :: a_num ! ancestor id
    integer :: cptr  ! Position in snode of the first row 
    ! matching a column of the current block column of anode.
    integer :: cptr2  ! Position in snode of the last row 
    ! matching a column of the current block column of anode.
    logical :: map_done
    integer :: i, ilast
    type(spllt_bc_type), pointer :: bc, a_bc ! node in the atree
    integer :: cb, jb
    integer :: jlast ! Last column in the cb-th block column of anode
    integer :: k
    integer(long) :: a_dblk, a_blk ! id of block in scol containing row 

    ! write(*,*) "---------- node ----------"
    ! write(*,*) 'snode: ', snode 
    
    num_nodes = keep%info%num_nodes
    snum = snode%num 

    node => snode%node
    
    ! node priority
    ! prio = (num_nodes - snum + 1)*4
    prio = 0

    ! initialize node data
    sa = node%sa
    en = node%en
    numcol = en - sa + 1
    numrow = size(node%index)

    ! s_nb is the size of the blocks
    s_nb = node%nb

    nc = (numcol-1) / s_nb + 1 
    nr = (numrow-1) / s_nb + 1 

    ! write(*,'("numrow,numcol = ",i5,i5)')numrow, numcol
    ! write(*,'("nr,nc = ",i5,i5)')nr, nc

    ! first block in node
    dblk = node%blk_sa

    do kk = 1, nc
          
       ! A_kk
       bc_kk => fdata%bc(dblk)
       call spllt_factorize_block_task(fdata, fdata%nodes(snum), bc_kk, keep%lfact, prio+3)

       ! loop over the row blocks (that is, loop over blocks in block col)
       do ii = kk+1,nr
          ! do ii = nr,kk+1,-1             
          ! A_mk
          blk = dblk+ii-kk
          ! bc_ik => keep%blocks(blk)
          bc_ik => fdata%bc(blk)

          call spllt_solve_block_task(fdata, bc_kk, bc_ik, keep%lfact,prio+2)
       end do
       
       
       do jj = kk+1,nc

          ! L_jk
          blk2 = dblk+jj-kk
          bc_jk => fdata%bc(blk2)

          do ii = jj,nr

             ! L_ik
             blk1 = dblk+ii-kk                
             bc_ik => fdata%bc(blk1)

             ! A_ij
             ! blk = get_dest_block(keep%blocks(blk1), keep%blocks(blk2))
             blk = get_dest_block(keep%blocks(blk2), keep%blocks(blk1))
             bc_ij => fdata%bc(blk)
             call spllt_update_block_task(fdata, bc_ik, bc_jk, bc_ij, keep%lfact, prio+1)

          end do
       end do

       ! move to next block column in snode
       dblk = keep%blocks(dblk)%last_blk + 1
       ! numrow = numrow - s_nb
    end do

    ! update between
    
    !  Loop over ancestors of snode
    a_num = node%parent  
    !  Initialize cptr to correspond to the first row of the rectangular part of 
    !  the snode matrix.
    cptr = 1 + numcol
    ! write(*,*)"numrow: ", numrow
    do while(a_num.gt.0)
       anode => keep%nodes(a_num) 
       a_node => fdata%nodes(a_num)
       ! Skip columns that come from other children
       do cptr = cptr, numrow
          if(node%index(cptr).ge.anode%sa) exit
       end do
       if(cptr.gt.numrow) exit ! finished with node

       map_done = .false. ! We will only build a map when we need it

       ! Loop over affected block columns of anode
       bcols: do
          if(cptr.gt.numrow) exit
          if(node%index(cptr).gt.anode%en) exit

          ! compute local index of block column in anode and find the id of 
          ! its diagonal block
          cb = (node%index(cptr) - anode%sa)/anode%nb + 1
          a_dblk = anode%blk_sa
          do jb = 2, cb
             a_dblk = keep%blocks(a_dblk)%last_blk + 1
          end do

          ! Find cptr2
          jlast = min(anode%sa + cb*anode%nb - 1, anode%en)
          do cptr2 = cptr,numrow
             if(node%index(cptr2) > jlast) exit
          end do
          cptr2 = cptr2 - 1
          ! write(*,*)"j: ", nodes(snode)%index(cptr), ", jlast: ", nodes(snode)%index(cptr2)
          ! Set info for source block csrc (hold start location
          ! and number of entries)

          ! csrc(1) = 1 + (cptr-(kk-1)*node%nb-1)*blocks(bc_kk%blk%id)%blkn
          ! csrc(2) = (cptr2 - cptr + 1)*blocks(bc_kk%blk%id)%blkn

          ! write(*,*)"csrc(1): ", csrc(1), ", csrc(2): ", csrc(2)
          ! Build a map of anode's blocks if this is first for anode
          if(.not.map_done) call spllt_build_rowmap(anode, map) 
          ! write(*,*)"csrc(1): ", csrc(1), ", csrc(2): ", csrc(2)

          ! Loop over the blocks of snode
          ii = map(node%index(cptr)) 
          ! ii = -1 
          ilast = cptr ! Set start of current block

          do i = cptr, numrow
             k = map(node%index(i))

             if(k.ne.ii) then
                a_blk = a_dblk + ii - cb
                ! a_bc => keep%blocks(a_blk)
                a_bc => fdata%bc(a_blk)

                ! rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
                ! rsrc(2) = (i - ilast)*blocks(blk)%blkn

                dblk = node%blk_sa
                ! Loop over the block columns in node. 
                do kk = 1, nc

                   bc_kk => fdata%bc(dblk)

                   call spllt_update_between_task( &
                        & fdata, &
                                ! & bc, &
                        & bc_kk, &
                        & node, a_bc, a_node, &
                                ! & csrc, rsrc, &
                        & cptr, cptr2, ilast, i-1, &
                        & fdata%row_list, fdata%col_list, fdata%workspace, &
                        & keep%lfact, keep%blocks, fdata%bc, &
                        & control, prio)

                   dblk = keep%blocks(dblk)%last_blk + 1
                end do

                ii = k
                ilast = i ! Update start of current block
             end if
          end do
          ! #if defined(SPLLT_USE_OMP)
          ! !$omp taskwait
          ! #endif
          ! i = min(i,numrow)
          a_blk = a_dblk + ii - cb
          ! a_bc => keep%blocks(a_blk)
          a_bc => fdata%bc(a_blk)

          ! rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
          ! rsrc(2) = (i - ilast)*blocks(blk)%blkn

          dblk = node%blk_sa
          ! Loop over the block columns in node. 
          do kk = 1, nc

             bc_kk => fdata%bc(dblk)

             call spllt_update_between_task( &
                  & fdata, &
                                ! & bc, &
                  & bc_kk, &
                  & node, a_bc, a_node, &
                                ! & csrc, rsrc, &
                  & cptr, cptr2, ilast, i-1, &
                  & fdata%row_list, fdata%col_list, fdata%workspace, &
                  & keep%lfact, keep%blocks, fdata%bc, &
                  & control, prio)

             dblk = keep%blocks(dblk)%last_blk + 1
          end do

          ! Move cptr down, ready for next block column of anode
          cptr = cptr2 + 1
       end do bcols

       a_num = anode%parent
    end do

    return
  end subroutine spllt_factorize_node

#if defined(SPLLT_USE_STARPU)
  ! deinitialize factorization
  ! StarPU: unregister data handles (block handles) in StarPU
  subroutine spllt_deinit_task(keep, pbl)
    use hsl_ma87_double
    use  starpu_f_mod
    implicit none

    type(MA87_keep), target, intent(inout) :: keep 
    type(spllt_data_type), intent(inout) :: pbl

    integer :: i
    integer :: snode, num_nodes
    type(node_type), pointer :: node ! node in the atree    
    integer(long) :: blk, dblk
    integer :: l_nb, sz, sa, en

    num_nodes = keep%info%num_nodes
    blk = 1

    do snode = 1, num_nodes

       node => keep%nodes(snode)

       l_nb = node%nb
       sz = (size(node%index) - 1) / l_nb + 1
       sa = node%sa
       en = node%en

       do i = sa, en, l_nb
          dblk = blk
          do blk = dblk, dblk+sz-1
             ! write(*,*) 'unregister blk: ', blk
             call starpu_f_data_unregister_submit(pbl%bc(blk)%hdl)
          end do
          sz = sz - 1
       end do
    end do

  end subroutine spllt_deinit_task

  subroutine spllt_data_unregister(keep, pbl)
    use hsl_ma87_double
    use  starpu_f_mod
    implicit none

    type(MA87_keep), target, intent(inout) :: keep 
    type(spllt_data_type), intent(inout) :: pbl

    integer :: i
    integer :: snode, num_nodes
    type(node_type), pointer :: node ! node in the atree    
    integer(long) :: blk, dblk
    integer :: l_nb, sz, sa, en

    num_nodes = keep%info%num_nodes
    blk = 1

    do snode = 1, num_nodes

       node => keep%nodes(snode)

       l_nb = node%nb
       sz = (size(node%index) - 1) / l_nb + 1
       sa = node%sa
       en = node%en

       do i = sa, en, l_nb
          dblk = blk
          do blk = dblk, dblk+sz-1
             ! write(*,*) 'unregister blk: ', blk
             call starpu_f_data_unregister(pbl%bc(blk)%hdl)
          end do
          sz = sz - 1
       end do
    end do

  end subroutine spllt_data_unregister
#endif

  ! factorize block 
  ! _potrf
  subroutine spllt_factorize_block_task(fdata, node, bc, lfact, prio)
    use spllt_data_mod
    use hsl_ma87_double
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
    
    type(spllt_data_type), target, intent(in)  :: fdata
    type(spllt_node_type), intent(in)          :: node
#if defined(SPLLT_USE_OMP)
    type(spllt_bc_type), pointer, intent(inout) :: bc ! block to be factorized    
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
#else
    type(spllt_bc_type), target, intent(inout) :: bc ! block to be factorized    
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    integer, optional :: prio

    integer :: m, n, bcol, sa
    integer(long) :: id
    type(block_type), pointer :: blk ! block to be factorized
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
    node_hdl = c_null_ptr
    if (node%node%blk_sa .eq. bc%blk%id) then
       node_hdl = node%hdl
       ! write(*,*)'Test'
    end if
    call spllt_starpu_insert_factorize_block_c(bc%hdl, node_hdl, p)

    ! call spllt_starpu_insert_factorize_block(bc, p)
#elif defined(SPLLT_USE_OMP)

    blk => bc%blk
    bcol = blk%bcol
    lcol => lfact(bcol)%lcol

    m    = blk%blkm
    n    = blk%blkn
    sa   = blk%sa
    id   = blk%id
    bc_c => fdata%bc(id)%c

! !$omp taskwait

!$omp task firstprivate(m, n, sa) &
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(fac_blk_id) &
#endif
!$omp    & firstprivate(blk, bc_c, bcol, lcol) &
! !$omp    & depend(inout:lcol(sa:sa+n*m-1))
! !$omp    & depend(inout:fdata%bc(id)%c(:))
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

    blk => bc%blk

    m    = blk%blkm
    n    = blk%blkn
    id   = blk%id 
    bcol = blk%bcol

    sa   = blk%sa

    ! factorize_block
    ! call factor_diag_block(n, m, id, &
    !      & lfact(bcol)%lcol(sa:sa+n*m-1),   &
    !      & control, flag, detlog(0))

    call spllt_factor_diag_block(m, n, lfact(bcol)%lcol(sa:sa+n*m-1))
#endif
    
    return
  end subroutine spllt_factorize_block_task

  ! _trsm
  subroutine spllt_solve_block_task(fdata, bc_kk, bc_ik, lfact, prio)
    use spllt_data_mod
    use hsl_ma87_double
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

    type(spllt_data_type), target, intent(inout)  :: fdata
#if defined(SPLLT_USE_OMP)
    type(spllt_bc_type), pointer, intent(inout) :: bc_kk, bc_ik ! block to be factorized    
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
    type(block_type), pointer :: blk_kk, blk_ik ! block to be factorized
#else
    type(spllt_bc_type), intent(inout) :: bc_kk, bc_ik ! block to be factorized    
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
    

    blk_kk => bc_kk%blk
    blk_ik => bc_ik%blk
    ! bcol is block column that blk and dblk belong to
    bcol = blk_kk%bcol    
    ! write(*,*)"bcol: ", bcol
    lcol => lfact(bcol)%lcol

    ! bc_kk
    d_m  = blk_kk%blkm
    d_n  = blk_kk%blkn
    d_sa = blk_kk%sa
    d_id = blk_kk%id

    ! bc_ik
    n  = blk_ik%blkn
    m  = blk_ik%blkm
    sa = blk_ik%sa
    id = blk_ik%id
    ! write(*,*)"solve id: ", id

    bc_kk_c => fdata%bc(d_id)%c
    bc_ik_c => fdata%bc(id)%c

    snode = blk_kk%node
    blk_sa = fdata%nodes(snode)%node%blk_sa
    ! write(*,*)"blk_sa: ", blk_sa
! !$omp taskwait

!$omp task private(th_id) &
!$omp    & firstprivate(m, n, sa, d_m, d_n, d_sa) &
!$omp    & firstprivate(lcol, bcol) &
!$omp    & firstprivate(blk_kk, blk_ik) &
!$omp    & firstprivate(bc_kk_c, bc_ik_c) &
!$omp    & firstprivate(d_id, id) &
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(slv_blk_id) &
#endif
!$omp    & depend(in:bc_kk_c(1)) &
!$omp    & depend(inout:bc_ik_c(1))

! !$omp    & depend(in:fdata%bc(d_id)%c(1)) &
! !$omp    & depend(inout:fdata%bc(id)%c)

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
    d_m  = bc_kk%blk%blkm
    d_n  = bc_kk%blk%blkn
    d_sa = bc_kk%blk%sa
    d_id = bc_kk%blk%id

    ! bc_ik
    n  = bc_ik%blk%blkn
    m  = bc_ik%blk%blkm
    sa = bc_ik%blk%sa
    id = bc_ik%blk%id
    
    ! bcol is block column that blk and dblk belong to
    bcol = bc_kk%blk%bcol    

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
  subroutine spllt_update_block_task(fdata, bc_ik, bc_jk, bc_ij, lfact, prio)
    use spllt_data_mod
    use hsl_ma87_double
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
    
    type(spllt_data_type), target, intent(in)  :: fdata
    ! type(block_type), intent(inout) :: bc_ik, bc_jk, bc_ij ! block to be updated    
#if defined(SPLLT_USE_OMP)
    type(spllt_bc_type), pointer, intent(inout) :: bc_ik, bc_jk, bc_ij
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
#else
    type(spllt_bc_type), intent(inout) :: bc_ik, bc_jk, bc_ij
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    integer, optional :: prio 

    integer :: n1, m1, sa1, n2, m2, sa2, n, m, sa
    integer :: bcol1, bcol, bcol2
    integer :: p, d
#if defined(SPLLT_USE_OMP)
    logical :: is_diag
    real(wp), dimension(:), pointer :: lcol1, lcol2, lcol
    real(wp), dimension(:), pointer :: bc_ik_c, bc_jk_c, bc_ij_c
    type(block_type), pointer :: blk_ik, blk_jk, blk_ij ! block to be factorized
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
    if (bc_ij%blk%dblk .eq. bc_ij%blk%id) then
       d = 1
    else
       d = 0
    end if

    call spllt_starpu_insert_update_block_c(bc_ik%hdl, bc_jk%hdl, bc_ij%hdl, d, p)
#elif defined(SPLLT_USE_OMP)

    blk_ik => bc_ik%blk
    blk_jk => bc_jk%blk
    blk_ij => bc_ij%blk

    id_ik = blk_ik%id
    id_jk = blk_jk%id
    id_ij = blk_ij%id

    bc_ik_c => fdata%bc(id_ik)%c
    bc_jk_c => fdata%bc(id_jk)%c
    bc_ij_c => fdata%bc(id_ij)%c

    bcol1 = bc_ik%blk%bcol
    lcol1 => lfact(bcol1)%lcol

    bcol2 = bc_jk%blk%bcol
    lcol2 => lfact(bcol2)%lcol

    bcol = bc_ij%blk%bcol
    lcol => lfact(bcol)%lcol

    snode = blk_ij%node
    blk_sa = fdata%nodes(snode)%node%blk_sa
    
! !$omp taskwait

!$omp task private(n1, m1, sa1, n2, m2, sa2, &
!$omp    & n, m, sa, is_diag, th_id) &
!$omp    & firstprivate(bcol2, lcol2, bcol1, lcol1, bcol, lcol) &
!$omp    & firstprivate(bc_ik, bc_jk, bc_ij) &
!$omp    & firstprivate(bc_ik_c, bc_jk_c, bc_ij_c) &
!$omp    & firstprivate(blk_ik, blk_jk, blk_ij) &
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(upd_blk_id) &
#endif
!$omp    & depend(in:bc_ik_c(1), bc_jk_c(1)) &
!$omp    & depend(inout: bc_ij_c(1))

! !$omp    & depend(in:fdata%bc(id_ik)%c, fdata%bc(id_jk)%c) &
! !$omp    & depend(inout: fdata%bc(id_ij)%c)

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
    n1  = blk_ik%blkn
    m1  = blk_ik%blkm
    sa1 = blk_ik%sa

    ! bc_jk
    n2  = blk_jk%blkn
    m2  = blk_jk%blkm
    sa2 = blk_jk%sa

    ! bc_ij
    n  = blk_ij%blkn
    m  = blk_ij%blkm
    sa = blk_ij%sa

    is_diag = blk_ij%dblk.eq.blk_ij%id 

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
    n1  = bc_ik%blk%blkn
    m1  = bc_ik%blk%blkm
    sa1 = bc_ik%blk%sa
    bcol1 = bc_ik%blk%bcol
    
    ! bc_jk
    n2  = bc_jk%blk%blkn
    m2  = bc_jk%blk%blkm
    sa2 = bc_jk%blk%sa
    bcol2 = bc_jk%blk%bcol
    
    ! bc_ij
    n  = bc_ij%blk%blkn
    m  = bc_ij%blk%blkm
    sa = bc_ij%blk%sa
    
    bcol = bc_ij%blk%bcol

    ! write(*,*)"bc_ik id: ", bc_ik%id, ", m: ", m1, ", n: ", n1
    ! write(*,*)"bc_jk id: ", bc_jk%id, ", m: ", m2, ", n: ", n2
    ! write(*,*)"bc_ij id: ", bc_ij%id, ", m: ", m, ", n: ", n
    
    ! write(*,*) "size lfact(bcol1)%lcol: ", size(lfact(bcol1)%lcol)
    ! write(*,*) "sa2+n1*m2-1: ", sa2+n1*m2-1

    call spllt_update_block(m, n, &
         & lfact(bcol)%lcol(sa:sa+n*m-1), &
         & bc_ij%blk%dblk.eq.bc_ij%blk%id, n1, &
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
!     type(node_type), intent(in) :: dnode ! Node to which blk belongs
!     integer :: n1 ! number of columns in source block column
!     ! integer(long), intent(in) :: src  ! identifier of block in source block col
!     integer, intent(in) :: scol ! index of block column that src belongs to in snode
!     type(node_type), intent(in) :: snode ! Node to which src belongs
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
  subroutine spllt_update_between_task(fdata, &
       & dbc, snode, a_bc, anode, &
       ! & csrc, rsrc, &
       & cptr, cptr2, rptr, rptr2, &
       & row_list, col_list, workspace, &
       & lfact, blocks, bcs, &
       & control, prio)
    use spllt_data_mod
    use hsl_ma87_double
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

    type(spllt_data_type), target, intent(in)       :: fdata
#if defined(SPLLT_USE_OMP)
    type(spllt_bc_type), pointer, intent(inout)     :: a_bc ! dest block
    type(spllt_bc_type), pointer, intent(inout)     :: dbc ! diag block in source node

    type(spllt_bc_type), allocatable, target        :: workspace(:)
    type(spllt_bc_type), target                     :: bcs(:) ! block info.
    type(spllt_workspace_i), allocatable, target    :: row_list(:), col_list(:)
#else
    type(spllt_bc_type), intent(inout)              :: a_bc ! dest block
    type(spllt_bc_type), intent(in)                 :: dbc ! diag block in source node

    type(spllt_bc_type)                             :: workspace
    type(spllt_bc_type)                             :: bcs(:) ! block info.
    type(spllt_workspace_i)                         :: row_list, col_list
#endif

#if defined(SPLLT_USE_STARPU)
    type(node_type), target                         :: snode ! src node
    type(spllt_node_type), target                   :: anode ! dest node
#elif defined(SPLLT_USE_OMP)
    type(node_type), pointer, intent(in)            :: snode ! src node
    type(spllt_node_type), pointer, intent(in)      :: anode ! dest node
#else
    type(node_type)                                 :: snode ! src node
    type(spllt_node_type)                           :: anode ! dest node
#endif
    integer                                         :: cptr, cptr2, rptr, rptr2 
!    integer :: csrc(2), rsrc(2) ! used for update_between tasks to


    ! real(wp), dimension(:), allocatable :: buffer ! update_buffer workspace
    type(block_type), dimension(:)                  :: blocks ! block info. 

#if defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)
    type(lfactor), allocatable, target, intent(inout)  :: lfact(:)
#else
    type(lfactor), allocatable, intent(inout)          :: lfact(:)
#endif
    type(MA87_control), intent(in) :: control     
    ! integer, intent(out) :: st ! TODO error managment
    integer, optional :: prio 

    ! integer :: info ! TODO error managment
    
    integer :: p
    integer :: m, n, n1, sa
    integer :: bcol, bcol1
    integer(long) :: id, id1
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
    type(c_ptr) :: ljk_hdl, lik_hdl
    type(c_ptr) :: snode_c, anode_c, a_bc_c
#endif

#if defined(SPLLT_USE_OMP)
    real(wp), dimension(:), pointer :: lcol1, lcol
    integer :: th_id
    integer :: csrc_sa, rsrc_sa, csrc_en, rsrc_en
    ! type(spllt_bc_type), pointer :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    real(wp), dimension(:), pointer :: a_bc_c, dbc_c
    ! real(wp), dimension(:), pointer :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    ! type(spllt_bc_type), pointer    :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    real(wp), dimension(:), pointer    :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    type(block_type), pointer :: blk_kk, a_blk
    type(spllt_bc_type), pointer :: p_workspace(:) => null()
    type(spllt_workspace_i), pointer :: p_rlst(:) => null(), p_clst(:) => null()
    ! real(wp), dimension(:), allocatable :: work
    real(wp), dimension(:), pointer :: work
    integer, dimension(:), pointer :: rlst, clst
    integer :: min_width_blas
    type(node_type), pointer                :: p_snode ! src node
    type(spllt_node_type), pointer          :: p_anode ! dest node
    integer(long) :: id_ik_sa, id_ik_en, id_jk_sa, id_jk_en, id_ij
    real(wp), dimension(:), pointer :: lcol_ik, lcol_jk
    type(lfactor), pointer :: p_lfact(:)
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

    s_nb = snode%nb    ! block size in source node
    n1  = dbc%blk%blkn ! width of column

    bcol1 = dbc%blk%bcol
    scol = bcol1 - blocks(snode%blk_sa)%bcol + 1

    bcol = a_bc%blk%bcol
    dcol = bcol - blocks(anode%node%blk_sa)%bcol + 1
    ! dcol = a_bc%blk%bcol ! blocks(anode%node%blk_sa)%bcol + 1

#if defined(SPLLT_USE_STARPU)
! #if 0
    ! write(*,*)"associated blk: ", associated(a_bc%blk)
    dblk = dbc%blk%id

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
       ljk_m = ljk_m + bcs(blk)%blk%blkm
    end do
    ljk_n = dbc%blk%blkn

#if defined(SPLLT_USE_GPU)
    ! create temporary handle for gathering blocks in Ljk factor
    call starpu_matrix_data_register(ljk_hdl, fdata%bc(blk_sa)%mem_node, &
         & c_loc(fdata%bc(blk_sa)%c(1)), ljk_m, ljk_m, ljk_n, &
         & int(wp,kind=c_size_t))
#endif

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
       lik_m = lik_m + bcs(blk)%blk%blkm
    end do

#if defined(SPLLT_USE_GPU)
    ! create temporary handle for gathering blocks in Ljk factor
    call starpu_matrix_data_register(lik_hdl, fdata%bc(blk_sa)%mem_node, &
         & c_loc(fdata%bc(blk_sa)%c(1)), lik_m, lik_m, ljk_n, &
         & int(wp,kind=c_size_t))
#endif
    
    csrc  = 1 + (mod(cptr-1, s_nb))*n1
    csrc2 = (cptr2 - cptr + 1)*n1

    rsrc  = 1 + (mod(rptr-1, s_nb))*n1
    rsrc2 = (rptr2 - rptr + 1)*n1

    snode_c = c_loc(snode)
    anode_c = c_loc(anode%node)
    a_bc_c  = c_loc(a_bc%blk)
    
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
         & control%min_width_blas, &
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
         & control%min_width_blas, &
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
    a_blk => a_bc%blk

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
    blk_kk => dbc%blk
    dblk = blk_kk%id

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
    bc_jk_sa => fdata%bc(id_jk_sa)%c ! bcs(blk_sa)%c
    bc_jk_en => fdata%bc(id_jk_en)%c ! bcs(blk_en)%c
    ! write(*,*) "nb_blk: ", nb_blk
    ! write(*,*) "nb_blk:", nb_blk, "blk_sa: ", blk_sa, ", blk_en: ", blk_en
    csrc_sa = blocks(blk_sa)%sa
    csrc_en = blocks(blk_en)%sa

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

    bc_ik_sa => fdata%bc(id_ik_sa)%c ! bcs(blk_sa)%c
    bc_ik_en => fdata%bc(id_ik_en)%c ! bcs(blk_en)%c

    rsrc_sa = blocks(blk_sa)%sa
    rsrc_en = blocks(blk_en)%sa

    id_ij = a_blk%id
    a_bc_c => fdata%bc(id_ij)%c

    dbc_c => dbc%c

    bcol1 = blk_kk%bcol
    bcol  = a_bc%blk%bcol

    lcol1 => lfact(bcol1)%lcol
    lcol  => lfact(bcol)%lcol

    p_workspace => workspace

    p_rlst => row_list
    p_clst => col_list

    min_width_blas = control%min_width_blas

    ! p_snode => snode
    ! p_anode => anode

    p_snode => fdata%nodes(blk_kk%node)%node
    p_anode => fdata%nodes(a_blk%node)

    s_nb = p_snode%nb    ! block size in source node
    n1  = blk_kk%blkn ! width of column

    csrc  = 1 + (cptr-(scol-1)*s_nb-1)*n1
    csrc2 = (cptr2 - cptr + 1)*n1

    rsrc  = 1 + (rptr-(scol-1)*s_nb-1)*n1
    rsrc2 = (rptr2 - rptr + 1)*n1

    lcol_ik => lcol1(rsrc:rsrc+rsrc2-1)
    lcol_jk => lcol1(csrc:csrc+csrc2-1)

    m  = a_blk%blkm
    n  = a_blk%blkn
    sa = a_blk%sa

    blk_sa = snode%blk_sa

! !$omp taskwait    

!$omp task firstprivate(m, n, a_blk, dcol, p_anode, n1, scol, p_snode, &
!$omp    & lcol, lcol1, min_width_blas, sa, csrc, csrc2, rsrc, rsrc2) &
!$omp    & private(rlst, clst, work) & 
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(upd_btw_id) &
#endif
!$omp    & firstprivate(bc_ik_sa, bc_ik_en, bc_jk_sa, bc_jk_en, a_bc_c) &
!$omp    & depend(in:bc_ik_sa(1), bc_ik_en(1)) &
!$omp    & depend(in:bc_jk_sa(1), bc_jk_en(1)) &
!$omp    & depend(inout: a_bc_c(1))

! !$omp    & depend(in:fdata%bc(id_ik_sa)%c, fdata%bc(id_ik_en)%c) &
! !$omp    & depend(in:fdata%bc(id_jk_sa)%c, fdata%bc(id_jk_en)%c) &
! !$omp    & depend(inout:fdata%bc(id_ij)%c)

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

    ! call spllt_omp_update_between_cpu_func(m, n, a_blk, dcol, p_anode%node, &
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

    call spllt_update_between(m, n, a_blk, dcol, p_anode%node, &
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
    
    m  = a_bc%blk%blkm
    n  = a_bc%blk%blkn
    id = a_bc%blk%id
    sa = a_bc%blk%sa
    
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

    call spllt_update_between(m, n, a_bc%blk, dcol, anode%node, &
         & n1, scol, snode, &
         & lfact(bcol)%lcol(sa:sa+m*n-1), &
         & lfact(bcol1)%lcol(csrc:csrc+csrc2-1), &
         & lfact(bcol1)%lcol(rsrc:rsrc+rsrc2-1), &
         & row_list%c, col_list%c, workspace%c, &
         & control%min_width_blas)

#endif

    return
  end subroutine spllt_update_between_task

  ! init node
  subroutine spllt_init_node_task(fdata, node, val, keep, prio)
    use spllt_data_mod
    use hsl_ma87_double
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

    type(spllt_data_type), target, intent(in)  :: fdata    
    type(spllt_node_type), intent(in) :: node
    real(wp), dimension(:), target, intent(in) :: val ! user's matrix values

     ! so that, if variable (row) i is involved in node,
     ! map(i) is set to its local row index
    type(MA87_keep), target, intent(inout) :: keep ! on exit, matrix a copied
     ! into relevant part of keep%lfact
    integer, optional :: prio 

    integer :: p
#if defined(SPLLT_USE_STARPU)
    type(c_ptr) :: val_c, keep_c
    integer(c_int) :: nval
#endif

#if defined(SPLLT_USE_OMP)
    integer :: snum
    integer :: th_id
    ! type(MA87_keep), intent(inout) :: keep ! on exit, matrix a copied
    type(MA87_keep), pointer :: p_keep 
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
    keep_c = c_loc(keep)

    call spllt_insert_init_node_task_c(node%hdl, &
         & node%num, val_c, nval, keep_c, p)
#elif defined(SPLLT_USE_OMP)

    snum = node%num
    p_keep => keep
    p_val => val

!$omp task firstprivate(snum) firstprivate(p_keep, p_val) &
#if defined(SPLLT_OMP_TRACE)
!$omp    & shared(ini_nde_id) &
#endif
!$omp    & private(th_id)

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(ini_nde_id, th_id)
#endif

    call spllt_init_node(snum, p_val, p_keep)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(ini_nde_id, th_id)
#endif

!$omp end task

#else
    call spllt_init_node(node%num, val, keep)
#endif
    
    return
  end subroutine spllt_init_node_task

end module spllt_factorization_task_mod
