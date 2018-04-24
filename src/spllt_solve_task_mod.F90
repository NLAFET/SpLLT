module spllt_solve_task_mod
contains
  
  !*************************************************  
  !
  ! Forward solve with block on diagonal
  !
  subroutine spllt_solve_fwd_block_task_with_ftask(dblk, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id, task_manager)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    use spllt_solve_dep_mod
    use task_manager_omp_mod
    use timer_mod
    implicit none
    
    integer,                    intent(in)    :: dblk !Index of diagonal block
    integer,                    intent(in)    :: nrhs !Number of RHS
    integer,                    intent(in)    :: ldr  !Leading dimension of RHS
    real(wp), target,           intent(inout) :: upd(:, :)
    real(wp), target,           intent(inout) :: rhs(ldr * nrhs)
    real(wp), target,           intent(inout) :: xlocal(:, :)
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: trace_id
    type(task_manager_omp_t),  intent(inout)  :: task_manager
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: node
    integer                     :: dep
    integer                     :: i, j, r
    integer                     :: ndep
    integer                     :: nthread, threadID
    integer,           pointer  :: p_index(:)
    real(wp),          pointer  :: p_lcol(:)
    real(wp),          pointer  :: p_lcol_update(:)
    type(spllt_block), pointer  :: p_blk_dep_update

    real(wp),          pointer  :: p_upd(:,:)
    real(wp),          pointer  :: p_xlocal(:,:)
    real(wp),          pointer  :: p_rhs(:)
    integer,           pointer  :: p_dep(:)
    integer                     :: nftask

    type(spllt_block), pointer  :: p_bc(:)
    type(spllt_timer_t), save   :: timer
        
    nthread   = omp_get_num_threads()
    call spllt_open_timer(task_manager%workerID, &
      "spllt_solve_fwd_block_task_with_ftask", timer)

    ! Get block info
    node      = fkeep%bc(dblk)%node
    m         = fkeep%bc(dblk)%blkm
    n         = fkeep%bc(dblk)%blkn
    sa        = fkeep%bc(dblk)%sa
    bcol      = fkeep%bc(dblk)%bcol ! Current block column
    dcol      = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col       = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
    p_bc      => fkeep%bc
    
    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs
    
    nftask      = 0

    p_dep => fkeep%bc(dblk)%fwd_update_dep
    ndep  = size(p_dep)

    do dep = 1, ndep - 1

      !$omp task firstprivate(dep)                       &
      !$omp firstprivate(p_bc, p_dep)                    &
      !$omp private(threadID)                            &
      !$omp depend(in: p_bc(p_dep(dep)))                 &
      !$omp depend(inout: p_dep(1))

      threadID = omp_get_thread_num()
#if defined(SPLLT_VERBOSE)
      print '(a, i3, a, i3)', "SLV Fake Task dep ", dblk, " [in : ", &
        p_dep(dep)
#endif

      !$omp end task

    end do
    nftask = merge(ndep - 1,0, ndep .gt. 1)

    !$omp task firstprivate(m, n, col, offset)                      &
    !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)              &
    !$omp firstprivate(dblk)                                        &
    !$omp firstprivate(nthread)                                     &
    !$omp firstprivate(p_blk_dep_update, p_lcol_update)             &
    !$omp firstprivate(p_dep, dep, p_bc)                            &
    !$omp firstprivate(p_upd, p_rhs, p_xlocal)                      &
    !$omp private(threadID, r, j, i)                                &
    !$omp depend(in: p_bc(p_dep(dep)))                              &
    !$omp depend(in: p_dep(1))                                      &
    !$omp depend(inout: p_bc(dblk))                                  

#if defined(SPLLT_PROFILING_TRACE)
    call trace_event_start(trace_id, omp_get_thread_num())
#endif

    threadID  = omp_get_thread_num()
#if defined(SPLLT_VERBOSE)
    print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
    print *, p_dep
#endif

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          p_rhs(i)    = p_rhs(i) + p_upd(i, j)
          p_upd(i,j)  = zero ! Reset in case of bwd solve
        end do
      end do
    end do


    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa:sa+n*n-1),    &
         'Transpose    ', 'Non-unit', nrhs, p_rhs, ldr)
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m .gt. 0) then
       sa = sa + n * n
       call slv_fwd_update(m, n, col, offset, p_index,                &
            p_lcol(sa : sa + n * m - 1), n, nrhs,                     &
            p_upd(:, threadID + 1), ldr, p_rhs,                       &
            ldr, p_xlocal(:, threadID + 1))
    endif

#if defined(SPLLT_PROFILING_TRACE)
    call trace_event_stop (trace_id, omp_get_thread_num())
#endif

    !$omp end task

    call spllt_close_timer(task_manager%workerID, timer)
    call task_manager%ntask_submitted(1, nftask)
  end subroutine spllt_solve_fwd_block_task_with_ftask



  subroutine spllt_solve_fwd_update_task_with_ftask(blk, node, nrhs, upd, rhs, &
      ldr, xlocal, fkeep, trace_id, task_manager)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num
    use trace_mod
    use spllt_solve_dep_mod
   !use utils_mod, ONLY : spllt_update_omp_task_info
    use task_manager_omp_mod
    use timer_mod
    implicit none
    
    integer,                    intent(in)    :: blk  ! Index of block
    integer,                    intent(in)    :: node
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp), target,           intent(inout) :: upd(:,:)        
    real(wp), target,           intent(in)    :: rhs(ldr*nrhs)
    real(wp), target,           intent(out)   :: xlocal(:,:)
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer, intent(in)                       :: trace_id
    type(task_manager_omp_t),   intent(inout) :: task_manager

    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID
    integer                     :: dep, ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    integer, pointer            :: p_dep(:)
    integer                     :: blk_dep_solve
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer                     :: nftask

    type(spllt_timer_t), save   :: timer
        
    call spllt_open_timer(task_manager%workerID, &
      "spllt_solve_fwd_update_task_with_ftask", timer)

    ! Establish variables describing block
    n         = fkeep%bc(blk)%blkn
    m         = fkeep%bc(blk)%blkm
    blk_sa    = fkeep%bc(blk)%sa
    bcol      = fkeep%bc(blk)%bcol
    dcol      = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col       = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
    offset    = offset + (blk-fkeep%bc(blk)%dblk) &
      * fkeep%nodes(node)%nb ! this blk
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
    p_bc      => fkeep%bc

    p_upd         => upd
    p_xlocal      => xlocal
    p_rhs         => rhs
    p_dep         => fkeep%bc(blk)%fwd_update_dep
    blk_dep_solve = fkeep%bc(blk)%fwd_solve_dep
    ndep          = size(p_dep)
    nftask        = 0

    do dep = 1, ndep - 1

      !$omp task firstprivate(p_dep, p_bc, dep)                 &
      !$omp private(threadID)                                   &
      !$omp depend(in: p_bc(p_dep(dep)))                        &
      !$omp depend(inout: p_dep(1))

      threadID = omp_get_thread_num()
#if defined(SPLLT_VERBOSE)
      print '(a, i3, a, i3)', "UPD Fake Task dep of ", blk, &
        " [in : ", p_dep(dep)
#endif

      !$omp end task

    end do
    nftask = merge(ndep - 1, 0, ndep .gt. 1)

    !$omp task firstprivate(m, n, col, offset)                    &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(blk)                                       &
    !$omp firstprivate(p_bc, p_dep, dep, blk_dep_solve)           &
    !$omp private(threadID)                                       &
    !$omp depend(in: p_bc(p_dep(dep)))                            &
    !$omp depend(in: p_bc(blk_dep_solve))                         &
    !$omp depend(in: p_dep(1))                                    & !fixme ?
    !$omp depend(inout: p_bc(blk))

#if defined(SPLLT_PROFILING_TRACE)
    call trace_event_start(trace_id, omp_get_thread_num())
#endif
    
    threadID  = omp_get_thread_num()
#if defined(SPLLT_VERBOSE)
    print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
    print *, p_dep
    print *, blk_dep_solve
#endif 

    call slv_fwd_update(m, n, col, offset, p_index,         &
      p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
      p_upd(:, threadID + 1), ldr, p_rhs,                   &
      ldr, p_xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_TRACE)
    call trace_event_stop (trace_id, omp_get_thread_num())
#endif

    !$omp end task

    call spllt_close_timer(task_manager%workerID, timer)
    call task_manager%ntask_submitted(1, nftask)

  end subroutine spllt_solve_fwd_update_task_with_ftask

  !*************************************************  
  !
  ! Backward solve with block on diagonal
  !         
  subroutine spllt_solve_bwd_block_task_with_ftask(dblk, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id, task_manager)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    use spllt_solve_dep_mod
    use task_manager_omp_mod
    use timer_mod
    use utils_mod
    implicit none

    integer, intent(in)                       :: dblk ! Index of diagonal block
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:, :)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(task_manager_omp_t), intent(inout)   :: task_manager
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: bcol, col
    integer                     :: offset
    integer                     :: node
    integer                     :: dep
    integer                     :: i, j, r, chunkth
    integer                     :: threadID, nthread
    integer                     :: ndep
    integer                     :: blk_dep_update
    integer,            pointer :: p_index(:)
    real(wp),           pointer :: p_lcol(:)
    real(wp),           pointer :: p_lcol_update(:)
    real(wp),           pointer :: p_lcol_solve(:)
    type(spllt_block),  pointer :: p_blk_dep_update
    type(spllt_block),  pointer :: p_blk_dep_solve

    real(wp),           pointer :: p_upd(:,:)
    real(wp),           pointer :: p_xlocal(:,:)
    real(wp),           pointer :: p_rhs(:)
    integer,            pointer :: p_dep_solve(:)
    integer,            pointer :: p_dep(:)
    integer                     :: nftask

    type(spllt_block), pointer  :: p_bc(:)
    type(spllt_timer_t), save   :: timer

    nthread = omp_get_num_threads()
    call spllt_open_timer(task_manager%workerID, &
      "spllt_solve_bwd_block_task_with_ftask", timer)

    ! Get block info
    node      = fkeep%bc(dblk)%node
    n         = fkeep%bc(dblk)%blkn
    m         = fkeep%bc(dblk)%blkm
    sa        = fkeep%bc(dblk)%sa
    bcol      = fkeep%bc(dblk)%bcol ! Current block column
    col       = calc_col(fkeep%nodes(node), fkeep%bc(dblk)) ! current bcol
    col       = fkeep%nodes(node)%sa + (col-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
    p_bc      => fkeep%bc

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

    nftask      = 0

    p_dep_solve => fkeep%bc(dblk)%bwd_solve_dep
    ndep        = size(p_dep_solve)

    do dep = 1, ndep - 1

      p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol
      p_blk_dep_solve => fkeep%bc(p_dep_solve(dep))

      !$omp task firstprivate(p_lcol_solve, p_blk_dep_solve)        &
      !$omp firstprivate(p_dep_solve)                               &
      !$omp firstprivate(dblk, dep)                                 &
      !$omp private(threadID)                                       &
      !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
      !$omp depend(inout: p_dep_solve(1))

      threadID = omp_get_thread_num()
  !   print '(a, i3, a, i3)', "SLV Fake Task dep ", dblk, " [in : ",&
  !     p_dep_solve(dep)

      !$omp end task

    end do
    nftask = merge(ndep - 1, 0, ndep .gt. 1)

    blk_dep_update  = fkeep%bc(dblk)%bwd_update_dep
    p_lcol_update   => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
    p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol

    p_blk_dep_update  => fkeep%bc(blk_dep_update)   
    p_blk_dep_solve   => fkeep%bc(p_dep_solve(dep))

    !$omp task firstprivate(m, n, col, offset, node, bcol)        &
    !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)            &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(p_dep_solve, dblk)                         &
    !$omp firstprivate(p_lcol_update, p_blk_dep_update)           &
    !$omp firstprivate(p_lcol_solve, p_blk_dep_solve)             &
    !$omp firstprivate(blk_dep_update)                            &
    !$omp firstprivate(nthread)                                   &
    !$omp private(threadID, r, j, i)                              &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
    !$omp depend(in: p_dep_solve(1))                              &
    !$omp depend(inout: p_lcol(sa))                                  

#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(trace_id, omp_get_thread_num())
#endif

    threadID = omp_get_thread_num()
#if defined(SPLLT_VERBOSE)
    print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
    print *, p_dep_solve
    print *, blk_dep_update
#endif

    ! Perform retangular update from diagonal block
    if(m .gt. n) then
       call slv_bwd_update(m - n, n, col, offset + n, p_index,      &
            p_lcol(sa + n * n : sa + n * m - 1), n, nrhs, p_rhs,    &
            p_upd(:, threadID + 1), ldr, p_xlocal(:, threadID + 1))
    endif

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          p_rhs(i)    = p_rhs(i) + p_upd(i, j)
        end do
      end do
    end do

    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa : sa + n * n - 1), &
         'Non-Transpose', 'Non-unit', nrhs, p_rhs, ldr)
    
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop (trace_id, omp_get_thread_num())
#endif
    !$omp end task

    call spllt_close_timer(task_manager%workerID, timer)
    call task_manager%ntask_submitted(1, nftask)
  end subroutine spllt_solve_bwd_block_task_with_ftask



  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         



  subroutine spllt_solve_bwd_update_task_with_ftask(blk, node, nrhs, upd, rhs, &
      ldr, xlocal, fkeep, trace_id, task_manager)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    use spllt_solve_dep_mod
    use task_manager_omp_mod
    use timer_mod
    implicit none

    integer, intent(in)                       :: blk  ! Index of block 
    integer, intent(in)                       :: node 
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:,:)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(task_manager_omp_t), intent(inout)  :: task_manager
    
    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID
    integer                     :: dep, ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    real(wp), pointer           :: p_lcol_solve(:)
    integer, pointer            :: p_dep_solve(:)
    integer, pointer            :: p_dep(:)
    integer                     :: blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer :: j
    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted
    integer :: nftask         ! #fake tasks inserted into the runtime

    type(spllt_timer_t), save   :: timer

    call spllt_open_timer(task_manager%workerID, &
      "spllt_solve_bwd_update_task_with_ftask", timer)

    ! Establish variables describing block
    n       = fkeep%bc(blk)%blkn
    m       = fkeep%bc(blk)%blkm
    blk_sa  = fkeep%bc(blk)%sa
    bcol    = fkeep%bc(blk)%bcol
    dcol    = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col     = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset  = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
    offset  = offset + (blk-fkeep%bc(blk)%dblk) * fkeep%nodes(node)%nb !this blk
    p_index => fkeep%nodes(node)%index
    p_lcol  => fkeep%lfact(bcol)%lcol
    p_bc    => fkeep%bc

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

    nftask      = 0

    blk_dep_update  = fkeep%bc(blk)%bwd_update_dep
    p_dep_solve     => fkeep%bc(blk)%bwd_solve_dep
    ndep            = size(p_dep_solve)

    do dep = 1, ndep - 1

      p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol
      p_blk_dep_solve => fkeep%bc(p_dep_solve(dep))

      !$omp task firstprivate(p_lcol_solve, p_blk_dep_solve)        &
      !$omp firstprivate(p_dep_solve)                               &
      !$omp firstprivate(blk, dep)                                  &
      !$omp private(threadID)                                       &
      !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
      !$omp depend(inout: p_dep_solve(1))

      threadID = omp_get_thread_num()
!     print '(a, i3, a, i3)', "UPD Fake Task dep ", blk, " [in : ", &
!       p_dep_solve(dep)

      !$omp end task

    end do
    nftask = merge(ndep - 1, 0, ndep .gt. 1)

    p_lcol_update => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
    p_lcol_solve  => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol

    p_blk_dep_update  => fkeep%bc(blk_dep_update)
    p_blk_dep_solve   => fkeep%bc(p_dep_solve(dep))

    !$omp task firstprivate(m, n, col, offset, node, bcol)        &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(p_blk_dep_update, p_blk_dep_solve)         &
    !$omp firstprivate(p_lcol_update, p_lcol_solve)               &
    !$omp firstprivate(blk)                                       &
    !$omp private(threadID)                                       &
    !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(in: p_dep_solve(1))                              &
    !$omp depend(inout: p_lcol(blk_sa))                            

    call trace_event_start(trace_id, omp_get_thread_num())

    threadID  = omp_get_thread_num()
!   print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!   print *, blk_dep_update
!   print *, p_dep_solve

    call  slv_bwd_update(m, n, col, offset, p_index,      &
          p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,   &
          p_rhs, p_upd(:, threadID + 1), ldr,             &
          xlocal(:, omp_get_thread_num() + 1))

    call trace_event_stop (trace_id, omp_get_thread_num())

    !$omp end task

    call spllt_close_timer(task_manager%workerID, timer)
    call task_manager%ntask_submitted(1, nftask)
  end subroutine spllt_solve_bwd_update_task_with_ftask
end module spllt_solve_task_mod
