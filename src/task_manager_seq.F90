module task_manager_seq_mod
  use task_manager_mod
  use worker_info_mod
  implicit none


  type, extends(task_manager_base) :: task_manager_seq_t
   !integer           :: ntask_run
   !double precision  :: nflop                  ! # flop performed
   !integer           :: narray_allocated       ! #allocation
    type(worker_info_t), pointer :: p_worker

    contains
      procedure :: init       => task_manager_seq_init
      procedure :: reset      => task_manager_seq_reset
      procedure :: print      => task_manager_seq_print
      procedure :: incr_alloc => task_manager_seq_incr_alloc
      procedure :: deallocate => task_manager_seq_deallocate

      procedure :: refresh_master
      procedure :: refresh_worker

      procedure :: incr_nrun
      procedure :: nflop_performed
      procedure :: get_nflop_performed
      procedure :: nflop_reset

      procedure :: solve_fwd_block_task   => solve_fwd_block_task_worker
      procedure :: solve_fwd_update_task  => solve_fwd_update_task_worker
      procedure :: solve_bwd_block_task   => solve_bwd_block_task_worker
      procedure :: solve_bwd_update_task  => solve_bwd_update_task_worker
      procedure :: solve_fwd_subtree_task => solve_fwd_subtree
      procedure :: solve_bwd_subtree_task => solve_bwd_subtree

  end type task_manager_seq_t

 contains



  subroutine task_manager_seq_incr_alloc(self, stat)
    implicit none
    class(task_manager_seq_t),  intent(inout) :: self
    integer,                    intent(in)    :: stat

    if(stat .ne. 0) then
      write(0,'(a)') "[>Allocation::Error] can not allocate the memory"
    else
      self%p_worker%narray_allocated = self%p_worker%narray_allocated + 1
    end if

  end subroutine task_manager_seq_incr_alloc



  subroutine nflop_reset(self)
    implicit none
    class(task_manager_seq_t),  intent(inout) :: self

    self%p_worker%nflop           = 0.0
    self%p_worker%nflop_performed = 0.0

  end subroutine nflop_reset



  subroutine nflop_performed(self, nflop)
    implicit none
    class(task_manager_seq_t),  intent(inout) :: self
    double precision,           intent(in)    :: nflop
   !integer(long),              intent(in)    :: nflop

    self%p_worker%nflop_performed = self%p_worker%nflop_performed + nflop

  end subroutine nflop_performed



  subroutine get_nflop_performed(self, nflop, thn)
    implicit none
    class(task_manager_seq_t),  intent(inout) :: self
    double precision,           intent(out)   :: nflop
    integer, optional,          intent(in)    :: thn

    nflop = self%p_worker%nflop_performed

  end subroutine get_nflop_performed



  subroutine incr_nrun(self)
    implicit none
    class(task_manager_seq_t), intent(inout)  :: self

    self%p_worker%ntask_run = self%p_worker%ntask_run + 1

  end subroutine incr_nrun



  subroutine refresh_master(self)
    implicit none
    class(task_manager_seq_t), intent(inout)  :: self

  end subroutine refresh_master



  subroutine refresh_worker(self)
 !$ use omp_lib, only : omp_get_thread_num
    implicit none
    class(task_manager_seq_t), intent(inout)  :: self

    self%workerID = 0
 !$ self%workerID = omp_get_thread_num()
    self%p_worker => self%worker_info(self%workerID)

  end subroutine refresh_worker



  subroutine task_manager_seq_reset(self)
 !$ use omp_lib, only : omp_get_thread_num, omp_get_num_threads
    implicit none
    class(task_manager_seq_t), intent(inout) :: self

    ! Parent class variables
    self%nworker             = 1
 !$ self%nworker             = omp_get_num_threads()
    self%workerID            = 0
 !$ self%workerID            = omp_get_thread_num()
    self%masterWorker        = self%workerID
    call self%p_worker%reset()

  end subroutine task_manager_seq_reset



  subroutine task_manager_seq_init(self, trace_names, stat)
    use trace_mod, only : trace_create_event
    implicit none
    class(task_manager_seq_t), target,  intent(inout) :: self
    character(len=*), optional,         intent(in)    :: trace_names(:)
    integer,          optional,         intent(out)   :: stat

    integer                           :: st0
    integer                           :: st1
    integer                           :: st2
    integer                           :: i
    integer                           :: ntotal_trace

    ! Initialize variables
    st0                 = 0
    st1                 = 0
    st2                 = 0
    self%trace_ids      => null()
    self%worker_info    => null()

    allocate(self%p_worker, stat = st0)
    call self%reset()

    ! Get ids of trace
    if(present(trace_names)) then
      ntotal_trace = ntrace_id + size(trace_names)
    else
      ntotal_trace = ntrace_id
    end if

    allocate(self%trace_ids(ntrace_id), stat=st2)

    if(st2 .eq. 0) then
#if defined(SPLLT_OMP_TRACE)
      do i = 1, ntrace_id
        call trace_create_event(task_manager_trace_names(i), self%trace_ids(i))
!       print *, "Create id ", self%trace_ids(i), " for step ", &
!         trace_names(i)
      end do
#else
      self%trace_ids(:) = 0
#endif
      if(present(trace_names)) then
        do i = ntrace_id + 1, ntotal_trace
          call trace_create_event(trace_names(i), self%trace_ids(i))
!         print *, "Create id ", self%trace_ids(i), " for step ", &
!           trace_names(i)
        end do
      end if
    end if

    ini_nde_id = self%trace_ids(trace_init_node_pos)
    fac_blk_id = self%trace_ids(trace_facto_blk_pos) 
    slv_blk_id = self%trace_ids(trace_solve_blk_pos) 
    upd_blk_id = self%trace_ids(trace_update_blk_pos) 
    upd_btw_id = self%trace_ids(trace_update_btw_pos) 

    if(present(stat)) then
      stat = st1 + st2
    end if

  end subroutine task_manager_seq_init



  subroutine task_manager_seq_print(self, msg, option)
    implicit none
    class(task_manager_seq_t), intent(in) :: self
    character(len=*), optional            :: msg
    integer, optional                     :: option
    
    if(present(msg)) then
      print *, msg
    end if

    print *, "Task manager env ::"
    print '(a, i3)', "workerID        : ", self%workerID
    print '(a, i3)', "masterWorker    : ", self%masterWorker
    print '(a, i3)', "nworker         : ", self%nworker
    call self%p_worker%print("Seq worker info", self%workerID)

  end subroutine task_manager_seq_print



  subroutine task_manager_seq_deallocate(self)
    implicit none
    class(task_manager_seq_t), intent(inout) :: self

    if(associated(self%trace_ids)) then
      deallocate(self%trace_ids)
    end if

  end subroutine task_manager_seq_deallocate
  


  subroutine solve_fwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, &
      ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    implicit none
    
    class(task_manager_seq_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: dblk !Index of diagonal block
    integer,                    intent(in)    :: nrhs !Number of RHS
    integer,                    intent(in)    :: ldr  !Leading dimension of RHS
    real(wp), target,           intent(inout) :: upd(:, :)
    real(wp), target,           intent(inout) :: rhs(ldr * nrhs)
    real(wp), target,           intent(inout) :: xlocal(:, :)
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer, optional,          intent(in)    :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: node
    integer                     :: nthread, threadID
    integer,           pointer  :: p_index(:)
    real(wp),          pointer  :: p_lcol(:)

    real(wp),          pointer  :: p_upd(:,:)
    real(wp),          pointer  :: p_xlocal(:,:)
    real(wp),          pointer  :: p_rhs(:)
    double precision            :: flops
    integer                     :: traceID

    type(spllt_timer_t), save   :: timer
        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_block_task_worker", timer)
#endif

    nthread   = task_manager%nworker
    nthread   = 1
    threadID  = task_manager%workerID

    if(present(trace_id)) then 
      traceID = trace_id
    else
      traceID = task_manager%trace_ids(trace_fwd_block_pos)
    end if

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
    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(traceID, threadID)
#endif

    call solve_fwd_block_work(m, n, col, offset, p_index, p_lcol, sa, nrhs, &
      p_upd, p_rhs, threadID, nthread, ldr, p_xlocal, flops)

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(flops)
#endif
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(traceID, threadID)
#endif
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%incr_nrun()

  end subroutine solve_fwd_block_task_worker



  subroutine solve_fwd_update_task_worker(task_manager, blk, node, nrhs, upd, &
      rhs, ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    implicit none
    
    class(task_manager_seq_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: blk  ! Index of block
    integer,                    intent(in)    :: node
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp), target,           intent(inout) :: upd(:,:)        
    real(wp), target,           intent(in)    :: rhs(ldr*nrhs)
    real(wp), target,           intent(out)   :: xlocal(:,:)
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer, optional,          intent(in)    :: trace_id

    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    double precision            :: flops
    integer                     :: traceID

    type(spllt_timer_t), save   :: timer
        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_update_task_worker", timer)
#endif

    if(present(trace_id)) then
      traceID = trace_id
    else
      traceID = task_manager%trace_ids(trace_fwd_update_pos)
    end if

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

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(traceID, task_manager%workerID)
#endif

    call solve_fwd_update_work(m, n, col, offset, p_index, p_lcol, blk_sa, &
      nrhs, p_upd, p_rhs, task_manager%workerID, ldr, p_xlocal, flops)

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(flops)
#endif
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(traceID, task_manager%workerID)
#endif

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%incr_nrun()

  end subroutine solve_fwd_update_task_worker



  subroutine solve_bwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, &
      ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    use utils_mod
    implicit none

    class(task_manager_seq_t), intent(inout)  :: task_manager
    integer, intent(in)                       :: dblk ! Index of diagonal block
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:, :)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, optional,          intent(in)    :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: bcol, col
    integer                     :: offset
    integer                     :: node
    integer                     :: threadID, nthread
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)

    real(wp), pointer :: p_upd(:,:)
    real(wp), pointer :: p_xlocal(:,:)
    real(wp), pointer :: p_rhs(:)

    double precision  :: flops
    integer           :: traceID

    type(spllt_block), pointer  :: p_bc(:)
    type(spllt_timer_t), save   :: timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_block_task_worker", timer)
#endif

    if(present(trace_id)) then
      traceID = trace_id
    else
      traceID = task_manager%trace_ids(trace_bwd_block_pos)
    end if

    threadID  = task_manager%workerID
    nthread   = task_manager%nworker
    nthread   = 1
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

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(traceID, threadID)
#endif

    call solve_bwd_block_work(m, n, col, offset, p_index, p_lcol, sa, nrhs, &
      p_upd, p_rhs, threadID, nthread, ldr, p_xlocal, flops)

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(flops)
#endif
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(traceID, threadID)
#endif
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer_flop(task_manager%workerID, timer, flops)
#endif
    call task_manager%incr_nrun()
  end subroutine solve_bwd_block_task_worker

  

  subroutine solve_bwd_update_task_worker(task_manager, blk, node, nrhs, upd, &
      rhs, ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    implicit none

    class(task_manager_seq_t), intent(inout)  :: task_manager
    integer, intent(in)                       :: blk  ! Index of block 
    integer, intent(in)                       :: node 
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:,:)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, optional,          intent(in)    :: trace_id
    
    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    integer, pointer            :: p_dep(:)
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    double precision            :: flops
    integer                     :: traceID

    type(spllt_timer_t), save   :: timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_update_task_worker", timer)
#endif

    if(present(trace_id)) then
      traceID = trace_id
    else
      traceID  = task_manager%trace_ids(trace_bwd_update_pos)
    end if

    threadID  = task_manager%workerID
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

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(traceID, threadID)
#endif

    call solve_bwd_update_work(m, n, col, offset, p_index, p_lcol, blk_sa, &
      nrhs, p_upd, p_rhs, threadID, ldr, p_xlocal, flops)

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(flops)
#endif
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(traceID, threadID)
#endif
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%incr_nrun()
  end subroutine solve_bwd_update_task_worker



  subroutine solve_fwd_subtree(task_manager, nrhs, rhs, ldr, fkeep, &
      tree, xlocal, rhs_local)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    implicit none

    class(task_manager_seq_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
    real(wp),                   intent(inout) :: xlocal(:,:)
    real(wp),                   intent(inout) :: rhs_local(:,:)
  
    integer :: i
    type(spllt_timer_t), save :: timer

    print *, "Treatment of subtree ", tree%num
    call spllt_print_subtree(tree)

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_fwd_subtree", timer)
#endif

    do i = tree%node_sa, tree%node_en
      call solve_fwd_node(nrhs, rhs, ldr, fkeep, i, xlocal, &
        rhs_local, task_manager)
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_fwd_subtree



  subroutine solve_bwd_subtree(task_manager, nrhs, rhs, ldr, fkeep, &
      tree, xlocal, rhs_local)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    implicit none

    class(task_manager_seq_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
    real(wp),                   intent(inout) :: xlocal(:,:)
    real(wp),                   intent(inout) :: rhs_local(:,:)
  
    integer :: i
    type(spllt_timer_t), save :: timer

    print *, "Treatment of subtree ", tree%num
    call spllt_print_subtree(tree)

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_bwd_subtree", timer)
#endif

    do i = tree%node_en, tree%node_sa, -1
      call solve_bwd_node(nrhs, rhs, ldr, fkeep, i, xlocal, &
        rhs_local, task_manager)
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_bwd_subtree

end module task_manager_seq_mod
