!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
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
  


  subroutine solve_fwd_block_task_worker(task_manager, dblk, nrhs, &
      n, rhs, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    implicit none
    
    class(task_manager_seq_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: dblk !Index of diagonal block
    integer,                    intent(in)    :: nrhs !Number of RHS
    integer,                    intent(in)    :: n
    real(wp), target,           intent(inout) :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, optional,          intent(in)    :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: blkm, blkn ! Block dimension
    integer                     :: bcol, node
    integer                     :: threadID
    integer                     :: i, ldy, nwdep
    double precision            :: flops
    integer                     :: traceID
    type(spllt_timer_t), save   :: timer

    real(wp),          pointer  :: p_y(:,:)
    real(wp),          pointer  :: p_rhs(:,:)
    integer,           pointer  :: p_wdep(:)
    integer,           pointer  :: p_index(:)
    integer,           pointer  :: p_order(:)
    real(wp),          pointer  :: p_lcol(:)

        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_block_task_worker", timer)
#endif

    threadID  = task_manager%workerID

    if(present(trace_id)) then 
      traceID = trace_id
    else
      traceID = task_manager%trace_ids(trace_fwd_block_pos)
    end if

    ! Get block info
    blkm      = fkeep%sbc(dblk)%blkm
    blkn      = fkeep%sbc(dblk)%blkn
    node      = fkeep%sbc(dblk)%node
    sa        = fkeep%sbc(dblk)%sa
    bcol      = fkeep%sbc(dblk)%bcol ! Current block column
    ldy       = fkeep%sbc(dblk)%ldu

    p_y         => fkeep%sbc(dblk)%p_upd(1 : ldy, 1 : nrhs)
    p_rhs       => rhs(1 : n, 1 : nrhs)
    p_index     => fkeep%sbc(dblk)%p_index
    p_order     => fkeep%porder
    p_lcol      => fkeep%lfact(bcol)%lcol(sa : sa + blkn * blkn - 1)
    p_wdep      => fkeep%sbc(dblk)%fwd_wdep

    nwdep       = size(p_wdep)
#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(traceID, threadID)
#endif

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd get RHS", 1, threadID, timer)
#endif

  !Gather part of RHS into y
  if(fkeep%nodes(node)%sblk_sa .eq. dblk) then
    do i = 1, blkn
     !print *, "p_y(", i, ",:) =", p_y(i,:)
     !print *, "p_rhs(", p_order(p_index(i)), ",:) =", p_rhs(p_order(p_index(i)), :)
      p_y(i, :) = p_rhs(p_order(p_index(i)), :)
    end do
  else
    do i = 1, blkn
      p_y(i, :) = p_rhs(p_order(p_index(i)), :) + p_y(i, :)
    end do
  end if
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(1, threadID, timer)
#endif
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd update UPD", 2, threadID, timer)
#endif
  !Compute the reduction if blk \in L_{:,1}
  if(fkeep%nodes(node)%sblk_sa .eq. dblk) then
    do i = 1, nwdep
      ! reduce children dep into p_y
      call fwd_update_upd(fkeep, dblk, p_wdep(i))
    end do
  end if
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(2, threadID, timer)
#endif

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd solve KERNEL", 3, threadID, timer)
#endif
  call slv_solve(blkm, blkn, p_lcol, 'Transpose    ', 'Non-unit', &
    nrhs, p_y, ldy)
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(3, threadID, timer)
#endif

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(flops)
#endif
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(traceID, threadID)
#endif
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_fwd_block_task_worker



  subroutine solve_fwd_update_task_worker(task_manager, blk, node, nrhs, &
      n, rhs, fkeep, trace_id)
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
    integer,                    intent(in)    :: n
    real(wp), target,           intent(in)    :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, optional,          intent(in)    :: trace_id

    ! Block info
    integer                     :: blkm, blkn         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol
    integer                     :: ldy, ldx
    integer                     :: nwdep
    integer                     :: i, dblk
    double precision            :: flops
    integer                     :: traceID
    logical                     :: reduction

    real(wp), pointer             :: p_lcol(:)
    real(wp)         , pointer    :: p_y(:,:)
    real(wp)         , pointer    :: p_xlocal(:,:)
    type(spllt_sblock_t), pointer :: p_bc(:)
    integer, pointer              :: p_wdep(:)

    type(spllt_timer_t), save   :: timer
        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_update_task_worker", timer)
#endif

    if(present(trace_id)) then
      traceID = trace_id
    else
      traceID = task_manager%trace_ids(trace_fwd_update_pos)
    end if

    ! Establish variables describing block
    blkn      = fkeep%sbc(blk)%blkn
    blkm      = fkeep%sbc(blk)%blkm
    dblk      = fkeep%sbc(blk)%dblk
    blk_sa    = fkeep%sbc(blk)%sa
    bcol      = fkeep%sbc(blk)%bcol
    ldy       = fkeep%sbc(dblk)%ldu
    ldx       = fkeep%sbc(blk)%ldu

    p_lcol    => fkeep%lfact(bcol)%lcol(blk_sa : blk_sa + blkn * blkm - 1)
    p_y       => fkeep%sbc(dblk)%p_upd
    p_xlocal  => fkeep%sbc(blk)%p_upd
    p_bc      => fkeep%sbc
    p_wdep    => fkeep%sbc(blk)%fwd_wdep

    nwdep     = size(p_wdep)
    reduction = (fkeep%nodes(node)%sblk_sa .eq. fkeep%sbc(blk)%dblk)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(traceID, task_manager%workerID)
#endif

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd update KERNEL", 2, task_manager%workerID, timer)
#endif
  call slv_fwd_update(blkm, blkn, p_lcol, blkn, n, nrhs, p_y, ldy, &
    p_xlocal, ldx, reduction)
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(2, task_manager%workerID, timer)
#endif

  !Compute the reduction if blk \in L_{:,1}
  if(reduction) then
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd update UPD", 1, task_manager%workerID, timer)
#endif
    do i = 1, nwdep
      ! reduce children dep into p_y
      call fwd_update_upd(fkeep, blk, p_wdep(i))
    end do
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(1, task_manager%workerID, timer)
#endif
  end if

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(flops)
#endif
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(traceID, task_manager%workerID)
#endif

#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_fwd_update_task_worker



  subroutine solve_bwd_block_task_worker(task_manager, dblk, nrhs, &
      n, rhs, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    implicit none
    
    class(task_manager_seq_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: dblk !Index of diagonal block
    integer,                    intent(in)    :: nrhs !Number of RHS
    integer,                    intent(in)    :: n
    real(wp), target,           intent(inout) :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, optional,          intent(in)    :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: blkm, blkn ! Block dimension
    integer                     :: bcol
    integer                     :: threadID
    integer                     :: i, ldy
    double precision            :: flops
    integer                     :: traceID
    type(spllt_timer_t), save   :: timer

    real(wp),          pointer  :: p_y(:,:)
    real(wp),          pointer  :: p_rhs(:,:)
    integer,           pointer  :: p_index(:)
    integer,           pointer  :: p_order(:)
    real(wp),          pointer  :: p_lcol(:)

        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_block_task_worker", timer)
#endif

    threadID  = task_manager%workerID

    if(present(trace_id)) then 
      traceID = trace_id
    else
      traceID = task_manager%trace_ids(trace_fwd_block_pos)
    end if

    ! Get block info
    blkm      = fkeep%sbc(dblk)%blkm
    blkn      = fkeep%sbc(dblk)%blkn
    sa        = fkeep%sbc(dblk)%sa
    bcol      = fkeep%sbc(dblk)%bcol ! Current block column
    ldy       = fkeep%sbc(dblk)%ldu

    p_y         => fkeep%sbc(dblk)%p_upd(1 : ldy, 1 : nrhs)
    p_rhs       => rhs(1 : n, 1 : nrhs)
    p_index     => fkeep%sbc(dblk)%p_index
    p_order     => fkeep%porder
    p_lcol      => fkeep%lfact(bcol)%lcol(sa : sa + blkn * blkn - 1)
#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(traceID, threadID)
#endif

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd solve KERNEL", 1, threadID, timer)
#endif
  call slv_solve(blkn, blkm, p_lcol, 'Non-Transpose', 'Non-unit', &
    nrhs, p_y, ldy)
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(1, threadID, timer)
#endif

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd set X", 2, threadID, timer)
#endif
  !Scatter the result into the solution vector
  do i = 1, blkn
    p_rhs(p_order(p_index(i)), :) = p_y(i, :)
  end do
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(2, threadID, timer)
#endif

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(flops)
#endif
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(traceID, threadID)
#endif
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(threadID, timer)
#endif

  end subroutine solve_bwd_block_task_worker



  subroutine solve_bwd_update_task_worker(task_manager, blk, node, nrhs, &
      n, rhs, fkeep, trace_id)
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
    integer,                    intent(in)    :: n
    real(wp), target,           intent(in)    :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, optional,          intent(in)    :: trace_id

    ! Block info
    integer                     :: blkm, blkn         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol
    integer                     :: ldy, ldx
    integer                     :: nwdep
    integer                     :: i, dblk
    double precision            :: flops
    integer                     :: traceID

    real(wp), pointer             :: p_lcol(:)
    real(wp)         , pointer    :: p_y(:,:)
    real(wp)         , pointer    :: p_xlocal(:,:)
    type(spllt_sblock_t), pointer :: p_bc(:)
    integer, pointer              :: p_wdep(:)

    type(spllt_timer_t), save   :: timer
        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_update_task_worker", timer)
#endif

    if(present(trace_id)) then
      traceID = trace_id
    else
      traceID = task_manager%trace_ids(trace_fwd_update_pos)
    end if

    ! Establish variables describing block
    blkn      = fkeep%sbc(blk)%blkn
    blkm      = fkeep%sbc(blk)%blkm
    dblk      = fkeep%sbc(blk)%dblk
    blk_sa    = fkeep%sbc(blk)%sa
    bcol      = fkeep%sbc(blk)%bcol
    ldy       = fkeep%sbc(dblk)%ldu
    ldx       = fkeep%sbc(blk)%ldu

    p_lcol    => fkeep%lfact(bcol)%lcol(blk_sa : blk_sa + blkn * blkm - 1)
    p_y       => fkeep%sbc(dblk)%p_upd
    p_xlocal  => fkeep%sbc(blk)%p_upd
    p_bc      => fkeep%sbc
    p_wdep    => fkeep%sbc(blk)%bwd_wdep

    nwdep        = size(p_wdep)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(traceID, task_manager%workerID)
#endif

  !Gather data into xlocal if belongs to last block column of L
  if(fkeep%nodes(node)%sblk_en .eq. fkeep%sbc(blk)%last_blk) then
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd update UPD", 1, task_manager%workerID, timer)
#endif
    do i = 1, nwdep
      ! reduce ancestor dep into p_y
      call bwd_update_upd(fkeep, blk, p_wdep(i))
    end do
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(1, task_manager%workerID, timer)
#endif
  end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd update KERNEL", 2, task_manager%workerID, timer)
#endif
  call slv_bwd_update(blkm, blkn, p_lcol, blkn, n, nrhs, &
    p_xlocal, ldx, p_y, ldy)
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(2, task_manager%workerID, timer)
#endif

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(flops)
#endif
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(traceID, task_manager%workerID)
#endif

#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_bwd_update_task_worker



  subroutine solve_fwd_subtree(task_manager, nrhs, rhs, n, fkeep, tree)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    implicit none

    class(task_manager_seq_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),                   intent(inout) :: rhs(n,nrhs)
    integer,                    intent(in)    :: n
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
  
    integer :: i
    type(spllt_timer_t), save :: timer

    print *, "Treatment of subtree ", tree%num
    call spllt_print_subtree(tree)

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_fwd_subtree", timer)
#endif

    do i = tree%node_sa, tree%node_en
      call solve_fwd_node(nrhs, rhs, n, fkeep, i, task_manager)
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_fwd_subtree



  subroutine solve_bwd_subtree(task_manager, nrhs, rhs, n, fkeep, tree)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    implicit none

    class(task_manager_seq_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    integer,                    intent(in)    :: n
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
  
    integer :: i
    type(spllt_timer_t), save :: timer

    print *, "Treatment of subtree ", tree%num
    call spllt_print_subtree(tree)

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_bwd_subtree", timer)
#endif

    do i = tree%node_en, tree%node_sa, -1
      call solve_bwd_node(nrhs, rhs, n, fkeep, i, task_manager)
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_bwd_subtree

end module task_manager_seq_mod
