module task_manager_omp_mod
  use task_manager_mod
  implicit none

  type, extends(task_manager_base) :: task_manager_omp_t
    integer                       :: nthread_max
    integer                       :: nfork

    contains
      procedure :: init       => task_manager_omp_init
      procedure :: copy       => task_manager_omp_copy
      procedure :: reset      => task_manager_omp_reset
      procedure :: print      => task_manager_omp_print
      procedure :: incr_alloc => task_manager_omp_incr_alloc
      procedure :: deallocate => task_manager_omp_deallocate
      procedure :: semifork   => task_manager_omp_semifork

      procedure :: refresh_master
      procedure :: refresh_worker

      procedure :: incr_nrun
      procedure :: ntask_submitted
      procedure :: nflop_performed
      procedure :: get_nflop_performed
      procedure :: nflop_reset

      procedure :: solve_fwd_block_task
      procedure :: solve_fwd_update_task
      procedure :: solve_bwd_block_task
      procedure :: solve_bwd_update_task
      procedure :: solve_fwd_subtree_task
      procedure :: solve_bwd_subtree_task

  end type task_manager_omp_t

 contains



  subroutine task_manager_omp_incr_alloc(self, stat)
    implicit none
    class(task_manager_omp_t),  intent(inout) :: self
    integer,                    intent(in)    :: stat

    if(stat .ne. 0) then
      write(0,'(a)') "[>Allocation::Error] can not allocate the memory"
    else
      self%worker_info(self%workerID)%narray_allocated = 1 + &
        self%worker_info(self%workerID)%narray_allocated
    end if

  end subroutine task_manager_omp_incr_alloc



  subroutine nflop_reset(self, thn)
    implicit none
    class(task_manager_omp_t),  intent(inout) :: self
    integer, optional,          intent(in)    :: thn

    integer :: i

    if(present(thn)) then
      self%worker_info(thn)%nflop_performed = 0.0
    else
      do i = 0, self%nworker - 1
        self%worker_info(i)%nflop_performed = 0.0
      end do
    end if

  end subroutine nflop_reset



  subroutine nflop_performed(self, nflop)
    implicit none
    class(task_manager_omp_t),  intent(inout) :: self
    double precision,           intent(in)    :: nflop
   !integer(long),              intent(in)    :: nflop

    self%worker_info(self%workerID)%nflop_performed = self%worker_info&
      &(self%workerID)%nflop_performed + nflop

  end subroutine nflop_performed



  subroutine get_nflop_performed(self, nflop, thn)
    implicit none
    class(task_manager_omp_t),  intent(inout) :: self
    double precision,           intent(out)   :: nflop
    integer, optional,          intent(in)    :: thn

    integer :: i

    nflop = 0.0
    do i = 0, self%nworker - 1
      nflop = nflop + self%worker_info(i)%nflop_performed
    end do

  end subroutine get_nflop_performed



  subroutine incr_nrun(self)
    implicit none
    class(task_manager_omp_t), intent(inout)  :: self

    self%worker_info(self%workerID)%ntask_run = self%worker_info&
      &(self%workerID)%ntask_run + 1

  end subroutine incr_nrun



  subroutine refresh_master(self)
 !$ use omp_lib, ONLY : omp_get_thread_num
    implicit none
    class(task_manager_omp_t), intent(inout)  :: self

    self%masterWorker = omp_get_thread_num()
    self%workerID     = self%masterWorker
!   print *, "Master thread is reset to ", self%masterWorker

  end subroutine refresh_master



  subroutine refresh_worker(self)
 !$ use omp_lib, ONLY : omp_get_thread_num
    implicit none
    class(task_manager_omp_t), intent(inout)  :: self

    self%workerID = omp_get_thread_num()

  end subroutine refresh_worker



  subroutine ntask_submitted(self, ntask, nftask)
    implicit none
    class(task_manager_omp_t), target,  intent(inout) :: self
    integer,                            intent(in)    :: ntask  ! #task insert
    integer,                            intent(in)    :: nftask ! #fake task

    type(worker_info_t), pointer :: p_th_task_info

    p_th_task_info => self%worker_info(self%workerID)

    p_th_task_info%nfake_task_insert = p_th_task_info%nfake_task_insert + nftask
    p_th_task_info%ntask_insert      = p_th_task_info%ntask_insert + ntask

  end subroutine ntask_submitted



  subroutine task_manager_omp_reset(self)
 !$ use omp_lib, only : omp_get_num_threads, omp_get_thread_num
    implicit none
    class(task_manager_omp_t), intent(inout) :: self

    ! Parent class variables
    call self%refresh_master()
    self%nworker                 = 1
 !$ self%nworker                 = omp_get_num_threads()

    ! Additional variables
    self%nthread_max             = self%nworker
    self%nfork                   = 0

  end subroutine task_manager_omp_reset



  subroutine task_manager_omp_worker_reset(self)
 !$ use omp_lib, only : omp_get_num_threads, omp_get_thread_num
    implicit none
    class(task_manager_omp_t), intent(inout) :: self

    integer :: i
    do i = lbound(self%worker_info, 1), ubound(self%worker_info, 1)
      call self%worker_info(i)%reset()
    end do

  end subroutine task_manager_omp_worker_reset



  subroutine task_manager_omp_init(self, trace_names, stat)
    use trace_mod, only : trace_create_event
    implicit none
    class(task_manager_omp_t), target,  intent(inout) :: self
    character(len=*), optional,         intent(in)    :: trace_names(:)
    integer,          optional,         intent(out)   :: stat

    integer                           :: st1
    integer                           :: st2
    integer                           :: i
    integer                           :: ntotal_trace

    ! Initialize variables
    st1 = 0
    st2 = 0
    self%trace_ids => null()
    call self%reset()

    ! Build task info space
    allocate(self%worker_info(0 : self%nthread_max - 1), stat=st1)
    if(st1 .eq. 0) then
      do i = lbound(self%worker_info, 1), ubound(self%worker_info, 1)
        call self%worker_info(i)%init()
      end do
    end if

    ! Get ids of trace
    if(present(trace_names)) then
      ntotal_trace = ntrace_id + size(trace_names)
    else
      ntotal_trace = ntrace_id
    end if

    allocate(self%trace_ids(ntotal_trace), stat=st2)

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

!   print *, "Init of the task_manager"
!   call self%print()
  end subroutine task_manager_omp_init



  subroutine task_manager_omp_semifork(self, task_manager)
    implicit none
    class(task_manager_omp_t), target,  intent(inout) :: self
    class(task_manager_base),           intent(inout) :: task_manager

    ! Initialize variables
    call self%copy(task_manager)
    self%nfork = self%nfork + 1

  end subroutine task_manager_omp_semifork



  subroutine task_manager_omp_copy(self, task_manager)
    implicit none
    class(task_manager_omp_t), target,  intent(in)    :: self
    class(task_manager_base),           intent(inout) :: task_manager

    task_manager%trace_ids    => self%trace_ids
    task_manager%worker_info  => self%worker_info
    task_manager%nworker      =  self%nworker
    task_manager%masterWorker =  self%masterWorker
    task_manager%workerID     =  self%workerID

  end subroutine task_manager_omp_copy



  subroutine task_manager_omp_print(self, msg, option)
    implicit none
    class(task_manager_omp_t), intent(in) :: self
    character(len=*), optional            :: msg
    integer, optional                     :: option
    
    integer :: i
    integer :: print_full

    print_full = 1
    if(present(option)) then
      print_full = option
    end if

    if(present(msg)) then
      print *, msg
    end if

    print *, "Task manager env ::"
    print '(a, i6)', "|workerID      :", self%workerID
    print '(a, i6)', "|masterWorker  :", self%masterWorker
    print '(a, i6)', "|nworker       :", self%nworker
    print '(a, i6)', "|nthread_max   :", self%nthread_max
    print '(a, i6)', "|nfork         :", self%nfork
    if(print_full .eq. 1) then
      if(associated(self%worker_info)) then
        do i = 0, self%nworker - 1
          call self%worker_info(i)%print("OMP task information", i)
        end do
      end if
    end if

  end subroutine task_manager_omp_print



  subroutine task_manager_omp_deallocate(self)
    implicit none
    class(task_manager_omp_t), intent(inout) :: self

    if(associated(self%worker_info)) then
      deallocate(self%worker_info)
    end if

  end subroutine task_manager_omp_deallocate
  

  subroutine solve_fwd_block_task(task_manager, dblk, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id)
    use spllt_data_mod
    implicit none
    
    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: dblk !Index of diagonal block
    integer,                    intent(in)    :: nrhs !Number of RHS
    integer,                    intent(in)    :: ldr  !Leading dimension of RHS
    real(wp), target,           intent(inout) :: upd(:, :)
    real(wp), target,           intent(inout) :: rhs(ldr * nrhs)
    real(wp), target,           intent(inout) :: xlocal(:, :)
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer, optional,          intent(in)    :: trace_id

    if(present(trace_id)) then
      call solve_fwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, ldr,&
        xlocal, fkeep, trace_id)
    else
      call solve_fwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, ldr,&
        xlocal, fkeep, task_manager%trace_ids(trace_fwd_block_pos))
    end if
  end subroutine solve_fwd_block_task

  subroutine solve_fwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, &
      ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use omp_lib, ONLY : omp_get_thread_num
    use timer_mod
    implicit none
    
    type (task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: dblk !Index of diagonal block
    integer,                    intent(in)    :: nrhs !Number of RHS
    integer,                    intent(in)    :: ldr  !Leading dimension of RHS
    real(wp), target,           intent(inout) :: upd(:, :)
    real(wp), target,           intent(inout) :: rhs(ldr * nrhs)
    real(wp), target,           intent(inout) :: xlocal(:, :)
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: node
   !integer                     :: dep
    integer                     :: i, j, r, chunkth
    integer                     :: ndep
    integer                     :: nthread, threadID
    integer,           pointer  :: p_index(:)
    real(wp),          pointer  :: p_lcol(:)

    real(wp),          pointer  :: p_upd(:,:)
    real(wp),          pointer  :: p_xlocal(:,:)
    real(wp),          pointer  :: p_rhs(:)
    integer,           pointer  :: p_dep(:)
    integer                     :: chunk, chunk_size
    integer                     :: ndep_lvl, lvl, nchunk
    integer                     :: beta, alpha
    logical                     :: all_task_submitted
    integer                     :: nftask
    double precision            :: flops

    type(spllt_block), pointer  :: p_bc(:)

    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer
        
   !nthread   = omp_get_num_threads()
    nthread = task_manager%nworker
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_block_task_worker", timer)
#endif

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
    p_dep     => fkeep%bc(dblk)%fwd_dep
    ndep      = size(p_dep)
    nftask    = 0
    chunk     = 0
    ndep_lvl  = 0

    if(ndep .eq. 0) then
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("CASE(0)", 1, task_manager%workerID, timer)
#endif
      !$omp task                                &
      include 'include/spllt_solve_fwd_block_omp_decl.F90.inc'

#include "include/spllt_solve_fwd_block_worker.F90.inc"

      !$omp end task
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(1, task_manager%workerID, timer)
#endif
    else
      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      lvl   = 1
      alpha = 1
      ndep_lvl = ndep ! #dep local to the lvl
      all_task_submitted = .false.

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
#endif
      do while(.not. all_task_submitted)
      
        nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

        beta = 1 - alpha

        do chunkth = 1, nchunk
          chunk_size = merge(ndep_lvl - (chunkth - 1) * chunk, chunk, &
            chunkth * chunk .gt. ndep_lvl)

          select case(chunk_size)

            case(0)
              print *, "No dep ??"

            !
            !This file contains the remaining cases that are generated through a script
            !
#include "include/spllt_fwd_block_cases.F90.inc"

          end select

          beta = beta + chunk_size
          if(ndep_lvl .le. chunk) then
            all_task_submitted = .true.
          else
            nftask = nftask + 1
          end if

        end do
        if(ndep_lvl .gt. chunk) then
          ndep_lvl = nchunk
          lvl = lvl + 1
          alpha = alpha * chunk
        end if
      end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(12, task_manager%workerID, timer)
#endif
    end if

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%ntask_submitted(1, nftask)

  end subroutine solve_fwd_block_task_worker


  subroutine solve_fwd_update_task(task_manager, blk, node, nrhs, upd, rhs, &
      ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    implicit none

    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: blk  ! Index of block
    integer,                    intent(in)    :: node
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp), target,           intent(inout) :: upd(:,:)        
    real(wp), target,           intent(in)    :: rhs(ldr*nrhs)
    real(wp), target,           intent(out)   :: xlocal(:,:)
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer, optional,          intent(in)    :: trace_id

    if(present(trace_id)) then
      call solve_fwd_update_task_worker(task_manager, blk, node, nrhs, upd, &
        rhs, ldr, xlocal, fkeep, trace_id)
    else
      call solve_fwd_update_task_worker(task_manager, blk, node, nrhs, upd, &
        rhs, ldr, xlocal, fkeep, task_manager%trace_ids(trace_fwd_update_pos))
    end if
  end subroutine solve_fwd_update_task

  subroutine solve_fwd_update_task_worker(task_manager, blk, node, nrhs, upd, &
      rhs, ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    use omp_lib, ONLY : omp_get_thread_num
    implicit none
    
    type (task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: blk  ! Index of block
    integer,                    intent(in)    :: node
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp), target,           intent(inout) :: upd(:,:)        
    real(wp), target,           intent(in)    :: rhs(ldr*nrhs)
    real(wp), target,           intent(out)   :: xlocal(:,:)
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer, intent(in)                       :: trace_id

    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID
    integer                     :: ndep
   !integer                     :: dep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    integer, pointer            :: p_dep(:)
   !integer                     :: blk_dep_solve
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer                     :: j
    integer                     :: chunk, chunk_size
    integer                     :: ndep_lvl, lvl, nchunk
    integer                     :: beta, alpha
    logical                     :: all_task_submitted
    integer                     :: nftask
    double precision            :: flops

   !type(spllt_timer_t), save   :: timer
    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer
        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_update_task_worker", timer)
#endif

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

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs
    p_dep     => fkeep%bc(blk)%fwd_dep

    nftask    = 0
    ndep      = size(p_dep)
    chunk     = 0
    ndep_lvl  = 0

    if(ndep .eq. 0) then
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("CASE(0)", 1, task_manager%workerID, timer)
#endif
      !$omp task                                &
      include 'include/spllt_solve_fwd_update_omp_decl.F90.inc'

#include "include/spllt_solve_fwd_update_worker.F90.inc"

      !$omp end task
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(1, task_manager%workerID, timer)
#endif
    else

      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      lvl   = 1
      alpha = 1
      ndep_lvl = ndep ! #dep local to the lvl
      all_task_submitted = .false.

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
#endif
      do while(.not. all_task_submitted)
      
        nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

        beta = 1 - alpha

        do j = 1, nchunk
          chunk_size = merge(ndep_lvl - (j - 1) * chunk, chunk, &
            j * chunk .gt. ndep_lvl)
          select case(chunk_size)

            case(0)
              print *, "No dep ?? "

            !
            !This file contains the remaining cases that are generated through a script
            !
#include "include/spllt_fwd_update_cases.F90.inc"

          end select
          beta = beta + chunk_size
          if(ndep_lvl .le. chunk) then
            all_task_submitted = .true.
          else
            nftask = nftask + 1
          end if

        end do
        if(ndep_lvl .gt. chunk) then
          ndep_lvl = nchunk
          lvl = lvl + 1
          alpha = alpha * chunk
        end if
      end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(12, task_manager%workerID, timer)
#endif
    end if

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%ntask_submitted(1, nftask)

  end subroutine solve_fwd_update_task_worker



  subroutine solve_bwd_block_task(task_manager, dblk, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id)
    use spllt_data_mod
    implicit none

    class(task_manager_omp_t), intent(inout)  :: task_manager
    integer, intent(in)                       :: dblk ! Index of diagonal block
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:, :)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, optional,          intent(in)    :: trace_id
    
    if(present(trace_id)) then
      call solve_bwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, &
        ldr, xlocal, fkeep, trace_id)
    else
      call solve_bwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, &
        ldr, xlocal, fkeep, task_manager%trace_ids(trace_bwd_block_pos))
    end if

  end subroutine solve_bwd_block_task

  subroutine solve_bwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, &
      ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    use utils_mod
    use omp_lib, ONLY : omp_get_thread_num
    implicit none

    type (task_manager_omp_t), intent(inout)  :: task_manager
    integer, intent(in)                       :: dblk ! Index of diagonal block
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:, :)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: bcol, col
    integer                     :: offset
    integer                     :: node
    integer                     :: i, j, r, chunkth
    integer                     :: threadID, nthread
    integer                     :: ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)

    real(wp), pointer :: p_upd(:,:)
    real(wp), pointer :: p_xlocal(:,:)
    real(wp), pointer :: p_rhs(:)
    integer,  pointer :: p_dep(:)

    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted
    integer :: nftask
    double precision :: flops

    type(spllt_block), pointer  :: p_bc(:)
   !type(spllt_timer_t), save   :: timer
    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer

   !nthread = omp_get_num_threads()
    nthread = task_manager%nworker
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_block_task_worker", timer)
#endif

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

    nftask    = 0
    p_dep     => fkeep%bc(dblk)%bwd_dep
    ndep      = size(p_dep)
    chunk     = 0
    ndep_lvl  = 0

    if(ndep .eq. 0) then
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("CASE(0)", 1, task_manager%workerID, timer)
#endif
      !$omp task                                &
      include 'include/spllt_solve_bwd_block_omp_decl.F90.inc'

#include "include/spllt_solve_bwd_block_worker.F90.inc"

      !$omp end task
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(1, task_manager%workerID, timer)
#endif
    else
      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      lvl   = 1
      alpha = 1
      ndep_lvl = ndep ! #dep local to the lvl
      all_task_submitted = .false.

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
#endif
      do while(.not. all_task_submitted)
      
        nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

        beta = 1 - alpha

        do chunkth = 1, nchunk
          chunk_size = merge(ndep_lvl - (chunkth - 1) * chunk, chunk, &
            chunkth * chunk .gt. ndep_lvl)

          select case(chunk_size)

            case(0)
              print *, "No dep ??"

            !
            !This file contains the remaining cases that are generated through a script
            !
#include "include/spllt_bwd_block_cases.F90.inc"

          end select

          beta = beta + chunk_size
          if(ndep_lvl .le. chunk) then
            all_task_submitted = .true.
          else
            nftask = nftask + 1
          end if

        end do
        if(ndep_lvl .gt. chunk) then
          ndep_lvl = nchunk
          lvl = lvl + 1
          alpha = alpha * chunk
        end if
      end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(12, task_manager%workerID, timer)
#endif
    end if
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%ntask_submitted(1, nftask)
  end subroutine solve_bwd_block_task_worker

  

  subroutine solve_bwd_update_task(task_manager, blk, node, nrhs, upd, rhs, &
      ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    implicit none

    class(task_manager_omp_t), intent(inout)  :: task_manager
    integer, intent(in)                       :: blk  ! Index of block 
    integer, intent(in)                       :: node 
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:,:)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, optional,          intent(in)    :: trace_id

    if(present(trace_id)) then
      call solve_bwd_update_task_worker(task_manager, blk, node, nrhs, upd, &
        rhs, ldr, xlocal, fkeep, trace_id)
    else
      call solve_bwd_update_task_worker(task_manager, blk, node, nrhs, upd, &
        rhs, ldr, xlocal, fkeep, task_manager%trace_ids(trace_bwd_update_pos))
    end if
  end subroutine solve_bwd_update_task
    
  subroutine solve_bwd_update_task_worker(task_manager, blk, node, nrhs, upd, &
      rhs, ldr, xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    use omp_lib, ONLY : omp_get_thread_num
    implicit none

    type (task_manager_omp_t), intent(inout)  :: task_manager
    integer, intent(in)                       :: blk  ! Index of block 
    integer, intent(in)                       :: node 
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:,:)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    
    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID
    integer                     :: ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    integer, pointer            :: p_dep(:)
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer                     :: j
    integer                     :: chunk, chunk_size
    integer                     :: ndep_lvl, lvl, nchunk
    integer                     :: beta, alpha
    logical                     :: all_task_submitted
    integer                     :: nftask   ! #fake tasks inserted 
                                            ! into the runtime
    double precision            :: flops

   !type(spllt_timer_t), save   :: timer
    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_update_task_worker", timer)
#endif

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

    nftask    = 0
    p_dep     => fkeep%bc(blk)%bwd_dep
    ndep      = size(p_dep)
    chunk     = 0
    ndep_lvl  = 0

    if(ndep .eq. 0) then
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("CASE(0)", 1, task_manager%workerID, timer)
#endif
      !$omp task                                &
      include 'include/spllt_solve_bwd_update_omp_decl.F90.inc'

#include "include/spllt_solve_bwd_update_worker.F90.inc"

      !$omp end task
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(1, task_manager%workerID, timer)
#endif
    else

      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      lvl   = 1
      alpha = 1
      ndep_lvl = ndep ! #dep local to the lvl
      all_task_submitted = .false.

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
#endif
      do while(.not. all_task_submitted)
      
        nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

        beta = 1 - alpha

        do j = 1, nchunk
          chunk_size = merge(ndep_lvl - (j - 1) * chunk, chunk, &
            j * chunk .gt. ndep_lvl)

          select case(chunk_size)

            case(0)
              print *, "No dep ?? "

            !
            !This file contains the remaining cases that are generated through a script
            !
#include "include/spllt_bwd_update_cases.F90.inc"

          end select

          beta = beta + chunk_size
          if(ndep_lvl .le. chunk) then
            all_task_submitted = .true.
          else
            nftask = nftask + 1
          end if

        end do
        if(ndep_lvl .gt. chunk) then
          ndep_lvl = nchunk
          lvl = lvl + 1
          alpha = alpha * chunk
        end if
      end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(12, task_manager%workerID, timer)
#endif
    end if
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%ntask_submitted(1, nftask)
  end subroutine solve_bwd_update_task_worker



  subroutine solve_fwd_subtree_task(task_manager, nrhs, rhs, ldr, fkeep, &
      tree, xlocal, rhs_local)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_seq_mod
    implicit none

    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
    real(wp),                   intent(inout) :: xlocal(:,:)
    real(wp),                   intent(inout) :: rhs_local(:,:)

    call solve_fwd_subtree_task_worker(task_manager, nrhs, rhs, ldr, &
      fkeep, tree, xlocal, rhs_local)
  end subroutine solve_fwd_subtree_task



  subroutine solve_fwd_subtree_task_worker(task_manager, nrhs, rhs, ldr, &
      fkeep, tree, xlocal, rhs_local)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_seq_mod
    implicit none

    type (task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),          target,  intent(inout) :: rhs(ldr*nrhs)
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
    real(wp),          target,  intent(inout) :: xlocal(:,:)
    real(wp),          target,  intent(inout) :: rhs_local(:,:)
  
    integer                     :: trace_id
    real(wp), pointer           :: p_rhs_local(:,:)
    real(wp), pointer           :: p_xlocal(:,:)
    real(wp), pointer           :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer                     :: i
    integer                     :: sa, en, blk_en
    type(task_manager_seq_t)    :: sub_task_manager
   !type(spllt_timer_t), save   :: timer
    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_subtree", timer)
#endif

    p_rhs_local => rhs_local
    p_xlocal    => xlocal
    p_rhs       => rhs
    p_bc        => fkeep%bc
    sa          = tree%node_sa
    en          = tree%node_en
    blk_en      = fkeep%nodes(en)%blk_en
    call task_manager%semifork(sub_task_manager)

    !$omp task                                        &
!   !$omp default(none)                               &
    !$omp depend(inout: p_bc(blk_en))                 &
    !$omp firstprivate(p_rhs, p_xlocal, p_rhs_local)  &
    !$omp firstprivate(nrhs, ldr)                     &
    !$omp firstprivate(sub_task_manager)              &
    !$omp firstprivate(p_bc, blk_en)                  &
    !$omp firstprivate(en, sa)                        &
    !$omp shared(fkeep)                               &
    !$omp private(i, trace_id)

    call sub_task_manager%refresh_worker()

#if defined(SPLLT_OMP_TRACE)
    trace_id = task_manager%trace_ids(trace_fwd_subtree_pos)
    call trace_event_start(trace_id, sub_task_manager%workerID)
#endif

#if defined(SPLLT_VERBOSE)
  print '(a, i3, a)', "SLV_FWD  Task dep of ", blk_en, " [in : "
#endif

    do i = sa, en
      call solve_fwd_node(nrhs, p_rhs, ldr, fkeep, i, p_xlocal, &
        p_rhs_local, sub_task_manager, no_trace)
    end do

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(trace_id, sub_task_manager%workerID)
#endif

    !$omp end task

    call task_manager%ntask_submitted(1, 0)
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_fwd_subtree_task_worker



  subroutine solve_bwd_subtree_task(task_manager, nrhs, rhs, ldr, fkeep, &
      tree, xlocal, rhs_local)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_seq_mod
    implicit none

    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
    real(wp),                   intent(inout) :: xlocal(:,:)
    real(wp),                   intent(inout) :: rhs_local(:,:)

    call solve_bwd_subtree_task_worker(task_manager, nrhs, rhs, ldr, &
      fkeep, tree, xlocal, rhs_local)
  end subroutine solve_bwd_subtree_task



  subroutine solve_bwd_subtree_task_worker(task_manager, nrhs, rhs, ldr, &
      fkeep, tree, xlocal, rhs_local)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_seq_mod
    implicit none

    type (task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),          target,  intent(inout) :: rhs(ldr*nrhs)
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    type(spllt_fkeep), target,  intent(in)    :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
    real(wp),          target,  intent(inout) :: xlocal(:,:)
    real(wp),          target,  intent(inout) :: rhs_local(:,:)
  
    integer,  pointer           :: p_dep(:)
    real(wp), pointer           :: p_rhs_local(:,:)
    real(wp), pointer           :: p_xlocal(:,:)
    real(wp), pointer           :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer                     :: ndep
    integer                     :: i
    integer                     :: sa, en, blk_en
    integer                     :: j
    integer                     :: chunk, chunk_size
    integer                     :: ndep_lvl, lvl, nchunk
    integer                     :: beta, alpha
    integer                     :: nftask   ! #fake tasks inserted 
                                            ! into the runtime
    logical                     :: all_task_submitted
    type(task_manager_seq_t)    :: sub_task_manager
   !type(spllt_timer_t), save   :: timer
    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    integer                     :: trace_id
    p_timer => timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_bwd_subtree", timer)
#endif

    p_rhs_local => rhs_local
    p_xlocal    => xlocal
    p_rhs       => rhs
    p_bc        => fkeep%bc
    sa          = tree%node_sa
    en          = tree%node_en
    blk_en      = fkeep%nodes(en)%blk_en
    call task_manager%semifork(sub_task_manager)

    nftask      = 0
    p_dep       => fkeep%bc(blk_en)%bwd_dep
    ndep        = size(p_dep)

    if(ndep .eq. 0) then
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("CASE(0)", 1, task_manager%workerID, timer)
#endif
      !$omp task                                &
      include 'include/spllt_solve_bwd_node_omp_decl.F90.inc'

#include "include/spllt_solve_bwd_node_worker.F90.inc"

      !$omp end task
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(1, task_manager%workerID, timer)
#endif
    else

      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      lvl   = 1
      alpha = 1
      ndep_lvl = ndep ! #dep local to the lvl
      all_task_submitted = .false.

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
#endif
      do while(.not. all_task_submitted)
      
        nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

        beta = 1 - alpha

        do j = 1, nchunk
          chunk_size = merge(ndep_lvl - (j - 1) * chunk, chunk, &
            j * chunk .gt. ndep_lvl)

          select case(chunk_size)

            case(0)
              print *, "No dep ?? "

            !
            !This file contains the remaining cases that are generated through a script
            !
#include "include/spllt_bwd_node_cases.F90.inc"

          end select

          beta = beta + chunk_size
          if(ndep_lvl .le. chunk) then
            all_task_submitted = .true.
          else
            nftask = nftask + 1
          end if

        end do
        if(ndep_lvl .gt. chunk) then
          ndep_lvl = nchunk
          lvl = lvl + 1
          alpha = alpha * chunk
        end if
      end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(12, task_manager%workerID, timer)
#endif
    end if

    call task_manager%ntask_submitted(1, nftask)
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif

  end subroutine solve_bwd_subtree_task_worker

end module task_manager_omp_mod
