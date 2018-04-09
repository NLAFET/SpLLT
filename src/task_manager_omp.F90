module task_manager_omp_mod
  use task_manager_mod
  implicit none

  type thread_task_info_t
    double precision  :: nflop                  ! # flop performed
    integer           :: ntask_run              ! #task run by this thread
    integer           :: ntask_insert           ! #task insert to the runtim
    integer           :: nfake_task_insert      ! #fake task insert
    integer           :: narray_allocated       ! #allocation
  end type thread_task_info_t

  type, extends(task_manager_base) :: task_manager_omp_t
    integer                           :: nthread_max
    type(thread_task_info_t), pointer :: thread_task_info(:)

    contains
      procedure :: init       => task_manager_omp_init
      procedure :: reset      => task_manager_omp_reset
      procedure :: print      => task_manager_omp_print
      procedure :: incr_alloc => task_manager_omp_incr_alloc
      procedure :: deallocate => task_manager_omp_deallocate

      procedure :: refresh_master
      procedure :: refresh_worker

      procedure :: incr_nrun
      procedure :: ntask_submitted
      procedure :: nflop_performed
      procedure :: nflop_reset
      procedure :: nflop_reset_all

      procedure :: solve_fwd_block_task
      procedure :: solve_fwd_update_task
      procedure :: solve_bwd_block_task
      procedure :: solve_bwd_update_task

  end type task_manager_omp_t

 contains



  subroutine task_manager_omp_incr_alloc(self, stat)
    implicit none
    class(task_manager_omp_t),  intent(inout) :: self
    integer,                    intent(in)    :: stat

    if(stat .ne. 0) then
      write(0,'(a)') "[>Allocation::Error] can not allocate the memory"
    else
      self%thread_task_info(self%workerID)%narray_allocated = 1 + &
        self%thread_task_info(self%workerID)%narray_allocated
    end if

  end subroutine task_manager_omp_incr_alloc



  subroutine nflop_reset_all(self)
    implicit none
    class(task_manager_omp_t),  intent(inout) :: self

    integer :: i

    do i = 0, self%nworker - 1
      self%thread_task_info(i)%nflop = 0.0
    end do

  end subroutine nflop_reset_all



  subroutine nflop_reset(self)
    implicit none
    class(task_manager_omp_t),  intent(inout) :: self

    self%thread_task_info(self%workerID)%nflop = 0.0

  end subroutine nflop_reset



  subroutine nflop_performed(self, nflop)
    implicit none
    class(task_manager_omp_t),  intent(inout) :: self
    double precision,           intent(in)    :: nflop
   !integer(long),              intent(in)    :: nflop

    self%thread_task_info(self%workerID)%nflop = self%thread_task_info&
      &(self%workerID)%nflop + nflop

  end subroutine nflop_performed



  subroutine incr_nrun(self)
    implicit none
    class(task_manager_omp_t), intent(inout)  :: self

    self%thread_task_info(self%workerID)%ntask_run = self%thread_task_info&
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
!   print *, "Current thread is reset to ", self%workerID

  end subroutine refresh_worker



  subroutine ntask_submitted(self, ntask, nftask)
    implicit none
    class(task_manager_omp_t), target,  intent(inout) :: self
    integer,                            intent(in)    :: ntask  ! #task insert
    integer,                            intent(in)    :: nftask ! #fake task

    type(thread_task_info_t), pointer :: p_th_task_info

    p_th_task_info => self%thread_task_info(self%workerID)

    p_th_task_info%nfake_task_insert = p_th_task_info%nfake_task_insert + nftask
    p_th_task_info%ntask_insert      = p_th_task_info%ntask_insert + ntask

  end subroutine ntask_submitted



  subroutine print_thread_task_info(msg, thread_id, thread_task_info)
    implicit none
    character(len=*),         intent(in)  :: msg
    integer,                  intent(in)  :: thread_id
    type(thread_task_info_t), intent(in)  :: thread_task_info

    print *, msg, " : ", thread_id
    print '(a, i9)', "#task insert        : ", thread_task_info%ntask_insert
    print '(a, i9)', "#fake task insert   : ", thread_task_info%&
      &nfake_task_insert
    print '(a, i9)', "#task run           : ", thread_task_info%ntask_run
    print '(a, i9)', "#array allocate     : ", thread_task_info%narray_allocated
    print '(a, es10.2)', "#flop               : ", thread_task_info%nflop

  end subroutine print_thread_task_info



  subroutine task_manager_omp_init_thread_task_info(thread_task_info)
    implicit none
    type(thread_task_info_t), intent(out) :: thread_task_info

    thread_task_info%nflop               = 0.0
    thread_task_info%ntask_run           = 0
    thread_task_info%ntask_insert        = 0
    thread_task_info%nfake_task_insert   = 0
    thread_task_info%narray_allocated    = 0
    
  end subroutine task_manager_omp_init_thread_task_info



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

   !print *, "Reset of the task_manager"
   !call self%print()
  end subroutine task_manager_omp_reset



  subroutine task_manager_omp_init(self, trace_names, stat)
    use trace_mod, only : trace_create_event
 !$ use omp_lib, only : omp_get_num_threads, omp_get_thread_num
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
    allocate(self%thread_task_info(0 : self%nthread_max - 1), stat=st1)
    if(st1 .eq. 0) then
      do i = lbound(self%thread_task_info, 1), ubound(self%thread_task_info, 1)
        call task_manager_omp_init_thread_task_info(self%thread_task_info(i))
      end do
    end if

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

   !print *, "Init of the task_manager"
   !call self%print()
  end subroutine task_manager_omp_init



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
    print '(a, i3)', "workerID      :", self%workerID
    print '(a, i3)', "masterWorker  :", self%masterWorker
    print '(a, i3)', "nworker       :", self%nworker
    print '(a, i3)', "nthread_max   :", self%nthread_max
    if(print_full .eq. 1) then
      if(associated(self%thread_task_info)) then
        do i = 0, self%nworker - 1
          call print_thread_task_info("OMP task information", i, &
            self%thread_task_info(i))
        end do
      end if
    end if

  end subroutine task_manager_omp_print



  subroutine task_manager_omp_deallocate(self)
    implicit none
    class(task_manager_omp_t), intent(inout) :: self

    if(associated(self%thread_task_info)) then
      deallocate(self%thread_task_info)
    end if

  end subroutine task_manager_omp_deallocate
  

  subroutine solve_fwd_block_task(task_manager, dblk, nrhs, upd, rhs, ldr, &
      xlocal, fkeep)
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

    call solve_fwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, task_manager%trace_ids(trace_fwd_block_pos))
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

    type(spllt_block), pointer  :: p_bc(:)
    type(spllt_timer), save     :: timer
        
   !nthread   = omp_get_num_threads()
    nthread = task_manager%nworker
    call spllt_open_timer(task_manager%nworker, task_manager%workerID, &
      "solve_fwd_block_task_worker", timer)

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

      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
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
      call spllt_tac(12, task_manager%workerID, timer)
    end if

    call spllt_close_timer(task_manager%workerID, timer)
    call task_manager%ntask_submitted(1, nftask)

  end subroutine solve_fwd_block_task_worker


  subroutine solve_fwd_update_task(task_manager, blk, node, nrhs, upd, rhs, &
      ldr, xlocal, fkeep)
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

    call solve_fwd_update_task_worker(task_manager, blk, node, nrhs, upd, rhs,&
      ldr, xlocal, fkeep, task_manager%trace_ids(trace_fwd_update_pos))
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

    type(spllt_timer), save     :: timer
        
    call spllt_open_timer(task_manager%nworker, task_manager%workerID, &
      "solve_fwd_update_task_worker", timer)

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

      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
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
      call spllt_tac(12, task_manager%workerID, timer)
    end if

    call spllt_close_timer(task_manager%workerID, timer)
    call task_manager%ntask_submitted(1, nftask)

  end subroutine solve_fwd_update_task_worker



  subroutine solve_bwd_block_task(task_manager, dblk, nrhs, upd, rhs, ldr, &
      xlocal, fkeep)
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
    
    call solve_bwd_block_task_worker(task_manager, dblk, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, task_manager%trace_ids(trace_bwd_block_pos))

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

    type(spllt_block), pointer  :: p_bc(:)
    type(spllt_timer), save     :: timer

   !nthread = omp_get_num_threads()
    nthread = task_manager%nworker
    call spllt_open_timer(task_manager%nworker, task_manager%workerID, &
      "solve_bwd_block_task_worker", timer)

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

      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
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
      call spllt_tac(12, task_manager%workerID, timer)
    end if
    call spllt_close_timer(task_manager%workerID, timer)
    call task_manager%ntask_submitted(1, nftask)
  end subroutine solve_bwd_block_task_worker

  

  subroutine solve_bwd_update_task(task_manager, blk, node, nrhs, upd, rhs, &
      ldr, xlocal, fkeep)
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

    call solve_bwd_update_task_worker(task_manager, blk, node, nrhs, upd, rhs,&
      ldr, xlocal, fkeep, task_manager%trace_ids(trace_bwd_update_pos))
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
    integer :: j
    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted
    integer :: nftask         ! #fake tasks inserted into the runtime

    type(spllt_timer), save     :: timer

    call spllt_open_timer(task_manager%nworker, task_manager%workerID, &
      "solve_bwd_update_task_worker", timer)

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

      call spllt_tic("Submit k-ary tree", 12, task_manager%workerID, timer)
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
      call spllt_tac(12, task_manager%workerID, timer)
    end if
    call spllt_close_timer(task_manager%workerID, timer)
    call task_manager%ntask_submitted(1, nftask)
  end subroutine solve_bwd_update_task_worker


end module task_manager_omp_mod