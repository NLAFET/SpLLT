!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
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
  


  subroutine solve_fwd_block_task(task_manager, dblk, nrhs, &
      n, rhs, fkeep, trace_id)
    use spllt_data_mod
    implicit none
    
    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: dblk !Index of diagonal block
    integer,                    intent(in)    :: nrhs !Number of RHS
    integer,                    intent(in)    :: n
    real(wp), target,           intent(inout) :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, optional,          intent(in)    :: trace_id

    integer :: id

    id = task_manager%trace_ids(trace_fwd_block_pos)
    if(present(trace_id)) then
      id = trace_id
    end if

    call solve_fwd_block_task_worker(task_manager, dblk, nrhs, &
      n, rhs, fkeep, id)


  end subroutine solve_fwd_block_task

  subroutine solve_fwd_block_task_worker(task_manager, dblk, nrhs, n, rhs,&
      fkeep, trace_id)
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
    integer,                    intent(in)    :: n
    real(wp), target,           intent(inout) :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer,                    intent(in)    :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: blkm, blkn ! Block dimension
    integer                     :: bcol
    integer                     :: node
    integer                     :: i, j, r, chunkth
    integer                     :: ndep, nwdep
    integer                     :: threadID
    integer                     :: ldy

    real(wp),          pointer    :: p_y(:,:)
    real(wp),          pointer    :: p_rhs(:,:)
    integer,           pointer    :: p_index(:)
    integer,           pointer    :: p_order(:)
    real(wp),          pointer    :: p_lcol(:)
    type(spllt_sblock_t), pointer :: p_bc(:)
    integer,           pointer    :: p_dep(:)
    integer,           pointer    :: p_wdep(:)

    integer                     :: chunk, chunk_size
    integer                     :: ndep_lvl, lvl, nchunk
    integer                     :: beta, alpha
    logical                     :: all_task_submitted
    integer                     :: nftask
    double precision            :: flops

    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer
        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_block_task_worker", timer)
#endif

    ! Get block info
    node      = fkeep%sbc(dblk)%node
    blkm      = fkeep%sbc(dblk)%blkm
    blkn      = fkeep%sbc(dblk)%blkn
    sa        = fkeep%sbc(dblk)%sa
    bcol      = fkeep%sbc(dblk)%bcol ! Current block column
    ldy       = fkeep%sbc(dblk)%ldu

    p_y         => fkeep%sbc(dblk)%p_upd(1 : ldy, 1 : nrhs)
    p_rhs       => rhs(1 : n, 1 : nrhs)
    p_index     => fkeep%sbc(dblk)%p_index
    p_order     => fkeep%porder
   !print *, "from ", sa, 'to', sa + blkn * blkn - 1, "FOR bcol", bcol
   !print *, "lbound lfact", lbound(fkeep%lfact), "ubound", ubound(fkeep%lfact)
   !print *, "lbound lcol ", lbound(fkeep%lfact(bcol)%lcol), "ubound", ubound(fkeep%lfact(bcol)%lcol)
    p_lcol      => fkeep%lfact(bcol)%lcol(sa : sa + blkn * blkn - 1)
    p_bc        => fkeep%sbc
    p_dep       => fkeep%sbc(dblk)%fwd_dep
    p_wdep      => fkeep%sbc(dblk)%fwd_wdep

    ndep      = size(p_dep)
    nwdep     = size(p_wdep)
    nftask    = 0
    chunk     = 0
    ndep_lvl  = 0

    if(ndep .eq. 0) then
#include "spllt_fwd_block_nodep.F90"
    else
      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      chunk = fkeep%chunk
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
#include "include/spllt_fwd_block_cases.F90"

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

#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%ntask_submitted(1, nftask)

  end subroutine solve_fwd_block_task_worker



  subroutine solve_fwd_update_task(task_manager, blk, node, nrhs, &
      n, rhs, fkeep, trace_id)
    use spllt_data_mod
    implicit none

    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: blk  ! Index of block
    integer,                    intent(in)    :: node
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n
    real(wp), target,           intent(in)    :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, optional,          intent(in)    :: trace_id

    integer :: id

    id = task_manager%trace_ids(trace_fwd_update_pos)
    if(present(trace_id)) then
      id = trace_id
    end if

    call solve_fwd_update_task_worker(task_manager, blk, node, nrhs, & 
      n, rhs, fkeep, id)

  end subroutine solve_fwd_update_task

  subroutine solve_fwd_update_task_worker(task_manager, blk, node, nrhs, &
      n, rhs, fkeep, trace_id)
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
    integer,                    intent(in)    :: n
    real(wp), target,           intent(in)    :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, intent(in)                       :: trace_id

    ! Block info
    integer                     :: blkm, blkn         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol
    integer                     :: dblk
    integer                     :: ldy, ldx
    integer                     :: threadID
    integer                     :: ndep, nwdep

    real(wp), pointer             :: p_lcol(:)
    real(wp)         , pointer    :: p_y(:,:)
    real(wp)         , pointer    :: p_xlocal(:,:)
    type(spllt_sblock_t), pointer :: p_bc(:)
    integer, pointer              :: p_dep(:)
    integer, pointer              :: p_wdep(:)

    integer                     :: i
    integer                     :: j
    integer                     :: chunk, chunk_size
    integer                     :: ndep_lvl, lvl, nchunk
    integer                     :: beta, alpha
    logical                     :: all_task_submitted
    integer                     :: nftask
    double precision            :: flops

    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer
        
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, &
      "solve_fwd_update_task_worker", timer)
#endif

    ! Establish variables describing block
    blkm        = fkeep%sbc(blk)%blkm
    blkn        = fkeep%sbc(blk)%blkn
    blk_sa      = fkeep%sbc(blk)%sa
    dblk        = fkeep%sbc(blk)%dblk
    bcol        = fkeep%sbc(blk)%bcol
    ldy         = fkeep%sbc(dblk)%ldu
    ldx         = fkeep%sbc(blk)%ldu

    p_lcol      => fkeep%lfact(bcol)%lcol(blk_sa : blk_sa + blkn * blkm - 1)
    p_y         => fkeep%sbc(dblk)%p_upd(1 : ldy, 1 : nrhs)
    p_xlocal    => fkeep%sbc(blk)%p_upd(1 : ldx, 1 : nrhs)
    p_bc        => fkeep%sbc
    p_dep       => fkeep%sbc(blk)%fwd_dep
    p_wdep      => fkeep%sbc(blk)%fwd_wdep

    nftask      = 0
    ndep        = size(p_dep)
    nwdep       = size(p_wdep)
    chunk       = 0
    ndep_lvl    = 0

    if(ndep .eq. 0) then
#include "spllt_fwd_update_nodep.F90"
    else

      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      chunk = fkeep%chunk
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
#include "include/spllt_fwd_update_cases.F90"

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

#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%ntask_submitted(1, nftask)

  end subroutine solve_fwd_update_task_worker



  subroutine solve_bwd_block_task(task_manager, dblk, nrhs, &
      n, rhs, fkeep, trace_id)
    use spllt_data_mod
    implicit none

    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: dblk ! Index of diagonal block
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n
    real(wp), target,           intent(inout) :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, optional,          intent(in)    :: trace_id

    integer :: id
    
    id = task_manager%trace_ids(trace_bwd_block_pos)
    if(present(trace_id)) then
      id = trace_id
    end if

    call solve_bwd_block_task_worker(task_manager, dblk, nrhs, &
      n, rhs, fkeep, id)

  end subroutine solve_bwd_block_task

  subroutine solve_bwd_block_task_worker(task_manager, dblk, nrhs, &
      n, rhs, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use spllt_solve_dep_mod
    use timer_mod
    use omp_lib, ONLY : omp_get_thread_num
    implicit none

    type (task_manager_omp_t), intent(inout)  :: task_manager
    integer, intent(in)                       :: dblk ! Index of diagonal block
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: n
    real(wp), target, intent(inout)           :: rhs(n, nrhs)
    type(spllt_fkeep), target, intent(inout)  :: fkeep
    integer, intent(in)                       :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: blkm, blkn ! Block dimension
    integer                     :: bcol
    integer                     :: node
    integer                     :: i, j, r, chunkth
    integer                     :: ndep
    integer                     :: threadID
    integer                     :: ldy

    real(wp),             pointer :: p_y(:,:)
    real(wp),             pointer :: p_rhs(:,:)
    integer,              pointer :: p_index(:)
    integer,              pointer :: p_order(:)
    real(wp),             pointer :: p_lcol(:)
    type(spllt_sblock_t), pointer :: p_bc(:)
    integer,              pointer :: p_dep(:)

    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted
    integer :: nftask
    double precision :: flops

    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_block_task_worker", timer)
#endif

    ! Get block info
    node      = fkeep%sbc(dblk)%node
    blkn      = fkeep%sbc(dblk)%blkn
    blkm      = fkeep%sbc(dblk)%blkm
    sa        = fkeep%sbc(dblk)%sa
    bcol      = fkeep%sbc(dblk)%bcol ! Current block column
    ldy       = fkeep%sbc(dblk)%ldu
    
    p_index     => fkeep%sbc(dblk)%p_index
    p_order     => fkeep%porder
   !print *, "Current block column", bcol
   !print *, "sa = ", sa, "to sa + n * n - 1 = ", sa + blkn*blkn-1
    p_lcol      => fkeep%lfact(bcol)%lcol(sa : sa + blkn * blkn - 1)
    p_y         => fkeep%sbc(dblk)%p_upd(1 : ldy, 1 : nrhs)
    p_rhs       => rhs(1 : n, 1 : nrhs)
    p_bc        => fkeep%sbc
    p_dep       => fkeep%sbc(dblk)%bwd_dep

    nftask    = 0
    ndep      = size(p_dep)
    chunk     = 0
    ndep_lvl  = 0

    if(ndep .eq. 0) then
#include "spllt_bwd_block_nodep.F90"
    else
      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      chunk = fkeep%chunk
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
#include "include/spllt_bwd_block_cases.F90"

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
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%ntask_submitted(1, nftask)
  end subroutine solve_bwd_block_task_worker

  

  subroutine solve_bwd_update_task(task_manager, blk, node, nrhs, & 
      n, rhs, fkeep, trace_id)
    use spllt_data_mod
    implicit none

    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: blk  ! Index of block 
    integer,                    intent(in)    :: node 
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n
    real(wp), target,           intent(in)    :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, optional,          intent(in)    :: trace_id

    integer :: id

    id = task_manager%trace_ids(trace_bwd_update_pos)

    if(present(trace_id)) then
      id = trace_id
    end if

    call solve_bwd_update_task_worker(task_manager, blk, node, nrhs, &
      n, rhs, fkeep, id)

  end subroutine solve_bwd_update_task
    
  subroutine solve_bwd_update_task_worker(task_manager, blk, node, nrhs, &
      n, rhs, fkeep, trace_id)
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
    integer,                    intent(in)    :: n
    real(wp), target,           intent(in)    :: rhs(n, nrhs)
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer, intent(in)                       :: trace_id
    
    ! Block info
    integer                     :: blkm, blkn         ! Block dimension
    integer                     :: blk_sa, blk_en
    integer                     :: sa
    integer                     :: bcol, nbcol
    integer                     :: threadID
    integer                     :: ndep, nwdep
    integer                     :: dblk
    integer                     :: j, i
    integer                     :: chunk, chunk_size
    integer                     :: ndep_lvl, lvl, nchunk
    integer                     :: beta, alpha
    integer                     :: ldy, ldx
    logical                     :: all_task_submitted
    integer                     :: nftask   ! #fake tasks inserted 
                                            ! into the runtime
    double precision            :: flops

    real(wp),             pointer :: p_lcol(:)
    integer,              pointer :: p_dep(:)
    integer,              pointer :: p_wdep(:)
    real(wp)         ,    pointer :: p_xlocal(:,:)
    real(wp)         ,    pointer :: p_y(:,:)
    type(spllt_sblock_t), pointer :: p_bc(:)

   !type(spllt_timer_t), save   :: timer
    type(spllt_timer_t), target, save :: timer
    type(spllt_timer_t), pointer      :: p_timer
    p_timer => timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_update_task_worker", timer)
#endif

    ! Establish variables describing block
    blkn    = fkeep%sbc(blk)%blkn
    blkm    = fkeep%sbc(blk)%blkm
    blk_sa  = fkeep%nodes(node)%sblk_sa
    blk_en  = fkeep%nodes(node)%sblk_en
    sa      = fkeep%sbc(blk)%sa
    bcol    = fkeep%sbc(blk)%bcol
    dblk    = fkeep%sbc(blk)%dblk
    ldy     = fkeep%sbc(dblk)%ldu
    ldx     = fkeep%sbc(blk)%ldu

    p_lcol      => fkeep%lfact(bcol)%lcol(sa : sa + blkn * blkm - 1)
    p_bc        => fkeep%sbc
    p_y         => fkeep%sbc(dblk)%p_upd(1 : ldy, 1 : nrhs)
    p_xlocal    => fkeep%sbc(blk)%p_upd(1 : ldx, 1 : nrhs)
    p_dep       => fkeep%sbc(blk)%bwd_dep
    p_wdep      => fkeep%sbc(blk)%bwd_wdep

    nftask    = 0
    ndep      = size(p_dep)
    nwdep     = size(p_wdep)
    chunk     = 0
    ndep_lvl  = 0

    nbcol     = fkeep%sbc(blk_en)%bcol - fkeep%sbc(blk_sa)%bcol + 1

    if(ndep .eq. 0) then
#include "spllt_bwd_update_nodep.F90"
    else

      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      chunk = fkeep%chunk
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
#include "include/spllt_bwd_update_cases.F90"

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
#if defined(SPLLT_TIMER_TASKS_SUBMISSION) || defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
    call task_manager%ntask_submitted(1, nftask)
  end subroutine solve_bwd_update_task_worker



  subroutine solve_fwd_subtree_task(task_manager, nrhs, rhs, n, fkeep, tree)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_seq_mod
    implicit none

    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    integer,                    intent(in)    :: n
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree

    call solve_fwd_subtree_task_worker(task_manager, nrhs, rhs, n, &
      fkeep, tree)
  end subroutine solve_fwd_subtree_task



  subroutine solve_fwd_subtree_task_worker(task_manager, nrhs, rhs, n, &
      fkeep, tree)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_seq_mod
    implicit none

    type (task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),          target,  intent(inout) :: rhs(n, nrhs)
    integer,                    intent(in)    :: n
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
  
    integer                     :: trace_id
    real(wp), pointer           :: p_rhs(:,:)
    type(spllt_sblock_t), pointer  :: p_bc(:)
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

    p_rhs       => rhs(1 : n, 1 : nrhs)
    p_bc        => fkeep%sbc
    sa          = tree%node_sa
    en          = tree%node_en
    blk_en      = fkeep%nodes(en)%sblk_en
    call task_manager%semifork(sub_task_manager)

    !$omp task                                        &
!   !$omp default(none)                               &
    !$omp depend(inout: p_bc(blk_en))                 &
    !$omp firstprivate(p_rhs)                         &
    !$omp firstprivate(nrhs)                          &
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

    call sub_task_manager%incr_nrun()

    do i = sa, en
      call solve_fwd_node(nrhs, p_rhs, n, fkeep, i, &
        sub_task_manager, no_trace)
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



  subroutine solve_bwd_subtree_task(task_manager, nrhs, rhs, n, fkeep, tree)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_seq_mod
    implicit none

    class(task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    integer,                    intent(in)    :: n
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree

    call solve_bwd_subtree_task_worker(task_manager, nrhs, rhs, n, &
      fkeep, tree)
  end subroutine solve_bwd_subtree_task



  subroutine solve_bwd_subtree_task_worker(task_manager, nrhs, rhs, n, &
      fkeep, tree)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_seq_mod
    implicit none

    type (task_manager_omp_t),  intent(inout) :: task_manager
    integer,                    intent(in)    :: nrhs ! Number of RHS
    real(wp),          target,  intent(inout) :: rhs(n, nrhs)
    integer,                    intent(in)    :: n
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    type(spllt_tree_t),         intent(in)    :: tree
  
    integer,  pointer           :: p_dep(:)
    real(wp), pointer           :: p_rhs(:,:)
    type(spllt_sblock_t), pointer  :: p_bc(:)
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
    call spllt_open_timer(task_manager%workerID, &
      "solve_bwd_subtree_task_worker", timer)
#endif

    p_rhs       => rhs(1 : n, 1 : nrhs)
    p_bc        => fkeep%sbc
    sa          = tree%node_sa
    en          = tree%node_en
    blk_en      = fkeep%nodes(en)%sblk_en
    call task_manager%semifork(sub_task_manager)

    nftask      = 0
    p_dep       => fkeep%sbc(blk_en)%bwd_dep
    ndep        = size(p_dep)

    chunk = fkeep%chunk
    ndep_lvl = ndep ! #dep local to the lvl

    if(ndep .eq. 0) then
#include "spllt_bwd_node_nodep.F90"
    else

      chunk = 10 ! Do not use chunk = 1 ; a non-sence
      chunk = fkeep%chunk
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
#include "include/spllt_bwd_node_cases.F90"

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
