module spllt_solve_mod
  !$ use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads

   interface spllt_solve
      module procedure spllt_solve_one_double
      module procedure spllt_solve_mult_double
      module procedure spllt_solve_mult_double_worker
   end interface

contains
  
  !*************************************************
  !
  !
  ! Solve phase. simplified interface for a single rhs
  !
  subroutine spllt_solve_one_double(fkeep, options, order, x, info, job, &
      task_manager)
    use spllt_data_mod
    use utils_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep),    intent(in)  :: fkeep          ! Factorization data
    type(spllt_options),  intent(in)  :: options        ! User-supplied options
    integer,              intent(in)  :: order(fkeep%n) ! pivot order
    ! For details of fkeep, control, info : see derived type description        
    real(wp),           intent(inout) :: x(fkeep%n)     ! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i) is the corresponding component of the right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i) holds solution for variable i.
    type(spllt_inform),   intent(out) :: info
    integer, optional,    intent(in)  :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)
    class(task_manager_base), intent(inout), optional, target  :: task_manager

    integer                        :: n            ! #rows
    integer                        :: solve_step   ! Selector step
    integer                        :: st           ! stat parameter
    real(wp),             pointer  :: work(:)
    class(task_manager_base), pointer  :: p_ltask_manager

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    p_ltask_manager => null()
    n               = fkeep%n

    !!!!!!!!!!!!!!!!!!!!!
    ! Check optional parameters
    !
    if(present(job)) then
      solve_step = job
    else
      solve_step = 0
    end if

    if(.not.present(task_manager)) then
     !allocate(p_ltask_manager, stat = st)
      
      if(st .ne. 0) then
        write (0,*) "[>Allocation::error] Can not create a task manager"
        stop
      end if

      call p_ltask_manager%init(stat = st)
    else
      p_ltask_manager => task_manager
    end if
    
    !!!!!!!!!!!!!!!!!!!!!
    ! worker
    !
    allocate(work(n + (fkeep%maxmn + n) * p_ltask_manager%nworker), stat = st)
    call p_ltask_manager%incr_alloc(st)

    call spllt_solve_mult_double_worker(fkeep, options, order, 1, x, info, &
      solve_step, work, p_ltask_manager)

    !!!!!!!!!!!!!!!!!!!!!
    ! Desallocation
    !
    deallocate(work)

    if(.not.present(task_manager)) then
      call p_ltask_manager%deallocate()
      deallocate(p_ltask_manager)
    end if
    
  end subroutine spllt_solve_one_double

  



  subroutine spllt_solve_mult_double(fkeep, options, order, nrhs, x, info, &
      job, task_manager)
    use spllt_data_mod
    use utils_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep),    intent(in)  :: fkeep          ! Factorization data
    type(spllt_options),  intent(in)  :: options        ! User-supplied options
    integer,              intent(in)  :: order(fkeep%n) ! pivot order
    integer,              intent(in)  :: nrhs           ! Number of RHS
    ! For details of fkeep, control, info : see derived type description        
    real(wp),           intent(inout) :: x(fkeep%n,nrhs)! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i, j) is the corresponding component of the j-th right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i, j) holds solution for variable i of rhs-th right-hand side
    type(spllt_inform),   intent(out) :: info
    integer, optional,    intent(in)  :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)
    class(task_manager_base), intent(inout), optional, target  :: task_manager
  
    integer                       :: n            ! #rows
    integer                       :: solve_step   ! Selector step
    integer                       :: st           ! stat parameter
    integer(long)                 :: worksize
    real(wp),             pointer :: work(:)
    class(task_manager_base), pointer :: p_ltask_manager

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    p_ltask_manager => null()
    n               = fkeep%n

    !!!!!!!!!!!!!!!!!!!!!
    ! Check optional parameters
    !
    if(present(job)) then
      solve_step = job
    else
      solve_step = 0
    end if

    if(.not.present(task_manager)) then
     !allocate(p_ltask_manager)

      if(st .ne. 0) then
        write (0,*) "[>Allocation::error] Can not create a task manager"
        stop
      end if

      call p_ltask_manager%init(stat = st)
    else
      p_ltask_manager => task_manager
    end if
    
    !!!!!!!!!!!!!!!!!!!!!
    ! worker
    !

   !allocate(work(n*nrhs + (fkeep%maxmn + n) * nrhs * &
   !  p_ltask_manager%nworker), stat = st)
    call spllt_solve_workspace_size(fkeep, p_ltask_manager%nworker, &
      nrhs, worksize)
    allocate(work(worksize), stat = st)
    call p_ltask_manager%incr_alloc(st)
    work = 0.0
    
    call spllt_solve_mult_double_worker(fkeep, options, order, nrhs, x, &
      info, solve_step, work, p_ltask_manager)

    !!!!!!!!!!!!!!!!!!!!!
    ! Desallocation
    !
    deallocate(work)

    if(.not.present(task_manager)) then
      call p_ltask_manager%deallocate()
      deallocate(p_ltask_manager)
    end if
  end subroutine spllt_solve_mult_double




!#define SPLLT_ILEAVE
  subroutine spllt_solve_mult_double_worker(fkeep, options, order, nrhs, x, &
      info, job, workspace, task_manager)
    use spllt_data_mod
    use spllt_solve_dep_mod, ONLY : spllt_compute_rhs_block
    use utils_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep),    intent(in)  :: fkeep          ! Factorization data
    type(spllt_options),  intent(in)  :: options        ! User-supplied options
    integer,              intent(in)  :: order(fkeep%n) ! pivot order
    integer,              intent(in)  :: nrhs           ! Number of RHS
    ! For details of fkeep, control, info : see derived type description        
    real(wp),           intent(inout) :: x(fkeep%n,nrhs)! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i, j) is the corresponding component of the j-th right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i, j) holds solution for variable i of rhs-th right-hand side
    type(spllt_inform),   intent(out) :: info
    integer, intent(in)  :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)
!   real(wp), dimension(fkeep%n, nrhs), intent(out)  :: work
    real(wp),             intent(out),   target :: workspace(:)
    class(task_manager_base), intent(inout), target :: task_manager

    integer           :: r
    integer           :: j        ! Iterator
    integer           :: n        ! Order of the system
    character(len=30) :: context  ! Name of the subroutine
    real(wp), pointer :: x_tmp(:,:)
    real(wp), pointer :: work(:)
    integer(long)     :: size_nrhs, size_work
#if defined(SPLLT_ILEAVE)
    integer               :: ldr, bdr
    integer               :: st
    real(wp), allocatable :: xil(:)
    real(wp), allocatable :: norm(:)
    integer,  allocatable :: rhsPtr(:)
#endif

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    n = fkeep%n
    context = 'spllt_solve_mult_double_worker'

   !work(1 : n, 1 : nrhs) => workspace(1 : n * nrhs)
   !work2(1 : (fkeep%maxmn + n) * nrhs * task_manager%nworker) =>      &
   !  workspace(n * nrhs + 1 : nrhs * (n + (fkeep%maxmn + n) *      &
   !  task_manager%nworker))

    size_nrhs = int(n, long) * nrhs
    size_work = int(fkeep%maxmn + n, long) * int(nrhs, long)

    x_tmp(1 : n, 1 : nrhs) => workspace(1 : size_nrhs)
    work(1 : size_work * task_manager%nworker) =>                        &
      workspace(size_nrhs + 1 : nrhs * (n + int(fkeep%maxmn + n, long) *  &
      task_manager%nworker))

    select case(job)
      case(0)
       !
       ! Reorder x
       !
       !TODO Add subroutine to interleave and order at the same time
        do j = 1, n
          !work(order(j),:) = x(j,:)
           x_tmp(order(j),:) = x(j,:)
        end do
 !!     do r = 1, nrhs
 !!       call print_darray("Original RHS", n, x_tmp(:, r), 1)
 !!     end do

#if defined(SPLLT_ILEAVE)

      ! bdr = fkeep%nodes(1)%nb
      ! ldr = nrhs * bdr
!     ! print *, "Allocate xil of size ", size_nrhs
      ! allocate(xil(size_nrhs))
      ! allocate(norm(nrhs))

      ! call pack_rhs(nrhs, x_tmp, n, ldr, bdr, xil, st)

      ! call print_darray("Interleaved RHS", int(size_nrhs), xil, 1)

      ! call unpack_rhs(nrhs, xil, n, ldr, bdr, x, st)

      ! call vector_norm_2(nrhs, x_tmp - x, norm)
      ! print *, "Norm ", norm

      ! deallocate(xil)
      ! deallocate(norm)
        allocate(xil(size_nrhs))
       !allocate(norm(nrhs))
        call pack_rhs(nrhs, x_tmp, n, fkeep%rhsPtr, xil, st)

       !call print_darray("Interleaved RHS", int(size_nrhs), xil, 1)

       !call unpack_rhs(nrhs, xil, n, fkeep%rhsPtr, x, st)

       !do r = 1, nrhs
       !  call print_darray("Unpack RHS", n, x(:, r), 1)
       !end do

       !call vector_norm_2(nrhs, x_tmp - x, norm)
       !print *, "Norm ", norm
       !stop

 !!     print *, "RHSPTR", fkeep%rhsPtr

        ! Forward solve
        call solve_fwd_ileave(nrhs, xil, n, ldr, bdr, fkeep, work, task_manager)

!       call print_darray("Interleaved returned solution", int(size_nrhs), &
!         xil, 1)

       !call unpack_rhs(nrhs, xil, n, ldr, bdr, x_tmp, st)
        call unpack_rhs(nrhs, xil, n, fkeep%rhsPtr, x_tmp, st)
#else
        ! Forward solve
        call solve_fwd(nrhs, x_tmp, n, fkeep, work, task_manager)

#endif
       !!$omp taskwait
 !!     do r = 1, nrhs
 !!       call print_darray("Solution returned by fwd", n, x_tmp(:, r), 1)
 !!     end do
    
       !stop

        ! Backward solve
        call solve_bwd(nrhs, x_tmp, n, fkeep, work, task_manager)

       !do r = 1, nrhs
       !  call print_darray("Solution returned by bwd", n, x_tmp(:, r), 0)
       !end do
       !
       ! Reorder soln
       !
      !$omp taskwait
        do j = 1, n
           x(j,:) = x_tmp(order(j),:)
        end do

     !case(1)
     !  x_tmp = x
     !  do j = 1, n
     !     x(order(j),:) = x_tmp(j,:)
     !  end do

     !  call solve_fwd(nrhs, x, n, fkeep, work, task_manager)
     !!$omp taskwait
     !
     !case(2)

     !  call solve_bwd(nrhs, x, n, fkeep, work, task_manager)
     !!$omp taskwait
     !  x_tmp = x
     !  do j = 1, n
     !     x(j,:) = x_tmp(order(j),:)
     !  end do

      case default
        info%flag = SPLLT_WARNING_PARAM_VALUE
        write (0, '(a, i2, a, i4)') "Unknown requested job = ", job, &
          " returned code : ", info%flag
        return
    end select

  end subroutine spllt_solve_mult_double_worker





  !*************************************************
  !
  ! Forward solve routine
  subroutine solve_fwd(nrhs, rhs, ldr, fkeep, workspace, task_manager)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use timer_mod
    use task_manager_mod
 !$ use omp_lib, ONLY : omp_lock_kind
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    real(wp), target,           intent(out)   :: workspace(:)
    class(task_manager_base ),  intent(inout) :: task_manager

    integer                     :: node
    integer                     :: fwd_submit_tree_id
    integer                     :: fwd_submit_id
    integer                     :: nworker
    real(wp), pointer           :: xlocal(:,:)    ! update_buffer workspace
    real(wp), pointer           :: rhs_local(:,:) ! update_buffer workspace
    type(spllt_block), pointer  :: p_bc(:)
    type(spllt_timer_t), save   :: timer
 !$ integer(kind=omp_lock_kind) :: lock
    integer                     :: tree_num
    integer(long)               :: size_rhs_local, size_xlocal
    double precision            :: nflop_sa, nflop_en

#if defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, "solve_fwd", timer)
#endif
    call task_manager%get_nflop_performed(nflop_sa)


    nworker       = task_manager%nworker
#if defined(SPLLT_OMP_TRACE)
    fwd_submit_id = task_manager%trace_ids(trace_fwd_submit_pos)
    call trace_event_start(fwd_submit_id, -1)
#endif

   !xlocal(1 : fkeep%maxmn * nrhs, 0 : nworker - 1) => workspace(1 : &
   !  fkeep%maxmn * nrhs * nworker)
   !rhs_local(1 : ldr * nrhs, 0 : nworker - 1) => workspace(fkeep%maxmn * nrhs &
   !  * nworker + 1 : (fkeep%maxmn + ldr) * nrhs * nworker)

    size_xlocal     = int(fkeep%maxmn, long) * nrhs
    size_rhs_local  = int(fkeep%maxmn + ldr, long) * int(nrhs, long)

    xlocal(1 : size_xlocal, 0 : nworker - 1) => workspace(1 : &
      size_xlocal * nworker)
    rhs_local(1 : int(ldr, long) * nrhs, 0 : nworker - 1) =>  &
      workspace(size_xlocal * nworker + 1 : size_rhs_local * nworker)

#if defined(SOLVE_TASK_LOCKED)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Lock the execution of the tasks that will then be submitted by the master
    p_bc => fkeep%bc(:)
 !$ call omp_init_lock(lock)
 !$ call omp_set_lock(lock)

    !$omp task depend(inout: p_bc(1))   &
    !$omp shared(lock)                  &
    !$omp firstprivate(p_bc)

 !$ call omp_set_lock(lock) 
    print *, "Lock locked in task"
 !$ call omp_unset_lock(lock)

    !$omp end task
#endif

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tic("Submit tasks", 4, task_manager%workerID, timer)
#endif

#if defined(SPLLT_OMP_TRACE)
    fwd_submit_tree_id = task_manager%trace_ids(trace_fwd_submit_tree_pos)
    call trace_event_start(fwd_submit_tree_id, -2)
#endif

    do tree_num = 1, size(fkeep%trees)
      call task_manager%solve_fwd_subtree_task(nrhs, rhs, ldr, fkeep, &
        fkeep%trees(tree_num), xlocal, rhs_local)
    end do

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(fwd_submit_tree_id, -2)
#endif

    do node = 1, fkeep%info%num_nodes

      if(fkeep%small(node) .ne. 0) then
        cycle
      end if

      call solve_fwd_node(nrhs, rhs, ldr, fkeep, node, xlocal, rhs_local, &
        task_manager)
    end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tac(4, task_manager%workerID, timer)
#endif

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(fwd_submit_id, -1)
#endif

#if defined(SOLVE_TASK_LOCKED)
call task_manager%print("fwd end of submitted task", 0)
 !$ call omp_unset_lock(lock)
 !$ call omp_destroy_lock(lock)
#endif
    
    call task_manager%print("fwd end of submitted task", 0)

    call task_manager%get_nflop_performed(nflop_en)
#if defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer, nflop_en - nflop_sa)
#endif

  end subroutine solve_fwd



  subroutine solve_bwd(nrhs, rhs, ldr, fkeep, workspace, task_manager)
    use spllt_data_mod
    use spllt_solve_task_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
 !$ use omp_lib, ONLY : omp_lock_kind
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp),                   intent(inout) :: rhs(ldr, nrhs)
    real(wp), target,           intent(out)   :: workspace(:)
    class(task_manager_base),   intent(inout) :: task_manager

    ! Node info
    integer                     :: node
    integer                     :: bwd_submit_tree_id
    integer                     :: bwd_submit_id
    integer                     :: nworker
    real(wp), pointer           :: xlocal(:,:)    ! update_buffer workspace
    real(wp), pointer           :: rhs_local(:,:) ! update_buffer workspace
    type(spllt_block), pointer  :: p_bc(:)
    integer                     :: tree_num
    type(spllt_timer_t), save   :: timer
    integer(long)               :: size_rhs_local, size_xlocal
    double precision            :: nflop_sa, nflop_en
 !$ integer(kind=omp_lock_kind) :: lock

#if defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, "solve_bwd", timer)
#endif
    call task_manager%get_nflop_performed(nflop_sa)

    nworker       = task_manager%nworker

#if defined(SPLLT_OMP_TRACE)
    bwd_submit_id = task_manager%trace_ids(trace_bwd_submit_pos)
    call trace_event_start(bwd_submit_id, -1)
#endif

   !xlocal(1 : fkeep%maxmn * nrhs, 0 : nworker - 1) => workspace(1 : &
   !  fkeep%maxmn * nrhs * nworker)
   !rhs_local(1 : ldr * nrhs, 0 : nworker - 1) => workspace(fkeep%maxmn * nrhs &
   !  * nworker + 1 : (fkeep%maxmn + ldr) * nrhs * nworker)
    size_xlocal     = int(fkeep%maxmn, long) * nrhs
    size_rhs_local  = int(fkeep%maxmn + ldr, long) * int(nrhs, long)

    xlocal(1 : size_xlocal, 0 : nworker - 1) => workspace(1 : &
      size_xlocal * nworker)
    rhs_local(1 : int(ldr, long) * nrhs, 0 : nworker - 1) =>  &
      workspace(size_xlocal * nworker + 1 : size_rhs_local * nworker)


#if defined(SOLVE_TASK_LOCKED)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Lock the execution of the tasks that will then be submitted by the master
    p_bc => fkeep%bc(:)
 !$ call omp_init_lock(lock)
 !$ call omp_set_lock(lock)

    !$omp task depend(inout: p_bc(num_node))  &
    !$omp shared(lock)                        &
    !$omp firstprivate(p_bc, num_node)

 !$ call omp_set_lock(lock) 
    print *, "[BWD] Lock locked in task"
 !$ call omp_unset_lock(lock)

    !$omp end task
#endif

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tic("Submit tasks", 4, task_manager%workerID, timer)
#endif
    do node = fkeep%info%num_nodes, 1, -1

      if(fkeep%small(node) .eq. 0) then
!       print *, "Node ", node, " is not part of a tree"
        call solve_bwd_node(nrhs, rhs, ldr, fkeep, node, xlocal, rhs_local, &
          task_manager)
      else if(fkeep%small(node) .eq. 1) then
        tree_num = fkeep%assoc_tree(node)

        call task_manager%solve_bwd_subtree_task(nrhs, rhs, ldr, fkeep, &
          fkeep%trees(tree_num), xlocal, rhs_local)
      else
        cycle
      end if

    end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tac(4, task_manager%workerID, timer)
#endif


#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(bwd_submit_id, -1)
#endif
#if defined(SOLVE_TASK_LOCKED)
  call task_manager%print("bwd end of submitted task", 0)
 !$ call omp_unset_lock(lock)
 !$ call omp_destroy_lock(lock)
#endif

    call task_manager%print("bwd end of execution task", 0)

    call task_manager%get_nflop_performed(nflop_en)
#if defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer, nflop_en - nflop_sa)
#endif

  end subroutine solve_bwd
  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! InterLeave subroutines
  !
  subroutine solve_fwd_ileave(nrhs, rhs, n, ldr, bdr, fkeep, workspace, task_manager)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use timer_mod
    use utils_mod, ONLY : print_node, print_blk_index
    use task_manager_mod
 !$ use omp_lib, ONLY : omp_lock_kind
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n    ! Order of the system
    integer,                    intent(in)    :: ldr  ! Move block by block in rhs
    integer,                    intent(in)    :: bdr  ! Move RHS by RHS in rhs
    real(wp),                   intent(inout) :: rhs(n*nrhs)
    real(wp), target,           intent(out)   :: workspace(:)
    class(task_manager_base ),  intent(inout) :: task_manager

    integer                     :: nb
    integer                     :: ldu, bdu, tdu
    integer                     :: node
    integer                     :: fwd_submit_tree_id
    integer                     :: fwd_submit_id
    integer                     :: nworker
    real(wp), pointer           :: xlocal(:,:)    ! update_buffer workspace
    real(wp), pointer           :: rhs_local(:) ! update_buffer workspace
    type(spllt_block), pointer  :: p_bc(:)
    type(spllt_timer_t), save   :: timer
 !$ integer(kind=omp_lock_kind) :: lock
    integer                     :: tree_num
    integer(long)               :: size_rhs_local, size_xlocal
    double precision            :: nflop_sa, nflop_en

#if defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, "solve_fwd_ileave", timer)
#endif
    call task_manager%get_nflop_performed(nflop_sa)


    nworker       = task_manager%nworker
#if defined(SPLLT_OMP_TRACE)
    fwd_submit_id = task_manager%trace_ids(trace_fwd_submit_pos)
    call trace_event_start(fwd_submit_id, -1)
#endif

    size_xlocal     = int(fkeep%maxmn, long) * nrhs
    size_rhs_local  = int(fkeep%maxmn + n, long) * int(nrhs, long) ! Shifted

    xlocal(1 : size_xlocal, 0 : nworker - 1) => workspace(1 : &
      size_xlocal * nworker)
    rhs_local => workspace(size_xlocal * nworker + 1 : size_rhs_local * nworker)
    nb  = fkeep%nodes(1)%nb
    tdu = nworker
   !bdu = nb * tdu
   !ldu = nrhs * bdu

#if defined(SOLVE_TASK_LOCKED)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Lock the execution of the tasks that will then be submitted by the master
    p_bc => fkeep%bc(:)
 !$ call omp_init_lock(lock)
 !$ call omp_set_lock(lock)

    !$omp task depend(inout: p_bc(1))   &
    !$omp shared(lock)                  &
    !$omp firstprivate(p_bc)

 !$ call omp_set_lock(lock) 
    print *, "Lock locked in task"
 !$ call omp_unset_lock(lock)

    !$omp end task
#endif

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tic("Submit tasks", 4, task_manager%workerID, timer)
#endif

#if defined(SPLLT_OMP_TRACE)
    fwd_submit_tree_id = task_manager%trace_ids(trace_fwd_submit_tree_pos)
    call trace_event_start(fwd_submit_tree_id, -2)
#endif

   !do tree_num = 1, size(fkeep%trees)
   !  call task_manager%solve_fwd_subtree_task(nrhs, rhs, ldr, fkeep, &
   !    fkeep%trees(tree_num), xlocal, rhs_local)
   !end do

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(fwd_submit_tree_id, -2)
#endif

    do node = 1, fkeep%info%num_nodes

     !if(fkeep%small(node) .ne. 0) then
     !  cycle
     !end if

 !!   print *, "Submit node ", node
 !!   call print_node(fkeep, node)
 !!   call print_blk_index('', size(fkeep%nodes(node)%index), &
 !!       fkeep%nodes(node)%index, 1)
      call solve_fwd_node_ileave(nrhs, rhs, n, ldr, bdr, fkeep, node, xlocal, &
        rhs_local, ldu, bdu, tdu, task_manager)
    end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tac(4, task_manager%workerID, timer)
#endif

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(fwd_submit_id, -1)
#endif

#if defined(SOLVE_TASK_LOCKED)
call task_manager%print("fwd end of submitted task", 0)
 !$ call omp_unset_lock(lock)
 !$ call omp_destroy_lock(lock)
#endif
    
    call task_manager%print("fwd end of submitted task", 0)

    call task_manager%get_nflop_performed(nflop_en)
#if defined(SPLLT_TIMER_TASKS)
    call spllt_close_timer(task_manager%workerID, timer, nflop_en - nflop_sa)
#endif

  end subroutine solve_fwd_ileave



end module spllt_solve_mod
