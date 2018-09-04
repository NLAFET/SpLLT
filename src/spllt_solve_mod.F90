module spllt_solve_mod

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
  subroutine spllt_solve_one_double(fkeep, options, x, job, info)
    use spllt_mod
    use spllt_data_mod
    use utils_mod
    use task_manager_omp_mod
    implicit none

    type(spllt_fkeep),    intent(in)  :: fkeep          ! Factorization data
    type(spllt_options),  intent(in)  :: options        ! User-supplied options
    ! For details of fkeep, control, info : see derived type description        
    real(wp),           intent(inout) :: x(fkeep%n)     ! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i) is the corresponding component of the right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i) holds solution for variable i.
    integer, optional,    intent(in)  :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)
    type(spllt_inform),   intent(out) :: info

    integer                   :: st
    integer                   :: solve_step   ! Selector step
    type(task_manager_omp_t)  :: task_manager

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    !!!!!!!!!!!!!!!!!!!!!
    ! Check optional parameters
    !
    if(present(job)) then
      solve_step = job
    else
      solve_step = 0
    end if

    call task_manager%init(stat = st)

    call spllt_solve_mult_double_worker(fkeep, options, 1, x, &
      solve_step, task_manager, info)

    call spllt_wait()

    !!!!!!!!!!!!!!!!!!!!!
    ! Desallocation
    !
    call task_manager%deallocate()
    
  end subroutine spllt_solve_one_double



  subroutine spllt_solve_mult_double(fkeep, options, nrhs, x, job, info)
    use spllt_mod
    use spllt_data_mod
    use utils_mod
    use task_manager_omp_mod
    implicit none

    type(spllt_fkeep),    intent(in)  :: fkeep          ! Factorization data
    type(spllt_options),  intent(in)  :: options        ! User-supplied options
    integer,              intent(in)  :: nrhs           ! Number of RHS
    ! For details of fkeep, control, info : see derived type description        
    real(wp),           intent(inout) :: x(fkeep%n,nrhs)! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i, j) is the corresponding component of the j-th right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i, j) holds solution for variable i of rhs-th right-hand side
    integer, optional,    intent(in)  :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)
    type(spllt_inform),   intent(out) :: info
  
    integer                   :: st
    integer                   :: solve_step   ! Selector step
    type(task_manager_omp_t)  :: task_manager

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    !!!!!!!!!!!!!!!!!!!!!
    ! Check optional parameters
    !
    if(present(job)) then
      solve_step = job
    else
      solve_step = 0
    end if

    call task_manager%init(stat = st)
    
    call spllt_solve_mult_double_worker(fkeep, options, nrhs, x, &
      solve_step, task_manager, info)

    call spllt_wait()

    !!!!!!!!!!!!!!!!!!!!!
    ! Desallocation
    !
    call task_manager%deallocate()

  end subroutine spllt_solve_mult_double




  subroutine spllt_solve_mult_double_worker(fkeep, options, nrhs, x, &
      job, task_manager, info)
    use spllt_data_mod
    use timer_mod
    use utils_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep)                 :: fkeep          ! Factorization data
    type(spllt_options),  intent(in)  :: options        ! User-supplied options
    integer,              intent(in)  :: nrhs           ! Number of RHS
    ! For details of fkeep, control, info : see derived type description        
    real(wp), target,   intent(inout) :: x(fkeep%n,nrhs)! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i, j) is the corresponding component of the j-th right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i, j) holds solution for variable i of rhs-th right-hand side
    integer, intent(in)  :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)
    class(task_manager_base), intent(inout), target :: task_manager
    type(spllt_inform),   intent(out) :: info

    character(len=30)         :: context  ! Name of the subroutine
    type(spllt_timer_t), save :: timer

    context = 'spllt_solve_mult_double_worker'
   !call spllt_open_timer(task_manager%workerID, &
   !  "spllt_solve_mult_double_worker", timer)

    ! immediate return if n = 0
    if (fkeep%n == 0) return


    select case(job)

      case(0)
        call solve_fwd(nrhs, x, fkeep%n, fkeep, task_manager, info)

        call solve_bwd(nrhs, x, fkeep%n, fkeep, task_manager, info)

      case(1)
        call solve_fwd(nrhs, x, fkeep%n, fkeep, task_manager, info)

      case(2)
        call solve_bwd(nrhs, x, fkeep%n, fkeep, task_manager, info)

      case default
        info%flag = SPLLT_WARNING_PARAM_VALUE
        write (0, '(a, i2, a, i4)') "Unknown requested job = ", job, &
          " returned code : ", info%flag
        return
    end select

   !call spllt_close_timer(task_manager%workerID, timer)
  end subroutine spllt_solve_mult_double_worker





  !*************************************************
  !
  ! Forward solve routine
  subroutine solve_fwd(nrhs, rhs, n, fkeep, task_manager, info)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use timer_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n    ! Order of the system
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    class(task_manager_base ),  intent(inout) :: task_manager
    type(spllt_inform),         intent(out)   :: info

    integer                     :: node
    integer                     :: tree_num
    integer                     :: fwd_submit_tree_id
    integer                     :: fwd_submit_id
    double precision            :: nflop_sa, nflop_en
    character(len=30)           :: context = 'solve_fwd'
#if defined(SPLLT_TIMER_TASKS)
    type(spllt_timer_t), save   :: timer
#endif

!   context = 'solve_fwd'
#if defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, "solve_fwd", timer)
    call task_manager%get_nflop_performed(nflop_sa)
#endif
#if defined(SPLLT_OMP_TRACE)
    fwd_submit_id = task_manager%trace_ids(trace_fwd_submit_pos)
    call trace_event_start(fwd_submit_id, -1)
#endif

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tic("Submit tasks", 4, task_manager%workerID, timer)
#endif


    !!!!!!!!!!!!!!!!!!!!!!
    ! Submit subtrees to the runtime system
    !
#if defined(SPLLT_OMP_TRACE)
    fwd_submit_tree_id = task_manager%trace_ids(trace_fwd_submit_tree_pos)
    call trace_event_start(fwd_submit_tree_id, -2)
#endif
    do tree_num = 1, size(fkeep%trees)
      call task_manager%solve_fwd_subtree_task(nrhs, rhs, n, &
        fkeep, fkeep%trees(tree_num))
    end do
#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(fwd_submit_tree_id, -2)
#endif


    !!!!!!!!!!!!!!!!!!!!!
    ! Submit the treatment of the remaining nodes
    do node = 1, fkeep%info%num_nodes

      if(fkeep%small(node) .ne. 0) then
        cycle
      end if

      call solve_fwd_node(nrhs, rhs, n, fkeep, node, task_manager)

    end do
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tac(4, task_manager%workerID, timer)
#endif

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(fwd_submit_id, -1)
#endif
    
   !call task_manager%print("fwd end of submitted task", 0)

#if defined(SPLLT_TIMER_TASKS)
    call task_manager%get_nflop_performed(nflop_en)
    call spllt_close_timer(task_manager%workerID, timer, nflop_en - nflop_sa)
#endif

  end subroutine solve_fwd



  subroutine solve_bwd(nrhs, rhs, n, fkeep, task_manager, info)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
 !$ use omp_lib, ONLY : omp_lock_kind
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n    ! Order of the system
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    class(task_manager_base),   intent(inout) :: task_manager
    type(spllt_inform),         intent(out)   :: info

    integer                     :: node
    integer                     :: tree_num
    integer                     :: bwd_submit_tree_id
    integer                     :: bwd_submit_id
    double precision            :: nflop_sa, nflop_en
    character(len=30)           :: context = 'solve_bwd'
#if defined(SPLLT_TIMER_TASKS)
    type(spllt_timer_t), save   :: timer
#endif

#if defined(SPLLT_TIMER_TASKS)
    call spllt_open_timer(task_manager%workerID, "solve_bwd", timer)
    call task_manager%get_nflop_performed(nflop_sa)
#endif

#if defined(SPLLT_OMP_TRACE)
    bwd_submit_id = task_manager%trace_ids(trace_bwd_submit_pos)
    call trace_event_start(bwd_submit_id, -1)
#endif

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_tic("Submit tasks", 4, task_manager%workerID, timer)
#endif
    do node = fkeep%info%num_nodes, 1, -1
     !call print_node_solve(fkeep, node)

      if(fkeep%small(node) .eq. 0) then
        call solve_bwd_node(nrhs, rhs, n, fkeep, node, task_manager)
      else if(fkeep%small(node) .eq. 1) then
        tree_num = fkeep%assoc_tree(node)

        call task_manager%solve_bwd_subtree_task(nrhs, rhs, n, fkeep, &
          fkeep%trees(tree_num))
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

   !call task_manager%print("bwd end of execution task", 0)

#if defined(SPLLT_TIMER_TASKS)
    call task_manager%get_nflop_performed(nflop_en)
    call spllt_close_timer(task_manager%workerID, timer, nflop_en - nflop_sa)
#endif

  end subroutine solve_bwd


  !______________________________________________________
  ! Locked versions of the forward and backward solve
  ! Its purpose is to ensure that the submission of the tasks
  ! is not altered by the execution of some by the submitting
  ! master worker
  !
  subroutine solve_fwd_lock(nrhs, rhs, n, fkeep, task_manager, info)
    use spllt_data_mod
    use task_manager_mod
 !$ use omp_lib, ONLY : omp_lock_kind
    implicit none

    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n    ! Order of the system
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    class(task_manager_base ),  intent(inout) :: task_manager
    type(spllt_inform),         intent(out)   :: info

    type(spllt_sblock_t), pointer :: p_bc(:)
 !$ integer(kind=omp_lock_kind)   :: lock

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Lock the execution of the tasks until all the tasks
    ! have been submitted by the master
    p_bc => fkeep%sbc(:)
    !$ call omp_init_lock(lock)
    !$ call omp_set_lock(lock)

    !$omp task depend(inout: p_bc(1))   &
    !$omp shared(lock)                  &
    !$omp firstprivate(p_bc)

    !$ call omp_set_lock(lock) 
    print *, "Lock locked in task"
    !$ call omp_unset_lock(lock)

    !$omp end task

    call solve_fwd(nrhs, rhs, n, fkeep, task_manager, info)

    call task_manager%print("Release of lock in fwd", 0)
    !$ call omp_unset_lock(lock)
    !$ call omp_destroy_lock(lock)

  end subroutine solve_fwd_lock



  subroutine solve_bwd_lock(nrhs, rhs, n, fkeep, task_manager, info)
    use spllt_data_mod
    use task_manager_mod
 !$ use omp_lib, ONLY : omp_lock_kind
    implicit none

    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n    ! Order of the system
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    class(task_manager_base),   intent(inout) :: task_manager
    type(spllt_inform),         intent(out)   :: info

    integer                       :: num_nodes
    type(spllt_sblock_t), pointer :: p_bc(:)
 !$ integer(kind=omp_lock_kind) :: lock

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Lock the execution of the tasks until all the tasks
    ! have been submitted by the master
    p_bc => fkeep%sbc(:)
    num_nodes = fkeep%info%num_nodes
    !$ call omp_init_lock(lock)
    !$ call omp_set_lock(lock)

    !$omp task depend(inout: p_bc(num_nodes))   &
    !$omp shared(lock)                          &
    !$omp firstprivate(p_bc, num_nodes)

    !$ call omp_set_lock(lock) 
    print *, "[BWD] Lock locked in task"
    !$ call omp_unset_lock(lock)

    !$omp end task

    call solve_bwd(nrhs, rhs, n, fkeep, task_manager, info)

    call task_manager%print("Release of lock in bwd", 0)
    !$ call omp_unset_lock(lock)
    !$ call omp_destroy_lock(lock)

  end subroutine solve_bwd_lock


end module spllt_solve_mod
