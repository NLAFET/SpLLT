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

    allocate(work(n*nrhs + (fkeep%maxmn + n) * nrhs * &
      p_ltask_manager%nworker), stat = st)
    call p_ltask_manager%incr_alloc(st)
    
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





  subroutine spllt_solve_mult_double_worker(fkeep, options, order, nrhs, x, &
      info, job, workspace, task_manager)
    use spllt_data_mod
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

    integer           :: j        ! Iterator
    integer           :: n        ! Order of the system
    character(len=30) :: context  ! Name of the subroutine
    real(wp), pointer :: work(:,:)
    real(wp), pointer :: work2(:)

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    n = fkeep%n
    context = 'spllt_solve_mult_double_worker'

    work(1 : n, 1 : nrhs) => workspace(1 : n * nrhs)
    work2(1 : (fkeep%maxmn + n) * nrhs * task_manager%nworker) =>      &
      workspace(n * nrhs + 1 : nrhs * (n + (fkeep%maxmn + n) *      &
      task_manager%nworker))

    select case(job)
      case(0)
       !
       ! Reorder x
       !
        do j = 1, n
           work(order(j),:) = x(j,:)
        end do

        ! Forward solve
        call solve_fwd(nrhs, work, n, fkeep, work2, task_manager)

        ! Backward solve
        call solve_bwd(nrhs, work, n, fkeep, work2, task_manager)

       !
       ! Reorder soln
       !
        do j = 1, n
           x(j,:) = work(order(j),:)
        end do

      case(1)
        do j = 1, n
           work(order(j),:) = x(j,:)
        end do

        call solve_fwd(nrhs, work, n, fkeep, work2, task_manager)

        x = work
      
      case(2)
        work = x

        call solve_bwd(nrhs, work, n, fkeep, work2, task_manager)

        do j = 1, n
           x(j,:) = work(order(j),:)
        end do

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
    use spllt_solve_task_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_mod
    use task_manager_omp_mod
 !$ use omp_lib, ONLY : omp_lock_kind
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    real(wp), target,           intent(out)   :: workspace(:)
    class(task_manager_base ),  intent(inout) :: task_manager

    integer                 :: num_node
    integer                 :: node
    integer                 :: sa, en
    integer                 :: numcol, numrow ! #column/row in node 
    integer                 :: nc, nr         ! #block-column/block-row in node
    integer                 :: jj, ii
    integer                 :: dblk           ! Diagonal index 
    integer                 :: s_nb           ! Block size in node
    integer                 :: blk            ! Block index
    integer                 :: fwd_update_id
    integer                 :: fwd_block_id
    integer                 :: fwd_submit_id
    integer                 :: nworker
    real(wp), pointer       :: xlocal(:,:)    ! update_buffer workspace
    real(wp), pointer       :: rhs_local(:,:) ! update_buffer workspace
    type(spllt_block), pointer :: p_bc(:)
    type(spllt_timer), save :: timer
 !$ integer(kind=omp_lock_kind) :: lock

    call spllt_open_timer(task_manager%nworker, task_manager%workerID, &
      "solve_fwd", timer)

    nworker       = task_manager%nworker
#if defined(SPLLT_OMP_TRACE)
    fwd_submit_id = task_manager%trace_ids(trace_fwd_submit_pos)
    call trace_event_start(fwd_submit_id, -1)
#endif

    xlocal(1 : fkeep%maxmn * nrhs, 0 : nworker - 1) => workspace(1 : &
      fkeep%maxmn * nrhs * nworker)
    rhs_local(1 : ldr * nrhs, 0 : nworker - 1) => workspace(fkeep%maxmn * nrhs &
      * nworker + 1 : (fkeep%maxmn + ldr) * nrhs * nworker)

!   ! initialise rhs_local
!   rhs_local(:,:) = zero

!   call spllt_tic("Reset workspace", 1, task_manager%workerID, timer)
!   rhs_local(:,:) = zero
!   call spllt_tac(1, task_manager%workerID, timer)

!   call spllt_tic("Reset workspace in parallel", 5, task_manager%workerID,&
!     timer)
!   !$omp parallel
!   rhs_local(1 : ldr * nrhs, omp_get_thread_num()) = zero
!   !$omp end parallel
!   call spllt_tac(5, task_manager%workerID, timer)

    num_node = fkeep%info%num_nodes

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

    call spllt_tic("Submit tasks", 4, task_manager%workerID, timer)
    do node = 1, num_node

       ! Get node info
       s_nb   = fkeep%nodes(node)%nb
       sa     = fkeep%nodes(node)%sa
       en     = fkeep%nodes(node)%en
       numcol = en - sa + 1
       numrow = size(fkeep%nodes(node)%index)
       nc     = (numcol-1) / s_nb + 1
       nr     = (numrow-1) / s_nb + 1 
       
       ! Get first diag block in node
       dblk = fkeep%nodes(node)%blk_sa

       ! Loop over block columns
       do jj = 1, nc
          
          !
          ! Forward solve with block on diagoanl
          !
          call spllt_tic("submit fwd block", 2, task_manager%workerID, timer)
         !call spllt_solve_fwd_block_task(dblk, nrhs, rhs_local, rhs, ldr, &
         !  xlocal, fkeep, fwd_block_id, task_manager)
            call task_manager%solve_fwd_block_task(dblk, nrhs, rhs_local, rhs, &
              ldr, xlocal, fkeep)
          call spllt_tac(2, task_manager%workerID, timer)

          do ii = jj+1, nr

            blk = dblk+ii-jj

            !
            ! Forward update with off-diagonal
            !
            call spllt_tic("submit fwd update", 3, task_manager%workerID, timer)
           !call spllt_solve_fwd_update_task(blk, node, nrhs, rhs_local, &
           !  rhs, ldr, xlocal, fkeep, fwd_update_id, task_manager)
            call task_manager%solve_fwd_update_task(blk, node, nrhs, rhs_local,&
              rhs, ldr, xlocal, fkeep)
            call spllt_tac(3, task_manager%workerID, timer)

          end do
          
          ! Update diag block in node          
          dblk = fkeep%bc(dblk)%last_blk + 1
       end do
    end do
    call spllt_tac(4, task_manager%workerID, timer)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(fwd_submit_id, -1)
#endif
#if defined(SOLVE_TASK_LOCKED)
call task_manager%print("fwd end of submitted task", 0)
 !$ call omp_unset_lock(lock)
 !$ call omp_destroy_lock(lock)
#endif
    
    !$omp taskwait
call task_manager%print("fwd end of execution task", 0)

    call spllt_close_timer(task_manager%workerID, timer)

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

    integer                 :: num_node
    ! Node info
    integer                 :: node
    integer                 :: sa, en
    integer                 :: numcol, numrow ! #column/row in node 
    integer                 :: nc, nr         ! #block-column/block-row in node
    integer                 :: s_nb           ! Block size in node
    integer                 :: jj, ii
    integer                 :: dblk
    ! Block info
    integer                 :: blk            ! Block index
    integer                 :: bwd_submit_id
    integer                 :: nworker
    real(wp), pointer       :: xlocal(:,:)    ! update_buffer workspace
    real(wp), pointer       :: rhs_local(:,:) ! update_buffer workspace
    type(spllt_block), pointer :: p_bc(:)
    type(spllt_timer), save :: timer
 !$ integer(kind=omp_lock_kind) :: lock

    call spllt_open_timer(task_manager%nworker, task_manager%workerID, &
      "solve_bwd", timer)

    nworker       = task_manager%nworker

#if defined(SPLLT_OMP_TRACE)
    bwd_submit_id = task_manager%trace_ids(trace_bwd_submit_pos)
    call trace_event_start(bwd_submit_id, -1)
#endif

    xlocal(1 : fkeep%maxmn * nrhs, 0 : nworker - 1) => workspace(1 : &
      fkeep%maxmn * nrhs * nworker)
    rhs_local(1 : ldr * nrhs, 0 : nworker - 1) => workspace(fkeep%maxmn * nrhs &
      * nworker + 1 : (fkeep%maxmn + ldr) * nrhs * nworker)

    num_node = fkeep%info%num_nodes

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
    
    call spllt_tic("Submit tasks", 4, task_manager%workerID, timer)
    do node = num_node, 1, -1

      ! Get node info
      s_nb   = fkeep%nodes(node)%nb
      sa     = fkeep%nodes(node)%sa
      en     = fkeep%nodes(node)%en
      numcol = en - sa + 1
      numrow = size(fkeep%nodes(node)%index)
      nc     = (numcol-1) / s_nb + 1
      nr     = (numrow-1) / s_nb + 1 

      ! Get first diag block in node
      dblk = fkeep%bc(fkeep%nodes(node)%blk_en)%dblk

      ! Loop over block columns
      do jj = nc, 1, -1

        do ii = nr, jj+1, -1
          
          blk = dblk+ii-jj ! Block index

          call spllt_tic("submit bwd update", 2, task_manager%workerID, timer)
         !call spllt_solve_bwd_update_task(blk, node, nrhs, rhs_local,   &
         !  rhs, ldr, xlocal, fkeep, bwd_update_id, task_manager)
          call task_manager%solve_bwd_update_task(blk, node, nrhs, rhs_local, &
            rhs, ldr, xlocal, fkeep)
          call spllt_tac(2, task_manager%workerID, timer)

           !
           ! Backward update with block on diagoanl
           !

        end do

        !
        ! Backward solve with block on diagoanl
        !
        call spllt_tic("submit bwd block", 3, task_manager%workerID, timer)
       !call spllt_solve_bwd_block_task(dblk, nrhs, rhs_local, rhs, ldr, &
       !  xlocal, fkeep, bwd_block_id, task_manager)
        call task_manager%solve_bwd_block_task(dblk, nrhs, rhs_local, rhs, ldr,&
          xlocal, fkeep)
        call spllt_tac(3, task_manager%workerID, timer)
       
       ! Update diag block in node       
       if (jj .gt. 1) dblk = fkeep%bc(dblk-1)%dblk
      end do
      
    end do
    call spllt_tac(4, task_manager%workerID, timer)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(bwd_submit_id, -1)
#endif
#if defined(SOLVE_TASK_LOCKED)
  call task_manager%print("bwd end of submitted task", 0)
 !$ call omp_unset_lock(lock)
 !$ call omp_destroy_lock(lock)
#endif

    !$omp taskwait
    call task_manager%print("bwd end of execution task", 0)

    call spllt_close_timer(task_manager%workerID, timer)

  end subroutine solve_bwd
  
end module spllt_solve_mod
