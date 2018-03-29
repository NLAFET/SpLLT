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
      scheduler)
    use spllt_data_mod
    use utils_mod
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
    type(spllt_omp_scheduler), intent(inout), optional, target  :: scheduler

    integer                             :: n            ! #rows
    integer                             :: solve_step   ! Selector step
    integer                             :: st           ! stat parameter
    real(wp),                  pointer  :: work(:)
    type(spllt_omp_scheduler), pointer  :: sched

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    n = fkeep%n

    !!!!!!!!!!!!!!!!!!!!!
    ! Check optional parameters
    !
    if(present(job)) then
      solve_step = job
    else
      solve_step = 0
    end if

    if(.not.present(scheduler)) then
      allocate(sched)
      call spllt_omp_init_scheduler(sched, stat=st)
      if(st .ne. 0) then
        print *, "Error in creation of the scheduler"
      end if
    else
      sched => scheduler
    end if
    
    !!!!!!!!!!!!!!!!!!!!!
    ! worker
    !
   !allocate(work(n), stat = st)
   !call spllt_scheduler_alloc(sched, st)
    allocate(work(n + (fkeep%maxmn + n) * sched%nworker), stat = st)
    call spllt_scheduler_alloc(sched, st)

    call spllt_solve_mult_double_worker(fkeep, options, order, 1, x, info, &
      solve_step, work, sched)

    !!!!!!!!!!!!!!!!!!!!!
    ! Desallocation
    !
    deallocate(work)

    if(.not.present(scheduler)) then
      deallocate(sched%task_info)
      deallocate(sched)
    end if
    
  end subroutine spllt_solve_one_double

  



  subroutine spllt_solve_mult_double(fkeep, options, order, nrhs, x, info, &
      job, scheduler)
    use spllt_data_mod
    use utils_mod
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
    type(spllt_omp_scheduler), intent(inout), optional, target  :: scheduler
  
    integer                             :: n            ! #rows
    integer                             :: solve_step   ! Selector step
    integer                             :: st           ! stat parameter
    real(wp),                  pointer  :: work(:)
    type(spllt_omp_scheduler), pointer  :: sched

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    n = fkeep%n

    !!!!!!!!!!!!!!!!!!!!!
    ! Check optional parameters
    !
    if(present(job)) then
      solve_step = job
    else
      solve_step = 0
    end if

    if(.not.present(scheduler)) then
      allocate(sched)
      call spllt_omp_init_scheduler(sched, stat=st)
      if(st .ne. 0) then
        print *, "Error in creation of the scheduler"
      end if
    else
      sched => scheduler
    end if
    
    !!!!!!!!!!!!!!!!!!!!!
    ! worker
    !

   !allocate(work(n, nrhs), stat = st)
   !call spllt_scheduler_alloc(sched, st)
   !allocate(work2((fkeep%maxmn+n)*nrhs, sched%nworker), stat = st)
   !call spllt_scheduler_alloc(sched, st)
    allocate(work(n*nrhs + (fkeep%maxmn + n) * nrhs * sched%nworker), stat = st)
    call spllt_scheduler_alloc(sched, st)
    
    call spllt_solve_mult_double_worker(fkeep, options, order, nrhs, x, &
      info, solve_step, work, sched)

    !!!!!!!!!!!!!!!!!!!!!
    ! Desallocation
    !
    deallocate(work)

    if(.not.present(scheduler)) then
      deallocate(sched%task_info)
      deallocate(sched)
    end if
  end subroutine spllt_solve_mult_double





  subroutine spllt_solve_mult_double_worker(fkeep, options, order, nrhs, x, &
      info, job, workspace, scheduler)
    use spllt_data_mod
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
    real(wp), intent(out), target                    :: workspace(:)
    type(spllt_omp_scheduler), intent(inout), target :: scheduler

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
    work2(1 : (fkeep%maxmn + n) * nrhs * scheduler%nworker) =>      &
      workspace(n * nrhs + 1 : nrhs * (n + (fkeep%maxmn + n) *      &
      scheduler%nworker))

    select case(job)
      case(0)
       !
       ! Reorder x
       !
        do j = 1, n
           work(order(j),:) = x(j,:)
        end do

        ! Forward solve
        call solve_fwd(nrhs, work, n, fkeep, work2, scheduler)

        ! Backward solve
        call solve_bwd(nrhs, work, n, fkeep, work2, scheduler)

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

        call solve_fwd(nrhs, work, n, fkeep, work2, scheduler)

        x = work
      
      case(2)
        work = x

        call solve_bwd(nrhs, work, n, fkeep, work2, scheduler)

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
  subroutine solve_fwd(nrhs, rhs, ldr, fkeep, workspace, scheduler)
    use spllt_data_mod
    use spllt_solve_task_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use omp_lib, ONLY : omp_lock_kind
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    real(wp), target,           intent(out)   :: workspace(:)
    type(spllt_omp_scheduler),  intent(inout) :: scheduler

    integer                 :: num_node
    integer                 :: node
    integer                 :: sa, en
    integer                 :: numcol, numrow ! #column/row in node 
    integer                 :: nc, nr         ! #block-column/block-row in node
    integer                 :: jj, ii
    integer                 :: dblk           ! Diagonal index 
    integer                 :: s_nb           ! Block size in node
    integer                 :: blk            ! Block index
    integer                 :: fwd_update_id, fwd_block_id, fwd_solve_id
    integer                 :: nworker
    real(wp), pointer       :: xlocal(:,:)    ! update_buffer workspace
    real(wp), pointer       :: rhs_local(:,:) ! update_buffer workspace
    type(spllt_block), pointer :: p_bc(:)
    type(spllt_timer), save :: timer
 !$ integer(kind=omp_lock_kind) :: lock

    call spllt_open_timer(scheduler%nworker, scheduler%workerID, "solve_fwd", &
      timer)

    fwd_update_id = 0
    fwd_block_id  = 0
    fwd_solve_id  = 0
    nworker       = scheduler%nworker

#if defined(SPLLT_OMP_TRACE)
    if(associated(scheduler%trace_ids)) then
      fwd_update_id = scheduler%trace_ids(1)
      fwd_block_id  = scheduler%trace_ids(2)
      fwd_solve_id  = scheduler%trace_ids(5)
    else
      call trace_create_event("fwd_update", fwd_update_id)
      call trace_create_event("fwd_block", fwd_block_id)
      call trace_create_event("fwd_submit", fwd_solve_id)
    end if
    call trace_event_start(fwd_solve_id, -1)
#endif

    xlocal(1 : fkeep%maxmn * nrhs, 0 : nworker - 1) => workspace(1 : &
      fkeep%maxmn * nrhs * nworker)
    rhs_local(1 : ldr * nrhs, 0 : nworker - 1) => workspace(fkeep%maxmn * nrhs &
      * nworker + 1 : (fkeep%maxmn + ldr) * nrhs * nworker)

    call spllt_tic("Reset workspace", 1, scheduler%workerID, timer)
   !! initialise rhs_local
   !!$omp parallel
   !rhs_local(1 : ldr * nrhs, omp_get_thread_num()) = zero
   !!$omp end parallel
   !call spllt_tac(1, scheduler%workerID, timer)

    num_node = fkeep%info%num_nodes
    
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

    call spllt_tic("Submit tasks", 4, scheduler%workerID, timer)
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
          call spllt_tic("submit fwd block", 2, scheduler%workerID, timer)
          call spllt_solve_fwd_block_task(dblk, nrhs, rhs_local, rhs, ldr, &
            xlocal, fkeep, fwd_block_id, scheduler)
          call spllt_tac(2, scheduler%workerID, timer)

          do ii = jj+1, nr

            blk = dblk+ii-jj

            !
            ! Forward update with off-diagonal
            !
            call spllt_tic("submit fwd update", 3, scheduler%workerID, timer)
            call spllt_solve_fwd_update_task(blk, node, nrhs, rhs_local, &
              rhs, ldr, xlocal, fkeep, fwd_update_id, scheduler)
            call spllt_tac(3, scheduler%workerID, timer)

          end do
          
          ! Update diag block in node          
          dblk = fkeep%bc(dblk)%last_blk + 1
       end do
    end do
    call spllt_tac(4, scheduler%workerID, timer)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(fwd_solve_id, -1)
#endif
 !$ call omp_unset_lock(lock)
 !$ call omp_destroy_lock(lock)
    
    !$omp taskwait

    call spllt_close_timer(scheduler%workerID, timer)

  end subroutine solve_fwd

  subroutine solve_bwd(nrhs, rhs, ldr, fkeep, workspace, scheduler)
    use spllt_data_mod
    use spllt_solve_task_mod
    use spllt_solve_kernels_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use omp_lib, ONLY : omp_lock_kind
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    real(wp),                   intent(inout) :: rhs(ldr, nrhs)
    real(wp), target,           intent(out)   :: workspace(:)
    type(spllt_omp_scheduler),  intent(inout) :: scheduler

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
    integer                 :: bwd_update_id, bwd_block_id, bwd_solve_id
    integer                 :: nworker
    real(wp), pointer       :: xlocal(:,:)    ! update_buffer workspace
    real(wp), pointer       :: rhs_local(:,:) ! update_buffer workspace
    type(spllt_block), pointer :: p_bc(:)
    type(spllt_timer), save :: timer
 !$ integer(kind=omp_lock_kind) :: lock

    call spllt_open_timer(scheduler%nworker, scheduler%workerID, "solve_bwd", &
      timer)

    bwd_update_id = 0
    bwd_block_id  = 0
    bwd_solve_id  = 0
    nworker       = scheduler%nworker

#if defined(SPLLT_OMP_TRACE)
    if(associated(scheduler%trace_ids)) then
      bwd_update_id = scheduler%trace_ids(3)
      bwd_block_id  = scheduler%trace_ids(4)
      bwd_solve_id  = scheduler%trace_ids(6)
    else
      call trace_create_event("bwd_update", bwd_update_id)
      call trace_create_event("bwd_block", bwd_block_id)
      call trace_create_event("bwd_submit", bwd_solve_id)
    end if
    call trace_event_start(bwd_solve_id, -1)
#endif

    xlocal(1 : fkeep%maxmn * nrhs, 0 : nworker - 1) => workspace(1 : &
      fkeep%maxmn * nrhs * nworker)
    rhs_local(1 : ldr * nrhs, 0 : nworker - 1) => workspace(fkeep%maxmn * nrhs &
      * nworker + 1 : (fkeep%maxmn + ldr) * nrhs * nworker)

!   call spllt_tic("Reset workspace", 1, scheduler%workerID, timer)
!   ! initialise rhs_local
!  !rhs_local = zero
!   call spllt_tac(1, scheduler%workerID, timer)

    num_node = fkeep%info%num_nodes

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
    
    call spllt_tic("Submit tasks", 4, scheduler%workerID, timer)
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

          call spllt_tic("submit bwd update", 2, scheduler%workerID, timer)
          call spllt_solve_bwd_udpate_task(blk, node, nrhs, rhs_local,   &
            rhs, ldr, xlocal, fkeep, bwd_update_id, scheduler)
          call spllt_tac(2, scheduler%workerID, timer)

           !
           ! Backward update with block on diagoanl
           !

        end do

        !
        ! Backward solve with block on diagoanl
        !
        call spllt_tic("submit bwd block", 3, scheduler%workerID, timer)
        call spllt_solve_bwd_block_task(dblk, nrhs, rhs_local, rhs, ldr, &
          xlocal, fkeep, bwd_block_id, scheduler)
        call spllt_tac(3, scheduler%workerID, timer)
       
       ! Update diag block in node       
       if (jj .gt. 1) dblk = fkeep%bc(dblk-1)%dblk
      end do
      
    end do
    call spllt_tac(4, scheduler%workerID, timer)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(bwd_solve_id, -1)
#endif
 !$ call omp_unset_lock(lock)
 !$ call omp_destroy_lock(lock)

    !$omp taskwait

    call spllt_close_timer(scheduler%workerID, timer)

  end subroutine solve_bwd
  
end module spllt_solve_mod
