!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
program test_solve_phasis
  use spral_rutherford_boeing
  use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_PSDEF
  use spllt_mod
  use spllt_solve_mod
  use utils_mod
  use spllt_solve_dep_mod
  use trace_mod
  use timer_mod
  use task_manager_omp_mod
  implicit none

  type(spllt_options)                 :: options ! User-supplied options 
  type(spllt_inform)                  :: info
  integer                             :: st

  ! matrix reader options (Rutherford Boeing)
  type(rb_read_options)               :: rb_options
  integer                             :: rb_flag
  ! Matrix description (Rutherford Boeing)
  character(len=200)                  :: matfile = ''
  integer                             :: m, n
  integer, dimension(:), allocatable  :: ptr, row
  real(wp), dimension(:), allocatable :: val
  ! Matrix reader options (Matrix Market)
  integer                             :: mm_flag
  integer                             :: nnz
  ! Matrix description (Rutherford Boeing)
  integer, dimension(:), allocatable  :: indx, jndx
  real(wp), dimension(:), allocatable :: val_in

  ! right-hand side and solution
  integer                             :: nrhs     = 1 ! # of right-hand side
  double precision, allocatable       :: rhs(:,:)
  double precision, allocatable       :: sol(:,:) 
  double precision, allocatable       :: sol_computed(:,:)
  integer                             :: i, j, k, r
  integer                             :: flag, more
  type(spllt_akeep)                   :: akeep ! Symbolic factorization data
  type(spllt_fkeep), target           :: fkeep ! Factorization data

  ! Solving
! type(spllt_sblock_t), target, allocatable :: sbc(:)

  ! stats
  integer, target,   allocatable      :: order(:)     ! Matrix permutation array
! integer, target,   allocatable      :: porder(:)    ! Matrix permutation array
  real(wp),          allocatable      :: workspace(:) ! Workspace
  real(wp), target,  allocatable      :: y(:) ! Workspace
  integer(long)                       :: worksize
  character(len=1024)                 :: header
  character(len=10)                   :: time
  character(len=8)                    :: date
  integer                             :: trace_chkerr_id
  integer                             :: ppos
  integer                             :: check

  ! runtime
  type(task_manager_omp_t)            :: task_manager
 !type(task_manager_seq_t)            :: task_manager
  type(spllt_timer_t),        save    :: timer

  call spllt_parse_args(options, matfile, nrhs)

  ! If input matrix is not specified then use matrix.rb file. 
  if (matfile .eq. '') matfile = 'matrix.rb'
  check = 0

  ! Print user-supplied options
  print "(a, a)",   "Matrix file                  = ", matfile
  print "(a, a)",   "Matrix format                = ", options%fmt
  print "(a, i4)",  "Number of CPUs               = ", options%ncpu
  print "(a, i4)",  "Block size                   = ", options%nb
  print "(a, i4)",  "Block size for solve         = ", options%snb
  print "(a, i4)",  "Supernode amalgamation nemin = ", options%nemin

  ! If not provided by the user, set the block size for the solve to 
  ! the block size of the Factorization
  if(options%snb .eq. -1) options%snb = options%nb

  ! Read in a matrix
  write (*, "(a)") "Reading..."

  ! Rutherford boeing format
  rb_options%values = 3 ! Force diagonal dominance
  call rb_read(matfile, m, n, ptr, row, val, rb_options, rb_flag)
  if(rb_flag.ne.0) then
    print *, "Rutherford-Boeing read failed with error ", rb_flag
    stop
  endif

  !!!!!!!!!!!!!!!!!!!!!
  ! Init env of test
  !

  !$omp parallel 
  !$omp single

#if defined(SPLLT_OMP_TRACE)
  call trace_init(omp_get_num_threads())

  call task_manager%init(stat=st)
  trace_chkerr_id = task_manager%trace_ids(trace_chk_err_pos)
#else
  trace_chkerr_id = 0 ! be sure the trace module fails
  call task_manager%init(stat=st)
#endif

  write (header, '(a, a, a, a, i4, a,&
    &a, i4, a, a, i4, a,&
    &a, i4, a, a, i4, a,&
    &a, i4)')  &
    '# matrix   = ', trim(matfile), ACHAR(10),                                &
    '# nrhs     = ', nrhs, ACHAR(10),                                         &
    '# nb       = ', options%nb, ACHAR(10),                                   &
    '# snb      = ', options%snb, ACHAR(10),                                  &
    '# nworker  = ', task_manager%nworker, ACHAR(10),                         &
    '# ncpu     = ', options%ncpu, ACHAR(10),                                 &
    '# nemin    = ', options%nemin

  print *, trim(header)

  !$omp end single
  !$omp end parallel
  
  call spllt_init_timer(st)

  !!!!!!!!!!!!!!!!!!!!
  ! Create RHS
  !
  ! Make up a rhs associated with the solution x = 1.0
  allocate(sol(n, nrhs), stat=st)
  allocate(sol_computed(n, nrhs), stat=st)
  allocate(rhs(n, nrhs), stat=st)

  rhs = 0
  ! Set up solution to constant vector equals to its index in sol array
  do r = 1, nrhs
    sol(:,r) = r
  end do

  do r = 1, nrhs
     do i = 1, n
        do j = ptr(i), ptr(i+1)-1
           k = row(j)
           rhs(k, r) = rhs(k, r) + val(j) * sol(k, r)
           if(i.eq.k) cycle
           rhs(i, r) = rhs(i, r) + val(j) * sol(k, r)
        end do
     end do
  end do

  ! check matrix format is correct
  call cscl_verify(6, SPRAL_MATRIX_REAL_SYM_PSDEF, n, n, &
       ptr, row, flag, more)
  if(flag.ne.0) then
     print *, "CSCL_VERIFY failed: ", flag, more
     stop
  endif
  
  allocate(order(n), stat = st)

  !$omp parallel
  !$omp single
  
  call task_manager%reset()
  call spllt_open_timer(task_manager%workerID, "MAIN", timer)

  !!!!!!!!!!!!!!!!!!!!
  ! Analyse SpLLT
  !
  call spllt_tic("Symbolic", 1, task_manager%workerID, timer)

  call spllt_analyse(akeep, fkeep, options, n, ptr, row, info, order)

  call spllt_tac(1, task_manager%workerID, timer)
  if(info%flag .lt. spllt_success) then
     write(*, "(a)") "error detected during analysis"
     stop
  endif

  !!!!!!!!!!!!!!!!!!!!
  ! Numerical Factorization
  !
  call spllt_tic("Factorize", 2, task_manager%workerID, timer)

  call spllt_factor(akeep, fkeep, options, val, info)
  call spllt_wait()

  !!!!!!!!!!!!!!!!!!!!
  ! Compute dependencies of each blk
  !
  call spllt_tic("Compute Dep", 3, task_manager%workerID, timer)

  call spllt_create_subtree(akeep, fkeep)

  call get_solve_blocks(fkeep, options%nb, nrhs, worksize, fkeep%sbc)

  call spllt_compute_solve_dep(fkeep, stat = st)

  call spllt_tac(3, task_manager%workerID, timer)

  call spllt_tac(2, task_manager%workerID, timer)

  allocate(workspace(worksize), stat = st)
  print '(a, es10.2)', "Allocation of a workspace of size ", real(worksize)
  if(st.ne.0) then
    write(0,fmt='(a, es10.2)') "Can not allocate the workspace of size ", &
      real(worksize)
    stop
  end if
  call task_manager%incr_alloc(st)

  allocate(y(fkeep%n * nrhs), stat=st)
  print '(a, es10.2)', "Allocation of y vector of size ", real(fkeep%n * nrhs)
  if(st.ne.0) then
    write(0,fmt='(a, es10.2)') "Can not allocate y vector of size ", &
      real(fkeep%n * nrhs)
    stop
  end if
  call task_manager%incr_alloc(st)

  call sblock_assoc_mem(fkeep, options%nb, nrhs, y, workspace, fkeep%sbc)

  !$omp end single
  !$omp end parallel

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !     TEST

  !$omp parallel num_threads(task_manager%nworker)
  !$omp single

  !======================================================
  !!!!!!!!!!!!!!!!!!!!!
  !   Synchronous call
  sol_computed = rhs

  call spllt_solve(fkeep, options, nrhs, sol_computed, 0, info)
  call check_backward_error(n, ptr, row, val, nrhs, sol_computed, rhs)
  check = check + 1

  sol_computed = rhs
  call spllt_solve(fkeep, options, nrhs, sol_computed, 1, info)
  call spllt_solve(fkeep, options, nrhs, sol_computed, 2, info)
  call check_backward_error(n, ptr, row, val, nrhs, sol_computed, rhs)
  check = check + 1

  do i = 1, nrhs
    sol_computed(:,i) = rhs(:,i)
    call spllt_solve(fkeep, options, sol_computed(:,i), 0, info)
    call check_backward_error(n, ptr, row, val, &
      sol_computed(:,i), rhs(:,i))
    check = check + 1

    sol_computed(:,i) = rhs(:,i)
    call spllt_solve(fkeep, options, sol_computed(:,i), 1, info)
    call spllt_solve(fkeep, options, sol_computed(:,i), 2, info)
    call check_backward_error(n, ptr, row, val, &
      sol_computed(:,i), rhs(:,i))
    check = check + 1
  end do

  !======================================================
  !!!!!!!!!!!!!!!!!!!!!
  !   Asynchronous call
  !
  sol_computed = rhs

  !!!!!!!!!!!!!!!!!!!!
  ! Forward substitution, then Backward substitution
  !
  call spllt_solve(fkeep, options, nrhs, sol_computed, 1, &
    task_manager, info)
  call spllt_wait()

  call spllt_solve(fkeep, options, nrhs, sol_computed, 2, &
    task_manager, info)
  call spllt_wait()

  !Compute the residual for each rhs
  call check_backward_error(n, ptr, row, val, nrhs, sol_computed, rhs)
  check = check + 1

  !!!!!!!!!!!!!!!!!!!!
  ! Full Solve
  !
  sol_computed = rhs

  call spllt_solve(fkeep, options, nrhs, sol_computed, &
    0, task_manager, info)
  call spllt_wait()

  !Compute the residual for each rhs
  call check_backward_error(n, ptr, row, val, nrhs, sol_computed, rhs)
  check = check + 1

  !!!!!!!!!!!!!!!!!!!!
  ! fwd & bwd 
  !

  sol_computed = rhs
  call spllt_solve(fkeep, options, nrhs, sol_computed, &
    1, task_manager, info)
  call spllt_solve(fkeep, options, nrhs, sol_computed, &
    2, task_manager, info)
  call spllt_wait()

  !Compute the residual for each rhs
  call check_backward_error(n, ptr, row, val, nrhs, sol_computed, rhs)
  check = check + 1

  !$omp end single
  !$omp end parallel

  write (0, '(a, i3)') "[INFO] NTEST : ", check

  call spllt_deallocate_akeep(akeep, st)
  deallocate(workspace, stat=st)
  deallocate(y, stat=st)
  call spllt_deallocate_fkeep(fkeep, st)

  call spllt_close_timer(task_manager%workerID, timer)
  call spllt_print_timers()

#if defined(SPLLT_OMP_TRACE)
  ppos = scan(trim(matfile),"/", BACK= .true.)
  if ( ppos > 0 ) matfile = matfile(ppos+1:)
  call trace_log_dump_paje('trace_spllt_'//trim(matfile)//'.out_'&
      &//date//'-'//time)
#endif

  call task_manager%print()

  deallocate(sol, sol_computed, rhs, order, ptr, row, val)
end program test_solve_phasis
