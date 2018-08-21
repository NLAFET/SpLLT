program spllt_test
  use spral_rutherford_boeing
  use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_PSDEF
  use spllt_mod
  use spllt_solve_mod
  use utils_mod
  use spllt_solve_dep_mod
  use trace_mod
  use timer_mod
  use task_manager_omp_mod
  use task_manager_seq_mod
  use ISO_Fortran_env, only: stdout => OUTPUT_UNIT, &
    compiler_version, compiler_options
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

  ! runtime
  type(task_manager_omp_t)            :: task_manager
 !type(task_manager_seq_t)            :: task_manager
  type(spllt_timer_t),        save    :: timer

  call date_and_time(DATE=date, TIME=time)

  write( stdout, '(/8a/)') ' This output was compiled using ', &
    compiler_version(), ' with the options ', compiler_options(), &
    ' and executed ', date, ' at ', time

  call spllt_parse_args(options, matfile, nrhs)

  ! If input matrix is not specified then use matrix.rb file. 
  if (matfile .eq. '') matfile = 'matrix.rb'
  
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

  select case(options%fmt)
    case('csc')
       ! Rutherford boeing format
       rb_options%values = 3 ! Force diagonal dominance
       call rb_read(matfile, m, n, ptr, row, val, rb_options, rb_flag)
       if(rb_flag.ne.0) then
          print *, "Rutherford-Boeing read failed with error ", rb_flag
          stop
       endif
    case('coo')
       ! Matrix Market format
       call mm_read(matfile, m, n, nnz, indx, jndx, val_in, mm_flag)
       if(mm_flag.ne.0) then
          print *, "Matrix Market read failed with error ", mm_flag
          stop
       endif

       ! convert to csc format
       call coo_to_csc_double(m, n, nnz, indx, jndx, val_in, & 
            row, val, ptr, mm_flag)
       if(mm_flag.ne.0) then
          print *, "COO to CSC convertion failed with error ", mm_flag
          stop
       endif

       deallocate(indx, jndx, val_in)

    case default
       print *, "Matrix format not suported"
       stop
  end select
  write(*, "(a)") "ok"

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

  write (header, '(a, a, a, a, i4, a, a, i4, a, a, i4, a, a, i4, a, a, i4)')  &
    '# matrix   = ', trim(matfile), ACHAR(10),                                &
    '# nrhs     = ', nrhs, ACHAR(10),                                         &
    '# nb       = ', options%nb, ACHAR(10),                                   &
    '# snb      = ', options%snb, ACHAR(10),                                  &
    '# nworker  = ', task_manager%nworker, ACHAR(10),                         &
    '# ncpu     = ', options%ncpu, ACHAR(10),                                 &
    '# nemin    = ', options%nemin

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
 !call spllt_print_atree(akeep, fkeep, options)

 !allocate(porder(n))
 !do i = 1, n
 !  porder(order(i)) = i
 !end do
 !fkeep%p_porder => porder
 !print *, "Order : ", order

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

 !call get_solve_blocks(fkeep, options%nb, nrhs, rhs, workspace, sbc)
 !call get_solve_blocks(fkeep, options%nb, nrhs, worksize, sbc)
  call get_solve_blocks(fkeep, options%nb, nrhs, worksize, fkeep%sbc)
 !fkeep%sbc => sbc

  call spllt_compute_solve_dep(fkeep, stat = st)

  call spllt_tac(3, task_manager%workerID, timer)

 !call spllt_compute_rhs_block(fkeep, st)

  call spllt_tac(2, task_manager%workerID, timer)

 !call spllt_solve_workspace_size(fkeep, task_manager%nworker, nrhs, worksize)
 !worksize = fkeep%n
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

  call spllt_tic("Reset workspace", 6, task_manager%workerID, timer)
  y         = zero
  workspace = zero
  call spllt_tac(6, task_manager%workerID, timer)

  call sblock_assoc_mem(fkeep, options%nb, nrhs, y, workspace, fkeep%sbc)
! fkeep%p_y => y

  ! Init the computed solution with the rhs that is further updated by
  ! the subroutine
  sol_computed = rhs

  !!!!!!!!!!!!!!!!!!!!
  ! Forward substitution
  !
! call task_manager%nflop_reset()
! call spllt_tic("Forward", 4, task_manager%workerID, timer)

! call spllt_solve(fkeep, options, order, nrhs, sol_computed, info, job=7, &
!   workspace=workspace, task_manager=task_manager)
! call spllt_wait()

! call spllt_tac(4, task_manager%workerID, timer)


  !!!!!!!!!!!!!!!!!!!!
  ! Backward substitution
  !
! call task_manager%nflop_reset()
! call spllt_tic("Backward", 5, task_manager%workerID, timer)

! call spllt_solve(fkeep, options, order, nrhs, sol_computed, info, job=8, &
!   workspace=workspace, task_manager=task_manager)
! call spllt_wait()

! call spllt_tac(5, task_manager%workerID, timer)

  !!!!!!!!!!!!!!!!!!!!
  ! Solve
  !
! sol_computed = rhs
! call task_manager%nflop_reset()
  call spllt_tic("Solving", 7, task_manager%workerID, timer)

 !print *, "Order : ", order
  call spllt_solve(fkeep, options, order, nrhs, sol_computed, info, &
   !job=merge(3,0, options%ileave_solve), &
   !job=merge(3,6, options%ileave_solve), &
    job=6,                                                          &
    workspace=workspace, task_manager=task_manager)
  call spllt_wait()

  call spllt_tac(7, task_manager%workerID, timer)


  !!!!!!!!!!!!!!!!!!!!
  ! STATISTICS
  !
  !Compute the residual for each rhs
  call check_backward_error(n, ptr, row, val, nrhs, sol_computed, rhs)

  call spllt_deallocate_akeep(akeep, st)
  deallocate(workspace, stat=st)
  call spllt_deallocate_fkeep(fkeep, st)

  !$omp end single
  !$omp end parallel

  call spllt_close_timer(task_manager%workerID, timer)
  call spllt_print_timers()

#if defined(SPLLT_OMP_TRACE)
  ppos = scan(trim(matfile),"/", BACK= .true.)
  if ( ppos > 0 ) matfile = matfile(ppos+1:)
  call trace_log_dump_paje('trace_spllt_'//trim(matfile)//'.out_'&
      &//date//'-'//time)
#endif

  call task_manager%print()

  deallocate(order, rhs, sol, sol_computed, ptr, row, val)

end program spllt_test
