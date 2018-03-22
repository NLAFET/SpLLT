program spllt_omp
  use spral_rutherford_boeing
  use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_PSDEF
  use spllt_mod
  use spllt_solve_mod
  use utils_mod
  use spllt_solve_dep_mod
  use trace_mod
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
  integer                             :: nrhs_max = 1 ! # of right-hand side
  double precision, allocatable       :: rhs(:,:)
  double precision, allocatable       :: sol(:,:) 
  double precision, allocatable       :: sol_computed(:,:)
  double precision, allocatable       :: res(:,:)
  integer                             :: i, j, k, r, th
  integer                             :: flag, more
  logical                             :: bwd_error_ok
  type(spllt_akeep)                   :: akeep ! Symbolic factorization data
  type(spllt_fkeep), target           :: fkeep ! Factorization data

  ! timing
  integer                             :: start_t, stop_t, rate_t
  double precision, allocatable       :: fwd_timer(:,:)
  double precision, allocatable       :: bwd_timer(:,:)
  double precision, allocatable       :: analyse_timer(:)
  double precision, allocatable       :: facto_timer(:)

  ! stats
  integer,           allocatable      :: order(:)     ! Matrix permutation array
  real(wp),          allocatable      :: workspace(:) ! Workspace
  double precision,  allocatable      :: normRes(:), normRHS(:)
  double precision,  allocatable      :: errNorm(:), solNorm(:)
  integer                             :: nnrhs, nnb, nb_i
  integer,           allocatable      :: nrhs_list(:), nb_list(:)
  character(len=10), allocatable      :: trace_names(:)
  character(len=1024)                 :: header
  character(len=10)                   :: time
  character(len=8)                    :: date
  double precision, allocatable       :: fwd_flops(:,:), bwd_flops(:,:)

  ! runtime
  type(spllt_omp_scheduler)           :: scheduler

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
  print "(a, i4)",  "Supernode amalgamation nemin = ", options%nemin

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

  allocate(trace_names(4))
  trace_names = [character(len=10) :: "fwd_update", "fwd_block", &
    "bwd_update", "bwd_block" ]

  call spllt_omp_init_scheduler(scheduler, trace_names, st)
#else
  call spllt_omp_init_scheduler(scheduler, stat=st)
#endif

  call compute_range(options%nb_min, options%nb_max,      &
    options%nb_linear_comp, nb_list)
  call compute_range(options%nrhs_min, options%nrhs_max,  &
    options%nrhs_linear_comp, nrhs_list)

  nnrhs     = size(nrhs_list)
  nnb       = size(nb_list)
  
  write (header, '(a, a, a, a, i4, a, i4, a, a, a, i4, a, i4, a, a, a, i4,    &
    &a, a, i4, a, a, i4)')                                                    &
    '# matrix   =  ', trim(matfile), ACHAR(10),                               &
    '# nrhs     = [', nrhs_list(1), ',', nrhs_list(nnrhs), ']', ACHAR(10),    &
    '# nb       = [', nb_list(1), ',', nb_list(nnb), ']', ACHAR(10),          &
    '# nworker  = ', scheduler%nworker, ACHAR(10),                            &
    '# ncpu     = ', options%ncpu, ACHAR(10),                                 &
    '# nemin    = ', options%nemin

  nrhs_max  = nrhs_list(nnrhs)

  !$omp end single
  !$omp end parallel
  


  !!!!!!!!!!!!!!!!!!!!
  ! Create RHS
  !
  ! Make up a rhs associated with the solution x = 1.0
  allocate(sol(n, nrhs_max), stat=st)
  allocate(sol_computed(n, nrhs_max), stat=st)
  allocate(rhs(n, nrhs_max), stat=st)

  rhs = 0
  ! Set up solution to constant vector equals to its index in sol array
  do r = 1, nrhs_max
    sol(:,r) = r
  end do

  do r = 1, nrhs_max
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
  allocate(analyse_timer(nnb), facto_timer(nnb))
  allocate(res(n, nrhs_max), normRes(nrhs_max), normRHS(nrhs_max), &
    errNorm(nrhs_max), solNorm(nrhs_max))
  allocate(fwd_timer(nnrhs, nnb), bwd_timer(nnrhs, nnb))
  allocate(fwd_flops(nnrhs, nnb), bwd_flops(nnrhs, nnb))


  !$omp parallel
  !$omp single
  do nb_i=1, nnb

    options%nb = nb_list(nb_i)
    !!!!!!!!!!!!!!!!!!!!
    ! Analyse SpLLT
    !
    call system_clock(start_t, rate_t)
    call spllt_analyse(akeep, fkeep, options, n, ptr, row, info, order)
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "error detected during analysis"
       stop
    endif
    call system_clock(stop_t)
    analyse_timer(nb_i) = (stop_t - start_t)/real(rate_t)
   !print "(a,es10.2)", "Predict nfact = ", real(info%ssids_inform%num_factor)
   !print "(a,es10.2)", "Predict nflop = ", real(info%ssids_inform%num_flops)

    !!!!!!!!!!!!!!!!!!!!
    ! Numerical Factorization
    !
   !!$omp parallel
   !!$omp single
    call system_clock(start_t, rate_t)
    call spllt_factor(akeep, fkeep, options, val, info)
    call spllt_wait()
    call system_clock(stop_t)
   !!$omp end single
   !!$omp end parallel
    facto_timer(nb_i) = (stop_t - start_t)/real(rate_t)

    allocate(workspace(n * nrhs_max + (fkeep%maxmn + n) * &
      nrhs_max * scheduler%nworker), stat = st)
    call spllt_scheduler_alloc(scheduler, st)
    

    !!!!!!!!!!!!!!!!!!!!
    ! Compute dependencies of each blk
    !
    call spllt_compute_solve_dep(fkeep)

   !!$omp parallel
   !!$omp single
    do j=1, nnrhs

      nrhs = nrhs_list(j)
      ! Init the computed solution with the rhs that is further updated by
      ! the subroutine
      sol_computed = rhs

      !!!!!!!!!!!!!!!!!!!!
      ! Forward substitution
      !
      call system_clock(start_t, rate_t)
      call spllt_solve(fkeep, options, order, nrhs, sol_computed, info, job=1, &
        workspace=workspace, scheduler=scheduler)
      call system_clock(stop_t)
      fwd_timer(j,nb_i) = (stop_t - start_t)/real(rate_t)
#if defined(SPLLT_PROFILING_FLOP)
      fwd_flops(j, nb_i) =  0.0
      do th = lbound(scheduler%task_info, 1), ubound(scheduler%task_info, 1)
        fwd_flops(j,nb_i) = fwd_flops(j,nb_i) + scheduler%task_info(th)%nflop
        scheduler%task_info(th)%nflop = 0.0
      end do
#endif

      !!!!!!!!!!!!!!!!!!!!
      ! Backward substitution
      !
      call system_clock(start_t, rate_t)
      call spllt_solve(fkeep, options, order, nrhs, sol_computed, info, job=2, &
        workspace=workspace, scheduler=scheduler)
      call system_clock(stop_t)
      bwd_timer(j,nb_i) = (stop_t - start_t)/real(rate_t)
#if defined(SPLLT_PROFILING_FLOP)
      bwd_flops(j, nb_i) =  0.0
      do th = lbound(scheduler%task_info, 1), ubound(scheduler%task_info, 1)
        bwd_flops(j,nb_i) = bwd_flops(j,nb_i) + scheduler%task_info(th)%nflop
        scheduler%task_info(th)%nflop = 0.0
      end do
#endif


#if defined(SPLLT_DRIVER_CHECK_ERROR)
      !!!!!!!!!!!!!!!!!!!!
      ! STATISTICS
      !
      !Compute the residual for each rhs
      call compute_residual(n, ptr, row, val, nrhs, &
        sol_computed, rhs, res)

      call vector_norm_2(n, res, normRes)
      call vector_norm_2(n, rhs, normRHS)
      call vector_norm_2(n, abs(sol_computed - sol), errNorm)
      call vector_norm_2(n, sol, solNorm)

      bwd_error_ok = .true.
      do i = 1, nrhs
        if(normRes(i) / normRHS(i) .gt. 1e-14) then
          write(0, "(a, i4, a, i4)") "Wrong Bwd error for ", i, "/", nrhs
          bwd_error_ok = .false.
        end if
      end do
      if(bwd_error_ok) then
        write(0, "(a)") "Backward error... ok"
      end if
#endif
    end do
   !!$omp end single
   !!$omp end parallel

    call spllt_deallocate_akeep(akeep, st)
    deallocate(workspace, stat=st)
    call spllt_deallocate_fkeep(fkeep, st)
  end do
  !$omp end single
  !$omp end parallel


  do i = 1, scheduler%nworker
    call print_omp_task_stat("SCHEDULER STAT", i, scheduler%task_info(i))
  end do

 !print *, "==========="
 !print *, "SOLVE TIMER"
 !do nb_i = 1, nnb
 !  do j = 1, nnrhs
 !    print *, "fwd ", fwd_timer(j,nb_i), " bwd ", bwd_timer(j,nb_i)
 !  end do
 !end do

  call timer_log_dump("fwd step"//ACHAR(10)//header, fwd_timer, &
    'fwd_time_'//trim(matfile)//'.out_'//date//'-'//time)
  call timer_log_dump("bwd step"//ACHAR(10)//header, bwd_timer, &
    'bwd_time_'//trim(matfile)//'.out_'//date//'-'//time)
  call timer_log_dump("analyse step"//ACHAR(10)//header, analyse_timer, &
    'analyse_time_'//trim(matfile)//'.out_'//date//'-'//time)
  call timer_log_dump("facto step"//ACHAR(10)//header, facto_timer, &
    'facto_time_'//trim(matfile)//'.out_'//date//'-'//time)

#if defined(SPLLT_PROFILING_FLOP)
  call flop_log_dump("fwd step"//ACHAR(10)//header, fwd_flops, &
    'fwd_time_'//trim(matfile)//'.out_'//date//'-'//time)
  call flop_log_dump("bwd step"//ACHAR(10)//header, bwd_flops, &
    'bwd_time_'//trim(matfile)//'.out_'//date//'-'//time)
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_log_dump_paje('trace_bench_full_'//trim(matfile)//'.out_'&
    &//date//'-'//time)
#endif

  deallocate(order, rhs, sol, sol_computed, ptr, row, val)
  deallocate(fwd_timer, bwd_timer, nrhs_list, facto_timer, analyse_timer)
  if(allocated(trace_names)) then
    deallocate(trace_names)
  end if

end program spllt_omp