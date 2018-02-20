
program spllt_omp
  use spral_rutherford_boeing
  use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_PSDEF
  use spllt_mod
  use spllt_solve_mod
! use trace_mod
  implicit none

! integer, external :: omp_get_thread_num, omp_get_num_threads

  type(spllt_options) :: options ! User-supplied options 
  type(spllt_inform) :: info

  ! matrix reader options (Rutherford Boeing)
  type(rb_read_options) :: rb_options
  integer :: rb_flag
  ! Matrix description (Rutherford Boeing)
  character(len=200) :: matfile = ''
  integer :: m, n
  integer, dimension(:), allocatable :: ptr, row
  real(wp), dimension(:), allocatable :: val
  ! Matrix reader options (Matrix Market)
  integer :: mm_flag
  integer :: nnz
  ! Matrix description (Rutherford Boeing)
  integer, dimension(:), allocatable :: indx, jndx
  real(wp), dimension(:), allocatable :: val_in

  ! right-hand side and solution
  integer :: nrhs = 1 ! Number of right-hand side
  double precision, dimension(:, :), allocatable :: rhs, sol 
  double precision, dimension(:, :), allocatable :: sol_computed
  double precision, dimension(:, :), allocatable :: res
  integer :: i, j, k, r
  integer :: flag, more
  type(spllt_akeep) :: akeep ! Symbolic factorization data
  type(spllt_fkeep), target  :: fkeep ! Factorization data

  ! timing
  integer :: start_t, stop_t, rate_t

  ! stats
! real                                :: smanal, smfact, smaflop, smafact
  integer, dimension(:), allocatable    :: order ! Matrix permutation array
  real(wp), dimension(:,:), allocatable :: work  ! Workspace
  double precision                      :: normInfA

! integer                             :: nthreads, threadId

! !$omp parallel 
! !$omp single
! nthreads = omp_get_num_threads()
! threadId = omp_get_thread_num()
! print *, "Total number of threads ", nthreads
! !$omp end single
! !$omp end parallel

! call trace_init(nthreads)

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
  end select
  write(*, "(a)") "ok"

  ! Make up a rhs associated with the solution x = 1.0
  allocate(sol(n, nrhs), sol_computed(n, nrhs))
  allocate(rhs(n, nrhs))
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
  
  allocate(order(n), work(n, nrhs))

  ! Analyse SpLLT
  write(*, "(a)") "Analyse..."
  call system_clock(start_t, rate_t)
  call spllt_analyse(akeep, fkeep, options, n, ptr, row, info, order)
  if(info%flag .lt. spllt_success) then
     write(*, "(a)") "error detected during analysis"
     stop
  endif
  call system_clock(stop_t)
  write(*, "(a)") "ok"
  print *, "Analyse took ", (stop_t - start_t)/real(rate_t)
  print "(a,es10.2)", "Predict nfact = ", real(info%ssids_inform%num_factor)
  print "(a,es10.2)", "Predict nflop = ", real(info%ssids_inform%num_flops)
! smaflop = real(info%ssids_inform%num_flops)
! smafact = real(info%ssids_inform%num_factor)
  ! Print elimination tree
! call spllt_print_atree(akeep, fkeep, options)

  !Numerical Factorization
  !$omp parallel
  !$omp single
  write(*, "(a)") "Numerical factorization..."
  call system_clock(start_t, rate_t)
  call spllt_factor(akeep, fkeep, options, val, info)
  call spllt_wait()
  call system_clock(stop_t)
  !$omp end single
  !$omp end parallel
  write(*, "(a)") "ok"
  print *, "Numerical factorization took ", (stop_t - start_t)/real(rate_t)

  ! Init the computed solution with the rhs that is further updated by
  ! the subroutine
  sol_computed = rhs

  !Forward substitution
  !$omp parallel
  !$omp single
  print '(a)', "Forward substitution..."
  call system_clock(start_t, rate_t)
  call spllt_solve(fkeep, options, order, nrhs, sol_computed, info, job=1, work=work)
  call system_clock(stop_t)
  print '(a)', "ok"
  print *, "Forward substitution took ", (stop_t - start_t)/real(rate_t)
  !$omp end single
  !$omp end parallel
! if (info%flag .ne. SPLLT_SUCCESS) then
!   write (0, '(a)') 'Execution aborted'
!   return
! end if


  !Backward substitution
  !$omp parallel
  !$omp single
  write(*, "(a)") "Backward substitution..."
  call system_clock(start_t, rate_t)
  call spllt_solve(fkeep, options, order, nrhs, sol_computed, info, job=2)
  call system_clock(stop_t)
  write(*, "(a)") "ok"
  print *, "Backward substitution took ", (stop_t - start_t)/real(rate_t)
  !$omp end single
  !$omp end parallel
! if (info%flag .ne. SPLLT_SUCCESS) then
!   write (0, '(a)') 'Execution aborted'
!   return
! end if

  !
  ! STATISTICS
  !

  ! Compute infinite matrix norm
  call matrix_norm_inf(n, ptr, row, val, normInfA)
  print '(a, es10.2)', "Computed norm infinite of A ", normInfA

  !Compute the residual for each rhs
  allocate(res(n, nrhs))
  call compute_residual(n, ptr, row, val, nrhs, &
    sol_computed, rhs, res)

  print '(a, 3Xa, 5Xa)', "rhs_id", "||x_comp - x_sol ||_inf", &
    "||Ax_comp - rhs ||_inf / ||A||_inf"
  do i = 1, nrhs
    print '(i3, 5Xes10.2, 18Xes10.2)', i,           &
      maxval(abs(sol_computed(:, i) - sol(:, i))),  &
      maxval(abs(res(:, i))) / normInfA
  end do
! call trace_log_dump_paje('trace.out')

  deallocate(order, rhs, sol, sol_computed, ptr, row, val, work)

  stop

contains

  !Compute res = b - Ax 
  subroutine compute_residual(n, ptr, row, val, nrhs, x, b, res)
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in) :: nrhs
    real(wp), dimension(n,nrhs), intent(in) :: x
    real(wp), dimension(n,nrhs), intent(in) :: b
    real(wp), dimension(n,nrhs), intent(inout) :: res

    integer :: i, j, k, r
    ! Find the residual
    !allocate(res(n,nrhs))
    res(:,:) = 0
    call compute_Ax(n, ptr, row, val, nrhs, x, res)
    res = b - res
  end subroutine compute_residual
  
  !Compute Ax
  subroutine compute_Ax(n, ptr, row, val, nrhs, x, res)
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in) :: nrhs
    real(wp), dimension(n,nrhs), intent(in) :: x
    real(wp), dimension(n,nrhs), intent(inout) :: res

    integer :: i, j, k, r
    res = 0
    do i = 1, n
      do j = ptr(i), ptr(i+1)-1
        r = row(j)
        do k = 1, nrhs
          res(r, k) = res(r, k) + val(j) * x(i, k)
          if(r .eq. i) cycle
          res(i, k) = res(i, k) + val(j) * x(r, k)
        end do
      end do
    end do
  end subroutine compute_Ax

  subroutine matrix_norm_inf(n, ptr, row, val, norm)
!   use spllt_data_mod
    integer, intent(in)                         :: n
    integer, dimension(n+1), intent(in)         :: ptr
    integer, dimension(ptr(n+1)-1), intent(in)  :: row    ! Unused
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    real(wp), intent(out)                       :: norm

    integer :: i

    norm = 0
    do i = 1, ptr(n+1)-1
    norm = merge(norm, abs(val(i)), norm > abs(val(i)))
    end do
  end subroutine matrix_norm_inf
  
! subroutine permute_darray(n, val, perm, val_perm, trans)
!   integer,                intent(in)      :: n
!   real(wp), dimension(n), intent(in)      :: val
!   integer,  dimension(n), intent(in)      :: perm
!   real(wp), dimension(n), intent(out)     :: val_perm
!   character(len=1), optional, intent(in)  :: trans

!   integer           :: i
!   character(len=1)  :: permute_type

!   if(.not. present(trans)) then
!     permute_type = 'N'
!   else
!     permute_type = trans
!   end if

!   if(permute_type == 'T') then
!     do i = 1, n
!       val_perm(perm(i)) = val(i)
!     end do
!   else
!     do i = 1, n
!       val_perm(i) = val(perm(i))
!     end do
!   end if
! end subroutine permute_darray

! subroutine permute_array(n, val, perm, val_perm)
!   integer,                intent(in)    :: n
!   real(wp), dimension(n), intent(in)    :: val
!   integer,  dimension(n), intent(in)    :: perm
!   real(wp), dimension(n), intent(out)   :: val_perm

!   integer :: i

!   do i = 1, n
!     val_perm(i) = val(perm(i))
!   end do
! end subroutine permute_array

! subroutine print_darray(array_name, n, val)
!   character(len=*)                      :: array_name
!   integer,                intent(in)    :: n
!   real(wp), dimension(n), intent(in)    :: val

!   integer :: i
!   print '(a)', array_name
!   do i = 1, n
!     print '(f20.8)', val(i)
!   end do
! end subroutine print_darray

! subroutine print_iarray(array_name, n, val)
!   character(len=*)                      :: array_name
!   integer,                intent(in)    :: n
!   integer, dimension(n),  intent(in)    :: val

!   integer :: i
!   print '(a)', array_name
!   do i = 1, n
!     print '(i9)', val(i)
!   end do
! end subroutine print_iarray

end program spllt_omp
