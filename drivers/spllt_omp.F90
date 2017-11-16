program spllt_omp
  use spral_rutherford_boeing
  use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_PSDEF
  use spllt_mod
  implicit none

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
  integer :: nrhs = 1 ! Numebr of right-hand side
  double precision, dimension(:, :), allocatable :: rhs, soln 
  double precision, dimension(:), allocatable :: res
  integer i, j, k, r
  integer :: flag, more
  type(spllt_akeep) :: akeep ! Symbolic factorization data
  type(spllt_fkeep), target  :: fkeep ! Factorization data

  ! timing
  integer :: start_t, stop_t, rate_t

  ! stats
  real :: smanal, smfact, smaflop, smafact
  integer, dimension(:), allocatable :: order ! Matrix permutation array

  call spllt_parse_args(options, matfile, nrhs)

  ! If input matrix is not specified then use matrix.rb file. 
  if (matfile .eq. '') matfile = 'matrix.rb'
  
  ! Print user-supplied options
  print "(a, a)", "Matrix file = ", matfile
  print "(a, a)", "Matrix format = ", options%fmt
  print "(a, i4)", "Number of CPUs = ", options%ncpu
  print "(a, i4)", 'Block size = ', options%nb
  print "(a, i4)", 'Supernode amalgamation nemin = ', options%nemin

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
  allocate(rhs(n, nrhs), soln(n, nrhs))
  rhs = 0
  do r = 1, nrhs
     do i = 1, n
        do j = ptr(i), ptr(i+1)-1
           k = row(j)
           rhs(k, r) = rhs(k, r) + val(j)
           if(i.eq.k) cycle
           rhs(i, r) = rhs(i, r) + val(j)
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
  
  allocate(order(n))

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
  smaflop = real(info%ssids_inform%num_flops)
  smafact = real(info%ssids_inform%num_factor)

  deallocate(order)
  
end program spllt_omp
