program run_prob
    use spral_rutherford_boeing
    use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_PSDEF
    use spral_metis_wrapper, only : metis_order
   ! use hsl_mc56_double
   ! use hsl_mc68_integer
   ! use hsl_mc69_double
   ! use hsl_zd11_double
   use hsl_ma87_double
   use spllt_mod
   use spllt_data_mod
   use spral_ssids
   implicit none

   integer, parameter :: dp = kind(0d0)

   ! type(mc56_control) :: control56
   ! integer :: info56
   ! type(mc68_control) :: order_control
   ! type(mc68_info) :: order_info

   ! type(ZD11_type) :: matrix
   ! matrix reader options (Rutherford Boeing)
   type(rb_read_options) :: rb_options
   integer :: rb_flag
   ! Matrix description (Rutherford Boeing)
   character(len=200) :: matfile = ''
   integer :: m, n
   integer, dimension(:), allocatable :: ptr, row
   real(wp), dimension(:), allocatable :: val
   ! matrix reader options (Matrix Market)
   integer :: mm_flag
   integer :: nnz
   ! Matrix description (Rutherford Boeing)
   integer, dimension(:), allocatable :: indx, jndx
   real(wp), dimension(:), allocatable :: val_in

   type(ma87_info) :: info
   type(ma87_keep) :: keep
   type(ma87_control) :: control
   double precision, dimension(:), allocatable :: work, scaling
   integer, dimension(:), allocatable :: order, invp
   ! integer, dimension(:), allocatable :: perm
   integer, allocatable :: iwork(:)

   ! right-hand side and solution
   integer :: nrhs
   double precision, dimension(:, :), allocatable :: rhs, soln 
   double precision, dimension(:), allocatable :: res

   integer :: i, j, k, r

   integer :: start_t, stop_t, rate_t, tsf_sa, tsf_en
   integer :: flag, more, stat

   integer, parameter :: unit_rhs = 14

   real :: smanal, smfact, smaflop, smafact

   integer, parameter :: nfact = 1
   !integer, parameter :: nfact = 5
   !integer, parameter :: nfact = 100

   integer, parameter :: nslv = 1
   !integer, parameter :: nslv = 10
   !integer, parameter :: nslv = 100
   type(spllt_options) :: options

   ! ssids options 
   type(ssids_options) :: ssids_opt

   ! ssids structures
   type(ssids_inform) :: inform ! stats
   type(ssids_akeep) :: akeep   ! analysis data

   ! Set nrhs
   nrhs = 1

   call spllt_parse_args(options)

   control%nb   = options%nb
   ! control%ncpu = options%ncpu

   if (options%mat .ne. '') then
      matfile = options%mat
   else
      matfile = 'matrix.rb'
   end if

    if (options%nemin .gt. 0) then
       control%nemin = options%nemin
    end if

    write(*,*) '  mat: ', matfile
    write(*,*) '   nb: ', control%nb    
    write(*,*) ' ncpu: ', options%ncpu
    write(*,*) 'nemin: ', control%nemin


!$ call omp_set_num_threads(options%ncpu)

   !control%factor_min = 0
   !control%u = 1e-8
   !control%small = 1e-16
   !control%nemin=1

   ! Read in a matrix
   write(*, "(a)") "Reading..."
   ! select case(matrix_type)
   ! case(HSL_MATRIX_REAL_SYM_PSDEF)
   !    control56%values = 3 ! make up values if necessary (posdef)
   ! case(HSL_MATRIX_REAL_SYM_INDEF)
   !    control56%values = 2 ! make up values if necessary (indef)
   ! end select
   ! call mc56_read(matfile, matrix, control56, info56)
   ! if(info56.ne.0) then
   !    print *, "mc56 failed with error ", info56
   !    stop
   ! endif
   if (options%fmt .eq. 'csc') then
      ! Rutherford boeing format

      ! DEBUG ensure matrix is diag dominant
      rb_options%values = 3 ! Force diagonal dominance
      call rb_read(matfile, m, n, ptr, row, val, rb_options, rb_flag)
      if(rb_flag.ne.0) then
         print *, "Rutherford-Boeing read failed with error ", rb_flag
         stop
      endif

   else if (options%fmt .eq. 'coo') then
      ! Matrix Market format

      ! read matrix
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

      ! print *, "m,n,nnz: ", m,n,nnz
      ! print *, "size ptr: ", size(ptr)
      ! print *, "size ptr: ", size(val)
      ! print *, "mm_flag: ", mm_flag
   else
      write(*, "(a)") "Matrix format not suported"
      stop
   end if
   write(*, "(a)") "ok"

   print *, "n = ", n

   ! check matrix format is correct
   call cscl_verify(6, SPRAL_MATRIX_REAL_SYM_PSDEF, n, n, &
        ptr, row, flag, more)
   if(flag.ne.0) then
      print *, "CSCL_VERIFY failed: ", flag, more
      stop
   endif

   write(*, "(a)") "Ordering..."
   call system_clock(start_t, rate_t)
   allocate(order(n), invp(n))
   call metis_order(n, ptr, row, order, invp, &
        flag, stat)
   ! call mc68_order(3, matrix%n, matrix%ptr, matrix%row, order, &
   !    order_control, order_info)
   ! !call mm_order(matrix%n, matrix%ptr, matrix%row, matrix%val, order)
   ! if(order_info%flag.ne.0) then
   !    print *, "mc68 failed with error ", order_info%flag
   !    stop
   ! endif
   ! call system_clock(stop_t)
   ! print *, "Ordering took ", (stop_t - start_t)/real(rate_t)
   ! !print *, "perm = ", order

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

   ! call mc69_verify(6, matrix_type, matrix%n, matrix%n, &
   !    matrix%ptr, matrix%row, flag, more)
   ! if(flag.ne.0) then
   !    print *, "MC69_VERIFY failed: ", flag, more
   !    stop
   ! endif

   ! ! Set options for analysis
   ! ssids_opt%ordering = 1 ! use Metis ordering
   ! ssids_opt%scaling = 0 ! no scaling

   ! ! Analyse SSIDS
   ! write(*, "(a)") "Analyse..."
   ! call system_clock(start_t, rate_t)
   ! call ssids_analyse(.false., n, ptr, row, akeep, &
   !      ssids_opt, inform, order, val=val)
   ! call system_clock(stop_t)
   ! print *, "Used order ", ssids_opt%ordering
   ! if (inform%flag < 0) then
   !    print *, "oops on analyse ", inform%flag
   !    stop
   ! endif
   ! write(*, "(a)") "ok"
   ! print *, "Analyse took ", (stop_t - start_t)/real(rate_t)
   ! !print *, "Used maximum memory of ", inform%maxmem
   ! smanal = (stop_t - start_t)/real(rate_t)
   ! print "(a,es10.2)", "Predict nfact = ", real(inform%num_factor)
   ! print "(a,es10.2)", "Predict nflop = ", real(inform%num_flops)
   ! print "(a6, i10)", "nparts", inform%nparts
   ! print "(a6, es10.2)", "cpu_flops", real(inform%cpu_flops)
   ! ! print "(a6, es10.2)", "gpu_flops", real(inform%gpu_flops)
   ! smaflop = real(inform%num_flops)
   ! smafact = real(inform%num_factor)

   ! Analyse and factor
   call system_clock(start_t, rate_t)
   call ma87_analyse(n, ptr, row, order, keep, control, &
      info)
   call system_clock(stop_t)
   if (info%flag < 0) then
      print *, "oops on analyse ", info%flag
      stop
   endif
   write(*, "(a)") "ok"
   print *, "Analyse took ", (stop_t - start_t)/real(rate_t)
   !print *, "Used maximum memory of ", info%maxmem
   smanal = (stop_t - start_t)/real(rate_t)
   print "(a,es10.2)", "Predict nfact = ", real(info%num_factor)
   print "(a,es10.2)", "Predict nflop = ", real(info%num_flops)
   smaflop = real(info%num_flops)
   smafact = real(info%num_factor)


   write(*, "(a)") "Factorize..."
   call system_clock(start_t, rate_t)
   do i = 1, nfact
      print *, "fact ", i
      call ma87_factor(n, ptr, row, val, order, &
         keep, control, info)
      !if(i.ne.nfact) call ma87_free(fkeep=fkeep)
   end do
   call system_clock(stop_t)
   call system_clock(tsf_en)
   if (info%flag < 0) then
      print *, "oops on factorize ", info%flag
      stop
   endif
   write(*, "(a)") "ok"
   print *, "Factor took ", (stop_t - start_t)/real(rate_t)
   smfact = (stop_t - start_t)/real(rate_t)

   ! Solve
   write(*, "(a)") "Solve..."
   call system_clock(start_t, rate_t)
   do i = 1, nslv
      soln = rhs
      !print *, "rhs = ", soln(1:matrix%n,1)
      call ma87_solve(soln(:,1), order, keep, control, info)
   end do
   call system_clock(stop_t)
   if (info%flag < 0) then
      print *, "oops on solve ", info%flag
      stop
   endif
   write(*, "(a)") "ok"
   print *, "Solve took ", (stop_t - start_t)/real(rate_t)

   !print *, "soln = ", soln(1:matrix%n,1)

   print *, "number bad cmp = ", count(abs(soln(1:n,1)-1.0).ge.1e-6)
   print *, "fwd error || ||_inf = ", maxval(abs(soln(1:n,1)-1.0))
   allocate(res(1))
    call internal_calc_norm(n, ptr, row, val, soln, rhs, nrhs, res)
    print *, "bwd error scaled = ", res

   call ma87_finalise(keep, control)

   print "(a6, a10)", "cmp:","SMFCT"
   print "(a6, f10.2)", "anal:", smanal
   print "(a6, f10.2)", "fact:", smfact
   print "(a6, es10.2)", "afact:", smafact
   print "(a6, es10.2)", "aflop:", smaflop
   print "(a6, es10.2)", "nfact:", real(info%num_factor)
   print "(a6, es10.2)", "nflop:", real(info%num_flops)

contains

  subroutine internal_calc_norm(n, ptr, row, val, x_vec, b_vec, nrhs, res)
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in) :: nrhs
    real(wp), dimension(nrhs*n), intent(in) :: x_vec
    real(wp), dimension(nrhs*n), intent(in) :: b_vec
    real(wp), dimension(nrhs), intent(out) :: res

    integer :: i, j, k, r
    double precision, allocatable, dimension(:) :: x_norm
    real(wp), dimension(:), allocatable :: res_vec
    double precision :: temp
    double precision :: normA

    ! Find the residual
    allocate(res_vec(n*nrhs), x_norm(nrhs))
    res_vec = 0
    do i = 1, n
       do j = ptr(i), ptr(i+1)-1
          r = row(j)
          do k = 0, nrhs-1
             res_vec(i+k*n) = res_vec(i+k*n) + &
                  val(j) * x_vec(r+k*n)
          end do
          if(r.eq.i) cycle
          do k = 0, nrhs-1
             res_vec(r+k*n) = res_vec(r+k*n) + &
                  val(j) * x_vec(i+k*n)
          end do
       end do
    end do
    res_vec(:) = res_vec(:) - b_vec(:)

    ! Find matrix norm
    call matrix_inf_norm(n, ptr, row, val, normA)

    ! Find x norm
    do i = 1, nrhs
       x_norm(i) = 0
       do j =1, n
          x_norm(i) = max(x_norm(i), abs(x_vec((i-1)*n+j)))
          if(x_vec((i-1)*n+j).ne.x_vec((i-1)*n+j)) then ! Tests for NaN
             x_norm(i) = x_vec((i-1)*n+j)
             exit
          endif
       end do
    end do

    ! Scaled residual = ||r|| / ( ||A|| ||x|| + ||b|| )
    do i = 1, nrhs
       temp = normA * x_norm(i) + &
            maxval(abs(b_vec((i-1)*n+1:i*n)))
       if(temp .eq. 0) then
          res(i) = maxval(abs(res_vec((i-1)*n+1:i*n)))
       else
          res(i) = maxval(abs(res_vec((i-1)*n+1:i*n))) / temp
       endif
    end do
  end subroutine internal_calc_norm

  subroutine matrix_inf_norm(n, ptr, row, val, norm)
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    real(wp), intent(out) :: norm

    real(wp), allocatable, dimension(:) :: row_norm
    integer :: i

    allocate(row_norm(n))

    row_norm = 0
    do i = 1, ptr(n+1)-1
       row_norm(row(i)) = row_norm(row(i)) + abs(val(i))
    end do

    norm = maxval(row_norm) 

  end subroutine matrix_inf_norm

end program
