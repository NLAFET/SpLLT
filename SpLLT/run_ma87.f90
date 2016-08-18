program run_prob
   use hsl_mc56_double
   use hsl_mc68_integer
   use hsl_mc69_double
   use hsl_zd11_double
   use hsl_ma87_double
   use spllt_mod
   use spllt_data_mod
   implicit none

   integer, parameter :: dp = kind(0d0)

   type(mc56_control) :: control56
   integer :: info56
   type(mc68_control) :: order_control
   type(mc68_info) :: order_info

   type(ZD11_type) :: matrix

   type(ma87_info) :: info
   type(ma87_keep) :: keep
   type(ma87_control) :: control
   double precision, dimension(:, :), allocatable :: rhs, soln
   double precision, dimension(:), allocatable :: work, res, scaling
   integer, dimension(:), allocatable :: order, perm
   integer, allocatable :: iwork(:)

   integer :: matrix_type
   integer :: i, j, k

   integer :: start_t, stop_t, rate_t, tsf_sa, tsf_en
   integer :: flag, more

   integer, parameter :: unit_rhs = 14

   real :: smanal, smfact, smaflop, smafact

   integer, parameter :: nfact = 1
   !integer, parameter :: nfact = 5
   !integer, parameter :: nfact = 100

   integer, parameter :: nslv = 1
   !integer, parameter :: nslv = 10
   !integer, parameter :: nslv = 100
   type(spllt_options) :: options
   character(len=200) :: matfile    

   matrix_type = HSL_MATRIX_REAL_SYM_PSDEF

   call splllt_parse_args(options)

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
   select case(matrix_type)
   case(HSL_MATRIX_REAL_SYM_PSDEF)
      control56%values = 3 ! make up values if necessary (posdef)
   case(HSL_MATRIX_REAL_SYM_INDEF)
      control56%values = 2 ! make up values if necessary (indef)
   end select
   call mc56_read(matfile, matrix, control56, info56)
   if(info56.ne.0) then
      print *, "mc56 failed with error ", info56
      stop
   endif
   write(*, "(a)") "ok"

   print *, "n = ", matrix%n

   call mc69_cscl_clean(matrix_type, matrix%n, matrix%n, &
      matrix%ptr, matrix%row, flag, val=matrix%val)

   write(*, "(a)") "Ordering..."
   call system_clock(start_t, rate_t)
   allocate(order(matrix%n))

   call mc68_order(3, matrix%n, matrix%ptr, matrix%row, order, &
      order_control, order_info)
   !call mm_order(matrix%n, matrix%ptr, matrix%row, matrix%val, order)
   if(order_info%flag.ne.0) then
      print *, "mc68 failed with error ", order_info%flag
      stop
   endif
   call system_clock(stop_t)
   print *, "Ordering took ", (stop_t - start_t)/real(rate_t)
   !print *, "perm = ", order

   ! Make up a rhs
   allocate(rhs(matrix%n, 1), soln(matrix%n, 1))
   rhs(:,:) = 0
   do i = 1, matrix%n
      do j = matrix%ptr(i), matrix%ptr(i+1)-1
         k = matrix%row(j)
         rhs(k, 1) = rhs(k, 1) + matrix%val(j)
         if(i.eq.k) cycle
         rhs(i, 1) = rhs(i, 1) + matrix%val(j)
      end do
   end do

   call mc69_verify(6, matrix_type, matrix%n, matrix%n, &
      matrix%ptr, matrix%row, flag, more)
   if(flag.ne.0) then
      print *, "MC69_VERIFY failed: ", flag, more
      stop
   endif

   ! Analyse and factor
   allocate(scaling(matrix%n))
   call system_clock(start_t, rate_t)
   call ma87_analyse(matrix%n, matrix%ptr, matrix%row, order, keep, control, &
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
      call ma87_factor(matrix%n, matrix%ptr, matrix%row, matrix%val, order, &
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

   print *, "number bad cmp = ", count(abs(soln(1:matrix%n,1)-1.0).ge.1e-6)
   print *, "fwd error || ||_inf = ", maxval(abs(soln(1:matrix%n,1)-1.0))
   allocate(res(1))
   call internal_calc_norm(matrix, soln, rhs, 1, res)
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
   subroutine internal_calc_norm(matrix, x_vec, b_vec, nrhs, res)
      type(ZD11_type), intent(in) :: matrix
      integer, intent(in) :: nrhs
      real(dp), dimension(nrhs*matrix%n), intent(in) :: x_vec
      real(dp), dimension(nrhs*matrix%n), intent(in) :: b_vec
      real(dp), dimension(nrhs), intent(out) :: res

      integer :: i, j, k, row
      double precision, allocatable, dimension(:) :: x_norm
      real(dp), dimension(:), allocatable :: res_vec
      double precision :: temp
      double precision :: normA

      ! Find the residual
      allocate(res_vec(matrix%n*nrhs), x_norm(nrhs))
      res_vec = 0
      do i = 1, matrix%n
         do j = matrix%ptr(i), matrix%ptr(i+1)-1
            row = matrix%row(j)
            do k = 0, nrhs-1
               res_vec(i+k*matrix%n) = res_vec(i+k*matrix%n) + &
                  matrix%val(j) * x_vec(row+k*matrix%n)
            end do
            if(row.eq.i) cycle
            do k = 0, nrhs-1
               res_vec(row+k*matrix%n) = res_vec(row+k*matrix%n) + &
                  matrix%val(j) * x_vec(i+k*matrix%n)
            end do
         end do
      end do
      res_vec(:) = res_vec(:) - b_vec(:)

      ! Find matrix norm
      call matrix_inf_norm(matrix, normA)

      ! Find x norm
      do i = 1, nrhs
         x_norm(i) = 0
         do j =1, matrix%n
            x_norm(i) = max(x_norm(i), abs(x_vec((i-1)*matrix%n+j)))
            if(x_vec((i-1)*matrix%n+j).ne.x_vec((i-1)*matrix%n+j)) then ! Tests for NaN
               x_norm(i) = x_vec((i-1)*matrix%n+j)
               exit
            endif
         end do
      end do

      print *, "||r|| = ", maxval(abs(res_vec(1:matrix%n)))
      print *, "||A|| = ", normA
      print *, "||x|| = ", x_norm(1)
      print *, "||b|| = ", maxval(abs(b_vec(1:matrix%n)))

      ! Scaled residual = ||r|| / ( ||A|| ||x|| + ||b|| )
      do i = 1, nrhs
         temp = normA * x_norm(i) + &
            maxval(abs(b_vec((i-1)*matrix%n+1:i*matrix%n)))
         if(temp .eq. 0) then
            res(i) = maxval(abs(res_vec((i-1)*matrix%n+1:i*matrix%n)))
         else
            res(i) = maxval(abs(res_vec((i-1)*matrix%n+1:i*matrix%n))) / temp
         endif
      end do
   end subroutine internal_calc_norm

   subroutine matrix_inf_norm(matrix, norm)
      type(ZD11_type), intent(in) :: matrix
      real(dp), intent(out) :: norm

      real(dp), allocatable, dimension(:) :: row_norm
      integer :: i

      allocate(row_norm(matrix%n))

      row_norm = 0
      do i = 1, matrix%ptr(matrix%n+1)-1
         row_norm(matrix%row(i)) = row_norm(matrix%row(i)) + &
            abs(matrix%val(i))
      end do

      norm = maxval(row_norm) 
   end subroutine matrix_inf_norm
end program
