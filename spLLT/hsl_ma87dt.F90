program main

! This program generates two files (we.out and dl.out) in addition to stdout
! The output in hsl_ma87dt.output is generated as follows:
!    ./hsl_ma87dt > temp1
!    cat dl.out we.out temp1 > hsl_ma87dt.output


!! Note: output is different if code run on single thread since
!        we set the max problem size to be smaller.


! To convert from double to single:
! * Change wp
! * Change _double
! * Change HSL names: fa14id, mc47ad, mc47id, ym11ad, ym11id
! * Change err_tol to 1e-6

!$ use omp_lib
   use hsl_fa14_double
   use hsl_ma87_double
   use hsl_zd11_double
   implicit none

   integer,  parameter :: wp = kind(0d0)
   real(wp), parameter :: err_tol = 1e-12
   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: zero = 0.0_wp

   integer :: errors

   integer, parameter :: we_unit = 11
   character(len=6)   :: we_file = "we.out"
   integer, parameter :: dl_unit = 12
   character(len=6)   :: dl_file = "dl.out"

   ! Error flags
   integer, parameter :: MA87_SUCCESS               = 0
   integer, parameter :: MA87_ERROR_ALLOCATION      = -1
   integer, parameter :: MA87_ERROR_ORDER           = -2
   integer, parameter :: MA87_ERROR_NOT_POSDEF      = -3
   integer, parameter :: MA87_ERROR_X_SIZE          = -4
   integer, parameter :: MA87_ERROR_INFINITY        = -5
   integer, parameter :: MA87_ERROR_JOB_OOR         = -6
   integer, parameter :: MA87_ERROR_NBI_OOR         = -7
   integer, parameter :: MA87_ERROR_UNKNOWN         = -99

   ! warning flags
   integer, parameter :: MA87_WARNING_POOL_SMALL    = 1


   if(we_unit.gt.6) open(unit=we_unit,file=we_file,status="replace")
   if(dl_unit.gt.6) open(unit=dl_unit,file=dl_file,status="replace")

   errors = 0

   call test_bad_args
   call test_not_posdef
   call test_warnings
   call test_random
   call test_random_sparse

   write(*, "(/a)") "=========================="
   write(*, "(a,i4)") "Total number of errors = ", errors

   if(we_unit.gt.6) close(we_unit)
   if(dl_unit.gt.6) close(dl_unit)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_random
   type(ma87_keep) :: keep
   type(ma87_control) :: control
   type(ma87_info) :: info


   integer :: maxn = 1000
   integer :: maxnb = 160
   integer :: maxnz =  1000000
   integer, parameter :: maxnrhs = 10
   integer, parameter :: nprob = 100

   type(fa14_seed) :: pseed

   type(ZD11_type) :: a
   real(wp), allocatable, dimension(:, :) :: rhs,x
   real(wp), allocatable, dimension(:) :: rhs1d,x1
   real(wp), allocatable, dimension(:, :) :: res

   integer :: iseed, nza, prblm, i, j, k, nrhs, st
   integer :: threads ! Number of threads
   integer, dimension(:), allocatable :: order
   real(wp) :: num_flops

   write(*, "(a)")
   write(*, "(a)") "======================="
   write(*, "(a)") "Testing random matrices"
   write(*, "(a)") "======================="

   control%unit_error = we_unit
   control%unit_warning = we_unit

   call fa14id(iseed)

   threads = 1
!$ threads = omp_get_max_threads()
   if (threads==1) then
      maxn = 100
      maxnz = 10000
      maxnb = 8
      control%cache_tq_sz = 10
   end if

   deallocate(a%ptr,a%row,a%val,a%col,order,stat=st)
   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz), a%col(2*maxnz))
   allocate(order(maxn))
   deallocate(rhs,rhs1d,x,x1,stat=st)
   allocate(rhs(maxn, maxnrhs),x(maxn, maxnrhs))
   allocate(x1(maxn),rhs1d(maxn))

   do prblm = 1, nprob

      ! Generate parameters
      call FA14_random_integer(pseed, maxn, a%n)
      if (prblm < 12) a%n = prblm ! check very small problems
      control%unit_diagnostics = dl_unit
      control%diagnostics_level = 1
      if (prblm == 11) then ! full printing for a small problem
        control%diagnostics_level = 3
      end if
      if (prblm == 12) a%n = 0
      i = a%n**2/2 - a%n
      i = max(0,i)
      call FA14_random_integer(pseed, i, nza)
      nza = nza + a%n

      call FA14_random_integer(pseed, max(maxnb,a%n/4), control%nb)
      control%nb = max(control%nb,2)
      if(prblm.ge.nprob-2) then
         ! Last two problems, enable lots of logging and debugging stuff
         control%unit_diagnostics = dl_unit
         control%diagnostics_level = 1
         if(prblm.eq.nprob) control%diagnostics_level = 2
         control%nb = 4
         control%cache_layout = 2 ! CACHE_SCATTER
         if(prblm.eq.nprob) control%cache_layout = 3 ! CACHE_IDENTITY
         a%n = min(maxn,a%n + 400)
         nza = min(nza + 400,maxnz)
      endif

      if(nza.gt.maxnz .or. a%n.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5)") "n = ", a%n, " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      write(*, "(a, i3, a, i5, a, i8, a)",advance="no") " * number ", &
         prblm,  &
         " n = ", a%n, " nza = ", nza, "..."


      call gen_random_posdef(a, nza, iseed)

      ! Generate a pivot order
      call amd_order(a, order)

      ! Peform analyse
      call MA87_analyse(a%n, a%ptr, a%row, order, keep, control, info)
      num_flops = info%num_flops
      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a)") "fail on analyse"
         errors = errors + 1
         cycle
      endif
      !print *, "nfact =", info%nfactor

      call FA14_random_integer(pseed, maxnrhs, nrhs)

      ! Generate rhs assuming x(k) = k/maxn. remember we have only
      ! half matrix held.
      rhs(1:a%n, 1:nrhs) = zero
      do k = 1, a%n
         do j = a%ptr(k), a%ptr(k+1)-1
            i = a%row(j)
            rhs(i, 1:nrhs) = rhs(i, 1:nrhs) + a%val(j)*real(k)/real(maxn)
            if(i.eq.k) cycle
            rhs(k, 1:nrhs) = rhs(k, 1:nrhs) + a%val(j)*real(i)/real(maxn)
         end do
      end do
      rhs1d(1:a%n) = rhs(1:a%n, 1)
      
      ! Perform straight forward factor then solve

      ! set right-hand side
      x(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
      x1(1:a%n) = rhs1d(1:a%n)

      call MA87_factor(a%n, a%ptr, a%row, a%val, order, keep, control, info)

      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a)") "fail on factor"
         errors = errors + 1
         cycle
      endif
      write(*,'(a,f6.1,1x)',advance="no") ' num_flops:',num_flops*1e-6


      ! Perform solve
      if (mod(a%n,2).eq.0) then
         if (a%n > maxn/2) then
            ! Even numbered problem, 50% of time at random
            ! Vanilla call

            call MA87_solve(x1, order, keep, control, info)
            if(info%flag .lt. MA87_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job absent", &
                  info%flag
               errors = errors + 1
               cycle
            endif

            call MA87_solve(nrhs, maxn, x, order, keep, &
               control, info)
            if(info%flag .lt. MA87_SUCCESS) then
               write(*, "(a,i4)") "fail on solve with job absent", info%flag
               errors = errors + 1
               cycle
            endif

         else
            ! Even numbered problem, 50% of time at random
            ! combined fwd+back with job = 0

            call MA87_solve(x1, order, keep,control, info, job=0)
            if(info%flag .lt. MA87_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job = 0", info%flag
               errors = errors + 1
               cycle
            endif

            call MA87_solve(nrhs, maxn, x, order, keep,control, info, job=0)
            if(info%flag .ne. MA87_SUCCESS) then
               write(*, "(a,i4)") "fail on solve with job = 0", info%flag
               errors = errors + 1
               cycle
            endif
         endif

      else
         ! Odd numbered problems
         ! seperate fwd and bwd calls

         call MA87_solve(x1, order, keep, control, info, job=1)
         if(info%flag .lt. MA87_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 1", info%flag
            errors = errors + 1
            cycle
         endif

         call MA87_solve(nrhs, maxn, x, order, keep, &
            control, info, job=1)
         if(info%flag .lt. MA87_SUCCESS) then
            write(*, "(a,i4)") "fail on solve with job = 1", info%flag
            errors = errors + 1
            cycle
         endif

         call MA87_solve(x1, order, keep, control, info, job=2)
         if(info%flag .lt. MA87_SUCCESS) then
             write(*, "(a,i4)") " fail on 1d solve with job = 2", info%flag
             errors = errors + 1
             cycle
         endif

         call MA87_solve(nrhs, maxn, x, order, keep, &
               control, info, job=2)
         if(info%flag .lt. MA87_SUCCESS) then
            write(*, "(a,i4)") "fail on solve with job = 2", info%flag
            errors = errors + 1
            cycle
         endif

      end if ! odd or even problem

      ! write (6,'(6es12.4)') x(1:a%n,1:1)

      ! Check residuals
      call compute_resid(1,a,x1,maxn,rhs1d,maxn,res)
      if(maxval(abs(res(1:a%n,1))) < err_tol) then
         write(*, "(a)",advance="no") "ok..."
      else
         write(*, "(a,es12.4)") " fail residual 1d = ", &
            maxval(abs(res(1:a%n,1)))
         errors = errors + 1
      endif

      call compute_resid(nrhs,a,x,maxn,rhs,maxn,res)

      if(maxval(abs(res(1:a%n,1:nrhs))) < err_tol) then
         write(*, "(a)",advance="no") "ok..."
      else
         write(*, "(a,es12.4)") "fail residual = ", &
            maxval(abs(res(1:a%n,1:nrhs)))
         errors = errors + 1
      endif

      !
      ! Use factor_solve
      !

      ! set right-hand side
      x(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
      x1(1:a%n) = rhs1d(1:a%n)

      ! perform factorization and solve (single rhs)
      call ma87_factor_solve(a%n, a%ptr, a%row, a%val, order, keep, control, &
         info, x1)

      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a)") " fail on factor_solve (single rhs)"
         errors = errors + 1
         cycle
      endif

      call compute_resid(1,a,x1,maxn,rhs1d,maxn,res)
      if(maxval(abs(res(1:a%n,1))) < err_tol) then
         write(*, "(a)",advance="no") "ok..."
      else
         write(*, "(a,es12.4)") " 1d factor_solve fail residual = ", &
            maxval(abs(res(1:a%n,1)))
         errors = errors + 1
      endif

      ! perform factorization and solve (multiple rhs)
      call ma87_factor_solve(a%n, a%ptr, a%row, a%val, order, keep, control,  &
         info, nrhs, maxn, x)
      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a)") " fail on factor_solve"
         errors = errors + 1
         cycle
      endif
      call compute_resid(nrhs,a,x,maxn,rhs,maxn,res)
      if(maxval(abs(res(1:a%n,1:nrhs))) < err_tol) then
         write(*, "(a)") "ok"
      else
         write(*, "(a,es12.4)") " factor_solve fail residual = ", &
            maxval(abs(res(1:a%n,1:nrhs)))
         errors = errors + 1
      endif

      call MA87_finalise(keep, control)

   end do

end subroutine test_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_random_sparse

   ! tests with sparse right-hand side.
   ! use ma87_sparse_fwd_solve for forward substitution with sparse rhs
   ! (single rhs only)
   type(ma87_keep) :: keep
   type(ma87_control) :: control
   type(ma87_info) :: info


   integer :: maxn = 1000
   integer :: maxnb = 32
   integer :: maxnz =  1000000
   integer, parameter :: nprob = 20

   type(fa14_seed) :: pseed

   type(ZD11_type) :: a
   real(wp), allocatable, dimension(:) :: rhs,x,w
   real(wp), allocatable, dimension(:,:) :: res

   integer :: iseed, nza, nbi, nxi, prblm, i, j, k, st
   integer :: threads ! Number of threads
   integer, dimension(:), allocatable :: order,invp,index,bindex
   real(wp) :: num_flops

   write(*, "(a)")
   write(*, "(a)") "======================="
   write(*, "(a)") "Testing sparse rhs"
   write(*, "(a)") "======================="

   control%unit_error = we_unit
   control%unit_warning = we_unit

   call fa14id(iseed)

   threads = 1
!$ threads = omp_get_max_threads()
   if (threads==1) then
      maxn = 100
      maxnz = 10000
      maxnb = 8
      control%cache_tq_sz = 10
   end if

   deallocate(a%ptr,a%row,a%val,a%col,order,index,bindex,stat=st)
   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz), a%col(2*maxnz))
   allocate(order(maxn),invp(maxn),bindex(maxn),index(maxn))
   deallocate(rhs,x,w,stat=st)
   allocate(x(maxn),rhs(maxn),w(maxn))

   do prblm = 1, nprob

      ! Generate parameters
      call FA14_random_integer(pseed, maxn, a%n)
      if (prblm < 10) a%n = prblm ! check very small problems
      control%unit_diagnostics = dl_unit
      control%diagnostics_level = 1
      i = a%n**2/2 - a%n
      i = max(0,i)
      call FA14_random_integer(pseed, i, nza)
      nza = nza + a%n

      call FA14_random_integer(pseed, maxnb, control%nb)
      control%nb = max(control%nb,2)

      if(nza.gt.maxnz .or. a%n.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5)") "n = ", a%n, " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      write(*, "(a, i3, a, i5, a, i8, a)",advance="no") " * number ", &
         prblm,  &
         " n = ", a%n, " nza = ", nza, "..."

      call gen_random_posdef(a, nza, iseed)

      ! Generate a pivot order
      call amd_order(a, order)

      if (prblm .eq. nprob) then
        do i = 1,a%n
          order(i) = i
        end do
      end if

      ! Peform analyse

      call MA87_analyse(a%n, a%ptr, a%row, order, keep, control, info)
      num_flops = info%num_flops
      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a)") "fail on analyse"
         errors = errors + 1
         cycle
      endif
      !print *, "nfact =", info%nfactor
      ! write (*,*) order(1:a%n)

      call MA87_factor(a%n, a%ptr, a%row, a%val, order, keep, control, info)

      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a)") "fail on factor"
         errors = errors + 1
         cycle
      endif
 !     write(*,'(a,f6.1,1x)',advance="no") ' num_flops:',num_flops*1e-6

      ! Generate sparse rhs. nbi is number of non-zeros in rhs.
      rhs(1:a%n) = zero
      call FA14_random_integer(pseed, a%n/12, nbi)
      nbi = max(1,nbi)
      if (prblm .eq. nprob) then
         ! test with single non zero entry in last position.
         nbi = 1
         bindex(1) = a%n
         call FA14_random_real(pseed, .false., rhs(a%n))
      else
         j = 0
         do while (j < nbi)
            call FA14_random_integer(pseed, a%n, k)
            if (rhs(k).eq.zero) then
               j = j + 1
               bindex(j) = k
               call FA14_random_real(pseed, .false., rhs(k))
            end if
         end do
      end if

      ! Perform solve. Sparse forward subsitution.

      ! inverse permutation
      do i = 1,a%n
        j = order(i)
        invp(j) = i
      end do

      x(1:a%n) = zero

         call MA87_sparse_fwd_solve(nbi,bindex,rhs,order,invp,nxi,index,x,w, &
            keep,control,info)
         if(info%flag .lt. MA87_SUCCESS) then
            write(*, "(a,i4)") " fail on sparse_fwd_solve", info%flag
            errors = errors + 1
            cycle
         endif
         write(*,'(a,i3,1x,i4,1x)',advance="no") ' nbi,nxi = ',nbi,nxi
 

         ! back substitution
         call MA87_solve(x, order, keep, control, info, job=2)
         if(info%flag .lt. MA87_SUCCESS) then
             write(*, "(a,i4)") " fail on solve with job = 2", info%flag
             errors = errors + 1
             cycle
         endif


      ! write (6,'(6es12.4)') x(1:a%n)

      ! Check residuals
      call compute_resid(1,a,x,maxn,rhs,maxn,res)
      if(maxval(abs(res(1:a%n,1))) < err_tol) then
         write(*, "(a)") "ok"
      else
         write(*, "(a,es12.4)") " fail residual = ", &
            maxval(abs(res(1:a%n,1)))
         errors = errors + 1
      endif

      call MA87_finalise(keep, control)

   end do

end subroutine test_random_sparse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine amd_order(a,order)
   type(zd11_type), intent(in) :: a
   integer, dimension(:), allocatable :: order

   logical :: realloc_flag
   integer :: i, st
   integer, dimension(10) :: icntl, info
   real(wp), dimension(10) :: rinfo
   integer, dimension(:), allocatable :: work, ptr

   realloc_flag = .true.
   if(allocated(order)) realloc_flag = size(order).lt.a%n

   if(realloc_flag) then
      deallocate(order,stat=st)
      allocate(order(a%n))
   endif

   ! Initilise control
   call mc47id(icntl)
   icntl(1:2) = -1 ! Supress warnings and errors
   icntl(5) = huge(0) ! Largest integer

   ! Copy ptr data to work array
   allocate(ptr(a%n+1))
   ptr(:) = a%ptr(1:a%n+1)
   ! Copy row data to work array
   allocate(work(2*a%ptr(a%n+1) + 10*a%n))
   work(1:a%ptr(a%n+1)-1) = &
      a%row(1:a%ptr(a%n+1)-1)

   ! Perform AMD
   call mc47ad(a%n, a%ptr(a%n+1)-1, ptr, work, &
      size(work), icntl, info, rinfo)
   if(info(1).lt.0) then
      ! Failed for some reason
      do i = 1, a%n
         order(i) = i
      end do
      return
   endif

   ! Extract ordering
   order(1:a%n) = work(size(work)-a%n+1:size(work))

end subroutine amd_order


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_random_posdef(matrix, nza, iseed)
   type(zd11_type), intent(inout) :: matrix
   integer, intent(in) :: nza
   integer, intent(inout) :: iseed

   integer, dimension(10) :: icntl
   character(len=8) :: key
   integer, dimension(:), allocatable :: work
   integer :: i, j, k
   real(wp) :: tempv

   ! Setup stuff for ym11
   key = 'nothing '
   call ym11id(icntl, i)
   icntl(3) = 0 ! Symmetric
   allocate(work(2*matrix%n))

   ! Generate matrix
   call ym11ad(matrix%n, matrix%n, nza, matrix%ne, matrix%row, &
      matrix%val, matrix%ptr, work, icntl, key, iseed)

   ! Make matrix diagonally dominant, observing first entry in column
   ! is always the diagonal when matrix generated by ym11
   do k = 1, matrix%n
      tempv = zero
      do j = matrix%ptr(k)+1, matrix%ptr(k+1)-1
         tempv = tempv + abs(matrix%val(j))
         i = matrix%ptr(matrix%row(j))
         matrix%val(i) = matrix%val(i) + abs(matrix%val(j))
      end do
      i = matrix%ptr(k)
      matrix%val(i) = one + matrix%val(i) + tempv
   end do
end subroutine gen_random_posdef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_bad_args

   type(zd11_type) :: a
   type(MA87_control) :: control
   type(MA87_keep) :: keep
   type(MA87_info) :: info

   integer :: i
   integer :: lx
   integer :: nbi,nrhs,nxi
   integer :: st
   integer, dimension(:), allocatable :: order,invp,index,bindex
   real(wp), dimension(:,:), allocatable :: x
   real(wp), dimension(:), allocatable :: b,x1,w

   control%unit_error = we_unit
   control%unit_warning = we_unit

   write(*,*)
   write(*,"(a)") "======================"
   write(*,"(a)") "Testing bad arguments:"
   write(*,"(a)") "======================"

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call simple_mat(a)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! tests on call to ma87_analyse

   write(*,"(a)",advance="no") " * Testing order too short..................."
   deallocate(order,stat=st)
   allocate(order(a%n-1))
   order(1:a%n-1) = 1
   call MA87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_ORDER)

   write(*,"(a)",advance="no") " * Recall MA87_analyse after error..........."
   ! recall MA87_analyse so that we test returning after previous error
   call MA87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_ORDER)

   deallocate(order,stat=st)
   allocate(order(a%n))
   do i = 1, a%n
      order(i) = i
   end do

   write(*,"(a)",advance="no") " * Testing order with repeats................"
   order(2) = order(1)
   call MA87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_ORDER)
   order(2) = 2

   write(*,"(a)",advance="no") " * Testing order with oor below.............."
   order(2) = 0
   call MA87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_ORDER)
   order(2) = 2

   write(*,"(a)",advance="no") " * Testing order with oor above.............."
   order(2) = a%n
   call MA87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_ORDER)
   order(2) = 2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! tests on call to ma87_factor_solve

   ! Run analyse properly
   call MA87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   if(info%flag.ne.0) then
      write(*, "(a,i4)") &
     "Unexpected error during setup in analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif
   nrhs = 2
   lx = a%n
   deallocate(x,stat=st)
   allocate(x(lx,nrhs))

   write(*,"(a)",advance="no") " * Testing factor_solve with nrhs=0.........."
   call MA87_factor_solve(a%n,a%ptr,a%row,a%val,order,keep,control,info, &
        0,lx,x)
   call print_result(info%flag, MA87_ERROR_X_SIZE)

   write(*,"(a)",advance="no") " * Testing factor_solve with lx=a%n-1........"
   deallocate(x,stat=st)
   lx = a%n-1
   allocate(x(lx,nrhs))
   call MA87_factor_solve(a%n,a%ptr,a%row,a%val,order,keep,control,info, &
         nrhs,lx,x)
   call print_result(info%flag, MA87_ERROR_X_SIZE)
   deallocate(x)
   lx = a%n
   allocate(x(lx,nrhs))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! tests on call to ma87_solve

   ! Run analyse and factorise properly
   call MA87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   if(info%flag.ne.0) then
      write(*, "(a,i4)") &
     "Unexpected error during setup in analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif
   call MA87_factor(a%n,a%ptr,a%row,a%val,order,keep,control,info)
   if(info%flag.ne.0) then
      write(*, "(a,i4)") &
     "Unexpected error during setup in factor. flag = ", info%flag
      errors = errors + 1
      return
   endif

   nrhs = 2
   lx = a%n
   deallocate(x)
   allocate(x(lx,nrhs))

   write(*,"(a)",advance="no") " * Testing solve with nrhs=0................."
   call MA87_solve(0,lx,x,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_X_SIZE)

   write(*,"(a)",advance="no") " * Testing solve with lx=a%n-1..............."
   deallocate(x,stat=st)
   lx = a%n-1
   allocate(x(lx,nrhs))
   call MA87_solve(nrhs,lx,x,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_X_SIZE)

   deallocate(x,stat=st)
   lx = a%n
   allocate(x(lx,nrhs))

   write(*,"(a)",advance="no") " * Testing solve with job=-1................."
   call MA87_solve(nrhs,lx,x,order,keep,control,info,job=-1)
   call print_result(info%flag, MA87_ERROR_JOB_OOR)

   write(*,"(a)",advance="no") " * Testing solve with job=5.................."
   call MA87_solve(nrhs,lx,x,order,keep,control,info,job=5)
   call print_result(info%flag, MA87_ERROR_JOB_OOR)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! tests on call to ma87_sparse_fwd_solve

   deallocate(x1,stat=st)
   deallocate(w,stat=st)
   deallocate(b,stat=st)
   deallocate(bindex,stat=st)
   deallocate(index,stat=st)
   deallocate(invp,stat=st)
   allocate(b(a%n),w(a%n),x1(a%n),index(a%n),invp(a%n),bindex(a%n))

   write(*,"(a)",advance="no") " * Testing sparse_fwd_solve with nbi=0......."
   nbi = 0
   call MA87_sparse_fwd_solve(nbi,bindex,b,order,invp,nxi,index,x1,w,keep, &
      control,info)
   call print_result(info%flag, MA87_ERROR_NBI_OOR)

end subroutine test_bad_args

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_not_posdef

   type(zd11_type) :: a
   type(ma87_control) :: control
   type(ma87_keep) :: keep
   type(ma87_info) :: info

   integer :: i, lx, nrhs
   integer, dimension(:), allocatable :: order
   real(wp), dimension(:,:), allocatable :: x

   control%unit_error = we_unit

   write(*,"(a)")
   write(*,"(a)") "============================="
   write(*,"(a)") "Testing not positive definite"
   write(*,"(a)") "============================="

   call simple_mat(a)

   allocate(order(a%n))
   do i = 1, a%n
      order(i) = i
   end do
   call ma87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   if(info%flag.ne.0) then
      write(*, "(a,i4)") &
     "Unexpected error during setup in analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif
   do i = 1, a%ptr(a%n+1)-1
      a%val(i) = -a%val(i)
   end do
   call ma87_factor(a%n,a%ptr,a%row,a%val,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_NOT_POSDEF)

   write(*,"(a)",advance="no") " * Recall ma87_factor after error......"
   do i = 1, a%ptr(a%n+1)-1
      a%val(i) = -a%val(i)
   end do
   call ma87_factor(a%n,a%ptr,a%row,a%val,order,keep,control,info)
   call print_result(info%flag, MA87_SUCCESS)

   write(*,"(a)",advance="no") " * Call ma87_solve after error........."
   do i = 1, a%ptr(a%n+1)-1
      a%val(i) = -a%val(i)
   end do
   call ma87_factor(a%n,a%ptr,a%row,a%val,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_NOT_POSDEF)

   lx = a%n
   nrhs = 1
   allocate(x(lx,nrhs))
   call ma87_solve(nrhs,lx,x,order,keep,control,info)
   call print_result(info%flag, MA87_ERROR_NOT_POSDEF)
   deallocate(x)

   call ma87_finalise(keep,control)

end subroutine test_not_posdef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_warnings

   type(zd11_type) :: a
   type(MA87_control) :: control
   type(MA87_keep) :: keep
   type(MA87_info) :: info

   integer :: i
   integer :: temp, temp2
   real(wp), dimension(:,:), allocatable :: res
   real(wp), dimension(:,:), allocatable :: rhs
   real(wp), dimension(:,:), allocatable :: x
   integer, dimension(:), allocatable :: order

   write(*,"(a)")
   write(*,"(a)") "================"
   write(*,"(a)") "Testing warnings"
   write(*,"(a)") "================"

   control%unit_warning = we_unit
   control%unit_diagnostics = dl_unit


   write(*,"(a)",advance="no") " * Testing task pool overflow............"
   control%diagnostics_level = 1

   a%n = 3000
   temp = (a%n * 7)/2
   allocate(a%ptr(a%n+1), a%row(temp*2), a%val(temp*2))
   call fa14id(temp2)
   call gen_random_posdef(a, temp, temp2)
   if (allocated(rhs)) deallocate(rhs)
   allocate(rhs(a%n,1))
   call gen_rhs(a,rhs,x)
   control%pool_size = 1
   control%cache_tq_sz = 2
   control%nb = 160

   allocate(order(a%n))
   do i = 1, a%n
      order(i) = i
   end do
   call MA87_analyse(a%n,a%ptr,a%row,order,keep,control,info)
   if(info%flag .ne. MA87_SUCCESS) then
      write(*, "(a,i4)") "fail on analyse", info%flag
      errors = errors + 1
      if (info%flag < 0) return
   endif

   call MA87_factor(a%n,a%ptr,a%row,a%val,order,keep,control,info)
   call print_result(info%flag,MA87_WARNING_POOL_SMALL)
   if(info%flag .lt. MA87_SUCCESS) then
      write(*, "(a,i4)") "fail on factor", info%flag
      errors = errors + 1
      return
   endif

   write(*,"(a)",advance="no") " *    checking answer...................."
   call MA87_solve(1, a%n, x, order, keep, control, info)
   if(info%flag .lt. MA87_SUCCESS) then
      write(*, "(a,i4)") "fail on solve", info%flag
      errors = errors + 1
      return
   endif
   call compute_resid(1,a,x,a%n,rhs,a%n,res)
   if(maxval(abs(res(1:a%n,1))) < err_tol) then
      write(*, "(a)") "ok"
   else
      write(*, "(a,es12.4)") "fail residual = ", &
         maxval(abs(res(1:a%n,1)))
      errors = errors + 1
   endif

   call MA87_finalise(keep,control)
end subroutine test_warnings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_rhs(a, rhs, x)
   type(zd11_type), intent(inout) :: a
   real(wp), dimension(:,:), allocatable, intent(inout) :: rhs
   real(wp), dimension(:,:), allocatable, intent(inout) :: x

   integer :: j, k

   if (allocated(rhs)) deallocate(rhs)
   allocate(rhs(a%n, 1))

   if (allocated(x)) deallocate(x)
   allocate(x(a%n, 1))

   ! Generate rhs assuming x = 1
   rhs = zero
   do k = 1, a%n
      do j = a%ptr(k), a%ptr(k+1)-1
         if(a%row(j).lt.k .or. a%row(j).gt.a%n) cycle
         rhs(k, 1) = rhs(k, 1) + a%val(j)
         if(a%row(j).eq.k) cycle
         rhs(a%row(j), 1) = rhs(a%row(j), 1) + a%val(j)
      end do
   end do
   x = rhs
end subroutine gen_rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chk_answer(a, keep, control, rhs, x, fs)

   type(zd11_type), intent(inout) :: a
   type(MA87_keep), intent(inout) :: keep
   type(MA87_control), intent(in) :: control
   real(wp), dimension(:,:), intent(inout) :: rhs
   real(wp), dimension(:,:), intent(inout) :: x
   logical, optional, intent(in) :: fs

   type(MA87_info) :: info
   integer, dimension(:), allocatable :: order
   real(wp), dimension(:,:), allocatable :: res
   integer :: i

   write(*,"(a)",advance="no") " *    checking answer...................."

   allocate(order(a%n))
   do i = 1, a%n
      order(i) = i
   end do

   call MA87_analyse(a%n, a%ptr, a%row, order,keep,control,info)
   if(info%flag .ne. MA87_SUCCESS) then
      write(*, "(a)") "fail on analyse"
      errors = errors + 1
      return
   endif

   if(present(fs)) then
      call MA87_factor_solve(a%n, a%ptr, a%row, a%val, &
           order,keep,control,info, 1, a%n, x)
      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a)") "fail on factor_solve"
         errors = errors + 1
         return
      endif
   else
      call MA87_factor(a%n, a%ptr, a%row, a%val, &
           order,keep,control,info)
      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a)") "fail on factor"
         errors = errors + 1
         return
      endif

      call MA87_solve(1, a%n, x, order, keep, control, info)
      if(info%flag .lt. MA87_SUCCESS) then
         write(*, "(a,i4)") "fail on solve", info%flag
         errors = errors + 1
         return
      endif
   endif

 ! Check residual
   call compute_resid(1,a,x,a%n,rhs,a%n,res)

   if(maxval(abs(res(1:a%n,1))) < err_tol) then
      write(*, "(a)") "ok"
   else
      write(*, "(a,es12.4)") "fail residual = ", &
         maxval(abs(res(1:a%n,1)))
      errors = errors + 1
   endif
end subroutine chk_answer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_mat(a,extra)
   type(zd11_type), intent(inout) :: a
   integer, optional, intent(in) :: extra

   integer :: myextra,st

   myextra = 0
   if(present(extra)) myextra = extra

   !
   ! Create the simple sparse matrix: store lower triangle
   !
   ! 10.0  2.0       3.0
   !  2.0 10.0
   !           10.0  4.0
   !  3.0       4.0 10.0
   !

   a%n = 4
   deallocate(a%ptr, a%row, a%val, stat=st)
   allocate(a%ptr(a%n+1))
   a%ptr = (/ 1, 4, 5, 7, 8 /)
   allocate(a%row(2*(a%ptr(a%n+1)-1)+myextra))
   allocate(a%val(2*(a%ptr(a%n+1)-1)+myextra))

   a%row(1:7) = (/ 1, 2, 4,     2,    3, 4,    4 /)
   a%val(1:7) = (/   10.0, 2.0, 3.0, &
                     10.0, &
                     10.0, 4.0, &
                     10.0 /)

end subroutine simple_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_mat_zero_diag(a,extra)
   type(zd11_type), intent(inout) :: a
   integer, optional, intent(in) :: extra

   integer :: myextra,st

   myextra = 0
   if(present(extra)) myextra = extra

   !
   ! Create the simple sparse matrix:
   !
   !  0.0  1.0  2.0
   !  1.0  0.0  1.0
   !  2.0  1.0  0.0

   a%n = 3
   deallocate(a%ptr, a%row, a%val, stat=st)
   allocate(a%ptr(a%n+1))
   a%ptr = (/ 1, 3, 4, 4 /)
   allocate(a%row(2*(a%ptr(a%n+1)-1)+a%n+myextra))
   allocate(a%val(2*(a%ptr(a%n+1)-1)+a%n+myextra))

   a%row(1:3) = (/ 2, 3, 3 /)
   a%val(1:3) = (/   1.0, 2.0, &
                     1.0 /)

end subroutine simple_mat_zero_diag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_sing_mat(a)
   type(zd11_type), intent(inout) :: a
   integer :: st

   !
   ! Create the simple singular sparse matrix:
   !
   !  0.0  2.0 
   !  2.0  0.0  1.0
   !       1.0  0.0
   !
   ! we will not enter diagonal entries explicitly

   a%n = 3
   deallocate(a%ptr, a%row, a%val, stat=st)
   allocate(a%ptr(a%n+1))
   a%ptr = (/ 1, 2, 3, 3 /)
   allocate(a%row(9))
   allocate(a%val(9))
   a%row(1:2) = (/ 2, 3 /)
   a%val(1:2) = (/   2.0, 1.0 /)

end subroutine simple_sing_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_result(actual, expected)
   integer :: actual
   integer :: expected

   if(actual.eq.expected) then
      write(*,"(a)") "ok"
      return
   endif

   write(*,"(a)") "fail"
   write(*,"(2(a,i4))") "returned ", actual, ", expected ", expected
   errors = errors + 1
end subroutine print_result

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine compute_resid(nrhs,a,x,lx,rhs,lrhs,res)

   integer, intent(in) :: nrhs, lrhs, lx
   type(zd11_type), intent(in) :: a
   real(wp), intent(in) :: rhs(lrhs,nrhs)
   real(wp), intent(in) :: x(lx,nrhs)
   real(wp), dimension(:,:), allocatable, intent(inout) :: res
   real(wp), dimension(:), allocatable :: work

   integer :: i, j, k
   real(wp) :: anorm, atemp, bnorm(1:nrhs), xnorm(1:nrhs)

   if (allocated(res)) deallocate(res)
   allocate(res(a%n,nrhs),work(a%n))

       anorm = zero
       bnorm = zero
       xnorm = zero
       work = zero

      ! Check residual
       res(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
       do k = 1, a%n
          do j = a%ptr(k), a%ptr(k+1)-1
             i = a%row(j)
             atemp = a%val(j)
             res(i, 1:nrhs) = res(i, 1:nrhs) - atemp*x(k,1:nrhs)
             work(i) = work(i) + abs(atemp)
             if(i.eq.k) cycle
             res(k, 1:nrhs) = res(k, 1:nrhs) - atemp*x(i,1:nrhs)
             work(k) = work(k) + abs(atemp)
          end do
       end do

       do k = 1, a%n
          anorm = max(anorm,work(k))
          do i = 1,nrhs
             bnorm(i) = max(bnorm(i),abs(rhs(k,i)))
             xnorm(i) = max(xnorm(i),abs(x(k,i)))
          end do
       end do

       do k = 1,a%n
          do i = 1,nrhs
             res(k,i) = res(k,i)/(anorm*xnorm(i) + bnorm(i))
          end do
       end do

end subroutine compute_resid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program
