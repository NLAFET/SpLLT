#include "pastix_fortran.h"

program run_prob
   use hsl_mc56_double
   use hsl_zd11_double
   use hsl_fa14_double
   ! use mpi
!$ use omp_lib
   implicit none

   integer, parameter :: dp = kind(0d0)

   integer :: matrix_type
   character(len=200) :: matfile    
   type(mc56_control) :: control56
   integer :: info56

   type(ZD11_type) :: matrix

   type(fa14_seed) :: seed

   double precision, dimension(:, :), allocatable :: rhs, soln
   double precision, dimension(:), allocatable :: work, res, scaling
   integer, dimension(:), allocatable :: order, invp
   integer :: icntl57(20), info57(40), icntl77(10)
   double precision :: rinfo57(20), cntl77(10), cntl57(5)
   integer, allocatable :: iwork(:), keep57(:)

   integer :: i, j, k

   integer :: start_t, stop_t, rate_t

   pastix_data_ptr_t :: pastix_data
   integer :: pastix_comm
   pastix_int_t :: icontrol(64)
   double precision :: dcontrol(64)
   integer :: nrhs
   integer :: required, provided, StatInfo

   interface
      subroutine pastix_fortran(pastix_data, pastix_comm, n, ia, ja, avals, &
            perm, invp, b, rhs, iparm, dparm)
         pastix_data_ptr_t :: pastix_data
         integer :: pastix_comm
         pastix_int_t :: n, rhs, ia(*), ja(*)
         pastix_float_t :: avals(*), b(*)
         pastix_int_t :: perm(*), invp(*), iparm(64)
         real*8 :: dparm(64)
      end subroutine pastix_fortran
   end interface

   ! pastix_comm = 0
   ! required=MPI_THREAD_MULTIPLE
   ! call MPI_Init_thread(required,provided,StatInfo)
   ! pastix_comm = MPI_COMM_WORLD

   ! Read in a matrix
   matfile = 'matrix.rb'
      
   ! matrix_type = HSL_MATRIX_REAL_SYM_PSDEF

   write(*, "(a)") "Reading..."
   ! select case(matrix_type)
   ! case(HSL_MATRIX_REAL_SYM_PSDEF)
   control56%values = 3 ! make up values if necessary (posdef)
   ! case(HSL_MATRIX_REAL_SYM_INDEF)
      ! control56%values = 2 ! make up values if necessary (indef)
   ! end select
   call mc56_read(matfile, matrix, control56, info56)
   if(info56.ne.0) then
      print *, "mc56 failed with error ", info56
      stop
   endif
   write(*, "(a)") "ok"

   ! write(*, "(a)") "Reading..."
   ! call MC56_init(file_control)
   ! file_control%iunit = 11
   ! open(unit=file_control%iunit, file="matrix.rb", status="old")
   ! call MC56_read(file_control, dattyp, matrix, seed, finfo)
   ! close(unit=file_control%iunit)
   ! write(*, "(a)") "ok"

   ! Use MA57 for ordering
   write(*, "(a)") "Ordering..."
   allocate(matrix%col(matrix%ne*2))
   do i = 1, matrix%n
      matrix%col(matrix%ptr(i):matrix%ptr(i+1)-1) = i
   end do
   allocate(order(matrix%n), invp(matrix%n))
   allocate(iwork(5*matrix%n), keep57(7*matrix%n+2*matrix%ne+42))
   call ma57id(cntl57, icntl57)
   !icntl57(6) = 2 ! Use MC47
   icntl57(6) = 4 ! Use MeTiS
   call ma57ad(matrix%n, matrix%ne, matrix%row, matrix%col, size(keep57), &
      keep57, iwork, icntl57, info57, rinfo57)
   order = keep57(1:matrix%n)
   deallocate(iwork,keep57,matrix%col)
   select case(info57(36))
   case(0)
      print *, "Used MC47, no dense rows"
   case(2)
      print *, "Used MC47"
   case(4)
      print *, "Used MeTiS"
   case default
      print *, "Unknown ordering ", info57(36)
   end select
   !print *, "perm = ", order

   do i = 1, matrix%n
      invp(order(i)) = i
   end do

   !! Swap order and invp meanings
   !invp(:) = order(:)
   !do i = 1, matrix%n
   !   order(invp(i)) = i
   !end do

   !! Use MC77_inf for scaling
   !write(*, "(a)") "Scaling..."
   !call mc77id(icntl77, cntl77)
   !icntl77(4) = 1 ! No data checking
   !icntl77(6) = 1 ! Symmetric
   !allocate(iwork(2*matrix%n), work(matrix%ne+2*matrix%n))
   !call mc77ad(1, matrix%n, matrix%n, matrix%ne, matrix%ptr, matrix%row, &
   !   matrix%val, iwork, size(iwork), work, size(work), icntl77, cntl77, &
   !   info57, rinfo57)
   !deallocate(iwork)
   !allocate(scaling(matrix%n))
   !scaling = 1/work(1:matrix%n)
   !deallocate(work)
   !write(*, "(a)") "ok"

   ! Make up a rhs
   allocate(rhs(matrix%n, 1), soln(matrix%n, 1))
   rhs = 0
   do i = 1, matrix%n
      do j = matrix%ptr(i), matrix%ptr(i+1)-1
         k = matrix%row(j)
         rhs(k, 1) = rhs(k, 1) + matrix%val(j)
         if(k.ne.i) rhs(i, 1) = rhs(i, 1) + matrix%val(j)
      end do
   end do

   !
   ! Initialize PaStiX parameters
   !
   pastix_data = 0
   nrhs        = 1
   icontrol(IPARM_MODIFY_PARAMETER)    = API_NO
   icontrol(IPARM_START_TASK)          = API_TASK_INIT
   icontrol(IPARM_END_TASK)            = API_TASK_INIT

   call pastix_fortran(pastix_data, pastix_comm, matrix%n, matrix%ptr, &
      matrix%row, matrix%val, order, invp, rhs, nrhs, icontrol, dcontrol)


   write(*, "(a)") "Analyse..."
   ! Analyse and factor
   icontrol(IPARM_START_TASK)    = API_TASK_ORDERING
   icontrol(IPARM_END_TASK)      = API_TASK_ANALYSE
   !icontrol(IPARM_MATRIX_VERIFICATION) = API_YES
   !icontrol(IPARM_ORDERING)      = API_ORDER_PERSONAL
   icontrol(IPARM_FACTORIZATION) = API_FACT_LLT
   !icontrol(IPARM_ESP)           = API_YES
   !icontrol(IPARM_VERBOSE)       = API_VERBOSE_YES
   icontrol(IPARM_THREAD_NBR)    = 1
!$ icontrol(IPARM_THREAD_NBR)    = omp_get_max_threads()
   call pastix_fortran(pastix_data, pastix_comm, matrix%n, matrix%ptr, &
      matrix%row, matrix%val, order, invp, rhs, nrhs, icontrol, dcontrol)
   write(*, "(a)") "ok"

   write(*, "(a)") "Factorize..."
   icontrol(IPARM_START_TASK)    = API_TASK_NUMFACT
   icontrol(IPARM_END_TASK)      = API_TASK_NUMFACT
   call system_clock(start_t, rate_t)
   call pastix_fortran(pastix_data, pastix_comm, matrix%n, matrix%ptr, &
      matrix%row, matrix%val, order, invp, rhs, nrhs, icontrol, dcontrol)
   call system_clock(stop_t)
   print *, "Factor took ", (stop_t - start_t)/real(rate_t)


   ! Solve
   write(*, "(a)") "Solve..."
   soln = rhs
   icontrol(IPARM_START_TASK)    = API_TASK_SOLVE
   icontrol(IPARM_END_TASK)      = API_TASK_SOLVE
   icontrol(IPARM_RHS_MAKING)    = API_RHS_B
   call system_clock(start_t, rate_t)
   call pastix_fortran(pastix_data, pastix_comm, matrix%n, matrix%ptr, &
      matrix%row, matrix%val, order, invp, soln, nrhs, icontrol, dcontrol)
   call system_clock(stop_t)
   print *, "Solve took ", (stop_t - start_t)/real(rate_t)

   print *, "fwd error || ||_inf = ", maxval(abs(soln(1:matrix%n,1)-1.0))
   allocate(res(1))
   call internal_calc_norm(matrix, soln, rhs, 1, res)
   print *, "bwd error scaled = ", res

   ! Finish
   icontrol(IPARM_START_TASK)    = API_TASK_CLEAN
   icontrol(IPARM_END_TASK)      = API_TASK_CLEAN
   call pastix_fortran(pastix_data, pastix_comm, matrix%n, matrix%ptr, &
      matrix%row, matrix%val, order, invp, rhs, nrhs, icontrol, dcontrol)

   ! call MPI_FINALIZE(StatInfo)

contains

   subroutine internal_calc_norm(matrix, x_vec, b_vec, nrhs, res)
      type(ZD11_type), intent(in) :: matrix
      integer, intent(in) :: nrhs
      real(dp), dimension(nrhs*matrix%n), intent(in) :: x_vec
      real(dp), dimension(nrhs*matrix%n), intent(in) :: b_vec
      real(dp), dimension(nrhs), intent(out) :: res

      integer :: i, j, k
      double precision, allocatable, dimension(:) :: x_norm
      real(dp), dimension(:), allocatable :: res_vec
      double precision :: temp
      double precision :: normA

      ! Find the residual
      allocate(res_vec(matrix%n*nrhs), x_norm(nrhs))
      res_vec = 0
      do i = 1, matrix%n
         do j = matrix%ptr(i), matrix%ptr(i+1)-1
            do k = 0, nrhs-1
               res_vec(i+k*matrix%n) = res_vec(i+k*matrix%n) + &
                  matrix%val(j) * x_vec(matrix%row(j)+k*matrix%n)
               if(i.ne.matrix%row(j)) res_vec(matrix%row(j)+k*matrix%n) = &
                  res_vec(matrix%row(j)+k*matrix%n) + &
                  matrix%val(j) * x_vec(matrix%row(j)+k*matrix%n)
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
      integer :: i, j, k

      allocate(row_norm(matrix%n))

      !FIXME:  This is wrong
      row_norm = 0
      do i = 1, matrix%n
         do j = matrix%ptr(i), matrix%ptr(i+1)-1
            k = matrix%row(j)
            row_norm(k) = row_norm(k) + abs(matrix%val(j))
            if(i.ne.k) row_norm(i) = row_norm(i) + abs(matrix%val(j))
         end do
      end do

      norm = maxval(row_norm) 
   end subroutine matrix_inf_norm
end program
