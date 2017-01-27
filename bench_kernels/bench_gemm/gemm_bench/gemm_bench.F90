program gemm_bench

  use omp_lib

  implicit none

  integer, parameter :: al = 16, seed(2) = (/ 0, 1072693248 /)

#ifdef INT_64
  integer, parameter :: ip = 8
#else
  integer, parameter :: ip = 4
#endif
  integer, parameter :: izero = 0_ip

#ifdef GEMM_S
  integer, parameter :: rp = 4, ic = 1
  real(rp), parameter :: zero = 0.0e0_rp
  character(len=5), parameter :: routine = 'SGEMM'
  real(rp), allocatable :: A(:,:), B(:,:), C(:,:), D(:,:)
  real(rp) :: alpha, beta
  external :: SGEMM
#endif
#ifdef GEMM_D
  integer, parameter :: rp = 8, ic = 1
  real(rp), parameter :: zero = 0.0e0_rp
  character(len=5), parameter :: routine = 'DGEMM'
  real(rp), allocatable :: A(:,:), B(:,:), C(:,:), D(:,:)
  real(rp) :: alpha, beta
  external :: DGEMM
#endif
#ifdef GEMM_C
  integer, parameter :: rp = 4, ic = 2
  complex(rp), parameter :: zero = (0.0e0_rp,0.0e0_rp)
  character(len=5), parameter :: routine = 'CGEMM'
  complex(rp), allocatable :: A(:,:), B(:,:), C(:,:), D(:,:)
  complex(rp) :: alpha, beta
  external :: CGEMM
#endif
#ifdef GEMM_Z
  integer, parameter :: rp = 8, ic = 2
  complex(rp), parameter :: zero = (0.0e0_rp,0.0e0_rp)
  character(len=5), parameter :: routine = 'ZGEMM'
  complex(rp), allocatable :: A(:,:), B(:,:), C(:,:), D(:,:)
  complex(rp) :: alpha, beta
  external :: ZGEMM
#endif

  character :: tranA, tranB
  integer(ip) :: m, n, k, rA, kA, rB, kB, ldA, ldB, ldC, i, j
  integer :: argc, runs, length, status
  character(len=53) :: argv ! (25,25)
  real(rp) :: harvest(2)
  double precision :: time_start, time_stop, times(4)

  integer, intrinsic :: command_argument_count, iand, ior, size
  integer(ip), intrinsic :: max
  intrinsic :: get_command_argument, random_number, random_seed

  argc = command_argument_count()
  if ((argc .lt. 6) .or. (argc .gt. 8)) stop 'gemm_bench.exe #runs tranA tranB m n k [ alpha beta ]'

  call random_seed(size=length)
  if (length .ne. size(seed)) stop 'random_seed(size) <> size(seed)'
  call random_seed(put=seed)

  call get_command_argument(1, argv, length, status)
  if (status .ne. 0) stop '#runs'
  read (argv,*) runs

  call get_command_argument(2, tranA, length, status)
  if (status .ne. 0) stop 'transA'

  call get_command_argument(3, tranB, length, status)
  if (status .ne. 0) stop 'transB'

  call get_command_argument(4, argv, length, status)
  if (status .ne. 0) stop 'm'
  read (argv,*) m

  call get_command_argument(5, argv, length, status)
  if (status .ne. 0) stop 'n'
  read (argv,*) n

  call get_command_argument(6, argv, length, status)
  if (status .ne. 0) stop 'k'
  read (argv,*) k

  if (argc .ge. 7) then
     call get_command_argument(7, argv, length, status)
     if (status .ne. 0) stop 'alpha'
     read (argv,*) alpha
  else if (ic .eq. 1) then
     harvest(1) = zero
     do while (harvest(1) .eq. zero)
        call random_number(harvest(1))
     end do
     alpha = harvest(1)
  else
     harvest = 0.0e0_rp
     do while (harvest(1) .eq. 0.0e0_rp)
        call random_number(harvest(1))
     end do
     do while (harvest(2) .eq. 0.0e0_rp)
        call random_number(harvest(2))
     end do
     alpha = cmplx(harvest(1), harvest(2), rp)
  end if

  if (argc .ge. 8) then
     call get_command_argument(8, argv, length, status)
     if (status .ne. 0) stop 'beta'
     read (argv,*) beta
  else if (ic .eq. 1) then
     harvest(1) = zero
     do while (harvest(1) .eq. zero)
        call random_number(harvest(1))
     end do
     beta = harvest(1)
  else
     harvest = 0.0e0_rp
     do while (harvest(1) .eq. 0.0e0_rp)
        call random_number(harvest(1))
     end do
     do while (harvest(2) .eq. 0.0e0_rp)
        call random_number(harvest(2))
     end do
     beta = cmplx(harvest(1), harvest(2), rp)
  end if

  if (runs .lt. 0) stop '#runs < 0'

  if (tranA .eq. 'n') then
     tranA = 'N'
  else if (tranA .eq. 't') then
     tranA = 'T'
  else if (tranA .eq. 'c') then
     tranA = 'C'
  end if

  if (tranA .eq. 'N') then
     rA = m
     kA = k
  else if (tranA .eq. 'T') then
     rA = k
     kA = m
  else if (tranA .eq. 'C') then
     rA = k
     kA = m
  else
     stop 'tranA not in { N, T, C }'
  end if

  if (tranB .eq. 'n') then
     tranB = 'N'
  else if (tranB .eq. 't') then
     tranB = 'T'
  else if (tranB .eq. 'c') then
     tranB = 'C'
  end if

  if (tranB .eq. 'N') then
     rB = k
     kB = n
  else if (tranB .eq. 'T') then
     rB = n
     kB = k
  else if (tranB .eq. 'C') then
     rB = n
     kB = k
  else
     stop 'tranB not in { N, T, C }'
  end if

  if (m .lt. izero) stop 'm < 0'
  if (n .lt. izero) stop 'n < 0'
  if (k .lt. izero) stop 'k < 0'

  ldA = align_me(rA)
  ldB = align_me(rB)
  ldC = align_me(m)

  status = 0
  if (alpha .ne. zero) status = ior(status, 1)
  if (beta .ne. zero) status = ior(status, 2)

  time_start = omp_get_wtime()
  allocate(A(ldA, kA))
  allocate(B(ldB, kB))
  if (iand(status, 1) .ne. 0) then
     if (ic .eq. 1) then
        do j = 1, kA
           do i = 1, rA
              call random_number(harvest)
              A(i,j) = harvest(1)
           end do
        end do
        do j = 1, kB
           do i = 1, rB
              call random_number(harvest)
              A(i,j) = harvest(1)
           end do
        end do
     else
        do j = 1, kA
           do i = 1, rA
              call random_number(harvest)
              A(i,j) = cmplx(harvest(1), harvest(2), rp)
           end do
        end do
        do j = 1, kB
           do i = 1, rB
              call random_number(harvest)
              A(i,j) = cmplx(harvest(1), harvest(2), rp)
           end do
        end do
     end if
  end if
  allocate(C(ldC, n))
  if (iand(status, 2) .ne. 0) then
     if (ic .eq. 1) then
        do j = 1, n
           do i = 1, m
              call random_number(harvest)
              C(i,j) = harvest(1)
           end do
        end do
     else
        do j = 1, n
           do i = 1, m
              call random_number(harvest)
              C(i,j) = cmplx(harvest(1), harvest(2), rp)
           end do
        end do
     end if
  end if
  if (runs .gt. 1) then
     allocate(D(ldC, n))
     D = C
  end if
  times(1) = omp_get_wtime() - time_start

  do i = 1, runs
     time_start = omp_get_wtime()
#ifdef GEMM_S
     call SGEMM(tranA, tranB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC)
#endif
#ifdef GEMM_D
     call DGEMM(tranA, tranB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC)
#endif
#ifdef GEMM_C
     call CGEMM(tranA, tranB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC)
#endif
#ifdef GEMM_Z
     call ZGEMM(tranA, tranB, m, n, k, alpha, A, ldA, B, ldB, beta, C, ldC)
#endif
     time_stop = omp_get_wtime() - time_start
     if (i .eq. 1) then
        times(2) = time_stop
        times(3) = time_stop / runs
        times(4) = time_stop
     else
        times(2) = min(times(2), time_stop)
        times(3) = times(3) + time_stop / runs
        times(4) = max(times(4), time_stop)
     end if
     if (i .lt. runs) C = D
  end do

  write (*,9) routine, tranA, tranB, m, n, k, times(1), times(2), times(3), times(4)

  if (runs .gt. 1) deallocate(D)
  deallocate(C)
  deallocate(B)
  deallocate(A)

9 format(A,',',A,',',A,',',I7,',',I7,',',I7,',',F9.2,',',F9.2,',',F9.2,',',F9.2)

contains

  elemental integer function align_me(n)

    implicit none

    integer, intent(in) :: n
    integer :: p, q
    integer, intrinsic :: mod

    p = mod(al, rp)
    if (p .ne. 0) then
       align_me = n
    else
       q = al / rp
       p = mod(n, q)
       if (p .eq. 0) then
          align_me = n
       else
          align_me = n + (q - p)
       end if
    end if

  end function align_me

end program gemm_bench
