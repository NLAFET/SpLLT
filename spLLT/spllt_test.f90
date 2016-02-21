program spllt_test
  use hsl_ma87_double
  use spllt_mod
  implicit none

  integer n, nnz
  
  write(*,'("[spllt test]")')

  n   = 256
  nnz = 10000

  call spllt_test_rand(n, nnz)

  stop

contains

  subroutine spllt_test_rand(n, nnz)
    use spllt_stf_mod
    implicit none
    
    integer :: n, nnz

    ! ma87
    type(ma87_keep) :: keep
    type(ma87_control) :: control
    type(ma87_info) :: info


    ! type(matrix_type) :: a
    integer :: iseed
    type(zd11_type) :: a
    real(wp), dimension(:), allocatable :: b, x
    integer :: i, nrhs
    integer, dimension(:), allocatable :: order
    real(wp) :: num_flops, resid

    a%m = n
    a%n = n
    a%ne = nnz
    
    allocate(a%ptr(n+1))
    allocate(a%row(2*a%ne), a%col(2*a%ne), a%val(2*a%ne))

    write(*,'("[spllt test rand] generate random matrix")')

    control%nb = 32
    ! control%nb = 64
    ! control%nb = 80
    ! control%nb = 120
    ! control%nb = 256

    call fa14id(iseed)

    ! generate matrix
    ! with spral
!!    call random_matrix_generate(state, SPRAL_MATRIX_REAL_SYM_PSDEF, & 
!!         & a%n, a%n, nnz, &
!!         & a%ptr, a%row, &
!!         & flag, &
!!         & val=a%val)

    call gen_random_posdef(a, a%ne, iseed)    

    write(*,'("[>] generate ordering")')

    ! ordering
    allocate(order(a%n))

    ! natural order
    do i = 1,a%n
       order(i) = i
    end do

    write(*,'("[>] analyse")')

    ! analysis
    call MA87_analyse(a%n, a%ptr, a%row, order, keep, control, info)
    num_flops = info%num_flops
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "error detected during analysis"
       goto 9999
    endif

    write(*,'("[>] [analysis] num flops: ", es10.3)') num_flops    
    write(*,'("[>] [analysis] num nodes: ", i10)') info%num_nodes    

    call spllt_print_atree(keep)
    
    ! generate rhs
    nrhs = 1
    allocate(b(a%m), x(a%n))
    call random_number(b)

    write(*,'("[>] factorize")')
    ! factorize matrix
    call spllt_stf_factorize(a%n, a%ptr, a%row, a%val, order, keep, control, info)
    ! call MA87_factor(a%n, a%ptr, a%row, a%val, order, keep, control, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "fail on factor"
    endif

! goto 9999 ! no solve

    write(*,'("[>] solve")')

    x = b
    ! solve
    call MA87_solve(x, order, keep, control, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a,i4)") " fail on 1d solve", &
            info%flag
    endif
    
    call spllt_bwerr(a, x, b, resid)
    
    write(*,'("[>] [solve] bwderr ||Ax-b|| / (||A||||x|| + ||b||): ", es10.3)') resid    

    call MA87_finalise(keep, control)
    
9999 continue

    return
  end subroutine spllt_test_rand
  
end program spllt_test
