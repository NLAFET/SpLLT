program spllt_test
  use hsl_ma87_double
  use spllt_mod
  implicit none

  type(spllt_cntl) :: cntl
  integer n, nnz
  
  write(*,'("[spllt test]")')

  cntl%ncpu = 4
  cntl%nb   = 8
  
  ! n   = 256
  ! nnz = 1000

  ! n   = 256
  ! nnz = 1000

  ! call spllt_test_rand(n, nnz, cntl)

  call spllt_test_mat("mesh2e1.rb", cntl)
  
  stop

contains

  subroutine spllt_test_mat(matfile, cntl)
    use spllt_stf_mod
    use spral_rutherford_boeing
    implicit none
    
    character(len=*), intent(in) :: matfile
    type(spllt_cntl) :: cntl

    type(rb_reader_options) :: rb_options
    integer :: flag
    integer :: m, n

    type(zd11_type) :: a
    integer :: i, nrhs
    integer, dimension(:), allocatable :: order
    real(wp) :: num_flops, resid
    real(wp), dimension(:), allocatable :: b, x
    type(spllt_data_type) :: pbl

    ! ma87
    type(ma87_keep) :: keep
    type(ma87_control) :: control
    type(ma87_info) :: info

    write(*,'("[spllt test mat] read matrix")')
    
    rb_options%values = 2 ! Make up values if necessary
    rb_options%format = 1 ! Coordinate
    call rb_read(matfile, a%m, a%n, a%ptr, a%row, a%col, &
         a%val, rb_options, flag)

    m = a%m
    n = a%n

    write(*,'("[>] generate ordering")')

    ! ordering
    allocate(order(a%n))

    ! amd ordering
    call amd_order(a, order)

    ! natural order
    ! do i = 1,n
    !    order(i) = i
    ! end do

    write(*,'("[>] perform analysis")')

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
    allocate(b(m), x(n))
    call random_number(b)

    write(*,'("[>] factorize")')
    control%nb = cntl%nb
    ! factorize matrix
    call spllt_stf_factorize(a%n, a%ptr, a%row, a%val, order, keep, control, info, pbl, cntl)
    ! call MA87_factor(a%n, a%ptr, a%row, a%val, order, keep, control, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "failed factorization"
    endif

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
  end subroutine spllt_test_mat

  subroutine spllt_test_rand(n, nnz, cntl)
    use spllt_stf_mod
    implicit none
    
    integer :: n, nnz
    type(spllt_cntl) :: cntl

    type(spllt_data_type) :: pbl
    
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

    ! control%nb = 10
    control%nb = 8
    ! control%nb = 39
    ! control%nb = 64
    ! control%nb = 119
    ! control%nb = 120
    ! control%nb = 256

    ! control%nemin = 10

    ! control%nb = 40

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
    call spllt_stf_factorize(a%n, a%ptr, a%row, a%val, order, keep, control, info, pbl, cntl)
    ! call MA87_factor(a%n, a%ptr, a%row, a%val, order, keep, control, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "failed factorization"
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
