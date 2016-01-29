program spllt_test
  use hsl_ma87_double
  implicit none

  integer, parameter :: wp = kind(0d0)

  integer n, nnz

  type :: matrix_type
     integer :: n, ne
     integer, dimension(:), allocatable :: ptr, row, col
     real(wp), dimension(:), allocatable :: val
  end type matrix_type
  
  write(*,*) "[spllt test]"

  n   = 100
  nnz = 1000

  call spllt_test_rand(n, nnz)

  stop

contains

  subroutine spllt_test_rand(n, nnz)
    use spral_random
    use spral_random_matrix, only : random_matrix_generate
    use spral_matrix_util, only : SPRAL_MATRIX_REAL_SYM_PSDEF
    implicit none
    
    integer :: n, nnz

    ! ma87
    type(ma87_keep) :: keep
    type(ma87_control) :: control
    type(ma87_info) :: info

    type(random_state) :: state
    type(matrix_type) :: a
    integer :: flag
    integer :: i
    integer, dimension(:), allocatable :: order
    real(wp) :: num_flops


    a%n = n
    a%ne = nnz
    
    allocate(a%ptr(n+1))
    allocate(a%row(a%ne), a%col(a%ne), a%val(a%ne))

    write(*,*) "[spllt test rand] generate random matrix"

    ! generate matrix
    ! with spral
    call random_matrix_generate(state, SPRAL_MATRIX_REAL_SYM_PSDEF, & 
         & a%n, a%n, nnz, &
         & a%ptr, a%row, &
         & flag, &
         & val=a%val)

    ! ordering
    allocate(order(a%n))

    ! natural order
    do i = 1,a%n
       order(i) = i
    end do

    ! analysis
    call MA87_analyse(a%n, a%ptr, a%row, order, keep, control, info)
    num_flops = info%num_flops
    if(info%flag .lt. MA87_SUCCESS) then
       write(*, "(a)") "error detected during analysis"
       goto 9999
    endif

9999 continue

    return
  end subroutine spllt_test_rand

end program spllt_test
