module spllt_factorization_mod
  use spllt_mod
  implicit none
  
contains

  ! factorize block 
  ! _potrf
  subroutine spllt_factorize_block_task(bc, lfact, & 
       & control, flag, detlog)
    use hsl_ma87_double
    implicit none
    
    type(block_type), intent(inout) :: bc ! block to be factorized    
    type(lfactor), dimension(:), allocatable, intent(inout) :: lfact
    type(MA87_control), intent(in) :: control 
    integer, intent(out) :: flag ! Error flag
    real(wp), dimension(:), allocatable, intent(inout) :: detlog ! per thread sum of log pivot

    integer :: m, n, bcol, sa
    integer(long) :: id
    
    m    = bc%blkm
    n    = bc%blkn
    id   = bc%id 
    bcol = bc%bcol

    sa   = bc%sa

    ! factorize_block
    call factor_diag_block(n, m, id, &
         & lfact(bcol)%lcol(sa:sa+n*m-1),   &
         & control, flag, detlog(0))
    
    return
  end subroutine spllt_factorize_block_task

  ! _trsm
  subroutine spllt_solve_block_task(bc_kk, bc_ik, lfact, control)
    use hsl_ma87_double
    implicit none

    type(block_type), intent(inout) :: bc_kk, bc_ik ! block to be factorized    
    type(lfactor), dimension(:), allocatable, intent(inout) :: lfact
    type(MA87_control), intent(in) :: control     
    
    integer :: m, n, bcol, sa

    ! bc_kk
    
    ! bc_ik
    n  = bc_ik%blkn
    m  = bc_ik%blkm
    sa    = bc_ik%sa
    id    = bc_ik%id
    
    ! bcol is block column that blk and dblk belong to
    bcol = bc%bcol    

    ! solve_block task
    call solv_col_block(m, n, id, & 
         & lfact(bcol)%lcol(sa:sa+n*m-1), &
         & bc_kk%id, keep%lfact(bcol)%lcol(bc_kk%sa:bc_kk%sa+blkn*bc_kk%blkm), &
         & control)

    return
  end subroutine spllt_solve_block_task

end module spllt_factorization_mod
