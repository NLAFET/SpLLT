module spllt_solve_task_mod
contains
  
  !*************************************************  
  ! Solve forward block task
  !
  subroutine spllt_solve_fwd_block_task(dblk, nrhs, rhs, ldr, xlocal, keep)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    implicit none
    
    integer, intent(in) :: dblk ! Index of block on the diagonal
    integer, intent(in) :: nrhs ! Number of RHS
    integer, intent(in) :: ldr ! Leading dimension of RHS
    real(wp), intent(inout) :: rhs(ldr)
    real(wp), dimension(:), intent(inout) :: xlocal
    type(spllt_keep), intent(inout) :: keep
    
    ! Node info
    integer :: sa
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa
    integer :: bcol, dcol, col
    integer :: offset
    integer :: node

        
    ! Get block info
    node = keep%blocks(dblk)%node
    sa = keep%nodes(node)%sa
    m = keep%blocks(dblk)%blkm
    n = keep%blocks(dblk)%blkn
    blk_sa = keep%blocks(dblk)%sa
    bcol = keep%blocks(dblk)%bcol ! Current block column
    dcol     = bcol - keep%blocks(keep%nodes(node)%blk_sa)%bcol + 1
    col      = keep%nodes(node)%sa + (dcol-1)*keep%nodes(node)%nb
    offset   = col - keep%nodes(node)%sa + 1
    
    !
    ! Forward solve with block on diagoanl
    !
    ! Perform triangular solve
    call slv_solve(n, n, col, keep%lfact(bcol)%lcol(sa:sa+n*n-1), &
         'Transpose    ', 'Non-unit', nrhs, rhs, ldr)
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m.gt.0) then
       sa = sa + n*n
       call slv_fwd_update(m, n, col, offset, keep%nodes(node)%index, &
            keep%lfact(bcol)%lcol(sa:sa+n*m-1), n, nrhs, &
            rhs, ldr, rhs, ldr, xlocal)
    endif

  end subroutine spllt_solve_fwd_block_task

  !*************************************************  
  ! Solve backward block task
  !
  subroutine spllt_solve_bwd_block_task(dblk, nrhs, rhs, ldr, xlocal, keep)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    implicit none

    integer, intent(in) :: dblk ! Index of block on the diagonal
    integer, intent(in) :: nrhs ! Number of RHS
    integer, intent(in) :: ldr ! Leading dimension of RHS
    real(wp), intent(inout) :: rhs(ldr)
    real(wp), dimension(:), intent(inout) :: xlocal
    type(spllt_keep), intent(inout) :: keep
    
    ! Node info
    integer :: sa
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa
    integer :: bcol, dcol, col
    integer :: offset
    integer :: node

    node = keep%blocks(dblk)%node
    ! Get node info
    sa = keep%nodes(node)%sa
    ! Get block info
    n      = keep%blocks(dblk)%blkn
    m      = keep%blocks(dblk)%blkm
    blk_sa = keep%blocks(dblk)%sa
    bcol   = keep%blocks(dblk)%bcol ! Current block column
    bcol   = keep%blocks(dblk)%bcol
    col    = calc_col(keep%nodes(node), keep%blocks(dblk)) ! current bcol
    col    = keep%nodes(node)%sa + (col-1)*keep%nodes(node)%nb
    offset = col - keep%nodes(node)%sa + 1

              !
    ! Backward solve with block on diagoanl
    !         
    
    ! Perform and retangular update from diagonal block
    if(m.gt.n) then
       call slv_bwd_update(m-n, n, col, offset+n, keep%nodes(node)%index, &
            keep%lfact(bcol)%lcol(sa+n*n:sa+n*m-1), n, nrhs, rhs, &
            rhs, ldr, xlocal)
    endif

    ! Perform triangular solve
    call slv_solve(n, n, col, keep%lfact(bcol)%lcol(sa:sa+n*n-1), &
         'Non-Transpose', 'Non-unit', nrhs, rhs, ldr)


  end subroutine spllt_solve_bwd_block_task

  !*************************************************
  !
  ! This function calculates column of a node we are on
  
  integer function calc_col(node, block)
    use spllt_data_mod    
    implicit none

    type(node_type), intent(in) :: node
    type(block_type), intent(in) :: block

    calc_col = (size(node%index)-1)/node%nb + 1 ! no. row blks for node

    calc_col = calc_col - (block%last_blk - block%dblk + 1) + 1 ! column of node

  end function calc_col

end module spllt_solve_task_mod
