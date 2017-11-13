module spllt_solve_task_mod
contains
  
  !*************************************************  
  !
  ! Forward solve with block on diagoanl
  !
  subroutine spllt_solve_fwd_block_task(dblk, nrhs, rhs, ldr, xlocal, fdata)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    implicit none
    
    integer, intent(in) :: dblk ! Index of block on the diagonal
    integer, intent(in) :: nrhs ! Number of RHS
    integer, intent(in) :: ldr ! Leading dimension of RHS
    real(wp), intent(inout) :: rhs(ldr)
    real(wp), dimension(:), intent(inout) :: xlocal
    type(spllt_fdata_type), intent(inout) :: fdata
    
    ! Node info
    integer :: sa
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa
    integer :: bcol, dcol, col
    integer :: offset
    integer :: node
        
    ! Get block info
    node = fdata%bc(dblk)%node
    m = fdata%bc(dblk)%blkm
    n = fdata%bc(dblk)%blkn
    sa = fdata%bc(dblk)%sa
    bcol = fdata%bc(dblk)%bcol ! Current block column
    dcol     = bcol - fdata%bc(fdata%nodes(node)%blk_sa)%bcol + 1
    col      = fdata%nodes(node)%sa + (dcol-1)*fdata%nodes(node)%nb
    offset   = col - fdata%nodes(node)%sa + 1
    
    ! Perform triangular solve
    call slv_solve(n, n, col, fdata%lfact(bcol)%lcol(sa:sa+n*n-1), &
         'Transpose    ', 'Non-unit', nrhs, rhs, ldr)
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m.gt.0) then
       sa = sa + n*n
       call slv_fwd_update(m, n, col, offset, fdata%nodes(node)%index, &
            fdata%lfact(bcol)%lcol(sa:sa+n*m-1), n, nrhs, &
            rhs, ldr, rhs, ldr, xlocal)
    endif

  end subroutine spllt_solve_fwd_block_task

  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         
  subroutine spllt_solve_bwd_block_task(dblk, nrhs, rhs, ldr, xlocal, fdata)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    implicit none

    integer, intent(in) :: dblk ! Index of block on the diagonal
    integer, intent(in) :: nrhs ! Number of RHS
    integer, intent(in) :: ldr ! Leading dimension of RHS
    real(wp), intent(inout) :: rhs(ldr)
    real(wp), dimension(:), intent(inout) :: xlocal
    type(spllt_fdata_type), intent(inout) :: fdata
    
    ! Node info
    integer :: sa
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa
    integer :: bcol, dcol, col
    integer :: offset
    integer :: node

    node = fdata%bc(dblk)%node

    ! print *, "[spllt_solve_bwd_block_task] node = ", node

    ! Get block info
    n      = fdata%bc(dblk)%blkn
    m      = fdata%bc(dblk)%blkm
    sa = fdata%bc(dblk)%sa
    bcol   = fdata%bc(dblk)%bcol ! Current block column
    col    = calc_col(fdata%nodes(node), fdata%bc(dblk)) ! current bcol
    col    = fdata%nodes(node)%sa + (col-1)*fdata%nodes(node)%nb
    offset = col - fdata%nodes(node)%sa + 1

    ! print *, "m = ", m, ", n = ", n
    ! print *, "blk_sa = ", blk_sa
    ! print *, "bcol = ", bcol
    ! print *, "col = ", col

    ! Perform and retangular update from diagonal block
    if(m.gt.n) then
       call slv_bwd_update(m-n, n, col, offset+n, fdata%nodes(node)%index, &
            fdata%lfact(bcol)%lcol(sa+n*n:sa+n*m-1), n, nrhs, rhs, &
            rhs, ldr, xlocal)
    endif

    ! Perform triangular solve
    call slv_solve(n, n, col, fdata%lfact(bcol)%lcol(sa:sa+n*n-1), &
         'Non-Transpose', 'Non-unit', nrhs, rhs, ldr)


  end subroutine spllt_solve_bwd_block_task

  !*************************************************
  !
  ! This function calculates column of a node we are on  
  integer function calc_col(node, bc)
    use spllt_data_mod    
    implicit none

    type(spllt_node), intent(in) :: node
    type(spllt_bc_type), intent(in) :: bc

    calc_col = (size(node%index)-1)/node%nb + 1 ! no. row blks for node

    calc_col = calc_col - (bc%last_blk - bc%dblk + 1) + 1 ! column of node

  end function calc_col

end module spllt_solve_task_mod
