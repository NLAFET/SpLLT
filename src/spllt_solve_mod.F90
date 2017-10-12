module spllt_solve_mod
  use hsl_ma87_double

   interface spllt_solve
      module procedure spllt_solve_one_double
   end interface

contains
  
  !*************************************************
  !
  !
  ! Solve phase. simplified interface for a single rhs
  !
  subroutine spllt_solve_one_double(x, order, keep, cntl, info, job)
    use spllt_data_mod
    ! use hsl_ma87_double
    implicit none

    type(spllt_keep), intent(inout) :: keep
    real(wp), intent(inout) :: x(keep%n) ! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i) is the corresponding component of the right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i) holds solution for variable i.
    integer, intent(in) :: order(:) ! pivot order. must be unchanged
    ! For details of keep, control, info : see derived type description
    type(spllt_cntl), intent(in) :: cntl
    type(spllt_info), intent(out) :: info
    integer, optional, intent(in) :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)

    integer :: j ! Iterator
    integer :: n ! Order of the system/First dimension of RHS
    integer :: nrhs = 1 ! Number of RHS
    real(wp), dimension(:), allocatable :: soln ! allocated to have
    ! size n*nrhs.  used to hold reordered rhs and then overwritten by
    ! reorder solution.
    integer :: st ! stat parameter

    n = keep%n

    ! immediate return if n = 0
    if (n == 0) return

    ! ma87
    ! type(ma87_control) :: ma_control
    ! type(ma87_keep) :: ma_keep
    ! type(ma87_info) :: ma_info
    
    ! Set control for HSL_MA87
    ! ma_control%nb = cntl%nb

    ! Set keep for HSL_MA87
    ! ma_keep%n = keep%n
    ! allocate(keep%nodes(-1:info%num_nodes+1))

    ! ma_keep = keep

    ! call spllt_solve_mult_double(1, keep%n, x, order, keep, &
    !      control, info, job)
    ! TODO Use ssids solve ?
    ! call MA87_solve(x, order, ma_keep, ma_control, ma_info, job)
    ! call MA87_solve(nrhs, n, soln, order, ma_keep, ma_control, ma_info)

    !
    ! Reorder rhs
    !
    deallocate(soln,stat=st)
    allocate(soln(n*nrhs),stat=st)
    
    ! do i = 1, nrhs
    !    do j = 1, n
    !       soln((i-1)*n + order(j)) = x(j, i)
    !    end do
    ! end do

    do j = 1, n
       soln(order(j)) = x(j)
    end do

    ! Forward solve
    !
    call solve_fwd(nrhs, soln, n, keep)
    
    ! Backward solve
    call solve_bwd(nrhs, soln, n, keep)

   !
   ! Reorder soln
   !
    do j = 1, n
       x(j) = soln(order(j))
    end do

  end subroutine spllt_solve_one_double

  !*************************************************
  !
  ! Forward solve routine
  subroutine solve_fwd(nrhs, rhs, ldr, keep)
    use spllt_data_mod
    use spllt_solve_task_mod
    use spllt_solve_kernels_mod
    implicit none

    type(spllt_keep), intent(inout) :: keep
    integer, intent(in) :: nrhs ! Number of RHS
    integer, intent(in) :: ldr ! Leading dimension of RHS
    real(wp), intent(inout) :: rhs(ldr)

    ! real(wp) :: xlocal(keep%n)
    integer :: num_node
    integer :: node
    integer :: i, j
    integer :: sa, en, blk_sa
    integer :: numcol, numrow ! Number of column/row in node 
    integer :: nc, nr ! Number of block-column/block-row in node
    integer :: jj, ii
    integer :: bcol, dcol, col, offset
    integer :: dblk ! Diagonal index 
    integer :: s_nb ! Block size in node
    integer :: m, n ! Block dimension 
    integer :: blk ! Block index
    real(wp), dimension(:), allocatable :: xlocal ! update_buffer workspace
    integer :: st ! Stat parameter

    print *, "[spllt_solve_mod] solve_fwd"

    ! Allocate workspace
    allocate(xlocal(keep%maxmn*nrhs), stat=st)

    num_node = keep%info%num_nodes

    do node = 1, num_node

       ! Get node info
       s_nb = keep%nodes(node)%nb
       sa = keep%nodes(node)%sa
       en = keep%nodes(node)%en
       numcol = en - sa + 1
       numrow = size(keep%nodes(node)%index)
       nc = (numcol-1) / s_nb + 1
       nr = (numrow-1) / s_nb + 1 
       
       ! Get first diag block in node
       dblk = keep%nodes(node)%blk_sa
       ! Loop over block columns
       do jj = 1, nc
          
          !
          ! Forward solve with block on diagoanl
          !
          call spllt_solve_fwd_block_task(dblk, nrhs, rhs, ldr, xlocal, keep)
          
          do ii = jj+1, nr
             blk = dblk+ii-jj

             !
             ! Forward update with off-diagonal
             !

             ! Establish variables describing block
             n        = keep%blocks(blk)%blkn
             m        = keep%blocks(blk)%blkm
             blk_sa       = keep%blocks(blk)%sa
             bcol     = keep%blocks(blk)%bcol
             dcol     = bcol - keep%blocks(keep%nodes(node)%blk_sa)%bcol + 1
             col      = keep%nodes(node)%sa + (dcol-1)*keep%nodes(node)%nb

             offset   = col - keep%nodes(node)%sa + 1 ! diagonal blk
             offset   = offset + (blk-keep%blocks(blk)%dblk) * keep%nodes(node)%nb ! this blk
             
             call slv_fwd_update(m, n, col, offset, keep%nodes(node)%index, &
                  keep%lfact(bcol)%lcol(blk_sa:blk_sa+n*m-1), n, nrhs, rhs, &
                  ldr, rhs, ldr, xlocal)
             
          end do
          
          ! Update diag block in node          
          dblk = keep%blocks(dblk)%last_blk + 1
       end do
              
    end do

    ! Deallocate workspace
    deallocate(xlocal)

  end subroutine solve_fwd

  subroutine solve_bwd(nrhs, rhs, ldr, keep)
    use spllt_data_mod
    use spllt_solve_task_mod
    use spllt_solve_kernels_mod
    implicit none

    type(spllt_keep), intent(inout) :: keep
    integer, intent(in) :: nrhs ! Number of RHS
    integer, intent(in) :: ldr ! Leading dimension of RHS
    real(wp), intent(inout) :: rhs(ldr)

    integer :: num_node
    ! Node info
    integer :: node
    integer :: sa, en
    integer :: numcol, numrow ! Number of column/row in node 
    integer :: nc, nr ! Number of block-column/block-row in node
    integer :: s_nb ! Block size in node
    integer :: jj, ii
    integer :: dblk
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa 
    integer :: bcol ! Global block-column index
    integer :: dcol, col, offset
    integer :: blk ! Block index
    real(wp), dimension(:), allocatable :: xlocal ! update_buffer workspace
    integer :: st ! Stat parameter
    
    print *, "[spllt_solve_mod] solve_bwd"

    ! Allocate workspace
    allocate(xlocal(keep%maxmn*nrhs), stat=st)

    num_node = keep%info%num_nodes
    
    do node = num_node, 1, -1

       ! Get node info
       s_nb = keep%nodes(node)%nb
       sa = keep%nodes(node)%sa
       en = keep%nodes(node)%en
       numcol = en - sa + 1
       numrow = size(keep%nodes(node)%index)
       nc = (numcol-1) / s_nb + 1
       nr = (numrow-1) / s_nb + 1 

       ! Get first diag block in node
       dblk = keep%blocks(keep%nodes(node)%blk_en)%dblk

       ! Loop over block columns
       do jj = nc, 1, -1

          do ii = nr, jj+1, -1
             
             blk = dblk+ii-jj ! Block index

             !
             ! Backward update with block on diagoanl
             !

             ! Establish variables describing block
             n        = keep%blocks(blk)%blkn
             m        = keep%blocks(blk)%blkm
             blk_sa       = keep%blocks(blk)%sa
             bcol     = keep%blocks(blk)%bcol
             ! node     = keep%blocks(blk)%node
             dcol     = bcol - keep%blocks(keep%nodes(node)%blk_sa)%bcol + 1
             col      = keep%nodes(node)%sa + (dcol-1)*keep%nodes(node)%nb

             offset   = col - keep%nodes(node)%sa + 1 ! diagonal blk
             offset   = offset + (blk-keep%blocks(blk)%dblk) * keep%nodes(node)%nb ! this blk

             call slv_bwd_update(m, n, col, offset, keep%nodes(node)%index, &
                  keep%lfact(bcol)%lcol(blk_sa:blk_sa+n*m-1), n, nrhs, rhs, &
                  rhs, ldr, xlocal)
             
          end do

          !
          ! Backward solve with block on diagoanl
          !
          call spllt_solve_bwd_block_task(dblk, nrhs, rhs, ldr, xlocal, keep)
          
          ! Update diag block in node          
          if (jj .gt. 1) dblk = keep%blocks(dblk-1)%dblk          
       end do
       
    end do

    ! Deallocate workspace
    deallocate(xlocal)

  end subroutine solve_bwd
  
end module spllt_solve_mod
