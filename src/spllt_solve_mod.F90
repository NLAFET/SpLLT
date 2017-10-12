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
    integer :: numcol, nc, numrow
    integer :: kk
    integer :: bcol, dcol, col, offset
    integer :: dblk ! Diagonal index 
    integer :: s_nb ! Block size in node
    integer :: m, n ! Block dimension 
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
       
       ! Get first diag block in node
       dblk = keep%nodes(node)%blk_sa
       ! Loop over block columns
       do kk = 1, nc
          bcol = keep%blocks(dblk)%bcol ! Current block column
          
          ! Forward solve with diagonal block
          m = keep%blocks(dblk)%blkm
          n = keep%blocks(dblk)%blkn
          blk_sa = keep%blocks(dblk)%sa
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
          
          ! Update diag block in node          
          dblk = keep%blocks(dblk)%last_blk + 1
       end do
              
    end do

    ! Deallocate workspace
    deallocate(xlocal)

  end subroutine solve_fwd

  subroutine solve_bwd(nrhs, rhs, ldr, keep)
    use spllt_data_mod
    implicit none

    type(spllt_keep), intent(inout) :: keep
    integer, intent(in) :: nrhs ! Number of RHS
    integer, intent(in) :: ldr ! Leading dimension of RHS
    real(wp), intent(inout) :: rhs(ldr)

    integer :: num_node
    ! Node info
    integer :: node
    integer :: sa, en
    integer :: numcol, nc, numrow
    integer :: s_nb ! Block size in node
    integer :: kk
    integer :: dblk
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa 
    integer :: bcol ! Global block-column index
    integer :: col, offset
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

       ! Get first diag block in node
       dblk = keep%blocks(keep%nodes(node)%blk_en)%dblk

       ! Loop over block columns
       do kk = nc, 1, -1
          bcol = keep%blocks(dblk)%bcol ! Current block column

          ! Establish variables describing block column
          n      = keep%blocks(dblk)%blkn
          m      = keep%blocks(dblk)%blkm
          blk_sa = keep%blocks(dblk)%sa
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


          ! Update diag block in node          
          if (kk .gt. 1) dblk = keep%blocks(dblk-1)%dblk          
       end do
       
    end do

    ! Deallocate workspace
    deallocate(xlocal)

  end subroutine solve_bwd
  
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

  !*************************************************
  !
  ! TASK_BSOLV, TASK_FSOLV
  ! B_j <- L_jj^-1 B_j
  ! B_j <- L_jj^-T B_j
  !
  ! Note: While diagonal blocks may be trapezoidal, this is handled at the
  ! level calling this subroutine through a call to slv_fwd_update or
  ! slv_bwd_update
  
  subroutine slv_solve(n, nelim, col, dest, trans, unit,  &
       nrhs, rhs, ldr)
    use spllt_data_mod
    implicit none

    integer, intent(in) :: n ! leading dimension of diag. block
    integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
    integer, intent(in) :: col ! start of block column variables in rhs
    real(wp), dimension(*), intent(in) :: dest ! holds destination block
    character(len=13), intent(in) :: trans ! set to 
    ! 'Transpose    ' for forward substitution and to 
    ! 'Non-Transpose' for back substitution
    character(len=8), intent(in) :: unit ! set to 
    ! 'Non-Unit' for positive-definite case
    integer, intent(in) :: nrhs ! number of right-hand sides
    integer, intent(in) :: ldr ! leading extent of rhs
    real(wp), intent(inout) :: rhs(ldr*nrhs)

    !%%%   integer :: this_thread
    !%%%   integer :: t_start, t_end

    if(nelim.eq.0) return

    !%%%   if(control%unit_log.gt.0) call system_clock(t_start)

    if(nrhs.eq.1) then
       call dtrsv('Upper', trans, unit, nelim, dest, n, rhs(col), 1)
    else
       call dtrsm('Left', 'Upper', trans, unit, nelim, nrhs, one, &
            dest, n, rhs(col), ldr)
    endif

  end subroutine slv_solve

  !*************************************************
  !
  ! TASK_FUPD
  ! B_j <- B_j - L_ij B_i
  !
  subroutine slv_fwd_update(m, nelim, col, offset, index, dest, ldd, nrhs, &
       upd, ldu, rhs, ldr, xlocal)
    use spllt_data_mod
    implicit none

    integer, intent(in) :: m ! number of rows in block
    integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
    integer, intent(in) :: col ! start of block column variables in rhs
    integer, intent(in) :: offset ! offset into index we start at
    integer, dimension(*), intent(in) :: index
    integer, intent(in) :: ldd ! leading dimension of block
    real(wp), dimension(m*ldd), intent(in) :: dest ! holds destination block
    integer, intent(in) :: nrhs
    integer, intent(in) :: ldu  ! leading extent of upd
    real(wp), intent(inout) :: upd(ldu*nrhs) ! vector to update
    integer, intent(in) :: ldr  ! leading extent of rhs
    real(wp), intent(in) :: rhs(ldr*nrhs) ! rhs vector
    real(wp), dimension(*), intent(out) :: xlocal

    integer :: i
    integer :: j
    integer :: k
    integer :: r ! right hand side loop variable
    real(wp) :: w ! temporary work value
    !%%%  integer :: t_start, t_end, this_thread

    if(nelim.eq.0) return

    ! forward substitution
    if(nrhs.eq.1) then
       if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single rhs, BLAS 2

          call dgemv('T', nelim, m, -one, dest, ldd, rhs(col), 1, zero, &
               xlocal, 1)

          ! Copy xlocal out
          j = 1
          do i = offset, offset+m-1
             upd(index(i)) = upd(index(i)) + xlocal(j)
             j = j + 1
          end do
       else
!!! Single rhs, direct update
          j = 1
          do i = offset, offset+m-1
             w = zero
             do k = col, col+nelim-1
                w = w - dest(j)*rhs(k)
                j = j + 1
             end do
             j = j + (ldd-nelim)
             upd(index(i)) = upd(index(i)) + w
          end do
       endif
    else
!!! Multiple rhs, BLAS 3
       call dgemm('T', 'N', m, nrhs, nelim, -one, dest, ldd, rhs(col), ldr, &
            zero, xlocal, m)

       ! Copy xlocal out
       j = 1
       do i = offset, offset+m-1
          do r = 0, nrhs-1
             upd(index(i)+r*ldu) = upd(index(i)+r*ldu) + xlocal(j+r*m)
          end do
          j = j + 1
       end do
    endif

  end subroutine slv_fwd_update

  !*************************************************
  !
  ! TASK_BUPD
  ! B_i <- B_i - L_ij^-T B_j
  !
  subroutine slv_bwd_update(m, nelim, col, offset, index, dest, ldd, nrhs, rhs, &
       upd, ldr, xlocal)
    use spllt_data_mod
    implicit none

    integer, intent(in) :: m ! number of rows in block
    integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
    integer, intent(in) :: col ! start of block column variables in rhs
    integer, intent(in) :: offset ! offset into index we start at
    integer, dimension(*), intent(in) :: index
    integer, intent(in) :: ldd ! leading dimension of block
    real(wp), dimension(m*ldd), intent(in) :: dest ! holds block
    integer, intent(in) :: nrhs
    integer, intent(in) :: ldr  ! leading extent of rhs
    real(wp), intent(inout) :: rhs(ldr*nrhs)
    real(wp), intent(inout) :: upd(ldr*nrhs)
    real(wp), dimension(*), intent(out) :: xlocal

    integer :: i
    integer :: j
    integer :: k
    integer :: r ! right hand side loop variable
    real(wp) :: w ! temporary work variable
    !%%%  integer :: t_start, t_end, this_thread

    if(nelim.eq.0) return

    !%%%  if(control%unit_log.gt.0) call system_clock(t_start)

    ! backward substitution
    if(nrhs.eq.1) then
       if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single right-hand side, BLAS 2

          ! Copy xlocal in
          j = 1
          do i = offset, offset+m-1
             xlocal(j) = rhs(index(i))
             j = j + 1
          end do

          call dgemv('N', nelim, m, -one, dest, ldd, xlocal, 1, one, &
               upd(col), 1)
       else
!!! Single right-hand side, direct update
          j = 1
          do i = offset, offset+m-1
             w = rhs(index(i))
             do k = col, col + nelim - 1
                upd(k) = upd(k) - dest(j)*w
                j = j + 1
             end do
             j = j + (ldd-nelim)
          end do
       endif
    else
!!! Multiple RHS, BLAS 3

       ! Copy xlocal in
       j = 1
       do i = offset, offset+m-1
          do r = 0, nrhs-1
             xlocal(j+r*m) = rhs(index(i)+r*ldr)
          end do
          j = j + 1
       end do

       call dgemm('N', 'N', nelim, nrhs, m, -one, dest, ldd, xlocal, m, &
            one, upd(col), ldr)
    endif

  end subroutine slv_bwd_update

end module spllt_solve_mod
