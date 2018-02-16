module spllt_solve_task_mod
contains
  
  !*************************************************  
  !
  ! Forward solve with block on diagoanl
  !
  subroutine spllt_solve_fwd_block_task(dblk, nrhs, upd, rhs, ldr, xlocal, &
      fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    implicit none
    
    integer, intent(in)                     :: dblk ! Index of diagonal block
    integer, intent(in)                     :: nrhs ! Number of RHS
    integer, intent(in)                     :: ldr  ! Leading dimension of RHS
    real(wp), intent(inout)                 :: upd(:,:)
    real(wp), intent(inout)                 :: rhs(ldr * nrhs)
    real(wp), dimension(:,:), intent(inout) :: xlocal
    type(spllt_fkeep), target, intent(in)   :: fkeep
    integer, intent(in)                     :: trace_id
    
    ! Node info
    integer :: sa
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa
    integer :: bcol, dcol, col
    integer :: offset
    integer :: node
    integer :: i, j, r
    integer :: threadID, nthread
    integer, dimension(:), pointer  :: p_index
    real(wp), dimension(:), pointer :: p_lcol
    real(wp), dimension(:), pointer :: p_lcol_update
    integer                         :: blk_dep_update
    type(spllt_block), pointer      :: p_blk_dep_update
    type(spllt_block), pointer      :: p_dblk   
        
    nthread   = omp_get_num_threads()
    threadID  = omp_get_thread_num()

    ! Get block info
    node      = fkeep%bc(dblk)%node
    m         = fkeep%bc(dblk)%blkm
    n         = fkeep%bc(dblk)%blkn
    sa        = fkeep%bc(dblk)%sa
    bcol      = fkeep%bc(dblk)%bcol ! Current block column
    dcol      = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col       = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
    p_dblk    => fkeep%bc(dblk)   
    
    blk_dep_update = fwd_update_dependency(fkeep, dblk)

    p_lcol_update => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol

    p_blk_dep_update    => fkeep%bc(blk_dep_update)   

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          rhs(i)    = rhs(i) + upd(i, j)
          upd(i,j)  = zero ! Reset in case of bwd solve
        end do
      end do
    end do

    !$omp task firstprivate(m, n, col, offset)                      &
    !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)              &
    !$omp firstprivate(blk_dep_update, dblk)                        &
    !$omp firstprivate(p_blk_dep_update   , p_dblk   )              &
    !$omp shared(fkeep, rhs, xlocal, upd)                           &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))            &
    !$omp depend(inout: p_lcol(sa))
!   print *, "[spllt_solve_fwd_block_task] solved by ", omp_get_thread_num()
!   call trace_event_start(trace_id, omp_get_thread_num())

    print "(a, i3, a, i4, a, i4)", "th :", threadID, " dblk :", dblk, &
      " has a dependency with ", blk_dep_update

    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa:sa+n*n-1),    &
         'Transpose    ', 'Non-unit', nrhs, rhs, ldr)
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m .gt. 0) then
       sa = sa + n * n
       call slv_fwd_update(m, n, col, offset, p_index,                &
            p_lcol(sa : sa + n * m - 1), n, nrhs,                     &
            upd(:, omp_get_thread_num() + 1), ldr, rhs,               &
            ldr, xlocal(:, omp_get_thread_num() + 1))
    endif

!   call trace_event_stop (trace_id, omp_get_thread_num())
    !$omp end task
  end subroutine spllt_solve_fwd_block_task


  subroutine spllt_solve_fwd_update_task(blk, node, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    implicit none
    
    integer, intent(in)                     :: blk  ! Index of block
    integer, intent(in)                     :: node
    integer, intent(in)                     :: nrhs ! Number of RHS
    integer, intent(in)                     :: ldr  ! Leading dimension of RHS
    real(wp), intent(inout)                 :: upd(:,:)        
    real(wp), intent(in)                    :: rhs(ldr*nrhs)
    real(wp), intent(out)                   :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)   :: fkeep
    integer, intent(in)                     :: trace_id

    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID, nthread
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    real(wp), pointer           :: p_lcol_solve(:)
    integer                     :: blk_dep_update
    integer                     :: blk_dep_solve
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve
    type(spllt_block), pointer  :: p_blk 

    ! Establish variables describing block
    n         = fkeep%bc(blk)%blkn
    m         = fkeep%bc(blk)%blkm
    blk_sa    = fkeep%bc(blk)%sa
    bcol      = fkeep%bc(blk)%bcol
    dcol      = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col       = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
    offset    = offset + (blk-fkeep%bc(blk)%dblk) &
      * fkeep%nodes(node)%nb ! this blk
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
    p_blk     => fkeep%bc(blk)

    blk_dep_update = fwd_update_dependency(fkeep, blk)
    blk_dep_solve  = fwd_solve_dependency(fkeep, blk)

    p_lcol_update => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
    p_lcol_solve  => fkeep%lfact(fkeep%bc(blk_dep_solve )%bcol)%lcol

    p_blk_dep_update  => fkeep%bc(blk_dep_update)
    p_blk_dep_solve   => fkeep%bc(blk_dep_solve)

    !$omp task firstprivate(m, n, col, offset)                    &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp firstprivate(p_blk_dep_update, p_blk_dep_solve)         &
    !$omp firstprivate(p_blk)                                     &
    !$omp firstprivate(blk, blk_dep_update)                       &
    !$omp shared(fkeep, xlocal, rhs, upd)                         &
    !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(inout: p_lcol(blk_sa))
!   print *, "[spllt_solve_fwd_update_task] solved by ", omp_get_thread_num()
!   call trace_event_start(trace_id, omp_get_thread_num())
    
    print '(a, i3, a, i4, a, i4, a, i4)', "th :", threadID, "  blk :", blk, &
      " has a dependency with ", blk_dep_update, " and ", fkeep%bc(blk)%dblk

    call slv_fwd_update(m, n, col, offset, p_index,         &
      p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
      upd(:, threadID + 1), ldr, rhs,           &
      ldr, xlocal(:, threadID + 1))

!   call trace_event_stop (trace_id, omp_get_thread_num())
    !$omp end task
  end subroutine spllt_solve_fwd_update_task

  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         
  subroutine spllt_solve_bwd_block_task(dblk, nrhs, rhs, ldr, xlocal, &
      fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num
    use trace_mod
    implicit none

    integer, intent(in)                     :: dblk ! Index of diagonal block
    integer, intent(in)                     :: nrhs ! Number of RHS
    integer, intent(in)                     :: ldr  ! Leading dimension of RHS
    real(wp), intent(inout)                 :: rhs(ldr, nrhs)
    real(wp), dimension(:,:), intent(inout) :: xlocal
    type(spllt_fkeep), target, intent(in)   :: fkeep
    integer, intent(in)                     :: trace_id
    
    ! Node info
    integer                         :: sa
    ! Block info
    integer                         :: m, n ! Block dimension
    integer                         :: blk_sa
    integer                         :: bcol, dcol, col
    integer                         :: offset
    integer                         :: node
    integer, dimension(:), pointer  :: p_index
    real(wp), dimension(:), pointer :: p_lcol

    node = fkeep%bc(dblk)%node

    ! print *, "[spllt_solve_bwd_block_task] node = ", node

    ! Get block info
    n       = fkeep%bc(dblk)%blkn
    m       = fkeep%bc(dblk)%blkm
    sa      = fkeep%bc(dblk)%sa
    bcol    = fkeep%bc(dblk)%bcol ! Current block column
    col     = calc_col(fkeep%nodes(node), fkeep%bc(dblk)) ! current bcol
    col     = fkeep%nodes(node)%sa + (col-1)*fkeep%nodes(node)%nb
    offset  = col - fkeep%nodes(node)%sa + 1
    p_index => fkeep%nodes(node)%index
    p_lcol  => fkeep%lfact(bcol)%lcol

    !$omp task                                                    &
    !$omp firstprivate(m, n, col, offset, node, bcol)             &
    !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)            &
    !$omp shared(fkeep, xlocal, rhs)
!   print *, "[spllt_solve_bwd_block_task] treated by ", omp_get_thread_num()
!   print *, "work on dblk ", dblk, "which belongs to block column ", & 
!     fkeep%bc(dblk)%bcol
!   call trace_event_start(trace_id, omp_get_thread_num())

    ! Perform retangular update from diagonal block
    if(m .gt. n) then
       call slv_bwd_update(m - n, n, col, offset + n, p_index, &
            p_lcol(sa + n * n : sa + n * m - 1), n, nrhs, rhs, &
            rhs, ldr, xlocal(:, omp_get_thread_num() + 1))
    endif


    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa : sa + n * n - 1), &
         'Non-Transpose', 'Non-unit', nrhs, rhs, ldr)
    
!   call trace_event_stop (trace_id, omp_get_thread_num())
    !$omp end task

  end subroutine spllt_solve_bwd_block_task

  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         
  subroutine spllt_solve_bwd_udpate_task(blk, node, nrhs, rhs, ldr, xlocal, &
      fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num
    use trace_mod
    implicit none

    integer, intent(in)                     :: blk  ! Index of block 
    integer, intent(in)                     :: node 
    integer, intent(in)                     :: nrhs ! Number of RHS
    integer, intent(in)                     :: ldr  ! Leading dimension of RHS
    real(wp), intent(inout)                 :: rhs(ldr, nrhs)
    real(wp), dimension(:,:), intent(inout) :: xlocal
    type(spllt_fkeep), target, intent(in)   :: fkeep
    integer, intent(in)                     :: trace_id
    
    ! Block info
    integer                         :: m, n         ! Block dimension
    integer                         :: blk_sa
    integer                         :: bcol, dcol, col
    integer                         :: offset
    integer, dimension(:), pointer  :: p_index
    real(wp), dimension(:), pointer :: p_lcol

    ! print *, "[spllt_solve_bwd_udpate_task] node = ", node

    !
    ! Backward update with block on diagoanl
    !

    ! Establish variables describing block
    n       = fkeep%bc(blk)%blkn
    m       = fkeep%bc(blk)%blkm
    blk_sa  = fkeep%bc(blk)%sa
    bcol    = fkeep%bc(blk)%bcol
    dcol    = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col     = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset  = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
    offset  = offset + (blk-fkeep%bc(blk)%dblk) * fkeep%nodes(node)%nb !this blk
    p_index => fkeep%nodes(node)%index
    p_lcol  => fkeep%lfact(bcol)%lcol

    !$omp task                                                    &
    !$omp firstprivate(m, n, col, offset, node, bcol)             &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp shared(fkeep, xlocal, rhs)
!   print *, "[spllt_solve_bwd_udpate_task] Thread id ", omp_get_thread_num()
!   print *, "work on blk  ", blk, " which belongs to ", fkeep%bc(blk)%bcol, &
!     ", and depends on ", fkeep%bc(blk)%dblk
!   print *, "last_blk ", fkeep%bc(blk)%last_blk
!   call trace_event_start(trace_id, omp_get_thread_num())
    call  slv_bwd_update(m, n, col, offset, p_index, &
          p_lcol(blk_sa:blk_sa+n*m-1), n, nrhs, rhs, &
          rhs, ldr, xlocal(:, omp_get_thread_num() + 1))
!   call trace_event_stop (trace_id, omp_get_thread_num())
    !$omp end task

  end subroutine spllt_solve_bwd_udpate_task

  !*************************************************
  !
  ! This function calculates column of a node we are on  
  integer function calc_col(node, bc)
    use spllt_data_mod    
    implicit none

    type(spllt_node), intent(in) :: node
    type(spllt_block), intent(in) :: bc

    calc_col = (size(node%index)-1)/node%nb + 1 ! no. row blks for node

    calc_col = calc_col - (bc%last_blk - bc%dblk + 1) + 1 ! column of node

  end function calc_col

  subroutine calc_col1(node, bc, col)
    use spllt_data_mod    
    implicit none

    type(spllt_node), intent(in) :: node
    type(spllt_block), intent(in) :: bc
    integer, intent(out) :: col

    col = (size(node%index)-1)/node%nb + 1 ! no. row blks for node

    col = col - (bc%last_blk - bc%dblk + 1) + 1 ! column of node

  end subroutine calc_col1

  integer function fwd_update_dependency(fkeep, blk)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)   :: fkeep
    integer, intent(in)             :: blk  ! Index of block 
    
    integer :: previous_dblk
    integer :: last_previous_dblk
    integer :: diff_bcol
    integer :: diff_previous_bcol

    fwd_update_dependency = blk

    if(fkeep%bc(blk)%bcol .gt. 1) then
      previous_dblk       = fkeep%bc(fkeep%bc(blk)%dblk - 1)%dblk
      last_previous_dblk  = fkeep%bc(previous_dblk)%last_blk
      diff_bcol           = blk - fkeep%bc(blk)%dblk
      diff_previous_bcol  = last_previous_dblk - previous_dblk

      if(last_previous_dblk .gt. (previous_dblk + diff_bcol)) then
        fwd_update_dependency = previous_dblk + diff_bcol + 1
      end if
    end if

  end function fwd_update_dependency

  integer function fwd_solve_dependency(fkeep, blk)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)   :: fkeep
    integer, intent(in)             :: blk  ! Index of block 

    if(fkeep%bc(blk)%dblk .ne. blk) then
      fwd_solve_dependency = fkeep%bc(blk)%dblk
    else
      fwd_solve_dependency = blk
    end if

  end function fwd_solve_dependency

end module spllt_solve_task_mod
