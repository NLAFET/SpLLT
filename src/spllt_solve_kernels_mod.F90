module spllt_solve_kernels_mod
  
contains


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

    integer,            intent(in)    :: n      ! leading dimension of 
                                                ! diag. block
    integer,            intent(in)    :: nelim  ! number eliminations 
                                                ! (immediate return if =0)
    integer,            intent(in)    :: col    ! start of block column 
                                                ! variables in rhs
    real(wp),           intent(in)    :: dest(*)! holds destination block
    character(len=13),  intent(in)    :: trans  ! set to :
                                                ! 'Transpose    ' 
                                                ! for forward substitution 
                                                ! 'Non-Transpose' 
                                                ! for back substitution
    character(len=8),   intent(in)    :: unit   ! set to 
                                                ! 'Non-Unit' 
                                                ! for positive-definite case
    integer,            intent(in)    :: nrhs   ! number of right-hand sides
    integer,            intent(in)    :: ldr    ! leading extent of rhs
    real(wp),           intent(inout) :: rhs(ldr*nrhs)

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

    integer,  intent(in)    :: m              ! number of rows in block
    integer,  intent(in)    :: nelim          ! # eliminations 
                                              ! (immediate return if =0)
    integer,  intent(in)    :: col            ! start of block column variables
                                              ! in rhs
    integer,  intent(in)    :: offset         ! offset into index we start at
    integer,  intent(in)    :: index(*)
    integer,  intent(in)    :: ldd            ! leading dimension of block
    real(wp), intent(in)    :: dest(m*ldd)    ! holds destination block
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: ldu            ! leading extent of upd
    real(wp), intent(inout) :: upd(ldu*nrhs)  ! vector to update
    integer,  intent(in)    :: ldr            ! leading extent of rhs
    real(wp), intent(in)    :: rhs(ldr*nrhs)  ! rhs vector
    real(wp), intent(out)   :: xlocal(*)

    integer   :: i
    integer   :: j
    integer   :: k
    integer   :: r ! right hand side loop variable
    real(wp)  :: w ! temporary work value
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
  subroutine slv_bwd_update(m, nelim, col, offset, index, dest, ldd, nrhs, &
      rhs, upd, ldr, xlocal)
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
    real(wp), dimension(:), intent(out) :: xlocal

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
!      do i = offset, offset+m-1
!         do r = 0, nrhs-1
!            xlocal(j+r*m) = rhs(index(i)+r*ldr)
!         end do
!         j = j + 1
!      end do
       !Test change order of loop
       do r = 0, nrhs-1
          j = 1
          do i = offset, offset+m-1
            xlocal(j+r*m) = rhs(index(i)+r*ldr)
            j = j + 1
          end do
       end do

       call dgemm('N', 'N', nelim, nrhs, m, -one, dest, ldd, xlocal, m, &
            one, upd(col), ldr)
    endif

  end subroutine slv_bwd_update


  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Kernel of solve steps
  !
  subroutine solve_bwd_block_work(m, n, col, offset, index, lcol, sa, nrhs, &
      upd, rhs, threadID, nthread, ldr, xlocal, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(in)    :: m
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: col
    integer,                  intent(in)    :: offset
    integer,                  intent(in)    :: sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: nthread
    integer,                  intent(in)    :: ldr
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:,:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)
    double precision,         intent(out)   :: flops

    integer :: i, r, j

    flops = zero

    ! Perform retangular update from diagonal block
    if(m .gt. n) then
      call slv_bwd_update(m - n, n, col, offset + n, index,     &
        lcol(sa + n * n : sa + n * m - 1), n, nrhs, rhs,      &
        upd(:, threadID + 1), ldr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (n * nrhs * (m-n))
#endif
    endif

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          rhs(i)    = rhs(i) + upd(i, j)
          upd(i,j)  = zero ! Reset in case of next fwd solve
        end do
      end do
    end do

    ! Perform triangular solve
    call slv_solve(n, n, col, lcol(sa : sa + n * n - 1), &
         'Non-Transpose', 'Non-unit', nrhs, rhs, ldr)

#if defined(SPLLT_PROFILING_FLOP)
    flops = flops + n * n * nrhs
#endif

  end subroutine solve_bwd_block_work



  subroutine solve_bwd_update_work(m, n, col, offset, index, lcol, blk_sa, &
      nrhs, upd, rhs, threadID, ldr, xlocal, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(in)    :: m
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: col
    integer,                  intent(in)    :: offset
    integer,                  intent(in)    :: blk_sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: ldr
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:,:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)
    double precision,         intent(out)   :: flops

    call  slv_bwd_update(m, n, col, offset, index,      &
          lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,   &
          rhs, upd(:, threadID + 1), ldr,               &
          xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (n * nrhs * m)
#endif

  end subroutine solve_bwd_update_work



  subroutine solve_fwd_update_work(m, n, col, offset, index, lcol, blk_sa, &
      nrhs, upd, rhs, threadID, ldr, xlocal, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(in)    :: m
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: col
    integer,                  intent(in)    :: offset
    integer,                  intent(in)    :: blk_sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: ldr
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:,:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)
    double precision,         intent(out)   :: flops

    call slv_fwd_update(m, n, col, offset, index,         &
      lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
      upd(:, threadID + 1), ldr, rhs,                   &
      ldr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (m * nrhs * n)
#endif

  end subroutine solve_fwd_update_work

  subroutine solve_fwd_block_work(m, n, col, offset, index, lcol, sa, nrhs, &
      upd, rhs, threadID, nthread, ldr, xlocal, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(inout) :: m
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: col
    integer,                  intent(inout) :: offset
    integer,                  intent(inout) :: sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: nthread
    integer,                  intent(in)    :: ldr
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:,:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)
    double precision,         intent(out)   :: flops

    integer :: i, r, j

    flops = zero
! ! Reset p_upd if this is the first diagonal block of the node
! if(dblk .eq. blk_sa) then
!!  print *, "RReset ", ubound(p_extra_row, 1) - lbound(p_extra_row, 1) + 1, &
!!    " elts of the workspace of blk ", dblk
!     do i = lbound(p_extra_row, 1), ubound(p_extra_row, 1)
!       if(.not. p_workspace_reset(p_extra_row(i))) then
!         !!$omp atomic write
!         p_workspace_reset(p_extra_row(i)) = .true. ! not safe for me ...
!         !!$omp end atomic
!         do r = 0, nrhs-1
!           do j = 1, nthread
!             p_upd(p_extra_row(i) + r*ldr, j)  = zero
!           end do
!         end do
!       end if
!     end do
! end if
! call spllt_tac(13, threadID, timer)

 !if(dblk .eq. blk_sa) then
 !  print *, "Check rhs_local for blk ", dblk
 !  do r = 0, nrhs - 1
 !    do j = 1, nthread
 !      do i = lbound(p_extra_row, 1), ubound(p_extra_row, 1)
 !        if(p_upd(p_extra_row(i) + r*ldr, j) .ne. 0) then
 !          print *, "Error ? blk ", dblk, " => row ", p_extra_row(i)
 !        end if
 !      end do
 !    end do
 !  end do
 !end if

    ! Sum contributions to rhs
    do r = 0, nrhs - 1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          rhs(i)    = rhs(i) + upd(i, j)
          upd(i,j)  = zero ! Reset in case of bwd solve
        end do
      end do
    end do


    ! Perform triangular solve
    call slv_solve(n, n, col, lcol(sa:sa+n*n-1),    &
      'Transpose    ', 'Non-unit', nrhs, rhs, ldr)
#if defined(SPLLT_PROFILING_FLOP)
    flops = n * n * nrhs
#endif
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m .gt. 0) then
      sa = sa + n * n
      call slv_fwd_update(m, n, col, offset, index,             &
        lcol(sa : sa + n * m - 1), n, nrhs,                     &
        upd(:, threadID + 1), ldr, rhs,                         &
        ldr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (n * nrhs * m)
#endif
    endif
  end subroutine solve_fwd_block_work
end module spllt_solve_kernels_mod
