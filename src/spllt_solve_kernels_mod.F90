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


#if 0
  subroutine solve_bwd_block_work(m, n, col, offset, index, lcol, sa, nrhs, &
      upd, threadID, ldr, xlocal, trace_id, task_manager)
    implicit none

#if defined(SPLLT_OMP_TRACE)
    call trace_event_start(trace_id, threadID)
#endif

#if defined(SPLLT_VERBOSE)
    print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
    print *, p_dep
#endif

    ! Perform retangular update from diagonal block
    if(m .gt. n) then
      call slv_bwd_update(m - n, n, col, offset + n, p_index,     &
        lcol(sa + n * n : sa + n * m - 1), n, nrhs, rhs,      &
        upd(:, threadID + 1), ldr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
      call task_manager%nflop_performed(2 * (n * nrhs * (m-n)) + zero)
#endif
    endif

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          p_rhs(i)    = p_rhs(i) + p_upd(i, j)
          p_upd(i,j)  = zero ! Reset in case of next fwd solve
        end do
      end do
    end do

    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa : sa + n * n - 1), &
         'Non-Transpose', 'Non-unit', nrhs, p_rhs, ldr)

#if defined(SPLLT_PROFILING_FLOP)
    call task_manager%nflop_performed(n * n * nrhs + zero)
#endif

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(trace_id, threadID)
#endif
  end subroutine solve_fwd_block_work
#endif

end module spllt_solve_kernels_mod
