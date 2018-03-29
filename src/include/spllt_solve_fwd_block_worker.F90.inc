threadID  = omp_get_thread_num()
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd block task", 13, threadID, timer)
#endif
scheduler%task_info(threadID)%ntask_run = scheduler%task_info(threadID)%&
  &ntask_run + 1

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
  call trace_event_start(trace_id, threadID)
#endif

#if defined(SPLLT_VERBOSE)
  print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
  print *, p_dep
#endif

! call spllt_tic("Reset workspace", 13, threadID, timer)
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
  do r = 0, nrhs-1
    do j = 1, nthread
      do i = col + r*ldr, col+n-1 + r*ldr
        p_rhs(i)    = p_rhs(i) + p_upd(i, j)
        p_upd(i,j)  = zero ! Reset in case of bwd solve
      end do
    end do
  end do


  ! Perform triangular solve
  call slv_solve(n, n, col, p_lcol(sa:sa+n*n-1),    &
    'Transpose    ', 'Non-unit', nrhs, p_rhs, ldr)
#if defined(SPLLT_PROFILING_FLOP)
  scheduler%task_info(threadID)%nflop = scheduler%task_info(threadID)%nflop + &
    (n * n * nrhs)
#endif
  offset = offset + n

  ! Deal with any left over trapezoidal part of diagonal block
  m = m - n
  if(m .gt. 0) then
    sa = sa + n * n
    call slv_fwd_update(m, n, col, offset, p_index,             &
      p_lcol(sa : sa + n * m - 1), n, nrhs,                     &
      p_upd(:, threadID + 1), ldr, p_rhs,                       &
      ldr, p_xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    scheduler%task_info(threadID)%nflop = scheduler%task_info(threadID)%nflop &
      + 2 * (n * nrhs * m)
#endif
  endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, timer)
#endif
