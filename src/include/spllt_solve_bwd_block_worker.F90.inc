threadID = omp_get_thread_num()
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd block task", 13, threadID, timer)
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

  ! Perform retangular update from diagonal block
  if(m .gt. n) then
    call slv_bwd_update(m - n, n, col, offset + n, p_index,     &
      p_lcol(sa + n * n : sa + n * m - 1), n, nrhs, p_rhs,      &
      p_upd(:, threadID + 1), ldr, p_xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    scheduler%task_info(threadID)%nflop = scheduler%task_info(threadID)%nflop &
      + 2 * (n * nrhs * (m-n))
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
  scheduler%task_info(threadID)%nflop = scheduler%task_info(threadID)%nflop + &
    (n * n * nrhs)
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop(trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, timer)
#endif
