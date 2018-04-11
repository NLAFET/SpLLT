 call task_manager%refresh_worker()
 threadID  = task_manager%workerID
!threadID = omp_get_thread_num()
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd block task", 13, threadID, timer)
#endif

call task_manager%incr_nrun()

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
  call trace_event_start(trace_id, threadID)
#endif

#if defined(SPLLT_VERBOSE)
  print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
  print *, p_dep
#endif


#if 1
  call solve_fwd_block_work(m, n, col, offset, p_index, p_lcol, sa, nrhs, &
      p_upd, p_rhs, threadID, nthread, ldr, p_xlocal, flops)

  call task_manager%nflop_performed(flops)
#else
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
  call task_manager%nflop_performed(n * n * nrhs + zero)
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
    call task_manager%nflop_performed(2 * (n * nrhs * m) + zero)
#endif
  endif
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, timer)
#endif
