threadID  = omp_get_thread_num()
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd update task", 13, threadID, timer)
#endif
scheduler%task_info(threadID)%ntask_run = scheduler%task_info(threadID)%&
  &ntask_run + 1

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
  call trace_event_start(trace_id, threadID)
#endif

#if defined(SPLLT_VERBOSE)
  print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : ]"
  print *, p_dep
#endif 
  call slv_fwd_update(m, n, col, offset, p_index,         &
    p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
    p_upd(:, threadID + 1), ldr, p_rhs,                   &
    ldr, p_xlocal(:, threadID + 1))
#if defined(SPLLT_PROFILING_FLOP)
  scheduler%task_info(threadID)%nflop = scheduler%task_info(threadID)%nflop+&
    2 * (m * nrhs * n)
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, timer)
#endif
