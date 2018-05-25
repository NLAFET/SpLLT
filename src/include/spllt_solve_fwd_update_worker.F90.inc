 call task_manager%refresh_worker()
 threadID  = task_manager%workerID
!threadID = omp_get_thread_num()
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd update task", 13, threadID, p_timer)
#endif
call task_manager%incr_nrun()

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
  call trace_event_start(trace_id, threadID)
#endif

#if defined(SPLLT_VERBOSE)
  print '(a, i3, a)', "UPD-FWD  Task dep of ", blk, " [in : ]"
  print *, p_dep
#endif 
#if 1
  call solve_fwd_update_work(m, n, col, offset, p_index, p_lcol, blk_sa, &
    nrhs, p_upd, p_rhs, threadID, ldr, p_xlocal, flops)

#if defined(SPLLT_PROFILING_FLOP)
  call task_manager%nflop_performed(flops)
#endif
#else
  call slv_fwd_update(m, n, col, offset, p_index,         &
    p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
    p_upd(:, threadID + 1), ldr, p_rhs,                   &
    ldr, p_xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
  call task_manager%nflop_performed(2 * (m * nrhs * n) + zero)
#endif
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, p_timer)
#endif
