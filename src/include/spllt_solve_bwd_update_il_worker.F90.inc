 call task_manager%refresh_worker()
 threadID  = task_manager%workerID
!threadID = omp_get_thread_num()
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd update task", 13, threadID, p_timer)
#endif
call task_manager%incr_nrun()

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
  call trace_event_start(trace_id, threadID)
#endif

#if defined(SPLLT_VERBOSE)
  print '(a, i3, a)', "UPD-BWD  Task dep of ", blk, " [in : "
  print *, p_dep
#endif

  call solve_bwd_update_work_ileave(blk, p_rhsPtr, p_indir_rhs, blkm,       &
    blkn, col, offset, p_index, p_lcol, blk_sa, nrhs, p_upd, tdu, p_rhs, n, &
    p_xlocal, threadID, flops)

#if defined(SPLLT_PROFILING_FLOP)
  call task_manager%nflop_performed(flops)
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, p_timer)
#endif
