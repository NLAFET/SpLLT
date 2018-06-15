call task_manager%refresh_worker()
threadID = task_manager%workerID
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd block task", 13, threadID, p_timer)
#endif
call task_manager%incr_nrun()

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
   call trace_event_start(trace_id, threadID)
#endif
 
#if defined(SPLLT_VERBOSE)
   print '(a, i3, a, i3)', "SLV-BWD  Task dep of ", dblk, " [in : "
   print *, p_dep
#endif
 
  call solve_bwd_block_work_ileave(dblk, p_rhsPtr, p_indir_rhs, blkm, blkn, col, &
    offset, p_index, p_lcol, sa, nrhs, p_upd, tdu, p_rhs, n, p_xlocal,    &
    threadID, nthread, flops)

#if defined(SPLLT_PROFILING_FLOP)
     call task_manager%nflop_performed(flops)
#endif
 
#if defined(SPLLT_OMP_TRACE)
   call trace_event_stop(trace_id, threadID)
#endif

end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, p_timer, flops)
#endif
