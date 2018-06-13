 call task_manager%refresh_worker()
 threadID  = task_manager%workerID
!threadID = omp_get_thread_num()
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd block task", 13, threadID, p_timer)
#endif

call task_manager%incr_nrun()

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
  call trace_event_start(trace_id, threadID)
#endif

#if defined(SPLLT_VERBOSE)
  print '(a, i3, a, i3)', "SLV-FWD  Task dep of ", dblk, " [in : "
  print *, p_dep
#endif

  call solve_fwd_block_work_ileave(dblk, p_rhsPtr, blkm, blkn, col, offset, p_index, p_lcol, sa, nrhs, &
      p_upd, ldu, bdu, tdu, p_rhs, n, ldr, bdr, p_xlocal, threadID, nthread, flops)

  call task_manager%nflop_performed(flops)

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, p_timer, flops)
#endif
