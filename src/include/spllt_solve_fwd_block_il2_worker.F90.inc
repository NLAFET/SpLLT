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

  !Gather part of RHS into y
  do i = 1, blkn
   !p_y(i, :) = rhs(p_order(p_index(i)), :)
   !p_y(i, :) = p_y(i, :) + p_rhs(p_index(p_order(i)), :)
    p_y(i, :) = p_y(i, :) + p_rhs(p_order(p_index(i)), :)
  end do

  call slv_solve_ileave2(blkm, blkn, p_lcol, 'Transpose    ', 'Non-unit', &
    nrhs, p_y, ldy)

#if defined(SPLLT_PROFILING_FLOP)
  flops = blkn * blkn * nrhs
  call task_manager%nflop_performed(flops)
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, p_timer, flops)
#endif
