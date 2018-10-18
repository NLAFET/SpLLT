!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
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

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd solve KERNEL", 15, threadID, p_timer)
#endif
  call slv_solve(blkn, blkm, p_lcol, 'Non-Transpose', 'Non-unit', &
    nrhs, p_y, ldy)
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(15, threadID, p_timer, flops)
#endif

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd set X", 14, threadID, p_timer)
#endif
  !Scatter the result into the solution vector
  do i = 1, blkn
    p_rhs(p_order(p_index(i)), :) = p_y(i, :)
  end do
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(14, threadID, p_timer, flops)
#endif

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
