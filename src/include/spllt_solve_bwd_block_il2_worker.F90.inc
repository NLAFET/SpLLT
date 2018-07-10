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

  !Gather data into xlocal if belongs to last block column of L
  if(fkeep%nodes(node)%sblk_en .eq. dblk) then
    do i = 1, ndep
      ! reduce ancestor dep into p_y
      call bwd_update_upd(fkeep, dblk, p_dep(i))
    end do
  end if
 
  call slv_solve_ileave2(blkn, blkm, p_lcol, 'Non-Transpose', 'Non-unit', &
    nrhs, p_y, ldy)

  !Scatter the result into the solution vector
  do i = 1, blkn
   !p_rhs(p_order(p_index(i)), :) = p_rhs(p_order(p_index(i)), :) + p_y(i, :)
    p_rhs(p_order(p_index(i)), :) = p_y(i, :)
   !print *, "Copy", p_y(i,:)
   !print *, "Into p_rhs(", p_order(p_index(i)), ") =", p_rhs(p_order(p_index(i)), :)
  end do

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
