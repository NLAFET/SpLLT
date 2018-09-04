call task_manager%refresh_worker()
threadID  = task_manager%workerID

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

  !Gather data into xlocal if belongs to last block column of L
  if(fkeep%nodes(node)%sblk_en .eq. fkeep%sbc(blk)%last_blk) then
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd update UPD", 14, threadID, p_timer)
#endif
   !print *, "Reduce in blk", blk
    do i = 1, nwdep
      ! reduce ancestor dep into p_y
   !  print *, "with ", p_dep(i)
      call bwd_update_upd(fkeep, blk, p_wdep(i))
    end do
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(14, threadID, p_timer)
#endif
  end if
 
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd update KERNEL", 15, threadID, p_timer)
#endif
  call slv_bwd_update(blkm, blkn, p_lcol, blkn, n, nrhs, &
    p_xlocal, ldx, p_y, ldy)
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(15, threadID, p_timer)
#endif

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
