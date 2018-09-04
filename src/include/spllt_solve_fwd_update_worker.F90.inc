call task_manager%refresh_worker()
threadID  = task_manager%workerID

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

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd update KERNEL", 15, threadID, p_timer)
#endif
  call slv_fwd_update(blkm, blkn, p_lcol, blkn, n, nrhs, p_y, ldy, &
    p_xlocal, ldx, fkeep%nodes(node)%sblk_sa .eq. fkeep%sbc(blk)%dblk)
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(15, threadID, p_timer)
#endif

  !Compute the reduction if blk \in L_{:,1}
  if(fkeep%nodes(node)%sblk_sa .eq. fkeep%sbc(blk)%dblk) then
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd update UPD", 14, threadID, p_timer)
#endif
    do i = 1, nwdep
      ! reduce children dep into p_y
      call fwd_update_upd(fkeep, blk, p_wdep(i))
    end do
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(14, threadID, p_timer)
#endif
  end if

#if defined(SPLLT_PROFILING_FLOP)
  flops = 2 * (blkm * nrhs * blkn)
  call task_manager%nflop_performed(flops)
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, p_timer)
#endif
