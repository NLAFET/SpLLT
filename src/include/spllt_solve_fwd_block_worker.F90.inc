!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
call task_manager%refresh_worker()
threadID  = task_manager%workerID

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

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd get RHS", 14, threadID, p_timer)
#endif
  !Gather part of RHS into y
  if(fkeep%nodes(node)%sblk_sa .eq. dblk) then
    do i = 1, blkn
      p_y(i, :) = p_rhs(p_order(p_index(i)), :)
    end do
  else
    do i = 1, blkn
      p_y(i, :) = p_rhs(p_order(p_index(i)), :) + p_y(i, :)
    end do
  end if
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(14, threadID, p_timer)
#endif
  !Compute the reduction if blk \in L_{:,1}
  if(fkeep%nodes(node)%sblk_sa .eq. dblk) then
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd update UPD", 15, threadID, p_timer)
#endif
    do i = 1, nwdep
      ! reduce children dep into p_y
      call fwd_update_upd(fkeep, dblk, p_wdep(i))
    end do
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(15, threadID, p_timer)
#endif
  end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd solve KERNEL", 16, threadID, p_timer)
#endif
  call slv_solve(blkm, blkn, p_lcol, 'Transpose    ', 'Non-unit', &
    nrhs, p_y, ldy)
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(16, threadID, p_timer)
#endif

#if defined(SPLLT_PROFILING_FLOP)
  flops = blkn * blkn * nrhs
  call task_manager%nflop_performed(flops)
#endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, p_timer) !, flops)
#endif
