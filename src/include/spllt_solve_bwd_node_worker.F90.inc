!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
call sub_task_manager%refresh_worker()
call sub_task_manager%incr_nrun()

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
    trace_id = sub_task_manager%trace_ids(trace_bwd_subtree_pos)
    call trace_event_start(trace_id, sub_task_manager%workerID)
#endif
#if defined(SPLLT_VERBOSE)
  print '(a, i3, a)', "SLV-BWD  Task dep of ", blk_en, " [in : "
  print *, p_dep
#endif
do i = en, sa, -1
  call solve_bwd_node(nrhs, p_rhs, n, fkeep, i, &
      sub_task_manager, no_trace)
end do

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(trace_id, sub_task_manager%workerID)
#endif
end if
