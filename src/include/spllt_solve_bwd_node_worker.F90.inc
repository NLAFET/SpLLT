call sub_task_manager%refresh_worker()

#if defined(SPLLT_OMP_TRACE)
    trace_id = task_manager%trace_ids(trace_bwd_subtree_pos)
    call trace_event_start(trace_id, sub_task_manager%workerID)
#endif

do i = en, sa, -1
  call solve_bwd_node(nrhs, p_rhs, ldr, fkeep, i, p_xlocal, &
    p_rhs_local, sub_task_manager, no_trace)
end do

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(trace_id, sub_task_manager%workerID)
#endif
