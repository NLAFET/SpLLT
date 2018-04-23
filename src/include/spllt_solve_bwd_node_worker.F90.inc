call sub_task_manager%refresh_worker()

do i = en, sa, -1
  call solve_bwd_node(nrhs, p_rhs, ldr, fkeep, i, p_xlocal, &
    p_rhs_local, sub_task_manager)
end do
