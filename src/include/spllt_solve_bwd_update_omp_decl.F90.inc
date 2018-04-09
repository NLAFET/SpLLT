!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))
