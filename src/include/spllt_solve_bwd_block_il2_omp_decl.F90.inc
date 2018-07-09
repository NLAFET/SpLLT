!$omp default(none)                                           &
!$omp firstprivate(blkm, blkn, nrhs, n, ldy)                  &
!$omp firstprivate(trace_id, p_timer)                         &
!$omp firstprivate(p_y, p_lcol, p_index, p_rhs, p_order)      &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp firstprivate(task_manager)                              &
!$omp private(threadID, i, flops)                             &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))
