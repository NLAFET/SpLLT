!$omp default(none)                                           &
!$omp firstprivate(blkm, blkn, nrhs, col, sa, offset)         &
!$omp firstprivate(n, tdu)                                    &
!$omp firstprivate(nthread, trace_id, p_timer)                &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_rhsPtr, p_indir_rhs)                     &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp firstprivate(task_manager)                              &
!$omp private(threadID, r, j, i, flops)                       &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))
