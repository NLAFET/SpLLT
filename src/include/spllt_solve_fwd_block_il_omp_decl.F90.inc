!$omp default(none)                                                       &
!$omp firstprivate(blkm, blkn, nrhs, col, sa, offset)                     &
!$omp firstprivate(n, ldr, bdr, ldu, bdu, tdu)                            &
!$omp firstprivate(nthread, trace_id, p_timer)                            &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal, p_rhsPtr)     &
!$omp firstprivate(p_bc, p_dep)                                           &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(i, j, r, threadID, flops)                                   &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))
