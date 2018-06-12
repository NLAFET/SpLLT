!$omp default(none)                                                       &
!$omp firstprivate(blkm, blkn, nrhs, col, blk_sa, offset)                 &
!$omp firstprivate(n, ldu, bdu, tdu, bdr)                                 &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                        &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))
