!$omp default(none)                                                       &
!$omp firstprivate(fkeep)                                                 &
!$omp firstprivate(blkm, blkn, nrhs, n, ldy, ldx, node, ndep)             &
!$omp firstprivate(p_lcol, p_y, p_xlocal)                                 &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                        &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, i, flops)                                         &
!$omp firstprivate(blk) depend(inout: p_bc(blk))
