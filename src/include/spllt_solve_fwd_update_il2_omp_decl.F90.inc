!$omp default(none)                                                       &
!$omp firstprivate(fkeep)                                                 &
!$omp firstprivate(blkm, blkn, n, nrhs)                                   &
!$omp firstprivate(ldy, ldx, node, blk, ndep)                             &
!$omp firstprivate(p_lcol, p_y, p_xlocal)                                 &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                        &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(i, threadID, flops)                                         &
!$omp depend(inout: p_bc(blk))
