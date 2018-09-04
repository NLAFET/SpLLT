!$omp default(none)                                                       &
!$omp shared(fkeep)                                                       &
!$omp firstprivate(blkm, blkn, n, nrhs)                                   &
!$omp firstprivate(ldy, ldx, node, blk, ndep, nwdep)                      &
!$omp firstprivate(p_lcol, p_y, p_xlocal)                                 &
!$omp firstprivate(p_bc, p_dep, p_wdep, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager, trace_id)                                &
!$omp private(i, threadID, flops)                                         &
!$omp depend(inout: p_bc(blk))
