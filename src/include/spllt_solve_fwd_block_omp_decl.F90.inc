!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)                      &
!$omp firstprivate(nthread, trace_id)                                     &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep)                                           &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(i, j, r, threadID)                                          &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))
