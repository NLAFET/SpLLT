!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))