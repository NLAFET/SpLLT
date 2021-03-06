!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols
!$omp default(none)                                           &
!$omp shared(fkeep)                                           &
!$omp firstprivate(blkm, blkn, nrhs, n, ldy, node, ndep)      &
!$omp firstprivate(trace_id, p_timer)                         &
!$omp firstprivate(p_y, p_lcol, p_index, p_rhs, p_order)      &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp firstprivate(task_manager)                              &
!$omp private(threadID, i, flops)                             &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))
