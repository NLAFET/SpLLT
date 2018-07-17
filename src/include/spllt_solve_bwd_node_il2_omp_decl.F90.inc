!!$omp default(none)                              &
!$omp firstprivate(p_rhs)                         &
!$omp firstprivate(nrhs, n)                       &
!$omp firstprivate(sub_task_manager)              &
!$omp firstprivate(p_bc, blk_en)                  &
!$omp firstprivate(en, sa)                        &
!$omp shared(fkeep)                               &
!$omp depend(inout: p_bc(blk_en))                 &
!$omp private(i, trace_id)
