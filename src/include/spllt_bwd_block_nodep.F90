



                                                                                




                                                                                   


                                                                                

                                                                                   

                                                                                   


                                                                                   


                                                                                   

                                                                                   

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(0)',1, task_manager%workerID, timer)
# endif
!$omp task &
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

# include "spllt_solve_bwd_block_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(1, task_manager%workerID, timer)
# endif

