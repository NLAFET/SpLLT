



                                                                                




                                                                                   


                                                                                

                                                                                   

                                                                                   


                                                                                   


                                                                                   

                                                                                   

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(0)',1, task_manager%workerID, timer)
# endif
!$omp task &
!!$omp default(none)                              &
!$omp firstprivate(p_rhs)                         &
!$omp firstprivate(nrhs, n)                       &
!$omp firstprivate(sub_task_manager)              &
!$omp firstprivate(p_bc, blk_en)                  &
!$omp firstprivate(en, sa, ndep_lvl, chunk)       &
!$omp shared(fkeep)                               &
!$omp depend(inout: p_bc(blk_en))                 &
!$omp private(i, trace_id)

# include "spllt_solve_bwd_node_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(1, task_manager%workerID, timer)
# endif

