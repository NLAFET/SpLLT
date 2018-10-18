!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols




                                                                                




                                                                                   


                                                                                

                                                                                   

                                                                                   


                                                                                   


                                                                                   

                                                                                   

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(0)',1, task_manager%workerID, timer)
# endif
!$omp task &
!$omp default(none)                                                       &
!$omp shared(fkeep)                                                       &
!$omp firstprivate(blkm, blkn, nrhs, ldy, ndep, nwdep, node)              &
!$omp firstprivate(trace_id, p_timer)                                     &
!$omp firstprivate(p_y, p_rhs, p_index, p_order, p_lcol)                  &
!$omp firstprivate(p_bc, p_dep, p_wdep)                                   &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(i, threadID, flops)                                         &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

# include "spllt_solve_fwd_block_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(1, task_manager%workerID, timer)
# endif

