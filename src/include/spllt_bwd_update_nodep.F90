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
!$omp firstprivate(blkm, blkn, nrhs, n, ldy, ldx, node, nwdep)            &
!$omp firstprivate(p_lcol, p_y, p_xlocal)                                 &
!$omp firstprivate(p_bc, p_dep, p_wdep, trace_id, p_timer)                &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(dblk, nbcol)                                           &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, i, flops)                                         &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(1, task_manager%workerID, timer)
# endif

