



                                                                                




                                                                                   


                                                                                

                                                                                   

                                                                                   


                                                                                   


                                                                                   

                                                                                   

case(1)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(1)',2, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)))   &
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(2, task_manager%workerID, timer)
# endif

case(2)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(2)',3, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(3, task_manager%workerID, timer)
# endif

case(3)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(3)',4, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)))&
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(4, task_manager%workerID, timer)
# endif

case(4)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(4)',5, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(5, task_manager%workerID, timer)
# endif

case(5)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(5)',6, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)))&
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(6, task_manager%workerID, timer)
# endif

case(6)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(6)',7, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(7, task_manager%workerID, timer)
# endif

case(7)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(7)',8, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)))&
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(8, task_manager%workerID, timer)
# endif

case(8)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(8)',9, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(9, task_manager%workerID, timer)
# endif

case(9)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(9)',10, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+8)))&
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(10, task_manager%workerID, timer)
# endif

case(10)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(10)',11, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+8)),p_bc(p_dep(alpha*2+beta+8)))&
!$omp default(none)                                                       &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id, p_timer)                          &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp firstprivate(task_manager)                                          &
!$omp private(threadID, flops)                                            &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tac(11, task_manager%workerID, timer)
# endif

