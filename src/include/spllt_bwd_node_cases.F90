!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Sebastien Cayrols




                                                                                




                                                                                   


                                                                                

                                                                                   

                                                                                   


                                                                                   


                                                                                   

                                                                                   

case(1)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(1)',2, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)))   &
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
call spllt_tac(2, task_manager%workerID, timer)
# endif

case(2)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(2)',3, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
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
call spllt_tac(3, task_manager%workerID, timer)
# endif

case(3)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(3)',4, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)))&
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
call spllt_tac(4, task_manager%workerID, timer)
# endif

case(4)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(4)',5, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
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
call spllt_tac(5, task_manager%workerID, timer)
# endif

case(5)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(5)',6, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)))&
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
call spllt_tac(6, task_manager%workerID, timer)
# endif

case(6)

# if defined(SPLLT_TIMER_TASKS_SUBMISSION)
call spllt_tic('CASE(6)',7, task_manager%workerID, timer)
# endif
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
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
call spllt_tac(11, task_manager%workerID, timer)
# endif

