



                                                                                




                                                                                   


                                                                                

                                                                                   

                                                                                   


                                                                                   


                                                                                   

                                                                                   

case(1)
call spllt_tic('CASE(1)',2, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)))   &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(2, scheduler%workerID, timer)
case(2)
call spllt_tic('CASE(2)',3, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(3, scheduler%workerID, timer)
case(3)
call spllt_tic('CASE(3)',4, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)))&
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(4, scheduler%workerID, timer)
case(4)
call spllt_tic('CASE(4)',5, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(5, scheduler%workerID, timer)
case(5)
call spllt_tic('CASE(5)',6, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)))&
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(6, scheduler%workerID, timer)
case(6)
call spllt_tic('CASE(6)',7, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(7, scheduler%workerID, timer)
case(7)
call spllt_tic('CASE(7)',8, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)))&
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(8, scheduler%workerID, timer)
case(8)
call spllt_tic('CASE(8)',9, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(9, scheduler%workerID, timer)
case(9)
call spllt_tic('CASE(9)',10, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+8)))&
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(10, scheduler%workerID, timer)
case(10)
call spllt_tic('CASE(10)',11, scheduler%workerID, timer)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+8)),p_bc(p_dep(alpha*2+beta+8)))&
!$omp firstprivate(m, n, nrhs, col, ldr, blk_sa, offset)                  &
!$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)               &
!$omp firstprivate(p_bc, p_dep, trace_id)                                 &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)                          &
!$omp private(threadID)                                                   &
!$omp firstprivate(blk) depend(inout: p_bc(blk))

# include "spllt_solve_bwd_update_worker.F90.inc"

!$omp end task
call spllt_tac(11, scheduler%workerID, timer)
