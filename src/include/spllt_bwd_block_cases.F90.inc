



                                                                                




                                                                                   


                                                                                

                                                                                   

                                                                                   


                                                                                   


                                                                                   

                                                                                   

case(1)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)))   &
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(2)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(3)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)))&
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(4)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(5)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)))&
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(6)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(7)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)))&
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(8)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(9)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+8)))&
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

case(10)
!$omp task depend(out: p_bc(p_dep(alpha*1+beta))) depend(in: p_bc(p_dep(alpha*1+beta)),p_bc(p_dep(alpha*2+beta)))   &
!$omp depend(in: p_bc(p_dep(alpha*1+beta+2)),p_bc(p_dep(alpha*2+beta+2)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+4)),p_bc(p_dep(alpha*2+beta+4)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+6)),p_bc(p_dep(alpha*2+beta+6)))&
!$omp depend(in: p_bc(p_dep(alpha*1+beta+8)),p_bc(p_dep(alpha*2+beta+8)))&
!$omp firstprivate(m, n, nrhs, col, ldr, sa, offset)          &
!$omp firstprivate(nthread, trace_id)                         &
!$omp firstprivate(p_index, p_lcol, p_xlocal, p_rhs, p_upd)   &
!$omp firstprivate(p_bc, p_dep)                               &
!$omp firstprivate(chunk, ndep_lvl, alpha, beta)              &
!$omp private(threadID, r, j, i)                              &
!$omp firstprivate(dblk) depend(inout: p_bc(dblk))

include '../src/include/spllt_solve_bwd_block_worker.F90.inc'

!$omp end task

