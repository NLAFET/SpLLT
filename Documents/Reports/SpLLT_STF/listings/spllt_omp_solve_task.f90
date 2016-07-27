!$omp task firstprivate(m, n)
!$omp    & depend(in:bc_kk%c) &
!$omp    & depend(inout:bc_ik%c)

call spllt_solve_block(m, n, bc_kk%c, bc_ik%c)

!$omp end task



