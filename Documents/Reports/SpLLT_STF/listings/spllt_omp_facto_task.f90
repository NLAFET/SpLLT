!$omp task firstprivate(m, n) &
!$omp    & depend(inout:bc%c(1))

call spllt_factorize_block(m, n, bc%c)

!$omp end task






