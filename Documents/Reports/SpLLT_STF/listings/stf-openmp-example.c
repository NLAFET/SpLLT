#pragma omp parallel
{
#pragma omp master
   {
      for (i = 1; i < N; i++) {
#pragma omp task depend(inout:x[i])
         x[i] = f(x[i]);
#pragma omp task depend(in:x[i], y[i-1]) depend(out:y[i])
         y[i] = g(x[i], y[i-1]);
      }
#pragma omp taskwait
   }
}
