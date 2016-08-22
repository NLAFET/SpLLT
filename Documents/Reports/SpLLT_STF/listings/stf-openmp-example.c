#pragma omp parallel
{
#pragma omp master
   {
      for (i = 1; i < N; i++) {
#pragma omp task depend(inout:x[i:1])
         x[i] = f(x[i:1]);
#pragma omp task depend(in:x[i], y[i-1:1]) depend(out:y[i:1])
         y[i] = g(x[i:1], y[i-1:1]);
      }
#pragma omp taskwait
   }
}
