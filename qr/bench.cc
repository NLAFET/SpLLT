#include <cstdio>
#include <stdlib.h>
#include <chrono>

#include "problem.hh"
#include "qr.hh"

int main(void) {

   int M, N, LDA;

   M = 24;
   N = 12;
   LDA = M;
   
   size_t sz = LDA*N;
   double *a = (double *)aligned_alloc(32, sz*sizeof(double));  
   
   Problem<double> *pbl = new Problem<double>(M, N, a, LDA);

   // pbl->print();
   
   printf("[bench] starting tests\n");

   printf("[bench] facto\n");
   
   double *r = (double *)aligned_alloc(32, N*N*sizeof(double));  
   pbl->r = r;
   pbl->ldr = N;

   auto begin = std::chrono::high_resolution_clock::now();   
   qr_mgs(pbl->m, pbl->n, pbl->a, pbl->lda, pbl->r, pbl->ldr);
   auto end = std::chrono::high_resolution_clock::now();
   long ttotal = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
   printf("[bench] facto took %ld ns\n", ttotal);

   pbl->q = a;
   pbl->ldq = pbl->lda;

   printf("[bench] check\n");
   
   double bwderr_lsq = pbl->bwderr_lsq();
   double bwderr = pbl->bwderr();

   printf("[bench] ||A^Tr|| / ||r||                = %e\n", bwderr_lsq);
   printf("[bench] ||Ax-b|| / (||A||||x|| + ||b||) = %e\n", bwderr);

   // if(bwderr_lsq >= 1e-14 || bwderr_lsq!=bwderr_lsq)
   //    printf("[bench] Failed bwderr_lsq = %e\n", bwderr_lsq);

   // if(bwderr >= 1e-14 || bwderr!=bwderr)
   //    printf("[bench] Failed bwderr = %e\n", bwderr);

 
   printf("[bench] done\n");

   return 0;
}
