#ifndef PROBLEM_H
#define PROBLEM_H

#include <algorithm>
#include <cstdlib>
#include <limits>

#include "cpu_wrappers.hh"

template <typename T>
class Problem {
public:
   int m;
   int n;
   T * a;
   int lda;
   
   // factors. dedicated class? 
   T *q;
   int ldq;
   T *r;
   int ldr;   

   Problem(int m, int n, T *a, int lda)
      : m(m), n(n), a(a), lda(lda),
        q(nullptr), ldq(0), r(nullptr), ldr(0),
        aorig_(nullptr), rhs_(nullptr),
        d_(nullptr){

      // array containing singular/eigen values
      d_ = new T[std::min(m,n)];  
      
      // generate random pbl
      // gen_rand_pbl();
      gen_rand_rr_pbl();

      // Take a copy for verficiation purposes
      ldaorig_ = m;
      aorig_ = new double[ldaorig_*n];
      cpy_orig();

      // Create a rhs
      rhs_ = new T[m];
      // corresponding to [0.1, 0.2, ... ]
      // for(int i=0; i<m; ++i) rhs_[i] = 0.1*(i+1);
      // rhs with random values
      gen_rand_rhs();

   }

   ~Problem()
   {
      delete[] rhs_;
      delete[] aorig_;
   }
   
   void gen_rand_rhs() {

      for(int i=0; i<m; ++i) rhs_[i] = 2 * ((T) std::rand()) / RAND_MAX - 1.0;
   }

   T bwderr() {
      double bwderr = std::numeric_limits<T>::quiet_NaN();

      // r and q should be allocated
      if (r && q) {
         
         double *x = new double[n];
      
         // y = Q^T.b
         gemv<double>('T', m, n, 1.0, q, ldq, rhs_, 1, 0.0, x, 1);
         // printf("||Q^T.b||  : %f\n", max_absval(n, x));

         // solve Rx = y
         trsv<double>('U', 'N', 'N', n, r, ldr, x, 1);
         // printf("||x||  : %f\n", max_absval(n, x));

         double *resid = new T[m];
         for(int i=0; i<m; i++) resid[i] = rhs_[i];

         gemv<double>('N', m, n, 1.0, aorig_, ldaorig_, x, 1, -1.0, resid, 1);

         // ||r||
         double resid_inf = max_absval(m, resid);
         double rhs_inf = max_absval(m, rhs_);
         double A_inf = calc_ainf();
         double soln_inf =  max_absval(m, x);
         bwderr = resid_inf / (A_inf*soln_inf + rhs_inf);      

         delete[] resid;
         delete[] x;
      }
      
      return bwderr;
   }

   T bwderr_lsq() {
      
      double bwderr = std::numeric_limits<T>::quiet_NaN();

      // r and q should be allocated
      if (r && q) {
         
         double *x = new double[n];
      
         // y = Q^T.b
         gemv<double>('T', m, n, 1.0, q, ldq, rhs_, 1, 0.0, x, 1);
         printf("m  : %d, n: %d\n", m, n);
         printf("||Q||  : %e\n", calc_ainf());
         printf("||b||  : %e\n", max_absval(m, rhs_));
         printf("||Q^T.b||  : %e\n", max_absval(n, x));

         // solve Rx = y
         trsv<double>('U', 'N', 'N', n, r, ldr, x, 1);
         // printf("||x||  : %f\n", max_absval(n, x));

         double *resid = new T[m];
         for(int i=0; i<m; i++) resid[i] = rhs_[i];

         gemv<double>('N', m, n, 1.0, aorig_, ldaorig_, x, 1, -1.0, resid, 1);
         
         // ||r||
         double normresid = max_absval(m, resid);
         // ||A^Tr||
         double *proj = new double[n];
         gemv<double>('T', m, n, 1.0, a, lda, resid, 1, 0.0, proj, 1);         
         double normproj = max_absval(n, proj);

         // printf("normproj  : %e\n", normproj);
         // printf("normresid : %e\n", normresid);
         
         bwderr = normproj/normresid;

         // free memory
         delete[] proj;
         delete[] resid;
         delete[] x;
      }

      return bwderr;
   }

   void gen_rand_pbl() {

      for(int col=0; col<n; ++col) {
         for(int row=0; row<m; ++row)
            a[col*lda+row] = 2 * ((T) std::rand()) / RAND_MAX - 1.0;
      }
   }

   void gen_rand_rr_pbl() {

      int k = std::min(m,n);
      // T *d = new T[k];
      T *work = new T[m+n];
      int iseed[4] = {0,0,0,1};
      int info;

      // first k-l values equal to 1.0
      for(int i=0; i<k; ++i) d_[i] = 1.0;

      // last l values equal to 1.0/1e12
      int l = 5;
      for(int i=0; i<l; ++i) d_[k-1-i] = 1e-9;
      
      
      // lagge<double>(m, n, m-1, n-1, d, a, lda, iseed, work, &info);
      int kl = m-1;
      int ku = n-1;

      dlagge_(&m, &n, &kl, &ku, d_, a, &lda, iseed, work, &info);
      // latms<double>(m, n, 'U', iseed, 'N', );
   }

   void print() {
      for(int i=0; i<m; ++i) {
         for(int j=0; j<n; ++j) {
            printf(" %f ", a[i+lda*j]);
         }
         printf("\n");
      }
   }
   
private:

   T calc_ainf() {
      T *rsum = new T[m];
      for(int row=0; row<m; ++row) rsum[row] = 0.0;
      for(int col=0; col<n; ++col) {
         rsum[col] += fabs(aorig_[col*m+col]);
         for(int row=col+1; row<m; ++row) {
            rsum[row] += fabs(aorig_[col*m+row]);
            rsum[col] += fabs(aorig_[col*m+row]);
         }
      }
      T best = 0.0;
      for(int row=0; row<m; ++row) best = std::max(best, rsum[row]);
      delete[] rsum;
      return best;
   }

   void cpy_orig() {
      
      for(int col=0; col<n; ++col) {
         for(int row=0; row<m; ++row)
            aorig_[col*ldaorig_+row] = a[col*lda+row];
      }
   }

   double max_absval(int n, double const* vec) {
      double v=0.0;
      for(int i=0; i<n; ++i)
         v = std::max(v, fabs(vec[i]));
      return v;
   }

   T *aorig_;
   int ldaorig_;
   T *rhs_;
   // singular/eigenvalues
   T *d_;

};

#endif // PROBLEM_H
