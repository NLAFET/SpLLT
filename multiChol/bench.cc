#include "AlignedCallocBlock.hpp"
#include "cpu_wrappers.hpp"
#include "multichol.hpp"
#include "SimdVec.hpp"

#include <algorithm>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <cstdio>

template <typename T>
void print_interleave_mat(int m, int n, T const* a, int lda, int ida) {
   printf("Printing %p %d\n", a, lda);
   for(int row=0; row<m; ++row) {
      printf("%d:", row);
      for(int col=0; col<n; ++col) {
         printf(" %e", a[col*lda+row*ida]);
      }
      printf("\n");
   }
}

template <typename T>
void print_mat(int m, int n, T const* a, int lda) {
   print_interleave_mat(m, n, a, lda, 1);
}

template <typename T>
void print_vec(char const* name, int n, T const* a) {
   if(name) printf("%s", name);
   for(int i=0; i<n; ++i)
      printf(" %e", a[i]);
   printf("\n");
}

void to_interleave(int m, int n, double const* src, int lds, double* dest, int ldd, int idd) {
   for(int col=0; col<n; ++col)
   for(int row=0; row<m; ++row)
      dest[col*ldd+row*idd] = src[col*lds+row];
}
void from_interleave(int m, int n, double const* src, int lds, int ids, double* dest, int ldd) {
   for(int col=0; col<n; ++col)
   for(int row=0; row<m; ++row)
      dest[col*ldd+row] = src[col*lds+row*ids];
}

double max_absval(int n, double const* vec) {
   double v=0.0;
   for(int i=0; i<n; ++i)
      v = std::max(v, fabs(vec[i]));
   return v;
}

template <typename T>
class Problem {
public:
   int const m;
   int const n;
   /* Interleved */
   T *const ai;
   int const ldai;
   int const idai;
   /* Standard */
   T *const as;
   int const ldas;

   Problem(int m, int n, T *ai, int ldai, int idai, T *as, int ldas)
   : m(m), n(n), ai(ai), ldai(ldai), idai(idai), as(as), ldas(ldas), aorig_(nullptr), rhs_(nullptr)
   {
      // Initialise a to random matrix (and upper triangle to NaN)
      for(int col=0; col<n; ++col) {
         for(int row=0; row<col; ++row)
            ai[col*ldai+row*idai] = std::numeric_limits<T>::quiet_NaN();
         for(int row=col; row<m; ++row)
            ai[col*ldai+row*idai] = 2 * ((T) std::rand()) / RAND_MAX - 1.0;
      }
      // Make diagonally dominant
      for(int col=0; col<n; ++col) {
         ai[col*(ldai+idai)] = 1.0 + fabs(ai[col*(ldai+idai)]);
         for(int row=0; row<col; ++row)
            ai[col*(ldai+idai)] += fabs( ai[row*ldai+col*idai] );
         for(int row=col+1; row<m; ++row)
            ai[col*(ldai+idai)] += fabs( ai[col*ldai+row*idai] );
      }
      // Copy ai to as
      from_interleave(m, n, ai, ldai, idai, as, ldas);
      // Take a copy for verficiation purposes
      aorig_ = new double[m*n];
      from_interleave(m, n, ai, ldai, idai, aorig_, m);

      // Create a rhs corresponding to [0.1, 0.2, ... ]
      T *x = new T[n];
      for(int i=0; i<n; ++i) x[i] = 0.1*(i+1);
      rhs_ = new T[m];
      spmv(x, rhs_);
      delete[] x;
   }
   ~Problem()
   {
      delete[] rhs_;
      delete[] aorig_;
   }

   /** Calculate backward error from pre-calculated rhs */
   T bwderr_interleve() {
      /* Perform solve */
      double *soln = new double[m];
      for(int i=0; i<m; ++i) soln[i] = rhs_[i];
      //print_vec("rhs =", m, rhs_);
      double *lval = new double[m*n];
      from_interleave(m, n, ai, ldai, idai, lval, m);
      trsv<double>('L', 'N', 'N', n, lval, m, soln, 1);
      trsv<double>('L', 'T', 'N', n, lval, m, soln, 1);
      //print_vec("soln =", m, soln);
      delete[] lval;
      
      /* Evaluate residual */
      double *resid = new T[m];
      spmv(soln, resid);
      for(int i=0; i<m; i++) resid[i] -= rhs_[i];

      /* Calculate scaled norm || Ax-b || / ( || A || || x || + || b || ) */
      double resid_inf = max_absval(m, resid);
      double rhs_inf = max_absval(m, rhs_);
      double A_inf = calc_ainf();
      double soln_inf =  max_absval(m, soln);
      double scaled_bwderr = resid_inf / (A_inf*soln_inf + rhs_inf);

      /* Free memory */
      delete[] resid;
      delete[] soln;

      return scaled_bwderr;
   }

   /** Calculate backward error from pre-calculated rhs */
   T bwderr_standard() {
      /* Perform solve */
      double *soln = new double[m];
      for(int i=0; i<m; ++i) soln[i] = rhs_[i];
      //print_vec("rhs =", m, rhs_);
      trsv<double>('L', 'N', 'N', n, as, ldas, soln, 1);
      trsv<double>('L', 'T', 'N', n, as, ldas, soln, 1);
      //print_vec("soln =", m, soln);
      
      /* Evaluate residual */
      double *resid = new T[m];
      spmv(soln, resid);
      for(int i=0; i<m; i++) resid[i] -= rhs_[i];

      /* Calculate scaled norm || Ax-b || / ( || A || || x || + || b || ) */
      double resid_inf = max_absval(m, resid);
      double rhs_inf = max_absval(m, rhs_);
      double A_inf = calc_ainf();
      double soln_inf =  max_absval(m, soln);
      double scaled_bwderr = resid_inf / (A_inf*soln_inf + rhs_inf);

      /* Free memory */
      delete[] resid;
      delete[] soln;

      return scaled_bwderr;
   }

private:
   void spmv(T const* x, T* y) {
      for(int row=0; row<m; ++row) y[row] = 0.0;
      for(int col=0; col<n; ++col) {
         y[col] += aorig_[col*m+col] * x[col];
         for(int row=col+1; row<m; ++row) {
            y[col] += aorig_[col*m+row] * x[row];
            y[row] += aorig_[col*m+row] * x[col];
         }
      }
   }

   /* Returns || A ||_inf, being the maxmimum row sum of abs values */
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

   T *aorig_;
   T *rhs_;
};

int main(void) {
   const int NCALL = 200000;
   const int M = 12;
   const int LDA = M;
   const int N = M;
   const int veclen = SimdVec<double>::vector_length;
   assert( M >= N );
   assert( M == N ); // FIXME: Remove and fix up code to check soln
   assert( LDA >= M );

   //
   // Generate sufficient diag dominant matrices
   //
   // Allocate mem
   size_t sz = NCALL*veclen*LDA*N;
   AlignedCallocBlock<double> amemi(sz);
   AlignedCallocBlock<double> amems(sz);
   double* const ai = amemi.get_ptr();
   double* const as = amems.get_ptr();
   int* const info = new int[NCALL*veclen];
   std::vector< Problem<double> > problems;
   problems.reserve(NCALL*veclen);
   // Create problems
   double* aptri = ai;
   double* aptrs = as;
   for(int i=0; i<NCALL; ++i) {
      for(int j=0; j<veclen; ++j) {
         problems.emplace_back( M, N, aptri+j, veclen*LDA, veclen, aptrs, LDA );
         aptrs += N*LDA;
      }
      aptri += N*veclen*LDA;
   }

   //
   // Interleave: Call code and time it
   //
   auto begin = std::chrono::high_resolution_clock::now();
#if 1
   multichol_interleave(NCALL, M, N, ai, LDA, info);
#else
   int idx=0;
   double *lval = new double[M*N];
   for(auto problem=problems.begin(); problem!=problems.end(); ++problem, ++idx) {
      //printf("A:\n");
      //print_interleave_mat<double>(problem->m, problem->n, problem->a, problem->lda, problem->ida);
      from_interleave(problem->m, problem->n, problem->a, problem->lda, problem->ida, lval, M);
      info[idx] = potrf<double>('L', problem->n, lval, M);
      to_interleave(problem->m, problem->n, lval, M, problem->a, problem->lda, problem->ida);
      //printf("L:\n");
      //print_interleave_mat<double>(problem->m, problem->n, problem->a, problem->lda, problem->ida);
      if(info[idx]!=0) {
         printf("info[%d] = %d\n", idx, info[idx]);
         exit(1);
      }
   }
   delete[] lval;
#endif
   auto end =  std::chrono::high_resolution_clock::now();
   long ttotal = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
   printf("Interleave took %e sec total for %d calls of veclen %d\n", 1e-9*ttotal, NCALL, veclen);
   printf("Interleave per matrix: %ld ns\n", ttotal/(NCALL*veclen));
   // Verify solutions
   for(auto problem=problems.begin(); problem!=problems.end(); ++problem) {
      double bwderr = problem->bwderr_interleve();
      if(bwderr >= 1e-14 || bwderr!=bwderr)
         printf("Failed bwderr = %e\n", bwderr);
   }

   //
   // Standard: Call code and time it
   //
   begin = std::chrono::high_resolution_clock::now();
#if 1
   multichol_standard(NCALL, M, N, as, LDA, info);
#else
   int idx=0;
   for(auto problem=problems.begin(); problem!=problems.end(); ++problem, ++idx) {
      info[idx] = potrf<double>('L', problem->n, problem->as, problem->ldas);
      /*if(info[idx]!=0) {
         printf("info[%d] = %d\n", idx, info[idx]);
         exit(1);
      }*/
   }
#endif
   end =  std::chrono::high_resolution_clock::now();
   ttotal = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
   printf("Standard took %e sec total for %d calls of veclen %d\n", 1e-9*ttotal, NCALL, veclen);
   printf("Standard per matrix: %ld ns\n", ttotal/(NCALL*veclen));
   // Verify solutions
   int idx = 0;
   for(auto problem=problems.begin(); problem!=problems.end(); ++problem, ++idx) {
      double bwderr = problem->bwderr_standard();
      if(bwderr >= 1e-14 || bwderr!=bwderr)
         printf("Failed bwderr for problem %d = %e\n", idx, bwderr);
   }

   //
   // Cleanup memory
   //
   delete[] info;

   return 0;
}
