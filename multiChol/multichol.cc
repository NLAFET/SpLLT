#include "multichol.hpp"
#include "SimdVec.hpp"

#include <algorithm>

void multichol_interleave(int num_vec, int m, int n, double* const a, int lda, int* const info) {
#if 1
   typedef double T;
   // AVX version
   const int veclen = SimdVec<T>::vector_length;
   const int sublda = veclen*lda;
   for(int set=0; set<num_vec; ++set) {
      T* const aptr = &a[set*veclen*lda*n];
#if 0
      if(set+1<num_vec) {
         T* const aptr2 = &a[(set+1)*veclen*lda*n];
         // Prefetch based on 8 doubles per cache line
         // NB L1 cache typically 32K, 8-way associative
         //    So can support approx two matrices of total size 512 each
         //    or about 20x20.
         // i.e. plenty of room???
         for(int i=0; i<m*n; i+=8)
            _mm_prefetch(&aptr2[i], _MM_HINT_T0);
      }
#endif
      for(int col=0; col<n; ++col) {
         // Take square root of diagonal (or error if non-positive)
         SimdVec<T> diag = SimdVec<T>::load_aligned(&aptr[col*(sublda+veclen)]);
         if(any(diag < SimdVec<T>::zero())) {
            // A diagonal entry is less than zero: flag in info, but apply ops
            // regardless
            for(int i=0; i<veclen; ++i) {
               int &flag = info[set*veclen+i];
               if(diag[i] < 0.0) flag = (flag!=0) ? flag : col;
            }
         }
         diag = sqrt(diag);
         SimdVec<T> rdiag = 1.0 / diag;
         diag.store_aligned(&aptr[col*(sublda+veclen)]);
         // Divide column by a
         for(int row=col+1; row<m; ++row) {
            SimdVec<T> rval = SimdVec<T>::load_aligned(
                  &aptr[col*sublda + row*veclen]
                  );
            rval *= rdiag;
            rval.store_aligned( &aptr[col*sublda + row*veclen] );
         }
         // Apply outer product update to remaining columns
         for(int upd=col+1; upd<n; ++upd) {
            // NB: Aim to do a_ij -= a_ki * a_kj
            //     i = row
            //     j = upd
            //     k = col
            SimdVec<T> a_kj = SimdVec<T>::load_aligned(
                  &aptr[col*sublda+upd*veclen]
                  );
            for(int row=upd; row<n; ++row) {
               SimdVec<T> a_ij = SimdVec<T>::load_aligned(
                     &aptr[upd*sublda + row*veclen]
                     );
               SimdVec<T> a_ki = SimdVec<T>::load_aligned(
                     &aptr[col*sublda+row*veclen]
                     );
               a_ij = fmsub(a_ij, a_ki, a_kj);
               a_ij.store_aligned(
                     &aptr[upd*sublda + row*veclen]
                     );
            }
         }
      }
   }
#else
   // Simplistic version
   const int veclen = SimdVec<double>::vector_length;
   const int sublda = veclen*lda;
   for(int set=0; set<num_vec; ++set) {
      for(int prblm=0; prblm<veclen; ++prblm) {
         double* const aprblm = &a[set*veclen*lda*n+prblm];
         for(int col=0; col<n; ++col) {
            // Take square root of diagonal (or error if non-positive)
            double adiag = aprblm[col*(sublda+veclen)];
            if(adiag < 0.0) {
               // Not positive definite
               info[prblm] = col;
               break;
            }
            adiag = sqrt(adiag);
            aprblm[col*(sublda+veclen)] = adiag; // Store back to a
            // Divide column by a
            for(int row=col+1; row<m; ++row)
               aprblm[col*sublda + row*veclen] /= adiag;
            // Apply outer product update to remaining columns
            for(int upd=col+1; upd<n; ++upd)
            for(int row=upd; row<n; ++row)
               aprblm[upd*sublda + row*veclen] -=
                  aprblm[col*sublda+upd*veclen] * aprblm[col*sublda+row*veclen];
         }
      }
   }
#endif
}

template<typename T>
void gather_chol(int n, double* const a, int ldmat, int lda, typename SimdVec<T>::simd_index_type const& gatheridx, int* info, int info_offset) {
   // NB: ldmat is lda for single mat, whilst lda is for block of mats (=N*lda)
   const int veclen = SimdVec<T>::vector_length;
   for(int col=0; col<n; ++col) {
      // Take square root of diagonal (or error if non-positive)
      SimdVec<T> diag = SimdVec<T>::gather(
               &a[col*(ldmat+1)], gatheridx, 1
               );
      if(any(diag < SimdVec<T>::zero())) {
         // A diagonal entry is less than zero: flag in info, but apply ops
         // regardless
         for(int i=0; i<veclen; ++i) {
            int &flag = info[i];
            if(diag[i] < 0.0) flag = (flag!=0) ? flag : info_offset + col;
         }
      }
      diag = sqrt(diag);
      SimdVec<T> rdiag = 1.0 / diag;
      diag.scatter(&a[col*(ldmat+1)], gatheridx, 1);
      // Divide column by a
      for(int row=col+1; row<n; ++row) {
         SimdVec<T> rval = SimdVec<T>::gather(
               &a[col*ldmat + row], gatheridx, 1
               );
         rval *= rdiag;
         rval.scatter(&a[col*ldmat+row], gatheridx, 1);
      }
      // Apply outer product update to remaining columns
      for(int upd=col+1; upd<n; ++upd) {
         // NB: Aim to do a_ij -= a_ik * a_jk
         //     i = row
         //     j = upd
         //     k = col
         SimdVec<T> a_kj = SimdVec<T>::gather(
               &a[col*ldmat+upd], gatheridx, 1
               );
         for(int row=upd; row<n; ++row) {
            SimdVec<T> a_ij = SimdVec<T>::gather(
               &a[upd*ldmat+row], gatheridx, 1
               );
            SimdVec<T> a_ki = SimdVec<T>::gather(
               &a[col*ldmat+row], gatheridx, 1
               );
            a_ij = fmsub(a_ij, a_ki, a_kj);
            a_ij.scatter(
               &a[upd*ldmat+row], gatheridx, 1
               );
         }
      }
   }
}

void print_mat(int n, double* const a, int lda) {
   for(int i=0; i<n; i++) {
      printf("%d:", i);
      for(int j=0; j<n; j++)
         printf(" %e", a[j*lda+i]);
      printf("\n");
   }
}

// #define BLOCKED

void multichol_standard(int num_vec, int m, int n, double* const a, int lda, int* const info) {
#if 1
   typedef double T;
   // AVX version
   const int veclen = SimdVec<T>::vector_length;
   SimdVec<T>::simd_index_type gatheridx;
   //for(int i=0; i<veclen; i++) gatheridx.set_element(i, i*n*lda);
   gatheridx.set_element(0, 0*n*lda);
   gatheridx.set_element(1, 1*n*lda);
   gatheridx.set_element(2, 2*n*lda);
   gatheridx.set_element(3, 3*n*lda);
   for(int set=0; set<num_vec; ++set) {
      double* const aptr = &a[set*veclen*lda*n];
#if defined(BLOCKED)
      for(int bcol=0; bcol<n; bcol+=veclen) {
         gather_chol<T>(veclen, &aptr[bcol*(lda+1)], lda, n*lda, gatheridx, &info[set*veclen], bcol);
         for(int col=bcol; col<bcol+veclen; ++col) {
            for(int prblm=set*veclen; prblm<set*veclen+veclen; ++prblm) {
               double* const aprblm = &a[prblm*lda*n];
               // Load diagonal
               SimdVec<T> diag( 1.0 / aprblm[col*(lda+1)] );
               // Divide column by a
               for(int row=bcol+veclen; row<m; row+=veclen) {
                  SimdVec<T> cvec = SimdVec<T>::load_aligned(
                        &aprblm[col*lda + row]
                        );
                  cvec *= diag;
                  cvec.store_aligned(
                        &aprblm[col*lda + row]
                        );
               }
               // Apply outer product update to remaining columns
               // First handle first veclen columns with rfrom=bcol+veclen
               for(int upd=col+1; upd<bcol+veclen; ++upd) {
                  // Update a_:j -= a_jk * a_:k
                  //     : = row range (represented by i in naming)
                  //     j = upd
                  //     k = col
                  SimdVec<T> a_jk( aprblm[col*lda + upd] );
                  // Start from first unfactorized row block
                  int rfrom = bcol+veclen;
                  for(int row=rfrom; row<n; row+=veclen) {
                     SimdVec<T> a_ij = SimdVec<T>::load_aligned(
                           &aprblm[upd*lda + row]
                           );
                     SimdVec<T> a_ik = SimdVec<T>::load_aligned(
                           &aprblm[col*lda + row]
                           );
                     a_ij = fmsub(a_ij, a_ik, a_jk);
                     a_ij.store_aligned(
                           &aprblm[upd*lda + row]
                           );
                  }
               }
               // Handle remaining columns with rfrom based on upd value
               for(int upd=bcol+veclen; upd<n; ++upd) {
                  // Update a_:j -= a_jk * a_:k
                  //     : = row range (represented by i in naming)
                  //     j = upd
                  //     k = col
                  SimdVec<T> a_jk( aprblm[col*lda + upd] );
                  // Start from aligned location at or just below upd
                  int rfrom = veclen * (upd/veclen);
                  for(int row=rfrom; row<n; row+=veclen) {
                     SimdVec<T> a_ij = SimdVec<T>::load_aligned(
                           &aprblm[upd*lda + row]
                           );
                     SimdVec<T> a_ik = SimdVec<T>::load_aligned(
                           &aprblm[col*lda + row]
                           );
                     a_ij = fmsub(a_ij, a_ik, a_jk);
                     a_ij.store_aligned(
                           &aprblm[upd*lda + row]
                           );
                  }
               }
            }
         }
      }
#else
   gather_chol<T>(n, aptr, lda, n*lda, gatheridx, &info[set*veclen], 0);
#endif
   }
#else
   // Simplistic version
   const int veclen = SimdVec<double>::vector_length;
   for(int prblm=0; prblm<num_vec*veclen; ++prblm) {
      double* const aprblm = &a[prblm*lda*n];
      for(int col=0; col<n; ++col) {
         // Take square root of diagonal (or error if non-positive)
         double adiag = aprblm[col*(lda+1)];
         if(adiag < 0.0) {
            // Not positive definite
            info[prblm] = col;
            break;
         }
         adiag = sqrt(adiag);
         aprblm[col*(lda+1)] = adiag; // Store back to a
         // Divide column by a
         for(int row=col+1; row<m; ++row)
            aprblm[col*lda + row] /= adiag;
         // Apply outer product update to remaining columns
         for(int upd=col+1; upd<n; ++upd)
         for(int row=upd; row<n; ++row)
            aprblm[upd*lda + row] -=
               aprblm[col*lda+upd] * aprblm[col*lda+row];
      }
   }
#endif
}
