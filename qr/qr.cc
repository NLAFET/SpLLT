#include "qr.hh"

double norm(int n, double const *a) {
   
   double nrm = 0.0;
   
   for (int i=0; i<n; ++i) {
      nrm += a[i]*a[i];
   }

   return sqrt(nrm);
}

// modified Gram-Schmidt
void qr_mgs(int m, int n, double *a, int lda, double *r, int ldr) {

   for (int j=0; j<n; ++j) {

      for (int i = 0; i<j; ++i) {
         
         for (int k=0; k<m; ++k) r[i+ldr*j] += a[k+lda*i]*a[k+lda*j];
         for (int k=0; k<m; ++k) a[k+lda*j] -= r[i+ldr*j]*a[k+lda*i];
      }
      
      r[j+ldr*j] = norm(m, a + j*lda);

      for (int k=0; k<m; ++k) a[k+lda*j] /= r[j+ldr*j];
   }
}

// modified Gram-Schmidt
void qr_mgs_vec(int m, int n, double *a, int lda, double *r, int ldr) {

   for (int j=0; j<n; ++j) {

      for (int i = 0; i<j; ++i) {
         
         for (int k=0; k<m; ++k) r[i+ldr*j] += a[k+lda*i]*a[k+lda*j];
         for (int k=0; k<m; ++k) a[k+lda*j] -= r[i+ldr*j]*a[k+lda*i];
      }
      
      r[j+ldr*j] = norm(m, a + j*lda);

      for (int k=0; k<m; ++k) a[k+lda*j] /= r[j+ldr*j];
   }
}

// void qr(int m, int n, double const *a, int lda) {

//    double r = 0.0;

//    for (int j=0; j<n; ++j) {
//       r = norm(m-j, a + j*lda);

      
//    }

// }
