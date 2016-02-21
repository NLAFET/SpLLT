#ifndef CPU_WRAPPERS_H
#define CPU_WRAPPERS_H

extern "C" {
   void dtrsv_(char *uplo, char *trans, char *diag, int *n, const double *a, int *lda, double *x, int *incx);
   void dgemv_(char *trans, int *m, int *n, double *alpha, const double *a, int *lda, const double *x, int *incx, double *beta, double *y, int *incy);
   void dlagge_(int *m, int *n, int *kl, int *ku, double* d, double* a, int *lda, int* iseed, double* work, int* info);
   void dlatms_(int *m, int *n, char *dist, int *iseed, char *sym, double *d, int *mode, double *cond, double *dmax, int *kl, int *ku, char *pack, double *a, int *lda, double *work, int *info);
}

/* _TRSV */
template <typename T>
void trsv(char uplo, char trans, char diag, int n, const T* a, int lda, T* x, int incx);
template <>
void trsv<double>(char uplo, char trans, char diag, int n, const double* a, int lda, double* x, int incx) {
   dtrsv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

/* _GEMV */
// y = alpha.A*x + beta.y
template <typename T>
void gemv(char trans, int m, int n, double alpha, const T* a, int lda, const T* x, int incx, double beta, T* y, int incy);
template <>
void gemv<double>(char trans, int m, int n, double alpha, const double* a, int lda, const double* x, int incx, double beta, double* y, int incy) {
   dgemv_(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

/* _DLAGGE */
template <typename T>
void lagge(int m, int n, int kl, int ku, T* d, T* a, int lda, int* iseed, T* work, int* info);
template <>
void lagge<double>(int m, int n, int kl, int ku, double* d, double* a, int lda, int* iseed, double* work, int* info) {
   dlagge_(&m, &n, &kl, &ku, d, a, &lda, iseed, work, info);
   // dlagge_(&m, &n, &kl, &ku, d, a, &lda, iseed, work, info);
}

/* _DLATMS */
template <typename T>
void latms(int m, int n, char dist, int *iseed, char sym, T *d, int mode, T cond, T dmax, int kl, int ku, char pack, T *a, int lda, T *work, int *info);
template <>
void latms<double>(int m, int n, char dist, int *iseed, char sym, double *d, int mode, double cond, double dmax, int kl, int ku, char pack, double *a, int lda, double *work, int *info) {
   dlatms_(&m, &n, &dist, iseed, &sym, d, &mode, &cond, &dmax, &kl, &ku, &pack, a, &lda, work, info);
}
#endif // CPU_WRAPPERS_H
