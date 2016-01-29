#ifndef CPU_WRAPPERS_H
#define CPU_WRAPPERS_H

extern "C" {
   void dtrsv_(char *uplo, char *trans, char *diag, int *n, const double *a, int *lda, double *x, int *incx);
   void dgemv_(char *trans, int *m, int *n, double *alpha, const double *a, int *lda, const double *x, int *incx, double *beta, double *y, int *incy);
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

#endif // CPU_WRAPPERS_H
