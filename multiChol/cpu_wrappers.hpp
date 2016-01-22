#pragma once

extern "C" {
   void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, const double* a, int* lda, const double* b, int* ldb, double *beta, double* c, int* ldc);
   void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
   void dsytrf_(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
   void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, const double *alpha, const double *a, int *lda, double *b, int *ldb);
   void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, const double *a, int *lda, double *beta, double *c, int *ldc);
   void dtrsv_(char *uplo, char *trans, char *diag, int *n, const double *a, int *lda, double *x, int *incx);
}

/* _GEMM */
template <typename T>
void gemm(char transa, char transb, int m, int n, int k, T alpha, const T* a, int lda, const T* b, int ldb, T beta, T* c, int ldc);
template <>
void gemm<double>(char transa, char transb, int m, int n, int k, double alpha, const double* a, int lda, const double* b, int ldb, double beta, double* c, int ldc) {
   dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

/* _POTRF */
template <typename T>
int potrf(char uplo, int n, T* a, int lda);
template<>
int potrf<double>(char uplo, int n, double* a, int lda) {
   int info;
   dpotrf_(&uplo, &n, a, &lda, &info);
   return info;
}

/* _SYTRF - Bunch-Kaufman factorization */
template <typename T>
int lapack_sytrf(char uplo, int n, T* a, int lda, int* ipiv, T* work, int lwork);
template<>
int lapack_sytrf<double>(char uplo, int n, double* a, int lda, int *ipiv, double* work, int lwork) {
   int info;
   dsytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info);
   return info;
}

/* _SYRK */
template <typename T>
void syrk(char uplo, char trans, int n, int k, T alpha, const T* a, int lda, T beta, T* c, int ldc);
template <>
void syrk<double>(char uplo, char trans, int n, int k, double alpha, const double* a, int lda, double beta, double* c, int ldc) {
   dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

/* _TRSV */
template <typename T>
void trsv(char uplo, char trans, char diag, int n, const T* a, int lda, T* x, int incx);
template <>
void trsv<double>(char uplo, char trans, char diag, int n, const double* a, int lda, double* x, int incx) {
   dtrsv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

/* _TRSM */
template <typename T>
void trsm(char side, char uplo, char transa, char diag, int m, int n, T alpha, const T* a, int lda, T* b, int ldb);
template <>
void trsm<double>(char side, char uplo, char transa, char diag, int m, int n, double alpha, const double* a, int lda, double* b, int ldb) {
   dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}
