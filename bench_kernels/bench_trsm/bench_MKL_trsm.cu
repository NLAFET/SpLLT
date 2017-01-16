/* xTRSM
 * SIDE   = L/R
 * UPLO   = L
 * TRANSA = N
 * DIAG   = N
 * ALPHA  = 1.0
 * M   == N
 * LDA == LDB
 */

#include "common.hpp"

// #include "mkl.h"

#ifdef USE_COMPLEX
#ifdef USE_FLOAT
#define dtype MKL_Complex8
#define btype float
#define Xtrsm CMPLX_BLAS(trsm)
#else // USE_DOUBLE
#define dtype MKL_Complex16
#define btype double
#define Xtrsm CMPLX_BLAS(trsm)
#endif // USE_FLOAT
#else // USE_REAL
#ifdef USE_FLOAT
#define dtype float
#define Xtrsm REAL_BLAS(trsm)
#else // USE_DOUBLE
#define dtype double
#define Xtrsm REAL_BLAS(trsm)
#endif // USE_FLOAT
#endif // USE_COMPLEX

static const char *const lin_fmt = "%d,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E\n";

static const int CACHE_LINE_BYTES = 64;
static const int CACHE_LINE_ELEMS = int(CACHE_LINE_BYTES / sizeof(dtype));

static int Nmin = 0, Nmax = 0, Nstep = 0, _samples = 0, lda = 0;
static dtype *A = (dtype*)NULL, *B = (dtype*)NULL, *X = (dtype*)NULL;

//

static double alloc_cpu_mtx()
{
  const double go = omp_get_wtime();
  const int rem = Nmax % CACHE_LINE_ELEMS;
  lda = (rem ? (Nmax + (CACHE_LINE_ELEMS - rem)) : Nmax);
  const size_t siz1 = size_t(lda) * Nmax;
  const size_t size = 3 * siz1 * sizeof(dtype);
  errno = posix_memalign((void**)&A, CACHE_LINE_BYTES, size); const int lin = __LINE__;
  switch (errno) {
  case 0:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
    break;
  default:
    (void)fprintf(stderr, "[%s@%s:%d,%d] ", __FUNCTION__, __FILE__, lin, errno);
    perror("posix_memalign");
    exit(errno);
  }
  B = A + siz1;
  X = B + siz1;
  const double end = (omp_get_wtime() - go);
  // don't time clearing the memory
  (void)memset(A, 0, size);
  return end;
}

static double init_cpu_mtx()
{
  const double go = omp_get_wtime();

  static const int idist = 1;
  int iseed[4] = { 0, 1, 2, 3 };
  const int k = Nmax - 1;
  int info = 0;

  dtype *const wrk = (dtype*)calloc(3 * Nmax, sizeof(dtype)); const int lin1 = __LINE__;
  if (!wrk) {
    (void)fprintf(stderr, "[%s@%s:%d,%d] ", __FUNCTION__, __FILE__, lin1, errno);
    perror("calloc");
    exit(errno);
  }

#ifdef USE_COMPLEX
  // Diagonal
  REAL_LAPACK(larnv)(&idist, iseed, &Nmax, (btype*)wrk);
  // A
  CMPLX_LAPACK(laghe)(&Nmax, &k, (btype*)wrk, A, &lda, iseed, wrk + Nmax, &info); const int lin2 = __LINE__;
#else // USE_REAL
  // Diagonal
  REAL_LAPACK(larnv)(&idist, iseed, &Nmax, wrk);
  // A
  REAL_LAPACK(lagsy)(&Nmax, &k, wrk, A, &lda, iseed, wrk + Nmax, &info); const int lin2 = __LINE__;
#endif // USE_COMPLEX
  if (info) {
    (void)fprintf(stderr, "[%s@%s:%d] INFO = %d\n", __FUNCTION__, __FILE__, lin2, info);
    exit(info);
  }

#ifdef USE_COMPLEX
  // Diagonal
  REAL_LAPACK(larnv)(&idist, iseed, &Nmax, (btype*)wrk);
  // A
  CMPLX_LAPACK(laghe)(&Nmax, &k, (btype*)wrk, B, &lda, iseed, wrk + Nmax, &info); const int lin3 = __LINE__;
#else // USE_REAL
  // Diagonal
  REAL_LAPACK(larnv)(&idist, iseed, &Nmax, wrk);
  // A
  REAL_LAPACK(lagsy)(&Nmax, &k, wrk, B, &lda, iseed, wrk + Nmax, &info); const int lin3 = __LINE__;
#endif // USE_COMPLEX
  if (info) {
    (void)fprintf(stderr, "[%s@%s:%d] INFO = %d\n", __FUNCTION__, __FILE__, lin3, info);
    exit(info);
  }

  free(wrk);
  return (omp_get_wtime() - go);
}

static double copy_mtx_cpu2cpu(const bool upper, const int n)
{
  const double go = omp_get_wtime();
#ifdef USE_COMPLEX
  CMPLX_LAPACK(lacpy)
#else // USE_REAL
  REAL_LAPACK(lacpy)
#endif // USE_COMPLEX
    (
#ifdef LACPY_ALL
     "A"
#else // U/L
     (upper ? "U" : "L")
#endif // LACPY_ALL
     , &n, &n, B, &lda, X, &lda);
  return (omp_get_wtime() - go);
}

static double trsm_cpu(const bool right, const int n)
{
  const dtype one = dtype(1.0);
  const double go = omp_get_wtime();

  if (n >= Nmin) {
    if (n > Nmax) {
      (void)fprintf(stderr, "[%s@%s] n == %d > Nmax == %d\n", __FUNCTION__, __FILE__, n, Nmax);
      exit(n);
    }
    Xtrsm((right ? "R" : "L"), "L", "N", "N", &n, &n, &one, A, &lda, X, &lda);
  }
  else {
    (void)fprintf(stderr, "[%s@%s] n == %d < Nmin == %d\n", __FUNCTION__, __FILE__, n, Nmin);
    exit(n);
  }

  return (omp_get_wtime() - go);
}

static double free_cpu_mtx()
{
  const double go = omp_get_wtime();

  if (A) {
    free(A);
    A = (dtype*)NULL;
  }

  return (omp_get_wtime() - go);
}

int main(int argc, char* argv[])
{
  if (argc != 5) {
    (void)fprintf(stderr, "%s Nmin Nmax Nstep #samples\n", *argv);
    return EXIT_FAILURE;
  }

  if ((Nmin = atoi(argv[1])) <= 0) {
    (void)fprintf(stderr, "Nmin == %d <= 0\n", Nmin);
    return EXIT_FAILURE;
  }

  if ((Nmax = atoi(argv[2])) <= 0) {
    (void)fprintf(stderr, "Nmax == %d <= 0\n", Nmax);
    return EXIT_FAILURE;
  }

  if (Nmax < Nmin) {
    (void)fprintf(stderr, "Nmax == %d < Nmin == %d\n", Nmax, Nmin);
    return EXIT_FAILURE;
  }

  if ((Nstep = atoi(argv[3])) <= 0) {
    (void)fprintf(stderr, "Nstep == %d <= 0\n", Nstep);
    return EXIT_FAILURE;
  }

  if ((_samples = atoi(argv[4])) <= 0) {
    (void)fprintf(stderr, "#samples == %d <= 0\n", _samples);
    return EXIT_FAILURE;
  }

  const char *const env_nthr = getenv("MKL_NUM_THREADS");
  if (!env_nthr) {
    (void)fprintf(stderr, "MKL_NUM_THREADS environment variable not set\n");
    return EXIT_FAILURE;
  }
  const int mkl_nthr = atoi(env_nthr);
  if (mkl_nthr <= 0) {
    (void)fprintf(stderr, "MKL_NUM_THREADS = %d <= 0\n", mkl_nthr);
    return EXIT_FAILURE;
  }

  const double resol = omp_get_wtick();
#ifndef NDEBUG
  (void)fprintf(stdout, "[omp_get_wtick] %#.17E s\n", resol);
#endif // !NDEBUG

  const double acpu_time = alloc_cpu_mtx();
  const double init_time = init_cpu_mtx();
#ifndef NDEBUG
  (void)fprintf(stdout, "[init_cpu_mtx] %#.17E s\n", init_time);
#endif // !NDEBUG

  (void)fprintf(stdout, "\"N\",\"COPY_H2H_MIN_s\",\"COPY_H2H_AVG_s\",\"COPY_H2H_MAX_s\",\"LTRSM_MIN_s\",\"LTRSM_AVG_s\",\"LTRSM_MAX_s\",\"RTRSM_MIN_s\",\"RTRSM_AVG_s\",\"RTRSM_MAX_s\"\n");
  (void)fflush(stdout);

  for (int n = Nmin; n <= Nmax; n += Nstep) {
    double Lcopy_times_min = INFINITY;
    double Lcopy_times_max = -0.0;
    double Lcopy_times_avg = -0.0;

    double Rcopy_times_min = INFINITY;
    double Rcopy_times_max = -0.0;
    double Rcopy_times_avg = -0.0;

    double Ltrsm_times_min = INFINITY;
    double Ltrsm_times_max = -0.0;
    double Ltrsm_times_avg = -0.0;

    double Rtrsm_times_min = INFINITY;
    double Rtrsm_times_max = -0.0;
    double Rtrsm_times_avg = -0.0;

    for (int sample = 0; sample < _samples; ++sample) {
      const double Lcopy_time = copy_mtx_cpu2cpu(false, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[copy_mtx_cpu2cpu(%d),%d,L] %#.17E s\n", n, sample, Lcopy_time);
#endif // !NDEBUG
      if (Lcopy_time < Lcopy_times_min)
        Lcopy_times_min = Lcopy_time;
      if (Lcopy_time > Lcopy_times_max)
        Lcopy_times_max = Lcopy_time;
      Lcopy_times_avg += Lcopy_time / _samples;

      const double Ltrsm_time = trsm_cpu(false, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[trsm_cpu(%d),%d,L] %#.17E s\n", n, sample, Ltrsm_time);
#endif // !NDEBUG
      if (Ltrsm_time < Ltrsm_times_min)
        Ltrsm_times_min = Ltrsm_time;
      if (Ltrsm_time > Ltrsm_times_max)
        Ltrsm_times_max = Ltrsm_time;
      Ltrsm_times_avg += Ltrsm_time / _samples;

      const double Rcopy_time = copy_mtx_cpu2cpu(true, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[copy_mtx_cpu2cpu(%d),%d,R] %#.17E s\n", n, sample, Rcopy_time);
#endif // !NDEBUG
      if (Rcopy_time < Rcopy_times_min)
        Rcopy_times_min = Rcopy_time;
      if (Rcopy_time > Rcopy_times_max)
        Rcopy_times_max = Rcopy_time;
      Rcopy_times_avg += Rcopy_time / _samples;

      const double Rtrsm_time = trsm_cpu(true, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[trsm_cpu(%d),%d,R] %#.17E s\n", n, sample, Rtrsm_time);
#endif // !NDEBUG
      if (Rtrsm_time < Rtrsm_times_min)
        Rtrsm_times_min = Rtrsm_time;
      if (Rtrsm_time > Rtrsm_times_max)
        Rtrsm_times_max = Rtrsm_time;
      Rtrsm_times_avg += Rtrsm_time / _samples;
    }

    const double copy_times_min = ((Lcopy_times_min <= Rcopy_times_min) ? Lcopy_times_min : Rcopy_times_min);
    const double copy_times_max = ((Lcopy_times_max >= Rcopy_times_max) ? Lcopy_times_max : Rcopy_times_max);
    const double copy_times_avg = (Lcopy_times_avg + Rcopy_times_avg) / 2;

    (void)fprintf(stdout, lin_fmt, n,
                  copy_times_min, copy_times_avg, copy_times_max,
                  Ltrsm_times_min, Ltrsm_times_avg, Ltrsm_times_max,
                  Rtrsm_times_min, Rtrsm_times_avg, Rtrsm_times_max);
    (void)fflush(stdout);
  }
  const double fcpu_time = free_cpu_mtx();
  const double agpu_time = -0.0;
  const double awrk_time = -0.0;
  const double fwrk_time = -0.0;
  const double fgpu_time = -0.0;

  (void)fprintf(stdout, lin_fmt, -_samples,
                resol, double(mkl_nthr), init_time,
                acpu_time, agpu_time, awrk_time,
                fcpu_time, fgpu_time, fwrk_time);
  (void)fflush(stdout);

  return EXIT_SUCCESS;
}
