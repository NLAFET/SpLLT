#include "common.hpp"

#ifdef USE_COMPLEX
#ifdef USE_FLOAT
#define dtype MKL_Complex8
#define Xpotrf CMPLX_LAPACK(potrf)
#else // USE_DOUBLE
#define dtype MKL_Complex16
#define Xpotrf CMPLX_LAPACK(potrf)
#endif // USE_FLOAT
#else // USE_REAL
#ifdef USE_FLOAT
#define dtype float
#define Xpotrf REAL_LAPACK(potrf)
#else // USE_DOUBLE
#define dtype double
#define Xpotrf REAL_LAPACK(potrf)
#endif // USE_FLOAT
#endif // USE_COMPLEX

static const int CACHE_LINE_BYTES = 64;
static const int CACHE_LINE_ELEMS = int(CACHE_LINE_BYTES / sizeof(dtype));

static int Nmin = 0, Nmax = 0, Nstep = 0, _samples = 0, lda = 0;
static dtype *A = (dtype*)NULL, *B = (dtype*)NULL;

static void alloc_cpu_mtx()
{
  const int rem = Nmax % CACHE_LINE_ELEMS;
  lda = (rem ? (Nmax + (CACHE_LINE_ELEMS - rem)) : Nmax);
  const size_t siz1 = size_t(lda) * Nmax;
  const size_t size = 2 * siz1 * sizeof(dtype);
  if (size > 0) {
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
    (void)memset(A, 0, size);
    B = A + siz1;
  }
}

// TODO: fix for COMPLEX
static double init_cpu_mtx()
{
  static const int idist = 1;
  int iseed[4] = { 0, 1, 2, 3 };
  const int k = Nmax - 1;
  int info = 0;

  const double go = omp_get_wtime();

  dtype *const wrk = (dtype*)calloc(3 * Nmax, sizeof(dtype)); const int lin1 = __LINE__;
  if (!wrk) {
    (void)fprintf(stderr, "[%s@%s:%d,%d] ", __FUNCTION__, __FILE__, lin1, errno);
    perror("calloc");
    exit(errno);
  }

  // diagonal
  REAL_LAPACK(larnv)(&idist, iseed, &Nmax, wrk);
  // A
  REAL_LAPACK(lagsy)(&Nmax, &k, wrk, A, &lda, iseed, wrk + Nmax, &info); const int lin2 = __LINE__;
  if (info) {
    (void)fprintf(stderr, "[%s@%s:%d] INFO = %d\n", __FUNCTION__, __FILE__, lin2, info);
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
     , &n, &n, A, &lda, B, &lda);
  return (omp_get_wtime() - go);
}

static double potrf_cpu(const bool upper, const int n)
{
  const double go = omp_get_wtime();
  if (n >= Nmin) {
    if (n > Nmax) {
      (void)fprintf(stderr, "[%s@%s] n == %d > Nmax == %d\n", __FUNCTION__, __FILE__, n, Nmax);
      exit(n);
    }
    int info = 0;
    Xpotrf((upper ? "U" : "L"), &n, A, &lda, &info); const int lin = __LINE__;
    if (info) {
      (void)fprintf(stderr, "[%s@%s:%d] INFO = %d\n", __FUNCTION__, __FILE__, lin, info);
      exit(info);
    }
  }
  else {
    (void)fprintf(stderr, "[%s@%s] n == %d < Nmin == %d\n", __FUNCTION__, __FILE__, n, Nmin);
    exit(n);
  }
  return (omp_get_wtime() - go);
}

static void free_cpu_mtx()
{
  if (A) {
    free(A);
    A = (dtype*)NULL;
  }
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

#ifndef NDEBUG
  (void)fprintf(stdout, "[omp_get_wtick] %#.17E s\n", omp_get_wtick());
#endif // !NDEBUG

  alloc_cpu_mtx();
  const double init_time = init_cpu_mtx();
#ifndef NDEBUG
  (void)fprintf(stdout, "[init_cpu_mtx] %#.17E s\n", init_time);
#endif // !NDEBUG

  (void)fprintf(stdout, "\"N\",\"COPY_H2H_MIN_s\",\"COPY_H2H_AVG_s\",\"COPY_H2H_MAX_s\",\"LPOTRF_MIN_s\",\"LPOTRF_AVG_s\",\"LPOTRF_MAX_s\",\"UPOTRF_MIN_s\",\"UPOTRF_AVG_s\",\"UPOTRF_MAX_s\"\n");
  for (int n = Nmin; n <= Nmax; n += Nstep) {
    double Lcopy_times_min = INFINITY;
    double Lcopy_times_max = -0.0;
    double Lcopy_times_avg = -0.0;

    double Ucopy_times_min = INFINITY;
    double Ucopy_times_max = -0.0;
    double Ucopy_times_avg = -0.0;

    double Lpotrf_times_min = INFINITY;
    double Lpotrf_times_max = -0.0;
    double Lpotrf_times_avg = -0.0;

    double Upotrf_times_min = INFINITY;
    double Upotrf_times_max = -0.0;
    double Upotrf_times_avg = -0.0;

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

      const double Lpotrf_time = potrf_cpu(false, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[potrf_cpu(%d),%d,L] %#.17E s\n", n, sample, Lpotrf_time);
#endif // !NDEBUG
      if (Lpotrf_time < Lpotrf_times_min)
        Lpotrf_times_min = Lpotrf_time;
      if (Lpotrf_time > Lpotrf_times_max)
        Lpotrf_times_max = Lpotrf_time;
      Lpotrf_times_avg += Lpotrf_time / _samples;

      const double Ucopy_time = copy_mtx_cpu2cpu(true, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[copy_mtx_cpu2cpu(%d),%d,U] %#.17E s\n", n, sample, Ucopy_time);
#endif // !NDEBUG
      if (Ucopy_time < Ucopy_times_min)
        Ucopy_times_min = Ucopy_time;
      if (Ucopy_time > Ucopy_times_max)
        Ucopy_times_max = Ucopy_time;
      Ucopy_times_avg += Ucopy_time / _samples;

      const double Upotrf_time = potrf_cpu(true, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[potrf_cpu(%d),%d,U] %#.17E s\n", n, sample, Upotrf_time);
#endif // !NDEBUG
      if (Upotrf_time < Upotrf_times_min)
        Upotrf_times_min = Upotrf_time;
      if (Upotrf_time > Upotrf_times_max)
        Upotrf_times_max = Upotrf_time;
      Upotrf_times_avg += Upotrf_time / _samples;
    }

    const double copy_times_min = ((Lcopy_times_min <= Ucopy_times_min) ? Lcopy_times_min : Ucopy_times_min);
    const double copy_times_max = ((Lcopy_times_max >= Ucopy_times_max) ? Lcopy_times_max : Ucopy_times_max);
    const double copy_times_avg = (Lcopy_times_avg + Ucopy_times_avg) / 2;

    (void)fprintf(stdout, "%d,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E\n", n,
                  copy_times_min, copy_times_avg, copy_times_max,
                  Lpotrf_times_min, Lpotrf_times_avg, Lpotrf_times_max,
                  Upotrf_times_min, Upotrf_times_avg, Upotrf_times_max);
  }
  free_cpu_mtx();

  return EXIT_SUCCESS;
}
