/* xTRSM
 * SIDE   = L/R
 * UPLO   = L
 * TRANSA = N
 * DIAG   = N
 * ALPHA  = 1.0
 * M   == N
 * LDA == LDB
 */

#ifndef USE_MKL
#define USE_MKL
#endif // USE_MKL
#include "common.hpp"

#include "cublas_v2.h"

#ifdef USE_COMPLEX
#ifdef USE_FLOAT
#define dtype cuComplex
#define btype float
#define Xtrsm cublasCtrsm
#else // USE_DOUBLE
#define dtype cuDoubleComplex
#define btype double
#define Xtrsm cublasZtrsm
#endif // USE_FLOAT
#else // USE_REAL
#ifdef USE_FLOAT
#define dtype float
#define Xtrsm cublasStrsm
#else // USE_DOUBLE
#define dtype double
#define Xtrsm cublasDtrsm
#endif // USE_FLOAT
#define btype dtype
#endif // USE_COMPLEX

static const char *const lin_fmt = "%d,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E\n";

__constant__ dtype one_gpu = 1.0;
static const dtype *alpha = (const dtype*)NULL;

static int Nmin = 0, Nmax = 0, Nstep = 0, _samples = 0, device_ = 0, _devices = 0, lda = 0;
static dtype *Agpu = (dtype*)NULL, *Xgpu = (dtype*)NULL, *Acpu = (dtype*)NULL, *Bcpu = (dtype*)NULL;

static cublasHandle_t handle;

static double device_count()
{
  const double go = omp_get_wtime();

  const cudaError_t error = cudaGetDeviceCount(&_devices); const int lin = __LINE__;
  switch (error) {
  case cudaSuccess:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
    break;
  case cudaErrorNoDevice:
    (void)fprintf(stderr, "[%s@%s:%d] NoDevice\n", __FUNCTION__, __FILE__, lin);
    exit(error);
  case cudaErrorInsufficientDriver:
    (void)fprintf(stderr, "[%s@%s:%d] InsufficientDriver\n", __FUNCTION__, __FILE__, lin);
    exit(error);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin, error);
    exit(error);
  }

  return (omp_get_wtime() - go);
}

static double set_device()
{
  const double go = omp_get_wtime();

  int device = 0;
  (void)cudaGetDevice(&device);
  if (device != device_) {
    const cudaError_t error = cudaSetDevice(device_); const int lin = __LINE__;
    switch (error) {
    case cudaSuccess:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
      break;
    case cudaErrorInvalidDevice:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidDevice\n", __FUNCTION__, __FILE__, lin);
      exit(error);
    case cudaErrorDeviceAlreadyInUse:
      (void)fprintf(stderr, "[%s@%s:%d] DeviceAlreadyInUse\n", __FUNCTION__, __FILE__, lin);
      exit(error);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin, error);
      exit(error);
    }
  }

  return (omp_get_wtime() - go);
}

static double create_handle()
{
  const double go = omp_get_wtime();

  cublasStatus_t status = cublasCreate(&handle); const int lin1 = __LINE__;
  switch (status) {
  case CUBLAS_STATUS_SUCCESS:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin1);
#endif // !NDEBUG
    break;
  case CUBLAS_STATUS_NOT_INITIALIZED:
    (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin1);
    exit(status);
  case CUBLAS_STATUS_ALLOC_FAILED:
    (void)fprintf(stderr, "[%s@%s:%d] ALLOC_FAILED\n", __FUNCTION__, __FILE__, lin1);
    exit(status);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin1, status);
    exit(status);
  }

  status = cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE); const int lin2 = __LINE__;
  switch (status) {
  case CUBLAS_STATUS_SUCCESS:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin2);
#endif // !NDEBUG
    break;
  case CUBLAS_STATUS_NOT_INITIALIZED:
    (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin2);
    exit(status);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin2, status);
    exit(status);
  }

  status = cublasSetAtomicsMode(handle, CUBLAS_ATOMICS_ALLOWED); const int lin3 = __LINE__;
  switch (status) {
  case CUBLAS_STATUS_SUCCESS:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin3);
#endif // !NDEBUG
    break;
  case CUBLAS_STATUS_NOT_INITIALIZED:
    (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin3);
    exit(status);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin3, status);
    exit(status);
  }

  const cudaError_t err = cudaGetSymbolAddress((void**)&alpha, one_gpu); const int lin4 = __LINE__;
  switch (err) {
  case cudaSuccess:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin4);
#endif // !NDEBUG
    break;
  case cudaErrorInvalidSymbol:
    (void)fprintf(stderr, "[%s@%s:%d] InvalidSymbol\n", __FUNCTION__, __FILE__, lin4);
    exit(err);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin4, err);
    exit(err);
  }

  return (omp_get_wtime() - go);
}

static double alloc_gpu_mtx()
{
  const double go = omp_get_wtime();

  size_t pitch = 0;
  const cudaError_t err1 = cudaMallocPitch(&Agpu, &pitch, Nmax * sizeof(dtype), 2 * Nmax); const int lin1 = __LINE__;
  switch (err1) {
  case cudaSuccess:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin1);
#endif // !NDEBUG
    break;
  case cudaErrorMemoryAllocation:
    (void)fprintf(stderr, "[%s@%s:%d] MemoryAllocation\n", __FUNCTION__, __FILE__, lin1);
    exit(err1);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin1, err1);
    exit(err1);
  }
  lda = int(pitch / sizeof(dtype));
#ifndef NDEBUG
  (void)fprintf(stdout, "lda = %d\n", lda);
#endif // !NDEBUG
  Xgpu = Agpu + lda * Nmax;
  const double end = (omp_get_wtime() - go);
  // don't time clearing the memory
  const cudaError_t err2 = cudaMemset2D(Agpu, pitch, 0, pitch, 2 * Nmax); const int lin2 = __LINE__;
  switch (err2) {
  case cudaSuccess:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin2);
#endif // !NDEBUG
    break;
  case cudaErrorInvalidValue:
    (void)fprintf(stderr, "[%s@%s:%d] InvalidValue\n", __FUNCTION__, __FILE__, lin2);
    exit(err2);
  case cudaErrorInvalidDevicePointer:
    (void)fprintf(stderr, "[%s@%s:%d] InvalidDevicePointer\n", __FUNCTION__, __FILE__, lin2);
    exit(err2);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin2, err2);
    exit(err2);
  }

  return end;
}

static double alloc_cpu_mtx()
{
  const double go = omp_get_wtime();

  const size_t size = size_t(lda) * 2 * Nmax * sizeof(dtype);
  const cudaError_t error = cudaMallocHost(&Acpu, size); const int lin = __LINE__;
  switch (error) {
  case cudaSuccess:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
    break;
  case cudaErrorMemoryAllocation:
    (void)fprintf(stderr, "[%s@%s:%d] MemoryAllocation\n", __FUNCTION__, __FILE__, lin);
    exit(error);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin, error);
    exit(error);
  }
  Bcpu = Acpu + lda * Nmax;
  const double end = (omp_get_wtime() - go);
  // don't time clearing the memory
  (void)memset(Acpu, 0, size);

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
  // Acpu
  CMPLX_LAPACK(laghe)(&Nmax, &k, (btype*)wrk, (MKL_Complex*)Acpu, &lda, iseed, (MKL_Complex*)(wrk + Nmax), &info); const int lin2 = __LINE__;
#else // USE_REAL
  // Diagonal
  REAL_LAPACK(larnv)(&idist, iseed, &Nmax, wrk);
  // Acpu
  REAL_LAPACK(lagsy)(&Nmax, &k, wrk, Acpu, &lda, iseed, wrk + Nmax, &info); const int lin2 = __LINE__;
#endif // USE_COMPLEX
  if (info) {
    (void)fprintf(stderr, "[%s@%s:%d] INFO = %d\n", __FUNCTION__, __FILE__, lin2, info);
    exit(info);
  }

#ifdef USE_COMPLEX
  // Diagonal
  REAL_LAPACK(larnv)(&idist, iseed, &Nmax, (btype*)wrk);
  // Acpu
  CMPLX_LAPACK(laghe)(&Nmax, &k, (btype*)wrk, (MKL_Complex*)Acpu, &lda, iseed, (MKL_Complex*)(wrk + Nmax), &info); const int lin3 = __LINE__;
#else // USE_REAL
  // Diagonal
  REAL_LAPACK(larnv)(&idist, iseed, &Nmax, wrk);
  // Acpu
  REAL_LAPACK(lagsy)(&Nmax, &k, wrk, Bcpu, &lda, iseed, wrk + Nmax, &info); const int lin3 = __LINE__;
#endif // USE_COMPLEX
  if (info) {
    (void)fprintf(stderr, "[%s@%s:%d] INFO = %d\n", __FUNCTION__, __FILE__, lin3, info);
    exit(info);
  }

  free(wrk);
  return (omp_get_wtime() - go);
}

static double copy_mtx_cpu2gpu(const int n)
{
  const double go = omp_get_wtime();

  if (n >= Nmin) {
    if (n > Nmax) {
      (void)fprintf(stderr, "[%s@%s] n == %d > Nmax == %d\n", __FUNCTION__, __FILE__, n, Nmax);
      exit(n);
    }
    const size_t pitch = lda * sizeof(dtype);
    const cudaError_t err1 = cudaMemcpy2D(Agpu, pitch, Acpu, pitch, n * sizeof(dtype), n, cudaMemcpyHostToDevice); const int lin1 = __LINE__;
    switch (err1) {
    case cudaSuccess:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin1);
#endif // !NDEBUG
      break;
    case cudaErrorInvalidValue:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidValue\n", __FUNCTION__, __FILE__, lin1);
      exit(err1);
    case cudaErrorInvalidPitchValue:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidPitchValue\n", __FUNCTION__, __FILE__, lin1);
      exit(err1);
    case cudaErrorInvalidDevicePointer:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidDevicePointer\n", __FUNCTION__, __FILE__, lin1);
      exit(err1);
    case cudaErrorInvalidMemcpyDirection:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidMemcpyDirection\n", __FUNCTION__, __FILE__, lin1);
      exit(err1);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin1, err1);
      exit(err1);
    }
    const cudaError_t err2 = cudaMemcpy2D(Xgpu, pitch, Bcpu, pitch, n * sizeof(dtype), n, cudaMemcpyHostToDevice); const int lin2 = __LINE__;
    switch (err2) {
    case cudaSuccess:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin2);
#endif // !NDEBUG
      break;
    case cudaErrorInvalidValue:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidValue\n", __FUNCTION__, __FILE__, lin2);
      exit(err2);
    case cudaErrorInvalidPitchValue:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidPitchValue\n", __FUNCTION__, __FILE__, lin2);
      exit(err2);
    case cudaErrorInvalidDevicePointer:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidDevicePointer\n", __FUNCTION__, __FILE__, lin2);
      exit(err2);
    case cudaErrorInvalidMemcpyDirection:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidMemcpyDirection\n", __FUNCTION__, __FILE__, lin2);
      exit(err2);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin2, err2);
      exit(err2);
    }
    // just to be sure...
    (void)cudaDeviceSynchronize();
  }
  else {
    (void)fprintf(stderr, "[%s@%s] n == %d < Nmin == %d\n", __FUNCTION__, __FILE__, n, Nmin);
    exit(n);
  }

  return (omp_get_wtime() - go);
}

static double trsm_gpu(const bool right, const int n)
{
  const double go = omp_get_wtime();

  if (n >= Nmin) {
    if (n > Nmax) {
      (void)fprintf(stderr, "[%s@%s] n == %d > Nmax == %d\n", __FUNCTION__, __FILE__, n, Nmax);
      exit(n);
    }
    const cublasSideMode_t side = (right ? CUBLAS_SIDE_RIGHT : CUBLAS_SIDE_LEFT);
    const cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    const cublasOperation_t trans = CUBLAS_OP_N;
    const cublasDiagType_t diag = CUBLAS_DIAG_NON_UNIT;
    
    const cublasStatus_t status = Xtrsm(handle, side, uplo, trans, diag, n, n, alpha, Agpu, lda, Xgpu, lda); const int lin = __LINE__;
    (void)cudaDeviceSynchronize();
    switch (status) {
    case CUBLAS_STATUS_SUCCESS:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
      break;
    case CUBLAS_STATUS_NOT_INITIALIZED:
      (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin);
      exit(status);
    case CUBLAS_STATUS_INVALID_VALUE:
      (void)fprintf(stderr, "[%s@%s:%d] INVALID_VALUE\n", __FUNCTION__, __FILE__, lin);
      exit(status);
    case CUBLAS_STATUS_ARCH_MISMATCH:
      (void)fprintf(stderr, "[%s@%s:%d] ARCH_MISMATCH\n", __FUNCTION__, __FILE__, lin);
      exit(status);
    case CUBLAS_STATUS_INTERNAL_ERROR:
      (void)fprintf(stderr, "[%s@%s:%d] INTERNAL_ERROR\n", __FUNCTION__, __FILE__, lin);
      exit(status);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin, status);
      exit(status);
    }
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

  if (Acpu) {
    const cudaError_t error = cudaFreeHost(Acpu); const int lin = __LINE__;
    switch (error) {
    case cudaSuccess:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
      break;
    case cudaErrorInitializationError:
      (void)fprintf(stderr, "[%s@%s:%d] InitializationError\n", __FUNCTION__, __FILE__, lin);
      exit(error);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin, error);
      exit(error);
    }
    Acpu = (dtype*)NULL;
  }

  return (omp_get_wtime() - go);
}

static double free_gpu_mtx()
{
  const double go = omp_get_wtime();

  const cudaError_t error = cudaFree(Agpu); const int lin = __LINE__;
  switch (error) {
  case cudaSuccess:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
    break;
  case cudaErrorInvalidDevicePointer:
    (void)fprintf(stderr, "[%s@%s:%d] InvalidDevicePointer\n", __FUNCTION__, __FILE__, lin);
    exit(error);
  case cudaErrorInitializationError:
    (void)fprintf(stderr, "[%s@%s:%d] InitializationError\n", __FUNCTION__, __FILE__, lin);
    exit(error);    
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin, error);
    exit(error);
  }
  Agpu = (dtype*)NULL;

  return (omp_get_wtime() - go);
}

static double destroy_handle()
{
  const double go = omp_get_wtime();

  const cublasStatus_t status = cublasDestroy(handle); const int lin = __LINE__;
  switch (status) {
  case CUBLAS_STATUS_SUCCESS:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
    break;
  case CUBLAS_STATUS_NOT_INITIALIZED:
    (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin);
    exit(status);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin, status);
    exit(status);
  }

  return (omp_get_wtime() - go);
}

int main(int argc, char* argv[])
{
  if ((argc < 5) || (argc > 6)) {
    (void)fprintf(stderr, "%s Nmin Nmax Nstep #samples [device#]\n", *argv);
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

  if ((argc == 6) && ((device_ = atoi(argv[5])) < 0)) {
    (void)fprintf(stderr, "device# == %d < 0\n", device_);
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

  (void)device_count();
  if (device_ > _devices) {
    (void)fprintf(stderr, "device# == %d > #devices == %d\n", device_, _devices);
    return EXIT_FAILURE;
  }
  (void)set_device();

  (void)create_handle();
  const double agpu_time = alloc_gpu_mtx();
  const double acpu_time = alloc_cpu_mtx();
  const double init_time = init_cpu_mtx();
#ifndef NDEBUG
  (void)fprintf(stdout, "[init_cpu_mtx] %#.17E s\n", init_time);
#endif // !NDEBUG

  (void)fprintf(stdout, "\"N\",\"COPY_H2D_MIN_s\",\"COPY_H2D_AVG_s\",\"COPY_H2D_MAX_s\",\"LTRSM_MIN_s\",\"LTRSM_AVG_s\",\"LTRSM_MAX_s\",\"RTRSM_MIN_s\",\"RTRSM_AVG_s\",\"RTRSM_MAX_s\"\n");
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
      const double Lcopy_time = copy_mtx_cpu2gpu(n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[copy_mtx_cpu2gpu(%d),%d,L] %#.17E s\n", n, sample, Lcopy_time);
#endif // !NDEBUG
      if (Lcopy_time < Lcopy_times_min)
        Lcopy_times_min = Lcopy_time;
      if (Lcopy_time > Lcopy_times_max)
        Lcopy_times_max = Lcopy_time;
      Lcopy_times_avg += Lcopy_time / _samples;

      const double Ltrsm_time = trsm_gpu(false, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[trsm_gpu(%d),%d,L] %#.17E s\n", n, sample, Ltrsm_time);
#endif // !NDEBUG
      if (Ltrsm_time < Ltrsm_times_min)
        Ltrsm_times_min = Ltrsm_time;
      if (Ltrsm_time > Ltrsm_times_max)
        Ltrsm_times_max = Ltrsm_time;
      Ltrsm_times_avg += Ltrsm_time / _samples;

      const double Rcopy_time = copy_mtx_cpu2gpu(n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[copy_mtx_cpu2gpu(%d),%d,R] %#.17E s\n", n, sample, Rcopy_time);
#endif // !NDEBUG
      if (Rcopy_time < Rcopy_times_min)
        Rcopy_times_min = Rcopy_time;
      if (Rcopy_time > Rcopy_times_max)
        Rcopy_times_max = Rcopy_time;
      Rcopy_times_avg += Rcopy_time / _samples;

      const double Rtrsm_time = trsm_gpu(true, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[trsm_gpu(%d),%d,R] %#.17E s\n", n, sample, Rtrsm_time);
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
  const double fgpu_time = free_gpu_mtx();
  (void)destroy_handle();
  const double awrk_time = -0.0;
  const double fwrk_time = -0.0;

  (void)fprintf(stdout, lin_fmt, -_samples,
                resol, double(mkl_nthr), init_time,
                acpu_time, agpu_time, awrk_time,
                fcpu_time, fgpu_time, fwrk_time);
  (void)fflush(stdout);

  return EXIT_SUCCESS;
}
