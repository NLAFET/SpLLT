#include "common.hpp"

#include "cusolverDn.h"

#ifdef USE_COMPLEX
#ifdef USE_FLOAT
#define dtype cuComplex
#define btype float
#define cusolverDnXpotrf_bufferSize cusolverDnCpotrf_bufferSize
#define cusolverDnXpotrf cusolverDnCpotrf
#else // USE_DOUBLE
#define dtype cuDoubleComplex
#define btype double
#define cusolverDnXpotrf_bufferSize cusolverDnZpotrf_bufferSize
#define cusolverDnXpotrf cusolverDnZpotrf
#endif // USE_FLOAT
#else // USE_REAL
#ifdef USE_FLOAT
#define dtype float
#define cusolverDnXpotrf_bufferSize cusolverDnSpotrf_bufferSize
#define cusolverDnXpotrf cusolverDnSpotrf
#else // USE_DOUBLE
#define dtype double
#define cusolverDnXpotrf_bufferSize cusolverDnDpotrf_bufferSize
#define cusolverDnXpotrf cusolverDnDpotrf
#endif // USE_FLOAT
#endif // USE_COMPLEX

static const char *const lin_fmt = "%d,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E,%#.17E\n";

static int Nmin = 0, Nmax = 0, Nstep = 0, _samples = 0, device_ = 0, _devices = 0, lda = 0, Lwork = 0;
static dtype *Agpu = (dtype*)NULL, *Workspace = (dtype*)NULL, *Acpu = (dtype*)NULL;

static cusolverDnHandle_t handle;

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

  const cusolverStatus_t status = cusolverDnCreate(&handle); const int lin = __LINE__;
  switch (status) {
  case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
    break;
  case CUSOLVER_STATUS_NOT_INITIALIZED:
    (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin);
    exit(status);
  case CUSOLVER_STATUS_ALLOC_FAILED:
    (void)fprintf(stderr, "[%s@%s:%d] ALLOC_FAILED\n", __FUNCTION__, __FILE__, lin);
    exit(status);
  case CUSOLVER_STATUS_ARCH_MISMATCH:
    (void)fprintf(stderr, "[%s@%s:%d] ARCH_MISMATCH\n", __FUNCTION__, __FILE__, lin);
    exit(status);
  default:
    (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin, status);
    exit(status);
  }

  return (omp_get_wtime() - go);
}

static double alloc_gpu_mtx()
{
  const double go = omp_get_wtime();

  size_t pitch = 0;
  const cudaError_t err1 = cudaMallocPitch(&Agpu, &pitch, Nmax * sizeof(dtype), Nmax); const int lin1 = __LINE__;
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
  const double end = (omp_get_wtime() - go);
  // don't time clearing the memory
  const cudaError_t err2 = cudaMemset2D(Agpu, pitch, 0, pitch, Nmax); const int lin2 = __LINE__;
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

static double find_lwork()
{
  const double go = omp_get_wtime();

  int LworkL = -1, maxLworkL = 0;
  int LworkU = -1, maxLworkU = 0;

  for (int n = Nmin; n <= Nmax; n += Nstep) {
    const cusolverStatus_t stat1 = cusolverDnXpotrf_bufferSize(handle, CUBLAS_FILL_MODE_LOWER, n, Agpu, lda, &LworkL); const int lin1 = __LINE__;
    switch (stat1) {
    case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin1);
#endif // !NDEBUG
      break;
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin1);
      exit(stat1);
    case CUSOLVER_STATUS_INVALID_VALUE:
      (void)fprintf(stderr, "[%s@%s:%d] INVALID_VALUE\n", __FUNCTION__, __FILE__, lin1);
      exit(stat1);
    case CUSOLVER_STATUS_ARCH_MISMATCH:
      (void)fprintf(stderr, "[%s@%s:%d] ARCH_MISMATCH\n", __FUNCTION__, __FILE__, lin1);
      exit(stat1);
    case CUSOLVER_STATUS_INTERNAL_ERROR:
      (void)fprintf(stderr, "[%s@%s:%d] INTERNAL_ERROR\n", __FUNCTION__, __FILE__, lin1);
      exit(stat1);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin1, stat1);
      exit(stat1);
    }
    if (LworkL > maxLworkL)
      maxLworkL = LworkL;
    const cusolverStatus_t stat2 = cusolverDnXpotrf_bufferSize(handle, CUBLAS_FILL_MODE_UPPER, n, Agpu, lda, &LworkU); const int lin2 = __LINE__;
    switch (stat2) {
    case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin2);
#endif // !NDEBUG
      break;
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin2);
      exit(stat2);
    case CUSOLVER_STATUS_INVALID_VALUE:
      (void)fprintf(stderr, "[%s@%s:%d] INVALID_VALUE\n", __FUNCTION__, __FILE__, lin2);
      exit(stat2);
    case CUSOLVER_STATUS_ARCH_MISMATCH:
      (void)fprintf(stderr, "[%s@%s:%d] ARCH_MISMATCH\n", __FUNCTION__, __FILE__, lin2);
      exit(stat2);
    case CUSOLVER_STATUS_INTERNAL_ERROR:
      (void)fprintf(stderr, "[%s@%s:%d] INTERNAL_ERROR\n", __FUNCTION__, __FILE__, lin2);
      exit(stat2);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin2, stat2);
      exit(stat2);
    }
    if (LworkU > maxLworkU)
      maxLworkU = LworkU;
  }

  Lwork = ((maxLworkL >= maxLworkU) ? maxLworkL : maxLworkU);
#ifndef NDEBUG
  (void)fprintf(stdout, "Lwork = %d, LworkL = %d, LworkU = %d\n", Lwork, maxLworkL, maxLworkU);
#endif // !NDEBUG

  return (omp_get_wtime() - go);
}

static double alloc_gpu_wrk()
{
  const double go = omp_get_wtime();

  if (Lwork > 0) {
    const cudaError_t err1 = cudaMalloc(&Workspace, (Lwork + 1) * sizeof(dtype)); const int lin1 = __LINE__;
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
    const cudaError_t err2 = cudaMemset(Workspace, 0, (Lwork + 1) * sizeof(dtype)); const int lin2 = __LINE__;
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
  }

  return (omp_get_wtime() - go);
}

static double alloc_cpu_mtx()
{
  const double go = omp_get_wtime();

  const size_t size = size_t(lda) * Nmax * sizeof(dtype);
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
    const cudaError_t error = cudaMemcpy2D(Agpu, pitch, Acpu, pitch, n * sizeof(dtype), n, cudaMemcpyHostToDevice); const int lin1 = __LINE__;
    switch (error) {
    case cudaSuccess:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin1);
#endif // !NDEBUG
      break;
    case cudaErrorInvalidValue:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidValue\n", __FUNCTION__, __FILE__, lin1);
      exit(error);
    case cudaErrorInvalidPitchValue:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidPitchValue\n", __FUNCTION__, __FILE__, lin1);
      exit(error);
    case cudaErrorInvalidDevicePointer:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidDevicePointer\n", __FUNCTION__, __FILE__, lin1);
      exit(error);
    case cudaErrorInvalidMemcpyDirection:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidMemcpyDirection\n", __FUNCTION__, __FILE__, lin1);
      exit(error);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin1, error);
      exit(error);
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

static double potrf_gpu(const bool upper, const int n)
{
  const double go = omp_get_wtime();

  if (n >= Nmin) {
    if (n > Nmax) {
      (void)fprintf(stderr, "[%s@%s] n == %d > Nmax == %d\n", __FUNCTION__, __FILE__, n, Nmax);
      exit(n);
    }
    const cublasFillMode_t uplo = (upper ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
    const cusolverStatus_t status = cusolverDnXpotrf(handle, uplo, n, Agpu, lda, Workspace, Lwork, (int*)(Workspace + Lwork)); const int lin1 = __LINE__;
    int devInfo = 0;
    (void)cudaDeviceSynchronize();
    const cudaError_t error = cudaMemcpy(&devInfo, Workspace + Lwork, sizeof(int), cudaMemcpyDeviceToHost); const int lin2 = __LINE__;
    (void)cudaDeviceSynchronize();
    switch (error) {
    case cudaSuccess:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d] Success\n", __FUNCTION__, __FILE__, lin2);
#endif // !NDEBUG
      break;
    case cudaErrorInvalidValue:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidValue\n", __FUNCTION__, __FILE__, lin2);
      exit(error);
    case cudaErrorInvalidDevicePointer:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidDevicePointer\n", __FUNCTION__, __FILE__, lin2);
      exit(error);
    case cudaErrorInvalidMemcpyDirection:
      (void)fprintf(stderr, "[%s@%s:%d] InvalidMemcpyDirection\n", __FUNCTION__, __FILE__, lin2);
      exit(error);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin2, error);
      exit(error);      
    }
    switch (status) {
    case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
      (void)fprintf(stdout, "[%s@%s:%d,%d] SUCCESS\n", __FUNCTION__, __FILE__, lin1, devInfo);
#endif // !NDEBUG
      break;
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      (void)fprintf(stderr, "[%s@%s:%d,%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin1, devInfo);
      exit(status);
    case CUSOLVER_STATUS_INVALID_VALUE:
      (void)fprintf(stderr, "[%s@%s:%d,%d] INVALID_VALUE\n", __FUNCTION__, __FILE__, lin1, devInfo);
      exit(status);
    case CUSOLVER_STATUS_ARCH_MISMATCH:
      (void)fprintf(stderr, "[%s@%s:%d,%d] ARCH_MISMATCH\n", __FUNCTION__, __FILE__, lin1, devInfo);
      exit(status);
    case CUSOLVER_STATUS_INTERNAL_ERROR:
      (void)fprintf(stderr, "[%s@%s:%d,%d] INTERNAL_ERROR\n", __FUNCTION__, __FILE__, lin1, devInfo);
      exit(status);
    default:
      (void)fprintf(stderr, "[%s@%s:%d,%d] unknown error %d\n", __FUNCTION__, __FILE__, lin1, devInfo, status);
      exit(status);
    }
    if (devInfo) {
      (void)fprintf(stderr, "[%s@%s:%d] INFO = %d\n", __FUNCTION__, __FILE__, lin1, devInfo);
      exit(devInfo);
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

static double free_gpu_wrk()
{
  const double go = omp_get_wtime();

  if (Workspace) {
    const cudaError_t error = cudaFree(Workspace); const int lin = __LINE__;
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
    Workspace = (dtype*)NULL;
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

  const cusolverStatus_t status = cusolverDnDestroy(handle); const int lin = __LINE__;
  switch (status) {
  case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin);
#endif // !NDEBUG
    break;
  case CUSOLVER_STATUS_NOT_INITIALIZED:
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
  const double lwrk_time = find_lwork();
  double awrk_time = alloc_gpu_wrk();
  // number of tests
  int ntst = (Nmax - Nmin) + 1;
  // at least one test
  ntst = ((ntst < Nstep) ? 1 : (ntst / Nstep));
  // add the average time to find Lwork
  awrk_time += lwrk_time / (2 * ntst);
  const double acpu_time = alloc_cpu_mtx();
  const double init_time = init_cpu_mtx();
#ifndef NDEBUG
  (void)fprintf(stdout, "[init_cpu_mtx] %#.17E s\n", init_time);
#endif // !NDEBUG

  (void)fprintf(stdout, "\"N\",\"COPY_H2D_MIN_s\",\"COPY_H2D_AVG_s\",\"COPY_H2D_MAX_s\",\"LPOTRF_MIN_s\",\"LPOTRF_AVG_s\",\"LPOTRF_MAX_s\",\"UPOTRF_MIN_s\",\"UPOTRF_AVG_s\",\"UPOTRF_MAX_s\"\n");
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
      const double Lcopy_time = copy_mtx_cpu2gpu(n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[copy_mtx_cpu2gpu(%d),%d,L] %#.17E s\n", n, sample, Lcopy_time);
#endif // !NDEBUG
      if (Lcopy_time < Lcopy_times_min)
        Lcopy_times_min = Lcopy_time;
      if (Lcopy_time > Lcopy_times_max)
        Lcopy_times_max = Lcopy_time;
      Lcopy_times_avg += Lcopy_time / _samples;

      const double Lpotrf_time = potrf_gpu(false, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[potrf_gpu(%d),%d,L] %#.17E s\n", n, sample, Lpotrf_time);
#endif // !NDEBUG
      if (Lpotrf_time < Lpotrf_times_min)
        Lpotrf_times_min = Lpotrf_time;
      if (Lpotrf_time > Lpotrf_times_max)
        Lpotrf_times_max = Lpotrf_time;
      Lpotrf_times_avg += Lpotrf_time / _samples;

      const double Ucopy_time = copy_mtx_cpu2gpu(n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[copy_mtx_cpu2gpu(%d),%d,U] %#.17E s\n", n, sample, Ucopy_time);
#endif // !NDEBUG
      if (Ucopy_time < Ucopy_times_min)
        Ucopy_times_min = Ucopy_time;
      if (Ucopy_time > Ucopy_times_max)
        Ucopy_times_max = Ucopy_time;
      Ucopy_times_avg += Ucopy_time / _samples;

      const double Upotrf_time = potrf_gpu(true, n);
#ifndef NDEBUG
      (void)fprintf(stdout, "[potrf_gpu(%d),%d,U] %#.17E s\n", n, sample, Upotrf_time);
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

    (void)fprintf(stdout, lin_fmt, n,
                  copy_times_min, copy_times_avg, copy_times_max,
                  Lpotrf_times_min, Lpotrf_times_avg, Lpotrf_times_max,
                  Upotrf_times_min, Upotrf_times_avg, Upotrf_times_max);
  }
  const double fcpu_time = free_cpu_mtx();
  const double fgpu_time = free_gpu_wrk();
  const double fwrk_time = free_gpu_mtx();
  (void)destroy_handle();

  (void)fprintf(stdout, lin_fmt, -_samples,
                resol, double(mkl_nthr), init_time,
                acpu_time, agpu_time, awrk_time,
                fcpu_time, fgpu_time, fwrk_time);

  return EXIT_SUCCESS;
}
