#include "common.hpp"

#include "cusolverDn.h"

#ifdef USE_COMPLEX
#ifdef USE_FLOAT
#define cusolverDnXpotrf_bufferSize cusolverDnCpotrf_bufferSize
#define cusolverDnXpotrf cusolverDnCpotrf
#define cusolverDnX cuComplex
#else // USE_DOUBLE
#define cusolverDnXpotrf_bufferSize cusolverDnZpotrf_bufferSize
#define cusolverDnXpotrf cusolverDnZpotrf
#define cusolverDnX cuDoubleComplex
#endif // USE_FLOAT
#else // USE_REAL
#ifdef USE_FLOAT
#define cusolverDnXpotrf_bufferSize cusolverDnSpotrf_bufferSize
#define cusolverDnXpotrf cusolverDnSpotrf
#define cusolverDnX float
#else // USE_DOUBLE
#define cusolverDnXpotrf_bufferSize cusolverDnDpotrf_bufferSize
#define cusolverDnXpotrf cusolverDnDpotrf
#define cusolverDnX double
#endif // USE_FLOAT
#endif // USE_COMPLEX

static int Nmin = 0, Nmax = 0, Nstep = 0, _samples = 0, device_ = 0, _devices = 0;
static int lda = 0, Lwork = 0;

static cusolverStatus_t status;
static cusolverDnHandle_t handle;

static cusolverDnX *Agpu = (cusolverDnX*)NULL, *Workspace = (cusolverDnX*)NULL, *Acpu = (cusolverDnX*)NULL;

static void device_count()
{
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
}

static void set_device()
{
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
}

static void create_handle()
{
  status = cusolverDnCreate(&handle); const int lin = __LINE__;
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
}

static void alloc_gpu_mtx()
{
  size_t pitch = 0;
  const cudaError_t error = cudaMallocPitch(&Agpu, &pitch, Nmax * sizeof(cusolverDnX), Nmax); const int lin = __LINE__;
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
  lda = int(pitch / sizeof(cusolverDnX));
#ifndef NDEBUG
  (void)fprintf(stdout, "lda = %d\n", lda);
#endif // !NDEBUG
}

static void find_lwork()
{
  int LworkL = -1, maxLworkL = 0;
  int LworkU = -1, maxLworkU = 0;

  for (int n = Nmin; n <= Nmax; n += Nstep) {
    status = cusolverDnXpotrf_bufferSize(handle, CUBLAS_FILL_MODE_LOWER, n, Agpu, lda, &LworkL); const int lin1 = __LINE__;
    switch (status) {
    case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
      // (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin1);
#endif // !NDEBUG
      break;
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin1);
      exit(status);
    case CUSOLVER_STATUS_INVALID_VALUE:
      (void)fprintf(stderr, "[%s@%s:%d] INVALID_VALUE\n", __FUNCTION__, __FILE__, lin1);
      exit(status);
    case CUSOLVER_STATUS_ARCH_MISMATCH:
      (void)fprintf(stderr, "[%s@%s:%d] ARCH_MISMATCH\n", __FUNCTION__, __FILE__, lin1);
      exit(status);
    case CUSOLVER_STATUS_INTERNAL_ERROR:
      (void)fprintf(stderr, "[%s@%s:%d] INTERNAL_ERROR\n", __FUNCTION__, __FILE__, lin1);
      exit(status);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin1, status);
      exit(status);
    }
    if (LworkL > maxLworkL)
      maxLworkL = LworkL;
    status = cusolverDnXpotrf_bufferSize(handle, CUBLAS_FILL_MODE_UPPER, n, Agpu, lda, &LworkU); const int lin2 = __LINE__;
    switch (status) {
    case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
      // (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin2);
#endif // !NDEBUG
      break;
    case CUSOLVER_STATUS_NOT_INITIALIZED:
      (void)fprintf(stderr, "[%s@%s:%d] NOT_INITIALIZED\n", __FUNCTION__, __FILE__, lin2);
      exit(status);
    case CUSOLVER_STATUS_INVALID_VALUE:
      (void)fprintf(stderr, "[%s@%s:%d] INVALID_VALUE\n", __FUNCTION__, __FILE__, lin2);
      exit(status);
    case CUSOLVER_STATUS_ARCH_MISMATCH:
      (void)fprintf(stderr, "[%s@%s:%d] ARCH_MISMATCH\n", __FUNCTION__, __FILE__, lin2);
      exit(status);
    case CUSOLVER_STATUS_INTERNAL_ERROR:
      (void)fprintf(stderr, "[%s@%s:%d] INTERNAL_ERROR\n", __FUNCTION__, __FILE__, lin2);
      exit(status);
    default:
      (void)fprintf(stderr, "[%s@%s:%d] unknown error %d\n", __FUNCTION__, __FILE__, lin2, status);
      exit(status);
    }
    if (LworkU > maxLworkU)
      maxLworkU = LworkU;
  }

  Lwork = ((maxLworkL >= maxLworkU) ? maxLworkL : maxLworkU);
#ifndef NDEBUG
  (void)fprintf(stdout, "Lwork = %d, LworkL = %d, LworkU = %d\n", Lwork, maxLworkL, maxLworkU);
#endif // !NDEBUG
}

static void alloc_gpu_wrk()
{
  if (Lwork > 0) {
    const cudaError_t error = cudaMalloc(&Workspace, Lwork * sizeof(cusolverDnX)); const int lin = __LINE__;
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
  }
}

static void alloc_cpu_mtx()
{
  const size_t size = size_t(lda) * Nmax * sizeof(cusolverDnX);
  if (size > 0) {
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
    (void)memset(Acpu, 0, size);
  }
}

static void free_cpu_mtx()
{
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
    Acpu = (cusolverDnX*)NULL;
  }
}

static void free_gpu_wrk()
{
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
    Workspace = (cusolverDnX*)NULL;
  }
}

static void free_gpu_mtx()
{
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
  Agpu = (cusolverDnX*)NULL;
}

static void destroy_handle()
{
  status = cusolverDnDestroy(handle); const int lin = __LINE__;
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

  device_count();
  if (device_ > _devices) {
    (void)fprintf(stderr, "device# == %d > #devices == %d\n", device_, _devices);
    return EXIT_FAILURE;
  }
  set_device();

  create_handle();
  alloc_gpu_mtx();
  find_lwork();
  alloc_gpu_wrk();
  alloc_cpu_mtx();
  free_cpu_mtx();
  free_gpu_wrk();
  free_gpu_mtx();
  destroy_handle();

  return EXIT_SUCCESS;
}
