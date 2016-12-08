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

static int Nmin, Nmax, Nstep, _samples;
static int lda, Lwork;

static cusolverStatus_t status;
static cusolverDnHandle_t handle;

static cusolverDnX *Adev, *Workspace;

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
  const cudaError_t error = cudaMallocPitch(&Adev, &pitch, Nmax * sizeof(cusolverDnX), Nmax); const int lin = __LINE__;
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
  int LworkL = -1;
  int LworkU = -1;

  status = cusolverDnXpotrf_bufferSize(handle, CUBLAS_FILL_MODE_LOWER, Nmax, Adev, lda, &LworkL); const int lin1 = __LINE__;
  switch (status) {
  case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin1);
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

  status = cusolverDnXpotrf_bufferSize(handle, CUBLAS_FILL_MODE_UPPER, Nmax, Adev, lda, &LworkU); const int lin2 = __LINE__;
  switch (status) {
  case CUSOLVER_STATUS_SUCCESS:
#ifndef NDEBUG
    (void)fprintf(stdout, "[%s@%s:%d] SUCCESS\n", __FUNCTION__, __FILE__, lin2);
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

  Lwork = ((LworkL >= LworkU) ? LworkL : LworkU);
#ifndef NDEBUG
  (void)fprintf(stdout, "Lwork = %d, LworkL = %d, LworkU = %d\n", Lwork, LworkL, LworkU);
#endif // !NDEBUG  
}

static void alloc_gpu_wrk()
{
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

static void free_gpu_wrk()
{
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

static void free_gpu_mtx()
{
  const cudaError_t error = cudaFree(Adev); const int lin = __LINE__;
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
  Adev = (cusolverDnX*)NULL;
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

  create_handle();
  alloc_gpu_mtx();
  find_lwork();
  alloc_gpu_wrk();
  free_gpu_wrk();
  free_gpu_mtx();
  destroy_handle();

  return EXIT_SUCCESS;
}
