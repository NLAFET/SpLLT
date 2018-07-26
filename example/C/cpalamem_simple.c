#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <spllt_iface.h>
#include <omp.h>

#include <cpalamem_macro.h>
#include <mat_csr.h>
#include <mat_load_mm.h>
#include <cpalamem_handler.h>
#include <cpalamem_instrumentation.h>

#define USAGE "Usage %s -m <matrixFileName>"
#define FILENAMESIZE 256

int main(int argc, char ** argv){

  void *akeep = NULL;
  void *fkeep = NULL;
  void *tm    = NULL;
  int *ptr    = NULL;
  int *row    = NULL;
  double *val = NULL;
  int *order  = NULL;
  double *x   = NULL;
  double *rhs = NULL;
  double *y   = NULL;
  double *workspace = NULL;
  int n, nnz, nrhs;
  int nb = 32;
  long worksize;
  spllt_inform_t info;
  spllt_options_t options = SPLLT_OPTIONS_NULL();
  int stat;
  CPLM_Mat_CSR_t A      = CPLM_MatCSRNULL();
  CPLM_Mat_CSR_t U      = CPLM_MatCSRNULL();
  int ierr      = 0;
  int rank      = 0;
  int matrixSet = 0;
  int ncpu      = 1;
  long size     = 0;
  char matrixFileName[FILENAMESIZE];

  CPLM_Init(&argc, &argv);
CPLM_PUSH
CPLM_BEGIN_TIME
CPLM_OPEN_TIMER

  //Handle parameters
  if(argc>0)
  {
    for(int i = 1; i < argc; i++)
    {
      if(!strcmp(argv[i],"-h"))
      {
        CPLM_Abort(USAGE,argv[0]);
      }
      else if(!strcmp(argv[i],"-m"))
      {
        i++;
        strcpy(matrixFileName, argv[i]);
        matrixSet = 1;
      }
      else{
        if(!rank)
        {
          CPLM_Abort(USAGE,argv[0]);
        }
      }
    }
  }

  options.ncpu  = ncpu;
  options.nb    = nb;

  if(!matrixSet)
  {
    CPLM_Abort("Error, you have to provide a matrix to factor");
  }

  ierr = CPLM_LoadMatrixMarket(matrixFileName, &A);CPLM_CHKERR(ierr);

  ierr = CPLM_MatCSRGetTriU(&A, &U);CPLM_CHKERR(ierr);
  CPLM_MatCSRConvertTo1BasedIndexing(&U);

  nnz = U.info.nnz;
  n   = U.info.n;
  ptr = U.rowPtr;
  row = U.colInd;
  val = U.val;

  order = malloc(n * sizeof(int));

  //Create RHS
  nrhs = 1;
  x     = malloc(n * nrhs * sizeof(double));
  rhs   = malloc(n * nrhs * sizeof(double));
  for(int i = 0; i < n; i++) rhs[i] = 1.0;
  memcpy(x, rhs, n * sizeof(double));

  #pragma omp parallel 
  #pragma omp single
  {
    spllt_task_manager_init(&tm);

    printf("Analysis\n");
    spllt_analyse(&akeep, &fkeep, &options, n, ptr,
        row, &info, order);

    printf("Factor\n");
    spllt_factor(akeep, fkeep, &options, nnz, val, &info);
    spllt_wait();

    printf("Prepare solve\n");
    spllt_prepare_solve(akeep, fkeep, nb, nrhs, &worksize, &info);
    printf("Need a workspace of size %ld\n", worksize);
    printf("Need a y vector of size %d\n", n * nrhs);

    y         = calloc( n * nrhs, sizeof(double));
    workspace = calloc( worksize, sizeof(double));

    spllt_set_mem_solve(akeep, fkeep, nb, nrhs, worksize, y, workspace, &info);
#if 0

    spllt_solve_worker(fkeep, &options, order, nrhs, x, &info, 1, workspace, size, tm);

    spllt_solve_worker(fkeep, &options, order, nrhs, x, &info, 2, workspace, size, tm);
#else
    spllt_solve(fkeep, &options, order, nrhs, x, &info, 7);
    spllt_solve(fkeep, &options, order, nrhs, x, &info, 8);
    spllt_wait();

  //spllt_solve(fkeep, &options, order, nrhs, x, &info, 2);
#endif
    spllt_chkerr(n, ptr, row, val, nrhs, x, rhs);
  }
  spllt_deallocate_akeep(&akeep, &stat);
  spllt_deallocate_fkeep(&fkeep, &stat);
  spllt_task_manager_deallocate(&tm, &stat);

  CPLM_MatCSRFree(&A);
  free(x);
  free(rhs);
  free(order);
  free(workspace);
CPLM_CLOSE_TIMER
CPLM_END_TIME
CPLM_POP
//CPLM_Finalize();
  return 0;
}
