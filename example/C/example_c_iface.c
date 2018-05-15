#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <spllt_iface.h>

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
  int *ptr    = NULL;
  int *row    = NULL;
  double *val = NULL;
  int *order  = NULL;
  double *x   = NULL;
  double *rhs = NULL;
  int n, nnz, nrhs;
  spllt_inform_t info;
  spllt_options_t options = SPLLT_OPTIONS_NULL();
  int stat;
  CPLM_Mat_CSR_t A      = CPLM_MatCSRNULL();
  CPLM_Mat_CSR_t U      = CPLM_MatCSRNULL();
  int ierr      = 0;
  int rank      = 0;
  int matrixSet = 0;
  int ncpu      = 1;
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

  options.ncpu = ncpu;

  if(!matrixSet)
  {
    CPLM_Abort("Error, you have to provide a matrix to factor");
  }

  ierr = CPLM_LoadMatrixMarket(matrixFileName, &A);CPLM_CHKERR(ierr);

  ierr = CPLM_MatCSRGetTriU(&A, &U);
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
  {
  spllt_analyse(&akeep, &fkeep, &options, n, ptr,
      row, &info, order);

  spllt_factor(akeep, fkeep, &options, nnz, val, &info);

  spllt_prepare_solve(akeep, fkeep, &info);

  spllt_solve(fkeep, &options, order, nrhs, x, &info, 1);

  spllt_solve(fkeep, &options, order, nrhs, x, &info, 2);

  spllt_chkerr(n, ptr, row, val, nrhs, x, rhs);
  }
  spllt_deallocate_akeep(&akeep, &stat);
  spllt_deallocate_fkeep(&fkeep, &stat);

  CPLM_MatCSRFree(&A);
  free(x);
  free(rhs);
  free(order);
CPLM_CLOSE_TIMER
CPLM_END_TIME
CPLM_POP
//CPLM_Finalize();
  return 0;
}
