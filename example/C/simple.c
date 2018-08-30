#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <spllt_iface.h>

int main(int argc, char ** argv){

  void *akeep = NULL;
  void *fkeep = NULL;
  int *ptr    = NULL;
  int *row    = NULL;
  double *val = NULL;
  int *order  = NULL;
  double *x   = NULL;
  double *rhs = NULL;
  double *y   = NULL;
  double *workspace = NULL;
  int n, nnz, nrhs, nb;
  long worksize;
  spllt_inform_t info;
  spllt_options_t options = SPLLT_OPTIONS_NULL();
  int stat;

    // Create the matrix 
  // [  2 -1  0 ]
  // [ -1  2 -1 ]
  // [  0 -1  2 ]
  n   = 3;
  nnz = 5;
  ptr = malloc((n+1) * sizeof(int));
  row = malloc(nnz * sizeof(int));
  val = malloc(nnz * sizeof(double));

  options.nb = nb = 4;

  ptr[0] = 1; ptr[1] = 3; ptr[2] = 5; ptr[3] = 6;
  row[0] = 1; row[1] = 2; row[2] = 2; row[3] = 3; row[4] = 3;
  for(int i = 0; i < nnz; i++) val[i] = 2.0;
  val[1] = - 1.0; val[3] = - 1.0;

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
  spllt_analyse(&akeep, &fkeep, &options, n, ptr,
      row, &info, order);

  spllt_factor(akeep, fkeep, &options, nnz, val, &info);
  spllt_wait();

  spllt_prepare_solve(akeep, fkeep, nb, nrhs, &worksize, &info);
  printf("Need a workspace of size %ld\n", worksize);

  y         = calloc( n * nrhs, sizeof(double));
  workspace = calloc( worksize, sizeof(double));

  spllt_set_mem_solve(akeep, fkeep, nb, nrhs, worksize, y, workspace, &info);

  spllt_solve(fkeep, &options, order, nrhs, x, &info, 6);
  spllt_wait();

//spllt_solve(fkeep, &options, order, nrhs, x, &info, 2);

  spllt_chkerr(n, ptr, row, val, nrhs, x, rhs);
  }
  spllt_deallocate_akeep(&akeep, &stat);
  spllt_deallocate_fkeep(&fkeep, &stat);

  return 0;
}
