#include <starpu.h>

#if defined(SPLLT_USE_GPU)
#include <starpu_cuda.h>
#include <magma.h>
/* #include <magmablas.h> */
#include <cublas.h>
/* #include <cusolverDn.h> */
#endif
/* void spllt_starpu_unpack_args_factorize_block(void *cl_arg, */
/*                                               int *m, int *n) { */
   
/*    starpu_codelet_unpack_args(cl_arg, */
/*                               m, n); */

/*   return; */
/* } */

// factorize block StarPU task insert

void spllt_starpu_factorize_block_cpu_func(void *buffers[], void *cl_arg);

#if defined(SPLLT_USE_GPU)
void spllt_starpu_factorize_block_cuda_func(void *buffers[], void *cl_arg) {

   int info;
   /* starpu_codelet_unpack_args(cl_arg, ); */
   /* printf("[spllt_starpu_factorize_block_cuda_func]\n"); */

   unsigned m = STARPU_MATRIX_GET_NX(buffers[0]);
   unsigned n = STARPU_MATRIX_GET_NY(buffers[0]);
   unsigned ld = STARPU_MATRIX_GET_LD(buffers[0]);

   double *dest = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);

   cudaStream_t local_stream = starpu_cuda_get_local_stream();

   /* cusolverDnHandle_t cu_hdl; */
   /* cusolverDnCreate(&cu_hdl); */
   /* cusolverDnDestroy(cu_hdl); */

   magmablasSetKernelStream(local_stream);

   /* magma_dpotrf_gpu(MagmaUpper, n, dest, n, &info); */
   magma_dpotf2_gpu(MagmaUpper, n, dest, n, &info);

   if (m > n) {
      /* printf("TET\n"); */
      magma_dtrsm(MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit,
                  n, m-n, 1.0, dest, n, dest + n*n, n);
   }

  /* cudaStreamSynchronize(local_stream); */
}
#endif


// factorize block task codelet
struct starpu_codelet cl_factorize_block = {
#if defined(SPLLT_USE_GPU)
  .where = STARPU_CUDA /* | STARPU_CPU */,
  /* .where = STARPU_CPU, */
  .cuda_flags = {STARPU_CUDA_ASYNC},
  .cuda_funcs = {spllt_starpu_factorize_block_cuda_func, NULL},
#else
  .where = STARPU_CPU,
#endif
  .cpu_funcs = {spllt_starpu_factorize_block_cpu_func, NULL},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .name = "FACTO_BLK"
};

void spllt_starpu_insert_factorize_block_c(starpu_data_handle_t l_handle,
                                           starpu_data_handle_t node_handle,
                                           int prio) {
   
   int ret;
   /* printf("[spllt_starpu_insert_factorize_block_c]\n"); */
   if (node_handle) {
      /* printf("Test\n"); */
      ret = starpu_task_insert(&cl_factorize_block,                           
                               STARPU_RW,	l_handle,
                               STARPU_R,	node_handle,
                               STARPU_PRIORITY, prio,
                               0);
   }
   else {
      ret = starpu_task_insert(&cl_factorize_block,                           
                               STARPU_RW,	l_handle,
                               STARPU_PRIORITY, prio,
                               0);
   }
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");
}

// solve block StarPU task insert

void spllt_starpu_solve_block_cpu_func(void *buffers[], void *cl_arg);

#if defined(SPLLT_USE_GPU)
void spllt_starpu_solve_block_cuda_func(void *buffers[], void *cl_arg) {

   /* printf("[spllt_starpu_solve_block_cuda_func]\n"); */

   double *bc_kk = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
 
   unsigned m = STARPU_MATRIX_GET_NX(buffers[1]);
   unsigned n = STARPU_MATRIX_GET_NY(buffers[1]);
   unsigned ld = STARPU_MATRIX_GET_LD(buffers[1]);
   double *bc_ik = (double *)STARPU_MATRIX_GET_PTR(buffers[1]);

   cudaStream_t local_stream = starpu_cuda_get_local_stream();

   magmablasSetKernelStream(local_stream);

   magma_dtrsm(MagmaLeft, MagmaUpper, MagmaTrans, MagmaNonUnit,
               n, m, 1.0, bc_kk, n, bc_ik, n);
  
}
#endif

// solve block task codelet
struct starpu_codelet cl_solve_block = {
#if defined(SPLLT_USE_GPU)
   /* .where = STARPU_CPU, */
   .where = STARPU_CUDA,
   .cuda_flags = {STARPU_CUDA_ASYNC},
   .cuda_funcs = {spllt_starpu_solve_block_cuda_func, NULL},
#else
   .where = STARPU_CPU,
#endif
   .cpu_funcs = {spllt_starpu_solve_block_cpu_func, NULL},
   .nbuffers = STARPU_VARIABLE_NBUFFERS,
   .name = "SOLVE_BLK"
};

void spllt_starpu_insert_solve_block_c(starpu_data_handle_t lkk_handle,
                                       starpu_data_handle_t lik_handle,
                                       int prio) {
   
   int ret;

   ret = starpu_task_insert(&cl_solve_block,                           
                            STARPU_R,	     lkk_handle,
                            STARPU_RW,	     lik_handle,
                            STARPU_PRIORITY, prio,
                            0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");
}

// update block StarPU task insert

void spllt_starpu_update_block_cpu_func(void *buffers[], void *cl_arg);

#if defined(SPLLT_USE_GPU)
void spllt_starpu_update_block_cuda_func(void *buffers[], void *cl_arg) {

   /* printf("[spllt_starpu_update_block_cuda_func]\n"); */
   
   int diag;

   starpu_codelet_unpack_args(cl_arg,
                              &diag);

   double *bc_ik = (double *)STARPU_MATRIX_GET_PTR(buffers[0]);
   unsigned n1 = STARPU_MATRIX_GET_NY(buffers[1]);
   double *bc_jk = (double *)STARPU_MATRIX_GET_PTR(buffers[1]);

   unsigned m = STARPU_MATRIX_GET_NX(buffers[2]);
   unsigned n = STARPU_MATRIX_GET_NY(buffers[2]);
   unsigned ld = STARPU_MATRIX_GET_LD(buffers[2]);
   double *bc_ij = (double *)STARPU_MATRIX_GET_PTR(buffers[2]);

   cudaStream_t local_stream = starpu_cuda_get_local_stream();

   magmablasSetKernelStream(local_stream);

   if (diag) {

      magma_dsyrk(MagmaUpper, MagmaTrans, 
                  n, n1, -1.0, bc_jk, n1, 1.0, bc_ij, n);

      if (m > n) {
         
         magma_dgemm(MagmaTrans, MagmaNoTrans, 
                     n, m-n, n1, -1.0, 
                     bc_jk, n1, bc_ik + n*n1, n1, 1.0, bc_ij + n*n, n);
      }
   }
   else {
      
         magma_dgemm(MagmaTrans, MagmaNoTrans, 
                     n, m, n1, -1.0, 
                     bc_jk, n1, bc_ik, n1, 1.0, bc_ij, n);      
   }   
}
#endif

void spllt_starpu_codelet_unpack_args_update_block(void *cl_arg,
                                                   int *diag) {
   
   starpu_codelet_unpack_args(cl_arg,
                              diag);
}

// update block task codelet
struct starpu_codelet cl_update_block = {
#if defined(SPLLT_USE_GPU)
   .where = STARPU_CUDA,
   /* .where = STARPU_CPU, */
   .cuda_flags = {STARPU_CUDA_ASYNC},
   .cuda_funcs = {spllt_starpu_update_block_cuda_func, NULL},
#else
   .where = STARPU_CPU,
#endif
   .cpu_funcs = {spllt_starpu_update_block_cpu_func, NULL},
   .nbuffers = STARPU_VARIABLE_NBUFFERS,
   .name = "UPDATE_BLK"
};

void spllt_starpu_insert_update_block_c(starpu_data_handle_t lik_handle,
                                        starpu_data_handle_t ljk_handle,
                                        starpu_data_handle_t lij_handle,
                                        int diag, int prio) {
   
   int ret;

   ret = starpu_task_insert(&cl_update_block,                           
                            STARPU_VALUE,    &diag, sizeof(int),
                            STARPU_R,	     lik_handle,
                            STARPU_R,	     ljk_handle,
                            STARPU_RW | STARPU_COMMUTE,	     lij_handle,
                            STARPU_PRIORITY, prio,
                            0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

}

// update between StarPU task insert

void spllt_starpu_update_between_cpu_func(void *buffers[], void *cl_arg);

#if defined(SPLLT_USE_GPU)
/* void spllt_starpu_update_between_cuda_func(void *buffers[], void *cl_arg); */

void spllt_starpu_update_between_cuda_func(void *buffers[], void *cl_arg) {

   void *snode;
   int scol;
   void *anode;
   void *a_blk;
   int dcol;
   int csrc, csrc2;
   int rsrc, rsrc2;
   int min_with_blas;

   /* printf("[spllt_starpu_update_between_cuda_func]\n"); */

   starpu_codelet_unpack_args(cl_arg,
                              &snode, &scol, &anode, &a_blk, &dcol,
                              &csrc, &csrc2,
                              &rsrc, &rsrc2,
                              &min_with_blas);

   /* /\* printf("[spllt_starpu_update_between_cuda_func] csrc: %d, rsrc: %d\n", csrc, rsrc); *\/ */

   // workspace pointer
   double *d_buffer = (double *)STARPU_VECTOR_GET_PTR(buffers[0]);
   unsigned b_sz    = STARPU_VECTOR_GET_NX(buffers[0]);

   // row_list pointer
   int *d_row_list = (int *)STARPU_VECTOR_GET_PTR(buffers[1]);
   unsigned d_row_list_sz    = STARPU_VECTOR_GET_NX(buffers[1]);

   // col_list pointer
   int *d_col_list = (int *)STARPU_VECTOR_GET_PTR(buffers[2]);
   unsigned d_col_list_sz = STARPU_VECTOR_GET_NX(buffers[2]);

   // A_ij pointer
   double *d_bc_ij = (double *)STARPU_MATRIX_GET_PTR(buffers[3]);
   unsigned m    = STARPU_MATRIX_GET_NX(buffers[3]);
   unsigned n    = STARPU_MATRIX_GET_NY(buffers[3]);
   unsigned ld   = STARPU_MATRIX_GET_LD(buffers[3]);

   // L_ik pointer
   double *d_bc_ik = (double *)STARPU_MATRIX_GET_PTR(buffers[4]);   
   /* unsigned n1    = STARPU_MATRIX_GET_NY(buffers[2]); */
   unsigned m2    = STARPU_MATRIX_GET_NX(buffers[4]);
      
   // L_jk pointer
   double *d_bc_jk = (double *)STARPU_MATRIX_GET_PTR(buffers[5]);
   unsigned n1    = STARPU_MATRIX_GET_NY(buffers[5]);
   unsigned m1    = STARPU_MATRIX_GET_NX(buffers[5]);

   int s1sa, s1en, s2sa, s2en;
   int *col_list, *row_list;
   int cls, rls;

   cudaStream_t local_stream = starpu_cuda_get_local_stream();
   int worker_id = starpu_worker_get_id();
   unsigned worker_node = starpu_worker_get_memory_node(worker_id);

   
   double *src1, *src2, *a, *buff;

   row_list = (int *) malloc(m*sizeof(int));
   col_list = (int *) malloc(n*sizeof(int));

   spllt_update_between_compute_map_c(a_blk, dcol, anode, scol, snode,
                                      row_list, col_list, &rls, &cls,
                                      &s1sa, &s1en, &s2sa, &s2en);
   
   /* a    = (double *)malloc(m*n*sizeof(double));   // A_ij */
   /* src1 = (double *)malloc(m1*n1*sizeof(double)); // L_ik */
   /* src2 = (double *)malloc(m2*n1*sizeof(double)); // L_jk */
   /* buff  = (double *)malloc(b_sz*sizeof(double)); // buffer */

   magmablasSetKernelStream(local_stream);

   int diag = is_blk_diag(a_blk);
   int ndiag = 0;

   if (diag) {
      ndiag = s1en-s1sa+1;
      magma_dsyrk(MagmaUpper, MagmaTrans,
                  ndiag, n1, -1.0, &d_bc_jk[csrc-1], n1,
                  0.0,
                  d_buffer, cls);

      if (s2en-s2sa+1-ndiag > 0) {
         magma_dgemm(MagmaTrans, MagmaNoTrans,
                     ndiag, s2en-s2sa+1-ndiag, n1, -1.0,
                     &d_bc_jk[csrc-1], n1,
                     &d_bc_ik[rsrc-1 + n1*ndiag], n1, 0.0,
                     &d_buffer[cls*ndiag], cls);
      }
   }
   else {
      
      magma_dgemm(MagmaTrans, MagmaNoTrans,
                  s1en-s1sa+1, s2en-s2sa+1, n1,
                  -1.0,
                  &d_bc_jk[csrc-1], n1,
                  &d_bc_ik[rsrc-1], n1,
                  0.0,
                  d_buffer, cls);
   }
   /* printf("diag: %d, ndiag: %d\n", diag, ndiag); */
   /* cudaStreamSynchronize(local_stream); */

   /* starpu_cuda_copy_async_sync(d_bc_ij, worker_node, a, 0, m*n*sizeof(double), local_stream, cudaMemcpyDeviceToHost); */
   /* starpu_cuda_copy_async_sync(d_bc_jk, worker_node, src1, 0, m1*n1*sizeof(double), local_stream, cudaMemcpyDeviceToHost); */
   /* starpu_cuda_copy_async_sync(d_bc_ik, worker_node, src2, 0, m2*n1*sizeof(double), local_stream, cudaMemcpyDeviceToHost); */
   
   /* spllt_update_between_c(m, n, a_blk, dcol, anode, n1, scol, snode, a, src1 + csrc-1, src2 + rsrc-1, buf, 0); */
   
   /* starpu_cuda_copy_async_sync(a, 0, d_bc_ij, worker_node, m*n*sizeof(double), local_stream, cudaMemcpyHostToDevice); */
   /* starpu_cuda_copy_async_sync(buf, 0, d_buffer, worker_node, b_sz*sizeof(double), local_stream, cudaMemcpyHostToDevice); */


   /* starpu_cuda_copy_async_sync(d_bc_ij, worker_node, a, 0, m*n*sizeof(double), local_stream, cudaMemcpyDeviceToHost); */
   /* starpu_cuda_copy_async_sync(d_buffer, worker_node, buff, 0, b_sz*sizeof(double), local_stream, cudaMemcpyDeviceToHost); */

   /* spllt_expand_buffer_c(a, m, n, row_list, rls, col_list, cls, ndiag, buff); */
   /* int i, j, imax; */
   
   /* for (i = 0; i < rls; i++) { */
   /*    int row = row_list[i]-1; */
   /*    imax = cls-1; */
   /*    if (i < ndiag) imax = i; */
   /*    for (j = 0; j <= imax; j++) { */
   /*       int col = col_list[j]-1; */
   /*       a[row*n + col] += buff[i*cls + j];          */
   /*    } */
   /* } */

   /* starpu_cuda_copy_async_sync(a, 0, d_bc_ij, worker_node, m*n*sizeof(double), local_stream, cudaMemcpyHostToDevice); */

   /* int *d_row_list = starpu_malloc_on_node(worker_node, rls*sizeof(int)); */
   /* int *d_col_list = starpu_malloc_on_node(worker_node, cls*sizeof(int)); */

   starpu_cuda_copy_async_sync(row_list, 0, d_row_list, worker_node, rls*sizeof(int), local_stream, cudaMemcpyHostToDevice);
   starpu_cuda_copy_async_sync(col_list, 0, d_col_list, worker_node, cls*sizeof(int), local_stream, cudaMemcpyHostToDevice);

   /* cudaStreamSynchronize(local_stream); */

   spllt_cu_expand_buffer(n, d_bc_ij, d_buffer,
                          d_row_list, rls, d_col_list, cls,
                          ndiag,
                          local_stream);

   /* cudaStreamSynchronize(local_stream); */

}
#endif

#if defined(SPLLT_USE_GPU)

void spllt_starpu_codelet_unpack_args_update_between(void *cl_arg,
                                                     void *snode, int *scol, 
                                                     void *anode, void *a_blk, int *dcol,
                                                     int *csrc, int *csrc2,
                                                     int *rsrc, int *rsrc2,
                                                     int *min_with_blas) {
   
   starpu_codelet_unpack_args(cl_arg,
                              snode, scol, anode, a_blk, dcol,
                              csrc, csrc2,
                              rsrc, rsrc2,
                              min_with_blas);

   return;
}

#else

void spllt_starpu_codelet_unpack_args_update_between(void *cl_arg,
                                                     void *snode, int *scol, 
                                                     void *anode, void *a_blk, int *dcol,
                                                     int *csrc, int *csrc2,
                                                     int *rsrc, int *rsrc2,
                                                     int *min_with_blas,
                                                     int *nhlik, int *nhljk) {
   
   starpu_codelet_unpack_args(cl_arg,
                              snode, scol, anode, a_blk, dcol,
                              csrc, csrc2, 
                              rsrc, rsrc2,
                              min_with_blas,
                              nhlik, nhljk);

   return;
}
#endif

// update block task codelet
struct starpu_codelet cl_update_between = {
#if defined(SPLLT_USE_GPU)
   .where =  STARPU_CUDA,
   /* .where =  /\* STARPU_CUDA | *\/ STARPU_CPU, */
   .cuda_flags = {STARPU_CUDA_ASYNC},
   .cuda_funcs = {spllt_starpu_update_between_cuda_func, NULL},
#else
   .where = STARPU_CPU,
#endif
   .cpu_funcs = {spllt_starpu_update_between_cpu_func, NULL},
   .nbuffers = STARPU_VARIABLE_NBUFFERS,
   .name = "UPDATE_BETWEEN"
};

#if defined(SPLLT_USE_GPU)

void spllt_starpu_insert_update_between_c(starpu_data_handle_t lik_handle,
                                          starpu_data_handle_t ljk_handle,
                                          starpu_data_handle_t lij_handle,
                                          void *snode, int scol,
                                          void *anode, void *a_blk, int dcol,
                                          int csrc, int csrc2, int rsrc, int rsrc2,
                                          int min_width_blas,
                                          starpu_data_handle_t workspace_handle,
                                          starpu_data_handle_t row_list_hdl,
                                          starpu_data_handle_t col_list_hdl,
                                          starpu_data_handle_t node_handle,
                                          int prio) {

   
   int ret;

   ret = starpu_task_insert(&cl_update_between,
                            STARPU_VALUE, &snode, sizeof(void *),
                            STARPU_VALUE, &scol, sizeof(int),
                            STARPU_VALUE, &anode, sizeof(void *),
                            STARPU_VALUE, &a_blk, sizeof(void *),
                            STARPU_VALUE, &dcol, sizeof(int),
                            STARPU_VALUE, &csrc, sizeof(int),
                            STARPU_VALUE, &csrc2, sizeof(int),
                            STARPU_VALUE, &rsrc, sizeof(int),
                            STARPU_VALUE, &rsrc2, sizeof(int),
                            STARPU_VALUE, &min_width_blas, sizeof(int),
                            STARPU_SCRATCH, workspace_handle,
                            STARPU_SCRATCH, row_list_hdl,
                            STARPU_SCRATCH, col_list_hdl,
                            STARPU_RW | STARPU_COMMUTE, lij_handle,
                            STARPU_R, lik_handle,
                            STARPU_R, ljk_handle,
                            STARPU_PRIORITY, prio,
                            0);

}

#else

void spllt_starpu_insert_update_between_c(starpu_data_handle_t *lik_handles, int nhlik,
                                          starpu_data_handle_t *ljk_handles, int nhljk,
                                          starpu_data_handle_t lij_handle,
                                          void *snode, int scol,
                                          void *anode, void *a_blk, int dcol,                                          
                                          int csrc, int csrc2, int rsrc, int rsrc2,
                                          int min_width_blas,
                                          starpu_data_handle_t workspace_handle,
                                          starpu_data_handle_t row_list_hdl,
                                          starpu_data_handle_t col_list_hdl,
                                          starpu_data_handle_t node_handle,
                                          int prio) {

   int ret, i, nh;
   struct starpu_data_descr *descrs;

   nh = 0;
   descrs = malloc((nhlik+nhljk+5) * sizeof(struct starpu_data_descr));

   // worksapce 
   descrs[nh].handle =  workspace_handle; descrs[nh].mode = STARPU_SCRATCH;
   nh = nh + 1;

   // row_list
   descrs[nh].handle =  row_list_hdl; descrs[nh].mode = STARPU_SCRATCH;
   nh = nh + 1;

   // col_list
   descrs[nh].handle =  col_list_hdl; descrs[nh].mode = STARPU_SCRATCH;
   nh = nh + 1;

   // A_ij 
   descrs[nh].handle =  lij_handle; descrs[nh].mode = STARPU_RW | STARPU_COMMUTE;
   nh = nh + 1;

   for(i=0; i<nhlik; i++) {
      descrs[nh].handle = lik_handles[i];  descrs[nh].mode = STARPU_R;
      nh = nh + 1;
   }

   for(i=0; i<nhljk; i++){
      descrs[nh].handle = ljk_handles[i];  descrs[nh].mode = STARPU_R;
      nh = nh + 1;
   }

   descrs[nh].handle = node_handle;  descrs[nh].mode = STARPU_R;
   nh = nh + 1;

   ret = starpu_task_insert(&cl_update_between,
                            STARPU_VALUE, &snode, sizeof(void *),
                            STARPU_VALUE, &scol, sizeof(int),
                            STARPU_VALUE, &anode, sizeof(void *),
                            STARPU_VALUE, &a_blk, sizeof(void *),
                            STARPU_VALUE, &dcol, sizeof(int),
                            STARPU_VALUE, &csrc, sizeof(int),
                            STARPU_VALUE, &csrc2, sizeof(int),
                            STARPU_VALUE, &rsrc, sizeof(int),
                            STARPU_VALUE, &rsrc2, sizeof(int),
                            STARPU_VALUE, &min_width_blas, sizeof(int),
                            STARPU_VALUE, &nhlik, sizeof(int),
                            STARPU_VALUE, &nhljk, sizeof(int),
                            STARPU_DATA_MODE_ARRAY, descrs,   nh,
                            STARPU_PRIORITY, prio,
                            0);
   
   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

   return;
}

#endif

// tests

void test_cpu_func(void *buffers[], void *cl_arg) {

   /* printf("TEST\n"); */

   return;
}

struct starpu_codelet cl_test = {
  .where = STARPU_CPU,
  .cpu_funcs = {test_cpu_func, NULL},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .name = "TEST"
};

void test_insert_c(starpu_data_handle_t *a_handles, int nah,
                   starpu_data_handle_t *b_handles, int nbh) {

   int ret, i, nh;
   struct starpu_data_descr *descrs;
   printf("test insert nah: %d, nbh: %d\n", nah, nbh);
   descrs = malloc((nah+nbh) * sizeof(struct starpu_data_descr));
   nh = 0;
   for(i=0; i<nah; i++){
      descrs[nh].handle = a_handles[i];  descrs[nh].mode = STARPU_R /* STARPU_RW */;
      nh = nh + 1;
   }

   for(i=0; i<nbh; i++){
      descrs[nh].handle = b_handles[i];  descrs[nh].mode = STARPU_R /* STARPU_RW */;
      nh = nh + 1;
   }

   ret = starpu_task_insert(&cl_test,
                            STARPU_DATA_MODE_ARRAY, descrs,   nh,
                            0);

   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

   return;   
}

// partitioning and unpartitioning data

void spllt_starpu_partition_cpu_func(void *buffers[], void *cl_arg){
  /* printf("partition\n"); */
  return;
}

struct starpu_codelet cl_partition = {
  .where = STARPU_CPU,
  .cpu_funcs = {spllt_starpu_partition_cpu_func, NULL},
  .nbuffers = STARPU_VARIABLE_NBUFFERS, 
  .name = "PARTITION",
/* #if defined(usebound) || defined(calibratecpus) || defined(calibrategpus) */
/*   .model = &partition_history_based_model */
/* #elif defined(useperfmodels) */
/*   .model = &partition_history_based_model */
/* #endif */
};


void spllt_insert_partition_task_c(starpu_data_handle_t bc_handle, 
                                   starpu_data_handle_t *in_bc_handles, int nh,
                                   int prio) {
   int ret;
   int i;
   struct starpu_data_descr *descrs;

   descrs[0].handle = bc_handle;        descrs[0].mode = STARPU_R;

   for(i=0; i<nh; i++){
      descrs[i+1].handle = in_bc_handles[i];  descrs[i+1].mode = STARPU_W /* STARPU_RW */;
   }

   ret = starpu_task_insert(&cl_partition,
                            STARPU_DATA_MODE_ARRAY, descrs,   nh+1,
                            STARPU_PRIORITY,        prio,
                            0);
   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");
} 

void spllt_insert_unpartition_task_c(starpu_data_handle_t bc_handle, 
                                     starpu_data_handle_t *in_bc_handles, int nh,
                                     int prio) {

   int ret, i;
   struct starpu_data_descr *descrs;

   descrs = malloc((nh+1) * sizeof(struct starpu_data_descr));
   
   descrs[0].handle = bc_handle;        descrs[0].mode = STARPU_W;
   /* printf("nh: %d\n", nh); */
   for(i=0; i<nh; i++){
      descrs[i+1].handle = in_bc_handles[i];  descrs[i+1].mode = STARPU_R;
   }

   ret = starpu_task_insert(&cl_partition/* &cl_unpartition_hierarchical */,
                            STARPU_DATA_MODE_ARRAY, descrs,   nh+1,
                            STARPU_PRIORITY,        prio,
                            0);
   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

   return;
}

void spllt_starpu_init_node_cpu_func(void *buffers[], void *cl_arg);

void spllt_starpu_codelet_unpack_args_init_node(void *cl_arg,
                                                int *snode,
                                                void *val, int *nval, 
                                                void *keep) {
   
   starpu_codelet_unpack_args(cl_arg,
                              snode,
                              val, nval, keep);
}

struct starpu_codelet cl_init_node = {
  .where = STARPU_CPU,
  .cpu_funcs = {spllt_starpu_init_node_cpu_func, NULL},
  .nbuffers = STARPU_VARIABLE_NBUFFERS, 
  .name = "INIT_NODE",
/* #if defined(usebound) || defined(calibratecpus) || defined(calibrategpus) */
/*   .model = &partition_history_based_model */
/* #elif defined(useperfmodels) */
/*   .model = &partition_history_based_model */
/* #endif */
};

void spllt_insert_init_node_task_c(starpu_data_handle_t node_handle, 
                                   int snode, void *val, int nval, void *keep, 
                                   int prio) {
   
   
   int ret;

   ret = starpu_task_insert(&cl_init_node,
                            STARPU_VALUE, &snode, sizeof(int),
                            STARPU_VALUE, &val, sizeof(void*),
                            STARPU_VALUE, &nval, sizeof(int),
                            STARPU_VALUE, &keep, sizeof(void*),
                            STARPU_RW, node_handle,
                            STARPU_PRIORITY, prio,
                            0);

   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");
}

void spllt_starpu_factorize_node_cpu_func(void *buffers[], void *cl_arg);

struct starpu_codelet cl_factorize_node = {
  .where = STARPU_CPU,
  .cpu_funcs = {spllt_starpu_factorize_node_cpu_func, NULL},
  .nbuffers = STARPU_VARIABLE_NBUFFERS, 
  .name = "FACTO_NODE",
/* #if defined(usebound) || defined(calibratecpus) || defined(calibrategpus) */
/*   .model = &partition_history_based_model */
/* #elif defined(useperfmodels) */
/*   .model = &partition_history_based_model */
/* #endif */
};

void spllt_starpu_codelet_unpack_args_factorize_node(void *cl_arg,
                                                     void *snode, void *fdata, 
                                                     void *keep, void *control) {
   
   starpu_codelet_unpack_args(cl_arg,
                              snode, fdata, 
                              keep, control);
}

void spllt_insert_factorize_node_task_c(starpu_data_handle_t node_hdl,
                                        starpu_data_handle_t *cnode_hdls, int nc,
                                        starpu_data_handle_t map_hdl,
                                        void *snode, void *fdata, 
                                        void *keep, void *control,
                                        int prio) {

   /* printf("[spllt_insert_factorize_node_task_c]\n"); */

   int ret;
   int i, nh;

   struct starpu_data_descr *descrs;

   descrs = malloc((nc+2) * sizeof(struct starpu_data_descr));

   nh = 0;
   
   descrs[nh].handle = map_hdl;  descrs[nh].mode = STARPU_SCRATCH;
   nh = nh + 1;
   
   descrs[nh].handle = node_hdl;  descrs[nh].mode = STARPU_RW;
   nh = nh + 1;
   /* printf("[spllt_insert_factorize_node_task_c] nc: %d\n", nc); */
   for(i=0; i<nc; i++){
      descrs[nh].handle = cnode_hdls[i];  descrs[nh].mode = STARPU_R;
      nh = nh + 1;
   }

   ret = starpu_task_insert(&cl_factorize_node,
                            STARPU_VALUE, &snode, sizeof(void*),
                            STARPU_VALUE, &fdata, sizeof(void*),
                            STARPU_VALUE, &keep, sizeof(void*),
                            STARPU_VALUE, &control, sizeof(void*),
                            STARPU_DATA_MODE_ARRAY, descrs, nh,
                            STARPU_PRIORITY, prio,
                            0);

   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");
}

// subtree_factorize StarPU task insert

struct starpu_codelet cl_subtree_factorize = {
  .where = STARPU_CPU,
  .cpu_funcs = {spllt_starpu_subtree_factorize_cpu_func, NULL},
  .nbuffers = STARPU_VARIABLE_NBUFFERS, 
  .name = "SUBTREE_FACTORIZE",
/* #if defined(usebound) || defined(calibratecpus) || defined(calibrategpus) */
/*   .model = &partition_history_based_model */
/* #elif defined(useperfmodels) */
/*   .model = &partition_history_based_model */
/* #endif */
};

void spllt_insert_subtree_factorize_task_c(int root, 
                                           void *val, void *keep,
                                           starpu_data_handle_t buffer_hdl,
                                           void *cntl,
                                           starpu_data_handle_t map_hdl,
                                           starpu_data_handle_t row_list_hdl,
                                           starpu_data_handle_t col_list_hdl,
                                           starpu_data_handle_t workspace_handle) {

   int ret;
   
   ret = starpu_task_insert(&cl_subtree_factorize,
                            STARPU_VALUE, &root, sizeof(int),
                            STARPU_VALUE, &val, sizeof(void *),
                            STARPU_VALUE, &keep, sizeof(void *),
                            STARPU_VALUE, &cntl, sizeof(void *));
   
   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");
}

void spllt_starpu_codelet_unpack_args_subtree_factorize(void *cl_arg,
                                                        int *root, void *val, 
                                                        void *keep, void *cntl) {
   
   starpu_codelet_unpack_args(cl_arg,
                              root, fdata, 
                              keep, control);
}
