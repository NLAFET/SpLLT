#include <starpu.h>

/* void spllt_starpu_unpack_args_factorize_block(void *cl_arg, */
/*                                               int *m, int *n) { */
   
/*    starpu_codelet_unpack_args(cl_arg, */
/*                               m, n); */

/*   return; */
/* } */

// factorize block StarPU task insert

void spllt_starpu_factorize_block_cpu_func(void *buffers[], void *cl_arg);

// factorize block task codelet
struct starpu_codelet cl_factorize_block = {
  .where = STARPU_CPU,
  .cpu_funcs = {spllt_starpu_factorize_block_cpu_func, NULL},
  .nbuffers = 1,
  .modes = {STARPU_RW},
  .name = "FACTO_BLK"
};

void spllt_starpu_insert_factorize_block_c(starpu_data_handle_t l_handle,
                                           int prio) {
   
   int ret;

   ret = starpu_task_insert(&cl_factorize_block,                           
                            STARPU_RW,	     l_handle,
                            STARPU_PRIORITY, prio,
                            0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;
}

// solve block StarPU task insert

void spllt_starpu_solve_block_cpu_func(void *buffers[], void *cl_arg);

// solve block task codelet
struct starpu_codelet cl_solve_block = {
  .where = STARPU_CPU,
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

  return;
}

// update block StarPU task insert

void spllt_starpu_update_block_cpu_func(void *buffers[], void *cl_arg);

void spllt_starpu_codelet_unpack_args_update_block(void *cl_arg,
                                                   int *diag) {
   
   starpu_codelet_unpack_args(cl_arg,
                              diag);

   return;
}

// update block task codelet
struct starpu_codelet cl_update_block = {
  .where = STARPU_CPU,
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
                            STARPU_RW,	     lij_handle,
                            STARPU_PRIORITY, prio,
                            0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;
}

// update between StarPU task insert

void spllt_starpu_update_between_cpu_func(void *buffers[], void *cl_arg);

void spllt_starpu_codelet_unpack_args_update_between(void *cl_arg,
                                                     void *snode, int *scol, 
                                                     void *anode, void *a_blk, int *dcol,
                                                     int *csrc_1, int *csrc_2,
                                                     int *rsrc_1, int *rsrc_2,
                                                     int *min_with_blas,
                                                     int *nhlik, int *nhljk) {
   
   starpu_codelet_unpack_args(cl_arg,
                              snode, scol, anode, a_blk, dcol,
                              csrc_1, csrc_2, 
                              rsrc_1, rsrc_2,
                              min_with_blas,
                              nhlik, nhljk);

   return;
}

// update block task codelet
struct starpu_codelet cl_update_between = {
  .where = STARPU_CPU,
  .cpu_funcs = {spllt_starpu_update_between_cpu_func, NULL},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .name = "UPDATE_BETWEEN"
};

void spllt_starpu_insert_update_between_c(starpu_data_handle_t *lik_handles, int nhlik,
                                          starpu_data_handle_t *ljk_handles, int nhljk,
                                          starpu_data_handle_t lij_handle,
                                          void *snode, int scol, 
                                          void *anode, void *a_blk, int dcol, 
                                          int csrc_1, int csrc_2, int rsrc_1, int rsrc_2,
                                          int min_width_blas,
                                          starpu_data_handle_t workspace_handle,
                                          int prio) {

   int ret, i, nh;
   struct starpu_data_descr *descrs;

   nh = 0;
   descrs = malloc((nhlik+nhljk+2) * sizeof(struct starpu_data_descr));

   descrs[nh].handle =  workspace_handle; descrs[nh].mode = STARPU_SCRATCH;
   nh = nh + 1;   

   descrs[nh].handle =  lij_handle; descrs[nh].mode = STARPU_RW;
   nh = nh + 1;

   for(i=0; i<nhlik; i++){
      descrs[i+nh].handle = lik_handles[i];  descrs[i+nh].mode = STARPU_R;
      nh = nh + 1;
   }  
   
   for(i=0; i<nhljk; i++){
      descrs[i+nh].handle = ljk_handles[i];  descrs[i+nh].mode = STARPU_R;
      nh = nh + 1;
   }

   ret = starpu_task_insert(&cl_update_between,
                            STARPU_VALUE, &snode, sizeof(void *),
                            STARPU_VALUE, &scol, sizeof(int),
                            STARPU_VALUE, &anode, sizeof(void *),
                            STARPU_VALUE, &a_blk, sizeof(void *),
                            STARPU_VALUE, &dcol, sizeof(int),
                            STARPU_VALUE, &csrc_1, sizeof(int),
                            STARPU_VALUE, &csrc_2, sizeof(int),
                            STARPU_VALUE, &rsrc_1, sizeof(int),
                            STARPU_VALUE, &rsrc_2, sizeof(int),
                            STARPU_VALUE, &min_width_blas, sizeof(int),
                            STARPU_VALUE, &nhlik, sizeof(int),
                            STARPU_VALUE, &nhljk, sizeof(int),
                            STARPU_DATA_MODE_ARRAY, descrs,   nh,
                            STARPU_PRIORITY, prio,
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


void spllt_insert_partition_task_c(starpu_data_handle_t *bc_handle, 
                                   starpu_data_handle_t **in_bc_handles, int nh,
                                   int prio) {
   int ret;
   int i;
   struct starpu_data_descr *descrs;

   descrs[0].handle = *bc_handle;        descrs[0].mode = STARPU_R;

   for(i=0; i<nh; i++){
      descrs[i+1].handle = *in_bc_handles[i];  descrs[i+1].mode = STARPU_W /* STARPU_RW */;
   }

   ret = starpu_task_insert(&cl_partition,
                            STARPU_DATA_MODE_ARRAY, descrs,   nh+1,
                            STARPU_PRIORITY,        prio,
                            0);
   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

   return;
} 

void spllt_insert_unpartition_hierarchical_task_c(starpu_data_handle_t *bc_handle, 
                                                  starpu_data_handle_t **in_bc_handles, int nh,
                                                  int prio) {

   int ret, i;
   struct starpu_data_descr *descrs;

   descrs = malloc((nh+1) * sizeof(struct starpu_data_descr));
   
   descrs[0].handle = *bc_handle;        descrs[0].mode = STARPU_W;
   /* printf("nh: %d\n", nh); */
   for(i=0; i<nh; i++){
      descrs[i+1].handle = *in_bc_handles[i];  descrs[i+1].mode = STARPU_R;
   }

   ret = starpu_task_insert(&cl_partition/* &cl_unpartition_hierarchical */,
                            STARPU_DATA_MODE_ARRAY, descrs,   nh+1,
                            STARPU_PRIORITY,        prio,
                            0);
   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

   return;
}
