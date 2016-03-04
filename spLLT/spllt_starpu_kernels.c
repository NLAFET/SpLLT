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
                                                     void *blocks,
                                                     void *snode, int *blk, 
                                                     void *anode, int *a_blk,
                                                     void *csrc, void *rsrc,
                                                     void *row_list, void *col_list,
                                                     int *min_with_blas,
                                                     int *nhlik, int *nhljk) {
   
   starpu_codelet_unpack_args(cl_arg,
                              blocks,
                              snode, blk, anode, a_blk,
                              csrc, rsrc, row_list, col_list,
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
                                          void *blocks,
                                          void *snode, long int blk, 
                                          void *anode, long int a_blk,
                                          int *csrc, int *rsrc,
                                          int *row_list, int *col_list,
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
                            STARPU_VALUE, &blocks, sizeof(void *),
                            STARPU_VALUE, &snode, sizeof(void *),
                            STARPU_VALUE, &blk, sizeof(long int),
                            STARPU_VALUE, &anode, sizeof(void *),
                            STARPU_VALUE, &a_blk, sizeof(long int),
                            STARPU_VALUE, &csrc, sizeof(int *),
                            STARPU_VALUE, &rsrc, sizeof(int *),
                            STARPU_VALUE, &row_list, sizeof(int *),
                            STARPU_VALUE, &col_list, sizeof(int *),
                            STARPU_VALUE, &min_width_blas, sizeof(int),
                            STARPU_VALUE, &nhlik, sizeof(int),
                            STARPU_VALUE, &nhljk, sizeof(int),
                            STARPU_DATA_MODE_ARRAY, descrs,   nh,
                            STARPU_PRIORITY, prio,
                            0);
   
   STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

   return;
}
