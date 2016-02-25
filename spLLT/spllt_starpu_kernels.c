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
                            STARPU_RW,	      l_handle,
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

void spllt_starpu_insert_solve_block_c(starpu_data_handle_t lik_handle,
                                           starpu_data_handle_t lkk_handle,
                                           int prio) {
   
   int ret;

   ret = starpu_task_insert(&cl_solve_block,                           
                            STARPU_RW,	      lik_handle,
                            STARPU_R,	      lkk_handle,
                            STARPU_PRIORITY, prio,
                            0);

  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_insert");

  return;
}
