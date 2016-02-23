#include <starpu.h>


void spllt_starpu_unpack_args_factorize_block(void *cl_arg,
                                              int *m, int *n) {
   
   starpu_codelet_unpack_args(cl_arg,
                              m, n);

  return;
}
