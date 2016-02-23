#include <starpu.h>


void _qrm_starpu_unpack_args_geqrt(void *cl_arg,
                                   int *m, int *n, int id) {
   
   starpu_codelet_unpack_args(cl_arg,
                              m, n, id);

  return;
}
