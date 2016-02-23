#include <starpu.h>
#include <starpu_cuda.h>
#include <starpu_profiling.h>
#include <limits.h>

int starpu_f_init_c(int ncpus) {

  int info;
  struct starpu_conf *conf = malloc(sizeof(struct starpu_conf));
  
  starpu_conf_init(conf);

  if(ncpus > 0)
    conf->ncpus = ncpus;  

  conf->sched_policy_name = "lws";

  info = starpu_init(conf);

  free(conf);
  return info;
}

void starpu_f_get_buffer(void *buffers[], int num, void **A, int *m, int *n, int *lda) {

  *A   = (void *)STARPU_MATRIX_GET_PTR(buffers[num]);
  *m   = (int)STARPU_MATRIX_GET_NX(buffers[num]);
  *n   = (int)STARPU_MATRIX_GET_NY(buffers[num]);
  *lda = (int)STARPU_MATRIX_GET_LD(buffers[num]);

  return;

}
