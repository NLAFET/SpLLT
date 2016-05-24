/* Codelet definition for kernel f */
struct starpu_codelet f_cl =
  {
    .where = STARPU_CPU,
    .cpu_funcs = { f_cpu_func },
    .nbuffers = 1
  };

/* Codelet definition for kernel g */
struct starpu_codelet g_cl =
  {
    .where = STARPU_CPU | STARPU_CUDA,
    .cpu_funcs = { g_cpu_func },
    .cuda_funcs = { g_cuda_func },
    .nbuffers = 3
  };

starpu_data_handle_t x_handle[N], y_handle[N];

starpu_init();

/* declaration of data handles */
for (i = 0; i < N; i++) {  >\label{code:ex1}>
  starpu_variable_data_register(&x_handle[i], STARPU_MAIN_RAM, 
                                (uintptr_t) &x[i], sizeof(double));
  starpu_variable_data_register(&y_handle[i], STARPU_MAIN_RAM, 
                                (uintptr_t) &y[i], sizeof(double));
 }

/* tasks submission */
for (i = 1; i < N; i++) {  >\label{code:ex2}>

  starpu_insert_task(&f_cl, 
                     STARPU_RW, x_handle[i],
                     0);

  starpu_insert_task(&g_cl, 
                     STARPU_R, x_handle[i],
                     STARPU_R, y_handle[i-1],
                     STARPU_W, y_handle[i],
                     0);
}

starpu_task_wait_for_all();  >\label{code:ex3}>

starpu_shutdown();
