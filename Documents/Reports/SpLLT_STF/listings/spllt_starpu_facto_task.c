struct starpu_codelet cl_factorize_block = {
  .where = STARPU_CPU,
  .cpu_funcs = {spllt_factorize_block, NULL},
  .nbuffers = STARPU_VARIABLE_NBUFFERS,
  .name = "FACTO_BLK"
};

starpu_task_insert(&cl_factorize_block,                           
                   STARPU_RW,	bc_handle,
                   STARPU_R,	node_handle,
                   STARPU_PRIORITY, prio,
                   0);
