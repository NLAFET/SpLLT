struct starpu_codelet cl_solve_block = {
   .where = STARPU_CPU,
   .cpu_funcs = {spllt_solve_block, NULL},
};

starpu_task_insert(&cl_solve_block,                     
                   STARPU_R,	     blk_kk_handle,
                   STARPU_RW,	     blk_ik_handle,
                   STARPU_PRIORITY, prio,
                   0);
