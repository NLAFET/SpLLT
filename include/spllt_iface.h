#ifndef SPLLT_IFACE_H
#define SPLLT_IFACE_H

typedef struct{
  void *akeep;
  void *fkeep;
  void *tm;
} spllt_data_t;

typedef struct{
  int print_level;
  int nrhs;
  int ncpu;
  int nb;
//char mat[100];
  int nemin;
  int prune_tree;
//char fmt[3];
  int min_width_blas;
  int nb_min;
  int nb_max;
  int nrhs_min;
  int nrhs_max;
  int nb_linear_comp;
  int nrhs_linear_comp;
  int chunk;
} spllt_options_t;

#define SPLLT_OPTIONS_NULL() {.print_level=0,     \
                              .nrhs=1,            \
                              .ncpu=1,            \
                              .nb=16,             \
                              .nemin=32,          \
                              .prune_tree=1,      \
                              .min_width_blas=8,  \
                              .nb_min=32,         \
                              .nb_max=32,         \
                              .nrhs_min=1,        \
                              .nrhs_max=1,        \
                              .nb_linear_comp=0,  \
                              .nrhs_linear_comp=0,\
                              .chunk=10           \
                              }

typedef struct{
//ssids_inform_t ssids_inform;
  int flag;
  int maxdepth;
  int num_factor;
  int num_flops;
  int num_nodes;
  int stat;
} spllt_inform_t;

extern void spllt_analyse(void            **akeep,
                          void            **fkeep,
                          spllt_options_t *options,
                          int             n,
                          int             *ptr,
                          int             *row,
                          spllt_inform_t  *info,
                          int             *order);

extern void spllt_factor( void            *akeep,
                          void            *fkeep,
                          spllt_options_t *options,
                          int             nnz,
                          double          *val,
                          spllt_inform_t  *info);

extern void spllt_prepare_solve(void            *akeep,
                                void            *fkeep,
                                int             nb,
                                int             nrhs,
                                long            *worksize,
                                spllt_inform_t  *info);

extern void spllt_set_mem_solve(void            *akeep,
                                void            *fkeep,
                                int             nb,
                                int             nrhs,
                                long            worksize,
                                double          *y,
                                double          *workspace,
                                spllt_inform_t  *info);

extern void spllt_solve_workspace_size( void  *fkeep,
                                        int   nworker,
                                        int   nrhs,
                                        long  *size);

extern void spllt_solve(void            *fkeep,
                        spllt_options_t *options,
                        int             *order,
                        int             nrhs,
                        double          *x,
                        spllt_inform_t  *info,
                        int             job);

extern void spllt_solve_worker( void            *fkeep,
                                spllt_options_t *options,
                                int             *order,
                                int             nrhs,
                                double          *x,
                                spllt_inform_t  *info,
                                int             job,
                                double          *workspace,
                                long            worksize,
                                void            *tm);

extern void spllt_wait();

extern void spllt_chkerr( int n,
                          int *ptr,
                          int *row,
                          double *val,
                          int nrhs,
                          double *x,
                          double *rhs);

extern void spllt_deallocate_fkeep( void  **fkeep,
                                    int   *stat);
                        
extern void spllt_deallocate_akeep( void  **akeep,
                                    int   *stat);

extern void spllt_task_manager_deallocate(void  **task_manager,
                                          int   *stat);

extern void spllt_task_manager_init(void  **task_manager);

extern void spllt_all(void **akeep,
                      void **fkeep,
                      spllt_options_t *options,
                      int             n,
                      int             nnz,
                      int             nrhs,
                      int             nb,
                      int             *ptr,
                      int             *row,
                      double          *val,
                      double          *x,
                      double          *rhs,
                      spllt_inform_t  *info);
#endif
