#ifndef SPLLT_IFACE_H
#define SPLLT_IFACE_H

typedef struct{
  int print_level;
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
} spllt_options_t;

#define SPLLT_OPTIONS_NULL() {.print_level=0,     \
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
                              .nrhs_linear_comp=0 \
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
                                spllt_inform_t  *info);


extern void spllt_solve(void            *fkeep,
                        spllt_options_t *options,
                        int             *order,
                        int             nrhs,
                        double          *x,
                        spllt_inform_t  *info,
                        int             job);

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
#if 0
typedef struct{
  double *lcol;
} spllt_lfactor_t;


typedef struct{
  int len_map;
  int **map;
} spllt_lmap_type_t;

typedef struct{
  int num;
  int nnode;
  int node_sa;
  int node_en;
  int nchild;
  int *child;
  int parent;
} spllt_tree_t;

typedef struct{
  int nnodes;
  int n;
  int num_factor;
  int num_flops;
  int *weight;
  int *small;
} spllt_akeep_t;

typedef struct{
  spllt_block_t *bc;
  spllt_block_t *workspace;
  spllt_node_t  *nodes;
  spllt_workspace_i_t *row_list;
  spllt_workspace_i_t *col_list;
  spllt_workspace_i_t *map;
  int *flag_array;
  int final_blk;
  spllt_inform_t info;
  int maxmn;
  int n;
  int nbcol;
  spllt_lfactor_t *lfactor;
  spllt_lmap_type_t *lmap;
  spllt_tree_t *trees;
  int *small;
  int *assoc_tree;
} spllt_fkeep_t;

typedef struct{
  double **c;
  int mem_node;
  int bcol;
  int blkm;
  int blkn;
  int dblk;
  int dep_initial;
  int id;
  int last_blk;
  int node;
  int sa;
  int dep;
  int *fwd_dep;
  int *fwd_update_dep;
  int *fwd_solve_dep;
  int *bwd_update_dep;
  int *bwd_solve_dep;
  int *bwd_dep;
} spllt_block_t;

typedef struct{
  int num;
  spllt_block_t *buffer;
  int blk_sa;
  int blk_en;
  int nb;
  int sa;
  int en;
  int *index;
  int nchild;
  int *child;
  int parent;
  int least_desc;
  int *extra_row;
} spllt_node_t;

typedef struct{
  int **c;
} spllt_workspace_i_t;

  
extern void spllt_analyse(spllt_akeep_t   *akeep,
                          spllt_fkeep_t   *fkeep,
                          spllt_options_t *options,
                          int             n,
                          int             *ptr,
                          int             *row,
                          spllt_inform_t  *info,
                          int             *order);

extern void spllt_factorize(spllt_akeep_t   *akeep,
                            spllt_fkeep_t   *fkeep,
                            spllt_options_t *options,
                            double          *val,
                            spllt_inform_t  *info);

extern void spllt_solve(spllt_fkeep_t   *fkeep,
                        spllt_options_t *options,
                        int             *order,
                        int             *nrhs,
                        double          *x,
                        spllt_inform_t  *info,
                        int             *job,
                        double          *workspace);

extern void spllt_deallocate_fkeep( spllt_fkeep_t *fkeep,
                                    int           *stat);
                        
extern void spllt_deallocate_akeep( spllt_akeep_t *akeep,
                                    int           *stat);
                        
#endif
#endif
