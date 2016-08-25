#include <dague.h>
#include <dague/data.h>
#include <dague/data_distribution.h>

#include "spllt_parsec_blk_data.h"

static inline uint32_t blk_key(int id) {

   uint32_t key;

   key = (uint32_t) id-1;

   return  key;
}

static uint32_t blk_data_key(dague_ddesc_t *desc, ... ) {

    va_list ap;
    /* sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat; */
    int id;

    va_start(ap, desc);
    id = va_arg(ap, int);
    va_end(ap);

    return blk_key(id);
}

static uint32_t blk_rank_of(dague_ddesc_t *desc, ... ) {
    (void)desc;
    return 0;
}

static uint32_t blk_rank_of_key(dague_ddesc_t *desc, dague_data_key_t key) {
(void)desc; (void)key;
return 0;
}

static int32_t blk_vpid_of(dague_ddesc_t *desc, ... ) {
    (void)desc;
    return 0;
}

static int32_t blk_vpid_of_key(dague_ddesc_t *desc, dague_data_key_t key) {
    (void)desc; (void)key;
    return 0;
}

static dague_data_t *blk_data_of(dague_ddesc_t *desc, ... ) {

    blk_desc_t *blk_desc = (blk_desc_t*)desc;

    va_list ap;
    int id;
    dague_data_key_t key;
    int pos;
    size_t size;
    void *bcs; // bloc list
    void *bc;  // bloc ptr

    va_start(ap, desc);
    id = va_arg(ap, int);
    va_end(ap);

    /* nbcol = blk_desc->nbcol; */
    bcs = blk_desc->bcs;
    bc  = (void *) get_blk_ptr(bcs, id);
    /* bc  = NULL; */

    key  = blk_key(id);
    pos  = key;
    size = get_blk_sze(bcs, id);

    /* printf("[blk_data_of] id: %d, key: %zu\n", id, key); */
    /* printf("[blk_data_of] key: %d, bc: %p, size: %zu\n", key, bc, size); */
    /* size = sizeof(double); */

    return dague_data_create(blk_desc->data_map + pos, desc, key, bc, size);
}

void spllt_parsec_blk_data_init(blk_desc_t *desc,
                                void *bcs, int nbc,
                                int nodes, int myrank) {
   
    dague_ddesc_t *o = (dague_ddesc_t*)desc;

    /* printf("[spllt_parsec_blk_data_init] myrank: %d\n", myrank); */
    
    /* Super setup */
    o->nodes     = nodes;
    o->myrank    = myrank;

    o->data_key      = blk_data_key;
/* #if defined(DAGUE_PROF_TRACE) */
/*     o->key_to_string = sparse_matrix_key_to_string; */
/* #endif */

    o->rank_of     = blk_rank_of;
    o->rank_of_key = blk_rank_of_key;
    o->vpid_of     = blk_vpid_of;
    o->vpid_of_key = blk_vpid_of_key;
    o->data_of     = blk_data_of;
    /* o->data_of_key = blk_data_of_key; */

    /* desc->typesize  = typesize; */
    desc->bcs       = bcs;
    desc->nbc       = nbc;
    desc->data_map  = (dague_data_t**)calloc( nbc, sizeof(dague_data_t*) );
}

blk_desc_t *spllt_alloc_blk_desc() {

   blk_desc_t *desc = malloc(sizeof(blk_desc_t));
                             
   return desc;
}