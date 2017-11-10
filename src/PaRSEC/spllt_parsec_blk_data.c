#include <parsec.h>
#include <parsec/data.h>
#include <parsec/data_distribution.h>

#include "spllt_parsec_blk_data.h"

static inline uint32_t blk_key(int id) {

   uint32_t key;

   key = (uint32_t) id-1;

   return  key;
}

static uint32_t blk_data_key(parsec_data_collection_t *desc, ... ) {

    va_list ap;
    /* sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat; */
    int id;

    va_start(ap, desc);
    id = va_arg(ap, int);
    va_end(ap);

    return blk_key(id);
}

static uint32_t blk_rank_of(parsec_data_collection_t *desc, ... ) {

    /* (void)desc; */
    /* return 0; */

    /* printf("[blk_rank_of]\n"); */

    int nodes = desc->nodes;

    /* printf("[blk_rank_of] nodes: %d\n", nodes); */

    va_list ap;
    int id;

    va_start(ap, desc);
    id = va_arg(ap, int);
    va_end(ap);

    blk_desc_t *blk_desc = (blk_desc_t*)desc;
    void *bcs = blk_desc->bcs;
    int bcol = get_bcol(bcs, id);
    /* printf("[blk_rank_of] bcol: %d\n", bcol); */

    /* printf("[blk_rank_of] id: %d\n", id); */
    /* int res = (bcol-1) % nodes; */
    /* int res = id % nodes; */
    int res = 1;

    /* printf("[blk_rank_of] rank_of: %d\n", res); */

    return res;
}

static uint32_t blk_rank_of_key(parsec_data_collection_t *desc, parsec_data_key_t key) {

   (void)desc; (void)key;
   return 0;

   /* int id = (int) key+1; */
   /* printf("[blk_rank_of_key] id: %d\n"); */
   /* return 0;    */
}

static int32_t blk_vpid_of(parsec_data_collection_t *desc, ... ) {
   
   /* printf("[blk_vpid_of] TETETETETETETETE\n"); */
   /* int nvp = vpmap_get_nb_vp(); */
   /* printf("[blk_vpid_of] nvp: %d\n", nvp); */
   
   (void)desc;
   return 0;
}

static int32_t blk_vpid_of_key(parsec_data_collection_t *desc, parsec_data_key_t key) {

   /* printf("[blk_vpid_of_key] TETETETETETETETE\n"); */

    (void)desc; (void)key;
    return 0;
}

static parsec_data_t *blk_data_of(parsec_data_collection_t *desc, ... ) {

    blk_desc_t *blk_desc = (blk_desc_t*)desc;

    va_list ap;
    int id;
    parsec_data_key_t key;
    int pos;
    size_t size;
    void *bcs; // bloc list
    void *bc;  // bloc ptr

    va_start(ap, desc);
    id = va_arg(ap, int);
    va_end(ap);

    /* nbcol = blk_desc->nbcol; */
    bcs = blk_desc->bcs;
    bc  = (void *)get_blk_ptr(bcs, id);
    /* bc  = NULL; */

    key  = blk_key(id);
    /* key  = 0; */
    pos  = key;
    size = get_blk_sze(bcs, id);
    /* size = 1; */

    /* printf("[blk_data_of] id: %d, key: %zu\n", id, key); */
    /* printf("[blk_data_of] key: %d, bc: %p, size: %zu\n", key, bc, size); */
    /* size = sizeof(double); */

    return parsec_data_create(blk_desc->data_map + pos, desc, key, bc, size);
}

static parsec_data_t *blk_data_of_key(parsec_data_collection_t *desc, parsec_data_key_t key) {
   /* printf("[blk_data_of_key] TETETETETETE\n"); */
    (void)desc; (void)key;
    return 0;   
}


void spllt_parsec_blk_data_init(blk_desc_t *desc,
                                void *bcs, int nbc,
                                int nodes, int myrank) {
   
    parsec_data_collection_t *o = (parsec_data_collection_t*)desc;

    /* printf("[spllt_parsec_blk_data_init] myrank: %d\n", myrank); */
    
    /* Super setup */
    parsec_data_collection_init( o, nodes, myrank );

    /* o->nodes     = nodes; */
    /* o->myrank    = myrank; */

    o->data_key      = blk_data_key;
/* #if defined(DAGUE_PROF_TRACE) */
/*     o->key_to_string = sparse_matrix_key_to_string; */
/* #endif */

    o->rank_of     = blk_rank_of;
    o->rank_of_key = blk_rank_of_key;
    o->vpid_of     = blk_vpid_of;
    o->vpid_of_key = blk_vpid_of_key;
    o->data_of     = blk_data_of;
    o->data_of_key = blk_data_of_key;

    /* desc->typesize  = sizeof(parsec_datatype_double_t); */
    desc->bcs       = bcs;
    desc->nbc       = nbc;
    desc->data_map  = (parsec_data_t**)calloc( nbc, sizeof(parsec_data_t*) );
}

blk_desc_t *spllt_alloc_blk_desc() {

   blk_desc_t *desc = malloc(sizeof(blk_desc_t));
                             
   return desc;
}
