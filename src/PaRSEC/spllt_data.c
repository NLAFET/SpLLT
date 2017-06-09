#include <parsec.h>
#include <parsec/data.h>
#include <parsec/data_distribution.h>

#include "spllt_data.h"

static inline uint32_t get_key(int id) {

   uint32_t key;

   key = (uint32_t) id-1;

   return  key;
}

static uint32_t data_key(parsec_ddesc_t *desc, ... ) {

    va_list ap;
    /* sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat; */
    int id;

    va_start(ap, desc);
    id = va_arg(ap, int);
    va_end(ap);

    return get_key(id);
}

static uint32_t rank_of(parsec_ddesc_t *desc, ... ) {
    (void)desc;
    return 0;
}

static uint32_t rank_of_key(parsec_ddesc_t *desc, parsec_data_key_t key) {

   (void)desc; (void)key;
   return 0;
}

static int32_t vpid_of(parsec_ddesc_t *desc, ... ) {
   
   /* printf("[blk_vpid_of] TETETETETETETETE\n"); */
   /* int nvp = vpmap_get_nb_vp(); */
   /* printf("[blk_vpid_of] nvp: %d\n", nvp); */
   
   (void)desc;
   return 0;
}

static int32_t vpid_of_key(parsec_ddesc_t *desc, parsec_data_key_t key) {

   /* printf("[blk_vpid_of_key] TETETETETETETETE\n"); */

    (void)desc; (void)key;
    return 0;
}

static parsec_data_t *data_of(parsec_ddesc_t *desc, ... ) {

    base_desc_t *base_desc = (base_desc_t*)desc;

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
    bcs = base_desc->bcs;
    bc  = (void *) get_blk_ptr(bcs, id);
    /* bc  = NULL; */

    key  = get_key(id);
    /* key  = 0; */
    pos  = key;
    size = get_blk_sze(bcs, id);
    /* size = 1; */

    /* printf("[data_of] id: %d, key: %zu\n", id, key); */
    printf("[data_of] key: %d, bc: %p, size: %zu\n", key, bc, size);
    /* size = sizeof(double); */

    return parsec_data_create(base_desc->data_map + pos, desc, key, bc, size);
}

static parsec_data_t *data_of_key(parsec_ddesc_t *desc, parsec_data_key_t key) {
   /* printf("[data_of_key] TETETETETETE\n"); */
   (void)desc; (void)key;
   return 0;   
}


void data_init(base_desc_t *desc,
               void *bcs, int nbc,
               int nodes, int myrank) {
   
    parsec_ddesc_t *o = (parsec_ddesc_t*)desc;

    /* printf("[spllt_parsec_blk_data_init] myrank: %d\n", myrank); */
    
    /* Super setup */
    parsec_ddesc_init(o, nodes, myrank);

    o->data_key      = data_key;
/* #if defined(DAGUE_PROF_TRACE) */
/*     o->key_to_string = key_to_string; */
/* #endif */

    o->rank_of     = rank_of;
    o->rank_of_key = rank_of_key;
    o->vpid_of     = vpid_of;
    o->vpid_of_key = vpid_of_key;
    o->data_of     = data_of;
    o->data_of_key = data_of_key;

    /* desc->typesize  = sizeof(double); */
    desc->bcs       = bcs;
    desc->nbc       = nbc;
    desc->data_map  = (parsec_data_t**) calloc(nbc, sizeof(parsec_data_t*));
}

base_desc_t *alloc_desc() {

   base_desc_t *desc = malloc(sizeof(base_desc_t));
                             
   return desc;
}
