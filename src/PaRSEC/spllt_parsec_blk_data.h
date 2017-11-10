#ifndef __SPLLT_PARSEC_BLK_DATA_H__
#define __SPLLT_PARSEC_BLK_DATA_H__

typedef struct blk_desc {
   parsec_data_collection_t super;
   parsec_data_t** data_map;   /**< map of the data */
   void *bcs;
   int nbc;
   int typesize;
} blk_desc_t;

void *get_blk_ptr(void *bcs, long int id);
size_t get_blk_sze(void *bcs, long int id);

void spllt_parsec_blk_data_init(blk_desc_t *desc, void *bcs, int nbc,
                                int nodes, int myrank);
#endif
