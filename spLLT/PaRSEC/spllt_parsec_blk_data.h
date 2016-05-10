#ifndef __SPLLT_PARSEC_BLK_DATA_H__
#define __SPLLT_PARSEC_BLK_DATA_H__

typedef struct blk_desc {
   dague_ddesc_t super;
   dague_data_t** data_map;   /**< map of the data */
   void *bcs;
   int nbc;
   int typesize;
} blk_desc_t;

#endif
