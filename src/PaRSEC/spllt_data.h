!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
#ifndef __SPLLT_DATA_H__
#define __SPLLT_DATA_H__

typedef struct base_desc {
   parsec_data_collection_t super;
   parsec_data_t** data_map;   /**< map of the data */
   void *bcs;
   int nbc;
   int typesize;
} base_desc_t;

void *get_blk_ptr(void *bcs, long int id);
size_t get_blk_sze(void *bcs, long int id);

void data_init(base_desc_t *desc, void *bcs, int nbc,
               int nodes, int myrank);
#endif
