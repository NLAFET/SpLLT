!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
#ifndef __SPLLT_PARSEC_DATA_H__
#define __SPLLT_PARSEC_DATA_H__

#include "parsec.h"
#include "parsec_internal.h"
#include "data_distribution.h"

#include "data_dist/matrix/matrix.h"
#include "data_dist/matrix/two_dim_tabular.h"

typedef struct spllt_parsec_data {
   parsec_data_collection_t super;
   int num_nodes;
   dague_data_t*  data;
} spllt_parsec_data_t;

#endif
