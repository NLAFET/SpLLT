!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
#include "parsec.h"
#include "parsec_internal.h"
#include "parsec/constants.h"

#include "spllt_parsec_data.h"

static inline parsec_data_t* get_or_create_data(spllt_parsec_data_t* dat, int snode) {

   parsec_data_t* data = dat->data;

   /* if (data == NULL) { */
      
   /*    dague_data_copy_t* data_copy = OBJ_NEW(dague_data_copy_t); */
   /*    data = OBJ_NEW(dague_data_t)/\* (dague_data_t *) malloc(sizeof(dague_data_t)); *\/; */
            
   /*    data_copy->coherency_state = DATA_COHERENCY_OWNED; */
   /*    data_copy->original = NULL /\* data *\/; */
   /*    data_copy->device_private = NULL; */

   /*    data->owner_device = 0; */
   /*    data->key = 0 /\* node *\/; */
   /*    data->nb_elts = 0; */
   /*    dague_data_copy_attach(data, data_copy, 0); */

   /* } */
    
   return parsec_data_create(&data, &(dat->super), 
                             0, NULL, 0);
}

static parsec_data_t* factorize_data_of(parsec_data_collection_t *desc, ...) {

  spllt_parsec_data_t* dat = (spllt_parsec_data_t*)desc;
  va_list ap;
  int snode;

  va_start(ap, desc);
  snode = va_arg(ap, int);
  va_end(ap);
    
  return get_or_create_data(dat, snode);
}

spllt_parsec_data_t* spllt_parsec_data_create(int num_nodes) {

   spllt_parsec_data_t *data = malloc(sizeof(spllt_parsec_data_t));

   data->num_nodes = num_nodes;
   data->data = NULL;
   memcpy(&(data->super), &parsec_static_local_data_ddesc, sizeof(parsec_data_collection_t));
   data->super.data_of = factorize_data_of;

   return data;
}  
