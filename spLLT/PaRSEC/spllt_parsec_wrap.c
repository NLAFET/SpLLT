#include "dague.h"

dague_context_t* parsec_init( int nb_cores ) {

   dague_context_t* context;

#if defined(DAGUE_HAVE_MPI)
#endif

   context = dague_init(nb_cores, NULL, NULL);

   return context;
}
