#include "dague.h"
#include "dague_internal.h"

dague_context_t* parsec_init( int nb_cores, int *nodes, int *rank) {

   dague_context_t* context;

#if defined(DAGUE_HAVE_MPI)
    {
        int provided;
        MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &provided);
    }
    MPI_Comm_size(MPI_COMM_WORLD, nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
#else
    *nodes = 1;
    *rank = 0;
#endif

   context = dague_init(nb_cores, NULL, NULL);

   return context;
}
