#include "parsec.h"
#include "parsec_internal.h"

parsec_context_t* spllt_parsec_init( int nb_cores, int *nodes, int *rank) {

   parsec_context_t* context;

#if defined(PARSEC_HAVE_MPI)
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

   context = parsec_init(nb_cores, NULL, NULL);

   return context;
}

void spllt_parsec_fini(parsec_context_t **ctx) {

   printf("[parsec fini]\n");
   
   parsec_fini(ctx);

#ifdef PARSEC_HAVE_MPI
   MPI_Finalize();
#endif

}
