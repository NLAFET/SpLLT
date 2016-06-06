module spllt_parsec_mod
  use dague_f08_interfaces

  interface
     function parsec_init(nb_cores, nodes, rank) bind(C)
       use iso_c_binding
       use dague_f08_interfaces
       integer(c_int), value   :: nb_cores
       integer(c_int)          :: nodes
       integer(c_int)          :: rank
       type(dague_context_t)   :: parsec_init
     end function parsec_init
  end interface

end module spllt_parsec_mod
