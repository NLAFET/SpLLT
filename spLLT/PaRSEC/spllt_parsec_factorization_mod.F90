module spllt_parsec_factorization_mod

  interface
     function spllt_parsec_factorize(num_nodes) bind(C)
       use iso_c_binding
       use dague_f08_interfaces
       integer(c_int), value   :: num_nodes
       type(dague_handle_t)    :: spllt_parsec_factorize
     end function spllt_parsec_factorize
  end interface

end module spllt_parsec_factorization_mod
