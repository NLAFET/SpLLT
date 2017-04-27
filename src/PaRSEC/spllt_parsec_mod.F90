module spllt_parsec_mod
  use iso_c_binding
  use parsec_f08_interfaces
  
  ! C binding funtions
  ! general Parsec function
  interface
     ! initialize Parsec
     function spllt_parsec_init(nb_cores, nodes, rank) bind(C)
       use iso_c_binding
       use parsec_f08_interfaces
       integer(c_int), value   :: nb_cores ! number of threads
       integer(c_int)          :: nodes ! number of nodes involved in the execution
       integer(c_int)          :: rank ! rank of the current process
       type(parsec_context_t)   :: spllt_parsec_init ! Parsec context i.e. context of execution 
     end function spllt_parsec_init
     ! finilize Parsec
     subroutine spllt_parsec_fini(ctx) bind(C)
       use iso_c_binding
       use parsec_f08_interfaces
       type(parsec_context_t)   :: ctx ! Parsec context i.e. context of execution 
     end subroutine spllt_parsec_fini
  end interface

  ! C binding funtion
  ! block data descriptor
  interface
     ! allocate block data descriptor
     function spllt_alloc_blk_desc() bind(C)
       use iso_c_binding
       type(c_ptr) :: spllt_alloc_blk_desc
     end function spllt_alloc_blk_desc
     ! initialize block data descriptor
     subroutine spllt_parsec_blk_data_init(blk_desc, bcs, nbc, nodes, myrank) bind(C)
       use iso_c_binding
       type(c_ptr), value      :: blk_desc ! data descriptor
       type(c_ptr), value      :: bcs ! blocks
       integer(c_int), value   :: nbc ! number of blocks
       integer(c_int), value   :: nodes ! number of nodes involved in the execution
       integer(c_int), value   :: myrank ! rank of the current process
     end subroutine spllt_parsec_blk_data_init
  end interface

  ! C binding funtion
  

end module spllt_parsec_mod
