module spllt_starpu_factorization_mod
  implicit none

  interface
     subroutine spllt_starpu_unpack_args_factorize_block(cl_arg, m, n) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg, m, n
     end subroutine spllt_starpu_unpack_args_factorize_block
  end interface

contains
  
  ! factorize block StarPU task 
  ! _potrf
  subroutine spllt_starpu_factorize_block_cpu_func(buffers, cl_arg) bind(C)
    use iso_c_binding
    use starpu_f_mod
    implicit none

    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers
    
    integer :: potrf_info
    integer :: wrk_id 
    integer, target :: m, n

    wrk_id = starpu_f_worker_get_id()

    call spllt_starpu_unpack_args_factorize_block(cl_arg, & 
         & c_loc(m), c_loc(n))

    ! factorize_block
    ! call factor_diag_block(n, m, id, &
    !      & lfact(bcol)%lcol(sa:sa+n*m-1),   &
    !      & control, potrf_info, detlog(0))

    return
  end subroutine spllt_starpu_factorize_block_cpu_func    

end module spllt_starpu_factorization_mod
