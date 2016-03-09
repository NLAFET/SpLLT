module spllt_starpu_factorization_mod
  implicit none

  ! interface
  !    subroutine spllt_starpu_unpack_args_factorize_block(cl_arg, m, n) bind(C)
  !      use iso_c_binding
  !      type(c_ptr), value :: cl_arg, m, n
  !    end subroutine spllt_starpu_unpack_args_factorize_block
  ! end interface

  ! factorize block StarPU task insert C
  interface
     subroutine spllt_starpu_insert_factorize_block_c(l_hdl, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value         :: l_hdl
       integer(c_int), value      :: prio
     end subroutine spllt_starpu_insert_factorize_block_c
  end interface

  ! solve block StarPU task insert C
  interface
     subroutine spllt_starpu_insert_solve_block_c(lkk_hdl, lik_hdl, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value         :: lik_hdl, lkk_hdl
       integer(c_int), value      :: prio
     end subroutine spllt_starpu_insert_solve_block_c
  end interface

  ! update block StarPU task insert C
  interface
     subroutine spllt_starpu_insert_update_block_c(lik_hdl, ljk_hdl, lij_hdl, &
          & diag, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value         :: lik_hdl, ljk_hdl, lij_hdl
       integer(c_int), value      :: diag, prio
     end subroutine spllt_starpu_insert_update_block_c
  end interface

  interface
     subroutine spllt_starpu_codelet_unpack_args_update_block(cl_arg, diag) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg, diag
     end subroutine spllt_starpu_codelet_unpack_args_update_block
  end interface

  ! update between StarPU task insert C
  interface
     subroutine spllt_starpu_insert_update_between_c(&
          & lik_hdl, nhlik, &
          & ljk_hdl, nhljk, &
          & lij_hdl, &
          & snode, scol, &
          & anode, a_bc, dcol, &
          & csrc, csrc2, rsrc, rsrc2, &
          & min_width_blas, &
          & workspace_hdl, &
          & prio) bind(C)
       use iso_c_binding
       type(c_ptr)            :: lik_hdl(*), ljk_hdl(*)
       type(c_ptr), value     :: snode, anode, a_bc
       type(c_ptr), value     :: lij_hdl, workspace_hdl
       integer(c_int), value  :: csrc, csrc2, rsrc, rsrc2
       integer(c_int), value  :: scol, dcol 
       integer(c_int), value  :: nhlik, nhljk
       integer(c_int), value  :: min_width_blas
       integer(c_int), value  :: prio
     end subroutine spllt_starpu_insert_update_between_c
  end interface

  interface
     subroutine spllt_starpu_codelet_unpack_args_update_between(cl_arg, &
          & snode, scol, a_node, a_blk, dcol, &
          & csrc, csrc2, rsrc, rsrc2, &
          & min_with_blas, &
          & nhlik, nhljk) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg 
       type(c_ptr), value :: snode, scol, a_node, a_blk, dcol
       type(c_ptr), value :: csrc, csrc2, rsrc, rsrc2
       type(c_ptr), value :: min_with_blas, nhlik, nhljk
     end subroutine spllt_starpu_codelet_unpack_args_update_between
  end interface

  ! data partitioning and unpartitioning task
  interface 
     subroutine spllt_insert_partition_task_c(hdl, in_hdls, nh, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value     :: hdl
       type(c_ptr)            :: in_hdls(*)
       integer(c_int), value  :: nh, prio
     end subroutine spllt_insert_partition_task_c

     subroutine spllt_insert_unpartition_task_c(hdl, in_hdls, nh, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value     :: hdl
       type(c_ptr)            :: in_hdls(*)
       integer(c_int), value  :: nh, prio
     end subroutine spllt_insert_unpartition_task_c
  end interface

contains

  ! factorize block StarPU task insert
  subroutine spllt_starpu_insert_factorize_block(bc, prio)
    use spllt_mod
    implicit none

    type(spllt_bc_type), intent(in) :: bc ! block to be factorized    
    integer, intent(in) :: prio

    call spllt_starpu_insert_factorize_block_c(bc%hdl, prio)
    
    return
  end subroutine spllt_starpu_insert_factorize_block

  ! factorize block StarPU task 
  ! _potrf
  subroutine spllt_starpu_factorize_block_cpu_func(buffers, cl_arg) bind(C)
    use spllt_mod
    use iso_c_binding
    use starpu_f_mod
    use spllt_kernels_mod
    implicit none

    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers
    
    integer, target          :: m, n, ld
    type(c_ptr), target      :: l
    real(wp), pointer        :: lkk(:)

    call starpu_f_get_buffer(buffers, 0, c_loc(l), c_loc(m), c_loc(n), c_loc(ld))
    ! call spllt_starpu_get_buffer(buffers, 0, c_loc(l), c_loc(ml), c_loc(nl), c_loc(ldl))
    call c_f_pointer(l, lkk,(/ld*n/))

    call spllt_factor_diag_block(m, n, lkk)

    ! factorize_block
    ! call factor_diag_block(n, m, id, &
    !      & lfact(bcol)%lcol(sa:sa+n*m-1),   &
    !      & control, potrf_info, detlog(0))

    return
  end subroutine spllt_starpu_factorize_block_cpu_func    

  ! solve block StarPU task insert
  subroutine spllt_starpu_insert_solve_block(bc_kk, bc_ik, prio)
    use spllt_mod
    implicit none
    
    ! block to be solved (bc_ik) wrt diag block (bc_kk)
    type(spllt_bc_type), intent(in) :: bc_kk, bc_ik 
    integer, intent(in) :: prio

    call spllt_starpu_insert_solve_block_c(bc_kk%hdl, bc_ik%hdl, prio)
    
    return
  end subroutine spllt_starpu_insert_solve_block

  ! solve block StarPU task 
  ! _potrf
  subroutine spllt_starpu_solve_block_cpu_func(buffers, cl_arg) bind(C)
    use spllt_mod
    use iso_c_binding
    use starpu_f_mod
    use spllt_kernels_mod
    implicit none

    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers

    type(c_ptr), target      :: a, b
    integer, target          :: m, n, ld, ma, na, lda
    real(wp), pointer        :: lkk(:), lik(:)

    call starpu_f_get_buffer(buffers, 0, c_loc(a), c_loc(ma), c_loc(na), c_loc(lda))
    call c_f_pointer(a, lkk,(/lda*na/))

    call starpu_f_get_buffer(buffers, 1, c_loc(b), c_loc(m), c_loc(n), c_loc(ld))
    call c_f_pointer(b, lik,(/ld*n/))
    
    call spllt_solve_block(m, n, lik, lkk)

    return
  end subroutine spllt_starpu_solve_block_cpu_func

  ! update block StarPU task 
  ! _syrk/_gemm
  subroutine spllt_starpu_update_block_cpu_func(buffers, cl_arg) bind(C)
    use spllt_mod
    use iso_c_binding
    use starpu_f_mod
    use spllt_kernels_mod
    implicit none

    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers
    
    type(c_ptr), target      :: dest, src1, src2
    integer, target          :: m, n, ld, m2, n2, ld2, m1, n1, ld1
    real(wp), pointer        :: lik(:), ljk(:), lij(:)
    integer, target          :: d
    logical                  :: diag

    call starpu_f_get_buffer(buffers, 0, c_loc(src2), c_loc(m2), c_loc(n2), c_loc(ld2))
    call c_f_pointer(src2, lik, (/ld2*n2/))

    call starpu_f_get_buffer(buffers, 1, c_loc(src1), c_loc(m1), c_loc(n1), c_loc(ld1))
    call c_f_pointer(src1, ljk, (/ld1*n1/))

    call starpu_f_get_buffer(buffers, 2, c_loc(dest), c_loc(m), c_loc(n), c_loc(ld))
    call c_f_pointer(dest, lij, (/ld*n/))

    call spllt_starpu_codelet_unpack_args_update_block(cl_arg, &
         & c_loc(d))
    
    if (d .eq. 1) then
       diag = .true.
    else
       diag = .false.
    end if

    call spllt_update_block(m, n, lij, diag, n1, ljk, lik)

    return
  end subroutine spllt_starpu_update_block_cpu_func

  ! update between StarPU task 
  ! _syrk/_gemm
  subroutine spllt_starpu_update_between_cpu_func(buffers, cl_arg) bind(C)
    use hsl_ma87_double
    use spllt_mod
    use spllt_kernels_mod
    use iso_c_binding
    implicit none
    
    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers

    type(c_ptr), target       :: snode_c, anode_c, a_bc_c
    integer, target      :: csrc, csrc2, rsrc, rsrc2 
    integer, target      :: scol, dcol
    integer, target      :: min_width_blas, nhlik, nhljk

    type(node_type), pointer :: snode, anode
    type(block_type), pointer :: a_bc
    integer, dimension(:), allocatable  :: row_list, col_list
    integer :: i, nh, cld2, cld1

    type(c_ptr), target      :: work_c, dest, src1, src2, ptr1, ptr2
    integer, target :: mw, nw, ldw, m, n, ld, m2, n2, ld2 
    integer, target :: m1, n1, ld1
    real(wp), pointer        :: work(:)
    ! real(wp), allocatable    :: work(:) ! TODO use StarPU scratch buffer
    real(wp), pointer        :: lij(:), lik(:), ljk(:)

    ! call spllt_starpu_codelet_unpack_args_update_between(cl_arg, &
    !      & c_loc(snode_c), c_loc(scol), &
    !      & c_loc(anode_c), c_loc(a_bc_c), c_loc(dcol), &
    !      & c_loc(csrc), c_loc(csrc2), &
    !      & c_loc(rsrc), c_loc(rsrc2), &
    !      & c_loc(min_width_blas), &
    !      & c_loc(nhlik), c_loc(nhljk))

    ! call c_f_pointer(snode_c, snode)
    ! call c_f_pointer(anode_c, anode)
    ! call c_f_pointer(a_bc_c, a_bc)

    ! nh = 0
    ! call starpu_f_get_buffer(buffers, nh, c_loc(work_c), c_loc(mw))
    ! ! TODO use scratch buffer
    ! ! call c_f_pointer(work_c, work,(/mw/))
    ! nh = nh + 1

    ! call starpu_f_get_buffer(buffers, nh, c_loc(dest), c_loc(m), c_loc(n), c_loc(ld))
    ! call c_f_pointer(dest, lij,(/ld*n/))
    ! nh = nh + 1

    ! cld2 = 0
    ! do i=1,nhlik       
    !    if (i.eq.1) then
    !       ! get pointer on first block in lik
    !       call starpu_f_get_buffer(buffers, nh, c_loc(src2), c_loc(m2), c_loc(n2), c_loc(ld2))
    !    else
    !       call starpu_f_get_buffer(buffers, nh, c_loc(ptr2), c_loc(m2), c_loc(n2), c_loc(ld2))
    !    end if
    !    cld2 = cld2 + ld2
    !    nh = nh + 1
    ! end do
    ! call c_f_pointer(src2, lik,(/cld2*n2/))    

    ! cld1 = 0
    ! do i=1,nhlik
    !    if (i.eq.1) then
    !       ! get pointer on first block in ljk
    !       call starpu_f_get_buffer(buffers, nh, c_loc(src1), c_loc(m1), c_loc(n1), c_loc(ld1))
    !    else
    !       call starpu_f_get_buffer(buffers, nh, c_loc(ptr1), c_loc(m1), c_loc(n1), c_loc(ld1))
    !    end if
    !    cld1 = cld1 + ld1
    !    nh = nh + 1
    ! end do
    ! call c_f_pointer(src1, ljk,(/cld1*n1/))    
    
    ! TODO use scratch memory
    ! allocate(col_list(1), row_list(1))
    ! allocate(work(m*n))

    ! call spllt_update_between(m, n, a_bc, dcol, anode, &
    !      & n1, scol, snode, &
    !      & lij, &
    !      & ljk(csrc:csrc+csrc2-1), &
    !      & lik(rsrc:rsrc+rsrc2-1), &
    !      & row_list, col_list, work, &
    !      & min_width_blas)

    ! deallocate(work)

  end subroutine spllt_starpu_update_between_cpu_func

end module spllt_starpu_factorization_mod
