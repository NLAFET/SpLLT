!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
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
     subroutine spllt_starpu_insert_factorize_block_c(l_hdl, node_hdl, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value         :: l_hdl, node_hdl
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
#if defined(SPLLT_USE_GPU)
     subroutine spllt_starpu_insert_update_between_c(&
          & lik_hdl, &
          & ljk_hdl, &
          & lij_hdl, &
          & snode, scol, &
          & anode, a_bc, dcol, &
          & csrc, csrc2, rsrc, rsrc2, &
          & min_width_blas, &
          & workspace_hdl, row_list_hdl, col_list_hdl, &
          & node_hdl, &
          & prio) bind(C)
       use iso_c_binding
       type(c_ptr), value     :: lik_hdl, ljk_hdl
       type(c_ptr), value     :: snode, anode, a_bc
       type(c_ptr), value     :: lij_hdl, workspace_hdl, node_hdl
       type(c_ptr), value     :: row_list_hdl, col_list_hdl 
       integer(c_int), value  :: csrc, csrc2, rsrc, rsrc2
       integer(c_int), value  :: scol, dcol
       integer(c_int), value  :: min_width_blas
       integer(c_int), value  :: prio
     end subroutine spllt_starpu_insert_update_between_c
#else
     subroutine spllt_starpu_insert_update_between_c(&
          & lik_hdl, nhlik, &
          & ljk_hdl, nhljk, &
          & lij_hdl, &
          & snode, scol, &
          & anode, a_bc, dcol, &
          & csrc, csrc2, rsrc, rsrc2, &
          & min_width_blas, &
          & workspace_hdl, row_list_hdl, col_list_hdl, &
          & node_hdl, &
          & prio) bind(C)
       use iso_c_binding
       type(c_ptr)            :: lik_hdl(*), ljk_hdl(*)
       type(c_ptr), value     :: snode, anode, a_bc
       type(c_ptr), value     :: lij_hdl, workspace_hdl, node_hdl
       type(c_ptr), value     :: row_list_hdl, col_list_hdl
       integer(c_int), value  :: csrc, csrc2, rsrc, rsrc2
       integer(c_int), value  :: scol, dcol 
       integer(c_int), value  :: nhlik, nhljk
       integer(c_int), value  :: min_width_blas
       integer(c_int), value  :: prio
     end subroutine spllt_starpu_insert_update_between_c
#endif
  end interface

  interface
#if defined(SPLLT_USE_GPU)
     subroutine spllt_starpu_codelet_unpack_args_update_between(cl_arg, &
          & snode, scol, a_node, a_blk, dcol, &
          & csrc, csrc2, rsrc, rsrc2, &
          & min_with_blas) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg 
       type(c_ptr), value :: snode, scol, a_node, a_blk, dcol
       type(c_ptr), value :: csrc, csrc2, rsrc, rsrc2
       type(c_ptr), value :: min_with_blas
     end subroutine spllt_starpu_codelet_unpack_args_update_between
#else
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
#endif
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
  
  ! tests
  interface
     subroutine test_insert_c(a_hdls, nah, b_hdls, nbh) bind(C)
       use iso_c_binding
       type(c_ptr)            :: a_hdls(*), b_hdls(*)
       integer(c_int), value  :: nah, nbh
     end subroutine test_insert_c
  end interface

  ! init node StarPU task insert
  interface
     subroutine spllt_insert_init_node_task_c(hdl, snode, val, nval, fkeep, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value     :: hdl
       integer(c_int), value  :: snode, prio, nval
       type(c_ptr), value     :: val, fkeep
     end subroutine spllt_insert_init_node_task_c
  end interface

  ! unpack args
  interface
     subroutine spllt_starpu_codelet_unpack_args_init_node(cl_arg, &
          & snode, val, nval, fkeep) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg
       type(c_ptr), value :: snode, val, fkeep, nval
     end subroutine spllt_starpu_codelet_unpack_args_init_node
  end interface

  ! subtree_factorize StarPU task insert
  interface
     subroutine spllt_insert_subtree_factorize_task_c(root, val, fkeep, &
          & buffer_hdl, cntl, map_hdl, rlst_hdl, clst_hdl, work_hdl) bind(C)
       use iso_c_binding
       integer(c_int), value  :: root ! root node of the subtree
       type(c_ptr), value     :: val, fkeep, cntl ! info
       type(c_ptr), value     :: buffer_hdl ! buffer memory containing accumulated update
       type(c_ptr), value     :: map_hdl, rlst_hdl, clst_hdl, work_hdl ! scratch memory
     end subroutine spllt_insert_subtree_factorize_task_c
  end interface

  ! unpack args
  interface
     subroutine spllt_starpu_codelet_unpack_args_subtree_factorize(cl_arg, root, &
          & val, fkeep, cntl) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg
       type(c_ptr), value :: root
       type(c_ptr), value :: val, fkeep, cntl
     end subroutine spllt_starpu_codelet_unpack_args_subtree_factorize
  end interface

  ! subtree_scatter_block StarPU task insert
  interface
     subroutine spllt_starpu_insert_subtree_scatter_block_task_c(rptr, rptr2, cptr, cptr2, &
          buffer_hdl, root, a_rptr, a_cptr, dest_hdl, anode) bind (C)
       use iso_c_binding
       integer(c_int), value  :: rptr, rptr2, cptr, cptr2
       type(c_ptr), value  :: buffer_hdl, root
       integer(c_int), value  :: a_rptr, a_cptr
       type(c_ptr), value  :: dest_hdl, anode
     end subroutine spllt_starpu_insert_subtree_scatter_block_task_c
  end interface

  ! unpack args
  interface
     subroutine spllt_starpu_codelet_unpack_args_subtree_scatter_block(cl_arg, & 
          rptr, rptr2, cptr, cptr2, root, a_rptr, a_cptr, anode) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg
       type(c_ptr), value :: rptr, rptr2, cptr, cptr2
       type(c_ptr), value :: root, anode
       type(c_ptr), value :: a_rptr, a_cptr
     end subroutine spllt_starpu_codelet_unpack_args_subtree_scatter_block
  end interface

contains

  ! factorize block StarPU task 
  ! _potrf
  subroutine spllt_starpu_factorize_block_cpu_func(buffers, cl_arg) bind(C)
    use spllt_data_mod
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
    use spllt_data_mod
    implicit none
    
    ! block to be solved (bc_ik) wrt diag block (bc_kk)
    type(spllt_block), intent(in) :: bc_kk, bc_ik 
    integer, intent(in) :: prio

    call spllt_starpu_insert_solve_block_c(bc_kk%hdl, bc_ik%hdl, prio)
    
    return
  end subroutine spllt_starpu_insert_solve_block

  ! solve block StarPU task 
  ! _potrf
  subroutine spllt_starpu_solve_block_cpu_func(buffers, cl_arg) bind(C)
    use spllt_data_mod
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
    use spllt_data_mod
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

#if defined(SPLLT_USE_GPU)

  ! subroutine spllt_starpu_update_between_cuda_func(buffers, cl_arg) bind(C)
  !   use iso_c_binding
  !   use hsl_ma87_double
  !   use spllt_data_mod
  !   use spllt_kernels_mod
  !   implicit none

  !   type(c_ptr), value        :: cl_arg
  !   type(c_ptr), value        :: buffers

  !   type(c_ptr), target :: snode_c, anode_c, a_bc_c
  !   integer, target     :: csrc, csrc2, rsrc, rsrc2 
  !   integer, target     :: scol, dcol
  !   integer, target     :: min_width_blas

  !   type(spllt_node), pointer  :: snode, anode
  !   type(block_type), pointer :: a_bc
  !   integer, dimension(:), allocatable  :: row_list, col_list
  !   integer :: i, nh, cld2, cld1

  !   type(c_ptr), target      :: work_c, dest, src1, src2, ptr1, ptr2
  !   integer, target :: mw, nw, ldw, m, n, ld, m2, n2, ld2 
  !   integer, target :: m1, n1, ld1
  !   real(wp), pointer        :: work(:)
  !   real(wp), pointer        :: lij(:), lik(:), ljk(:)

  !   call spllt_starpu_codelet_unpack_args_update_between(cl_arg, &
  !        & c_loc(snode_c), c_loc(scol), &
  !        & c_loc(anode_c), c_loc(a_bc_c), c_loc(dcol), &
  !        & c_loc(csrc), c_loc(csrc2), &
  !        & c_loc(rsrc), c_loc(rsrc2), &
  !        & c_loc(min_width_blas))

  !   call c_f_pointer(snode_c, snode)
  !   call c_f_pointer(anode_c, anode)
  !   call c_f_pointer(a_bc_c, a_bc)
      
  !   ! get workspace pointer
  !   call starpu_f_get_buffer(buffers, 0, c_loc(work_c), c_loc(mw))
  !   call c_f_pointer(work_c, work,(/mw/))

  !   ! Lij buffer
  !   call starpu_f_get_buffer(buffers, 1, c_loc(dest), c_loc(m), c_loc(n), c_loc(ld))
  !   call c_f_pointer(dest, lij,(/ld*n/))
    
  !   ! Lik buffer
  !   call starpu_f_get_buffer(buffers, 2, c_loc(src2), c_loc(m2), c_loc(n2), c_loc(ld2))
  !   call c_f_pointer(src2, lik,(/ld2*n2/))    
    
  !   ! Ljk buffer
  !   call starpu_f_get_buffer(buffers, 3, c_loc(src1), c_loc(m1), c_loc(n1), c_loc(ld1))
  !   call c_f_pointer(src1, ljk,(/ld1*n1/))    

  !   ! TODO use scratch memory
  !   allocate(col_list(1), row_list(1))

  !   call spllt_cuda_update_between(m, n, a_bc, dcol, anode, &
  !        & n1, scol, snode, &
  !        & lij, &
  !        & ljk(csrc:csrc+csrc2-1), &
  !        & lik(rsrc:rsrc+rsrc2-1), &
  !        & row_list, col_list, work, &
  !        & min_width_blas)

  !   ! TODO use scratch memory
  !   deallocate(col_list, row_list)

  ! end subroutine spllt_starpu_update_between_cuda_func

  subroutine spllt_starpu_update_between_cpu_func(buffers, cl_arg) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use spllt_kernels_mod
    implicit none
    
    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers

    type(c_ptr), target :: snode_c, anode_c, a_bc_c
    integer, target     :: csrc, csrc2, rsrc, rsrc2 
    integer, target     :: scol, dcol
    integer, target     :: min_width_blas

    type(spllt_node), pointer  :: snode, anode
    type(spllt_block), pointer :: a_bc
    integer :: i, nh, cld2, cld1

    type(c_ptr), target      :: rlst_c, clst_c, work_c, dest, src1, src2, ptr1, ptr2
    integer, target :: mw, nw, ldw, m, n, ld, m2, n2, ld2 
    integer, target :: m1, n1, ld1
    integer, target :: rls, cls
    integer, pointer  :: row_list(:), col_list(:)
    real(wp), pointer        :: work(:)
    real(wp), pointer        :: lij(:), lik(:), ljk(:)

    ! write(*,*)"update_between_cpu_func"

    call spllt_starpu_codelet_unpack_args_update_between(cl_arg, &
         & c_loc(snode_c), c_loc(scol), &
         & c_loc(anode_c), c_loc(a_bc_c), c_loc(dcol), &
         & c_loc(csrc), c_loc(csrc2), &
         & c_loc(rsrc), c_loc(rsrc2), &
         & c_loc(min_width_blas))

    call c_f_pointer(snode_c, snode)
    call c_f_pointer(anode_c, anode)
    call c_f_pointer(a_bc_c, a_bc)
      
    ! get workspace pointer
    call starpu_f_get_buffer(buffers, 0, c_loc(work_c), c_loc(mw))
    call c_f_pointer(work_c, work,(/mw/))

    ! get row_list pointer
    call starpu_f_get_buffer(buffers, 1, c_loc(rlst_c), c_loc(rls))
    call c_f_pointer(rlst_c, row_list,(/rls/))

    ! get col_list pointer
    call starpu_f_get_buffer(buffers, 2, c_loc(clst_c), c_loc(cls))
    call c_f_pointer(clst_c, col_list,(/cls/))

    ! Lij buffer
    call starpu_f_get_buffer(buffers, 3, c_loc(dest), c_loc(m), c_loc(n), c_loc(ld))
    call c_f_pointer(dest, lij,(/ld*n/))
    
    ! Lik buffer
    call starpu_f_get_buffer(buffers, 4, c_loc(src2), c_loc(m2), c_loc(n2), c_loc(ld2))
    call c_f_pointer(src2, lik,(/ld2*n2/))    
    
    ! Ljk buffer
    call starpu_f_get_buffer(buffers, 5, c_loc(src1), c_loc(m1), c_loc(n1), c_loc(ld1))
    call c_f_pointer(src1, ljk,(/ld1*n1/))    

    call spllt_update_between(m, n, a_bc, dcol, anode, &
         & n1, scol, snode, &
         & lij, &
         & ljk(csrc:csrc+csrc2-1), &
         & lik(rsrc:rsrc+rsrc2-1), &
         & row_list, col_list, work, &
         & min_width_blas)

  end subroutine spllt_starpu_update_between_cpu_func

#else

  subroutine spllt_starpu_update_between_cpu_func(buffers, cl_arg) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use spllt_kernels_mod
    implicit none
    
    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers

    type(c_ptr), target       :: snode_c, anode_c, a_bc_c
    integer, target           :: csrc, csrc2, rsrc, rsrc2 
    integer, target           :: scol, dcol
    integer, target           :: min_width_blas, nhlik, nhljk

    type(spllt_node), pointer  :: snode, anode
    type(spllt_block), pointer :: a_bc
    integer                   :: i, nh, cld2, cld1

    type(c_ptr), target      :: work_c, dest, src1, src2, ptr1, ptr2
    type(c_ptr), target      :: rlst_c, clst_c
    integer, target          :: mw, nw, ldw, m, n, ld, m2, n2, ld2 
    integer, target          :: m1, n1, ld1
    integer, target          :: rls, cls
    real(wp), pointer        :: work(:)
    ! real(wp), allocatable    :: work(:) ! TODO use StarPU scratch buffer
    real(wp), pointer        :: lij(:), lik(:), ljk(:)
    integer, pointer         :: row_list(:), col_list(:)
    
    ! write(*,*)"update_between_cpu_func"
    call spllt_starpu_codelet_unpack_args_update_between(cl_arg, &
         & c_loc(snode_c), c_loc(scol), &
         & c_loc(anode_c), c_loc(a_bc_c), c_loc(dcol), &
         & c_loc(csrc), c_loc(csrc2), &
         & c_loc(rsrc), c_loc(rsrc2), &
         & c_loc(min_width_blas), &
         & c_loc(nhlik), c_loc(nhljk))

    call c_f_pointer(snode_c, snode)
    call c_f_pointer(anode_c, anode)
    call c_f_pointer(a_bc_c, a_bc)
    ! write(*,*)"nhlik: ", nhlik, ", nhljk: ", nhljk
    ! write(*,*)"min_width_blas: ", min_width_blas
    nh = 0
    call starpu_f_get_buffer(buffers, nh, c_loc(work_c), c_loc(mw))
    call c_f_pointer(work_c, work,(/mw/))
    nh = nh + 1

    ! get row_list pointer
    call starpu_f_get_buffer(buffers, nh, c_loc(rlst_c), c_loc(rls))
    call c_f_pointer(rlst_c, row_list,(/rls/))
    nh = nh + 1

    ! get col_list pointer
    call starpu_f_get_buffer(buffers, nh, c_loc(clst_c), c_loc(cls))
    call c_f_pointer(clst_c, col_list,(/cls/))
    nh = nh + 1

    call starpu_f_get_buffer(buffers, nh, c_loc(dest), c_loc(m), c_loc(n), c_loc(ld))
    call c_f_pointer(dest, lij,(/ld*n/))
    nh = nh + 1

    ! Lik
    cld2 = 0
    do i=1,nhlik       
       if (i.eq.1) then
          ! get pointer on first block in lik
          call starpu_f_get_buffer(buffers, nh, c_loc(src2), c_loc(m2), c_loc(n2), c_loc(ld2))
       else
          call starpu_f_get_buffer(buffers, nh, c_loc(ptr2), c_loc(m2), c_loc(n2), c_loc(ld2))
       end if
       cld2 = cld2 + ld2
       nh = nh + 1
    end do
    call c_f_pointer(src2, lik,(/cld2*n2/))    

    ! Ljk
    cld1 = 0
    do i=1,nhljk
       if (i.eq.1) then
          ! get pointer on first block in ljk
          call starpu_f_get_buffer(buffers, nh, c_loc(src1), c_loc(m1), c_loc(n1), c_loc(ld1))
       else
          call starpu_f_get_buffer(buffers, nh, c_loc(ptr1), c_loc(m1), c_loc(n1), c_loc(ld1))
       end if
       cld1 = cld1 + ld1
       nh = nh + 1
    end do
    call c_f_pointer(src1, ljk,(/cld1*n1/))    
    
    ! TODO use scratch memory
    ! allocate(col_list(1), row_list(1))
    ! allocate(work(m*n))
    
    call spllt_update_between(m, n, a_bc, dcol, anode, &
         & n1, scol, snode, &
         & lij, &
         & ljk(csrc:csrc+csrc2-1), &
         & lik(rsrc:rsrc+rsrc2-1), &
         & row_list, col_list, work, &
         & min_width_blas)

    ! deallocate(col_list, row_list)
    ! deallocate(work)

  end subroutine spllt_starpu_update_between_cpu_func
#endif

  ! init node StarPU task 
  subroutine spllt_starpu_init_node_cpu_func(buffers, cl_arg) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use spllt_kernels_mod
    implicit none

    type(c_ptr), value :: cl_arg
    type(c_ptr), value :: buffers

    integer, target            :: nval, snode
    type(c_ptr), target        :: val_c, fkeep_c
    real(wp), pointer          :: val(:)
    type(spllt_fkeep), pointer  :: fkeep

    call spllt_starpu_codelet_unpack_args_init_node(cl_arg, &
         & c_loc(snode), &  
         & c_loc(val_c), c_loc(nval), c_loc(fkeep_c)) 
    
    call c_f_pointer(val_c, val, (/nval/))
    call c_f_pointer(fkeep_c, fkeep)    

    call spllt_init_node(snode, val, fkeep)
    
  end subroutine spllt_starpu_init_node_cpu_func

  ! init node StarPU CPU kernel  
  subroutine spllt_starpu_subtree_factorize_cpu_func(buffers, cl_arg) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use spllt_kernels_mod
    implicit none

    type(c_ptr), value :: cl_arg
    type(c_ptr), value :: buffers

    type(c_ptr), target       :: val_c, fkeep_c, cntl_c
    integer, target           :: root

    type(spllt_fkeep), pointer  :: fkeep
    integer :: n
    real(wp), pointer         :: val(:)
    type(spllt_options), pointer :: cntl

    type(c_ptr), target       :: buffer_c, map_c, rlst_c, clst_c, work_c
    real(wp), pointer         :: buffer(:), work(:)
    integer, pointer          :: map(:), rlst(:), clst(:)
    integer, target           :: mw, map_sz, rsize, csize, buf_m, buf_n, buf_ld

    call spllt_starpu_codelet_unpack_args_subtree_factorize(cl_arg, &
         & c_loc(root), c_loc(val_c), c_loc(fkeep_c), c_loc(cntl_c))

    call c_f_pointer(fkeep_c, fkeep)    
    n = fkeep%n
    call c_f_pointer(val_c, val, (/n*n/))
    
    call c_f_pointer(cntl_c, cntl)

    call starpu_f_get_buffer(buffers, 0, &
         & c_loc(buffer_c), c_loc(buf_m), c_loc(buf_n), c_loc(buf_ld))
    call c_f_pointer(buffer_c, buffer, (/buf_ld*buf_n/))

    ! get map
    call starpu_f_get_buffer(buffers, 1, c_loc(map_c), c_loc(map_sz))
    call c_f_pointer(map_c, map, (/map_sz/))
    
    ! get row_list pointer
    call starpu_f_get_buffer(buffers, 2, c_loc(rlst_c), c_loc(rsize))
    call c_f_pointer(rlst_c, rlst, (/rsize/))

    ! get col_list pointer
    call starpu_f_get_buffer(buffers, 3, c_loc(clst_c), c_loc(csize))
    call c_f_pointer(clst_c, clst, (/csize/))

    ! get workspace
    call starpu_f_get_buffer(buffers, 4, c_loc(work_c), c_loc(mw))
    call c_f_pointer(work_c, work, (/mw/))

    call spllt_subtree_factorize(root, val, fkeep, buffer, &
       & cntl, map, rlst, clst, work)

  end subroutine spllt_starpu_subtree_factorize_cpu_func

  ! init node StarPU CPU kernel  
  subroutine spllt_starpu_subtree_scatter_block_cpu_func(buffers, cl_arg) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use spllt_kernels_mod
    implicit none

    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers

    integer, target           :: rptr, rptr2, cptr, cptr2
    type(c_ptr), target       :: root_c, anode_c
    integer, target           :: a_rptr, a_cptr
    type(spllt_node), pointer  :: root, anode
    
    type(c_ptr), target       :: buffer_c, dest_c
    integer, target           :: bm, bn, bld 
    integer, target           :: blkm, blkn, blkld 
    real(wp), pointer         :: buffer(:), dest(:)

    integer :: rm, rn, b_sz
    integer :: m, n
    integer :: bsa, ben
    
    call spllt_starpu_codelet_unpack_args_subtree_scatter_block(cl_arg, &
         c_loc(rptr), c_loc(rptr2), c_loc(cptr), c_loc(cptr2), &
         c_loc(root_c), c_loc(a_rptr), c_loc(a_cptr), c_loc(anode_c))
    
    call c_f_pointer(root_c, root)
    call c_f_pointer(anode_c, anode)

    call starpu_f_get_buffer(buffers, 0, &
         & c_loc(buffer_c), c_loc(bm), c_loc(bn), c_loc(bld))
    call c_f_pointer(buffer_c, buffer, (/bld*bn/))

    call starpu_f_get_buffer(buffers, 1, &
         & c_loc(dest_c), c_loc(blkm), c_loc(blkn), c_loc(blkld))
    call c_f_pointer(dest_c, dest, (/blkld*blkn/))

    rm = size(root%index)
    rn = root%en - root%sa + 1
    b_sz = rm-rn    

    m = rptr2-rptr+1
    n = cptr2-cptr+1
    ! write(*,*)"m:", m, ", n:", n
    bsa = (rptr-rn-1)*b_sz+cptr-rn
    ben = (rptr2-rn-1)*b_sz+cptr2-rn

    call spllt_scatter_block(m, n, &
         root%index(rptr:rptr2), &
         root%index(cptr:cptr2), &
         buffer(bsa:ben), b_sz, &
         anode%index(a_rptr), &
         anode%index(a_cptr), &
         dest, &
         blkn)

  end subroutine spllt_starpu_subtree_scatter_block_cpu_func

end module spllt_starpu_factorization_mod
