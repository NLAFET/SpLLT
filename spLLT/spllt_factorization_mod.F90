module spllt_factorization_mod
  use spllt_mod
  implicit none
  
contains

  ! init node
  ! subroutine spllt_factorize_block_task(node)
  !   implicit none

  !   type(node_type), intent(inout) :: node ! node in the atree

  !   return
  ! end subroutine spllt_factorize_block_task

  ! factorize block 
  ! _potrf
  subroutine spllt_factorize_block_task(bc, lfact, prio)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#endif
    implicit none
    
    type(spllt_bc_type), target, intent(inout) :: bc ! block to be factorized    
    type(lfactor), dimension(:), allocatable, intent(inout) :: lfact
    integer, optional :: prio

    integer :: m, n, bcol, sa
    integer(long) :: id
    type(block_type), pointer :: blk ! block to be factorized
    integer :: p

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)
    call spllt_starpu_insert_factorize_block_c(bc%hdl, p)

    ! call spllt_starpu_insert_factorize_block(bc, p)
#else    

    blk => bc%blk

    m    = blk%blkm
    n    = blk%blkn
    id   = blk%id 
    bcol = blk%bcol

    sa   = blk%sa

    ! factorize_block
    ! call factor_diag_block(n, m, id, &
    !      & lfact(bcol)%lcol(sa:sa+n*m-1),   &
    !      & control, flag, detlog(0))

    call spllt_factor_diag_block(m, n, lfact(bcol)%lcol(sa:sa+n*m-1))
#endif
    
    return
  end subroutine spllt_factorize_block_task

  ! _trsm
  subroutine spllt_solve_block_task(bc_kk, bc_ik, lfact, prio)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#endif
    implicit none

    type(spllt_bc_type), intent(inout) :: bc_kk, bc_ik ! block to be factorized    
    type(lfactor), dimension(:), allocatable, intent(inout) :: lfact
    integer, optional :: prio 
    
    integer :: m, n, bcol, sa
    integer :: d_m, d_n, d_sa
    integer(long) :: id, d_id
    integer :: p

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)    
    call spllt_starpu_insert_solve_block_c(bc_kk%hdl, bc_ik%hdl, p)
    ! call spllt_starpu_insert_solve_block(bc_kk, bc_ik, 0)
#else    

    ! bc_kk
    d_m  = bc_kk%blkm
    d_n  = bc_kk%blkn
    d_sa = bc_kk%sa
    d_id = bc_kk%id

    ! bc_ik
    n  = bc_ik%blkn
    m  = bc_ik%blkm
    sa = bc_ik%sa
    id = bc_ik%id
    
    ! bcol is block column that blk and dblk belong to
    bcol = bc_kk%bcol    

    ! solve_block task
    ! call solv_col_block(m, n, id, & 
    !      & lfact(bcol)%lcol(sa:sa+n*m-1), &
    !      & d_id, lfact(bcol)%lcol(d_sa:d_sa+d_n*d_m), &
    !      & control)

    call spllt_solve_block(m, n, &
         & lfact(bcol)%lcol(sa:sa+n*m-1), & 
         & lfact(bcol)%lcol(d_sa:d_sa+d_n*d_m))

#endif

    return
  end subroutine spllt_solve_block_task

  ! syrk/gemm (same node)
  ! A_ij <- A_ij - A_ik A_jk^T
  subroutine spllt_update_block_task(bc_ik, bc_jk, bc_ij, lfact, prio)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#endif
    implicit none
    
    ! type(block_type), intent(inout) :: bc_ik, bc_jk, bc_ij ! block to be updated    
    type(spllt_bc_type), intent(inout) :: bc_ik, bc_jk, bc_ij
    type(lfactor), dimension(:), allocatable, intent(inout) :: lfact
    integer, optional :: prio 

    integer :: n1, m1, sa1, n2, m2, sa2, n, m, sa
    integer :: bcol1, bcol, bcol2
    integer :: p, d

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)
    if (bc_ij%blk%dblk .eq. bc_ij%blk%id) then
       d = 1
    else
       d = 0
    end if

    call spllt_starpu_insert_update_block_c(bc_ik%hdl, bc_jk%hdl, bc_ij%hdl, d, p)
#else

    ! bc_ik
    n1  = bc_ik%blk%blkn
    m1  = bc_ik%blk%blkm
    sa1 = bc_ik%blk%sa
    bcol1 = bc_ik%blk%bcol
    
    ! bc_jk
    n2  = bc_jk%blk%blkn
    m2  = bc_jk%blk%blkm
    sa2 = bc_jk%blk%sa
    bcol2 = bc_jk%blk%bcol
    
    ! bc_ij
    n  = bc_ij%blk%blkn
    m  = bc_ij%blk%blkm
    sa = bc_ij%blk%sa
    
    bcol = bc_ij%blk%bcol

    ! write(*,*)"bc_ik id: ", bc_ik%id, ", m: ", m1, ", n: ", n1
    ! write(*,*)"bc_jk id: ", bc_jk%id, ", m: ", m2, ", n: ", n2
    ! write(*,*)"bc_ij id: ", bc_ij%id, ", m: ", m, ", n: ", n
    
    ! write(*,*) "size lfact(bcol1)%lcol: ", size(lfact(bcol1)%lcol)
    ! write(*,*) "sa2+n1*m2-1: ", sa2+n1*m2-1

    call spllt_update_block(m, n, &
         & lfact(bcol)%lcol(sa:sa+n*m-1), &
         & bc_ij%blk%dblk.eq.bc_ij%blk%id, n1, &
         & lfact(bcol1)%lcol(sa2:sa2+n1*m2-1), &
         & lfact(bcol1)%lcol(sa1:sa1+n1*m1-1))

#endif

    return
  end subroutine spllt_update_block_task

  ! syrk/gemm (inter-node)
  subroutine spllt_update_between_task(bc, snode, a_bc, anode, &
       & csrc, rsrc, &
       & row_list, col_list, workspace, &
       & lfact, blocks, &
       & control)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#endif
    implicit none

    ! type(block_type), intent(inout)     :: a_bc ! dest block
    ! type(block_type), intent(in)        :: bc ! src block
    type(spllt_bc_type), intent(inout)     :: a_bc ! dest block
    type(spllt_bc_type), intent(in)        :: bc ! src block
    type(node_type)                     :: snode ! src node
    type(node_type)                     :: anode ! dest node
    integer :: csrc(2), rsrc(2) ! used for update_between tasks to
    integer, dimension(:), allocatable  :: row_list, col_list
    type(spllt_bc_type) :: workspace
    ! real(wp), dimension(:), allocatable :: buffer ! update_buffer workspace
    type(block_type), dimension(:)      :: blocks ! block info. 
    type(lfactor), dimension(:), allocatable, intent(inout) :: lfact
    type(MA87_control), intent(in) :: control     
    ! integer, intent(out) :: st ! TODO error managment

    ! integer :: info ! TODO error managment
    
    integer :: m, n, n1, sa
    integer :: bcol, bcol1
    integer(long) :: id, id1
    
    m  = a_bc%blk%blkm
    n  = a_bc%blk%blkn
    id = a_bc%blk%id
    sa = a_bc%blk%sa

    bcol = a_bc%blk%bcol
    
    n1  = bc%blk%blkn
    id1 = bc%blk%id

    bcol1 = bc%blk%bcol

    ! write(*,*) "update_between_task"
    ! call update_between(m, n, id, anode, &
    !      & n1, id1, snode, &
    !      & lfact(bcol)%lcol(sa:sa+m*n-1), &
    !      & lfact(bcol1)%lcol(csrc(1):csrc(1)+csrc(2)-1), &
    !      & lfact(bcol1)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1), &
    !      & blocks, row_list, col_list, buffer, &
    !      & control, info, st)

    call spllt_update_between(m, n, id, anode, &
         & n1, id1, snode, &
         & lfact(bcol)%lcol(sa:sa+m*n-1), &
         & lfact(bcol1)%lcol(csrc(1):csrc(1)+csrc(2)-1), &
         & lfact(bcol1)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1), &
         & blocks, row_list, col_list, workspace%c, &
         & control%min_width_blas)

    return
  end subroutine spllt_update_between_task

end module spllt_factorization_mod
