module spllt_kernels_mod
  implicit none

contains

  !********************************************************************  

  ! TASK_FACTORIZE_BLOCK (uses Lapack routine dpotrf and dtrsm)
  ! A_ii <- L_ii
  !
  subroutine spllt_factor_diag_block(m, n, dest)
    use spllt_mod
    implicit none
    
    integer, intent(in) :: m ! number of rows in dest
    integer, intent(in) :: n ! number of columns in dest
    real(wp), dimension(*), intent(inout) :: dest ! holds block
    ! on diagonal of factor L. It may not be square

    integer :: dpotrf_info ! error flag for dpotrf
    integer :: i, j ! Loop indices

    call dpotrf('Upper', n, dest, n, dpotrf_info)
    ! check for errors
    if(dpotrf_info.ne.0) return

    ! Do dtrsm with any remainder below diagonal block
    if(m.gt.n) then
       call dtrsm('Left', 'Upper', 'Transpose', 'Non-Unit', n, &
            m-n, one, dest, n, dest(1+n*n), n)
    endif

  end subroutine spllt_factor_diag_block  

  !********************************************************************  

  ! TASK_SOLVE_BLOCK
  ! Solve using factorization of diag. block (uses dtrsm)
  ! A_ij <- A_ij A_ii^-1
  ! dest <- dest diag^-1
  !
  subroutine spllt_solve_block(m, n, dest, diag)
    use spllt_mod
    implicit none

    integer, intent(in) :: m ! number of rows in dest
    integer, intent(in) :: n ! number of columns in dest
    real(wp), dimension(*), intent(inout) :: dest ! holds destination block
    real(wp), dimension(*), intent(in)    :: diag ! block
    
    call dtrsm('Left', 'Upper', 'Transpose', 'Non-Unit', n, m, &
         one, diag, n, dest, n)

  end subroutine spllt_solve_block
  
  !*************************************************  

  ! TASK_UPDATE_INTERNAL
  ! A_ik <- A_ik - A_ij A_kj^T
  ! dest <- dest - src2 src1^T
  ! Remember that the blocks are stored by rows.
  ! dest, src1 and src2 all belong to the same node.
  !
  subroutine spllt_update_block(m, n, dest, diag, n1, src1, src2)
    use spllt_mod
    implicit none

    integer, intent(in) :: m ! number of rows in dest
    integer, intent(in) :: n ! number of columns in dest
    real(wp), dimension(*), intent(inout) :: dest ! holds block in L
    ! that is to be updated. 
    logical :: diag ! set to true if dest is the diagonal block
    ! type(block_type), intent(inout) :: blk ! destination block  
    integer, intent(in) :: n1 ! number of columns in src1 and src2
    real(wp), dimension(*), intent(in) :: src1
    real(wp), dimension(*), intent(in) :: src2

    !%%%   integer :: t_start, t_end, this_thread

    ! diag = (blk%dblk.eq.blk%id)

    if(diag) then
       call dsyrk('U', 'T', n, n1, -one, src1, n1, one, dest, n)

       if(m.gt.n) then
          ! Do a dgemm on any remainder
          call dgemm('T', 'N', n, m-n, n1, -one,                   &
               src1, n1, src2(1+n*n1), n1, one, dest(1+n*n), n)
       endif
    else
       ! dest is an off-diagonal block
       call dgemm('T', 'N', n, m, n1, -one, src1, n1, src2, n1, one, dest, n)
    endif

  end subroutine spllt_update_block

  !*************************************************
  
  !   Given a destination block dest, update_between performs the update
  !                     L_dest <-- L_dest - L_rsrc (L_csrc)^T
  !   where L_dest is a submatrix of the block dest of an ancestor
  !   of the node snode and L_rsrc and L_csrc are submatrices of contiguous
  !   rows that belong to the same block column of snode as the block src
  !   (this block col. has local index scol).
  !   The first row of L_rsrc is the first row
  !   of the block column that corresponds to a row in the block dest and the 
  !   last row of L_rsrc is the last row of the block column that corresponds to 
  !   a row in the block dest. Similarly, the first/last row of L_csrc is the
  !   first/last row of the block column that corresponds to a column in the 
  !   block dest. The set of rows and columns of dest thus
  !   determine which two sets of contiguous rows in scol are involved.
  !   Unless the number of entries updated is very small, use BLAS 3 kernel 
  !   gemm or syrk by placing its result in a buffer from which we add the 
  !   update into the appropriate entries of the
  !   destination block dest.

  ! TODO error managment

  subroutine spllt_update_between(m, n, blk, dcol, dnode, n1, scol, snode, dest, & 
       & csrc, rsrc, row_list, col_list, buffer, min_width_blas)
    use spllt_mod
    use hsl_ma87_double
    implicit none

    integer, intent(in) :: m  ! number of rows in destination block
    integer, intent(in) :: n  ! number of columns in destination block
    ! integer(long), intent(in) :: blk ! identifier of destination block
    type(block_type), intent(in) :: blk ! destination block
    integer, intent(in) :: dcol ! index of block column that blk belongs to in dnode
    type(node_type), intent(in) :: dnode ! Node to which blk belongs
    integer :: n1 ! number of columns in source block column
    ! integer(long), intent(in) :: src  ! identifier of block in source block col
    integer, intent(in) :: scol ! index of block column that src belongs to in snode
    type(node_type), intent(in) :: snode ! Node to which src belongs
    real(wp), dimension(*), intent(inout) :: dest ! holds block in L
    ! that is to be updated.
    real(wp), dimension(*), intent(in) :: csrc ! holds csrc block
    real(wp), dimension(*), intent(in) :: rsrc ! holds rsrc block
    ! type(block_type), dimension(:), intent(inout) :: blocks
    ! real(wp), dimension(:), allocatable :: buffer
    real(wp), dimension(:), pointer :: buffer
    integer, dimension(:), allocatable :: row_list ! reallocated to min size m
    integer, dimension(:), allocatable :: col_list ! reallocated to min size n
    ! integer, dimension(:), pointer :: row_list ! reallocated to min size m
    ! integer, dimension(:), pointer :: col_list ! reallocated to min size n
    integer, intent(in) :: min_width_blas      ! Minimum width of source block
         ! before we use an indirect update_between    
    ! type(MA87_control), intent(in) :: control
    ! integer, intent(inout) :: info   
    ! integer, intent(inout) :: st
    

    ! Local scalars
    integer :: cptr ! used in determining the rows that belong to the
    ! source block csrc
    integer :: col_list_sz ! initialise to 0. then increment while
    ! rows involed in csrc are recorded, until
    ! holds the number of columns in blk (= number
    ! of rows in csrc)
    integer :: dcen ! index of end column in dcol
    ! integer :: dcol ! index of block column that blk belongs to in dnode
    integer :: dcsa ! index of first column in dcol
    logical :: diag ! set to true if blk is the diagonal block
    integer :: dptr
    integer :: dptr_sa
    integer :: drsa, dren ! point to first and last rows of destination
    ! block blk
    integer :: i

    integer :: ndiag ! set to int(s1en-s1sa+1) if blk is a
    ! block on diagonal and 0 ow. so is number of triangular rows of update
    integer :: row_list_sz ! initialise to 0. then increment while
    ! rows involed in rsrc are recorded, until
    ! holds the number of rows in blk (= number of rows in rsrc)
    integer :: rptr ! used in determining the rows that belong to the
    ! source block rsrc
    ! integer :: scol ! index of block column that src belongs to in snode
    integer :: s1sa, s1en ! point to the first and last rows of
    ! the block csrc within scol
    integer :: s2sa, s2en ! point to the first and last rows of
    ! the block rsrc within scol
    integer :: size_dnode ! size(dnode%index)
    integer :: size_snode ! size(snode%index)

    ! TODO error managment
    integer :: st
    
    ! set diag to true if blk is block on diagonal
    diag = (blk%dblk.eq.blk%id)
    ! diag = (blocks(blk)%dblk.eq.blocks(blk)%id)

    ! Make a list of incident csrc rows (ie. columns of blk)
    !
    ! Initialize lists
    ! TODO error managment
    if(size(col_list).lt.n) then
       deallocate(col_list, stat=st)
       allocate(col_list(n), stat=st)
       ! if (st.ne.0) go to 10
    endif
    if(size(row_list).lt.m) then
       deallocate(row_list, stat=st)
       allocate(row_list(m), stat=st)
       ! if (st.ne.0) go to 10
    endif

! 10  if (st.ne.0) then
!        info = MA87_ERROR_ALLOCATION
!        call MA87_print_flag(info, control, context='MA87_factor',st=st)
!        return
!     end if

    col_list_sz = 0
    row_list_sz = 0

    size_dnode = size(dnode%index)
    size_snode = size(snode%index)

    ! Find block column dcol of dnode that blk belongs to. The block
    ! cols are numbered locally within dnode as 1,2,3,...

    ! dcol = blk%bcol - blocks(dnode%blk_sa)%bcol + 1

    ! Set dcsa and dcen to hold indices
    ! of start and end columns in dcol (global column indices)
    dcsa = dnode%sa + (dcol-1)*dnode%nb                
    dcen = min(dnode%sa + dcol*dnode%nb-1, dnode%en)

    ! Find block column scol of snode that src belongs to. 
    ! scol = blocks(src)%bcol - blocks(snode%blk_sa)%bcol + 1

    ! Set cptr to point to the first row in csrc
    cptr = 1 + min(snode%en-snode%sa+1, (scol-1)*snode%nb)

    ! loop while row index within scol is less the index
    ! of the first column in blk

    do while(snode%index(cptr).lt.dcsa)
       cptr = cptr + 1
       if(cptr.gt.size_snode) return ! No incident columns
    end do

    ! Set s1sa to point to first row in csrc
    s1sa = cptr 

    ! Now record the rows in csrc. Local row numbers
    ! are held in col_list(1:slen-slsa+1)
    do while(snode%index(cptr).le.dcen)
       col_list_sz = col_list_sz + 1
       col_list(col_list_sz) = snode%index(cptr) - dcsa + 1
       cptr = cptr + 1
       if(cptr.gt.size_snode) exit ! No more rows
    end do

    ! Set slen to point to last row in csrc
    s1en = cptr - 1 

    ! Loop over rsrc rows, building row list. Identify required data, form
    ! outer product of it into buffer.

    ! Find first and last rows of destination block
    i = dcol + blk%id - blk%dblk ! block in snode
    drsa = dnode%index(1 + (i-1)*dnode%nb)
    dren = dnode%index(min(1 + i*dnode%nb - 1, size_dnode))

    ! Find first row in rsrc
    rptr = s1sa
    do while(snode%index(rptr).lt.drsa)
       rptr = rptr + 1
       if(rptr.gt.size_snode) return ! No incident row! Shouldn't happen.
    end do
    s2sa = rptr ! Points to first row in rsrc

    ! Find the first row of destination block
    i = blk%id - blk%dblk + 1 ! row block of blk column
    dptr_sa = 1 + (dcol-1 + i-1)*dnode%nb

    ! Now record the rows in rsrc. Local row numbers
    ! are held in row_list(1:s2en-s2sa+1)
    dptr = dptr_sa ! Pointer for destination block

    do rptr = s2sa, size_snode
       if(snode%index(rptr).gt.dren) exit
       do while(dnode%index(dptr).lt.snode%index(rptr))
          dptr = dptr + 1
       end do
       row_list_sz = row_list_sz + 1
       row_list(row_list_sz) = dptr - dptr_sa + 1
    end do
    s2en = rptr - 1 ! Points to last row in rsrc

    ! if(n1.ge.control%min_width_blas) then
    if(n1.ge.min_width_blas) then
       ! High flop/buffer sz ratio => perform operations into buffer with BLAS

       ! FIXME We probably dont need this as size(buffer) is equal to
       ! maxmn*maxmn which always more than row_list_sz*col_list_sz
       if(size(buffer).lt.row_list_sz*col_list_sz) then
          deallocate(buffer, stat=st)
          allocate(buffer(row_list_sz*col_list_sz), stat=st)
          ! if (st.ne.0) then
          !    info = MA87_ERROR_ALLOCATION
          !    call MA87_print_flag(info, control, context='MA87_factor',st=st)
          !    return
          ! end if
       endif

       if(diag) then
          ! blk is a block on diagonal
          ndiag = int(s1en-s1sa+1)
          call dsyrk('U', 'T', ndiag, n1, -one, csrc,                  &
               n1, zero, buffer, col_list_sz)

          if(s2en-s2sa+1-ndiag.gt.0) then
             call dgemm('T', 'N', ndiag, int(s2en-s2sa+1-ndiag),       &
                  n1, -one, csrc, n1, rsrc(1+n1*ndiag), n1, zero,        &
                  buffer(1+col_list_sz*ndiag), col_list_sz)
          endif
       else
          ! Off diagonal block
          ndiag = 0
          call dgemm('T', 'N', int(s1en-s1sa+1), int(s2en-s2sa+1), n1, &
               -one, csrc, n1, rsrc, n1, zero, buffer, col_list_sz)
       endif

       !
       ! Apply update
       !

       call expand_buffer(dest, n, row_list, row_list_sz, &
            col_list, col_list_sz, ndiag, buffer)

    else
       ! Low flop/buffer ratio => perform update operation directly
       ! set ndiag if blk is a diagonal block
       ndiag = 0
       if(diag) ndiag = int(s1en-s1sa+1)


       call update_direct(n, dest, n1, csrc, rsrc, row_list, row_list_sz, &
            col_list, col_list_sz, ndiag)

    endif

  end subroutine spllt_update_between

  ! init node
  subroutine spllt_init_node(snode, n, val, map, keep)
    use spllt_mod
    use hsl_ma87_double
    implicit none

    integer, intent(in) :: snode
    integer, intent(in) :: n      ! order of matrix 
    real(wp), dimension(*), intent(in) :: val ! user's matrix values
#if defined(SPLLT_USE_OMP)
    integer, pointer, intent(inout) :: map(:)  ! mapping array. Reset for each node
#else
    integer, intent(inout) :: map(n)  ! mapping array. Reset for each node
#endif
    ! so that, if variable (row) i is involved in node,
    ! map(i) is set to its local row index
    type(MA87_keep), intent(inout) :: keep ! on exit, matrix a copied
    ! into relevant part of keep%lfact

    integer(long) :: i, j ! Temporary variable   
    integer(long) :: dblk ! set to keep%nodes(snode)%blk_sa (first block
    ! in snode which is, of course, a diagonal block)
    integer :: l_nb ! set to keep%nodes(snode)%nb
    integer :: sa ! set keep%nodes(snode)%sa
    integer :: en ! set keep%nodes(snode)%en
    integer :: cb ! Temporary variable
    integer :: bcol ! block column
    integer :: sz ! set to size of keep%lfact(nbcol)%lcol
    integer(long) :: offset
    integer :: swidth ! set to keep%blocks(dblk)%blkn (number of columns
    ! in block column to which dblk belongs)
    ! write(*,*)"snode: ", snode
    ! write(*,*)"size map: ", size(map), ", snode: ", snode
    ! write(*,*)"size nodes: ", size(keep%nodes), ", snode: ", snode
    ! Build a map from global to local indices
    do j = 1, size(keep%nodes(snode)%index)
       i = keep%nodes(snode)%index(j)
       map(i) = j - 1
    end do

    ! Fill in keep%lfact by block columns
    dblk = keep%nodes(snode)%blk_sa

    l_nb = keep%nodes(snode)%nb
    sa = keep%nodes(snode)%sa
    en = keep%nodes(snode)%en

    do cb = sa, en, l_nb
       bcol = keep%blocks(dblk)%bcol
       sz = size(keep%lfact(bcol)%lcol)
       ! Zero the block column. 
       keep%lfact(bcol)%lcol(1:sz) = zero

       offset = keep%blocks(dblk)%sa - (cb-sa)*keep%blocks(dblk)%blkn
       swidth = keep%blocks(dblk)%blkn

       do i = 1, keep%lmap(bcol)%len_map
          keep%lfact(bcol)%lcol(keep%lmap(bcol)%map(1,i)) = &
               val(keep%lmap(bcol)%map(2,i))
       end do

       ! move to next block column in snode
       dblk = keep%blocks(dblk)%last_blk + 1
    end do

    return
  end subroutine spllt_init_node

  subroutine spllt_activate_node(snode, keep, fdata)
    use spllt_mod
    use hsl_ma87_double
#if defined(SPLLT_USE_STARPU)
    use  starpu_f_mod
#endif
    implicit none

    type(MA87_keep), target, intent(inout) :: keep 
    type(spllt_data_type), intent(inout) :: fdata
    integer :: snode

    type(node_type), pointer :: node ! node in the atree    
    integer(long) :: blk, dblk
    integer :: nbcol, l_nb, sz, sa, en
    integer :: blkm, blkn, size_bcol
    integer :: i
    integer :: st ! stat parameter
    integer :: ptr

    node => keep%nodes(snode)
    blk = node%blk_sa

    l_nb = node%nb
    sz = (size(node%index) - 1) / l_nb + 1
    sa = node%sa
    en = node%en
    
    size_bcol = 0
    do i = sa, en, l_nb
       ! nbcol = nbcol + 1
       size_bcol = 0
       dblk = blk
       nbcol = keep%blocks(dblk)%bcol
       ! loop over the row blocks
       do blk = dblk, dblk+sz-1
          blkm = keep%blocks(blk)%blkm
          blkn = keep%blocks(blk)%blkn
          size_bcol = size_bcol + blkm*blkn
       end do
       allocate (keep%lfact(nbcol)%lcol(size_bcol),stat=st)
       
       ptr = 1
       do blk = dblk, dblk+sz-1
          blkm = keep%blocks(blk)%blkm
          blkn = keep%blocks(blk)%blkn

          fdata%bc(blk)%blk => keep%blocks(blk)
#if defined(SPLLT_USE_STARPU)
          call starpu_matrix_data_register(fdata%bc(blk)%hdl, fdata%bc(blk)%mem_node, &
               & c_loc(keep%lfact(nbcol)%lcol(ptr)), blkm, blkm, blkn, &
               & int(wp,kind=c_size_t))
#endif
          fdata%bc(blk)%c => keep%lfact(nbcol)%lcol(ptr:ptr+blkm*blkn-1)
          ptr = ptr + blkm*blkn
       end do
       sz = sz - 1
    end do

    return
  end subroutine spllt_activate_node

  ! Build a map of node's blocks
  subroutine spllt_build_rowmap(node, rowmap)
    use hsl_ma87_double
    implicit none

    type(node_type), intent(in) :: node  ! current node in the atree
    integer, dimension(:), intent(out) :: rowmap ! Workarray to hold map from row 
    ! indices to block indices in ancestor node.

    integer :: a_nr ! number of rows in ancestor
    integer :: a_nb ! block size in ancestor
    integer :: rr  ! row index
    integer :: row, arow  ! row index
    integer :: i

    a_nr = size(node%index)
    a_nb = node%nb

    rr = 1
    do row = 1, a_nr, a_nb
       do i = row, min(row+a_nb-1, a_nr)
          arow = node%index(i)
          rowmap(arow) = rr
       end do
       rr = rr + 1
    end do

    return
  end subroutine spllt_build_rowmap
  
end module spllt_kernels_mod
