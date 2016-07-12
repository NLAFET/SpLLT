module spllt_kernels_mod
  implicit none

contains

  !********************************************************************  

  ! TASK_FACTORIZE_BLOCK (uses Lapack routine dpotrf and dtrsm)
  ! A_ii <- L_ii
  !
  subroutine spllt_factor_diag_block(m, n, dest)
    use spllt_data_mod
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

  ! C wrapper

  subroutine spllt_factor_diag_block_c(m, n, bc_c) bind(C)
    use iso_c_binding
    use spllt_data_mod
    implicit none

    integer(c_int), value :: m, n
    type(c_ptr), value :: bc_c
    
    real(wp), pointer :: bc(:) ! holds block

    call c_f_pointer(bc_c, bc,(/m*n/))
    
    call spllt_factor_diag_block(m, n, bc)

    return
  end subroutine spllt_factor_diag_block_c

  !********************************************************************  

  ! TASK_SOLVE_BLOCK
  ! Solve using factorization of diag. block (uses dtrsm)
  ! A_ij <- A_ij A_ii^-1
  ! dest <- dest diag^-1
  !
  subroutine spllt_solve_block(m, n, dest, diag)
    use spllt_data_mod
    implicit none

    integer, intent(in) :: m ! number of rows in dest
    integer, intent(in) :: n ! number of columns in dest
    real(wp), dimension(*), intent(inout) :: dest ! holds destination block
    real(wp), dimension(*), intent(in)    :: diag ! block
    
    call dtrsm('Left', 'Upper', 'Transpose', 'Non-Unit', n, m, &
         one, diag, n, dest, n)

  end subroutine spllt_solve_block

  ! C wrapper

  subroutine spllt_solve_block_c(m, n, bc_kk_c, bc_ik_c) bind(C)
    use iso_c_binding
    use spllt_data_mod
    implicit none

    integer(c_int), value :: m ! number of rows in dest
    integer(c_int), value :: n ! number of columns in dest
    type(c_ptr), value    :: bc_kk_c ! holds destination block
    type(c_ptr), value    :: bc_ik_c ! block
    
    real(wp), pointer :: bc_kk(:), bc_ik(:) 
    
    call c_f_pointer(bc_kk_c, bc_kk, (/n*n/))
    call c_f_pointer(bc_ik_c, bc_ik, (/m*n/))

    call spllt_solve_block(m, n, bc_ik, bc_kk)
    
    return
  end subroutine spllt_solve_block_c
  
  !*************************************************  

  ! TASK_UPDATE_INTERNAL
  ! A_ik <- A_ik - A_ij A_kj^T
  ! dest <- dest - src2 src1^T
  ! Remember that the blocks are stored by rows.
  ! dest, src1 and src2 all belong to the same node.
  !
  subroutine spllt_update_block(m, n, dest, diag, n1, src1, src2)
    use spllt_data_mod
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

  ! C wrapper

  subroutine spllt_update_block_c(m, n, dest_c, isDiag, n1, src1_c, src2_c) bind(C)
    use iso_c_binding
    use spllt_data_mod
    implicit none

    integer(c_int), value :: m ! number of rows in dest
    integer(c_int), value :: n ! number of columns in dest
    type(c_ptr),  value :: dest_c ! holds block in L
    integer(c_int), value :: isDiag ! set to true if dest is the diagonal block
    integer(c_int), value :: n1 ! number of columns in src1 and src2
    type(c_ptr), value :: src1_c ! holds block in L
    type(c_ptr), value :: src2_c ! holds block in L

    real(wp), pointer :: dest(:), src1(:), src2(:) 
    logical :: diag
    
    diag = .true.
    if (isDiag .le. 0) diag = .false.
    
    call c_f_pointer(dest_c, dest, (/m*n/))

    call c_f_pointer(src1_c, src1, (/n*n1/))
    call c_f_pointer(src2_c, src2, (/m*n1/))    

    call spllt_update_block(m, n, dest, diag, n1, src1, src2) 

    return
  end subroutine spllt_update_block_c

  subroutine spllt_update_between_block(m, n, dest, n1, csrc, rsrc, &
       & dcol, scol, &
       & blk_dest, dnode, blk_csrc, blk_rsrc, snode, &
       & min_width_blas, row_list, col_list, buffer)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer, intent(in) :: m  ! number of rows in destination block
    integer, intent(in) :: n  ! number of columns in destination block
    integer :: n1 ! number of columns in source block column
    real(wp), dimension(*), intent(inout) :: dest ! holds block in L
    ! that is to be updated.
    real(wp), dimension(*), intent(in) :: csrc ! holds csrc block
    real(wp), dimension(*), intent(in) :: rsrc ! holds rsrc block
    type(block_type), intent(in) :: blk_dest ! destination block
    type(block_type), intent(in) :: blk_csrc ! destination block
    type(block_type), intent(in) :: blk_rsrc ! destination block
    type(node_type), intent(in) :: dnode ! Node to which blk belongs
    type(node_type), intent(in) :: snode ! Node to which src belongs
    integer, intent(in) :: dcol ! index of block column that blk belongs to in dnode
    integer, intent(in) :: scol ! index of block column that blk belongs to in node
    integer, intent(in) :: min_width_blas      ! Minimum width of source block
         ! before we use an indirect update_between    
    integer, dimension(:), allocatable :: row_list ! reallocated to min size m
    integer, dimension(:), allocatable :: col_list ! reallocated to min size n
    real(wp), dimension(:), pointer :: buffer

    ! Local
    logical :: diag ! set to true if blk is the diagonal block
    integer :: size_dnode ! size(dnode%index)
    integer :: size_snode ! size(snode%index)
    integer :: col_list_sz ! initialise to 0. then increment while
    ! rows involed in csrc are recorded, until
    ! holds the number of columns in blk (= number
    ! of rows in csrc)
    integer :: row_list_sz ! initialise to 0. then increment while
    ! rows involed in rsrc are recorded, until
    ! holds the number of rows in blk (= number of rows in rsrc)
    integer :: dcen ! index of end column in dcol
    ! integer :: dcol ! index of block column that blk belongs to in dnode
    integer :: dcsa ! index of first column in dcol
    integer :: cptr ! used in determining the rows that belong to the
    ! source block csrc
    integer :: rptr ! used in determining the rows that belong to the
    ! source block rsrc
    integer :: drsa, dren ! point to first and last rows of destination
    ! block blk
    integer :: s1sa, s1en ! point to the first and last rows of
    ! the block csrc within scol
    integer :: s2sa, s2en ! point to the first and last rows of
    ! the block rsrc within scol
    integer :: drow, srow
    integer :: dptr_sa
    integer :: dptr
    integer :: ndiag ! set to int(s1en-s1sa+1) if blk is a
    ! block on diagonal and 0 ow. so is number of triangular rows of update
    integer :: bcptr, brptr

    ! TODO error managment
    integer :: st

    ! set diag to true if blk is block on diagonal
    diag = (blk_dest%dblk.eq.blk_dest%id)

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

    col_list_sz = 0
    row_list_sz = 0

    size_dnode = size(dnode%index)
    size_snode = size(snode%index)

    ! Set dcsa and dcen to hold indices
    ! of start and end columns in dcol (global column indices)
    dcsa = dnode%sa + (dcol-1)*dnode%nb                
    dcen = min(dnode%sa + dcol*dnode%nb-1, dnode%en)

    ! Find first and last rows of csrc block
    srow = scol + blk_csrc%id - blk_csrc%dblk ! block in snode
    cptr = 1 + (srow-1)*snode%nb

    ! loop while row index within scol is less the index
    ! of the first column in blk
    
    do while( snode%index(cptr) .lt. dcsa )
       cptr = cptr + 1
       if( cptr .gt. min(size_snode, srow*snode%nb) ) return ! No incident columns
    end do

    ! Set s1sa to point to first row in csrc
    s1sa = cptr
    ! first row within csrc block
    bcptr = 1 + (cptr - (srow-1)*snode%nb - 1)*n1

    ! Now record the rows in csrc. Local row numbers
    ! are held in col_list(1:slen-slsa+1)
    do while( snode%index(cptr) .le. dcen )
       col_list_sz = col_list_sz + 1
       col_list(col_list_sz) = snode%index(cptr) - dcsa + 1
       cptr = cptr + 1
       if(cptr .gt. min(size_snode, srow*snode%nb) ) exit ! No more rows
    end do

    ! Set slen to point to last row in csrc
    s1en = cptr - 1 

    ! Loop over rsrc rows, building row list. Identify required data, form
    ! outer product of it into buffer.
    
    ! Find first and last rows of destination block
    drow = dcol + blk_dest%id - blk_dest%dblk ! block in snode
    drsa = dnode%index(1 + (drow-1)*dnode%nb)
    dren = dnode%index(min(1 + drow*dnode%nb - 1, size_dnode))

    ! Find first row in rsrc
    srow = scol + blk_rsrc%id - blk_rsrc%dblk ! block in snode
    ! first row within rsrc block
    rptr = 1 + (srow-1)*snode%nb

    do while( snode%index(rptr) .lt. drsa)
       rptr = rptr + 1
       if( rptr .gt. min(size_snode, srow*snode%nb) ) return ! No incident row! Shouldn't happen.
    end do
    s2sa = rptr ! Points to first row in rsrc

    brptr = 1 + (rptr - (srow-1)*snode%nb - 1)*n1

    ! Find the first row of destination block
    drow = dcol + blk_dest%id - blk_dest%dblk ! row block of blk column
    dptr_sa = 1 + (drow-1)*dnode%nb

    ! Now record the rows in rsrc. Local row numbers
    ! are held in row_list(1:s2en-s2sa+1)
    dptr = dptr_sa ! Pointer for destination block

    do rptr = s2sa, min( size_snode, srow*snode%nb )
       if( snode%index(rptr) .gt. dren ) exit
       do while(dnode%index(dptr) .lt. snode%index(rptr))
          dptr = dptr + 1
       end do
       row_list_sz = row_list_sz + 1
       row_list(row_list_sz) = dptr - dptr_sa + 1
    end do
    s2en = rptr - 1 ! Points to last row in rsrc
    
    diag = diag .and. (blk_csrc%id .eq. blk_rsrc%id)

    if(n1.ge.min_width_blas) then

       ! write(*,*)'s_nb:',snode%nb, 'd_nb:',dnode%nb
       ! write(*,*)'m:', m, ', n:', n, ', n1:', n1, ', diag:', diag
       ! write(*,*)'blk_jk:', blk_csrc%id, 'blk_ik:', blk_rsrc%id, 'blk_ij:', blk_dest%id
       ! write(*,*)'col_list_sz:', col_list_sz

       if(diag) then
          ! blk is a block on diagonal
          ndiag = int(s1en-s1sa+1)
          ! write(*,*) 'cptr:', s1sa, 'cptr2:', s1en, 'rptr:', s2sa, 'rptr2:', s2en, 'ndiag:', ndiag
          ! write(*,*) 'bcptr:', bcptr, 'brptr:', brptr, 'ndiag:',ndiag

          call dsyrk('U', 'T', ndiag, n1, -one, csrc(bcptr), &
               n1, zero, buffer, col_list_sz)

          if( s2en-s2sa+1-ndiag .gt. 0 ) then
             call dgemm('T', 'N', ndiag, int(s2en-s2sa+1-ndiag), &
                  n1, -one, csrc(bcptr), n1, rsrc(brptr+n1*ndiag), n1, zero, &
                  buffer(1+col_list_sz*ndiag), col_list_sz)
          endif
       else
          ! Off diagonal block
          ndiag = 0
          call dgemm('T', 'N', int(s1en-s1sa+1), int(s2en-s2sa+1), n1, &
               -one, csrc(bcptr), n1, rsrc(brptr), n1, zero, buffer, col_list_sz)
       endif

       !
       ! Apply update
       !    

       ! ndiag = min(ndiag, row_list_sz)
       call expand_buffer(dest, n, row_list, row_list_sz, &
            col_list, col_list_sz, ndiag, buffer)

    else
       ! Low flop/buffer ratio => perform update operation directly
       ! set ndiag if blk is a diagonal block
       ndiag = 0
       if(diag) ndiag = int(s1en-s1sa+1)


       call update_direct(n, dest, n1, csrc(bcptr), rsrc(brptr), row_list, row_list_sz, &
            col_list, col_list_sz, ndiag)

    endif


    return
  end subroutine spllt_update_between_block

  ! C wrapper
  
  subroutine spllt_update_between_block_c(m, n, dest_c, n1, csrc_c, rsrc_c, &
       & dcol, scol, &
       & blk_dest_c, dnode_c, blk_csrc_c, blk_rsrc_c, snode_c, &
       & min_width_blas, buffer_c) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer(c_int), value  :: m ! number of rows in dest
    integer(c_int), value  :: n ! number of columns in dest    
    type(c_ptr), value     :: dest_c
    integer(c_int), value  :: n1 ! number of columns in src    
    type(c_ptr), value     :: csrc_c
    type(c_ptr), value     :: rsrc_c
    integer(c_int), value  :: dcol
    integer(c_int), value  :: scol
    type(c_ptr), value     :: blk_dest_c
    type(c_ptr), value     :: dnode_c
    type(c_ptr), value     :: blk_csrc_c
    type(c_ptr), value     :: blk_rsrc_c
    type(c_ptr), value     :: snode_c
    integer(c_int), value  :: min_width_blas
    type(c_ptr), value     :: buffer_c
    
    real(wp), pointer, dimension(:) :: dest, csrc, rsrc
    integer, dimension(:), allocatable :: row_list ! reallocated to min size m
    integer, dimension(:), allocatable :: col_list ! reallocated to min size n
    type(node_type), pointer :: dnode, snode
    type(block_type), pointer :: blk_dest, blk_csrc, blk_rsrc
    real(wp), pointer, dimension(:) :: buffer
    
    call c_f_pointer(dest_c, dest, (/m*n/))
    call c_f_pointer(csrc_c, csrc, (/n*n1/))
    call c_f_pointer(rsrc_c, rsrc, (/m*n1/))    

    call c_f_pointer(dnode_c, dnode)    
    call c_f_pointer(snode_c, snode)    

    call c_f_pointer(blk_dest_c, blk_dest)    
    call c_f_pointer(blk_csrc_c, blk_csrc)    
    call c_f_pointer(blk_rsrc_c, blk_rsrc)    
    
    call c_f_pointer(buffer_c, buffer, (/m*n/))    

    allocate(row_list(1), col_list(1))

    call spllt_update_between_block(m, n, dest, n1, csrc, rsrc, &
         & dcol, scol, &
         & blk_dest, dnode, blk_csrc, blk_rsrc, snode, &
         & min_width_blas, row_list, col_list, buffer)

    deallocate(row_list, col_list)

    return
  end subroutine spllt_update_between_block_c

  ! Compute mapping between blocks for updating
  subroutine spllt_update_between_compute_map(blk, dcol, dnode, scol, snode, &
       & row_list, col_list, row_list_sz, col_list_sz, &
       & s1sa, s1en, s2sa, s2en)
    use iso_c_binding
    use hsl_ma87_double
    implicit none

    type(block_type), intent(in) :: blk ! destination block
    integer(c_int), intent(in) :: dcol ! index of block column that blk belongs to in dnode
    type(node_type), intent(in) :: dnode ! destination node
    integer(c_int), intent(in) :: scol ! index of block column that src belongs to in snode
    type(node_type), intent(in) :: snode ! Node to which src belongs

    integer(c_int), dimension(:), pointer, intent(inout) :: row_list ! reallocated to min size m
    integer(c_int), dimension(:), pointer, intent(inout) :: col_list ! reallocated to min size n
    integer(c_int), intent(inout) :: row_list_sz
    integer(c_int), intent(inout) :: col_list_sz
    integer(c_int), intent(inout) :: s1sa, s1en ! point to the first and last rows of
    ! the block csrc within scol
    integer(c_int), intent(inout) :: s2sa, s2en ! point to the first and last rows of
    ! the block rsrc within scol

    ! Local scalars
    integer :: cptr ! used in determining the rows that belong to the
    ! source block csrc
    integer :: rptr ! used in determining the rows that belong to the
    ! source block rsrc
    integer :: dcen ! index of end column in dcol
    ! integer :: dcol ! index of block column that blk belongs to in dnode
    integer :: dcsa ! index of first column in dcol
    integer :: size_dnode ! size(dnode%index)
    integer :: size_snode ! size(snode%index)
    integer :: dptr
    integer :: dptr_sa
    integer :: drsa, dren ! point to first and last rows of destination
    ! block blk
    integer :: i

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

  end subroutine spllt_update_between_compute_map

  subroutine spllt_update_between_compute_map_c(blk_c, dcol, dnode_c, scol, snode_c, &
       & row_list_c, col_list_c, row_list_sz, col_list_sz, &
       & s1sa, s1en, s2sa, s2en) bind(C)
    use iso_c_binding
    use hsl_ma87_double
    implicit none

    type(c_ptr), value, intent(in) :: blk_c ! destination block C ptr
    integer(c_int), value, intent(in) :: dcol ! destination block-column
    type(c_ptr), value, intent(in) :: dnode_c ! destination node C ptr
    integer(c_int), value, intent(in) :: scol ! source block-column
    type(c_ptr), value, intent(in) :: snode_c ! source node C ptr
    type(c_ptr), value, intent(in) :: row_list_c
    type(c_ptr), value, intent(in) :: col_list_c
    integer(c_int), intent(inout) :: row_list_sz
    integer(c_int), intent(inout) :: col_list_sz
    integer(c_int), intent(inout) :: s1sa, s1en ! point to the first and last rows of
    ! the block csrc within scol
    integer(c_int), intent(inout) :: s2sa, s2en ! point to the first and last rows of
    ! the block rsrc within scol

    type(block_type), pointer :: blk ! destination block
    type(node_type), pointer  :: dnode ! destination node
    type(node_type), pointer  :: snode ! source node
    integer(c_int), dimension(:), pointer :: row_list ! reallocated to min size m
    integer(c_int), dimension(:), pointer :: col_list ! reallocated to min size n
    integer :: m, n
        
    call c_f_pointer(blk_c, blk)
    call c_f_pointer(dnode_c, dnode)
    call c_f_pointer(snode_c, snode)

    m = blk%blkm
    n = blk%blkn

    call c_f_pointer(row_list_c, row_list, (/m/))
    call c_f_pointer(col_list_c, col_list, (/n/))

    call spllt_update_between_compute_map(blk, dcol, dnode, scol, snode, &
         & row_list, col_list, row_list_sz, col_list_sz, &
         & s1sa, s1en, s2sa, s2en)

  end subroutine spllt_update_between_compute_map_c

#if defined(SPLLT_USE_GPU)

  subroutine spllt_cuda_update_between(m, n, blk, dcol, dnode, n1, scol, snode, dest, & 
       & csrc, rsrc, row_list, col_list, buffer, min_width_blas)
    use spllt_data_mod
    use hsl_ma87_double
    use magma
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
    integer, dimension(:), allocatable :: row_list ! reallocated to min size m
    integer, dimension(:), allocatable :: col_list ! reallocated to min size n
    ! integer, dimension(:), pointer :: row_list ! reallocated to min size m
    ! integer, dimension(:), pointer :: col_list ! reallocated to min size n
    real(wp), dimension(:), pointer :: buffer
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

    if(diag) then
       ! blk is a block on diagonal
       ndiag = int(s1en-s1sa+1)
       
       
    end if

    ! ! if(n1.ge.control%min_width_blas) then
    ! if(n1.ge.min_width_blas) then
    !    ! High flop/buffer sz ratio => perform operations into buffer with BLAS

    !    ! FIXME We probably dont need this as size(buffer) is equal to
    !    ! maxmn*maxmn which always more than row_list_sz*col_list_sz
    !    if(size(buffer).lt.row_list_sz*col_list_sz) then
    !       deallocate(buffer, stat=st)
    !       allocate(buffer(row_list_sz*col_list_sz), stat=st)
    !       ! if (st.ne.0) then
    !       !    info = MA87_ERROR_ALLOCATION
    !       !    call MA87_print_flag(info, control, context='MA87_factor',st=st)
    !       !    return
    !       ! end if
    !    endif

    !    if(diag) then
    !       ! blk is a block on diagonal
    !       ndiag = int(s1en-s1sa+1)
    !       call dsyrk('U', 'T', ndiag, n1, -one, csrc,                  &
    !            n1, zero, buffer, col_list_sz)

    !       if(s2en-s2sa+1-ndiag.gt.0) then
    !          call dgemm('T', 'N', ndiag, int(s2en-s2sa+1-ndiag),       &
    !               n1, -one, csrc, n1, rsrc(1+n1*ndiag), n1, zero,        &
    !               buffer(1+col_list_sz*ndiag), col_list_sz)
    !       endif
    !    else
    !       ! Off diagonal block
    !       ndiag = 0
    !       call dgemm('T', 'N', int(s1en-s1sa+1), int(s2en-s2sa+1), n1, &
    !            -one, csrc, n1, rsrc, n1, zero, buffer, col_list_sz)
    !    endif

    !    !
    !    ! Apply update
    !    !

    !    call expand_buffer(dest, n, row_list, row_list_sz, &
    !         col_list, col_list_sz, ndiag, buffer)

    ! else
    !    ! Low flop/buffer ratio => perform update operation directly
    !    ! set ndiag if blk is a diagonal block
    !    ndiag = 0
    !    if(diag) ndiag = int(s1en-s1sa+1)


    !    call update_direct(n, dest, n1, csrc, rsrc, row_list, row_list_sz, &
    !         col_list, col_list_sz, ndiag)

    ! endif

  end subroutine spllt_cuda_update_between
#endif

  ! fastest version of expand_buffer kernel 
  subroutine spllt_expand_buffer(a, blkn, row_list, rls, col_list, cls, ndiag, &
       & buffer)
    use spllt_data_mod
    implicit none    

    real(wp), dimension(*), intent(inout) :: a ! holds L
    integer, intent(in) :: blkn ! number of cols in destination block
    integer, intent(in) :: rls ! size of row_list
    integer, intent(in) :: row_list(rls)
    integer, intent(in) :: cls ! size of col_list
    integer, intent(in) :: col_list(cls)
    integer, intent(in) :: ndiag ! Number of triangular rows of update
    real(wp), intent(in) :: buffer(rls*cls)

    integer :: i, j, rptr, cptr, imax, cmax
    integer :: k0, k1, k2, k3

    ! write(*,*) 'TESTS'

    do j = 1, rls
       rptr = (j-1)*cls + 1
       cptr = (row_list(j)-1) * blkn
       imax = cls
       if(j.le.ndiag) imax = j
       cmax = 4*(imax/4)
       do i = 1, cmax, 4
          k0 = cptr + col_list(i+0)
          k1 = cptr + col_list(i+1)
          k2 = cptr + col_list(i+2)
          k3 = cptr + col_list(i+3)
          a(k0) = a(k0) + buffer(rptr+0)
          a(k1) = a(k1) + buffer(rptr+1)
          a(k2) = a(k2) + buffer(rptr+2)
          a(k3) = a(k3) + buffer(rptr+3)
          rptr = rptr + 4
       end do
       do i = cmax+1, imax
          k0 = cptr + col_list(i)
          a(k0) = a(k0) + buffer(rptr)
          rptr = rptr + 1
       end do
    end do

  end subroutine spllt_expand_buffer

  subroutine spllt_expand_buffer_c(a_c, blkm, blkn, row_list_c, rls, col_list_c, cls, ndiag, &
       & buffer_c) bind(C)
    use iso_c_binding
    use spllt_data_mod
    implicit none    
    
    type(c_ptr), value    :: a_c ! holds L
    integer(c_int), value :: blkm ! number of rows in destination block
    integer(c_int), value :: blkn ! number of cols in destination block
    integer(c_int), value :: rls ! size of row_list
    type(c_ptr), value    :: row_list_c
    integer(c_int), value :: cls ! size of col_list
    type(c_ptr), value    :: col_list_c
    integer(c_int), value :: ndiag ! Number of triangular rows of update
    type(c_ptr), value    :: buffer_c

    real(wp), dimension(:), pointer :: a ! holds L    
    integer(c_int), dimension(:), pointer :: row_list ! holds L    
    integer(c_int), dimension(:), pointer :: col_list ! holds L    
    real(wp), dimension(:), pointer :: buffer ! holds L    

    call c_f_pointer(a_c, a, (/blkm*blkn/))
    call c_f_pointer(row_list_c, row_list, (/rls/))
    call c_f_pointer(col_list_c, col_list, (/cls/))
    call c_f_pointer(buffer_c, buffer, (/rls*cls/))

    call spllt_expand_buffer(a, blkn, row_list, rls, col_list, cls, ndiag, &
       & buffer)

  end subroutine spllt_expand_buffer_c

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
    use spllt_data_mod
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
    integer, dimension(:), pointer :: row_list ! reallocated to min size m
    integer, dimension(:), pointer :: col_list ! reallocated to min size n
    ! integer, dimension(:), pointer :: row_list ! reallocated to min size m
    ! integer, dimension(:), pointer :: col_list ! reallocated to min size n
    real(wp), dimension(:), pointer :: buffer
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

    ! col_list_sz = 0
    ! row_list_sz = 0

    ! size_dnode = size(dnode%index)
    ! size_snode = size(snode%index)

    ! ! Find block column dcol of dnode that blk belongs to. The block
    ! ! cols are numbered locally within dnode as 1,2,3,...

    ! ! dcol = blk%bcol - blocks(dnode%blk_sa)%bcol + 1

    ! ! Set dcsa and dcen to hold indices
    ! ! of start and end columns in dcol (global column indices)
    ! dcsa = dnode%sa + (dcol-1)*dnode%nb                
    ! dcen = min(dnode%sa + dcol*dnode%nb-1, dnode%en)

    ! ! Find block column scol of snode that src belongs to. 
    ! ! scol = blocks(src)%bcol - blocks(snode%blk_sa)%bcol + 1

    ! ! Set cptr to point to the first row in csrc
    ! cptr = 1 + min(snode%en-snode%sa+1, (scol-1)*snode%nb)

    ! ! loop while row index within scol is less the index
    ! ! of the first column in blk

    ! do while(snode%index(cptr).lt.dcsa)
    !    cptr = cptr + 1
    !    if(cptr.gt.size_snode) return ! No incident columns
    ! end do

    ! ! Set s1sa to point to first row in csrc
    ! s1sa = cptr 

    ! ! Now record the rows in csrc. Local row numbers
    ! ! are held in col_list(1:slen-slsa+1)
    ! do while(snode%index(cptr).le.dcen)
    !    col_list_sz = col_list_sz + 1
    !    col_list(col_list_sz) = snode%index(cptr) - dcsa + 1
    !    cptr = cptr + 1
    !    if(cptr.gt.size_snode) exit ! No more rows
    ! end do

    ! ! Set slen to point to last row in csrc
    ! s1en = cptr - 1 

    ! ! Loop over rsrc rows, building row list. Identify required data, form
    ! ! outer product of it into buffer.

    ! ! Find first and last rows of destination block
    ! i = dcol + blk%id - blk%dblk ! block in snode
    ! drsa = dnode%index(1 + (i-1)*dnode%nb)
    ! dren = dnode%index(min(1 + i*dnode%nb - 1, size_dnode))

    ! ! Find first row in rsrc
    ! rptr = s1sa
    ! do while(snode%index(rptr).lt.drsa)
    !    rptr = rptr + 1
    !    if(rptr.gt.size_snode) return ! No incident row! Shouldn't happen.
    ! end do
    ! s2sa = rptr ! Points to first row in rsrc

    ! ! Find the first row of destination block
    ! i = blk%id - blk%dblk + 1 ! row block of blk column
    ! dptr_sa = 1 + (dcol-1 + i-1)*dnode%nb

    ! ! Now record the rows in rsrc. Local row numbers
    ! ! are held in row_list(1:s2en-s2sa+1)
    ! dptr = dptr_sa ! Pointer for destination block

    ! do rptr = s2sa, size_snode
    !    if(snode%index(rptr).gt.dren) exit
    !    do while(dnode%index(dptr).lt.snode%index(rptr))
    !       dptr = dptr + 1
    !    end do
    !    row_list_sz = row_list_sz + 1
    !    row_list(row_list_sz) = dptr - dptr_sa + 1
    ! end do
    ! s2en = rptr - 1 ! Points to last row in rsrc

    call spllt_update_between_compute_map(blk, dcol, dnode, scol, snode, &
         & row_list, col_list, row_list_sz, col_list_sz, &
         & s1sa, s1en, s2sa, s2en)

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

       ! call expand_buffer(dest, n, row_list, row_list_sz, &
       !      col_list, col_list_sz, ndiag, buffer)

       call spllt_expand_buffer(dest, n, row_list, row_list_sz, &
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

  ! C wrapper
  
  subroutine spllt_update_between_c(m, n, blk_c, &
       & dcol, dnode_c, n1, scol, snode_c, dest_c, & 
       & src1_c, src2_c, buffer_c, min_width_blas) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer(c_int), value :: m ! number of rows in dest
    integer(c_int), value :: n ! number of columns in dest
    type(c_ptr), value :: blk_c ! blocks array pointer
    integer(c_int), value :: dcol
    type(c_ptr), value :: dnode_c ! blocks array pointer
    integer(c_int), value :: n1 ! number of columns in source block column
    integer(c_int), value :: scol
    type(c_ptr), value :: snode_c
    type(c_ptr), value :: dest_c ! holds block in L
    type(c_ptr), value :: src1_c, src2_c ! holds block in L
    type(c_ptr), value :: buffer_c
    integer(c_int), value :: min_width_blas

    real(wp), pointer, dimension(:) :: dest, src1, src2
    integer, dimension(:), pointer :: row_list ! reallocated to min size m
    integer, dimension(:), pointer :: col_list ! reallocated to min size n
    real(wp), pointer :: buffer(:)
    type(spllt_bc_type), pointer :: bc(:) ! blocks
    type(node_type), pointer :: dnode ! destination node
    type(node_type), pointer :: snode ! source node
    type(block_type), pointer :: blk

    call c_f_pointer(dnode_c, dnode)
    call c_f_pointer(snode_c, snode)
    call c_f_pointer(blk_c, blk)

    call c_f_pointer(dest_c, dest, (/m*n/))
    call c_f_pointer(src1_c, src1, (/n*n1/))
    call c_f_pointer(src2_c, src2, (/m*n1/))    
    
    call c_f_pointer(buffer_c, buffer, (/m*n/))    
    ! allocate(work(m*n))
    allocate(row_list(1), col_list(1))

    ! write(*,*)"node: ", blk%node
    ! write(*,*)"snode: ", snode
    ! write(*,*)'m: ', m, ', n: ', n
    ! write(*,*)'min_width_blas: ', min_width_blas

    call spllt_update_between(m, n, blk, dcol, dnode, n1, scol, snode, dest, & 
         & src1, src2, row_list, col_list, buffer, min_width_blas)

    ! deallocate(work)
    deallocate(row_list, col_list)

    return
  end subroutine spllt_update_between_c

  ! init node
  ! copy matrix coefficients into lfact array within snode
  subroutine spllt_init_node(snode, val, keep)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer, intent(in) :: snode
    real(wp), dimension(*), intent(in) :: val ! user's matrix values
! #if defined(SPLLT_USE_OMP)
!     integer, pointer, intent(inout) :: map(:)  ! mapping array. Reset for each node
! #else
!     integer, intent(inout) :: map(n)  ! mapping array. Reset for each node
! #endif
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
    ! do j = 1, size(keep%nodes(snode)%index)
    !    i = keep%nodes(snode)%index(j)
    !    map(i) = j - 1
    ! end do

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

  ! C wrapper
  subroutine spllt_init_node_c(snode, val_c, nval, keep_c) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer(c_int), value  :: snode
    type(c_ptr), value     :: val_c
    integer(c_int), value  :: nval
    type(c_ptr), value     :: keep_c
    
    real(wp), pointer :: val(:) ! user's matrix values
    type(MA87_keep), pointer :: keep 

    call c_f_pointer(val_c, val, (/nval/))
    call c_f_pointer(keep_c, keep)    

    call spllt_init_node(snode, val, keep)

    return
  end subroutine spllt_init_node_c
  
  ! init blk
  ! copy matrix coefficicents into blk
  subroutine spllt_init_blk(id, val, keep)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer(long) :: id
    real(wp), dimension(*), intent(in) :: val ! user's matrix values
    type(MA87_keep), target, intent(inout) :: keep 

    type(block_type), pointer :: blk
    integer :: sa
    integer :: sz
    integer :: bcol
    integer :: i, j

    blk => keep%blocks(id)
    
    sa = blk%sa
    sz = blk%blkn*blk%blkm
    bcol = blk%bcol

    keep%lfact(bcol)%lcol(sa:sa+sz-1) = zero

    do i = 1, keep%lmap(bcol)%len_map
       j = keep%lmap(bcol)%map(1,i)
       if ((j .ge. sa) .and. (j .le. sa+sz-1)) then
          keep%lfact(bcol)%lcol(j) = &
               val(keep%lmap(bcol)%map(2,i))
       end if
    end do
    
    return
  end subroutine spllt_init_blk

  subroutine spllt_init_blk_c(id, val_c, nval, keep_c) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer(long), value  :: id
    type(c_ptr), value    :: val_c
    integer(c_int), value        :: nval
    type(c_ptr), value    :: keep_c

    real(wp), pointer :: val(:) ! user's matrix values
    type(MA87_keep), pointer :: keep 

    call c_f_pointer(val_c, val, (/nval/))
    call c_f_pointer(keep_c, keep)
    
    call spllt_init_blk(id, val, keep)

    return
  end subroutine spllt_init_blk_c

  subroutine spllt_activate_node(snode, keep, fdata)
    use iso_c_binding
    use spllt_data_mod
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
          ! write(*,'(z16)')c_loc(keep%lfact(nbcol)%lcol(ptr))
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
