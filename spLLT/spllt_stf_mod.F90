module spllt_stf_mod
  use spllt_mod
  implicit none

contains

  subroutine spllt_stf_factorize(n, ptr, row, val, order, keep, control, info, cntl)
    use hsl_ma87_double
    use spllt_factorization_mod
#if defined(SPLLT_USE_STARPU) 
    use starpu_f_mod
#endif
    implicit none
    
    integer, intent(in) :: n ! order of A
    integer, intent(in) :: row(:) ! row indices of lower triangular part
    integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
    real(wp), intent(in) :: val(:) ! matrix values
    integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
    ! since the analyse phase)

    type(MA87_keep), target, intent(inout) :: keep 
    type(MA87_control), intent(in) :: control 
    type(MA87_info), intent(out) :: info 

    type(spllt_cntl) :: cntl

    ! local arrays
    real(wp), dimension(:), allocatable :: detlog ! per thread sum of log pivot
    ! integer, dimension(:), allocatable ::  invp ! used to hold inverse ordering
    integer, dimension(:), allocatable ::  colmap, map ! allocated to have size n.
    ! used in copying entries of user's matrix a into factor storage 
    ! (keep%fact).

    ! shortcuts
    type(node_type), pointer :: node ! node in the atree    
    type(block_type), pointer :: bc_kk, bc_ik, bc_jk, bc_ij, bc ! node in the atree    
    type(lfactor), dimension(:), pointer :: lfact! block columns in L factor    
    type(block_type), dimension(:), pointer :: blocks ! block info. 
    type(node_type), dimension(:), pointer :: nodes 
 
    ! local scalars
    integer(long) :: blk, blk1, blk2 ! block identity
    integer :: blkm, blkn ! number of rows/cols in block
    integer(long) :: dblk ! diagonal block within block column
    integer :: en ! holds keep%nodes(snode)%en
    integer :: sa ! holds keep%nodes(snode)%sa
    integer(long) :: i
    integer :: j
    integer :: s_nb ! set to block size of snode (keep%nodes(snode)%nb)
    integer :: nc, nr ! number of block column/row
    integer :: snode, num_nodes, par
    integer :: st ! stat parameter
    integer :: numrow, numcol
    integer :: total_threads ! number of threads being used
    integer :: this_thread
    integer :: ii, jj, kk
    
    ! update between variables
    integer :: csrc(2), rsrc(2) ! used for update_between tasks to
    type(node_type), pointer :: anode ! ancestor node in the atree
    ! locate source blocks
    integer :: a_num ! ancestor id
    integer :: cptr  ! Position in snode of the first row 
    ! matching a column of the current block column of anode.
    integer :: cptr2  ! Position in snode of the last row 
    ! matching a column of the current block column of anode.
    logical :: map_done
    integer :: a_sa
    integer :: ilast
    type(block_type), pointer :: a_bc ! node in the atree    
    real(wp), dimension(:), allocatable :: buffer ! update_buffer workspace
    integer, dimension(:), allocatable :: row_list, col_list ! update_buffer workspace
    integer :: cb, jb
    integer :: jlast ! Last column in the cb-th block column of anode
    ! real(wp) :: soln(0)
    integer :: k
    integer(long) :: rblk, a_dblk, a_blk ! id of block in scol containing row 
    ! nodes(snode)%index(cptr).
    integer(long) :: rb ! Index of block row in snode
    integer :: iinfo, a_nb, k1, size_anode

#if defined(SPLLT_USE_STARPU) 
    ! when using StarPU
    integer(c_int) :: ret
#endif

    ! call factorize_posdef(n, val, order, keep, control, info, 0, 0, soln)

    write(*,*) 'control%nb: ', control%nb

    ! shortcut
    blocks => keep%blocks
    nodes  => keep%nodes

    num_nodes = keep%info%num_nodes
    write(*,*) 'num_nodes: ', num_nodes

    allocate(col_list(1), row_list(1), buffer(keep%maxmn), stat=st)
    if(st.ne.0) goto 10

    ! Set up inverse permutation
    ! deallocate (invp,stat=st)
    ! allocate (invp(n),stat=st)
    ! if(st.ne.0) go to 10

    ! do j = 1, n
    !    invp(order(j)) = j
    ! end do
    
    ! Allocate factor storage (in keep%lfact)
    call spllt_init_lfact(keep)

    total_threads = 1 ! sequential
    this_thread = 0
    
    ! Allocate information arrays
    deallocate(detlog,stat=st)
    allocate(detlog(0:total_threads-1),stat=st)

    !
    ! Copy matrix values across from a into keep%lfact
    !
    allocate(map(n),stat=st)
    if(st.ne.0) go to 10
    
    ! allocate array for colomn mapping beween snode and ancestors
    allocate(colmap(n),stat=st)
    if(st.ne.0) go to 10

    call copy_a_to_l(n,num_nodes,val,map,keep)

#if defined(SPLLT_USE_STARPU)
    ! initialize starpu
    ret = starpu_f_init(cntl%ncpu)
#endif

    ! factorize nodes
    ! blk = 1
    do snode = 1, num_nodes

       node => keep%nodes(snode)
       
       write(*,*) "---------- node ----------"
       write(*,*) 'snode: ', snode 

       ! initialize node data
       
       sa = node%sa
       en = node%en
       numcol = en - sa + 1
       numrow = size(node%index)

       ! s_nb is the size of the blocks
       s_nb = node%nb

       nc = (numcol-1) / s_nb + 1 
       nr = (numrow-1) / s_nb + 1 

       write(*,'("numrow,numcol = ",i5,i5)')numrow, numcol
       write(*,'("nr,nc = ",i5,i5)')nr, nc
       
       ! init column mapping
       colmap = 0
       
       ! first block in node
       dblk = node%blk_sa

       ! Loop over the block columns in node. 
       do kk = 1, nc

          ! A_kk          
          bc_kk => keep%blocks(dblk)

          call spllt_factorize_block_task(bc_kk, keep%lfact)

          ! loop over the row blocks (that is, loop over blocks in block col)
          do ii = kk+1,nr
          ! do blk = dblk+1, dblk+sz-1
             
             ! A_mk
             blk = dblk+ii-kk
             bc_ik => keep%blocks(blk)

             call spllt_solve_block_task(bc_kk, bc_ik, keep%lfact, control)
          end do

          do jj = kk+1,nc
             
             ! L_jk
             blk2 = dblk+jj-kk
             bc_jk => keep%blocks(blk2)

             do ii = jj,nr

                ! L_ik
                blk1 = dblk+ii-kk                
                bc_ik => keep%blocks(blk1)
                
                ! A_ij
                ! blk = get_dest_block(keep%blocks(blk1), keep%blocks(blk2))
                blk = get_dest_block(keep%blocks(blk2), keep%blocks(blk1))
                bc_ij => keep%blocks(blk)
                
                call spllt_update_block_task(bc_ik, bc_jk, bc_ij, keep%lfact, control)
             end do
          end do
          
          ! map update between current block column and ancestors
          ! rsrc = 0
          ! csrc = 0
          ! call map_update_between(keep%nodes, keep%blocks, &
          !      & snode, kk, dblk, csrc, map)

          ! update between

          !  Loop over ancestors of snode
          a_num = node%parent  
          !  Initialize cptr to correspond to the first row of the rectangular part of 
          !  the snode matrix.
          cptr = 1 + numcol
          ! write(*,*)"numrow: ", numrow
          do while(a_num.gt.0)
             anode => keep%nodes(a_num) 
             
             ! Skip columns that come from other children
             do cptr = cptr, numrow
                if(node%index(cptr).ge.anode%sa) exit
             end do
             if(cptr.gt.numrow) exit ! finished with node

             map_done = .false. ! We will only build a map when we need it
             
             ! Loop over affected block columns of anode
             bcols: do
                if(cptr.gt.numrow) exit
                if(node%index(cptr).gt.anode%en) exit
                
                ! compute local index of block column in anode and find the id of 
                ! its diagonal block
                cb = (node%index(cptr) - anode%sa)/anode%nb + 1
                a_dblk = anode%blk_sa
                do jb = 2, cb
                   a_dblk = blocks(a_dblk)%last_blk + 1
                end do

                ! Find cptr2
                jlast = min(anode%sa + cb*anode%nb - 1, anode%en)
                do cptr2 = cptr,numrow
                   if(nodes(snode)%index(cptr2) > jlast) exit
                end do
                cptr2 = cptr2 - 1
                ! write(*,*)"j: ", nodes(snode)%index(cptr), ", jlast: ", nodes(snode)%index(cptr2)
                ! Set info for source block csrc (hold start location
                ! and number of entries)
                csrc(1) = 1 + (cptr-(kk-1)*node%nb-1)*blocks(bc_kk%id)%blkn
                csrc(2) = (cptr2 - cptr + 1)*blocks(bc_kk%id)%blkn
                ! write(*,*)"csrc(1): ", csrc(1), ", csrc(2): ", csrc(2)
                ! Build a map of anode's blocks if this is first for anode
                if(.not.map_done) call spllt_build_rowmap(anode, map) 
                ! write(*,*)"csrc(1): ", csrc(1), ", csrc(2): ", csrc(2)

                ! write(*,*)"map ", map
                ! if(.not.map_done) then
                !    a_nb = anode%nb
                !    size_anode = size(anode%index)
                !    map_done = .true.
                !    jb = 1
                !    do i = 1, size_anode, a_nb
                !       do k = i, min(i+a_nb-1, size_anode) ! size_anode
                !          k1 = anode%index(k)
                !          map(k1) = jb
                !       end do
                !       jb = jb + 1 
                !    end do
                ! end if

                ! Loop over the blocks of snode
                ii = map(node%index(cptr)) 
                ! ii = -1 
                ilast = cptr ! Set start of current block

                do i = cptr, numrow
                   k = map(node%index(i))
                   ! write(*,*) "i:",i," k:",k
                   blk = bc_kk%id + (i-1)/s_nb - (kk-1)
                   bc => keep%blocks(blk) ! current block in node

                   if(k.ne.ii) then
                      a_blk = a_dblk + ii - cb
                      a_bc => keep%blocks(a_blk)
                      a_sa = a_bc%sa

                      rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
                      rsrc(2) = (i - ilast)*blocks(blk)%blkn
                      ! write(*,*) "kk: ", kk, ", dblk: ", dblk, ", blk: ", blk, ", a_num: ", a_num, ", a_blk: ", a_blk
                      ! write(*,*) "ilast: ", ilast, ", i: ", i, ", rsrc(2): ", rsrc(2)
                      call spllt_update_between_task(bc, node, a_bc, anode, &
                           & csrc, rsrc, &
                           & row_list, col_list, buffer, &
                           & keep%lfact, keep%blocks, &
                           & control, st)

                      ii = k
                      ilast = i ! Update start of current block
                   end if
                end do
                ! i = min(i,numrow)
                a_blk = a_dblk + ii - cb
                a_bc => keep%blocks(a_blk)

                rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
                rsrc(2) = (i - ilast)*blocks(blk)%blkn
                ! write(*,*)"i: ", i 
                ! write(*,*)"size lcol", size(keep%lfact(keep%blocks(blk)%bcol)%lcol), ", rsrc(1)+rsrc(2)-1: ", rsrc(1)+rsrc(2)-1
                ! write(*,*) "kk: ", kk, ", dblk: ", dblk, ", blk: ", blk, ", a_num: ", a_num, ", a_blk: ", a_blk
                ! write(*,*) "ilast: ", ilast, ", i: ", i, ", rsrc(2): ", rsrc(2)
                call spllt_update_between_task(bc, node, a_bc, anode, &
                     & csrc, rsrc, &
                     & row_list, col_list, buffer, &
                     & keep%lfact, keep%blocks, &
                     & control, st)
                
                ! find id of the block in current kk colums that contains the cptr-th row
                ! rblk = (cptr-1)/s_nb - (kk-1) + bc_kk%id
                ! i = cptr
                ! do blk = rblk, blocks(rblk)%last_blk
                !    bc => keep%blocks(blk) ! current block in node
                !    rb = blk - blocks(blk)%dblk + kk ! block index                   
                !    do i = i, min(numrow, 1+rb*node%nb-1)
                !       k = map(node%index(i))
                !       ! write(*,*)"k: ", k, ", ii: ", ii            
                !       if(k.ne.ii) then
                !          a_blk = a_dblk + ii - cb
                !          write(*, *) "a_dblk: ", a_dblk
                !          write(*, *) "ii: ", ii
                !          write(*, *) "cb: ", cb
                !          write(*, *) "a_blk: ", a_blk
                !          write(*, *) "size blocks: ", size(keep%blocks)
                !          a_bc => keep%blocks(a_blk)
                !          a_sa = a_bc%sa

                !          rsrc(1) = 1 + (ilast-(kk-1)*node%nb-1)*blocks(blk)%blkn
                !          rsrc(2) = (i - ilast)*blocks(blk)%blkn

                !          ! call spllt_update_between_task(bc, node, a_bc, anode, &
                !          !      & csrc, rsrc, &
                !          !      & row_list, col_list, buffer, &
                !          !      & keep%lfact, keep%blocks, &
                !          !      & control, st)
                         
                !          call update_between(a_bc%blkm, a_bc%blkn, a_bc%id, anode, &
                !               & bc%blkn, bc%id, node, & 
                !               & keep%lfact(a_bc%bcol)%lcol(a_sa:a_sa+a_bc%blkm*a_bc%blkn-1), &
                !               & keep%lfact(bc%bcol)%lcol(csrc(1):csrc(1)+csrc(2)-1), &
                !               & keep%lfact(bc%bcol)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1), &
                !               & keep%blocks, row_list, col_list, buffer, &
                !               & control, iinfo, st)
                         
                !          ii = k
                !          ilast = i ! Update start of current block
                !       end if
                !    end do
                ! end do
                
                ! blk = blk-1
                ! bc => keep%blocks(blk) ! current block in node

                ! a_blk = a_dblk + ii - cb
                ! rsrc(1) = 1 + (ilast-(kk-1)*node%nb-1)*blocks(blk)%blkn
                ! rsrc(2) = (i - ilast)*blocks(blk)%blkn

                ! a_bc => keep%blocks(a_blk)
                ! a_sa = a_bc%sa
                ! write(*, *) "rblk: ", rblk, ", blk: ", blk, ", a_dblk: ", a_dblk
                ! write(*, *) "kk: ", kk, ", ii: ", ii
                ! write(*, *) "cb: ", cb
                ! write(*, *) "a_blk: ", a_blk
                ! write(*, *) "i: ", i, ", ilast: ", ilast, ", blkn: ", blocks(blk)%blkn
                ! write(*, *) "rsrc(1): ", rsrc(1), ", rsrc(1)+rsrc(2)-1: ", rsrc(1)+rsrc(2)-1 
                ! write(*, *) "size lcol: ", size(keep%lfact(bc%bcol)%lcol)

                ! ! call spllt_update_between_task(bc, node, a_bc, anode, &
                ! !      & csrc, rsrc, &
                ! !      & row_list, col_list, buffer, &
                ! !      & keep%lfact, keep%blocks, &
                ! !      & control, st)

                ! call update_between(a_bc%blkm, a_bc%blkn, a_bc%id, anode, &
                !      & bc%blkn, bc%id, node, & 
                !      & keep%lfact(a_bc%bcol)%lcol(a_sa:a_sa+a_bc%blkm*a_bc%blkn-1), &
                !      & keep%lfact(bc%bcol)%lcol(csrc(1):csrc(1)+csrc(2)-1), &
                !      & keep%lfact(bc%bcol)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1), &
                !      & keep%blocks, row_list, col_list, buffer, &
                !      & control, iinfo, st)
                         
                ! Move cptr down, ready for next block column of anode
                cptr = cptr2 + 1
             end do bcols
          
             a_num = anode%parent             
          end do
          
          ! move to next block column in snode
          dblk = blocks(dblk)%last_blk + 1
          ! numrow = numrow - s_nb
       end do

    end do

#if defined(SPLLT_USE_STARPU)
    call starpu_f_task_wait_for_all()
#endif


#if defined(SPLLT_USE_STARPU)
    call starpu_f_shutdown()
#endif


10 if(st.ne.0) then
      info%flag = spllt_error_allocation
      info%stat = st
      call spllt_print_err(info%flag, control, context='spllt_stf_factorize',st=st)
      return
   endif

    return
  end subroutine spllt_stf_factorize

  ! TODO comments!

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

  subroutine spllt_build_colmap(node, anode, cptr, colmap, ncolmap)
    use hsl_ma87_double
    implicit none

    type(node_type), intent(in) :: node  ! current node in the atree
    type(node_type), intent(in) :: anode  ! ancestor node in the atree
    integer, dimension(:), intent(out) :: colmap ! Workarray to hold map from row 
    integer, intent(inout) :: cptr  ! Position in snode of the first row
    ! matching a column of the current block column of anode.    
    integer, intent(inout) :: ncolmap ! number of entries in colmap

    integer :: sa, en
    integer :: numrow, numcol ! number of row/col in node
    integer :: cptr2  ! Position in snode of the last row 
    ! matching a column of the current block column of anode.
    integer :: jlast ! Last column in the cb-th block column of anode
    integer :: cb
    integer :: a_nb
    
    sa = node%sa
    en = node%en
    numcol = en - sa + 1
    numrow = size(node%index)
    
    a_nb = anode%nb

    ncolmap = 0

    ! Skip columns that come from other children
    do cptr = cptr, numrow
       if(node%index(cptr).ge.anode%sa) exit
    end do
    if(cptr.gt.numrow) return ! finished with node

    ! Loop over affected block columns of anode
    acols: do
       if(cptr.gt.numrow) exit
       if(node%index(cptr).gt.anode%en) exit

       cb = (node%index(cptr) - anode%sa)/anode%nb + 1
       ! Find cptr2
       jlast = min(anode%sa + cb*a_nb - 1, anode%en)
       do cptr2 = cptr, numrow
          if(node%index(cptr2) > jlast) exit
       end do
       cptr2 = cptr2 - 1 
       
       ncolmap = ncolmap + 1 
       colmap(ncolmap) = cptr2

       ! Move cptr down, ready for next block column of anode
       cptr = cptr2 + 1
    end do acols
    
    return
  end subroutine spllt_build_colmap

  subroutine spllt_init_lfact(keep)
    use hsl_ma87_double
    implicit none

    type(MA87_keep), target, intent(inout) :: keep 

    type(node_type), pointer :: node ! node in the atree    
    integer(long) :: blk, dblk
    integer :: nbcol, l_nb, sz, sa, en
    integer :: blkm, blkn, size_bcol
    integer :: snode, num_nodes
    integer :: i
    integer :: st ! stat parameter

    num_nodes = keep%info%num_nodes

    deallocate (keep%lfact,stat=st)
    allocate (keep%lfact(keep%nbcol),stat=st)
    ! TODO trace error

    blk = 1
    nbcol = 0
    ! loop over the nodes
    
    do snode = 1, num_nodes
       ! Loop over the block columns in snode, allocating space 
       ! l_nb is the size of the blocks and sz is number of
       ! blocks in the current block column
       
       node => keep%nodes(snode)

       l_nb = node%nb
       sz = (size(node%index) - 1) / l_nb + 1
       sa = node%sa
       en = node%en

       size_bcol = 0
       do i = sa, en, l_nb
          nbcol = nbcol + 1
          size_bcol = 0
          dblk = blk
          ! loop over the row blocks
          do blk = dblk, dblk+sz-1
             blkm = keep%blocks(blk)%blkm
             blkn = keep%blocks(blk)%blkn
             size_bcol = size_bcol + blkm*blkn
          end do
          sz = sz - 1
          allocate (keep%lfact(nbcol)%lcol(size_bcol),stat=st)
          ! TODO trace error
       end do
    end do

  end subroutine spllt_init_lfact

end module spllt_stf_mod
