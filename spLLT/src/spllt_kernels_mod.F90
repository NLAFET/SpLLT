module spllt_kernels_mod
  implicit none

contains

  ! kernel for the factorization of a node within a subtree
  subroutine spllt_subtree_factorize_node(node, blocks, lfact)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    ! type(spllt_node_type), target, intent(inout)        :: snode ! node to factorize (spllt)    
    type(node_type), target, intent(inout)        :: node ! node to factorize (ma87)    
    type(block_type), dimension(*), intent(in)          :: blocks
    type(lfactor), dimension(*), intent(inout)          :: lfact

    integer :: prio ! task priority
    integer :: sa ! first column 
    integer :: en ! last column
    integer :: numcol ! number of columns in node 
    integer :: numrow ! number of rows in node 
    integer :: nc ! number of block columns
    integer :: nr ! number of block rows
    integer(long) :: dblk ! diagonal block
    integer :: s_nb ! block size
    integer :: ii, jj, kk ! indexes 
    type(spllt_bc_type), pointer :: bc_kk, bc_ik, bc_jk, bc_ij ! block pointers
    integer(long) :: blk, blk1, blk2 ! block id

    integer :: m, n, bcol 
    integer :: m1, n1, sa1, bcol1 
    integer :: m2, n2, sa2, bcol2 
    integer :: d_m, d_n, d_sa, d_bcol 
    logical :: is_diag

    ! node => snode%node

    ! node priority
    ! num_nodes = keep%info%num_nodes
    ! snum = snode%num 
    ! prio = (num_nodes - snum + 1)*4
    prio = 0

    ! initialize node data
    sa = node%sa
    en = node%en
    numcol = en - sa + 1
    numrow = size(node%index)

    ! s_nb is the size of the blocks
    s_nb = node%nb

    nc = (numcol-1) / s_nb + 1 
    nr = (numrow-1) / s_nb + 1 

    ! first block in node
    dblk = node%blk_sa

    do kk = 1, nc
       ! A_kk

       d_m    = blocks(dblk)%blkm
       d_n    = blocks(dblk)%blkn
       bcol   = blocks(dblk)%bcol
       
       d_sa   = blocks(dblk)%sa

       ! call spllt_factorize_block_task(fdata, snode, bc_kk, keep%lfact, prio+3)
       call spllt_factor_diag_block(d_m, d_n, lfact(bcol)%lcol(d_sa:d_sa+d_n*d_m-1))

       ! loop over the row blocks (that is, loop over blocks in block col)
       do ii = kk+1,nr
          ! do ii = nr,kk+1,-1             
          ! A_mk
          blk = dblk+ii-kk
          ! bc_ik => keep%blocks(blk)
          ! bc_ik => fdata%bc(blk)
          ! write(*,*)"ii: ", ii
          ! call spllt_solve_block_task(fdata, bc_kk, bc_ik, keep%lfact,prio+2)
          m  = blocks(blk)%blkm
          n  = blocks(blk)%blkn
          sa = blocks(blk)%sa

          call spllt_solve_block(m, n, &
               & lfact(bcol)%lcol(sa:sa+n*m-1), & 
               & lfact(bcol)%lcol(d_sa:d_sa+d_n*d_m))
       end do

       ! perform update operations within node
       do jj = kk+1,nc

          ! L_jk
          blk2 = dblk+jj-kk
          ! bc_jk => fdata%bc(blk2)
          n2  = blocks(blk2)%blkn
          m2  = blocks(blk2)%blkm
          sa2 = blocks(blk2)%sa

          do ii = jj,nr

             ! L_ik
             blk1 = dblk+ii-kk                
             ! bc_ik => fdata%bc(blk1)
             n1    = blocks(blk1)%blkn
             m1    = blocks(blk1)%blkm
             sa1   = blocks(blk1)%sa
             bcol1 = blocks(blk1)%bcol

             ! A_ij
             ! blk = get_dest_block(keep%blocks(blk1), keep%blocks(blk2))
             blk = get_dest_block(blocks(blk2), blocks(blk1))
             ! bc_ij => fdata%bc(blk)
             ! call spllt_update_block_task(fdata, bc_ik, bc_jk, bc_ij, keep%lfact, prio+1)
             n  = blocks(blk)%blkn
             m  = blocks(blk)%blkm
             sa = blocks(blk)%sa

             bcol = blocks(blk)%bcol
             
             call spllt_update_block(m, n, &
                  & lfact(bcol)%lcol(sa:sa+n*m-1), &
                  & blocks(blk)%dblk.eq.blocks(blk)%id, n1, &
                  & lfact(bcol1)%lcol(sa2:sa2+n1*m2-1), &
                  & lfact(bcol1)%lcol(sa1:sa1+n1*m1-1))

          end do
       end do

       ! move to next block column in snode
       dblk = blocks(dblk)%last_blk + 1
       ! numrow = numrow - s_nb
    end do

  end subroutine spllt_subtree_factorize_node

  ! kernel for applying update on nodes within a subtree
  subroutine spllt_subtree_apply_node(node, root, nodes, blocks, lfact, buffer, &
       & map, row_list, col_list, workspace, control)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    type(node_type), intent(in) :: node 
    integer, intent(in) :: root ! root node in the subtree
    type(node_type), dimension(-1:), target, intent(in) :: nodes
    type(block_type), dimension(*), intent(in)          :: blocks
    type(lfactor), dimension(*), intent(inout)          :: lfact
    real(wp), dimension(:), allocatable, intent(inout)  :: buffer ! data array used 
    ! to accumulate updates
    integer, pointer :: map(:)
    integer, pointer :: row_list(:), col_list(:) ! worskapce used for conputing indexes 
    ! when scatering elements in update_between
    real(wp), pointer :: workspace(:) ! workspace used for update between
    type(MA87_control), intent(in) :: control


    integer :: numcol, numrow
    integer :: s_nb ! block size of source node
    integer :: sa, en
    type(node_type), pointer :: anode ! ancestor node in the atree
    ! locate source blocks
    integer :: a_num ! ancestor id
    integer :: cptr  ! Position in snode of the first row 
    ! matching a column of the current block column of anode.
    integer :: cptr2  ! Position in snode of the last row 
    ! matching a column of the current block column of anode.
    logical :: map_done
    integer :: i, ilast
    ! type(spllt_bc_type), pointer :: bc, a_bc ! node in the atree
    integer :: cb, jb
    integer :: jlast ! Last column in the cb-th block column of anode
    integer :: k
    integer(long) :: a_dblk, a_blk ! id of block in scol containing row 
    integer :: m, n
    integer(long) :: id
    integer :: csrc, csrc2, rsrc, rsrc2
    integer :: n1
    integer :: bcol, bcol1
    integer(long) :: dblk
    integer :: dcol, scol
    integer :: ii, kk ! indexes
    integer :: nc, nr

    ! varable for the applying update on buffer
    integer :: acol, arow ! column and row pointer in the buffer
    integer :: am, an ! sizes of root node
    integer :: j
    real(wp) :: v, a_ik, a_kj
    integer :: b_sz ! sizes of the buffer
    integer :: p, q

    ! update between

    sa = node%sa
    en = node%en
    numcol = en - sa + 1
    numrow = size(node%index)

    ! s_nb is the size of the blocks
    s_nb = node%nb

    nc = (numcol-1) / s_nb + 1 
    nr = (numrow-1) / s_nb + 1 

    !  Loop over ancestors of snode
    a_num = node%parent  
    !  Initialize cptr to correspond to the first row of the rectangular part of 
    !  the snode matrix.
    cptr = 1 + numcol
    ! write(*,*)"numrow: ", numrow
    do while(a_num.gt.0)
       if (a_num .gt. root) exit ! make sure we stay in the subtree
       ! if (root .ne. 21 .and. a_num .gt. root) exit ! make sure we stay in the subtree
       anode => nodes(a_num) 
       ! a_node => fdata%nodes(a_num)
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
             if(node%index(cptr2) > jlast) exit
          end do
          cptr2 = cptr2 - 1
          ! write(*,*)"j: ", nodes(snode)%index(cptr), ", jlast: ", nodes(snode)%index(cptr2)
          ! Set info for source block csrc (hold start location
          ! and number of entries)

          ! csrc(1) = 1 + (cptr-(kk-1)*node%nb-1)*blocks(bc_kk%blk%id)%blkn
          ! csrc(2) = (cptr2 - cptr + 1)*blocks(bc_kk%blk%id)%blkn

          ! write(*,*)"csrc(1): ", csrc(1), ", csrc(2): ", csrc(2)
          ! Build a map of anode's blocks if this is first for anode
          if(.not.map_done) call spllt_build_rowmap(anode, map) 
          ! write(*,*)"csrc(1): ", csrc(1), ", csrc(2): ", csrc(2)

          ! Loop over the blocks of snode
          ii = map(node%index(cptr)) 
          ! ii = -1 
          ilast = cptr ! Set start of current block

          do i = cptr, numrow
             k = map(node%index(i))

             if(k.ne.ii) then
                a_blk = a_dblk + ii - cb
                ! a_bc => keep%blocks(a_blk)
                ! a_bc => fdata%bc(a_blk)

                m  = blocks(a_blk)%blkm
                n  = blocks(a_blk)%blkn
                id = blocks(a_blk)%id
                sa = blocks(a_blk)%sa

                ! rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
                ! rsrc(2) = (i - ilast)*blocks(blk)%blkn

                dblk = node%blk_sa
                ! Loop over the block columns in node. 
                do kk = 1, nc

                   ! bc_kk => fdata%bc(dblk)

                   n1 = blocks(dblk)%blkn
                   
                   csrc  = 1 + (cptr-(kk-1)*s_nb-1)*n1
                   csrc2 = (cptr2 - cptr + 1)*n1

                   rsrc  = 1 + (ilast-(kk-1)*s_nb-1)*n1
                   rsrc2 = (i - ilast)*n1

                   bcol1 = blocks(dblk)%bcol
                   scol = bcol1 - blocks(node%blk_sa)%bcol + 1
                   
                   bcol = blocks(a_blk)%bcol
                   dcol = bcol - blocks(anode%blk_sa)%bcol + 1

                   call spllt_update_between(m, n, blocks(a_blk), dcol, anode, &
                        & n1, scol, node, &
                        & lfact(bcol)%lcol(sa:sa+m*n-1), &
                        & lfact(bcol1)%lcol(csrc:csrc+csrc2-1), &
                        & lfact(bcol1)%lcol(rsrc:rsrc+rsrc2-1), &
                        & row_list, col_list, workspace, &
                        & control%min_width_blas)

                   dblk = blocks(dblk)%last_blk + 1
                end do

                ii = k
                ilast = i ! Update start of current block
             end if
          end do

          ! i = min(i,numrow)
          a_blk = a_dblk + ii - cb
          ! a_bc => keep%blocks(a_blk)
          ! a_bc => fdata%bc(a_blk
          m  = blocks(a_blk)%blkm
          n  = blocks(a_blk)%blkn
          id = blocks(a_blk)%id
          sa = blocks(a_blk)%sa

          ! rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
          ! rsrc(2) = (i - ilast)*blocks(blk)%blkn

          dblk = node%blk_sa
          ! Loop over the block columns in node. 
          do kk = 1, nc

             ! bc_kk => fdata%bc(dblk)

             n1 = blocks(dblk)%blkn

             csrc  = 1 + (cptr-(kk-1)*s_nb-1)*n1
             csrc2 = (cptr2 - cptr + 1)*n1

             rsrc  = 1 + (ilast-(kk-1)*s_nb-1)*n1
             rsrc2 = (i - ilast)*n1

             bcol1 = blocks(dblk)%bcol
             scol = bcol1 - blocks(node%blk_sa)%bcol + 1

             bcol = blocks(a_blk)%bcol
             dcol = bcol - blocks(anode%blk_sa)%bcol + 1

             call spllt_update_between(m, n, blocks(a_blk), dcol, anode, &
                  & n1, scol, node, &
                  & lfact(bcol)%lcol(sa:sa+m*n-1), &
                  & lfact(bcol1)%lcol(csrc:csrc+csrc2-1), &
                  & lfact(bcol1)%lcol(rsrc:rsrc+rsrc2-1), &
                  & row_list, col_list, workspace, &
                  & control%min_width_blas)

             dblk = blocks(dblk)%last_blk + 1
          end do

          ! Move cptr down, ready for next block column of anode
          cptr = cptr2 + 1
       end do bcols

       a_num = anode%parent
    end do
    

    ! print *, "root: ", root,", cptr: ", cptr, ", numrow: ", numrow, ", numcol: ", numcol, ", s_nb: ", s_nb
    ! print *, "nc: ", nc

    ! accumulate update in buffer for nodes above root
    anode => nodes(root)
    acol = 1
    arow = 1
    am = size(anode%index)
    an = anode%en - anode%sa + 1
    b_sz = am-an

    ! print *, "am: ", am, ", an: ", an, ", b_sz: ", b_sz

    map_done = .false. ! We will only build a map when we need it
    buff_bcols: do
       if(cptr.gt.numrow) exit

       ! compute local index of block column in anode and find the id of 
       ! its diagonal block
       cb = (node%index(cptr) - anode%sa)/anode%nb + 1
       print *, "cptr: ", cptr, "node%index(cptr): ", node%index(cptr), ", anode%sa: ", anode%sa 
       print *, "cb: ", cb, ", a_nb: ", anode%nb
       ! Find cptr2
       jlast = anode%sa + cb*anode%nb - 1
       do cptr2 = cptr,numrow
          if(node%index(cptr2) > jlast) exit
       end do
       cptr2 = cptr2 - 1

       m = cptr2 - cptr + 1 ! number of columns in update

       if(.not.map_done) call spllt_build_rowmap(anode, map) 

       ! Loop over the blocks of snode
       ii = map(node%index(cptr)) 
       ! ii = -1 
       ilast = cptr ! Set start of current block

       do i = cptr, numrow
          k = map(node%index(i))

          if(k.ne.ii) then

             ! n = i - ilast

             ! dblk = node%blk_sa
             ! ! Loop over block columns in node.
             ! ! accumulate update in workspace
             ! do kk = 1, nc

             !    n1 = blocks(dblk)%blkn

             !    csrc  = 1 + (cptr-(kk-1)*s_nb-1)*n1
             !    csrc2 = m*n1

             !    rsrc  = 1 + (ilast-(kk-1)*s_nb-1)*n1
             !    rsrc2 = n*n1

             !    bcol1 = blocks(dblk)%bcol

             !    ! print *, "kk: ",kk
             !    ! print *, "cptr: ", cptr, ", cptr2: ", cptr2
             !    ! print *, "ilast: ", ilast, ", i: ", i
             !    ! print *, "m: ", cptr2-cptr+1, ", n: ", i-ilast, ", n1: ", n1

             !    ! workspace += a_ik * a_kj
             !    call dgemm('T', 'N', int(m), int(n), n1, &
             !         & one, &
             !         & lfact(bcol1)%lcol(csrc:csrc+csrc2-1), n1, &
             !         & lfact(bcol1)%lcol(rsrc:rsrc+rsrc2-1), n1, &
             !         & one, workspace, int(m))

             ! end do

             ! ! expand workspace into buffer
             ! arow = 1
             ! acol = 1
             ! do p = cptr, cptr2
             !    ! column to update in buffer 
             !    do while(anode%index(acol).ne.node%index(p))
             !       acol = acol + 1
             !    end do

             !    arow = 1
             !    do q = ilast, i-1 
             !       ! rows to update in buffer
             !       do while(anode%index(arow).ne.node%index(q))
             !          arow = arow + 1
             !       end do

             !       ! print *, "arow: ", arow, "acol: ", acol
             !       ! print *, "workspace: ", workspace((q-ilast)*n + (p-cptr+1))
             !       buffer((arow-an-1)*b_sz + (acol-an)) = &
             !            buffer((arow-an-1)*b_sz + (acol-an)) + workspace((q-ilast)*n + (p-cptr+1))

             !    end do

             ! end do

             ! ii = k
             ! ilast = i ! Update start of current block
          end if
       end do

       n = i - ilast
       ! initialize workspace
       workspace = 0

       dblk = node%blk_sa
       ! Loop over the block columns in node. 
       do kk = 1, nc
          
          n1 = blocks(dblk)%blkn
          
          csrc  = 1 + (cptr-(kk-1)*s_nb-1)*n1
          csrc2 = m*n1

          rsrc  = 1 + (ilast-(kk-1)*s_nb-1)*n1
          rsrc2 = n*n1

          bcol1 = blocks(dblk)%bcol

          ! print *, "kk: ",kk
          ! print *, "cptr: ", cptr, ", cptr2: ", cptr2
          ! print *, "ilast: ", ilast, ", i: ", i
          ! print *, "m: ", cptr2-cptr+1, ", n: ", i-ilast, ", n1: ", n1
          
          ! workspace += a_ik * a_kj
          call dgemm('T', 'N', int(m), int(n), n1, &
               & one, &
               & lfact(bcol1)%lcol(csrc:csrc+csrc2-1), n1, &
               & lfact(bcol1)%lcol(rsrc:rsrc+rsrc2-1), n1, &
               & one, workspace, int(m))

       end do
       
       ! expand workspace into buffer
       arow = 1
       acol = 1
       do p = cptr, cptr2
          ! column to update in buffer 
          do while(anode%index(acol).ne.node%index(p))
             acol = acol + 1
          end do
          
          arow = 1
          do q = ilast, i-1 
             ! rows to update in buffer
             do while(anode%index(arow).ne.node%index(q))
                arow = arow + 1
             end do

             ! print *, "arow: ", arow, "acol: ", acol
             ! print *, "workspace: ", workspace((q-ilast)*n + (p-cptr+1))
             buffer((arow-an-1)*b_sz + (acol-an)) = &
                  buffer((arow-an-1)*b_sz + (acol-an)) + workspace((q-ilast)*n + (p-cptr+1))
             
          end do
          
       end do

       ! Move cptr down, ready for next block column of anode
       cptr = cptr2 + 1       
    end do buff_bcols

    ! ! start from curent position in node
    ! do i = cptr, numrow
    !    ! Find col of generated element
    !    do while(nodes(root)%index(acol) .ne. node%index(i))
    !       acol = acol + 1
    !    end do

    !    arow = acol
    !    do j = i, numrow ! rows to update
    !       ! Find row of ancestor
    !       do while(nodes(root)%index(arow) .ne. node%index(j))
    !          arow = arow + 1
    !       end do

    !       ! Do update a_ij -= sum_p a_ik a_kj
    !       v    = 0.0
    !       a_ik = 0.0
    !       a_kj = 0.0

    !       ! udpate buffer : non blocked version
    !       bcol = blocks(node%blk_sa)%bcol

    !       ! n = numcol
    !       ! do k = 1, n
    !       !    a_ik = lfact(bcol)%lcol((i-1)*n+k)
    !       !    a_kj = lfact(bcol)%lcol((j-1)*n+k)
    !       !    v = v + a_ik * a_kj
    !       ! end do
    !       ! buffer((arow-an-1)*b_sz + (acol-an)) = &
    !       !      buffer((arow-an-1)*b_sz + (acol-an)) + v

    !       ! udpate buffer : blocked version
    !       dblk = node%blk_sa
    !       ! Loop over the block columns in node. 
    !       do kk = 1, nc
             
    !          bcol = blocks(dblk)%bcol
    !          n = blocks(dblk)%blkn

    !          do k = 1, n
                
    !             csrc = (i-(kk-1)*s_nb-1)*n + k
    !             rsrc = (j-(kk-1)*s_nb-1)*n + k

    !             a_ik = lfact(bcol)%lcol(csrc)
    !             a_kj = lfact(bcol)%lcol(rsrc)
    !             v = v + a_ik * a_kj
    !          end do

    !          dblk = blocks(dblk)%last_blk + 1
    !       end do
    !       buffer((arow-an-1) * b_sz + (acol-an)) = &
    !            buffer((arow-an-1) * b_sz + (acol-an)) + v
          
    !    end do
    ! end do

  end subroutine spllt_subtree_apply_node

  subroutine spllt_subtree_factorize(root, keep, buffer, control)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none
    
    integer, intent(in) :: root ! root of subtree
    type(MA87_keep), target, intent(inout) :: keep ! on exit, matrix a copied
    real(wp), dimension(:), allocatable, intent(inout) :: buffer
    type(MA87_control), intent(in) :: control
    ! integer, intent(out) :: flag ! TODO error managment

    integer :: node, m, n
    type(node_type), pointer :: snode

    ! workspaces
    integer, pointer :: row_list(:), col_list(:)
    real(wp), pointer :: workspace(:)
    integer, pointer :: map(:)

    m = size(keep%nodes(root)%index)
    n = keep%nodes(root)%en - keep%nodes(root)%sa + 1
    buffer(1:(m-n)**2) = 0.0

    ! allocate workspaces
    ! TODO use workspaces that are initialized in init routine
    allocate(row_list(keep%maxmn))
    allocate(col_list(keep%maxmn))
    allocate(workspace(keep%maxmn*keep%maxmn))
    allocate(map(keep%n))

    ! Loop over nodes of tree in order
    do node = keep%nodes(root)%least_desc, root
       
       snode => keep%nodes(node)
       ! factorize supernode
       call spllt_subtree_factorize_node(snode, keep%blocks, keep%lfact)
       ! apply udpate on ancestor node (right-looking update)
       call spllt_subtree_apply_node(snode, root, keep%nodes, keep%blocks, keep%lfact, buffer, &
            & map, row_list, col_list, workspace, control)
    end do

    ! deallocate workspaces
    deallocate(row_list, col_list)
    deallocate(workspace)
    deallocate(map)

  end subroutine spllt_subtree_factorize
  
  subroutine spllt_factor_subtree(root, nodes, blocks, lfact, buffer)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none
    
    integer, intent(in) :: root ! root of subtree
    type(node_type), dimension(-1:), intent(in) :: nodes
    type(block_type), dimension(*), intent(in) :: blocks
    type(lfactor), dimension(*), intent(inout) :: lfact
    real(wp), dimension(:), allocatable, intent(inout) :: buffer
    ! integer, intent(out) :: flag ! TODO error managment

    integer(long) :: blk, sa, id
    integer :: node, m, n, bcol, abcol, am, an, anc, arow, acol, i, j, k, p
    real(wp) :: a_ip, a_kp, v

    !   print *, "ENTRY", nodes(root)%least_desc, root
    !   do node = nodes(root)%least_desc, root
    !      blk = nodes(node)%blk_sa
    !      m = blocks(blk)%blkm
    !      n = blocks(blk)%blkn
    !      sa = blocks(blk)%sa
    !      bcol = blocks(blk)%bcol
    !      id = blocks(blk)%id
    !
    !      print *, "node ", node, nodes(node)%sa, nodes(node)%en
    !      do i = 1, m
    !         print "(i4,a,10es10.2)", nodes(node)%index(i), ":", lfact(bcol)%lcol(sa+(i-1)*n:sa+i*n-1)
    !      end do
    !   end do

    ! print *, "factor_subtree, root: ", root
    
    ! Ensure buffer is sufficiently large, then zero out
    blk = nodes(root)%blk_sa
    m = blocks(blk)%blkm
    n = blocks(blk)%blkn
    buffer(1:(m-n)**2) = 0.0

    ! Loop over nodes of tree in order
    do node = nodes(root)%least_desc, root
       ! Gather stats
       blk = nodes(node)%blk_sa
       m = blocks(blk)%blkm
       n = blocks(blk)%blkn
       bcol = blocks(blk)%bcol
       sa = blocks(blk)%sa
       id = blocks(blk)%id
       ! print *, "node: ", node, ", blk: ", blk
       ! Factor node
       ! print *, "node: ", node, ", m: ", m, ", n: ", n
       call spllt_factor_diag_block(m, n, lfact(bcol)%lcol(sa:sa+n*m-1))

       ! Apply updates
       ! First, within subtree
       anc = nodes(node)%parent
       if(anc.le.root .and. anc.ne.-1) then
          acol = 1
          abcol = blocks(nodes(anc)%blk_sa)%bcol
          an = blocks(nodes(anc)%blk_sa)%blkn
       endif
       do i = n+1, m ! rows below diagonal block (=> col to update)
          j = nodes(node)%index(i)
          ! Find ancestor
          do while(j.gt.nodes(anc)%en)
             anc = nodes(anc)%parent
             if(anc.eq.-1) exit
             acol = 1
             abcol = blocks(nodes(anc)%blk_sa)%bcol
             an = blocks(nodes(anc)%blk_sa)%blkn
          end do
          if(anc.gt.root .or. anc.eq.-1) exit ! Not in this subtree
          ! Find col of ancestor
          do while(nodes(anc)%index(acol).ne.j)
             acol = acol + 1
          end do
          arow = acol
          do k = i, m ! rows to update
             ! Find row of ancestor
             do while(nodes(anc)%index(arow).ne.nodes(node)%index(k))
                arow = arow + 1
             end do
             ! Do update a_ki -= sum_p a_ip a_kp
             v = 0.0
             do p = 1, n
                a_ip = lfact(bcol)%lcol((i-1)*n+p)
                a_kp = lfact(bcol)%lcol((k-1)*n+p)
                v = v + a_ip * a_kp
             end do
             lfact(abcol)%lcol((arow-1)*an+acol) = &
                  lfact(abcol)%lcol((arow-1)*an+acol) - v
          end do
       end do
       ! Second, into buffer corresponding to root's generated element
       anc = root
       acol = 1
       abcol = blocks(nodes(root)%blk_sa)%bcol
       am = blocks(nodes(root)%blk_sa)%blkm
       an = blocks(nodes(root)%blk_sa)%blkn
       do i = i, m ! Remaining rows below diagonal (=> cols of gen. element)
          j = nodes(node)%index(i)
          ! Find col of generated element
          do while(nodes(anc)%index(acol).ne.j)
             acol = acol + 1
          end do
          arow = acol
          do k = i, m ! rows to update
             ! Find row of ancestor
             do while(nodes(anc)%index(arow).ne.nodes(node)%index(k))
                arow = arow + 1
             end do
             ! Do update a_ki -= sum_p a_ip a_kp
             v = 0.0
             do p = 1, n
                a_ip = lfact(bcol)%lcol((i-1)*n+p)
                a_kp = lfact(bcol)%lcol((k-1)*n+p)
                v = v + a_ip * a_kp
             end do
             buffer((arow-an-1)*(am-an) + (acol-an)) = &
                  buffer((arow-an-1)*(am-an) + (acol-an)) + v
          end do
       end do
    end do
    ! print *, "Finally, buffer = "
    ! do i = 1, m-n
    !   print "(i3,a,10es10.2)", i, ":", buffer((i-1)*(m-n)+1:i*(m-n))
    ! end do
  end subroutine spllt_factor_subtree

  ! Apply updates from subtree's generated element to its ancestors
  ! This routine aquires all locks as needed, and decrements dependencies
  ! upon release of relevant locks.
  ! NB: This is basically an immediate acting variant of add_between_updates()
  subroutine spllt_apply_subtree(root, buffer, nodes, blocks, lfact, map)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer, intent(in) :: root ! root of subtree
    real(wp), dimension(*), intent(in) :: buffer ! generated element
    type(node_type), dimension(-1:), intent(in) :: nodes
    type(block_type), dimension(*), intent(inout) :: blocks
    type(lfactor), dimension(*), intent(inout) :: lfact
    integer, dimension(:), intent(inout) :: map ! Workarray to hold map from row
    ! indices to block indices in ancestor node. 

    integer :: dest ! target block
    integer :: rsrc(2), csrc(2) ! specifices block of gen element to act on

    integer :: a_nb  ! Block size of anode
    integer :: anode ! Ancestor of snode
    integer(long) :: blk   ! Block id
    integer :: bsa, ben
    integer :: cb    ! Local index of column block in anode
    integer :: cptr  ! Position in snode of the first row 
    ! matching a column of the current block column of anode.
    integer :: cptr2  ! Position in snode of the last row 
    ! matching a column of the current block column of anode.
    integer(long) :: dblk ! id of diagonal block of anode
    integer :: i
    integer :: ilast
    integer :: jb ! Block index in anode
    integer :: jlast ! Last column in the cb-th block column of anode
    integer :: k
    integer :: k1
    integer :: lds ! leading dimension (row width) of generated element buffer
    logical :: map_done ! True if map has been built for anode.
    integer :: m ! set to blocks(src)%blkm
    integer :: n ! set to blocks(src)%blkn
    integer :: numcol ! number of cols in snode
    integer :: rb ! Index of block row in snode
    integer :: size_anode ! size(nodes(anode)%index)
    integer :: size_snode ! size(nodes(snode)%index)
    integer :: s_nb   ! Block size of snode

    ! cache some values in variables
    size_snode = size(nodes(root)%index)
    s_nb = nodes(root)%nb
    m = size(nodes(root)%index)
    n = nodes(root)%en - nodes(root)%sa + 1
    ! m = blocks(nodes(root)%blk_sa)%blkm
    ! n = blocks(nodes(root)%blk_sa)%blkn
    lds = size_snode-n

    if(m-n.eq.0) return ! was a root already

    ! Loop over ancestors of subtree root
    anode = nodes(root)%parent
    ! Initialize cptr to correspond to the first row of the rectangular part of
    ! the snode matrix.
    numcol = nodes(root)%en - nodes(root)%sa + 1
    cptr = 1 + numcol

    do while(anode.gt.0)
       ! Skip columns that come from other children
       do cptr = cptr, size_snode
          if(nodes(root)%index(cptr).ge.nodes(anode)%sa) exit
       end do
       if(cptr.gt.size_snode) exit ! finished with snode

       map_done = .false. ! We will only build a map when we need it
       a_nb = nodes(anode)%nb

       ! Loop over affected block columns of anode
       bcols: do
          if(cptr.gt.size_snode) exit
          if(nodes(root)%index(cptr).gt.nodes(anode)%en) exit

          ! compute local index of block column in anode and find the id of 
          ! its diagonal block
          cb = (nodes(root)%index(cptr) - nodes(anode)%sa)/a_nb + 1
          dblk = nodes(anode)%blk_sa
          do jb = 2, cb
             dblk = blocks(dblk)%last_blk + 1
          end do

          ! Find cptr2
          jlast = min(nodes(anode)%sa + cb*a_nb - 1, nodes(anode)%en)
          do cptr2 = cptr, size_snode
             if(nodes(root)%index(cptr2) > jlast) exit
          end do
          cptr2 = cptr2 - 1 

          ! Set info for source block csrc (hold start and end locations)
          csrc(1) = cptr
          csrc(2) = cptr2

          ! Build a map of anode's blocks if this is first for anode
          if(.not.map_done) then
             ! The indices for each row block in anode are mapped to a local row
             ! block index.
             size_anode = size(nodes(anode)%index)
             map_done = .true.
             jb = 1
             do i = 1, size_anode, a_nb
                do k = i, size_anode
                   k1 = nodes(anode)%index(k)
                   map(k1) = jb
                end do
                jb = jb + 1 
             end do
          endif

          ! Loop over the blocks of snode
          jb = map(nodes(root)%index(cptr))
          i = cptr
          ilast = i ! Set start of current block
          blk = nodes(root)%blk_sa
          rb = 1 ! block index 
          do i = i, size_snode
             k = map(nodes(root)%index(i))
             if(k.ne.jb) then
                ! Moved to a new block in anode
                dest = dblk + jb - cb
                ! Set info for source block rsrc (hold start and end locations)
                rsrc(1) = ilast
                rsrc(2) = i-1
                bsa = (rsrc(1)-n-1)*lds+csrc(1)-n
                ben = (rsrc(2)-n-1)*lds+csrc(2)-n

                call spllt_scatter_block(rsrc(2)-rsrc(1)+1, csrc(2)-csrc(1)+1, &
                     nodes(root)%index(rsrc(1):rsrc(2)), &
                     nodes(root)%index(csrc(1):csrc(2)), &
                     buffer(bsa:ben), lds, &
                     nodes(anode)%index((jb-1)*a_nb+1), &
                     nodes(anode)%index((cb-1)*a_nb+1), &
                     lfact(blocks(dest)%bcol)%lcol(blocks(dest)%sa:), &
                     blocks(dest)%blkn)
                jb = k
                ilast = i ! Update start of current block
             endif
          end do
          dest = dblk+jb-cb
          rsrc(1) = ilast
          rsrc(2) = i-1
          bsa = (rsrc(1)-n-1)*lds+csrc(1)-n
          ben = (rsrc(2)-n-1)*lds+csrc(2)-n

          call spllt_scatter_block(rsrc(2)-rsrc(1)+1, csrc(2)-csrc(1)+1, &
               nodes(root)%index(rsrc(1):rsrc(2)), &
               nodes(root)%index(csrc(1):csrc(2)), &
               buffer(bsa:ben), lds, &
               nodes(anode)%index((jb-1)*a_nb+1), &
               nodes(anode)%index((cb-1)*a_nb+1), &
               lfact(blocks(dest)%bcol)%lcol(blocks(dest)%sa:), &
               blocks(dest)%blkn)

          ! Move cptr down, ready for next block column of anode
          cptr = cptr2 + 1
       end do bcols

       ! Move up the tree
       anode = nodes(anode)%parent
    end do
  end subroutine spllt_apply_subtree

  ! Scatters buffer src across dest: dest <= dest - src
  subroutine spllt_scatter_block(m, n, rsrc_index, csrc_index, src, lds, &
       rdest_index, cdest_index, dest, ldd)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none

    integer, intent(in) :: m ! number of rows in source buffer
    integer, intent(in) :: n ! number of cols in source buffer
    integer, dimension(n), intent(in) :: rsrc_index ! source col indices
    integer, dimension(n), intent(in) :: csrc_index ! source col indices
    real(wp), dimension(*), intent(in) :: src ! source buffer
    integer, intent(in) :: lds ! leading dimension (row major) or source buffer
    integer, dimension(n), intent(in) :: rdest_index ! dest col indices
    integer, dimension(n), intent(in) :: cdest_index ! dest col indices
    real(wp), dimension(*), intent(inout) :: dest ! dest buffer
    integer, intent(in) :: ldd ! leading dimension (row major) or dest buffer

    integer :: sr, sc, srow, scol
    integer :: dr, dc

    dr = 1
    do sr = 1, m
       srow = rsrc_index(sr)
       do while(rdest_index(dr).ne.srow)
          dr = dr + 1
       end do
       dc = 1
       do sc = 1, n
          scol = csrc_index(sc)
          do while(cdest_index(dc).ne.scol)
             dc = dc + 1
          end do
          dest((dr-1)*ldd+dc) = dest((dr-1)*ldd+dc) - src((sr-1)*lds+sc)
       end do
    end do
  end subroutine spllt_scatter_block


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
