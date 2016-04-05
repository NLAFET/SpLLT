module spllt_analyse_mod
  use spllt_mod
  implicit none

contains

  ! analysis phase
  subroutine spllt_analyse(adata, n, ptr, row, order, keep, cntl, info)
    use hsl_mc78_integer
    use hsl_mc34_double
    use hsl_ma87_double
    implicit none

    type(spllt_adata_type) :: adata
    integer, intent(in) :: n ! order of A
    integer, intent(in) :: row(:) ! row indices of lower triangular part
    integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
    integer, intent(inout), dimension(:) :: order
    ! order(i) must hold position of i in the pivot sequence. 
    ! On exit, holds the pivot order to be used by MA87_factor.
    type(MA87_keep), intent(out) :: keep
    type(spllt_cntl), intent(in) :: cntl
    type(MA87_info), intent(inout) :: info

    integer, allocatable :: amap(:) ! map from user a to reordered a
    integer, allocatable :: aptr(:) ! column pointers of expanded matrix
    integer, allocatable :: arow(:) ! row pointers of expanded matrix
    integer, allocatable :: iw(:) ! work array
    integer, allocatable :: perm(:) ! inverse permutation.
    ! perm(i) holds the variable that is i-th in the pivot sequence.
    ! Also used for checking user-supplied permutation.
    integer, allocatable :: map(:) ! Allocated to have size n.
    ! used in computing dependency counts. For each row k in 
    ! j-th block column of a node, map(k1) is set to j
    
    integer :: num_nodes ! number of nodes in tree
    integer :: node ! a node in tree
    integer :: l_nb ! set to block size of snode (keep%nodes(snode)%nb)
    integer :: ne ! set to ptr(n+1) - 1
    integer :: nemin ! min. number of eliminations
    integer :: st ! stat parameter
    integer :: i, j, k, jj, k1 ! temporary variable
    integer(long) :: ii ! temporary variable
    integer :: sz ! number of blocks in a block column of node
    integer :: sa ! holds keep%nodes(snode)%sa
    integer :: en ! holds keep%nodes(snode)%en
    integer :: numcol ! number of cols in node (en-sa+1)
    integer :: numrow ! number of rows in a block column
    integer :: row_used ! used in computing number of rows in a block.
    integer :: col_used ! used in computing number of cols in block col.
    integer :: cb ! index of block column within node
    integer :: anode ! ancestoral node of snode in tree
    integer(long) :: blk ! temporary variable for holding block number
    integer(long) :: dblk ! diagonal block within block column
    integer :: blkn ! number of columns within block column
    integer :: ci ! do loop variable. current block col.

    type(mc78_control) :: control78
    integer :: par
    integer :: info78
    integer, dimension(:), allocatable :: sptr, sparent, rlist
    integer(long), dimension(:), allocatable :: rptr

    ! initialise
    info%flag = 0
    info%num_factor = 0_long
    info%num_flops = 0_long
    info%num_nodes = 0
    info%maxdepth = 0
    info%stat = 0

    keep%n = n
    ne = ptr(n+1) - 1

    ! immediate return if n = 0
    if (n == 0) return

    ! expand the matrix
    
    ! allocate space for expanded matrix (held in aptr,arow)
    allocate (arow(2*ne-n),aptr(n+3),iw(n+1),amap(ptr(n+1)-1),stat=st)
    if (st /= 0) go to 9999

    arow(1:ne) = row(1:ne)
    aptr(1:n+1) = ptr(1:n+1)
    call mc34_expand(n, arow, aptr, iw)
    deallocate(iw,stat=st)

    nemin = cntl%nemin
    ! Check nemin (a node is merged with its parent if both involve
    ! fewer than nemin eliminations). If out of range, use the default
    if (nemin < 1) nemin = spllt_nemin_default

    ! Check the user-supplied array order and set the inverse in perm.
    if (size(order).lt.n) then
       go to 9999
    end if

    deallocate (perm,stat=st)
    allocate (perm(n),stat=st)
    if (st /= 0) go to 9999
    perm(:) = 0
    k1 = 0
    do i = 1, n
       jj = order(i)
       if (jj < 1 .or. jj > n) exit
       if (perm(jj) /= 0) exit ! Duplicate found
       perm(jj) = i
    end do
    if (i-1 /= n) then
       go to 9999
    end if
    
    control78%nemin = nemin
    control78%sort = .true.
    control78%lopt = .true.
    call mc78_analyse(n, aptr, arow, order, num_nodes, &
         sptr, sparent, rptr, rlist, control78, info78, nfact=adata%num_factor, &
         nflops=adata%num_flops)


    info%num_nodes = num_nodes
    !**************************************
    ! Set up nodal data structures
    ! For each node, hold data in keep%nodes(node) 
    deallocate(keep%nodes, stat=st)

    allocate(keep%nodes(-1:num_nodes),stat=st)
    if (st /= 0) go to 9999

    keep%nodes(0)%blk_en = 0
    keep%nodes(1)%blk_sa = 1
    keep%nodes(1)%sa = 1

    ! loop over root nodes
    keep%nodes(:)%nchild = 0
    do node = 1, num_nodes
       keep%nodes(node)%sa = sptr(node)
       keep%nodes(node)%en = sptr(node+1)-1

       par = sparent(node)
       keep%nodes(node)%parent = par
       if(par .le. num_nodes) then
          keep%nodes(par)%nchild = keep%nodes(par)%nchild + 1
       else
          keep%nodes(node)%parent = -1
       endif

       ! determine and record the block size for node
       ! note we are careful in case l_nb**2 overflows (in fact 1+l_nb must
       ! not overflow at the end), and limit the answer to huge(l_nb)/2
       l_nb = cntl%nb
       if (l_nb < 1) l_nb = spllt_nb_default
       l_nb = min(huge(l_nb)/2_long, &
            (l_nb**2_long) / min(sptr(node+1)-sptr(node), l_nb) )
       l_nb = (l_nb-1) / 8 + 1
       l_nb = 8 * l_nb
       keep%nodes(node)%nb = l_nb

       ! Copy row list into keep%nodes
       allocate(keep%nodes(node)%index(rptr(node+1)-rptr(node)),stat=st)
       if (st /= 0) go to 9999
       j = 1
       do ii = rptr(node), rptr(node+1)-1
          keep%nodes(node)%index(j) = rlist(ii)
          j = j + 1
       end do

       ! Allocate space to store child nodes
       allocate(keep%nodes(node)%child(keep%nodes(node)%nchild), stat=st)
       if(st.ne.0) goto 9999

       ! Calculate number j of blocks in node and set
       ! keep%nodes(node)%blk_en
       sz = int(rptr(node+1)-rptr(node) - 1) / l_nb + 1
       j = 0
       do i = keep%nodes(node)%sa, keep%nodes(node)%en, l_nb
          j = j + sz
          sz = sz - 1
       end do
       keep%nodes(node)%blk_en = keep%nodes(node-1)%blk_en + j

       ! if node is not the final node, record first block
       ! for the next node (node+1)
       if (node < num_nodes)  &
            keep%nodes(node+1)%blk_sa = keep%nodes(node)%blk_en + 1

       ! if(cntl%prune_tree) then
       !    keep%nodes(node)%subtree_root = subtree_root(node)
       ! else
       !    keep%nodes(node)%subtree_root = -1
       ! endif
    end do

    ! set keep%final_blk to hold total number of blocks.
    keep%final_blk = keep%nodes(num_nodes)%blk_en

    ! Add children to nodes, use sptr as a counter as it has fufilled its purpose
    sptr(:) = 0
    do node = 1, num_nodes
       par = sparent(node)
       if(par.gt.num_nodes) cycle
       sptr(par) = sptr(par) + 1
       keep%nodes(par)%child(sptr(par)) = node
    end do

    ! Setup least descendants, to allow easy walk of subtrees
    do node = -1, num_nodes
       ! initialise least descendat to self
       keep%nodes(node)%least_desc = node
    end do
    do node = 1, num_nodes
       ! walk up tree from leaves. A parent's least descendant is either this
       ! nodes least descendant (itself in case of a leaf), or its existing
       ! one if that is smaller.
       anode = keep%nodes(node)%parent
       keep%nodes(anode)%least_desc = &
            min(keep%nodes(node)%least_desc, keep%nodes(anode)%least_desc)
    end do

    !**************************************   
    ! Fill out block information. 
    
    deallocate(keep%blocks,stat=st)
    allocate(keep%blocks(keep%final_blk),stat=st)
    if(st.ne.0) go to 9999

    ! Loop over the nodes. Number the blocks in the first node
    ! contiguously, then those in the second node, and so on.
    ! Each node has a number of block columns; the blocks within
    ! each block column are numbered contiguously.
    ! Also add up the number of block columns and store largest block dimension.
    blk = 1
    keep%nbcol = 0
    keep%maxmn = 0
    do node = 1, num_nodes

       sa = keep%nodes(node)%sa
       en = keep%nodes(node)%en
       numcol = en - sa + 1
       numrow = size(keep%nodes(node)%index)

       ! l_nb is the size of the blocks
       l_nb = keep%nodes(node)%nb

       ! sz is number of blocks in the current block column
       sz = (numrow - 1) / l_nb + 1

       ! cb is the index of the block col. within node
       cb = 0
       col_used = 0

       ! Loop over the block columns in node. 
       do ci = sa, en, l_nb
          k = 1 ! use k to hold position of block within block column
          ! increment count of block columns
          keep%nbcol = keep%nbcol + 1

          cb = cb + 1

          ! blkn is the number of columns in the block column.
          ! For all but the last block col. blkn = l_nb.
          blkn = min(l_nb, numcol-col_used)
          col_used = col_used + blkn

          dblk = blk

          ! loop over the row blocks (that is, loop over blocks in block col)
          row_used = 0 
          do blk = dblk, dblk+sz-1
             ! store identity of block
             keep%blocks(blk)%id       = blk
             !print *, "node ", node, "=> blk", blk

             ! store number of rows in the block.
             ! For all but the last block, the number of rows is l_nb
             keep%blocks(blk)%blkm     = min(l_nb, numrow-row_used)
             row_used = row_used + keep%blocks(blk)%blkm

             ! store number of columns in the block.
             keep%blocks(blk)%blkn     = blkn

             keep%maxmn = max(keep%maxmn, &
                  keep%blocks(blk)%blkm,  keep%blocks(blk)%blkn)

             ! store position of the first entry of the block within the
             ! block column of L
             keep%blocks(blk)%sa       = k

             ! store identity of diagonal block within current block column
             keep%blocks(blk)%dblk     = dblk

             ! store identity of last block within current block column
             keep%blocks(blk)%last_blk = dblk + sz - 1

             ! store node the blk belongs to
             keep%blocks(blk)%node     = node

             ! initialise dependency count
             keep%blocks(blk)%dep_initial = cb

             ! store identity of block column that blk belongs to
             keep%blocks(blk)%bcol     = keep%nbcol

             !$          call omp_init_lock(keep%blocks(blk)%lock)
             !$          call omp_init_lock(keep%blocks(blk)%alock)

             ! increment k by number of entries in block
             k = k + keep%blocks(blk)%blkm * keep%blocks(blk)%blkn

          end do

          ! Diagonal block has no dependency for factor(dblk)
          keep%blocks(dblk)%dep_initial = cb - 1 

          ! decrement number of row blocks and rows in next block column
          sz = sz - 1
          numrow = numrow - l_nb
       end do
    end do

    allocate (map(n),stat=st)
    if(st.ne.0) go to 9999

    ! Note: following two functions could probably be combined with mc34 call
    ! above, but have been left as is to ease maintenance
    allocate(keep%lmap(keep%nbcol),stat=st)
    if(st.ne.0) go to 9999
    call make_map(n, order, ptr, row, aptr, arow, amap)
    call lcol_map(aptr, arow, num_nodes, keep%nodes, keep%blocks, &
         keep%lmap, map, amap, st)
    if(st.ne.0) goto 9999

    ! before returning take copy of components of info set by MA87_analyse
    keep%info%flag         = info%flag
    keep%info%num_factor   = info%num_factor
    keep%info%num_flops    = info%num_flops
    keep%info%num_nodes    = info%num_nodes
    keep%info%maxdepth     = info%maxdepth
    keep%info%stat         = info%stat

    deallocate (arow,stat=st)
    deallocate (aptr,stat=st)
    deallocate (amap,stat=st)
    deallocate (iw,stat=st)
    deallocate (map,stat=st)
    deallocate (perm,stat=st)

    return

9999 continue

    ! TODO print error
    
    return
  end subroutine spllt_analyse

  ! tree pruning method. Inspired by the strategy employed in qr_mumps
  ! for pruning the atree
  subroutine spllt_prune_tree(num_nodes, sparent, nth)
    implicit none

    integer :: num_nodes ! number of nodes if the atree
    integer, dimension(:), intent(inout) :: sparent ! atree
    integer :: nth ! number of workers

    integer :: node, nlz
    integer, allocatable :: lzero(:)
    real(kind(1.d0)), allocatable   :: lzero_w(:), proc_w(:)

    ! initialize the l0 layer with the root nodes
    nlz = 0
    do node = 1,num_nodes
       if (sparent(node) .eq. 0) then
         nlz = nlz+1
         lzero(nlz) = node
       end if
    end do

    godown: do
       ! if (nth .eq. 1) exit ! serial execution process the whole tree as a subtree
       if(nlz .gt. nth*max(2.d0,(log(real(nth,kind(1.d0)))/log(2.d0))**2)) exit ! exit if already too many nodes in l0

       proc_w = 0.d0

    end do godown

    return
  end subroutine spllt_prune_tree

  ! subroutine spllt_node_flops()
  !   implicit none    

  !   return
  ! end subroutine spllt_node_flops

  subroutine spllt_symbolic(adata, nnodes, sptr, sparent, rptr)
    use spllt_mod
    implicit none

    type(spllt_adata_type), intent(inout) :: adata
    integer, intent(in) :: nnodes ! number of nodes if the atree
    integer, dimension(:), intent(in) :: sptr, sparent
    integer(long), dimension(:), intent(in) :: rptr

    integer :: node, parent
    integer(long) :: mm, m, n, j
    integer(long) :: nflops
    
    allocate(adata%weight(nnodes))

    do node = 1,nnodes

       parent = sparent(node)
       m = rptr(node+1) - rptr(node)
       n = sptr(node+1) - sptr(node)
       
       mm = m-n

       nflops = 0
       do j=1,n
          nflops = nflops + (mm+j)**2
       end do
       
       adata%weight(node) = adata%weight(node) + real(nflops, kind(1.d0)) 
       if (parent .gt. 0) then 
          adata%weight(parent) = adata%weight(parent) + adata%weight(node)
       end if

    end do

    return
  end subroutine spllt_symbolic

end module spllt_analyse_mod
