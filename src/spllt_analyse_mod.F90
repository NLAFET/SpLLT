!> \file
!> \copyright 2016-17 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
module spllt_analyse_mod
  implicit none

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Performs the analysis phase.
  !>
  !> Performs the symbolic factorization of the input matrix. It uses
  !> the ssids_analyse from SSIDSto compute the assembly tree.
  !> @param adata Symbolic factorization data.
  !> @param fdata Factorization data.
  !> @param n Order of the input matrix.
  !> @param ptr Column pointers of A.
  !> @param row Row indices of A.
  !> @param order Return ordering to user / allow user to supply order.
  !> @param options User-supplied options. 
  !> @param info Stats/information returned to user.
  subroutine spllt_analyse(adata, fdata, options, n, ptr, row, info, order)
    use spllt_data_mod
    use spral_ssids_datatypes, only: ssids_options
    use spral_ssids_akeep, only: ssids_akeep
    use spral_ssids_anal, only: expand_pattern
    use spral_ssids, only: ssids_analyse
    implicit none

    type(spllt_adata), intent(inout) :: adata ! data related to the analysis
    type(spllt_fdata), intent(inout) :: fdata ! data related to the factorization
    type(spllt_options), intent(in) :: options ! User-supplied options
    integer, intent(in) :: n ! order of A
    integer, intent(in) :: row(:) ! row indices of lower triangular part
    integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
    integer, dimension(:), optional, intent(inout) :: order
    ! order(i) must hold position of i in the pivot sequence. 
    ! On exit, holds the pivot order to be used by MA87_factor.
    ! type(spllt_keep), target, intent(out) :: keep
    type(spllt_inform), intent(inout) :: info

    integer, allocatable :: amap(:) ! map from user a to reordered a
    integer(long), allocatable :: aptr(:) ! column pointers of expanded matrix
    ! integer, allocatable :: aptr(:) ! column pointers of expanded matrix
    integer, allocatable :: arow(:) ! row pointers of expanded matrix
    integer, allocatable :: iw(:) ! work array
    integer, allocatable :: perm(:) ! inverse permutation.
    ! perm(i) holds the variable that is i-th in the pivot sequence.
    ! Also used for checking user-supplied permutation.
    integer, allocatable :: map(:) ! Allocated to have size n.
    ! used in computing dependency counts. For each row k in 
    ! j-th block column of a node, map(k1) is set to j
    integer(long), allocatable :: ptr64(:)
    
    integer :: num_nodes ! number of nodes in tree
    integer :: node ! a node in tree
    integer :: l_nb ! set to block size of snode (keep%nodes(snode)%nb)
    integer :: ne ! set to ptr(n+1) - 1
    integer(long) :: nz
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
    integer :: ci ! do loop variable. current block col
    integer :: root_node

    integer :: par
    integer, dimension(:), pointer :: sptr, sparent, rlist
    integer(long), dimension(:), pointer :: rptr
    integer, dimension(:), allocatable :: nchild ! pointer on the nth child in the current node

    ! SSIDS data
    type(ssids_options) :: ssids_opt ! SSIDS options
    ! type(ssids_inform) :: inform ! SSDIS stats
    type(ssids_akeep), target :: akeep   ! analysis data from SSIDS

    ! Setup options for analysis in SSIDS
    ssids_opt%ordering = 1 ! use Metis ordering
    ssids_opt%scaling = 0 ! no scaling

    ! Check nemin and set to default if out of range.
    nemin = options%nemin
    if (nemin .lt. 1) nemin = nemin_default

    ssids_opt%nemin = nemin

    ! Perform analysis with SSIDS.
    ! if (present(order)) then
    !    call ssids_analyse(.false., n, ptr, row, akeep, ssids_opt, info%ssids_inform, order)
    ! else
    !    call ssids_analyse(.false., n, ptr, row, akeep, ssids_opt, info%ssids_inform)
    ! end if
    call ssids_analyse(.false., n, ptr, row, akeep, ssids_opt, info%ssids_inform, order=order)
    if (info%ssids_inform%flag .lt. 0) then
       ! Problem occured in SSIDS routine.
       info%flag = SPLLT_ERROR_UNKNOWN
       goto 100
    end if

    ! initialise
    info%flag = 0
    info%num_factor = 0_long
    info%num_flops = 0_long
    info%num_nodes = 0
    info%maxdepth = 0
    info%stat = 0

    fdata%n = n
    adata%n = n
    ne = ptr(n+1) - 1
    nz = ptr(n+1) - 1
    num_nodes = akeep%nnodes

    ! immediate return if n = 0
    if (n == 0) return

    sptr => akeep%sptr
    sparent => akeep%sparent
    rlist => akeep%rlist 
    rptr => akeep%rptr
    
    ! expand the matrix
    
    ! allocate space for expanded matrix (held in aptr,arow)
    allocate(arow(2*nz), aptr(n+3), amap(ptr(n+1)-1), stat=st)
    if (st /= 0) goto 100

    allocate(ptr64(size(ptr,1)))
    ptr64 = ptr

    call expand_pattern(n, nz, ptr64, row, aptr, arow)

    ! allocate (arow(2*ne-n),aptr(n+3),iw(n+1),amap(ptr(n+1)-1),stat=st)
    ! arow(1:ne) = row(1:ne)
    ! aptr(1:n+1) = ptr(1:n+1)
    ! call mc34_expand(n, arow, aptr, iw)
    ! deallocate(iw,stat=st)

    ! nemin = options%nemin
    ! ! Check nemin (a node is merged with its parent if both involve
    ! ! fewer than nemin eliminations). If out of range, use the default
    ! if (nemin < 1) nemin = spllt_nemin_default

    ! ! Check the user-supplied array order and set the inverse in perm.
    ! if (size(order).lt.n) then
    !    go to 9999
    ! end if

    ! deallocate (perm,stat=st)
    ! allocate (perm(n),stat=st)
    ! if (st /= 0) go to 9999
    ! perm(:) = 0
    ! k1 = 0
    ! do i = 1, n
    !    jj = order(i)
    !    if (jj < 1 .or. jj > n) exit
    !    if (perm(jj) /= 0) exit ! Duplicate found
    !    perm(jj) = i
    ! end do
    ! if (i-1 /= n) then
    !    go to 9999
    ! end if
    
    ! control78%nemin = nemin
    ! control78%sort = .true.
    ! control78%lopt = .true.
    ! call mc78_analyse(n, aptr, arow, order, num_nodes, &
    !      sptr, sparent, rptr, rlist, control78, info78, nfact=adata%num_factor, &
    !      nflops=adata%num_flops)
    ! ! print *, "sz sparent:", size(sparent)
    
    adata%nnodes = num_nodes

    ! perform symbolic factorization
    call spllt_symbolic(adata, num_nodes, sptr, sparent, rptr)

    info%num_nodes = num_nodes
    !**************************************
    ! Set up nodal data structures
    ! For each node, hold data in keep%nodes(node) 
    deallocate(fdata%nodes, stat=st)    
    allocate(fdata%nodes(-1:num_nodes+1),stat=st)
    if (st /= 0) goto 100

    deallocate(fdata%nodes, stat=st)
    allocate(fdata%nodes(-1:num_nodes+1),stat=st)
    if (st /= 0) goto 100

    fdata%nodes(0)%blk_en = 0
    fdata%nodes(1)%blk_sa = 1
    fdata%nodes(1)%sa = 1
    
    ! loop over root nodes
    fdata%nodes(:)%nchild = 0
    ! setup nchild and child info
    ! print *, "sz sparent:", size(akeep%sparent)
    ! print *, "num_nodes:", num_nodes
    do node = 1, num_nodes+1

       ! print *, "node: ", node, ", nchild: ", fdata%nodes(node)%nchild
       if (node.le.num_nodes) then
          par = sparent(node)
          fdata%nodes(node)%parent = par 
          fdata%nodes(par)%nchild = fdata%nodes(par)%nchild + 1
       else
          fdata%nodes(node)%parent = -1 ! virutal root node: no parent node
       end if

       ! Allocate space to store child nodes
       ! print *, "node:", node, ", nchild:", fdata%nodes(node)%nchild
       allocate(fdata%nodes(node)%child(fdata%nodes(node)%nchild), stat=st)
       if(st.ne.0) goto 100
    end do

    allocate(nchild(num_nodes+1),stat=st)
    if (st /= 0) goto 100

    ! Add children list to nodes, use nchild as a counter
    nchild(:) = 0
    do node = 1, num_nodes+1
       ! par = sparent(node)
       par = fdata%nodes(node)%parent
       !    ! if(par.gt.num_nodes) cycle
       !    print *, "par: ", par
       if (par .gt. 0) then
          nchild(par) = nchild(par) + 1
          ! print *, "node:", node, "par:", par, "nchild(par):", nchild(par), "sz child:", size(fdata%nodes(par)%child)
          ! nchild(par) = 0
          fdata%nodes(par)%child(nchild(par)) = node
       end if
    end do

    ! goto 9999 ! DEBUG
    ! Setup least descendants, to allow easy walk of subtrees
    do node = -1, num_nodes
       ! initialise least descendat to self
       fdata%nodes(node)%least_desc = node
    end do
    do node = 1, num_nodes
       ! walk up tree from leaves. A parent's least descendant is either this
       ! nodes least descendant (itself in case of a leaf), or its existing
       ! one if that is smaller.
       anode = fdata%nodes(node)%parent
       fdata%nodes(anode)%least_desc = &
            min(fdata%nodes(node)%least_desc, fdata%nodes(anode)%least_desc)
    end do

    ! prune tree
    allocate(adata%small(adata%nnodes))
    if (st /= 0) goto 100
    adata%small = 0
    if (options%prune_tree) then
       call spllt_prune_tree(adata, sparent, options%ncpu, fdata)
       ! call spllt_prune_tree(adata, sparent, 1, fdata) ! TESTS sequential
       ! call spllt_prune_tree(adata, sparent, 2, fdata) ! TESTS
       ! call spllt_prune_tree(adata, sparent, 4, fdata) ! TESTS
       ! call spllt_prune_tree(adata, sparent, 16, fdata) ! TESTS
    end if

    ! set up blocking info
    do node = 1, num_nodes
       fdata%nodes(node)%sa = sptr(node)
       fdata%nodes(node)%en = sptr(node+1)-1
       ! set node id 
       ! fdata%nodes(node)%node => fdata%nodes(node)
       fdata%nodes(node)%num = node
       
       ! determine and record the block size for node
       ! note we are careful in case l_nb**2 overflows (in fact 1+l_nb must
       ! not overflow at the end), and limit the answer to huge(l_nb)/2
       ! if (adata%small(node) .eq. 0) then
       l_nb = options%nb
       if (l_nb < 1) l_nb = nb_default
       ! l_nb = min(huge(l_nb)/2_long, &
       !      (l_nb**2_long) / min(sptr(node+1)-sptr(node), l_nb) )

       ! l_nb = (l_nb-1) / 8 + 1
       ! l_nb = 8 * l_nb

       ! else
       ! l_nb = max(rptr(node+1)-rptr(node), sptr(node+1)-sptr(node))
       ! end if
       fdata%nodes(node)%nb = l_nb

       ! Copy row list into fdata%nodes
       allocate(fdata%nodes(node)%index(rptr(node+1)-rptr(node)),stat=st)
       if (st /= 0) goto 100
       j = 1
       do ii = rptr(node), rptr(node+1)-1
          fdata%nodes(node)%index(j) = rlist(ii)
          j = j + 1
       end do

       ! Calculate number j of blocks in node and set
       ! fdata%nodes(node)%blk_en
       sz = int(rptr(node+1)-rptr(node) - 1) / l_nb + 1
       j = 0
       do i = fdata%nodes(node)%sa, fdata%nodes(node)%en, l_nb
          j = j + sz
          sz = sz - 1
       end do
       fdata%nodes(node)%blk_en = fdata%nodes(node-1)%blk_en + j

       ! if node is not the final node, record first block
       ! for the next node (node+1)
       if (node < num_nodes)  &
            fdata%nodes(node+1)%blk_sa = fdata%nodes(node)%blk_en + 1

       ! if(options%prune_tree) then
       !    fdata%nodes(node)%subtree_root = subtree_root(node)
       ! else
       !    fdata%nodes(node)%subtree_root = -1
       ! endif
    end do

    ! set fdata%final_blk to hold total number of blocks.
    fdata%final_blk = fdata%nodes(num_nodes)%blk_en

    !**************************************   
    ! Fill out block information.    
    deallocate(fdata%bc,stat=st)
    allocate(fdata%bc(fdata%final_blk),stat=st)
    if(st.ne.0) goto 100
    ! if(st.ne.0) go to 9999

#if defined(SPLLT_USE_PARSEC)
    ! allocate diag in fdata
    deallocate(fdata%diags,stat=st)
    allocate(fdata%diags(fdata%final_blk),stat=st) ! large over-estimation of array size
#endif

    ! Loop over the nodes. Number the blocks in the first node
    ! contiguously, then those in the second node, and so on.
    ! Each node has a number of block columns; the blocks within
    ! each block column are numbered contiguously.
    ! Also add up the number of block columns and store largest block dimension.
    blk = 1
    fdata%nbcol = 0
    fdata%maxmn = 0
    do node = 1, num_nodes

       sa = fdata%nodes(node)%sa
       en = fdata%nodes(node)%en
       numcol = en - sa + 1
       numrow = size(fdata%nodes(node)%index)

       ! l_nb is the size of the blocks
       l_nb = fdata%nodes(node)%nb

       ! sz is number of blocks in the current block column
       sz = (numrow - 1) / l_nb + 1

       ! if (adata%small(node) .ge. 0) then

       ! cb is the index of the block col. within node
       cb = 0
       col_used = 0

       ! Loop over the block columns in node. 
       do ci = sa, en, l_nb
          k = 1 ! use k to hold position of block within block column
          ! increment count of block columns
          fdata%nbcol = fdata%nbcol + 1

          cb = cb + 1

          ! blkn is the number of columns in the block column.
          ! For all but the last block col. blkn = l_nb.
          blkn = min(l_nb, numcol-col_used)
          col_used = col_used + blkn

          dblk = blk

#if defined(SPLLT_USE_PARSEC)
          fdata%diags(fdata%nbcol) = dblk
#endif
          ! loop over the row blocks (that is, loop over blocks in block col)
          row_used = 0 
          do blk = dblk, dblk+sz-1
             ! store identity of block
             fdata%bc(blk)%id       = blk
             !print *, "node ", node, "=> blk", blk

             ! store number of rows in the block.
             ! For all but the last block, the number of rows is l_nb
             fdata%bc(blk)%blkm     = min(l_nb, numrow-row_used)
             row_used = row_used + fdata%bc(blk)%blkm

             ! store number of columns in the block.
             fdata%bc(blk)%blkn     = blkn

             fdata%maxmn = max(fdata%maxmn, &
                  fdata%bc(blk)%blkm,  fdata%bc(blk)%blkn)

             ! store position of the first entry of the block within the
             ! block column of L
             fdata%bc(blk)%sa       = k

             ! store identity of diagonal block within current block column
             fdata%bc(blk)%dblk     = dblk

             ! store identity of last block within current block column
             fdata%bc(blk)%last_blk = dblk + sz - 1

             ! store node the blk belongs to
             fdata%bc(blk)%node     = node

             ! initialise dependency count
             fdata%bc(blk)%dep_initial = cb

             ! store identity of block column that blk belongs to
             fdata%bc(blk)%bcol     = fdata%nbcol

             ! increment k by number of entries in block
             k = k + fdata%bc(blk)%blkm * fdata%bc(blk)%blkn

          end do

          ! Diagonal block has no dependency for factor(dblk)
          fdata%bc(dblk)%dep_initial = cb - 1 

          ! decrement number of row blocks and rows in next block column
          sz = sz - 1
          numrow = numrow - l_nb
       end do

       ! else
       !    ! In a subtree, put all entries for a node in a single block
       !    ! TODO error managment
       !    if(sz.ne.1) then
       !       print *, "ERROR: Subtree node should only ever be a single block!"
       !       stop -1
       !    endif

       !    if (adata%small(node) .eq. 1) then 
       !       root_node = node ! the root is the node himself when small equal to 1
       !    else
       !       root_node = -adata%small(node)
       !    end if

       !    if(node .eq. fdata%nodes(root_node)%least_desc) then
       !       ! First node of subtree, all to be in the same bcol array space
       !       k = 1 ! Position within lfact(bcol)%lcol array for start of node
       !    endif

       !    ! increment count of block columns
       !    k = 1 ! FIXME delete (+repoint at lcol, not a subsection, FIXME below)
       !    fdata%nbcol = fdata%nbcol + 1

       !    ! blkn is the number of columns in the block column.
       !    ! For all but the last block col. blkn = l_nb.
       !    blkn = numcol

       !    dblk = blk

       !    ! store identity of block
       !    fdata%blocks(blk)%id       = blk

       !    ! store number of rows in the block.
       !    ! For all but the last block, the number of rows is l_nb
       !    fdata%blocks(blk)%blkm     = numrow

       !    ! store number of columns in the block.
       !    fdata%blocks(blk)%blkn     = blkn

       !    fdata%maxmn = max(fdata%maxmn, &
       !         fdata%blocks(blk)%blkm,  fdata%blocks(blk)%blkn)

       !    ! store position of the first entry of the block within the
       !    ! block column of L
       !    fdata%blocks(blk)%sa       = k

       !    ! store identity of diagonal block within current block column
       !    fdata%blocks(blk)%dblk     = dblk

       !    ! store identity of last block within current block column
       !    fdata%blocks(blk)%last_blk = dblk + sz - 1

       !    ! store node the blk belongs to
       !    fdata%blocks(blk)%node     = node

       !    ! store identity of block column that blk belongs to
       !    fdata%blocks(blk)%bcol     = fdata%nbcol

       !    ! increment k by number of entries in block
       !    k = k + fdata%blocks(blk)%blkm * fdata%blocks(blk)%blkn

       !    blk = dblk + 1

       !    ! Diagonal block has no dependency for factor(dblk)
       !    fdata%blocks(dblk)%dep_initial = 0

       ! end if

    end do

    allocate (map(n),stat=st)
    if(st.ne.0) goto 100

    ! Note: following two functions could probably be combined with mc34 call
    ! above, but have been left as is to ease maintenance
    allocate(fdata%lmap(fdata%nbcol),stat=st)
    if(st.ne.0) goto 100
    call spllt_make_map(n, order, ptr64, row, aptr, arow, amap)
    ! call spllt_make_map(n, order, ptr, row, aptr, arow, amap)
    call spllt_lcol_map(aptr, arow, num_nodes, fdata%nodes, fdata%bc, &
         fdata%lmap, map, amap, st)
    if(st.ne.0) goto 100

#if defined(SPLLT_USE_PARSEC)
    call spllt_compute_dep(adata, fdata, fdata, map)
#endif

    ! before returning take copy of components of info set by MA87_analyse
    fdata%info%flag         = info%flag
    fdata%info%num_factor   = info%num_factor
    fdata%info%num_flops    = info%num_flops
    fdata%info%num_nodes    = info%num_nodes
    fdata%info%maxdepth     = info%maxdepth
    fdata%info%stat         = info%stat

    deallocate (ptr64,stat=st)
    deallocate (arow,stat=st)
    deallocate (aptr,stat=st)
    deallocate (amap,stat=st)
    deallocate (iw,stat=st)
    deallocate (map,stat=st)
    deallocate (perm,stat=st)

    return

100 continue

    ! TODO print error
    
    return
  end subroutine spllt_analyse

#if defined(SPLLT_USE_PARSEC)
  subroutine spllt_compute_upd_bet_add(fdata, node, cptr, cptr2, rptr, rptr2, dest_blk)
    use spllt_data_mod
    implicit none

    type(spllt_fdata), target, intent(inout) :: fdata ! data related to the factorization
    type(spllt_node)     :: node
    integer :: cptr, cptr2, rptr, rptr2
    integer(long)     :: dest_blk

    integer :: s_nb ! block size in node
    integer :: numcol ! number of column in node
    integer :: nc ! number of block column in node
    integer(long) :: dblk ! diag block
    integer :: c ! block column index
    integer :: ljk_sa, ljk_en ! start and end index for ljk contribution
    integer :: lik_sa, lik_en ! start and end index for lik contribution
    integer(long) :: id_ik, id_jk
    integer :: i, j ! index
    type(spllt_block), pointer :: bc_jk, bc_ik, bc_ij
    logical :: sync ! true if task is use for synchronization purpose
    integer :: n1 ! column width of blocks Ljk and Lik
    integer :: csrc, rsrc ! start of fist row respectively in Ljk and Lik block 
    integer :: p, p1, p2 ! dependency index

    ! numcol is number of columns in node
    numcol = node%en-node%sa + 1
    ! block size
    s_nb = node%nb
    ! number of block colum
    nc = (numcol-1) / s_nb + 1 

    ljk_sa = (cptr -1)/s_nb
    ljk_en = (cptr2-1)/s_nb
    
    lik_sa = (rptr -1)/s_nb
    lik_en = (rptr2-1)/s_nb    

    ! first (diag) block
    dblk = node%blk_sa
    
    do c = 1, nc

       sync = .true.

       do j = ljk_en,ljk_sa,-1
          do i = lik_en, max(lik_sa, j), -1
             
             ! if (i .lt. j) cycle

             ! compute current ljk block id
             id_jk = dblk + j - (c-1)
             ! compute current lik block id
             id_ik = dblk + i - (c-1)

             ! compute first row of Ljk block
             n1 = fdata%bc(dblk)%blkn
             csrc  = 1 + (mod(cptr-1, s_nb))*n1
             rsrc  = 1 + (mod(rptr-1, s_nb))*n1

             ! write(*,*)'[spllt_analyse_mod][compute_dep] add upd_bet, id_jk: ', &
             !      & id_jk, ', id_ik: ', id_ik, ', id_ij: ', dest_blk
             
             bc_jk => fdata%bc(id_jk)
             bc_ik => fdata%bc(id_ik)
             bc_ij => fdata%bc(dest_blk)

             if ((j .eq. ljk_sa) .and. (i .eq. lik_sa)) sync = .false.
             ! sync = .false.

             p  = spllt_dep_in_add(bc_ij%dep_in, id_jk, id_ik, csrc, rsrc, sync)

             p1 = spllt_dep_out_add(bc_jk%dep_out, dest_blk, 1)
             p2 = spllt_dep_out_add(bc_ik%dep_out, dest_blk, 2)
             
             bc_ij%dep_in(p)%p1 = p1
             bc_ij%dep_in(p)%p2 = p2

             bc_jk%dep_out(p1)%p = p
             bc_ik%dep_out(p2)%p = p
             
          end do
       end do

       dblk = fdata%bc(dblk)%last_blk + 1
    end do

    return
  end subroutine spllt_compute_upd_bet_add

  subroutine spllt_compute_dep(adata, fdata, map)
    use spllt_data_mod
    use spllt_kernels_mod
    implicit none

    type(spllt_adata), intent(inout) :: adata ! data related to the analysis
    type(spllt_fdata), target, intent(inout) :: fdata ! data related to the factorization
    ! type(spllt_keep), target, intent(in) :: keep
    integer, allocatable :: map(:)

    type(spllt_node), pointer     :: node, anode ! node in the atree
    integer :: num_nodes, snum, anum
    integer :: numcol, numrow
    integer :: s_nb
    integer :: cptr, cptr2
    integer :: n1 ! column width
    logical :: map_done
    integer :: jb, cb, jlast
    integer(long) :: blk, dblk, a_blk, a_dblk
    integer :: i, k, ii, ilast
    integer :: blk_jk_sa, blk_jk_en
    integer :: lik, ljk, ljk_sa, ljk_en
    integer :: nr, nc
    integer(long) :: id_ik, id_jk, id, id_jk_sa, id_jk_en
    type(spllt_block), pointer :: bc_jk, bc_ik, bc_ij

    write(*,*)'[spllt_analyse_mod][compute_dep]'

    num_nodes = adata%nnodes

    ! print *, "[compute_dep] num_nodes: ", num_nodes

    do snum = 1, num_nodes

       node => fdata%nodes(snum)

       ! parent node
       anum = node%parent
       numrow = size(node%index)
       ! numcol is number of columns in node
       numcol = node%en-node%sa + 1
       ! block size
       s_nb = node%nb
       ! number of block colum
       nc = (numcol-1) / s_nb + 1 
       ! number of block row
       nr = (numrow-1) / s_nb + 1 

       ! initialise cptr 
       cptr = 1 + numcol

       ! loop over ancestors of node
       do while(anum.gt.0)
          anode => fdata%nodes(anum)

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
                a_dblk = fdata%bc(a_dblk)%last_blk + 1
             end do

             ! Find cptr2
             jlast = min(anode%sa + cb*anode%nb - 1, anode%en)
             do cptr2 = cptr,numrow
                if(node%index(cptr2) > jlast) exit
             end do
             cptr2 = cptr2 - 1
             
             ljk_sa = (cptr -1)/s_nb
             ljk_en = (cptr2-1)/s_nb
             
             if(.not.map_done) call spllt_build_rowmap(anode, map) 

             ! Loop over the blocks of snode
             ii = map(node%index(cptr)) 
             ! ii = -1 
             ilast = cptr ! Set start of current block
             
             do i = cptr, numrow
                k = map(node%index(i))

                if(k.ne.ii) then

                   a_blk = a_dblk + ii - cb
                   ! add dep for updating a_blk block 
                   call spllt_compute_upd_bet_add(fdata, node, cptr, cptr2, ilast, i-1, a_blk)
                   
                   ii = k
                   ilast = i ! Update start of current block
                end if
             end do

             a_blk = a_dblk + ii - cb
             ! add dep for updating a_blk block 
             call spllt_compute_upd_bet_add(fdata, node, cptr, cptr2, ilast, i-1, a_blk)

             ! Move cptr down, ready for next block column of anode
             cptr = cptr2 + 1
          end do bcols
          
          ! Move up the tree to the parent of anode
          anum = fdata%nodes(anum)%parent
       end do

    end do

    return
  end subroutine spllt_compute_dep
#endif

  ! Tree pruning method. Inspired by the strategy employed in qr_mumps
  ! for pruning the atree.
  subroutine spllt_prune_tree(adata, sparent, nth, fdata)
    use spllt_data_mod
    use spllt_utils_mod
    implicit none
    
    type(spllt_adata), intent(inout) :: adata
    integer, dimension(:), intent(in) :: sparent ! atree
    integer :: nth ! number of workers
    type(spllt_fdata), target, intent(in) :: fdata

    integer                         :: i, j
    integer                         :: c
    integer :: nlz ! number of nodes in the lzero layer
    integer :: node , leaves, totleaves
    integer                         :: n ! node to be replaced by its direct descendents in lzero
    integer(long) :: p, totflops 
    real(kind(1.d0))                :: rm, smallth
    integer, allocatable       :: lzero(:)
    integer(long), allocatable :: lzero_w(:), proc_w(:)
    logical                         :: found

    write(*,*)'[analysis][prune_tree] nth: ', nth
    
    allocate(lzero_w (adata%nnodes))
    allocate(lzero   (adata%nnodes))
    allocate(proc_w  (nth))

    smallth = 0.01

10 continue
    totleaves = 0
    adata%small = 0
    
    totflops = adata%weight(adata%nnodes+1)
    ! write(*,*)'totflops: ', totflops
    ! write(*,*)'weights: ', adata%weight 
    ! write(*,*)'nnodes: ', adata%nnodes
    ! write(*,*)'root: ', sparent(1)
    ! initialize the l0 layer with the root nodes
    nlz = 0

    ! do node = 1, adata%nnodes+1
    !    write(*,*) 'node: ', node, ', weight: ', adata%weight(node)
    !    if (sparent(node) .gt. adata%nnodes) then
    !    ! if (sparent(node) .le. 0) then
    !       ! write(*,*) 'weight: ', real(adata%weight(node), kind(1.d0))
    !       ! write(*,*) 'thresh: ', smallth*real(totflops, kind(1.d0))
    !       if(real(adata%weight(node), kind(1.d0)) .gt. smallth*real(totflops, kind(1.d0))) then
    !          nlz = nlz+1
    !          lzero(nlz) = node
    !          lzero_w(nlz) = -adata%weight(node)
    !          write(*,*) 'node: ', node,', weight: ', adata%weight(node)
    !       else
    !          adata%small(node) = 1 ! node is too small; mark it
    !       end if
    !    end if
    !    if(keep%nodes(node)%nchild .eq. 0) totleaves = totleaves+1
    ! end do

    node = adata%nnodes+1 ! root node (symbolic)
    nlz = nlz+1
    lzero(nlz) = node
    lzero_w(nlz) = -adata%weight(node)
    ! count leaf nodes
    do node = 1, adata%nnodes+1
       if(fdata%nodes(node)%nchild .eq. 0) totleaves = totleaves+1       
    end do
    
    leaves = 0

    godown: do
       ! if (nth .eq. 1) exit ! serial execution process the whole tree as a subtree
       if (nlz .le. 0) exit ! only small nodes ! 
       if(nlz .gt. nth*max(2.d0,(log(real(nth,kind(1.d0)))/log(2.d0))**2)) exit ! exit if already too many nodes in l0

       proc_w = 0
       
       ! sort lzero_w into ascending order and apply the same order on
       ! lzero array
       call spllt_sort(lzero_w, nlz, map=lzero)
       ! write(*,*) 'lzero_w: ', lzero_w(1:nlz)
       ! map subtrees to threads round-robin 
       do node=1, nlz
          ! find the least loaded proc
          p = minloc(proc_w,1)
          proc_w(p) = proc_w(p) + abs(lzero_w(node))
       end do
       ! write(*,*)'nlz: ', nlz       
       ! all the subtrees have been mapped. Evaluate load balance
       rm = minval(proc_w)/maxval(proc_w)
       ! print *, "rm: ", rm
       if((rm .gt. 0.9) .and. (nlz .ge. 1*nth)) exit ! if balance is higher than 90%, we're happy

       ! if load is not balanced, replace heaviest node with its kids (if any)
       found = .false.
       findn: do
          if(leaves .eq. totleaves) exit godown ! reached the bottom of the tree
        
          if(leaves .eq. nlz) then
             if(nlz .ge. nth*max(2.d0,(log(real(nth,kind(1.d0)))/log(2.d0))**2)) then 
                exit godown ! all the nodes in l0 are leaves. nothing to do
             else
                smallth = smallth/2.d0
                if(smallth .lt. 1e-4) then
                   exit godown
                else
                   goto 10
                end if
             end if
          end if
          n = lzero(leaves+1) ! n is the node that must be replaced
          ! print *, "n:", n, ", nchild:", fdata%nodes(n)%nchild, ", sz:", size(fdata%nodes(n)%child)
          ! print *, "children:", fdata%nodes(n)%child
          ! append children of n
          do i=1, size(fdata%nodes(n)%child) ! fdata%nodes(n)%nchild
             c = fdata%nodes(n)%child(i)
             ! print *, "c:", c             
             if(real(adata%weight(c), kind(1.d0)) .gt. smallth*real(totflops, kind(1.d0))) then
                ! this child is big enough, add it
                found = .true.
                nlz = nlz+1
                lzero  (nlz) = c
                lzero_w(nlz) = -adata%weight(c)
             else
                ! print *, "TETET"
                adata%small(fdata%nodes(c)%least_desc:c) = -c
                adata%small(c) = 1 ! node is too smal; mark it
             end if

          end do
          if(found) exit findn ! if at least one child was added then we redo the mapping
          leaves = leaves+1          
       end do findn
       ! write(*,*) 'lzero: ', lzero(1:nlz)

       ! swap n with last element
       lzero  (leaves+1) = lzero  (nlz)
       lzero_w(leaves+1) = lzero_w(nlz)
       nlz = nlz-1

    end do godown

    ! write(*,*)'nlz: ', nlz
    ! write(*,*)'final lzero: ', lzero(1:nlz)

    ! DEBUG
    ! adata%small = 0
    ! nlz = 1
    ! lzero(1) = 1

    ! mark all the children of nodes in l0
    do i=1, nlz
       n = lzero(i)
       ! write(*,*)'n: ', n
       ! print *, "n:", n, ", nchild:", fdata%nodes(n)%nchild, ", sz:", size(fdata%nodes(n)%child) 
       do j=1, size(fdata%nodes(n)%child) ! fdata%nodes(n)%nchild
          c = fdata%nodes(n)%child(j)
          ! print *, "c: ", c
          ! write(*,*)'desc: ', fdata%nodes(c)%least_desc
          adata%small(fdata%nodes(c)%least_desc:c) = -c
          adata%small(c) = 1
       end do
    end do

    ! print *, "TETETET"

    ! DEBUG
    ! adata%small = 0
    ! c = 703
    ! adata%small(fdata%nodes(c)%least_desc:c) = -c
    ! adata%small(c) = 1

    ! c = 668
    ! adata%small(fdata%nodes(c)%least_desc:c) = -c
    ! adata%small(c) = 1

    deallocate(lzero_w)
    deallocate(lzero)
    deallocate(proc_w)

    return
  end subroutine spllt_prune_tree

  ! Performs the symbolic factorization of the assembly tree.
  subroutine spllt_symbolic(adata, nnodes, sptr, sparent, rptr)
    use spllt_data_mod
    implicit none

    type(spllt_adata), intent(inout) :: adata
    integer, intent(in) :: nnodes ! number of nodes if the atree
    integer, dimension(:), intent(in) :: sptr, sparent
    integer(long), dimension(:), intent(in) :: rptr

    integer :: node, parent
    integer(long) :: mm, m, n, j
    integer(long) :: nflops
    
    allocate(adata%weight(nnodes+1))

    adata%weight = 0 

    do node = 1,nnodes

       parent = sparent(node)
       m = rptr(node+1) - rptr(node)
       n = sptr(node+1) - sptr(node)
       
       mm = m-n

       nflops = 0
       do j=1,n
          nflops = nflops + (mm+j)**2
       end do
       ! write(*,*) 'nflops: ', nflops
       adata%weight(node) = adata%weight(node) + nflops
       adata%weight(parent) = adata%weight(parent) + adata%weight(node)

    end do
    ! write(*,*) 'symbolic totflops: ', adata%weight(nnodes+1)

    ! print *, "root weight: ", adata%weight(nnodes+1)

    return
  end subroutine spllt_symbolic

  ! Make a map from original A to reordered half matrix A
  ! The reordered half matrix's pattern is returned in nptr and nrow
  subroutine spllt_make_map(n, perm, optr, orow, nptr, nrow, map)
    use spllt_data_mod
    implicit none

    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: perm
    integer(long), dimension(n+1), intent(in) :: optr
    ! integer, dimension(n+1), intent(in) :: optr
    integer, dimension(optr(n+1)-1), intent(in) :: orow
    integer(long), dimension(n+3), intent(out) :: nptr ! extra space used for tricks
    ! integer, dimension(n+3), intent(out) :: nptr ! extra space used for tricks
    integer, dimension(optr(n+1)-1), intent(out) :: nrow
    integer, dimension(optr(n+1)-1), intent(out) :: map

    integer :: i, k, l
    integer(long) :: j

    nptr(:) = 0

    ! Count number of entries in each column of new matrix (at offset 2)
    do i = 1, n
       l = perm(i)
       do j = optr(i), optr(i+1)-1
          k = perm(orow(j))
          if(k<l) then
             nptr(k+2) = nptr(k+2) + 1
          else
             nptr(l+2) = nptr(l+2) + 1
          endif
       end do
    end do

    ! Calculate column starts (at offset 1)
    nptr(1:2) = 1
    do i = 2, n
       nptr(i+1) = nptr(i) + nptr(i+1)
    end do

    ! Now build map
    do i = 1, n
       l = perm(i)
       do j = optr(i), optr(i+1)-1
          k = perm(orow(j))
          if(k<l) then
             map(nptr(k+1)) = j
             nrow(nptr(k+1)) = l
             nptr(k+1) = nptr(k+1) + 1
          else
             map(nptr(l+1)) = j
             nrow(nptr(l+1)) = k
             nptr(l+1) = nptr(l+1) + 1
          endif
       end do
    end do
  end subroutine spllt_make_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Build mapping on per block column basis from user's val to block col's lcol
  ! This routine uses the reordered half matrix and map from the make_map routine
  subroutine spllt_lcol_map(aptr, arow, num_nodes, nodes, blocks, lmap, map, amap, st)
    use spllt_data_mod
    implicit none

    ! Reordered lower triangle of reordered matrix held using aptr and arow
    integer(long), dimension(:), intent(in) :: aptr
    ! integer, dimension(:), intent(in) :: aptr
    integer, dimension(:), intent(in) :: arow
    integer, intent(in) :: num_nodes
    type(spllt_node), dimension(-1:), intent(in) :: nodes ! Node info
    type(spllt_block), dimension(:), intent(in) :: blocks ! block info
    type(lmap_type), dimension(:), intent(out) :: lmap ! output lcol map
    integer, dimension(:), intent(out) :: map ! work array
    integer, dimension(:), intent(in) :: amap ! map set up by make_map
    integer, intent(out) :: st

    ! Local scalars
    integer :: bcol ! block column
    integer :: cb ! Temporary variable
    integer :: col ! do loop index
    integer(long) :: dblk ! set to fdata%nodes(snode)%blk_sa (first block
    ! in snode which is, of course, a diagonal block)
    integer :: en ! set fdata%nodes(snode)%en
    integer(long) :: i ! Temporary variable. global row index
    integer(long) :: j ! Temporary variable
    integer :: l_nb ! set to fdata%nodes(snode)%nb
    integer(long) :: offset
    integer :: sa ! set fdata%nodes(snode)%sa
    integer :: snode ! node
    integer :: swidth ! set to fdata%blocks(dblk)%blkn (number of columns
    ! in block column to which dblk belongs)

    integer :: k

    st = 0

    do snode = 1, num_nodes ! loop over nodes
       ! Build a map from global to local indices
       do j = 1, size(nodes(snode)%index)
          i = nodes(snode)%index(j)
          map(i) = j - 1
       end do

       ! Fill in lfact by block columns
       dblk = nodes(snode)%blk_sa

       l_nb = nodes(snode)%nb
       sa = nodes(snode)%sa
       en = nodes(snode)%en

       do cb = sa, en, l_nb
          bcol = blocks(dblk)%bcol

          offset = blocks(dblk)%sa - (cb-sa)*blocks(dblk)%blkn
          swidth = blocks(dblk)%blkn

          k = 1
          lmap(bcol)%len_map = aptr(min(cb+l_nb-1,en)+1) - aptr(cb)
          allocate(lmap(bcol)%map(2,lmap(bcol)%len_map), stat=st)
          if (st.ne.0) return

          ! Loop over columns in the block column bcol
          do col = cb, min(cb+l_nb-1, en)
             ! loop over rows in column col
             do j = aptr(col), aptr(col+1)-1
                i = arow(j)
                i = map(i)
                i = offset + i*swidth
                lmap(bcol)%map(1,k) = i ! destination in lfact(:)%lcol
                lmap(bcol)%map(2,k) = amap(j) ! source
                k = k + 1
             end do
             offset = offset + 1
          end do
          ! move to next block column in snode
          dblk = blocks(dblk)%last_blk + 1
       end do
    end do
  end subroutine spllt_lcol_map

end module spllt_analyse_mod
