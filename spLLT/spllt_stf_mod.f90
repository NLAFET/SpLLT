module spllt_stf_mod
  use spllt_mod
  implicit none

contains

  subroutine spllt_stf_factorize(n, ptr, row, val, order, keep, control, info)
    use hsl_ma87_double
    implicit none
    
    integer, intent(in) :: n ! order of A
    integer, intent(in) :: row(:) ! row indices of lower triangular part
    integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
    real(wp), intent(in) :: val(:) ! matrix values
    integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
    ! since the analyse phase)

    type(MA87_keep), intent(inout) :: keep 
    type(MA87_control), intent(in) :: control 
    type(MA87_info), intent(out) :: info 

    ! local arrays
    real(wp), dimension(:), allocatable :: detlog ! per thread sum of log pivot
    ! integer, dimension(:), allocatable ::  invp ! used to hold inverse ordering
    integer, dimension(:), allocatable ::  map ! allocated to have size n.
    ! used in copying entries of user's matrix a into factor storage 
    ! (keep%fact).
 
    ! local scalars
    integer(long) :: blk, blk1, blk2 ! block identity
    integer :: blkm, blkn, blkm1, blkn1, blkm2, blkn2 ! number of rows/cols in block
    integer :: dblkm,dblkn ! number of rows/cols in block
    integer(long) :: dblk ! diagonal block within block column
    integer :: en ! holds keep%nodes(snode)%en
    integer :: sa, diagsa, sa1, sa2 ! holds keep%nodes(snode)%sa
    integer(long) :: i, id, diagid
    integer :: j
    integer :: l_nb ! set to block size of snode (keep%nodes(snode)%nb)
    integer :: nbcol, nbrow ! number of block column/row
    integer :: size_bcol ! number of entries in the block column (sum over the
    ! row blocks in the column)
    integer :: sz ! number of blocks in a block column of snode
    integer :: snode, node, num_nodes, par 
    integer :: flag ! Error flag
    integer :: st ! stat parameter
    integer :: numrow, numcol
    integer :: bcol, dbcol, bcol_src
    integer :: total_threads ! number of threads being used
    integer :: this_thread
    integer :: mm, nn, kk
    
    ! update between variables
    integer :: anode ! Ancestor of node
    integer :: cptr  ! Position in snode of the first row 
    ! matching a column of the current block column of anode.
    logical :: map_done ! True if map has been built for anode.
    integer :: a_nb  ! Block size of anode
    integer(long) :: rblk ! id of block in scol containing row 
    ! nodes(snode)%index(cptr).
    
    ! real(wp) :: soln(0)

    ! call factorize_posdef(n, val, order, keep, control, info, 0, 0, soln)

    write(*,*) 'control%nb: ', control%nb

    ! Initialise
    flag = 0

    num_nodes = keep%info%num_nodes
    write(*,*) 'num_nodes: ', num_nodes

    ! Set up inverse permutation
    ! deallocate (invp,stat=st)
    ! allocate (invp(n),stat=st)
    ! if(st.ne.0) go to 10

    ! do j = 1, n
    !    invp(order(j)) = j
    ! end do
    
    ! Allocate factor storage (in keep%lfact)
    deallocate (keep%lfact,stat=st)
    allocate (keep%lfact(keep%nbcol),stat=st)
    if(st.ne.0) go to 10
    
    blk = 1
    nbcol = 0
    ! loop over the nodes
    do snode = 1, num_nodes
       ! Loop over the block columns in snode, allocating space 
       ! l_nb is the size of the blocks and sz is number of
       ! blocks in the current block column
       l_nb = keep%nodes(snode)%nb
       sz = (size(keep%nodes(snode)%index) - 1) / l_nb + 1
       sa = keep%nodes(snode)%sa
       en = keep%nodes(snode)%en

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
          if(st.ne.0) go to 10
       end do
    end do

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

    call copy_a_to_l(n,num_nodes,val,map,keep)

    ! factorize nodes
    ! blk = 1
    do node = 1, num_nodes
       
       write(*,*) 'snode: ', node, & 
            & ', parent: ', keep%nodes(node)%parent, &
            & ', nchild: ', keep%nodes(node)%nchild

       sa = keep%nodes(node)%sa
       en = keep%nodes(node)%en
       numcol = en - sa + 1
       numrow = size(keep%nodes(node)%index)

       nbcol = (numcol-1) / l_nb + 1 
       nbrow = (numrow-1) / l_nb + 1 

       ! l_nb is the size of the blocks
       l_nb = keep%nodes(node)%nb

       ! sz is number of blocks in the current block column
       sz = (numrow - 1) / l_nb + 1

       write(*,*) 'numrow: ', numrow, ', numcol: ', numcol, ', sz: ', sz
       ! write(*,*) 'l_nb: ', l_nb

       ! first block in node
       blk = keep%nodes(node)%blk_sa

       ! Loop over the block columns in node. 
       do kk = 1, nbcol

          ! A_kk          
          dblk = blk

          ! use n, m, and sa to point to where blk is stored in lfact
          dblkn  = keep%blocks(dblk)%blkn
          dblkm  = keep%blocks(dblk)%blkm
          diagsa = keep%blocks(dblk)%sa
          diagid = keep%blocks(dblk)%id

          ! bcol is block column that blk belongs to
          dbcol = keep%blocks(blk)%bcol

          ! factorize_block task
          call factor_diag_block(dblkn, dblkm, diagid, &
               & keep%lfact(dbcol)%lcol(diagsa:diagsa+dblkn*dblkm-1),   &
               & control, flag, detlog(this_thread))

          ! loop over the row blocks (that is, loop over blocks in block col)
          do mm = kk+1,nbrow
          ! do blk = dblk+1, dblk+sz-1
             
             ! A_mk
             blk = dblk+mm-kk

             ! write(*,*)'blk: ', blk

             blkn  = keep%blocks(blk)%blkn
             blkm  = keep%blocks(blk)%blkm
             sa = keep%blocks(blk)%sa
             id = keep%blocks(blk)%id

             ! bcol is block column that blk and dblk belong to
             bcol = keep%blocks(blk)%bcol

             ! solve_block task
             call solv_col_block(blkm, blkn, id, & 
                  & keep%lfact(bcol)%lcol(sa:sa+blkn*blkm-1), &
                  & diagid, keep%lfact(bcol)%lcol(diagsa:diagsa+blkn*dblkm), control)


          end do

          do nn = kk+1,nbcol
             do mm = nn,nbrow
                
                ! L_mk
                blk1 = dblk+mm-kk                
                blkn1 = keep%blocks(blk1)%blkn
                blkm1 = keep%blocks(blk1)%blkm
                sa1   = keep%blocks(blk1)%sa

                ! L_nk
                blk2 = dblk+nn-kk
                blkn2 = keep%blocks(blk2)%blkn
                blkm2 = keep%blocks(blk2)%blkm
                sa2   = keep%blocks(blk2)%sa
                
                ! A_mn
                blk = get_dest_block(keep%blocks(blk1), keep%blocks(blk2))
                blkn = keep%blocks(blk)%blkn
                blkm = keep%blocks(blk)%blkm
                sa   = keep%blocks(blk)%sa                
                
                bcol = keep%blocks(blk)%bcol
                bcol_src = keep%blocks(blk1)%bcol

                call update_block_block(blkm, blkn, &
                     & keep%lfact(bcol)%lcol(sa:sa+blkn*blkm-1),  &
                     & keep%blocks(blk), blkn1, &
                     & keep%lfact(bcol_src)%lcol(sa1:sa1+blkn1*blkm1-1), &
                     & keep%lfact(bcol_src)%lcol(sa2:sa2+blkn1*blkm2-1), control)

             end do
          end do

          ! update between
          

          ! decrement number of row blocks and rows in next block column
          sz = sz - 1
          numrow = numrow - l_nb
       end do

    end do

10 if(st.ne.0) then
      info%flag = spllt_error_allocation
      info%stat = st
      call spllt_print_err(info%flag, control, context='spllt_stf_factorize',st=st)
      return
   endif

    return
  end subroutine spllt_stf_factorize

end module spllt_stf_mod
