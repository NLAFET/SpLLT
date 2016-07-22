! node factorization
! Submit the DAG for the factorization of a node
subroutine spllt_factorize_node(snode, map, fdata, keep, control)
  use spllt_data_mod
  use hsl_ma87_double
  use spllt_kernels_mod
  use spllt_factorization_task_mod
  implicit none

  type(spllt_node_type), target, intent(inout)        :: snode ! node to factorize (spllt)    
  integer, dimension(:), pointer, intent(inout)       :: map
  type(spllt_data_type), target, intent(inout)        :: fdata
  type(MA87_keep), target, intent(inout)              :: keep 
  type(MA87_control), intent(in)                      :: control 

  type(node_type), pointer :: node ! node to factorize (hsl_ma87)
  integer :: num_nodes ! number of nodes in etree
  integer :: snum ! node idx
  integer :: prio ! task priority
  ! shortcuts
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

  ! update between
  type(node_type), pointer :: anode ! ancestor node in the atree
  type(spllt_node_type), pointer :: a_node ! ancestor node in the atree
  ! locate source blocks
  integer :: a_num ! ancestor id
  integer :: cptr  ! Position in snode of the first row 
  ! matching a column of the current block column of anode.
  integer :: cptr2  ! Position in snode of the last row 
  ! matching a column of the current block column of anode.
  logical :: map_done
  integer :: i, ilast
  type(spllt_bc_type), pointer :: bc, a_bc ! node in the atree
  integer :: cb, jb
  integer :: jlast ! Last column in the cb-th block column of anode
  integer :: k
  integer(long) :: a_dblk, a_blk ! id of block in scol containing row 

  ! write(*,*) "---------- node ----------"
  ! write(*,*) 'snode: ', snode 

  num_nodes = keep%info%num_nodes
  snum = snode%num 

  node => snode%node

  ! node priority
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

  ! write(*,'("numrow,numcol = ",i5,i5)')numrow, numcol
  ! write(*,'("nr,nc = ",i5,i5)')nr, nc

  ! first block in node
  dblk = node%blk_sa

  do kk = 1, nc

     ! A_kk
     bc_kk => fdata%bc(dblk)
     call spllt_factorize_block_task(fdata, fdata%nodes(snum), bc_kk, keep%lfact, prio+3)

     ! loop over the row blocks (that is, loop over blocks in block col)
     do ii = kk+1,nr
        ! do ii = nr,kk+1,-1             
        ! A_mk
        blk = dblk+ii-kk
        ! bc_ik => keep%blocks(blk)
        bc_ik => fdata%bc(blk)

        call spllt_solve_block_task(fdata, bc_kk, bc_ik, keep%lfact,prio+2)
     end do


     do jj = kk+1,nc

        ! L_jk
        blk2 = dblk+jj-kk
        bc_jk => fdata%bc(blk2)

        do ii = jj,nr

           ! L_ik
           blk1 = dblk+ii-kk                
           bc_ik => fdata%bc(blk1)

           ! A_ij
           ! blk = get_dest_block(keep%blocks(blk1), keep%blocks(blk2))
           blk = get_dest_block(keep%blocks(blk2), keep%blocks(blk1))
           bc_ij => fdata%bc(blk)
           call spllt_update_block_task(fdata, bc_ik, bc_jk, bc_ij, keep%lfact, prio+1)

        end do
     end do

     ! move to next block column in snode
     dblk = keep%blocks(dblk)%last_blk + 1
     ! numrow = numrow - s_nb
  end do

  ! update between

  !  Loop over ancestors of snode
  a_num = node%parent  
  !  Initialize cptr to correspond to the first row of the rectangular part of 
  !  the snode matrix.
  cptr = 1 + numcol
  ! write(*,*)"numrow: ", numrow
  do while(a_num.gt.0)
     anode => keep%nodes(a_num) 
     a_node => fdata%nodes(a_num)
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
           a_dblk = keep%blocks(a_dblk)%last_blk + 1
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
              a_bc => fdata%bc(a_blk)

              ! rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
              ! rsrc(2) = (i - ilast)*blocks(blk)%blkn

              dblk = node%blk_sa
              ! Loop over the block columns in node. 
              do kk = 1, nc

                 bc_kk => fdata%bc(dblk)

                 call spllt_update_between_task( &
                      & fdata, &
                                ! & bc, &
                      & bc_kk, &
                      & node, a_bc, a_node, &
                                ! & csrc, rsrc, &
                      & cptr, cptr2, ilast, i-1, &
                      & fdata%row_list, fdata%col_list, fdata%workspace, &
                      & keep%lfact, keep%blocks, fdata%bc, &
                      & control, prio)

                 dblk = keep%blocks(dblk)%last_blk + 1
              end do

              ii = k
              ilast = i ! Update start of current block
           end if
        end do
        ! #if defined(SPLLT_USE_OMP)
        ! !$omp taskwait
        ! #endif
        ! i = min(i,numrow)
        a_blk = a_dblk + ii - cb
        ! a_bc => keep%blocks(a_blk)
        a_bc => fdata%bc(a_blk)

        ! rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
        ! rsrc(2) = (i - ilast)*blocks(blk)%blkn

        dblk = node%blk_sa
        ! Loop over the block columns in node. 
        do kk = 1, nc

           bc_kk => fdata%bc(dblk)

           call spllt_update_between_task( &
                & fdata, &
                                ! & bc, &
                & bc_kk, &
                & node, a_bc, a_node, &
                                ! & csrc, rsrc, &
                & cptr, cptr2, ilast, i-1, &
                & fdata%row_list, fdata%col_list, fdata%workspace, &
                & keep%lfact, keep%blocks, fdata%bc, &
                & control, prio)

           dblk = keep%blocks(dblk)%last_blk + 1
        end do

        ! Move cptr down, ready for next block column of anode
        cptr = cptr2 + 1
     end do bcols

     a_num = anode%parent
  end do

  return
end subroutine spllt_factorize_node
