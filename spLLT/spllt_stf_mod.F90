module spllt_stf_mod
  use spllt_mod
  implicit none

contains

  subroutine spllt_stf_factorize(n, ptr, row, val, order, keep, control, info, fdata, cntl)
    use spllt_mod
    use hsl_ma87_double
    use spllt_factorization_mod
#if defined(SPLLT_USE_STARPU) 
    use iso_c_binding
    use starpu_f_mod
#endif
    use spllt_kernels_mod, only : spllt_activate_node
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

    type(spllt_data_type), target :: fdata
    type(spllt_cntl)      :: cntl

    ! local arrays
    ! integer, dimension(:), allocatable ::  invp ! used to hold inverse ordering
    integer, dimension(:), allocatable ::  colmap, map ! allocated to have size n.
    ! used in copying entries of user's matrix a into factor storage 
    ! (keep%fact).

    ! shortcuts
    type(node_type), pointer     :: node ! node in the atree    
    type(spllt_bc_type), pointer :: bc_kk, bc_ik, bc_jk, bc_ij
    type(block_type), dimension(:), pointer :: blocks ! block info. 
    type(node_type), dimension(:), pointer :: nodes 
 
    ! local scalars
    integer(long) :: blk, blk1, blk2 ! block identity
    integer(long) :: dblk ! diagonal block within block column
    integer :: en ! holds keep%nodes(snode)%en
    integer :: sa ! holds keep%nodes(snode)%sa
    integer :: s_nb ! set to block size of snode (keep%nodes(snode)%nb)
    integer :: nc, nr ! number of block column/row
    integer :: snode, num_nodes
    integer :: st ! stat parameter
    integer :: numrow, numcol
    integer :: ii, jj, kk
    integer :: prio
    
    ! update between variables
    ! integer :: csrc(2), rsrc(2) ! used for update_between tasks to
    type(node_type), pointer :: anode ! ancestor node in the atree
    ! locate source blocks
    integer :: a_num ! ancestor id
    integer :: cptr  ! Position in snode of the first row 
    ! matching a column of the current block column of anode.
    integer :: cptr2  ! Position in snode of the last row 
    ! matching a column of the current block column of anode.
    logical :: map_done
    integer :: i, ilast
    type(spllt_bc_type), pointer :: bc, a_bc ! node in the atree
    integer, dimension(:), allocatable :: row_list, col_list ! update_buffer workspace
    integer :: cb, jb
    integer :: jlast ! Last column in the cb-th block column of anode
    ! real(wp) :: soln(0)
    integer :: k
    integer(long) :: a_dblk, a_blk ! id of block in scol containing row 
    ! nodes(snode)%index(cptr).
    ! integer(long) :: rb ! Index of block row in snode

#if defined(SPLLT_USE_STARPU) 
    ! when using StarPU
    integer(c_int) :: ret
#endif
    ! timing
    integer :: start_t, stop_t, rate_t
    integer :: stf_start_t, stf_stop_t, stf_rate_t
    integer :: start_nosub_t, rate_nosub_t
    integer :: start_setup_t, stop_setup_t, rate_setup_t
    integer :: start_cpya2l_t, stop_cpya2l_t, rate_cpya2l_t
    ! call system_clock(start_t, rate_t)
    call system_clock(start_setup_t, rate_setup_t)

    ! call factorize_posdef(n, val, order, keep, control, info, 0, 0, soln)

    ! write(*,*) 'control%nb: ', control%nb

    ! shortcut
    blocks => keep%blocks
    nodes  => keep%nodes

    num_nodes = keep%info%num_nodes
    ! write(*,*) 'num_nodes: ', num_nodes

    allocate(col_list(1), row_list(1), stat=st)
    if(st.ne.0) goto 10

    ! Set up inverse permutation
    ! deallocate (invp,stat=st)
    ! allocate (invp(n),stat=st)
    ! if(st.ne.0) go to 10

    ! do j = 1, n
    !    invp(order(j)) = j
    ! end do
    
    ! TODO to be done at analyse?
    ! init factorization
    ! call spllt_factorization_init(keep, data)

    deallocate (keep%lfact,stat=st)
    allocate (keep%lfact(keep%nbcol),stat=st)
    if(st.ne.0) go to 10

    ! allocate blocks 
    deallocate (fdata%bc,stat=st)
    allocate(fdata%bc(keep%final_blk),stat=st)
    if(st.ne.0) go to 10

    ! call spllt_init_lfact(keep, fdata)

    !
    ! Copy matrix values across from a into keep%lfact
    !
    allocate(map(n),stat=st)
    if(st.ne.0) go to 10
    
    ! allocate array for colomn mapping beween snode and ancestors
    ! allocate(colmap(n),stat=st)
    ! if(st.ne.0) go to 10

    ! init facto    

#if defined(SPLLT_USE_STARPU)
    
    ! register workspace handle
    call starpu_f_vector_data_register(fdata%workspace%hdl, -1, c_null_ptr, &
         & int(keep%maxmn*keep%maxmn, kind=c_int), int(wp,kind=c_size_t))
#else
    allocate(fdata%workspace%c(keep%maxmn*keep%maxmn), stat=st)
    if(st.ne.0) goto 10
#endif

    ! call system_clock(start_cpya2l_t, rate_cpya2l_t)
    ! call copy_a_to_l(n,num_nodes,val,map,keep)
    do snode = 1, num_nodes ! loop over nodes
       call starpu_f_void_data_register(fdata%nodes(snode)%hdl)

       ! activate node: allocate factors, register handles
       call spllt_activate_node(snode, keep, fdata)

       ! init node
       call spllt_init_node_task(fdata%nodes(snode), n, val, map, keep, huge(1))

    end do
    ! call system_clock(stop_cpya2l_t)
! #if defined(SPLLT_USE_STARPU)
!     call starpu_f_task_wait_for_all()
! #endif

#if defined(SPLLT_USE_STARPU) && defined(SPLLT_STARPU_NOSUB)
    call starpu_f_pause()
#endif

    call system_clock(stop_setup_t)
    ! write(*,'("[>] [spllt_stf_factorize] cpy a2l time: ", es10.3, " s")') (stop_cpya2l_t - start_cpya2l_t)/real(rate_cpya2l_t)
    write(*,'("[>] [spllt_stf_factorize]   setup time: ", es10.3, " s")') (stop_setup_t - start_setup_t)/real(rate_setup_t)
    call system_clock(stf_start_t, stf_rate_t)

    ! factorize nodes
    ! blk = 1
    do snode = 1, num_nodes

       node => keep%nodes(snode)
       
       ! write(*,*) "---------- node ----------"
       ! write(*,*) 'snode: ', snode 

       ! task priority
       prio = (num_nodes - snode + 1)*4

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
       
       ! init column mapping
       ! colmap = 0
       
       ! first block in node
       dblk = node%blk_sa

       ! Loop over the block columns in node. 
       do kk = 1, nc

          ! A_kk          
          ! bc_kk => keep%blocks(dblk)
          bc_kk => fdata%bc(dblk)

          call spllt_factorize_block_task(fdata%nodes(snode), bc_kk, keep%lfact, prio+3)

! #if defined(SPLLT_USE_STARPU)
         ! call starpu_f_task_wait_for_all()
! #endif
          ! loop over the row blocks (that is, loop over blocks in block col)
          do ii = kk+1,nr
          ! do blk = dblk+1, dblk+sz-1
             
             ! A_mk
             blk = dblk+ii-kk
             ! bc_ik => keep%blocks(blk)
             bc_ik => fdata%bc(blk)

             call spllt_solve_block_task(bc_kk, bc_ik, keep%lfact,prio+2)
          end do

! #if defined(SPLLT_USE_STARPU)
!          call starpu_f_task_wait_for_all()
! #endif

          do jj = kk+1,nc
             
             ! L_jk
             blk2 = dblk+jj-kk
             bc_jk => fdata%bc(blk2)
             ! bc_jk => keep%blocks(blk2)

             do ii = jj,nr

                ! L_ik
                blk1 = dblk+ii-kk                
                ! bc_ik => keep%blocks(blk1)
                bc_ik => fdata%bc(blk1)
                
                ! A_ij
                ! blk = get_dest_block(keep%blocks(blk1), keep%blocks(blk2))
                blk = get_dest_block(keep%blocks(blk2), keep%blocks(blk1))
                bc_ij => fdata%bc(blk)
                ! bc_ij => keep%blocks(blk)
                
                call spllt_update_block_task(bc_ik, bc_jk, bc_ij, keep%lfact, prio+1)
             end do
          end do
          
! #if defined(SPLLT_USE_STARPU)
!           call starpu_f_task_wait_for_all()
! #endif

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
                   ! write(*,*) "i:",i," k:",k
                   blk = bc_kk%blk%id + (i-1)/s_nb - (kk-1)
                   ! bc => keep%blocks(blk) ! current block in node
                   bc => fdata%bc(blk) ! current block in node

                   if(k.ne.ii) then
                      a_blk = a_dblk + ii - cb
                      ! a_bc => keep%blocks(a_blk)
                      a_bc => fdata%bc(a_blk)

                      ! rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
                      ! rsrc(2) = (i - ilast)*blocks(blk)%blkn
                      ! write(*,*) "kk: ", kk, ", dblk: ", dblk, ", blk: ", blk, ", a_num: ", a_num, ", a_blk: ", a_blk
                      ! write(*,*) "ilast: ", ilast, ", i: ", i, ", rsrc(2): ", rsrc(2)

                      call spllt_update_between_task( &
                           ! & bc, &
                           & bc_kk, &
                           & node, a_bc, fdata%nodes(a_num), &
                           ! & csrc, rsrc, &
                           & cptr, cptr2, ilast, i-1, &
                           & row_list, col_list, fdata%workspace, &
                           & keep%lfact, keep%blocks, fdata%bc, &
                           & control, prio)

                      ii = k
                      ilast = i ! Update start of current block
                   end if
                end do
                ! i = min(i,numrow)
                a_blk = a_dblk + ii - cb
                ! a_bc => keep%blocks(a_blk)
                a_bc => fdata%bc(a_blk)

                ! rsrc(1) = 1 + (ilast-(kk-1)*s_nb-1)*blocks(blk)%blkn
                ! rsrc(2) = (i - ilast)*blocks(blk)%blkn
                ! write(*,*)"i: ", i 
                ! write(*,*)"size lcol", size(keep%lfact(keep%blocks(blk)%bcol)%lcol), ", rsrc(1)+rsrc(2)-1: ", rsrc(1)+rsrc(2)-1
                ! write(*,*) "kk: ", kk, ", dblk: ", dblk, ", blk: ", blk, ", a_num: ", a_num, ", a_blk: ", a_blk
                ! write(*,*) "ilast: ", ilast, ", i: ", i, ", rsrc(2): ", rsrc(2)

                call spllt_update_between_task( &
                     ! & bc, &
                     & bc_kk, &
                     & node, a_bc, fdata%nodes(a_num), &
                     ! & csrc, rsrc, &
                     & cptr, cptr2, ilast, i-1, &
                     & row_list, col_list, fdata%workspace, &
                     & keep%lfact, keep%blocks, fdata%bc, &
                     & control, prio)
                                         
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
    ! unregister workspace handle
    call starpu_f_data_unregister_submit(fdata%workspace%hdl)

    ! unregister data handles
    call spllt_deinit_task(keep, fdata)

    ! do snode = 1, num_nodes
    !    call starpu_f_void_unregister_submit(fdata%nodes(snode)%hdl)        
    ! end do

    call system_clock(stf_stop_t)
    write(*,'("[>] [spllt_stf_factorize] task insert time: ", es10.3, " s")') (stf_stop_t - stf_start_t)/real(stf_rate_t)

#if defined(SPLLT_STARPU_NOSUB)
    call system_clock(start_nosub_t, rate_nosub_t)
    call starpu_f_resume()
    ! wait for task completion
    call starpu_f_task_wait_for_all()
#endif
#endif

    call system_clock(stop_t)
    ! write(*,'("[>] [spllt_stf_factorize] time: ", es10.3, " s")') (stop_t - start_t)/real(rate_t)
#if defined(SPLLT_USE_STARPU) && defined(SPLLT_STARPU_NOSUB)
    write(*,'("[>] [spllt_stf_factorize] nosub time: ", es10.3, " s")') (stop_t - start_nosub_t)/real(rate_nosub_t)
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

  ! subroutine spllt_factorization_init(keep, data)
  !   use hsl_ma87_double
  !   implicit none

  !   type(MA87_keep), target, intent(inout) :: keep 
  !   type(spllt_data_type), intent(inout) :: data

  !   ! Allocate factor storage (in keep%lfact)
  !   ! TODO put block in data structure directly?    
  !   ! init lfactor
  !   call spllt_init_lfact(keep)
        
  !   return
  ! end subroutine spllt_factorization_init

  ! initialize (allocate) L factors
  ! StarPU: register data handles (block handles) in StarPU
  subroutine spllt_init_lfact(keep, pbl)
    use hsl_ma87_double
#if defined(SPLLT_USE_STARPU)
    use  starpu_f_mod
#endif
    implicit none

    type(MA87_keep), target, intent(inout) :: keep 
    type(spllt_data_type), intent(inout) :: pbl

    type(node_type), pointer :: node ! node in the atree    
    integer(long) :: blk, dblk
    integer :: nbcol, l_nb, sz, sa, en
    integer :: blkm, blkn, size_bcol
    integer :: snode, num_nodes
    integer :: i
    integer :: st ! stat parameter
    integer :: ptr

    num_nodes = keep%info%num_nodes

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
          allocate (keep%lfact(nbcol)%lcol(size_bcol),stat=st)
          ! TODO trace error
          ! TODO merge with previous loop?
          ! register blocks hanldes in StarPU

          ptr = 1
          do blk = dblk, dblk+sz-1
             blkm = keep%blocks(blk)%blkm
             blkn = keep%blocks(blk)%blkn

             pbl%bc(blk)%blk => keep%blocks(blk) 
#if defined(SPLLT_USE_STARPU)
             call starpu_matrix_data_register(pbl%bc(blk)%hdl, pbl%bc(blk)%mem_node, &
                  & c_loc(keep%lfact(nbcol)%lcol(ptr)), blkm, blkm, blkn, &
                  & int(wp,kind=c_size_t))
#endif
             pbl%bc(blk)%c => keep%lfact(nbcol)%lcol(ptr:ptr+blkm*blkn-1)
             ptr = ptr + blkm*blkn
          end do
          sz = sz - 1
       end do       
    end do

  end subroutine spllt_init_lfact

#if defined(SPLLT_USE_STARPU)
  ! deinitialize factorization
  ! StarPU: unregister data handles (block handles) in StarPU
  subroutine spllt_deinit_task(keep, pbl)
    use hsl_ma87_double
    use  starpu_f_mod
    implicit none

    type(MA87_keep), target, intent(inout) :: keep 
    type(spllt_data_type), intent(inout) :: pbl

    integer :: i
    integer :: snode, num_nodes
    type(node_type), pointer :: node ! node in the atree    
    integer(long) :: blk, dblk
    integer :: l_nb, sz, sa, en

    num_nodes = keep%info%num_nodes
    blk = 1

    do snode = 1, num_nodes

       node => keep%nodes(snode)

       l_nb = node%nb
       sz = (size(node%index) - 1) / l_nb + 1
       sa = node%sa
       en = node%en

       do i = sa, en, l_nb
          dblk = blk
          do blk = dblk, dblk+sz-1
             call starpu_f_data_unregister_submit(pbl%bc(blk)%hdl)
          end do
          sz = sz - 1
       end do
    end do

  end subroutine spllt_deinit_task
#endif

end module spllt_stf_mod
