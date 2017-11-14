module spllt_factorization_mod
  use spllt_data_mod
  implicit none

#if defined(SPLLT_USE_STARPU)
  ! factorize node StarPU task insert
  interface
     subroutine spllt_insert_factorize_node_task_c(node_hdl, cnode_hdls, nchild, &
          & map_hdl, snode, fdata, control, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value     :: node_hdl, map_hdl
       type(c_ptr)            :: cnode_hdls(*)
       type(c_ptr), value     :: snode, fdata, control
       integer(c_int), value  :: prio, nchild
     end subroutine spllt_insert_factorize_node_task_c
  end interface

  interface
     subroutine spllt_starpu_codelet_unpack_args_factorize_node(cl_arg, &
          & snode, fdata, control) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg 
       type(c_ptr), value :: snode, fdata, control
     end subroutine spllt_starpu_codelet_unpack_args_factorize_node
  end interface
#endif

contains

  ! Apply updates from subtree's generated element to its ancestors
  ! This routine aquires all locks as needed, and decrements dependencies
  ! upon release of relevant locks.
  ! NB: This is basically an immediate acting variant of add_between_updates()
  subroutine spllt_subtree_apply_buffer(root, buffer, nodes, blocks, lfact, map, &
       fdata)
    use spllt_data_mod
    use spllt_factorization_task_mod
    implicit none

    integer, intent(in) :: root ! root of subtree
    ! real(wp), dimension(*), intent(in) :: buffer ! generated element
    type(spllt_block), intent(in) :: buffer

    type(spllt_node), dimension(-1:), intent(in) :: nodes
    type(spllt_block), dimension(*), intent(inout) :: blocks
    type(lfactor), dimension(*), intent(inout) :: lfact
    integer, dimension(:), intent(inout) :: map ! Workarray to hold map from row
    ! indices to block indices in ancestor node. 
    type(spllt_fdata_type) :: fdata

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

                call spllt_scatter_block_task(ilast, i-1, cptr, cptr2,  buffer, nodes(root), &
                     & (jb-1)*a_nb+1, (cb-1)*a_nb+1, fdata%bc(dest), nodes(anode))

                jb = k
                ilast = i ! Update start of current block
             endif
          end do
          dest = dblk+jb-cb
          rsrc(1) = ilast
          rsrc(2) = i-1
          bsa = (rsrc(1)-n-1)*lds+csrc(1)-n
          ben = (rsrc(2)-n-1)*lds+csrc(2)-n

          call spllt_scatter_block_task(ilast, i-1, cptr, cptr2,  buffer, nodes(root), &
               & (jb-1)*a_nb+1, (cb-1)*a_nb+1, fdata%bc(dest), nodes(anode))

          ! Move cptr down, ready for next block column of anode
          cptr = cptr2 + 1
       end do bcols

       ! Move up the tree
       anode = nodes(anode)%parent
    end do
  end subroutine spllt_subtree_apply_buffer


  subroutine spllt_subtree_factorize_apply(root, fdata, val, cntl, map, buffer)
    use spllt_data_mod
    use spllt_kernels_mod
    use spllt_factorization_task_mod
    implicit none

    integer, intent(in) :: root ! Index of root node     
    type(spllt_fdata_type), target, intent(inout) :: fdata ! Factorize data
    real(wp), dimension(:), intent(in) :: val ! User's matrix values
    type(spllt_cntl), intent(inout) :: cntl
    integer, pointer, intent(inout) :: map(:)
    type(spllt_block), target, intent(inout) :: buffer ! update_buffer workspace

    type(spllt_node), pointer :: node ! node in the atree    
    integer :: m, n, b_sz
    type(spllt_block), pointer :: buf

    node => fdata%nodes(root)
    m = size(node%index)
    n = fdata%nodes(root)%en - fdata%nodes(root)%sa + 1
    b_sz = m-n

#if defined(SPLLT_USE_STARPU)

    ! Create a temporary buffer for the update operation. This is done
    ! by registering a StarPU handle associated to a NULL pointer
    ! which means that StarPU is in charge of allocating/deallocating
    ! memory for this buffer
    call starpu_f_matrix_data_register(fdata%nodes(root)%buffer%hdl, & 
         -1, c_null_ptr, int(b_sz, kind=c_int), int(b_sz, kind=c_int), &
         int(b_sz, kind=c_int), int(wp, kind=c_size_t))
    
#elif defined(SPLLT_USE_OMP)

    allocate(fdata%nodes(root)%buffer%c(b_sz**2)) ! FIXME deallocate memory after apply
    
#else
    ! FIXME use fdata%nodes(root)%buffer
    if (size(buffer%c).lt.(m-n)**2) then
       deallocate(buffer%c)
       allocate(buffer%c((m-n)**2))
    end if
    buf => buffer
#endif

    ! subtree factorization task

    ! call system_clock(subtree_start_t, subtree_rate_t)
    call spllt_subtree_factorize_task(root, fdata, val, fdata%nodes(root)%buffer, cntl, map)
    ! call system_clock(subtree_stop_t)
    ! write(*,'("[>] [spllt_stf_factorize] facto subtree: ", es10.3, " s")') &
         ! & (subtree_stop_t - subtree_start_t)/real(subtree_rate_t)

    ! Expand generated element out to ancestors
    ! call system_clock(subtree_start_t, subtree_rate_t)
    call spllt_subtree_apply_buffer(root, fdata%nodes(root)%buffer, fdata%nodes, fdata%bc, &
         fdata%lfact, map, fdata)


#if defined(SPLLT_USE_STARPU)   
    call starpu_f_data_unregister_submit(fdata%nodes(root)%buffer%hdl)
    ! deallocate(buf%c)
    ! deallocate(buf)
#endif

  end subroutine spllt_subtree_factorize_apply

#if defined(SPLLT_USE_STARPU)

  ! Factorize node task submission routine for StarPU
  subroutine spllt_factorize_apply_node_task(node, fdata, cntl, prio)
    use iso_c_binding
    use spllt_data_mod
    use starpu_f_mod
    use spllt_starpu_factorization_mod
    implicit none

    type(spllt_node), target, intent(inout) :: node ! supernode to be factorized 
    type(spllt_fdata_type), target, intent(inout) :: fdata ! factorization data
    type(spllt_cntl), target, intent(in) :: cntl ! options
    integer, intent(in) :: prio ! task priority for scheduler

    ! C pointers
    type(c_ptr) :: node_c
    type(c_ptr) :: fdata_c
    type(c_ptr) :: cntl_c

    type(spllt_node), pointer :: cnode ! child node 
    integer :: nchild, i, c
    type(c_ptr), allocatable, target :: cnode_handles(:) ! children node hanldes

    node_c = c_loc(node)
    fdata_c = c_loc(fdata)
    ! keep_c = c_loc(keep)
    cntl_c = c_loc(cntl)

    ! Gather the handles from the children nodes and put them in the
    ! cnode_handles array.
    nchild = node%nchild
    allocate(cnode_handles(nchild))
    do i = 1, nchild
       c = node%child(i)
       cnode => fdata%nodes(c)
       cnode_handles(i) = cnode%hdl2 
    end do

    call spllt_insert_factorize_node_task_c(node%hdl2, cnode_handles, &
         & int(nchild, kind=c_int), fdata%map%hdl, node_c, fdata_c, cntl_c, &
         & prio)

    deallocate(cnode_handles)

  end subroutine spllt_factorize_apply_node_task

  ! Factorize node CPU kernel for StarPU 
  subroutine spllt_starpu_factorize_node_cpu_func(buffers, cl_arg) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use spllt_kernels_mod
    implicit none

    type(c_ptr), value :: cl_arg
    type(c_ptr), value :: buffers

    type(c_ptr), target:: node_c, fdata_c, cntl_c
    type(c_ptr), target :: map_c
    type(spllt_node),pointer :: node
    type(spllt_fdata_type), pointer :: fdata
    type(spllt_cntl), pointer :: cntl 
    integer, pointer :: map(:)
    integer, target :: n

    call spllt_starpu_codelet_unpack_args_factorize_node(cl_arg, &
         & c_loc(node_c), c_loc(fdata_c), c_loc(cntl_c)) 

    call c_f_pointer(node_c, node)
    call c_f_pointer(fdata_c, fdata)
    call c_f_pointer(cntl_c, cntl)

    call starpu_f_get_buffer(buffers, 0, c_loc(map_c), c_loc(n))
    call c_f_pointer(map_c, map, (/n/))

    call spllt_factorize_apply_node(node, map, fdata, cntl)

  end subroutine spllt_starpu_factorize_node_cpu_func
#endif

  ! initialize factorization
  ! allocate map array
  ! allocate workspaces: fdata%workspace, fdata%row_list, fdata%col_list 
  subroutine spllt_factorization_init(fdata, map)
    use spllt_data_mod
#if defined(SPLLT_USE_OMP)
!$  use omp_lib
#endif
    implicit none

    type(spllt_fdata_type), target :: fdata
    integer, dimension(:), pointer ::  map
    ! type(spllt_keep), target, intent(inout) :: keep 

    integer :: n ! order of the system
    integer :: st ! stat parameter
#if defined(SPLLT_USE_OMP)
    integer :: i ! thread index
    integer :: nt ! num threads
#endif
   
    n = fdata%n

    deallocate (fdata%lfact,stat=st)
    allocate (fdata%lfact(fdata%nbcol),stat=st)

#ifndef SPLLT_USE_NESTED_STF
    ! when using the nested STF model we let StarPU manage the map
    ! array
    allocate(map(n),stat=st)
    ! TODO error managment
    ! if(st.ne.0) go to 10
#endif    

#if defined(SPLLT_USE_STARPU)

    ! register workspace handle
    call starpu_f_vector_data_register(fdata%workspace%hdl, -1, c_null_ptr, &
         & int(fdata%maxmn*fdata%maxmn, kind=c_int), int(wp, kind=c_size_t))

    ! register col_list and row_list workspaces
    ! register row_list handle
    call starpu_f_vector_data_register(fdata%row_list%hdl, -1, c_null_ptr, &
         & int(fdata%maxmn, kind=c_int), int(c_int, kind=c_size_t))

    ! register col_list handle
    call starpu_f_vector_data_register(fdata%col_list%hdl, -1, c_null_ptr, &
         & int(fdata%maxmn, kind=c_int), int(c_int, kind=c_size_t))

    call starpu_f_vector_data_register(fdata%map%hdl, -1, c_null_ptr, &
         & int(n, kind=c_int), int(c_int, kind=c_size_t))

#elif defined(SPLLT_USE_OMP)

    nt = 1
!$  nt = omp_get_num_threads()
    allocate(fdata%workspace(0:nt-1)) ! workspace
    allocate(fdata%row_list(0:nt-1)) ! row list workspace 
    allocate(fdata%col_list(0:nt-1)) ! col list workspace
    allocate(fdata%map(0:nt-1)) ! map workspace

    do i = 0, nt-1
       allocate(fdata%workspace(i)%c(fdata%maxmn*fdata%maxmn), stat=st)
       ! if(st.ne.0) goto 10
       allocate(fdata%row_list(i)%c(fdata%maxmn), stat=st)
       ! if(st.ne.0) goto 10
       allocate(fdata%col_list(i)%c(fdata%maxmn), stat=st)
       ! if(st.ne.0) goto 10
       allocate(fdata%map(i)%c(n), stat=st)
    end do
#else
    allocate(fdata%workspace%c(fdata%maxmn*fdata%maxmn), stat=st)
    ! if(st.ne.0) goto 10
    allocate(fdata%row_list%c(fdata%maxmn), stat=st)
    ! if(st.ne.0) goto 10
    allocate(fdata%col_list%c(fdata%maxmn), stat=st)
    ! if(st.ne.0) goto 10
#endif

  end subroutine spllt_factorization_init

! #if defined(SPLLT_USE_STARPU)
!   subroutine spllt_factorization_fini_task(fdata, map, keep)
!     implicit none

!     type(spllt_fdata_type), target :: fdata
!     integer, dimension(:), pointer ::  map
!     type(MA87_keep), target, intent(inout) :: keep 
    
!     type(c_ptr) :: fdata_c
!     type(c_ptr) :: keep_c

!     fdata_c = c_loc(fdata)
!     keep_c = c_loc(keep)

    

!   end subroutine spllt_factorization_fini_task
! #endif


  ! deinitialize factorization
  ! deallocate map array
  ! deallocate workspaces
  subroutine spllt_factorization_fini(fdata, map, adata)
    use spllt_data_mod
    use spllt_factorization_task_mod
    implicit none

    type(spllt_fdata_type), target :: fdata
    integer, dimension(:), pointer ::  map
    ! type(spllt_keep), target, intent(inout) :: keep 
    type(spllt_adata_type), target, intent(in) :: adata

    integer :: st ! stat parameter
    
#if defined(SPLLT_USE_STARPU)
    
    ! unregister workspace handle
    call starpu_f_data_unregister_submit(fdata%workspace%hdl)

    call starpu_f_data_unregister_submit(fdata%row_list%hdl)
    call starpu_f_data_unregister_submit(fdata%col_list%hdl)
    call starpu_f_data_unregister_submit(fdata%map%hdl)

#endif
    
    ! unregister data handle
#if defined(SPLLT_USE_STARPU)

#ifndef SPLLT_USE_GPU

    call spllt_data_unregister_task(fdata, adata)
#endif

    ! do snode = 1, num_nodes
    !    call starpu_f_void_unregister_submit(fdata%nodes(snode)%hdl)        
    ! end do
#endif

#ifndef SPLLT_USE_NESTED_STF
    ! when using the nested STF model we let StarPU manage the map
    ! array
    deallocate(map, stat=st)
    ! TODO error managment
    ! if(st.ne.0) go to 10
#endif    

  end subroutine spllt_factorization_fini

  subroutine spllt_factorize_node(node, fdata)
    use spllt_data_mod
    use spllt_factorization_task_mod
    implicit none

    type(spllt_node), target, intent(inout) :: node ! node to factorize (spllt)    
    type(spllt_fdata_type), target, intent(inout) :: fdata
    ! type(spllt_keep), target, intent(inout)              :: keep 

    ! type(spllt_node), pointer :: node ! node to factorize (hsl_ma87)
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
    type(spllt_block), pointer :: bc_kk, bc_ik, bc_jk, bc_ij ! block pointers
    integer(long) :: blk, blk1, blk2 ! block id
    
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
       ! write(*,*)"kk: ", kk
       ! A_kk
       bc_kk => fdata%bc(dblk)
       call spllt_factorize_block_task(fdata, node, bc_kk, fdata%lfact, prio+3)

       ! loop over the row blocks (that is, loop over blocks in block col)
       do ii = kk+1,nr
          ! do ii = nr,kk+1,-1             
          ! A_mk
          blk = dblk+ii-kk
          ! bc_ik => keep%blocks(blk)
          bc_ik => fdata%bc(blk)
          ! write(*,*)"ii: ", ii
          call spllt_solve_block_task(fdata, bc_kk, bc_ik, fdata%lfact,prio+2)
       end do

       ! perform update operations within node
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
             blk = get_dest_block(fdata%bc(blk2), fdata%bc(blk1))
             bc_ij => fdata%bc(blk)
             call spllt_update_block_task(fdata, bc_ik, bc_jk, bc_ij, fdata%lfact, prio+1)

          end do
       end do

       ! move to next block column in snode
       dblk = fdata%bc(dblk)%last_blk + 1
       ! numrow = numrow - s_nb
    end do

  end subroutine spllt_factorize_node

  ! node factorization
  ! Submit the DAG for the factorization of a node
  subroutine spllt_factorize_apply_node(node, map, fdata, cntl)
    use spllt_data_mod
    use spllt_kernels_mod
    use spllt_factorization_task_mod
    implicit none

    type(spllt_node), target, intent(inout) :: node ! node to factorize (spllt)    
    integer, dimension(:), pointer, intent(inout) :: map
    type(spllt_fdata_type), target, intent(inout) :: fdata
    ! type(spllt_keep), target, intent(inout) :: keep 
    type(spllt_cntl), intent(in) :: cntl

    ! type(spllt_node), pointer :: node ! node to factorize (hsl_ma87)
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
    type(spllt_block), pointer :: bc_kk, bc_ik, bc_jk, bc_ij ! block pointers
    integer(long) :: blk, blk1, blk2 ! block id

    ! update between
    type(spllt_node), pointer :: anode ! ancestor node in the atree
    ! type(spllt_node), pointer :: a_node ! ancestor node in the atree
    ! locate source blocks
    integer :: a_num ! ancestor id
    integer :: cptr  ! Position in snode of the first row 
    ! matching a column of the current block column of anode.
    integer :: cptr2  ! Position in snode of the last row 
    ! matching a column of the current block column of anode.
    logical :: map_done
    integer :: i, ilast
    type(spllt_block), pointer :: bc, a_bc ! node in the atree
    integer :: cb, jb
    integer :: jlast ! Last column in the cb-th block column of anode
    integer :: k
    integer(long) :: a_dblk, a_blk ! id of block in scol containing row 

    ! print *, 'factorize_apply_node' 
    ! write(*,*) "---------- node ----------"
    ! write(*,*) 'snode: ', snode 

    ! perform blocked cholesky factorizaton on node
    call spllt_factorize_node(node, fdata)
    ! return
    ! update between
    ! node => snode%node

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
       anode => fdata%nodes(a_num) 
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
             a_dblk = fdata%bc(a_dblk)%last_blk + 1
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
          if(.not.map_done) then
             call spllt_build_rowmap(anode, map)
             map_done = .true.
          end if
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
                        & node, a_bc, anode, &
                                ! & csrc, rsrc, &
                        & cptr, cptr2, ilast, i-1, &
                        & fdata%row_list, fdata%col_list, fdata%workspace, &
                        & fdata%lfact, fdata%bc, &
                        & cntl%min_width_blas, prio)

                   dblk = fdata%bc(dblk)%last_blk + 1
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
                  & node, a_bc, anode, &
                                ! & csrc, rsrc, &
                  & cptr, cptr2, ilast, i-1, &
                  & fdata%row_list, fdata%col_list, fdata%workspace, &
                  & fdata%lfact, fdata%bc, &
                  & cntl%min_width_blas, prio)

             dblk = fdata%bc(dblk)%last_blk + 1
          end do

          ! Move cptr down, ready for next block column of anode
          cptr = cptr2 + 1
       end do bcols

       a_num = anode%parent
    end do

    return
  end subroutine spllt_factorize_apply_node

end module spllt_factorization_mod
