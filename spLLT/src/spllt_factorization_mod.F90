module spllt_factorization_mod
  use spllt_data_mod
  implicit none

#if defined(SPLLT_USE_STARPU)
  ! factorize node StarPU task insert
  interface
     subroutine spllt_insert_factorize_node_task_c(node_hdl, cnode_hdls, nchild, &
          & map_hdl, snode, fdata, keep, control, prio) bind(C)
       use iso_c_binding
       type(c_ptr), value     :: node_hdl, map_hdl
       type(c_ptr)            :: cnode_hdls(*)
       type(c_ptr), value     :: snode, fdata, keep, control
       integer(c_int), value  :: prio, nchild
     end subroutine spllt_insert_factorize_node_task_c
  end interface

  interface
     subroutine spllt_starpu_codelet_unpack_args_factorize_node(cl_arg, &
          & snode, fdata, keep, control) bind(C)
       use iso_c_binding
       type(c_ptr), value :: cl_arg 
       type(c_ptr), value :: snode, fdata, keep, control
     end subroutine spllt_starpu_codelet_unpack_args_factorize_node
  end interface
#endif

contains

#if defined(SPLLT_USE_STARPU)
  subroutine spllt_factorize_apply_node_task(snode, fdata, keep, control, prio)
    use iso_c_binding
    use hsl_ma87_double
    use starpu_f_mod
    use spllt_starpu_factorization_mod
    implicit none

    type(spllt_node_type), target, intent(inout)      :: snode ! node to factorize (spllt)    
    type(spllt_data_type), target, intent(inout)      :: fdata
    type(MA87_keep), target, intent(inout)            :: keep 
    type(MA87_control), target, intent(in)            :: control 
    integer, intent(in)                               :: prio

    type(c_ptr) :: snode_c
    type(c_ptr) :: fdata_c
    type(c_ptr) :: keep_c
    type(c_ptr) :: control_c

    type(node_type), pointer        :: node
    type(spllt_node_type), pointer  :: cnode 
    integer :: nchild, i, c
    type(c_ptr), allocatable, target :: cnode_handles(:)

    snode_c = c_loc(snode)
    fdata_c = c_loc(fdata)
    keep_c = c_loc(keep)
    control_c = c_loc(control)

    node => snode%node
    nchild = node%nchild

    allocate(cnode_handles(nchild))
    do i = 1, nchild
       c = node%child(i)
       cnode => fdata%nodes(c)
       cnode_handles(i) = cnode%hdl2 
    end do

    call spllt_insert_factorize_node_task_c(snode%hdl2, cnode_handles, &
         & int(nchild, kind=c_int), fdata%map%hdl, snode_c, fdata_c, keep_c, control_c, &
         & prio)

    deallocate(cnode_handles)

  end subroutine spllt_factorize_apply_node_task

  ! factorize node StarPU task
  subroutine spllt_starpu_factorize_node_cpu_func(buffers, cl_arg) bind(C)
    use iso_c_binding
    use spllt_data_mod
    use hsl_ma87_double
    use spllt_kernels_mod
    implicit none

    type(c_ptr), value        :: cl_arg
    type(c_ptr), value        :: buffers

    type(c_ptr), target            :: snode_c, fdata_c, keep_c, control_c
    type(c_ptr), target            :: map_c
    type(spllt_node_type),pointer  :: snode
    type(spllt_data_type), pointer :: fdata
    type(ma87_keep), pointer       :: keep    
    type(MA87_control), pointer    :: control 
    integer, pointer               :: map(:)
    integer, target :: n

    call spllt_starpu_codelet_unpack_args_factorize_node(cl_arg, &
         & c_loc(snode_c), c_loc(fdata_c), &
         & c_loc(keep_c), c_loc(control_c)) 

    call c_f_pointer(snode_c, snode)    
    call c_f_pointer(fdata_c, fdata)    
    call c_f_pointer(keep_c, keep)    
    call c_f_pointer(control_c, control)    

    call starpu_f_get_buffer(buffers, 0, c_loc(map_c), c_loc(n))
    call c_f_pointer(map_c, map, (/n/))

    ! allocate(map(n))

    ! write(*,*)"num", snode%num

    call spllt_factorize_apply_node(snode, map, fdata, keep, control)

    ! deallocate(map)

    ! write(*,*)"map() :", map(n)
    ! write(*,*)"n :", n
    ! write(*,*)"n :", keep%n
    ! write(*,*)"final_blk :", keep%final_blk
    ! write(*,*)"min_width_blas", control%min_width_blas

  end subroutine spllt_starpu_factorize_node_cpu_func
#endif

  ! initialize factorization
  ! allocate map array
  ! allocate workspaces: fdata%workspace, fdata%row_list, fdata%col_list 
  subroutine spllt_factorization_init(fdata, map, keep)
    use hsl_ma87_double
    implicit none

    type(spllt_data_type), target :: fdata
    integer, dimension(:), pointer ::  map
    type(MA87_keep), target, intent(inout) :: keep 

    integer :: n ! order of the system
    integer :: st ! stat parameter
    
    n = keep%n

    deallocate (keep%lfact,stat=st)
    allocate (keep%lfact(keep%nbcol),stat=st)

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
         & int(keep%maxmn*keep%maxmn, kind=c_int), int(wp,kind=c_size_t))

    ! register col_list and row_list workspaces
    ! register row_list handle
    call starpu_f_vector_data_register(fdata%row_list%hdl, -1, c_null_ptr, &
         & int(keep%maxmn, kind=c_int), int(c_int, kind=c_size_t))

    ! register col_list handle
    call starpu_f_vector_data_register(fdata%col_list%hdl, -1, c_null_ptr, &
         & int(keep%maxmn, kind=c_int), int(c_int, kind=c_size_t))

#if defined(SPLLT_USE_NESTED_STF)
    call starpu_f_vector_data_register(fdata%map%hdl, -1, c_null_ptr, &
         & int(n, kind=c_int), int(c_int, kind=c_size_t))
#endif

#elif defined(SPLLT_USE_OMP)

    nt = 1
!$  nt = omp_get_num_threads()
    allocate(fdata%workspace(0:nt-1))
    allocate(fdata%row_list(0:nt-1))
    allocate(fdata%col_list(0:nt-1))

    do i=0,nt-1
       allocate(fdata%workspace(i)%c(keep%maxmn*keep%maxmn), stat=st)
       ! if(st.ne.0) goto 10
       allocate(fdata%row_list(i)%c(keep%maxmn), stat=st)
       ! if(st.ne.0) goto 10
       allocate(fdata%col_list(i)%c(keep%maxmn), stat=st)
       ! if(st.ne.0) goto 10
    end do
#else
    allocate(fdata%workspace%c(keep%maxmn*keep%maxmn), stat=st)
    ! if(st.ne.0) goto 10
    allocate(fdata%row_list%c(keep%maxmn), stat=st)
    ! if(st.ne.0) goto 10
    allocate(fdata%col_list%c(keep%maxmn), stat=st)
    ! if(st.ne.0) goto 10
#endif

  end subroutine spllt_factorization_init

! #if defined(SPLLT_USE_STARPU)
!   subroutine spllt_factorization_fini_task(fdata, map, keep)
!     implicit none

!     type(spllt_data_type), target :: fdata
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
  subroutine spllt_factorization_fini(fdata, map, keep)
    use hsl_ma87_double
    use spllt_factorization_task_mod
    implicit none

    type(spllt_data_type), target :: fdata
    integer, dimension(:), pointer ::  map
    type(MA87_keep), target, intent(inout) :: keep 

    integer :: st ! stat parameter
    
#if defined(SPLLT_USE_STARPU)
    
    ! unregister workspace handle
    call starpu_f_data_unregister_submit(fdata%workspace%hdl)

    call starpu_f_data_unregister_submit(fdata%row_list%hdl)
    call starpu_f_data_unregister_submit(fdata%col_list%hdl)

#if defined(SPLLT_USE_NESTED_STF)
    call starpu_f_data_unregister_submit(fdata%map%hdl)
#endif

#endif
    
    ! unregister data handle
#if defined(SPLLT_USE_STARPU)

#ifndef SPLLT_USE_GPU

    call spllt_data_unregister_task(keep, fdata)
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

  subroutine spllt_factorize_node(snode, fdata, keep)
    use spllt_data_mod
    use hsl_ma87_double
    use spllt_factorization_task_mod
    implicit none

    type(spllt_node_type), target, intent(inout)        :: snode ! node to factorize (spllt)    
    type(spllt_data_type), target, intent(inout)        :: fdata
    type(MA87_keep), target, intent(inout)              :: keep 

    type(node_type), pointer :: node ! node to factorize (hsl_ma87)
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
    
    node => snode%node

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
       call spllt_factorize_block_task(fdata, snode, bc_kk, keep%lfact, prio+3)

       ! loop over the row blocks (that is, loop over blocks in block col)
       do ii = kk+1,nr
          ! do ii = nr,kk+1,-1             
          ! A_mk
          blk = dblk+ii-kk
          ! bc_ik => keep%blocks(blk)
          bc_ik => fdata%bc(blk)
          ! write(*,*)"ii: ", ii
          call spllt_solve_block_task(fdata, bc_kk, bc_ik, keep%lfact,prio+2)
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
             blk = get_dest_block(keep%blocks(blk2), keep%blocks(blk1))
             bc_ij => fdata%bc(blk)
             call spllt_update_block_task(fdata, bc_ik, bc_jk, bc_ij, keep%lfact, prio+1)

          end do
       end do

       ! move to next block column in snode
       dblk = keep%blocks(dblk)%last_blk + 1
       ! numrow = numrow - s_nb
    end do

  end subroutine spllt_factorize_node

  ! node factorization
  ! Submit the DAG for the factorization of a node
  subroutine spllt_factorize_apply_node(snode, map, fdata, keep, control)
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

    ! perform blocked cholesky factorizaton on node
    call spllt_factorize_node(snode, fdata, keep)

    ! update between
    node => snode%node

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
  end subroutine spllt_factorize_apply_node

end module spllt_factorization_mod
