module spllt_factorization_mod
  use spllt_data_mod
  implicit none

  ! TODO put in a proper module file
  interface
     subroutine spllt_factorize_node(snode, map, fdata, keep, control)
       use spllt_data_mod
       use hsl_ma87_double
       implicit none
       type(spllt_node_type), target        :: snode ! node to factorize (spllt)    
       integer, dimension(:), pointer       :: map
       type(spllt_data_type), target        :: fdata
       type(MA87_keep), target              :: keep 
       type(MA87_control)                   :: control 
     end subroutine spllt_factorize_node
  end interface

contains

  subroutine spllt_stf_factorize(n, ptr, row, val, order, keep, control, info, fdata, cntl)
    use spllt_data_mod
    use spllt_error_mod
    use hsl_ma87_double
    use spllt_factorization_task_mod
#if defined(SPLLT_USE_STARPU) 
    use iso_c_binding
    use starpu_f_mod
#elif defined(SPLLT_USE_OMP)
!$ use omp_lib
#endif
    use spllt_kernels_mod
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
    ! integer, dimension(:), allocatable ::  map ! allocated to have size n.

    ! used in copying entries of user's matrix a into factor storage 
    ! (keep%fact).
    integer, dimension(:), pointer ::  tmpmap

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
    ! integer :: s_nb ! set to block size of snode (keep%nodes(snode)%nb)
    integer :: nc, nr ! number of block column/row
    integer :: snode, num_nodes
    integer :: st ! stat parameter
    integer :: numrow, numcol
    integer :: ii, jj, kk
    integer :: prio
    
    ! update between variables
    ! integer :: csrc(2), rsrc(2) ! used for update_between tasks to
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
    ! real(wp) :: soln(0)
    integer :: k
    integer(long) :: a_dblk, a_blk ! id of block in scol containing row 
    ! nodes(snode)%index(cptr).
    ! integer(long) :: rb ! Index of block row in snode

#if defined(SPLLT_USE_OMP)
    integer :: nt
#endif
    ! timing
    integer :: start_t, stop_t, rate_t
    integer :: stf_start_t, stf_stop_t, stf_rate_t
    integer :: start_nosub_t, rate_nosub_t
    integer :: start_setup_t, stop_setup_t, rate_setup_t
    ! integer :: start_cpya2l_t, stop_cpya2l_t, rate_cpya2l_t
    ! call system_clock(start_t, rate_t)
    call system_clock(start_setup_t, rate_setup_t)

    ! call factorize_posdef(n, val, order, keep, control, info, 0, 0, soln)

    ! write(*,*) 'control%nb: ', control%nb

    ! shortcut
    blocks => keep%blocks
    nodes  => keep%nodes

    num_nodes = keep%info%num_nodes
    ! write(*,*) 'num_nodes: ', num_nodes

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
    ! deallocate (fdata%bc,stat=st)
    ! allocate(fdata%bc(keep%final_blk),stat=st)
    ! if(st.ne.0) go to 10

    ! call spllt_init_lfact(keep, fdata)

    !
    ! Copy matrix values across from a into keep%lfact
    !
    ! allocate(map(n),stat=st)
    ! if(st.ne.0) go to 10

#ifndef SPLLT_USE_NESTED_STF
    ! write(*,*)"TETETET"
    allocate(tmpmap(n),stat=st)
    if(st.ne.0) go to 10
#endif
    
    ! init facto    

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
       if(st.ne.0) goto 10
       allocate(fdata%row_list(i)%c(keep%maxmn), stat=st)
       if(st.ne.0) goto 10
       allocate(fdata%col_list(i)%c(keep%maxmn), stat=st)
       if(st.ne.0) goto 10
    end do
#else
    allocate(fdata%workspace%c(keep%maxmn*keep%maxmn), stat=st)
    if(st.ne.0) goto 10
    allocate(fdata%row_list%c(keep%maxmn), stat=st)
    if(st.ne.0) goto 10
    allocate(fdata%col_list%c(keep%maxmn), stat=st)
    if(st.ne.0) goto 10
#endif

    call system_clock(stf_start_t, stf_rate_t)

    do snode = 1, num_nodes ! loop over nodes
#if defined(SPLLT_USE_STARPU)
       ! TODO put in activate routine
       call starpu_f_void_data_register(fdata%nodes(snode)%hdl)
#if defined(SPLLT_USE_NESTED_STF)
       call starpu_f_void_data_register(fdata%nodes(snode)%hdl2)
#endif
#endif
       ! activate node: allocate factors, register handles
       call spllt_activate_node(snode, keep, fdata)
    end do

    call system_clock(stop_setup_t)
    write(*,'("[>] [spllt_stf_factorize]   setup and activate nodes time: ", es10.3, " s")') &
         & (stop_setup_t - start_setup_t)/real(rate_setup_t)

! #if defined(SPLLT_USE_STARPU)
!     call starpu_f_task_wait_for_all()
! #endif    
    ! call system_clock(start_cpya2l_t, rate_cpya2l_t)
    ! call copy_a_to_l(n,num_nodes,val,map,keep)
    ! write(*,*)"num_nodes: ", num_nodes
    do snode = 1, num_nodes ! loop over nodes
       prio = 5 ! max priority 
       ! prio = huge(1)
       ! init node
       call spllt_init_node_task(fdata, fdata%nodes(snode), val, keep, prio)
    end do
    ! call system_clock(stop_cpya2l_t)

! #if defined(SPLLT_USE_STARPU)
!     call starpu_f_task_wait_for_all()
! #endif

#if defined(SPLLT_USE_OMP)
!$omp taskwait
#endif

#if defined(SPLLT_USE_STARPU) && defined(SPLLT_STARPU_NOSUB)
    call starpu_f_pause()
#endif

    ! factorize nodes
    do snode = 1, num_nodes
#if defined(SPLLT_USE_NESTED_STF)
       prio = 5 ! max priority 
       ! prio = huge(1)

       call spllt_factorize_node_task(fdata%nodes(snode), fdata, keep, control, prio)          
#else
       call spllt_factorize_node(fdata%nodes(snode), tmpmap, fdata, keep, control)          
#endif
    end do

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

    call spllt_deinit_task(keep, fdata)
#endif

    ! do snode = 1, num_nodes
    !    call starpu_f_void_unregister_submit(fdata%nodes(snode)%hdl)        
    ! end do
#endif

    call system_clock(stf_stop_t)
    write(*,'("[>] [spllt_stf_factorize] task insert time: ", es10.3, " s")') (stf_stop_t - stf_start_t)/real(stf_rate_t)

#if defined(SPLLT_USE_STARPU) && defined(SPLLT_STARPU_NOSUB)
    call system_clock(start_nosub_t, rate_nosub_t)
    call starpu_f_resume()
    ! wait for task completion
    call starpu_f_task_wait_for_all()
#endif

! #if defined(SPLLT_USE_OMP)
! !$omp taskwait
! #endif

    call system_clock(stop_t)
    ! write(*,'("[>] [spllt_stf_factorize] time: ", es10.3, " s")') (stop_t - start_t)/real(rate_t)
#if defined(SPLLT_USE_STARPU) && defined(SPLLT_STARPU_NOSUB)
    write(*,'("[>] [spllt_stf_factorize] nosub time: ", es10.3, " s")') (stop_t - start_nosub_t)/real(rate_nosub_t)
#endif

10 if(st.ne.0) then
      info%flag = spllt_error_allocation
      info%stat = st
      call spllt_print_err(info%flag, context='spllt_stf_factorize',st=st)
      return
   endif

    return
  end subroutine spllt_stf_factorize

  ! left looking variant of the STF algorithm
  subroutine spllt_stf_ll_factorize(n, ptr, row, val, order, keep, control, info, fdata, cntl)
    use spllt_data_mod
    use spllt_error_mod
    use hsl_ma87_double
    use spllt_factorization_task_mod
#if defined(SPLLT_USE_STARPU) 
    use iso_c_binding
    use starpu_f_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#endif
    use spllt_kernels_mod
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

    type(block_type), dimension(:), pointer :: blocks ! block info. 
    type(node_type), dimension(:), pointer :: nodes 
    integer :: snode, num_nodes
    type(node_type), pointer     :: node, dnode ! node in the atree    

    integer, dimension(:), allocatable ::  tmpmap

    integer :: st ! stat parameter
    integer :: least_desc, d_num, d_sa, d_en, d_numcol, d_numrow, d_nb, d_nc
    ! logical :: map_done
    integer :: cptr
    integer :: cptr2
    integer :: cb, jb
    integer(long) :: dblk, blk 
    integer :: jlast
    integer :: stf_start_t, stf_stop_t, stf_rate_t
    integer :: ii, i, ilast, k
    type(spllt_bc_type), pointer :: bc, d_bc_kk
    integer(long) :: d_dblk 
    integer :: kk

    integer :: prio
    integer :: sa, en
    integer :: numcol, numrow
    integer :: nc, nr
    integer :: s_nb
    integer(long) :: blk1, blk2 
    type(spllt_bc_type), pointer :: bc_kk, bc_ik, bc_jk, bc_ij
    integer :: jj

#if defined(SPLLT_USE_OMP)
    integer :: nt
#endif

    ! shortcut
    blocks => keep%blocks
    nodes  => keep%nodes
    num_nodes = keep%info%num_nodes

    ! allocate L factors
    deallocate (keep%lfact,stat=st)
    allocate (keep%lfact(keep%nbcol),stat=st)
    if(st.ne.0) go to 9999

    allocate(tmpmap(n),stat=st)
    if(st.ne.0) go to 9999
    
    ! init facto    

#if defined(SPLLT_USE_STARPU)
    
    ! register workspace handle
    call starpu_f_vector_data_register(fdata%workspace%hdl, -1, c_null_ptr, &
         & int(keep%maxmn*keep%maxmn, kind=c_int), int(wp,kind=c_size_t))
#elif defined(SPLLT_USE_OMP)

    nt = 1
!$  nt = omp_get_num_threads()
    allocate(fdata%workspace(0:nt-1))

    do i=0,nt-1
       allocate(fdata%workspace(i)%c(keep%maxmn*keep%maxmn), stat=st)
    end do
#else
    allocate(fdata%workspace%c(keep%maxmn*keep%maxmn), stat=st)
    if(st.ne.0) goto 9999
#endif

    call system_clock(stf_start_t, stf_rate_t)

    do snode = 1, num_nodes ! loop over nodes
#if defined(SPLLT_USE_STARPU)
       ! TODO put in activate routine
       call starpu_f_void_data_register(fdata%nodes(snode)%hdl)
#endif
       ! activate node: allocate factors, register handles
       call spllt_activate_node(snode, keep, fdata)
    end do

! #if defined(SPLLT_USE_STARPU)
!     call starpu_f_task_wait_for_all()
! #endif    
    ! write(*,*)"num_nodes: ", num_nodes
    do snode = 1, num_nodes ! loop over nodes
       ! init node
       call spllt_init_node_task(fdata, fdata%nodes(snode), val, keep, huge(1))
    end do
! #if defined(SPLLT_USE_STARPU)
!     call starpu_f_task_wait_for_all()
! #endif

#if defined(SPLLT_USE_OMP)
!$omp taskwait
#endif

    do snode = 1, num_nodes

       ! write(*,*)'snode: ', snode
       
       node => keep%nodes(snode)

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
       
       !  Loop over descendent of snode
       least_desc = node%least_desc  
       
       ! build a map between row and block index
       call spllt_build_rowmap(node, tmpmap) 
       
       do d_num=least_desc, snode-1
          
          ! write(*,*)'d_num: ', d_num

          dnode => keep%nodes(d_num)
          
          ! initialize node data
          d_sa = dnode%sa
          d_en = dnode%en
          d_numcol = d_en - d_sa + 1
          d_numrow = size(dnode%index)

          ! s_nb is the size of the blocks
          d_nb = dnode%nb
          d_nc = (d_numcol-1) / d_nb + 1 
          
          cptr = 1 + d_numcol

          do cptr = cptr, d_numrow
             if(dnode%index(cptr) .ge. node%sa) exit
          end do
          if(cptr .gt. d_numrow) cycle ! finished with dnode
          
          ! map_done = .false. ! We will only build a map when we need it

          ! Loop over affected block columns of anode
          bcols: do
             if(cptr .gt. d_numrow) exit
             if(dnode%index(cptr) .gt. node%en) exit

             ! compute local index of block column in anode and find the id of 
             ! its diagonal block
             cb = (dnode%index(cptr) - node%sa)/node%nb + 1
             dblk = node%blk_sa
             do jb = 2, cb
                dblk = blocks(dblk)%last_blk + 1
             end do

             ! Find cptr2
             jlast = min(node%sa + cb*node%nb - 1, node%en)
             do cptr2 = cptr, d_numrow
                if(dnode%index(cptr2) > jlast) exit
             end do
             cptr2 = cptr2 - 1
             ! write(*,*)'cptr: ', cptr, ', cptr2: ', cptr2
             ! if(.not.map_done) call spllt_build_rowmap(node, tmpmap) 
             ! Loop over the blocks of snode
             ii = tmpmap(dnode%index(cptr)) 
             ! ii = -1
             ilast = cptr ! Set start of current block

             do i = cptr, d_numrow
                k = tmpmap(dnode%index(i))

                if(k.ne.ii) then

                   blk = dblk + ii - cb
                   bc => fdata%bc(blk)
                   
                   d_dblk = dnode%blk_sa
                   ! Loop over the block columns in node. 
                   do kk = 1, d_nc

                      d_bc_kk => fdata%bc(d_dblk)

                      call spllt_update_between_task( &
                           & fdata, &
                           & d_bc_kk, dnode, bc, fdata%nodes(snode), &
                           & cptr, cptr2, ilast, i-1, &
                           & fdata%row_list, fdata%col_list, fdata%workspace, &
                           & keep%lfact, keep%blocks, fdata%bc, &
                           & control, prio)

                         d_dblk = blocks(d_dblk)%last_blk + 1
                   end do
             
                   ii = k
                   ilast = i ! Update start of current block      
                end if
             end do

             blk = dblk + ii - cb
             bc => fdata%bc(blk)

             d_dblk = dnode%blk_sa
             ! Loop over the block columns in node. 
             do kk = 1, d_nc

                d_bc_kk => fdata%bc(d_dblk)

                call spllt_update_between_task( &
                     & fdata, &
                     & d_bc_kk, dnode, bc, fdata%nodes(snode), &
                     & cptr, cptr2, ilast, i-1, &
                     & fdata%row_list, fdata%col_list, fdata%workspace, &
                     & keep%lfact, keep%blocks, fdata%bc, &
                     & control, prio)

                d_dblk = blocks(d_dblk)%last_blk + 1
             end do

             ! Move cptr down, ready for next block column of anode
             cptr = cptr2 + 1
          end do bcols

       end do

       ! first block in node
       dblk = node%blk_sa

       ! Loop over the block columns in node. 
       do kk = 1, nc

          ! #if defined(SPLLT_USE_OMP)
          ! !$omp taskwait
          ! #endif

          ! A_kk          

          bc_kk => fdata%bc(dblk)
          call spllt_factorize_block_task(fdata, fdata%nodes(snode), bc_kk, keep%lfact, prio+3)

          ! #if defined(SPLLT_USE_OMP)
          ! !$omp taskwait
          ! #endif

          ! #if defined(SPLLT_USE_STARPU)
          ! call starpu_f_task_wait_for_all()
          ! #endif
          ! loop over the row blocks (that is, loop over blocks in block col)
          do ii = kk+1,nr
             ! do ii = nr,kk+1,-1             
             ! A_mk
             blk = dblk+ii-kk
             ! bc_ik => keep%blocks(blk)
             bc_ik => fdata%bc(blk)

             call spllt_solve_block_task(fdata, bc_kk, bc_ik, keep%lfact,prio+2)
          end do

          ! #if defined(SPLLT_USE_OMP)
          ! !$omp taskwait
          ! #endif

          ! #if defined(SPLLT_USE_STARPU)
          !          call starpu_f_task_wait_for_all()
          ! #endif

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
#endif

    call system_clock(stf_stop_t)
    write(*,'("[>] [spllt_stf_factorize] task insert time: ", es10.3, " s")') (stf_stop_t - stf_start_t)/real(stf_rate_t)

    return

9999 if(st.ne.0) then
     info%flag = spllt_error_allocation
     info%stat = st
     call spllt_print_err(info%flag, context='spllt_stf_factorize',st=st)
     return
  endif

  end subroutine spllt_stf_ll_factorize

  ! TODO comments!

  ! subroutine spllt_build_colmap(node, anode, cptr, colmap, ncolmap)
  !   use hsl_ma87_double
  !   implicit none

  !   type(node_type), intent(in) :: node  ! current node in the atree
  !   type(node_type), intent(in) :: anode  ! ancestor node in the atree
  !   integer, dimension(:), intent(out) :: colmap ! Workarray to hold map from row 
  !   integer, intent(inout) :: cptr  ! Position in snode of the first row
  !   ! matching a column of the current block column of anode.    
  !   integer, intent(inout) :: ncolmap ! number of entries in colmap

  !   integer :: sa, en
  !   integer :: numrow, numcol ! number of row/col in node
  !   integer :: cptr2  ! Position in snode of the last row 
  !   ! matching a column of the current block column of anode.
  !   integer :: jlast ! Last column in the cb-th block column of anode
  !   integer :: cb
  !   integer :: a_nb
    
  !   sa = node%sa
  !   en = node%en
  !   numcol = en - sa + 1
  !   numrow = size(node%index)
    
  !   a_nb = anode%nb

  !   ncolmap = 0

  !   ! Skip columns that come from other children
  !   do cptr = cptr, numrow
  !      if(node%index(cptr).ge.anode%sa) exit
  !   end do
  !   if(cptr.gt.numrow) return ! finished with node

  !   ! Loop over affected block columns of anode
  !   acols: do
  !      if(cptr.gt.numrow) exit
  !      if(node%index(cptr).gt.anode%en) exit

  !      cb = (node%index(cptr) - anode%sa)/anode%nb + 1
  !      ! Find cptr2
  !      jlast = min(anode%sa + cb*a_nb - 1, anode%en)
  !      do cptr2 = cptr, numrow
  !         if(node%index(cptr2) > jlast) exit
  !      end do
  !      cptr2 = cptr2 - 1 
       
  !      ncolmap = ncolmap + 1 
  !      colmap(ncolmap) = cptr2

  !      ! Move cptr down, ready for next block column of anode
  !      cptr = cptr2 + 1
  !   end do acols
    
  !   return
  ! end subroutine spllt_build_colmap

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
!   subroutine spllt_init_lfact(keep, pbl)
!     use hsl_ma87_double
! #if defined(SPLLT_USE_STARPU)
!     use  starpu_f_mod
! #endif
!     implicit none

!     type(MA87_keep), target, intent(inout) :: keep 
!     type(spllt_data_type), intent(inout) :: pbl

!     type(node_type), pointer :: node ! node in the atree    
!     integer(long) :: blk, dblk
!     integer :: nbcol, l_nb, sz, sa, en
!     integer :: blkm, blkn, size_bcol
!     integer :: snode, num_nodes
!     integer :: i
!     integer :: st ! stat parameter
!     integer :: ptr

!     num_nodes = keep%info%num_nodes

!     blk = 1
!     nbcol = 0
!     ! loop over the nodes
    
!     do snode = 1, num_nodes
!        ! Loop over the block columns in snode, allocating space 
!        ! l_nb is the size of the blocks and sz is number of
!        ! blocks in the current block column
       
!        node => keep%nodes(snode)

!        l_nb = node%nb
!        sz = (size(node%index) - 1) / l_nb + 1
!        sa = node%sa
!        en = node%en

!        size_bcol = 0
!        do i = sa, en, l_nb
!           nbcol = nbcol + 1
!           size_bcol = 0
!           dblk = blk
!           ! loop over the row blocks
!           do blk = dblk, dblk+sz-1
!              blkm = keep%blocks(blk)%blkm
!              blkn = keep%blocks(blk)%blkn
!              size_bcol = size_bcol + blkm*blkn
!           end do
!           allocate (keep%lfact(nbcol)%lcol(size_bcol),stat=st)
!           ! TODO trace error
!           ! TODO merge with previous loop?
!           ! register blocks hanldes in StarPU

!           ptr = 1
!           do blk = dblk, dblk+sz-1
!              blkm = keep%blocks(blk)%blkm
!              blkn = keep%blocks(blk)%blkn

!              pbl%bc(blk)%blk => keep%blocks(blk) 
! #if defined(SPLLT_USE_STARPU)
!              call starpu_matrix_data_register(pbl%bc(blk)%hdl, pbl%bc(blk)%mem_node, &
!                   & c_loc(keep%lfact(nbcol)%lcol(ptr)), blkm, blkm, blkn, &
!                   & int(wp,kind=c_size_t))
! #endif
!              pbl%bc(blk)%c => keep%lfact(nbcol)%lcol(ptr:ptr+blkm*blkn-1)
!              ptr = ptr + blkm*blkn
!           end do
!           sz = sz - 1
!        end do       
!     end do

!   end subroutine spllt_init_lfact

end module spllt_factorization_mod
