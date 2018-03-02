module spllt_solve_task_mod
contains
  
  !*************************************************  
  !
  ! Forward solve with block on diagoanl
  !
  subroutine spllt_solve_fwd_block_task(dblk, nrhs, upd, rhs, ldr, xlocal, &
      fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    implicit none
    
    integer, intent(in)                       :: dblk !Index of diagonal block
    integer, intent(in)                       :: nrhs !Number of RHS
    integer, intent(in)                       :: ldr  !Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:, :)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target,  intent(inout)          :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(spllt_omp_scheduler), intent(inout)  :: scheduler
    
    ! Node info
    integer :: sa
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa
    integer :: bcol, dcol, col
    integer :: offset
    integer :: node
    integer :: dep
    integer :: i, j, r
    integer :: threadID, nthread
    integer :: ndep
    integer, dimension(:), pointer  :: p_index
    real(wp), dimension(:), pointer :: p_lcol
    real(wp), dimension(:), pointer :: p_lcol_update
!   integer, allocatable, target    :: dep_update(:)
    type(spllt_block), pointer      :: p_blk_dep_update
!   type(spllt_block), pointer      :: p_dblk   

    real(wp), pointer               :: p_upd(:,:)
    real(wp), pointer               :: p_xlocal(:,:)
    real(wp), pointer               :: p_rhs(:)
    integer, pointer                :: p_dep_update(:)
    integer          , pointer  :: p_child_blk_index(:)
        
    nthread   = omp_get_num_threads()

    ! Get block info
    node      = fkeep%bc(dblk)%node
    m         = fkeep%bc(dblk)%blkm
    n         = fkeep%bc(dblk)%blkn
    sa        = fkeep%bc(dblk)%sa
    bcol      = fkeep%bc(dblk)%bcol ! Current block column
    dcol      = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col       = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
!   p_dblk    => fkeep%bc(dblk)
    
    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs
    
!   print *, "Treat dep of ", dblk
!   blk_dep_update = fwd_update_dependency(fkeep, dblk)
    call fwd_update_dependency(fkeep, dblk, p_dep_update)

!   p_dep_update => dep_update
    ndep = size(p_dep_update)

!   print *, "Launch ", ndep, " tasks to solve the same block"

!   if(blk .eq. 39) then
!     print *, "========================== dblk : ", dblk 
!     call getPointerBlkIndex(fkeep, dblk, p_child_blk_index)
!     call print_iarray("dblk_index", size(p_child_blk_index), p_child_blk_index, 1)
!     call print_blk_index("dblk_index", size(p_child_blk_index), p_child_blk_index, 1)
!   end if

    do dep = 1, ndep - 1
!     print '(a, i3, a, i3, a, i3)', "Update Dep( ", dep, " ) of blk ", &
!       dblk, " is ", p_dep_update(dep)

      p_lcol_update     => fkeep%lfact(fkeep%bc(p_dep_update(dep))%bcol)%lcol
      p_blk_dep_update  => fkeep%bc(p_dep_update(dep))

      !$omp task firstprivate(p_lcol_update, p_blk_dep_update)      &
      !$omp firstprivate(p_dep_update, dep)                         &
      !$omp firstprivate(dblk)                                      &
      !$omp private(threadID)                                       &
      !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
      !$omp depend(inout: p_dep_update(  1))
     !print *, "Release of a fake task that treats ", p_dep_update(dep) 
      threadID = omp_get_thread_num()
!     print '(a, i3, a, i3)', "SLV Fake Task dep ", dblk, " [in : ", &
!       p_dep_update(dep)
      !$omp end task

      scheduler%nfake_task_insert = scheduler%nfake_task_insert + 1
    end do

    p_lcol_update     => fkeep%lfact(fkeep%bc(p_dep_update(dep))%bcol)%lcol
    p_blk_dep_update  => fkeep%bc(p_dep_update(dep))

    !$omp task firstprivate(m, n, col, offset)                      &
    !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)              &
!   !$omp firstprivate(dep_update, dblk)                            &
    !$omp firstprivate(dblk)                                        &
    !$omp firstprivate(nthread)                                     &
    !$omp firstprivate(p_blk_dep_update, p_lcol_update)             &
    !$omp firstprivate(p_dep_update)                                &
    !$omp firstprivate(p_upd, p_rhs, p_xlocal)                      &
    !$omp private(threadID, r, j, i)                                &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))            &
    !$omp depend(in: p_dep_update(1))                               &
    !$omp depend(inout: p_lcol(sa))                                  

!   print *, "[spllt_solve_fwd_block_task] solved by ", omp_get_thread_num()
    call trace_event_start(trace_id, omp_get_thread_num())

    threadID  = omp_get_thread_num()
!   print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
!   print *, p_dep_update
!   print "(a, i3, a, i4, a, i4, a, f20.2, a, f20.2)", "th :",            &
!     threadID, " dblk :", dblk,  " has a dependency with ",              &
!     dep_update, " using [in: ", p_lcol_update(p_blk_dep_update%sa), &
!     " and [inout: ", p_lcol(sa)

      ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
!       call print_darray("use upd ", n, &
!         p_upd(col + r * ldr : col + n - 1 + r * ldr, j))
        do i = col + r*ldr, col+n-1 + r*ldr
          p_rhs(i)    = p_rhs(i) + p_upd(i, j)
          p_upd(i,j)  = zero ! Reset in case of bwd solve
        end do
      end do
    end do


      ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa:sa+n*n-1),    &
         'Transpose    ', 'Non-unit', nrhs, p_rhs, ldr)
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m .gt. 0) then
       sa = sa + n * n
       call slv_fwd_update(m, n, col, offset, p_index,                &
            p_lcol(sa : sa + n * m - 1), n, nrhs,                     &
            p_upd(:, threadID + 1), ldr, p_rhs,                       &
            ldr, p_xlocal(:, threadID + 1))
    endif

    call trace_event_stop (trace_id, omp_get_thread_num())
    deallocate(p_dep_update)
    !$omp end task
!   !$omp taskwait
    scheduler%ntask_insert = scheduler%ntask_insert + 1
  end subroutine spllt_solve_fwd_block_task


  subroutine spllt_solve_fwd_update_task(blk, node, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num
    use trace_mod
    implicit none
    
    integer, intent(in)                       :: blk  ! Index of block
    integer, intent(in)                       :: node
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:,:)        
    real(wp), target, intent(in)              :: rhs(ldr*nrhs)
    real(wp), target, intent(out)             :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(spllt_omp_scheduler), intent(inout)  :: scheduler

    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID
    integer                     :: dep, ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    real(wp), pointer           :: p_lcol_solve(:)
    integer, pointer            :: p_dep_update(:)
    integer                     :: blk_dep_solve
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve
!   type(spllt_block), pointer  :: p_blk 
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    integer          , pointer  :: p_child_blk_index(:)
    integer                     :: tmp_blk

    ! Establish variables describing block
    n         = fkeep%bc(blk)%blkn
    m         = fkeep%bc(blk)%blkm
    blk_sa    = fkeep%bc(blk)%sa
    bcol      = fkeep%bc(blk)%bcol
    dcol      = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col       = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
    offset    = offset + (blk-fkeep%bc(blk)%dblk) &
      * fkeep%nodes(node)%nb ! this blk
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
!   p_blk     => fkeep%bc(blk)

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

  ! if(blk .eq. 39) then
  !   print *, "---------------------- INFO -------------------------"
  !   print *, "========================== blk : ", blk 
  !   call getPointerBlkIndex(fkeep, blk, p_child_blk_index)
! !   call print_iarray("blk_index", size(p_child_blk_index), p_child_blk_index, 1)
  !   call print_blk_index("blk_index", size(p_child_blk_index), p_child_blk_index, 1)

  !   do tmp_blk = 18, 20
  !     print *, "========================== blk : ", tmp_blk 
  !     call getPointerBlkIndex(fkeep, tmp_blk, p_child_blk_index)
  !     call print_blk_index("blk_index", size(p_child_blk_index), p_child_blk_index, 1)
  !   end do
  !   print *, "-----------------------------------------------"
  ! end if

!   print *, "Treat dep of ", blk
    call fwd_update_dependency(fkeep, blk, p_dep_update)
    blk_dep_solve  = fwd_solve_dependency(fkeep, blk)

    ndep = size(p_dep_update)
  ! if(ndep .gt. 1) then
! !   print *, "************************ CHECK update DEP"
  !   do dep = 1, ndep
  !     if(.not. contain(fkeep, blk, p_dep_update(dep))) then
  !   !   print '(a, i2, a, i2)', "SUCCESS !!! Child_blk ", p_dep_update(dep), &
  !   !     " intersect with ", blk
  !   ! else
  !       print '(a, i3, a, i3)', "ERROR   !!! Child_blk ", p_dep_update(dep), &
  !         " does not intersect with ", blk
  !       if(blk .eq. 39) then
  !         call getPointerBlkIndex(fkeep, p_dep_update(dep), p_child_blk_index)
! !         call print_iarray("child_index", size(p_child_blk_index), p_child_blk_index, 1)
  !         call print_blk_index("blk_index", size(p_child_blk_index), p_child_blk_index, 1)
  !       end if
  !     end if
  !   end do
! !   print *, "************************ END CHECK update DEP"
  ! end if

!   print *, "Launch ", ndep, " tasks to update the same block"

!   print '(a, i3, a, i3, a, i3)', "Solve  Dep( ", 1, " ) of blk ", blk, &
!     " is ", blk_dep_solve

!   if(ndep .gt. 1) then
!     print '(a, i5, a, i5)', "blk ", blk, " =>  Fake tasks : ", ndep - 1
!   end if
    do dep = 1, ndep - 1
!     print '(a, i3, a, i3, a, i3)', "Update Dep( ", dep, " ) of blk ", blk, &
!       " is ", p_dep_update(dep)

      p_lcol_update => fkeep%lfact(fkeep%bc(p_dep_update(dep))%bcol)%lcol
      p_blk_dep_update  => fkeep%bc(p_dep_update(dep))

      !$omp task firstprivate(p_lcol_update, p_blk_dep_update)      &
      !$omp firstprivate(p_dep_update, dep)                         &
      !$omp firstprivate(blk)                                       &
      !$omp private(threadID)                                       &
      !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
      !$omp depend(inout: p_dep_update(  1))
     !print *, "Release of a fake task that treats ", p_dep_update(dep) 
      threadID = omp_get_thread_num()
!     print '(a, i3, a, i3)', "UPD Fake Task dep of ", blk, &
!       " [in : ", p_dep_update(dep)
      !$omp end task

      scheduler%nfake_task_insert = scheduler%nfake_task_insert + 1
    end do

    p_lcol_update => fkeep%lfact(fkeep%bc(p_dep_update(dep))%bcol)%lcol
    p_lcol_solve  => fkeep%lfact(fkeep%bc(blk_dep_solve )%bcol)%lcol

    p_blk_dep_update  => fkeep%bc(p_dep_update(dep))
    p_blk_dep_solve   => fkeep%bc(blk_dep_solve)
    !$omp task firstprivate(m, n, col, offset)                    &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(p_blk_dep_update, p_blk_dep_solve)         &
    !$omp firstprivate(p_lcol_update, p_lcol_solve)               &
    !$omp firstprivate(blk)                                       &
!   !$omp firstprivate(p_blk)                                     &
!   !$omp firstprivate(blk, blk_dep_update, blk_dep_solve)        &
    !$omp private(threadID)                                       &
    !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(in: p_dep_update(1))                             &
    !$omp depend(inout: p_lcol(blk_sa))                            

!   print *, "[spllt_solve_fwd_update_task] solved by ", omp_get_thread_num()
    call trace_event_start(trace_id, omp_get_thread_num())
    
    threadID  = omp_get_thread_num()
!   print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!   print *, p_dep_update

!   print '(a, i3, a, i4, a, i4, a, i4, a, f20.2, a, f20.2, a, f20.2)',   &
!     "th :", threadID, "  blk :", blk, " has a dependency with ",        &
!     blk_dep_update, " and ", blk_dep_solve, " using [in: ",             &
!     p_lcol_solve(p_blk_dep_solve%sa), ", ",                             &
!     p_lcol_update(p_blk_dep_update%sa), " and [inout: ", p_lcol(blk_sa)

!   call print_darray("upd before update", ldr * nrhs, p_upd(:, threadID + 1))

    call slv_fwd_update(m, n, col, offset, p_index,         &
      p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
      p_upd(:, threadID + 1), ldr, p_rhs,                   &
      ldr, p_xlocal(:, threadID + 1))

!   call print_darray("upd after  update", ldr * nrhs, p_upd(:, threadID + 1))

    deallocate(p_dep_update)
    call trace_event_stop (trace_id, omp_get_thread_num())
    !$omp end task
!   !$omp taskwait
    scheduler%ntask_insert = scheduler%ntask_insert + 1
  end subroutine spllt_solve_fwd_update_task

  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         
  subroutine spllt_solve_bwd_block_task(dblk, nrhs, upd, rhs, ldr, xlocal, &
      fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    implicit none

    integer, intent(in)                     :: dblk ! Index of diagonal block
    integer, intent(in)                     :: nrhs ! Number of RHS
    integer, intent(in)                     :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)         :: upd(:, :)
    real(wp), target, intent(inout)         :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)         :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)   :: fkeep
    integer, intent(in)                     :: trace_id
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: node
    integer                     :: i, j, r
    integer                     :: threadID, nthread
    integer                     :: blk_dep_update
    real(wp), pointer           :: p_lcol_update(:)
    type(spllt_block), pointer  :: p_blk_dep_update
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_upd(:,:)
    real(wp), pointer           :: p_xlocal(:,:)
    real(wp), pointer           :: p_rhs(:)

    node = fkeep%bc(dblk)%node

    ! print *, "[spllt_solve_bwd_block_task] node = ", node

    nthread = omp_get_num_threads()

    ! Get block info
    n       = fkeep%bc(dblk)%blkn
    m       = fkeep%bc(dblk)%blkm
    sa      = fkeep%bc(dblk)%sa
    bcol    = fkeep%bc(dblk)%bcol ! Current block column
    col     = calc_col(fkeep%nodes(node), fkeep%bc(dblk)) ! current bcol
    col     = fkeep%nodes(node)%sa + (col-1)*fkeep%nodes(node)%nb
    offset  = col - fkeep%nodes(node)%sa + 1
    p_index => fkeep%nodes(node)%index
    p_lcol  => fkeep%lfact(bcol)%lcol

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

    blk_dep_update    = bwd_update_dependency(fkeep, dblk)
    p_lcol_update     => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
    p_blk_dep_update  => fkeep%bc(blk_dep_update)   

    !$omp task                                                    &
    !$omp firstprivate(m, n, col, offset, node, bcol)             &
    !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)            &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(p_lcol_update, p_blk_dep_update)           &
    !$omp firstprivate(blk_dep_update, nthread)                   &
    !$omp private(threadID, r, j, i)                              &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(inout: p_lcol(sa))                                  

!   print *, "[spllt_solve_bwd_block_task] treated by ", omp_get_thread_num()
!   print *, "work on dblk ", dblk, "which belongs to block column ", & 
!     fkeep%bc(dblk)%bcol
!   call trace_event_start(trace_id, omp_get_thread_num())

    threadID = omp_get_thread_num()
!   print "(a, i3, a, i4, a, i4)", "th :", threadID, " dblk :", dblk, &
!     " has a dependency with ", blk_dep_update

    ! Perform retangular update from diagonal block
    if(m .gt. n) then
       call slv_bwd_update(m - n, n, col, offset + n, p_index,      &
            p_lcol(sa + n * n : sa + n * m - 1), n, nrhs, p_rhs,    &
            p_upd(:, threadID + 1), ldr, p_xlocal(:, threadID + 1))
    endif

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
!       call print_darray("use upd ", n, &
!         p_upd(col + r * ldr : col + n - 1 + r * ldr, j))
        do i = col + r*ldr, col+n-1 + r*ldr
          p_rhs(i)    = p_rhs(i) + p_upd(i, j)
!         p_upd(i,j)  = zero ! Reset in case of bwd solve
        end do
      end do
    end do

    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa : sa + n * n - 1), &
         'Non-Transpose', 'Non-unit', nrhs, p_rhs, ldr)
    
!   call trace_event_stop (trace_id, omp_get_thread_num())
    !$omp end task

  end subroutine spllt_solve_bwd_block_task

  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         
  subroutine spllt_solve_bwd_udpate_task(blk, node, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    implicit none

    integer, intent(in)                     :: blk  ! Index of block 
    integer, intent(in)                     :: node 
    integer, intent(in)                     :: nrhs ! Number of RHS
    integer, intent(in)                     :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)         :: upd(:,:)
    real(wp), target, intent(inout)         :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)         :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)   :: fkeep
    integer, intent(in)                     :: trace_id
    
    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID, nthread
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    real(wp), pointer           :: p_lcol_solve(:)
    integer                     :: blk_dep_update
    integer                     :: blk_dep_solve
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve
    type(spllt_block), pointer  :: p_blk 
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)

    ! print *, "[spllt_solve_bwd_udpate_task] node = ", node

    nthread = omp_get_num_threads()

    !
    ! Backward update with block on diagoanl
    !

    ! Establish variables describing block
    n       = fkeep%bc(blk)%blkn
    m       = fkeep%bc(blk)%blkm
    blk_sa  = fkeep%bc(blk)%sa
    bcol    = fkeep%bc(blk)%bcol
    dcol    = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col     = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset  = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
    offset  = offset + (blk-fkeep%bc(blk)%dblk) * fkeep%nodes(node)%nb !this blk
    p_index => fkeep%nodes(node)%index
    p_lcol  => fkeep%lfact(bcol)%lcol

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

    blk_dep_update = bwd_update_dependency(fkeep, blk)
    blk_dep_solve  = bwd_solve_dependency(fkeep, blk)

    p_lcol_update => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
    p_lcol_solve  => fkeep%lfact(fkeep%bc(blk_dep_solve )%bcol)%lcol

    p_blk_dep_update  => fkeep%bc(blk_dep_update)
    p_blk_dep_solve   => fkeep%bc(blk_dep_solve)

    !$omp task                                                    &
    !$omp firstprivate(m, n, col, offset, node, bcol)             &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp firstprivate(p_blk_dep_update, p_blk_dep_solve)         &
    !$omp firstprivate(p_blk)                                     &
    !$omp firstprivate(blk, blk_dep_update, blk_dep_solve)        &
    !$omp firstprivate(p_lcol_update, p_lcol_solve)               &
    !$omp private(threadID)                                       &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(inout: p_lcol(blk_sa))                            

!   print *, "[spllt_solve_bwd_udpate_task] Thread id ", omp_get_thread_num()
!   call trace_event_start(trace_id, omp_get_thread_num())

    threadID  = omp_get_thread_num()
!   print '(a, i3, a, i4, a, i4, a, i4)', "th :", threadID, "  blk :", blk, &
!     " has a dependency with ", blk_dep_update, " and ", blk_dep_solve

!   call print_darray("upd before update", ldr * nrhs, p_upd(:, threadID + 1))

    call  slv_bwd_update(m, n, col, offset, p_index,      &
          p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,   &
          p_rhs, p_upd(:, threadID + 1), ldr,             &
          xlocal(:, omp_get_thread_num() + 1))

!   call print_darray("upd after  update", ldr * nrhs, p_upd(:, threadID + 1))

!   call trace_event_stop (trace_id, omp_get_thread_num())
    !$omp end task

  end subroutine spllt_solve_bwd_udpate_task

  !*************************************************
  !
  ! This function calculates column of a node we are on  
  integer function calc_col(node, bc)
    use spllt_data_mod    
    implicit none

    type(spllt_node), intent(in) :: node
    type(spllt_block), intent(in) :: bc

    calc_col = (size(node%index)-1)/node%nb + 1 ! no. row blks for node

    calc_col = calc_col - (bc%last_blk - bc%dblk + 1) + 1 ! column of node

  end function calc_col

  subroutine fwd_update_dependency(fkeep, blk, dep)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)     :: fkeep
    integer, intent(in)               :: blk    ! Index of block 
!   integer, allocatable, intent(out) :: dep(:) !List of dependencies
    integer, pointer,     intent(out) :: dep(:) !List of dependencies
    
    integer :: previous_dblk
    integer :: last_previous_dblk
    integer :: diff_bcol
    integer :: diff_previous_bcol
    integer :: node
    integer :: blk_sa, blk_en
    integer :: bcol_blk_sa, bcol
    integer :: blkm, nb, local_blk, offset
    integer :: child_node, i
    integer :: nblk_dep, nchild, nnode_child, nind
    integer :: new_method
    integer, allocatable :: nblk_child(:) ! Remove by using a workspace
    integer, allocatable :: node_child(:) ! Remove by using a workspace
    integer, allocatable :: node_index(:) ! Remove by using a workspace
    integer, allocatable :: node_child_bis(:) ! Remove by using a workspace

    integer :: target_blk1
    integer :: target_blk2
    integer :: target_node         
    integer :: target_blkm         
    integer :: local_target_blk   
    integer :: target_blk_sa       
    integer :: target_bcol_blk_sa  
    integer :: target_bcol         
    integer :: target_offset       

    node    = fkeep%bc(blk)%node
    blk_sa  = fkeep%nodes(node)%blk_sa
    blk_en  = fkeep%nodes(node)%blk_en
    bcol_blk_sa = fkeep%bc(blk_sa)%bcol
    bcol = fkeep%bc(blk)%bcol

    nb      = fkeep%nodes(node)%nb    
    blkm    = fkeep%bc(blk)%blkm
    local_blk = blk - fkeep%bc(blk)%dblk
    offset  = nb * ( local_blk + (bcol - bcol_blk_sa)) + 1
    
    new_method = 0
   !if(blk .eq. blk_sa) then
   !  print *, "********** Node ", node
   !  call print_node(fkeep, node)

   !end if

!   print '(a, i3, a, i3)', "Initial dep of blk ", blk,        &
!     " is ", fkeep%bc(blk)%dep_initial
!   print *, "Trying to print from ", offset, " to ", offset + blkm - 1
!   call print_iarray("block row", blkm,                &
!     fkeep%nodes(node)%index(offset: offset + blkm - 1))

!   if(blk .eq. blk_sa) then
!     print '(a, i3, a, i3, a, i3)', "In node ", node, " blk_sa = ",    &
!       blk_sa, " and blk_en ", blk_en
!     print *, "Node index ", fkeep%nodes(node)%index
!     print *, "Children nodes are ", fkeep%nodes(node)%child
!   end if
!       call get_update_dep(fkeep, child_node, node, pos)
!       print *, "This correponds in index of the child node as "
!       call print_iarray("rows in child_node", size(pos),    &
!         fkeep%nodes(node)%index(pos))
!       deallocate(pos)

!   print *, "Get dependencies of blk ", blk

    !If blk belongs to the first block column of the node,
    ! then we consider the node dependencies of the block
    if(fkeep%nodes(node)%blk_sa .eq. fkeep%bc(blk)%dblk) then

!     nchild  = fkeep%nodes(node)%nchild
      nchild  = node - fkeep%nodes(node)%least_desc
!     print *, "nnode to traverse ", nchild
      allocate(nblk_child(nchild + 1), node_child(nchild))
      nblk_child = zero

      nnode_child = 1
      cur_node    = node
!     print *, "Node ", node, " has ", nchild, " children"

      do i = 1, nchild
        node_child(nnode_child : nnode_child + fkeep%nodes(cur_node)%nchild - 1) = fkeep%nodes(cur_node)%child
        nnode_child = nnode_child + fkeep%nodes(cur_node)%nchild
        cur_node = node_child(i)
      end do
!     if(nchild .gt. 0) then
!       print *, "Children of node ", node, " are ", node_child
!     end if

      allocate(node_index(size(fkeep%nodes(node)%index)))
      node_index = fkeep%nodes(node)%index

!     allocate(dep(200))
!     do i = 1, nchild
!       child_node = node_child(i)
!       call reduce_ind_and_get_dep_blk(fkeep, node_index(offset : offset + blkm - 1), blkm, child_node, dep)
!     end do
      

      node_index  = fkeep%nodes(node)%index
!     target_blk1 = 17
!     target_blk2 = 14
!     if(blk .eq. target_blk1) then
!       print *, "index of blk ", target_blk1, " is ", node_index(offset: offset + blkm - 1)

!       target_node         = fkeep%bc(target_blk2)%node
!       target_blkm         = fkeep%bc(target_blk2)%blkm
!       local_target_blk    = target_blk2 - fkeep%bc(target_blk2)%dblk
!       target_blk_sa       = fkeep%nodes(target_node)%blk_sa
!       target_bcol_blk_sa  = fkeep%bc(target_blk_sa)%bcol
!       target_bcol         = fkeep%bc(target_blk2)%bcol
!       target_offset       = nb * ( local_target_blk + &
!         (target_bcol - target_bcol_blk_sa)) + 1

!       print *, "Index of blk ", target_blk2, " is ", &
!         fkeep%nodes(target_node)%index(target_offset : target_offset + target_blkm - 1)
!     end if
      if(new_method .eq. 1) then
        if(nchild .gt. 0) then

          allocate(node_child_bis(nchild + 1))
          node_child_bis    = zero
          node_child_bis(1) = one

          nind  = size(node_index(offset : offset + blkm - 1))
!         print *, "=============== GET number of update for block ", blk

          call getUpdateNDep(fkeep, node, node_index(offset: offset + blkm - 1),&
            nind, node_child_bis(2 : nchild+1))
!         call print_iarray("=====> #dep in NODE_CHILD ", nchild, &
!           node_child_bis(2 : nchild + 1), 1)

          !Restore the array
          node_index = fkeep%nodes(node)%index

          !PostTreatment of ndep to compute the accsum
          do i = 1, nchild
            node_child_bis(i + 1) = node_child_bis(i) + node_child_bis(i + 1)
          end do
!         call print_iarray("PostTreated node_child_bis ", nchild + 1, &
!           node_child_bis, 1)

          if(node_child_bis(nchild + 1) .gt. 1) then
            allocate(dep(node_child_bis(nchild + 1) - 1))

            call getUpdateDep(fkeep, node, node_index(offset : offset + blkm - 1),&
              nind, dep, node_child_bis)
!           print *, "=====> dep of ", blk, " in NODE_CHILD are ", dep

!           deallocate(dep)
          else
            allocate(dep(1))
            dep(1) = blk
          end if
        else
          allocate(dep(1))
          dep(1) = blk
        end if
      else
!     node_index = fkeep%nodes(node)%index

        do i = 1, nchild
  !       child_node = fkeep%nodes(fkeep%nodes(node)%least_desc + i - 1)%num
          child_node = node_child(i)

          call get_update_nblk(fkeep, child_node,                     &
  !         fkeep%nodes(node)%index(offset : offset + blkm - 1),       &
            node_index(offset : offset + blkm - 1),                    &
            nblk_child(i + 1))

!         print *, nblk_child(i + 1), " found in ", child_node

          if(nblk_child(i+1) .gt. 0) then
            nblk_child(i+1) = nblk_child(i+1) + nblk_child(i)
          else
            nblk_child(i+1) = nblk_child(i)
          end if
        end do

!       call print_iarray("nblk_child", nchild + 1, nblk_child)
        node_index = fkeep%nodes(node)%index

        if(nblk_child(nchild+1) .gt. 0) then
          allocate(dep(nblk_child(nchild + 1)))

          do i = 1, nchild
            if(nblk_child(i + 1) .ne. nblk_child(i)) then
  !           child_node = fkeep%nodes(fkeep%nodes(node)%least_desc + i - 1)%num
              child_node = node_child(i)
              call get_update_dep_blk(fkeep, child_node,            &
  !             fkeep%nodes(node)%index(offset : offset + blkm - 1), &
                node_index(offset : offset + blkm - 1), &
                dep(nblk_child(i) + 1 : nblk_child(i+1)))
            end if
          end do

!         print *, "For blk ", blk
!         call print_iarray("Node dependencies found ", size(dep), dep)
        else

          allocate(dep(1))
          dep(1) = blk

  !       call print_iarray("Simulated dependencies ",  &
  !         size(dep), dep)

        end if

        deallocate(nblk_child)
      end if
    else

      allocate(dep(1))

!     if(bcol .gt. bcol_blk_sa) then
      previous_dblk       = fkeep%bc(fkeep%bc(blk)%dblk - 1)%dblk
      last_previous_dblk  = fkeep%bc(previous_dblk)%last_blk
      diff_bcol           = blk - fkeep%bc(blk)%dblk
      diff_previous_bcol  = last_previous_dblk - previous_dblk

      if(last_previous_dblk .gt. (previous_dblk + diff_bcol)) then
        dep(1) = previous_dblk + diff_bcol + 1
      end if
!     end if
!     call print_iarray("Local dependencies ", size(dep), dep)
    end if

!   call print_iarray("Returned dep from fwd_update_dependency",  &
!     size(dep), dep)

  end subroutine fwd_update_dependency

  integer function bwd_update_dependency(fkeep, blk)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)   :: fkeep
    integer, intent(in)             :: blk  ! Index of block 
    
    bwd_update_dependency = blk

    if(fkeep%bc(blk)%last_blk .ne. blk) then
      bwd_update_dependency = blk + 1
    end if

  end function bwd_update_dependency

  integer function fwd_solve_dependency(fkeep, blk)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)   :: fkeep
    integer, intent(in)             :: blk  ! Index of block 

    if(fkeep%bc(blk)%dblk .ne. blk) then
      fwd_solve_dependency = fkeep%bc(blk)%dblk
    else
      fwd_solve_dependency = blk
    end if

  end function fwd_solve_dependency

  integer function bwd_solve_dependency(fkeep, blk)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)   :: fkeep
    integer, intent(in)             :: blk  ! Index of block 

    integer :: next_dblk    
    integer :: dist

    dist = blk - fkeep%bc(blk)%dblk
    next_dblk = fkeep%bc(blk)%last_blk + 1

    if(dist .gt. 0) then
      bwd_solve_dependency = next_dblk + dist - 1
    else
      bwd_solve_dependency = blk
    end if

  end function bwd_solve_dependency

end module spllt_solve_task_mod
