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
    use spllt_solve_dep_mod
    use utils_mod, ONLY : spllt_update_task_info
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
    integer :: bcol, dcol, col
    integer :: offset
    integer :: node
    integer :: dep
    integer :: i, j, r, chunkth
    integer :: threadID, nthread
    integer :: ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    type(spllt_block), pointer  :: p_blk_dep_update

    real(wp), pointer :: p_upd(:,:)
    real(wp), pointer :: p_xlocal(:,:)
    real(wp), pointer :: p_rhs(:)
    integer,  pointer :: p_dep(:)
    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted
    integer :: nftask

    type(spllt_block), pointer :: p_bc(:)
        
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
    p_bc      => fkeep%bc
    
    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs
    
    nftask      = 0

    if(use_omp_cases_method) then

      p_dep     => fkeep%bc(dblk)%fwd_dep
      ndep      = size(p_dep)
!     threadID  = omp_get_thread_num()
      chunk     = 0
      ndep_lvl  = 0

      if(ndep .eq. 0) then
        !$omp task                                &
        include 'spllt_solve_fwd_block_omp_decl.F90'

        include 'spllt_solve_fwd_block_worker.F90'

        !$omp end task
      else
        chunk = 10 ! Do not use chunk = 1 ; a non-sence
        lvl   = 1
        alpha = 1
        ndep_lvl = ndep ! #dep local to the lvl
        all_task_submitted = .false.

        do while(.not. all_task_submitted)
        
          nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

          beta = 1 - alpha

          do chunkth = 1, nchunk
            chunk_size = merge(ndep_lvl - (chunkth - 1) * chunk, chunk, &
              chunkth * chunk .gt. ndep_lvl)

            select case(chunk_size)

              case(0)
                print *, "No dep ??"

              !
              !This file contains the remaining cases that are generated through a script
              !
              include 'spllt_fwd_block_cases.F90'

            end select

            beta = beta + chunk_size
            if(ndep_lvl .le. chunk) then
              all_task_submitted = .true.
            else
              nftask = nftask + 1
            end if

          end do
          if(ndep_lvl .gt. chunk) then
            ndep_lvl = nchunk
            lvl = lvl + 1
            alpha = alpha * chunk
          end if
        end do
      end if

    else

      p_dep => fkeep%bc(dblk)%fwd_update_dep
      ndep  = size(p_dep)

      do dep = 1, ndep - 1

        !$omp task firstprivate(dep)                       &
        !$omp firstprivate(p_bc, p_dep)                    &
        !$omp private(threadID)                            &
        !$omp depend(in: p_bc(p_dep(dep)))                 &
        !$omp depend(inout: p_dep(1))

        threadID = omp_get_thread_num()
!       print '(a, i3, a, i3)', "SLV Fake Task dep ", dblk, " [in : ", &
!         p_dep(dep)

        !$omp end task

      end do
      nftask = merge(ndep - 1,0, ndep .gt. 1)

      !$omp task firstprivate(m, n, col, offset)                      &
      !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)              &
      !$omp firstprivate(dblk)                                        &
      !$omp firstprivate(nthread)                                     &
      !$omp firstprivate(p_blk_dep_update, p_lcol_update)             &
      !$omp firstprivate(p_dep, dep, p_bc)                            &
      !$omp firstprivate(p_upd, p_rhs, p_xlocal)                      &
      !$omp private(threadID, r, j, i)                                &
      !$omp depend(in: p_bc(p_dep(dep)))                              &
      !$omp depend(in: p_dep(1))                                      &
      !$omp depend(inout: p_bc(dblk))                                  

      call trace_event_start(trace_id, omp_get_thread_num())

      threadID  = omp_get_thread_num()
!     print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
!     print *, p_dep

      ! Sum contributions to rhs
      do r = 0, nrhs-1
        do j = 1, nthread
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

      !$omp end task
    end if

    call spllt_update_task_info(scheduler%task_info(scheduler%workerID), &
      ndep, 1, nftask)
  end subroutine spllt_solve_fwd_block_task


  subroutine spllt_solve_fwd_update_task(blk, node, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num
    use trace_mod
    use spllt_solve_dep_mod
    use utils_mod, ONLY : spllt_update_task_info
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
    integer, pointer            :: p_dep(:)
    integer                     :: blk_dep_solve
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer :: j
    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted
    integer :: nftask

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
    p_bc      => fkeep%bc

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

    nftask      = 0

    if(use_omp_cases_method) then

      p_dep     => fkeep%bc(blk)%fwd_dep
      ndep      = size(p_dep)
!     threadID  = omp_get_thread_num()
      chunk     = 0
      ndep_lvl  = 0

      if(ndep .eq. 0) then
        !$omp task                                &
        include 'spllt_solve_fwd_update_omp_decl.F90'

        include 'spllt_solve_fwd_update_worker.F90'

        !$omp end task
      else

        chunk = 10 ! Do not use chunk = 1 ; a non-sence
        lvl   = 1
        alpha = 1
        ndep_lvl = ndep ! #dep local to the lvl
        all_task_submitted = .false.

        do while(.not. all_task_submitted)
        
          nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

          beta = 1 - alpha

          do j = 1, nchunk
            chunk_size = merge(ndep_lvl - (j - 1) * chunk, chunk, &
              j * chunk .gt. ndep_lvl)
            select case(chunk_size)

              case(0)
                print *, "No dep ?? "

              !
              !This file contains the remaining cases that are generated through a script
              !
              include 'spllt_fwd_update_cases.F90'

            end select
            beta = beta + chunk_size
            if(ndep_lvl .le. chunk) then
              all_task_submitted = .true.
            else
              nftask = nftask + 1
            end if

          end do
          if(ndep_lvl .gt. chunk) then
            ndep_lvl = nchunk
            lvl = lvl + 1
            alpha = alpha * chunk
          end if
        end do
      end if
    else
    
    p_dep         => fkeep%bc(blk)%fwd_update_dep
    blk_dep_solve = fkeep%bc(blk)%fwd_solve_dep
    ndep          = size(p_dep)

    do dep = 1, ndep - 1

      !$omp task firstprivate(p_dep, p_bc, dep)                 &
      !$omp private(threadID)                                   &
      !$omp depend(in: p_bc(p_dep(dep)))                        &
      !$omp depend(inout: p_dep(1))

      threadID = omp_get_thread_num()
!     print '(a, i3, a, i3)', "UPD Fake Task dep of ", blk, &
!       " [in : ", p_dep(dep)

      !$omp end task

    end do
    nftask = merge(ndep - 1, 0, ndep .gt. 1)

    !$omp task firstprivate(m, n, col, offset)                    &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(blk)                                       &
    !$omp firstprivate(p_bc, p_dep, dep, blk_dep_solve)           &
    !$omp private(threadID)                                       &
    !$omp depend(in: p_bc(p_dep(dep)))                            &
    !$omp depend(in: p_bc(blk_dep_solve))                         &
    !$omp depend(in: p_dep(1))                                    &
    !$omp depend(inout: p_bc(blk))

    call trace_event_start(trace_id, omp_get_thread_num())
    
    threadID  = omp_get_thread_num()
!   print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!   print *, p_dep
!   print *, blk_dep_solve

    call slv_fwd_update(m, n, col, offset, p_index,         &
      p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
      p_upd(:, threadID + 1), ldr, p_rhs,                   &
      ldr, p_xlocal(:, threadID + 1))

    call trace_event_stop (trace_id, omp_get_thread_num())

    !$omp end task
  end if
  call spllt_update_task_info(scheduler%task_info(scheduler%workerID), &
    ndep, 1, nftask)
  end subroutine spllt_solve_fwd_update_task

  !*************************************************  
  !
  ! Backward solve with block on diagonal
  !         
  subroutine spllt_solve_bwd_block_task(dblk, nrhs, upd, rhs, ldr, xlocal, &
      fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    use spllt_solve_dep_mod
    use utils_mod, ONLY : spllt_update_task_info
    implicit none

    integer, intent(in)                       :: dblk ! Index of diagonal block
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:, :)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(spllt_omp_scheduler), intent(inout)  :: scheduler
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: bcol, col
    integer                     :: offset
    integer                     :: node
    integer                     :: dep
    integer                     :: i, j, r, chunkth
    integer                     :: threadID, nthread
    integer                     :: ndep
    integer                     :: blk_dep_update
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    real(wp), pointer           :: p_lcol_solve(:)
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve

    real(wp), pointer :: p_upd(:,:)
    real(wp), pointer :: p_xlocal(:,:)
    real(wp), pointer :: p_rhs(:)
    integer,  pointer :: p_dep_solve(:)
    integer,  pointer :: p_dep(:)

    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted
    integer :: nftask

    type(spllt_block), pointer :: p_bc(:)

    nthread = omp_get_num_threads()

    ! Get block info
    node      = fkeep%bc(dblk)%node
    n         = fkeep%bc(dblk)%blkn
    m         = fkeep%bc(dblk)%blkm
    sa        = fkeep%bc(dblk)%sa
    bcol      = fkeep%bc(dblk)%bcol ! Current block column
    col       = calc_col(fkeep%nodes(node), fkeep%bc(dblk)) ! current bcol
    col       = fkeep%nodes(node)%sa + (col-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
    p_bc      => fkeep%bc

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

    nftask      = 0

    if(use_omp_cases_method) then

      p_dep     => fkeep%bc(dblk)%bwd_dep
      ndep      = size(p_dep)
!     threadID  = omp_get_thread_num()
      chunk     = 0
      ndep_lvl  = 0

      if(ndep .eq. 0) then
        !$omp task                                &
        include 'spllt_solve_bwd_block_omp_decl.F90'

        include 'spllt_solve_bwd_block_worker.F90'

        !$omp end task
      else
        chunk = 10 ! Do not use chunk = 1 ; a non-sence
        lvl   = 1
        alpha = 1
        ndep_lvl = ndep ! #dep local to the lvl
        all_task_submitted = .false.

        do while(.not. all_task_submitted)
        
          nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

          beta = 1 - alpha

          do chunkth = 1, nchunk
            chunk_size = merge(ndep_lvl - (chunkth - 1) * chunk, chunk, &
              chunkth * chunk .gt. ndep_lvl)

            select case(chunk_size)

              case(0)
                print *, "No dep ??"

              !
              !This file contains the remaining cases that are generated through a script
              !
              include 'spllt_bwd_block_cases.F90'

            end select

            beta = beta + chunk_size
            if(ndep_lvl .le. chunk) then
              all_task_submitted = .true.
            else
              nftask = nftask + 1
            end if

          end do
          if(ndep_lvl .gt. chunk) then
            ndep_lvl = nchunk
            lvl = lvl + 1
            alpha = alpha * chunk
          end if
        end do
      end if
    else

      p_dep_solve => fkeep%bc(dblk)%bwd_solve_dep
      ndep        = size(p_dep_solve)

      do dep = 1, ndep - 1

        p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol
        p_blk_dep_solve => fkeep%bc(p_dep_solve(dep))

        !$omp task firstprivate(p_lcol_solve, p_blk_dep_solve)        &
        !$omp firstprivate(p_dep_solve)                               &
        !$omp firstprivate(dblk, dep)                                 &
        !$omp private(threadID)                                       &
        !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
        !$omp depend(inout: p_dep_solve(1))

        threadID = omp_get_thread_num()
  !     print '(a, i3, a, i3)', "SLV Fake Task dep ", dblk, " [in : ",&
  !       p_dep_solve(dep)

        !$omp end task

      end do
      nftask = merge(ndep - 1, 0, ndep .gt. 1)

      blk_dep_update  = fkeep%bc(dblk)%bwd_update_dep
      p_lcol_update   => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
      p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol

      p_blk_dep_update  => fkeep%bc(blk_dep_update)   
      p_blk_dep_solve   => fkeep%bc(p_dep_solve(dep))

      !$omp task firstprivate(m, n, col, offset, node, bcol)        &
      !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)            &
      !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
      !$omp firstprivate(p_dep_solve, dblk)                         &
      !$omp firstprivate(p_lcol_update, p_blk_dep_update)           &
      !$omp firstprivate(p_lcol_solve, p_blk_dep_solve)             &
      !$omp firstprivate(blk_dep_update)                            &
      !$omp firstprivate(nthread)                                   &
      !$omp private(threadID, r, j, i)                              &
      !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
      !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
      !$omp depend(in: p_dep_solve(1))                              &
      !$omp depend(inout: p_lcol(sa))                                  

      call trace_event_start(trace_id, omp_get_thread_num())

      threadID = omp_get_thread_num()
  !   print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
  !   print *, p_dep_solve
  !   print *, blk_dep_update

      ! Perform retangular update from diagonal block
      if(m .gt. n) then
         call slv_bwd_update(m - n, n, col, offset + n, p_index,      &
              p_lcol(sa + n * n : sa + n * m - 1), n, nrhs, p_rhs,    &
              p_upd(:, threadID + 1), ldr, p_xlocal(:, threadID + 1))
      endif

      ! Sum contributions to rhs
      do r = 0, nrhs-1
        do j = 1, nthread
          do i = col + r*ldr, col+n-1 + r*ldr
            p_rhs(i)    = p_rhs(i) + p_upd(i, j)
          end do
        end do
      end do

      ! Perform triangular solve
      call slv_solve(n, n, col, p_lcol(sa : sa + n * n - 1), &
           'Non-Transpose', 'Non-unit', nrhs, p_rhs, ldr)
      
      call trace_event_stop (trace_id, omp_get_thread_num())
      !$omp end task
    end if

    call spllt_update_task_info(scheduler%task_info(scheduler%workerID), &
      ndep, 1, nftask)
  end subroutine spllt_solve_bwd_block_task

  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         
  subroutine spllt_solve_bwd_udpate_task(blk, node, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    use spllt_solve_dep_mod
    use utils_mod, ONLY : spllt_update_task_info
    implicit none

    integer, intent(in)                       :: blk  ! Index of block 
    integer, intent(in)                       :: node 
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:,:)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:,:)
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
    integer, pointer            :: p_dep_solve(:)
    integer, pointer            :: p_dep(:)
    integer                     :: blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer :: j
    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted
    integer :: nftask         ! #fake tasks inserted into the runtime

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
    p_bc      => fkeep%bc

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

    nftask      = 0

    if(use_omp_cases_method) then

      p_dep     => fkeep%bc(blk)%bwd_dep
      ndep      = size(p_dep)
!     threadID  = omp_get_thread_num()
      chunk     = 0
      ndep_lvl  = 0

      if(ndep .eq. 0) then
        !$omp task                                &
        include 'spllt_solve_bwd_update_omp_decl.F90'

        include 'spllt_solve_bwd_update_worker.F90'

        !$omp end task
      else

        chunk = 10 ! Do not use chunk = 1 ; a non-sence
        lvl   = 1
        alpha = 1
        ndep_lvl = ndep ! #dep local to the lvl
        all_task_submitted = .false.

        do while(.not. all_task_submitted)
        
          nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

          beta = 1 - alpha

          do j = 1, nchunk
            chunk_size = merge(ndep_lvl - (j - 1) * chunk, chunk, &
              j * chunk .gt. ndep_lvl)
            select case(chunk_size)

              case(0)
                print *, "No dep ?? "

              !
              !This file contains the remaining cases that are generated through a script
              !
              include 'spllt_bwd_update_cases.F90'

            end select
            beta = beta + chunk_size
            if(ndep_lvl .le. chunk) then
              all_task_submitted = .true.
            else
              nftask = nftask + 1
            end if

          end do
          if(ndep_lvl .gt. chunk) then
            ndep_lvl = nchunk
            lvl = lvl + 1
            alpha = alpha * chunk
          end if
        end do
      end if
    else

      blk_dep_update  = fkeep%bc(blk)%bwd_update_dep
      p_dep_solve     => fkeep%bc(blk)%bwd_solve_dep
      ndep            = size(p_dep_solve)

      do dep = 1, ndep - 1

        p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol
        p_blk_dep_solve => fkeep%bc(p_dep_solve(dep))

        !$omp task firstprivate(p_lcol_solve, p_blk_dep_solve)        &
        !$omp firstprivate(p_dep_solve)                               &
        !$omp firstprivate(blk, dep)                                  &
        !$omp private(threadID)                                       &
        !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
        !$omp depend(inout: p_dep_solve(1))

        threadID = omp_get_thread_num()
!       print '(a, i3, a, i3)', "UPD Fake Task dep ", blk, " [in : ", &
!         p_dep_solve(dep)

        !$omp end task

      end do
      nftask = merge(ndep - 1, 0, ndep .gt. 1)

      p_lcol_update => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
      p_lcol_solve  => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol

      p_blk_dep_update  => fkeep%bc(blk_dep_update)
      p_blk_dep_solve   => fkeep%bc(p_dep_solve(dep))

      !$omp task firstprivate(m, n, col, offset, node, bcol)        &
      !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
      !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
      !$omp firstprivate(p_blk_dep_update, p_blk_dep_solve)         &
      !$omp firstprivate(p_lcol_update, p_lcol_solve)               &
      !$omp firstprivate(blk)                                       &
      !$omp private(threadID)                                       &
      !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
      !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
      !$omp depend(in: p_dep_solve(1))                              &
      !$omp depend(inout: p_lcol(blk_sa))                            

      call trace_event_start(trace_id, omp_get_thread_num())

      threadID  = omp_get_thread_num()
!     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!     print *, blk_dep_update
!     print *, p_dep_solve

      call  slv_bwd_update(m, n, col, offset, p_index,      &
            p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,   &
            p_rhs, p_upd(:, threadID + 1), ldr,             &
            xlocal(:, omp_get_thread_num() + 1))

      call trace_event_stop (trace_id, omp_get_thread_num())

      !$omp end task
    end if

    call spllt_update_task_info(scheduler%task_info(scheduler%workerID), &
      ndep, 1, nftask)
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
end module spllt_solve_task_mod
