module spllt_solve_kernels_mod
  
contains


  !*************************************************
  !
  ! TASK_BSOLV, TASK_FSOLV
  ! B_j <- L_jj^-1 B_j
  ! B_j <- L_jj^-T B_j
  !
  ! Note: While diagonal blocks may be trapezoidal, this is handled at the
  ! level calling this subroutine through a call to slv_fwd_update or
  ! slv_bwd_update

  subroutine slv_solve(n, nelim, col, dest, trans, unit,  &
       nrhs, rhs, ldr)
    use spllt_data_mod
    implicit none

    integer,            intent(in)    :: n      ! leading dimension of 
                                                ! diag. block
    integer,            intent(in)    :: nelim  ! number eliminations 
                                                ! (immediate return if =0)
    integer,            intent(in)    :: col    ! start of block column 
                                                ! variables in rhs
    real(wp),           intent(in)    :: dest(*)! holds destination block
    character(len=13),  intent(in)    :: trans  ! set to :
                                                ! 'Transpose    ' 
                                                ! for forward substitution 
                                                ! 'Non-Transpose' 
                                                ! for back substitution
    character(len=8),   intent(in)    :: unit   ! set to 
                                                ! 'Non-Unit' 
                                                ! for positive-definite case
    integer,            intent(in)    :: nrhs   ! number of right-hand sides
    integer,            intent(in)    :: ldr    ! leading extent of rhs
    real(wp),           intent(inout) :: rhs(ldr*nrhs)

    !%%%   integer :: this_thread
    !%%%   integer :: t_start, t_end
   !integer :: r

    if(nelim.eq.0) return

    !%%%   if(control%unit_log.gt.0) call system_clock(t_start)

    if(nrhs.eq.1) then
      call dtrsv('Upper', trans, unit, nelim, dest, n, rhs(col), 1)
    else
     !do r = 0, nrhs - 1
     !  print *, "SOLVE org RHS", rhs(col + r * ldr : col + n - 1 + r * ldr)
     !end do
      call dtrsm('Left', 'Upper', trans, unit, nelim, nrhs, one, &
           dest, n, rhs(col), ldr)
     !do r = 0, nrhs - 1
     !  print *, "SOLVE upd RHS", rhs(col + r * ldr : col + n - 1 + r * ldr)
     !end do
    endif

  end subroutine slv_solve
  !*************************************************
  !
  ! TASK_FUPD
  ! B_j <- B_j - L_ij B_i
  !
  subroutine slv_fwd_update(m, nelim, col, offset, index, dest, ldd, nrhs, &
       upd, ldu, rhs, ldr, xlocal)
    use spllt_data_mod
    use timer_mod
 !$ use omp_lib, ONLY : omp_get_thread_num
    implicit none

    integer,  intent(in)    :: m              ! number of rows in block
    integer,  intent(in)    :: nelim          ! # eliminations 
                                              ! (immediate return if =0)
    integer,  intent(in)    :: col            ! start of block column variables
                                              ! in rhs
    integer,  intent(in)    :: offset         ! offset into index we start at
    integer,  intent(in)    :: index(*)
    integer,  intent(in)    :: ldd            ! leading dimension of block
    real(wp), intent(in)    :: dest(m*ldd)    ! holds destination block
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: ldu            ! leading extent of upd
    real(wp), intent(inout) :: upd(ldu*nrhs)  ! vector to update
    integer,  intent(in)    :: ldr            ! leading extent of rhs
    real(wp), intent(in)    :: rhs(ldr*nrhs)  ! rhs vector
    real(wp), intent(out)   :: xlocal(*)

    integer   :: i
    integer   :: j
    integer   :: k
    integer   :: r ! right hand side loop variable
    real(wp)  :: w ! temporary work value
    !%%%  integer :: t_start, t_end, this_thread
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    type(spllt_timer_t), save :: timer
    integer                   :: threadID

    threadID = 0
 !$ threadID = omp_get_thread_num()

    call spllt_open_timer(threadID, "slv_fwd_update", timer)
#endif

    if(nelim.eq.0) return

    ! forward substitution
    if(nrhs.eq.1) then
       if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single rhs, BLAS 2

          call dgemv('T', nelim, m, -one, dest, ldd, rhs(col), 1, zero, &
               xlocal, 1)

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
          call spllt_tic("fwd::scatter", 1, threadID, timer)
#endif
          ! Copy xlocal out
          j = 1
          do i = offset, offset+m-1
             upd(index(i)) = upd(index(i)) + xlocal(j)
             j = j + 1
          end do
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
          call spllt_ftac(1, threadID, timer)
#endif
       else
!!! Single rhs, direct update
          j = 1
          do i = offset, offset+m-1
             w = zero
             do k = col, col+nelim-1
                w = w - dest(j)*rhs(k)
                j = j + 1
             end do
             j = j + (ldd-nelim)
             upd(index(i)) = upd(index(i)) + w
          end do
       endif
    else
!!! Multiple rhs, BLAS 3
      call dgemm('T', 'N', m, nrhs, nelim, -one, dest, ldd, rhs(col), ldr, &
           zero, xlocal, m)
     !do r = 0, nrhs - 1
     !  print *, "RHS", rhs(col + r * ldr: col + r * ldr + nelim - 1)
     !  print *, "xlocal", xlocal(1 + r * m : r * m + m)
     !end do

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
      call spllt_tic("fwd::scatter", 1, threadID, timer)
#endif
      ! Copy xlocal out
      j = 1
      do i = offset, offset+m-1
         do r = 0, nrhs-1
            upd(index(i)+r*ldu) = upd(index(i)+r*ldu) + xlocal(j+r*m)
         end do
         j = j + 1
      end do
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
      call spllt_ftac(1, threadID, timer)
#endif
    endif

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    call spllt_close_ftimer(threadID, timer)
#endif

  end subroutine slv_fwd_update



  !*************************************************
  !
  ! TASK_BUPD
  ! B_i <- B_i - L_ij^-T B_j
  !
  subroutine slv_bwd_update(m, nelim, col, offset, index, dest, ldd, nrhs, &
      rhs, upd, ldr, xlocal)
    use spllt_data_mod
    use timer_mod
 !$ use omp_lib, ONLY : omp_get_thread_num
    implicit none

    integer, intent(in) :: m ! number of rows in block
    integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
    integer, intent(in) :: col ! start of block column variables in rhs
    integer, intent(in) :: offset ! offset into index we start at
    integer, dimension(*), intent(in) :: index
    integer, intent(in) :: ldd ! leading dimension of block
    real(wp), dimension(m*ldd), intent(in) :: dest ! holds block
    integer, intent(in) :: nrhs
    integer, intent(in) :: ldr  ! leading extent of rhs
    real(wp), intent(inout) :: rhs(ldr*nrhs)
    real(wp), intent(inout) :: upd(ldr*nrhs)
    real(wp), dimension(:), intent(out) :: xlocal

    integer :: i
    integer :: j
    integer :: k
    integer :: r ! right hand side loop variable
    real(wp) :: w ! temporary work variable
    !%%%  integer :: t_start, t_end, this_thread
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer
    integer                   :: threadID

    threadID = 0
 !$ threadID = omp_get_thread_num()

    call spllt_open_timer(threadID, "slv_bwd_update", timer)
#endif

    if(nelim.eq.0) return

    !%%%  if(control%unit_log.gt.0) call system_clock(t_start)

    ! backward substitution
    if(nrhs.eq.1) then
      if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single right-hand side, BLAS 2

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
        call spllt_tic("bwd::gather", 1, threadID, timer)
#endif
         ! Copy xlocal in
        j = 1
        do i = offset, offset+m-1
           xlocal(j) = rhs(index(i))
           j = j + 1
        end do
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
        call spllt_ftac(1, threadID, timer)
#endif

        call dgemv('N', nelim, m, -one, dest, ldd, xlocal, 1, one, &
             upd(col), 1)
      else
!!! Single right-hand side, direct update
        j = 1
        do i = offset, offset+m-1
          w = rhs(index(i))
          do k = col, col + nelim - 1
            upd(k) = upd(k) - dest(j)*w
            j = j + 1
          end do
          j = j + (ldd-nelim)
        end do
      endif
    else
!!! Multiple RHS, BLAS 3

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
      call spllt_tic("bwd::gather", 1, threadID, timer)
#endif
      ! Copy xlocal in
      j = 1
!     do i = offset, offset+m-1
!        do r = 0, nrhs-1
!           xlocal(j+r*m) = rhs(index(i)+r*ldr)
!        end do
!        j = j + 1
!     end do
      !Test change order of loop
      do r = 0, nrhs-1
         j = 1
         do i = offset, offset+m-1
           xlocal(j+r*m) = rhs(index(i)+r*ldr)
           j = j + 1
         end do
      end do
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
      call spllt_ftac(1, threadID, timer)
#endif

      call dgemm('N', 'N', nelim, nrhs, m, -one, dest, ldd, xlocal, m, &
           one, upd(col), ldr)
    endif

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_ftimer(threadID, timer)
#endif

  end subroutine slv_bwd_update


  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Kernel of solve steps
  !
  subroutine solve_bwd_block_work(m, n, col, offset, index, lcol, sa, nrhs, &
      upd, rhs, threadID, nthread, ldr, xlocal, flops)
    use spllt_data_mod
    use timer_mod
    implicit none
    integer,                  intent(in)    :: m
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: col
    integer,                  intent(in)    :: offset
    integer,                  intent(in)    :: sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: nthread
    integer,                  intent(in)    :: ldr
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:,:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)
    double precision,         intent(out)   :: flops

    integer :: i, r, j
    double precision :: lflops
#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer

    call spllt_open_timer(threadID, "solve_bwd_block_work", timer)
#endif

    flops = zero

    ! Perform retangular update from diagonal block
    if(m .gt. n) then
      call slv_bwd_update(m - n, n, col, offset + n, index,     &
        lcol(sa + n * n : sa + n * m - 1), n, nrhs, rhs,      &
        upd(:, threadID + 1), ldr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (n * nrhs * (m-n))
#endif
    endif

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd block task reduction", 1, threadID, timer)
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_tic("bwd::reduction", 2, threadID, timer)
#endif
    ! Sum contributions to rhs
    if(nthread .eq. 1) then
!     print *, "Nthread given to solve_bwd_block_work : ", nthread
!     print *, "Lighter reduction"
      do r = 0, nrhs-1
        j = threadID + 1
        do i = col + r*ldr, col+n-1 + r*ldr
          rhs(i)    = rhs(i) + upd(i, j)
          upd(i,j)  = zero ! Reset in case of next fwd solve
        end do
      end do
#if defined(SPLLT_PROFILING_FLOP)
    lflops = n * nrhs
#endif
    else
!     print *, "Nthread given to solve_bwd_block_work : ", nthread
      do r = 0, nrhs-1
        do j = 1, nthread
          do i = col + r*ldr, col+n-1 + r*ldr
            rhs(i)    = rhs(i) + upd(i, j)
            upd(i,j)  = zero ! Reset in case of next fwd solve
          end do
        end do
      end do
#if defined(SPLLT_PROFILING_FLOP)
    lflops = n * nrhs * nthread
#endif
    end if
#if defined(SPLLT_PROFILING_FLOP)
    flops = flops + lflops
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_ftac(2, threadID, timer)
#endif
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(1, threadID, timer, lflops)
#endif

    ! Perform triangular solve
    call slv_solve(n, n, col, lcol(sa : sa + n * n - 1), &
         'Non-Transpose', 'Non-unit', nrhs, rhs, ldr)

#if defined(SPLLT_PROFILING_FLOP)
    flops = flops + n * n * nrhs
#endif

#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_timer(threadID, timer, flops)
#endif
  end subroutine solve_bwd_block_work



  subroutine solve_bwd_update_work(m, n, col, offset, index, lcol, blk_sa, &
      nrhs, upd, rhs, threadID, ldr, xlocal, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(in)    :: m
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: col
    integer,                  intent(in)    :: offset
    integer,                  intent(in)    :: blk_sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: ldr
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:,:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)
    double precision,         intent(out)   :: flops

    call  slv_bwd_update(m, n, col, offset, index,      &
          lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,   &
          rhs, upd(:, threadID + 1), ldr,               &
          xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (n * nrhs * m)
#endif

  end subroutine solve_bwd_update_work



  subroutine solve_fwd_update_work(m, n, col, offset, index, lcol, blk_sa, &
      nrhs, upd, rhs, threadID, ldr, xlocal, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(in)    :: m
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: col
    integer,                  intent(in)    :: offset
    integer,                  intent(in)    :: blk_sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: ldr
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:,:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)
    double precision,         intent(out)   :: flops

    call slv_fwd_update(m, n, col, offset, index,         &
      lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
      upd(:, threadID + 1), ldr, rhs,                   &
      ldr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (m * nrhs * n)
#endif

  end subroutine solve_fwd_update_work



  subroutine solve_fwd_block_work(m, n, col, offset, index, lcol, sa, nrhs, &
      upd, rhs, threadID, nthread, ldr, xlocal, flops)
    use spllt_data_mod
    use timer_mod
    implicit none
    integer,                  intent(inout) :: m
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: col
    integer,                  intent(inout) :: offset
    integer,                  intent(inout) :: sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: nthread
    integer,                  intent(in)    :: ldr
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:,:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)
    double precision,         intent(out)   :: flops

    integer :: i, r, j
    double precision :: lflops
#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer

    call spllt_open_timer(threadID, "solve_fwd_block_work", timer)
#endif

    flops = zero

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd block task reduction", 1, threadID, timer)
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_tic("fwd::reduction", 2, threadID, timer)
#endif
    ! Sum contributions to rhs
    if(nthread .eq. 1) then
      do r = 0, nrhs - 1
        j = threadID + 1
        do i = col + r*ldr, col+n-1 + r*ldr
          rhs(i)    = rhs(i) + upd(i, j)
          upd(i,j)  = zero ! Reset in case of bwd solve
        end do
      end do
#if defined(SPLLT_PROFILING_FLOP)
    lflops = n * nrhs
#endif
    else
      do r = 0, nrhs - 1
        do j = 1, nthread
          do i = col + r*ldr, col+n-1 + r*ldr
            rhs(i)    = rhs(i) + upd(i, j)
            upd(i,j)  = zero ! Reset in case of bwd solve
          end do
        end do
      end do
#if defined(SPLLT_PROFILING_FLOP)
    lflops = n * nrhs * nthread
#endif
    end if
#if defined(SPLLT_PROFILING_FLOP)
    flops = flops + lflops
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_ftac(2, threadID, timer)
#endif
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(1, threadID, timer, lflops)
#endif


    ! Perform triangular solve
    call slv_solve(n, n, col, lcol(sa:sa+n*n-1),    &
      'Transpose    ', 'Non-unit', nrhs, rhs, ldr)
#if defined(SPLLT_PROFILING_FLOP)
    flops = n * n * nrhs
#endif
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m .gt. 0) then
      sa = sa + n * n
      call slv_fwd_update(m, n, col, offset, index,             &
        lcol(sa : sa + n * m - 1), n, nrhs,                     &
        upd(:, threadID + 1), ldr, rhs,                         &
        ldr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (n * nrhs * m)
#endif
    endif

#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_timer(threadID, timer, flops)
#endif
  end subroutine solve_fwd_block_work



  subroutine solve_bwd_node(nrhs, rhs, ldr, fkeep, node, xlocal, rhs_local, &
      task_manager, trace_id)
    use spllt_data_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    integer,                    intent(in)    :: node 
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    real(wp),                   intent(inout) :: xlocal(:,:)
    real(wp),                   intent(inout) :: rhs_local(:,:)
    class(task_manager_base ),  intent(inout) :: task_manager
    integer, optional,          intent(in)    :: trace_id

    integer                 :: sa, en
    integer                 :: numcol, numrow ! #column/row in node 
    integer                 :: nc, nr         ! #block-column/block-row in node
    integer                 :: s_nb           ! Block size in node
    integer                 :: jj, ii
    integer                 :: dblk, blk
    integer                 :: bwd_block_id
    integer                 :: bwd_update_id
    type(spllt_timer_t), save :: timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_bwd_node", timer)
#endif

    if(present(trace_id)) then
      bwd_block_id  = trace_id
      bwd_update_id = trace_id
    else
      bwd_block_id  = task_manager%trace_ids(trace_bwd_block_pos)
      bwd_update_id = task_manager%trace_ids(trace_bwd_update_pos)
    end if

    ! Get node info
    s_nb   = fkeep%nodes(node)%nb
    sa     = fkeep%nodes(node)%sa
    en     = fkeep%nodes(node)%en
    numcol = en - sa + 1
    numrow = size(fkeep%nodes(node)%index)
    nc     = (numcol-1) / s_nb + 1
    nr     = (numrow-1) / s_nb + 1 

    ! Get first diag block in node
    dblk = fkeep%bc(fkeep%nodes(node)%blk_en)%dblk

    ! Loop over block columns
    do jj = nc, 1, -1

      do ii = nr, jj+1, -1
        
        blk = dblk+ii-jj ! Block index

        !
        ! Backward update with block on diagoanl
        !
        call spllt_tic("submit bwd update", 1, task_manager%workerID, timer)
        call task_manager%solve_bwd_update_task(blk, node, nrhs, rhs_local, &
          rhs, ldr, xlocal, fkeep, bwd_update_id)
        call spllt_tac(1, task_manager%workerID, timer)

      end do

      !
      ! Backward solve with block on diagoanl
      !
      call spllt_tic("submit bwd block", 2, task_manager%workerID, timer)
      call task_manager%solve_bwd_block_task(dblk, nrhs, rhs_local, rhs, ldr,&
        xlocal, fkeep, bwd_block_id)
      call spllt_tac(2, task_manager%workerID, timer)
     
      ! Update diag block in node       
      if (jj .gt. 1) dblk = fkeep%bc(dblk-1)%dblk
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_bwd_node



  subroutine solve_fwd_node(nrhs, rhs, ldr, fkeep, node, xlocal, rhs_local, &
      task_manager, trace_id)
    use spllt_data_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: ldr  ! Leading dimension of RHS
    integer,                    intent(in)    :: node 
    real(wp),                   intent(inout) :: rhs(ldr*nrhs)
    real(wp),                   intent(inout) :: xlocal(:,:)
    real(wp),                   intent(inout) :: rhs_local(:,:)
    class(task_manager_base ),  intent(inout) :: task_manager
    integer, optional,          intent(in)    :: trace_id

    integer                 :: sa, en
    integer                 :: numcol, numrow ! #column/row in node 
    integer                 :: nc, nr         ! #block-column/block-row in node
    integer                 :: jj, ii
    integer                 :: dblk           ! Diagonal index 
    integer                 :: s_nb           ! Block size in node
    integer                 :: blk            ! Block index
    integer                 :: fwd_block_id
    integer                 :: fwd_update_id
    type(spllt_timer_t), save :: timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_fwd_node", timer)
#endif

    if(present(trace_id)) then
      fwd_block_id  = trace_id
      fwd_update_id = trace_id
    else
      fwd_block_id  = task_manager%trace_ids(trace_fwd_block_pos)
      fwd_update_id = task_manager%trace_ids(trace_fwd_update_pos)
    end if

    ! Get node info
    s_nb   = fkeep%nodes(node)%nb
    sa     = fkeep%nodes(node)%sa
    en     = fkeep%nodes(node)%en
    numcol = en - sa + 1
    numrow = size(fkeep%nodes(node)%index)
    nc     = (numcol-1) / s_nb + 1
    nr     = (numrow-1) / s_nb + 1 
    
    ! Get first diag block in node
    dblk = fkeep%nodes(node)%blk_sa

    ! Loop over block columns
    do jj = 1, nc
       
      !
      ! Forward solve with block on diagoanl
      !
      call spllt_tic("submit fwd block", 1, task_manager%workerID, timer)
      call task_manager%solve_fwd_block_task(dblk, nrhs, rhs_local, rhs, &
        ldr, xlocal, fkeep, fwd_block_id)
      call spllt_tac(1, task_manager%workerID, timer)

      do ii = jj+1, nr

        blk = dblk+ii-jj

        !
        ! Forward update with off-diagonal
        !
        call spllt_tic("submit fwd update", 2, task_manager%workerID, timer)
        call task_manager%solve_fwd_update_task(blk, node, nrhs, rhs_local,&
          rhs, ldr, xlocal, fkeep, fwd_update_id)
        call spllt_tac(2, task_manager%workerID, timer)

      end do
      
      ! Update diag block in node          
      dblk = fkeep%bc(dblk)%last_blk + 1
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_fwd_node


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! InterLeave subroutines
!
!
!


  subroutine slv_solve_ileave(n, nelim, col, dest, trans, unit,  &
       nrhs, rhs)
    use spllt_data_mod
    implicit none

    integer,            intent(in)    :: n      ! leading dimension of 
                                                ! diag. block
    integer,            intent(in)    :: nelim  ! number eliminations 
                                                ! (immediate return if =0)
    integer,            intent(in)    :: col    ! start of block column 
                                                ! variables in rhs
    real(wp),           intent(in)    :: dest(*)! holds destination block
    character(len=13),  intent(in)    :: trans  ! set to :
                                                ! 'Transpose    ' 
                                                ! for forward substitution 
                                                ! 'Non-Transpose' 
                                                ! for back substitution
    character(len=8),   intent(in)    :: unit   ! set to 
                                                ! 'Non-Unit' 
                                                ! for positive-definite case
    integer,            intent(in)    :: nrhs   ! number of right-hand sides
    real(wp),           intent(inout) :: rhs(:)

    !%%%   integer :: this_thread
    !%%%   integer :: t_start, t_end
   !integer :: r

    if(nelim.eq.0) return

    if(nrhs.eq.1) then
     !print *, "SOLVE org RHS", rhs(col : col + n - 1)
      call dtrsv('Upper', trans, unit, nelim, dest, n, rhs(col), 1)
     !print *, "SOLVE upd RHS", rhs(col : col + n - 1)
    else
     !do r = 0, nrhs - 1
     !  print *, "SOLVE org RHS", rhs(col + r * nelim : col + n - 1 + r * nelim)
     !end do
      call dtrsm('Left', 'Upper', trans, unit, nelim, nrhs, one, &
           dest, n, rhs(col), nelim)
     !do r = 0, nrhs - 1
     !  print *, "SOLVE upd RHS", rhs(col + r * nelim : col + n - 1 + r * nelim)
     !end do
    endif

  end subroutine slv_solve_ileave



  subroutine slv_fwd_update_ileave(m, nelim, col, offset, index, dest, &
      ldd, nrhs, upd, tdu, rhs, n, rhsPtr, indir_rhs, xlocal)
    use spllt_data_mod
    use timer_mod
 !$ use omp_lib, ONLY : omp_get_thread_num
    implicit none

    integer,  intent(in)    :: m              ! number of rows in block
    integer,  intent(in)    :: nelim          ! # eliminations 
                                              ! (immediate return if =0)
    integer,  intent(in)    :: col            ! start of block column variables
                                              ! in rhs
    integer,  intent(in)    :: offset         ! offset into index we start at
    integer,  intent(in)    :: index(*)
    integer,  intent(in)    :: ldd            ! leading dimension of block
    real(wp), intent(in)    :: dest(m*ldd)    ! holds destination block
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: tdu            ! leading dim of upd
    integer,  intent(in)    :: n              ! Order of the system
    integer,  intent(in)    :: rhsPtr(:)
    integer,  intent(in)    :: indir_rhs(:)
    real(wp), intent(inout) :: upd(:)         ! vector to update
    real(wp), intent(in)    :: rhs(n*nrhs)    ! rhs vector
    real(wp), intent(out)   :: xlocal(*)

    integer   :: i, ind_i, ind_n
    integer   :: j
    integer   :: k
    integer   :: r ! right hand side loop variable
    integer   :: l
    integer   :: lrow, rowPtr, lcol, ndblk, pos
    integer   :: rhs_blk_size, rhs_thread_size
    real(wp)  :: w ! temporary work value
    integer                   :: threadID
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    type(spllt_timer_t), save :: timer
#endif

    threadID = 0
 !$ threadID = omp_get_thread_num()

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    call spllt_open_timer(threadID, "slv_fwd_update_ileave", timer)
#endif

    if(nelim.eq.0) return

    l = indir_rhs(col) 
    lcol = rhsPtr(l) * nrhs + 1

    ! forward substitution
    if(nrhs.eq.1) then
       if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single rhs, BLAS 2

          call dgemv('T', nelim, m, -one, dest, ldd, rhs(col), 1, zero, &
               xlocal, 1)

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
          call spllt_tic("fwd::scatter", 1, threadID, timer)
#endif
          ndblk   = size(rhsPtr) - 1

          ! Copy xlocal out
          l = 1
          i = offset
          do j = 1, m
            l       = indir_rhs(index(i))
            rowPtr  = rhsPtr(l)
            ind_n   = rhsPtr(l + 1) - rowPtr
            lrow  = index(i) - rowPtr - 1
            ind_i = rowPtr   * tdu  + &
                    threadID * ind_n + &
                    lrow + 1
            upd(ind_i) = upd(ind_i) + xlocal(j)

            i = i + 1
          end do
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
          call spllt_ftac(1, threadID, timer)
#endif
       else
!!! Single rhs, direct update
          j     = 1
          l     = 1
          ndblk = size(rhsPtr) - 1
          do i = offset, offset+m-1
            l = indir_rhs(index(i))
            rowPtr  = rhsPtr(l)
            ind_n   = rhsPtr(l + 1) - rowPtr
            w = zero
            do k = col, col + nelim - 1
               w = w - dest(j) * rhs(k)
               j = j + 1
            end do
            j = j + (ldd-nelim)

            lrow  = index(i) - rowPtr - 1
            ind_i = rowPtr   * tdu  + &
                    threadID * ind_n + &
                    lrow + 1
            upd(ind_i) = upd(ind_i) + w
          end do
       endif
    else
!!! Multiple rhs, BLAS 3
      call dgemm('T', 'N', m, nrhs, nelim, -one, dest, ldd, rhs(lcol), nelim, &
           zero, xlocal, m)

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
      call spllt_tic("fwd::scatter", 1, threadID, timer)
#endif

      ! Copy xlocal out
      j = 1
      rhs_blk_size    = nrhs * tdu - 1
      rhs_thread_size = threadID * nrhs
      do i = offset, offset + m - 1
        l       = indir_rhs(index(i))
        rowPtr  = rhsPtr(l)
        ind_n   = rhsPtr(l + 1) - rowPtr
        pos     = rowPtr * rhs_blk_size + ind_n * rhs_thread_size
        lrow    = pos + index(i)
        do r = 0, nrhs-1
          ind_i = r * ind_n + lrow
          upd(ind_i) = upd(ind_i) + xlocal(j + r * m)
        end do
        j = j + 1
      end do

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
      call spllt_ftac(1, threadID, timer)
#endif
    endif

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    call spllt_close_ftimer(threadID, timer)
#endif

  end subroutine slv_fwd_update_ileave



  subroutine slv_bwd_update_ileave(m, nelim, col, offset, index, dest, &
      ldd, nrhs, upd, tdu, rhs, n, rhsPtr, indir_rhs, xlocal)
    use spllt_data_mod
    use timer_mod
 !$ use omp_lib, ONLY : omp_get_thread_num
    implicit none

    integer,  intent(in)    :: m            ! number of rows in block
    integer,  intent(in)    :: nelim        ! # eliminations
                                            ! (immediate return if =0)
    integer,  intent(in)    :: col          ! start of block column variables
                                            ! in rhs
    integer,  intent(in)    :: offset       ! offset into index we start at
    integer,  intent(in)    :: index(*)
    integer,  intent(in)    :: ldd          ! leading dimension of block
    real(wp), intent(in)    :: dest(m*ldd)  ! holds block
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: tdu
    integer,  intent(in)    :: n
    integer,  intent(in)    :: rhsPtr(:)
    integer,  intent(in)    :: indir_rhs(:)
    real(wp), intent(inout) :: upd(:)
    real(wp), intent(inout) :: rhs(n*nrhs)
    real(wp), intent(out)   :: xlocal(:)

    integer   :: i, ind_i, ind_n
    integer   :: j
    integer   :: k
    integer   :: l
    integer   :: r ! right hand side loop variable
    integer   :: lrow, rowPtr, lcol, ndblk, pos
    integer   :: rhs_blk_size, rhs_thread_size
    real(wp)  :: w ! temporary work variable
    integer   :: threadID
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer
#endif

    threadID = 0
 !$ threadID = omp_get_thread_num()

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_open_timer(threadID, "slv_bwd_update_ileave", timer)
#endif

    if(nelim.eq.0) return

    l             = indir_rhs(col) 
    lcol          = rhsPtr(l) * nrhs * tdu + 1
    lcol          = lcol + threadID * nelim * nrhs
    rhs_blk_size  = nrhs - 1
    j             = 1

    ! backward substitution
    if(nrhs.eq.1) then
      if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single right-hand side, BLAS 2

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
        call spllt_tic("bwd::gather", 1, threadID, timer)
#endif
         ! Copy xlocal in
       !j = 1
       !do i = offset, offset+m-1
       !   xlocal(j) = rhs(index(i))
       !   j = j + 1
       !end do
        do i = offset, offset + m - 1
          l       = indir_rhs(index(i))
          rowPtr  = rhsPtr(l)
          pos     = rowPtr * rhs_blk_size
          ind_i   = pos + index(i)
   !!     print *, "Gather rhs(", ind_i, ") into xlocal(", j + r * m
          xlocal(j) = rhs(ind_i)
          j = j + 1
        end do
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
        call spllt_ftac(1, threadID, timer)
#endif

        call dgemv('N', nelim, m, -one, dest, ldd, xlocal, 1, one, &
             upd(lcol), 1)
      else
!!! Single right-hand side, direct update
        do i = offset, offset+m-1
          l       = indir_rhs(index(i))
          rowPtr  = rhsPtr(l)
          pos     = rowPtr * rhs_blk_size
          ind_i   = pos + index(i)

          w = rhs(ind_i)
          do k = lcol, lcol + nelim - 1
            upd(k) = upd(k) - dest(j)*w
            j = j + 1
          end do
          j = j + (ldd-nelim)
        end do
      endif
    else
!!! Multiple RHS, BLAS 3

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
      call spllt_tic("bwd::gather", 1, threadID, timer)
#endif
#if 0
      ! Copy xlocal in
      do i = offset, offset + m - 1
        l       = indir_rhs(index(i))
        rowPtr  = rhsPtr(l)
        ind_n   = rhsPtr(l + 1) - rowPtr
        pos     = rowPtr * rhs_blk_size
        lrow    = pos + index(i)
        do r = 0, nrhs-1
          ind_i = r * ind_n + lrow
 !!       print *, "Gather rhs(", ind_i, ") into xlocal(", j + r * m
         !upd(ind_i) = upd(ind_i) + xlocal(j + r * m)
          xlocal(j + r * m) = rhs(ind_i)
        end do
        j = j + 1
      end do
#else
      do r = 0, nrhs - 1
        j = 1
        do i = offset, offset + m - 1
          l       = indir_rhs(index(i))
          rowPtr  = rhsPtr(l)
          ind_n   = rhsPtr(l + 1) - rowPtr
          pos     = rowPtr * rhs_blk_size
          lrow    = pos + index(i)

          xlocal(j + r * m) = rhs(lrow + r * ind_n)
          j = j + 1
        end do
      end do
#endif
!     do i = offset, offset+m-1
!        do r = 0, nrhs-1
!           xlocal(j+r*m) = rhs(index(i)+r*ldr)
!        end do
!        j = j + 1
!     end do
      !Test change order of loop
     !do r = 0, nrhs-1
     !   j = 1
     !   do i = offset, offset+m-1
     !     xlocal(j+r*m) = rhs(index(i)+r*ldr)
     !     j = j + 1
     !   end do
     !end do

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
      call spllt_ftac(1, threadID, timer)
#endif

 !!   print *, "Call DGEMM on upd(", lcol
      call dgemm('N', 'N', nelim, nrhs, m, -one, dest, ldd, xlocal, m, &
           one, upd(lcol), nelim)
    endif

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_ftimer(threadID, timer)
#endif

  end subroutine slv_bwd_update_ileave
  


  subroutine solve_fwd_block_work_ileave(blk, rhsPtr, indir_rhs, blkm,    &
      blkn, col, offset, index, lcol, sa, nrhs, upd, tdu, rhs, n, xlocal, &
      threadID, nthread, flops)
    use spllt_data_mod
    use timer_mod
    implicit none
    integer,                  intent(in)    :: blk
    integer,                  intent(inout) :: blkm
    integer,                  intent(in)    :: blkn
    integer,                  intent(in)    :: col
    integer,                  intent(inout) :: offset
    integer,                  intent(inout) :: sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: tdu
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: nthread
    double precision,         intent(out)   :: flops
    integer,                  intent(in)    :: rhsPtr(:)
    integer,                  intent(in)    :: indir_rhs(:)
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),  target,        intent(inout) :: upd(:)
    real(wp),  target,        intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)

    integer           :: i, r, j, l
    integer           :: rowPtr_rhs, rowPtr_upd
    real(wp), pointer :: p_upd(:,:)
    real(wp), pointer :: p_rhs(:)
    double precision  :: lflops
    integer           :: nlrhs_val
#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer

    call spllt_open_timer(threadID, "solve_fwd_block_work_ileave", timer)
#endif

    flops = zero
    l     = indir_rhs(col)

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd block task reduction", 1, threadID, timer)
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_tic("fwd::reduction", 2, threadID, timer)
#endif

    ! Sum contributions to rhs
    if(nthread .eq. 1) then

      rowPtr_rhs  = rhsPtr(l) * nrhs
      rowPtr_upd  = rhsPtr(l) * nrhs * tdu + threadID * nrhs * blkn
      nlrhs_val   = blkn * nrhs
      p_rhs(1 : nlrhs_val) => rhs(rowPtr_rhs + 1 : rowPtr_rhs + nlrhs_val)

      p_rhs = p_rhs + upd(rowPtr_upd + 1 : rowPtr_upd + blkn * nrhs)
      upd(rowPtr_upd + 1 : rowPtr_upd + blkn * nrhs) = zero

#if defined(SPLLT_PROFILING_FLOP)
    lflops = blkn * nrhs
#endif
    else
      rowPtr_rhs = rhsPtr(l) * nrhs
      rowPtr_upd = rhsPtr(l) * nrhs * tdu
      nlrhs_val = blkn * nrhs
      p_rhs(1 : nlrhs_val)                  => rhs(rowPtr_rhs + 1 : &
        rowPtr_rhs + nlrhs_val)
      p_upd(1 : nlrhs_val, 0 : nthread - 1) => upd(rowPtr_upd + 1 : &
        rowPtr_upd + nlrhs_val * tdu)

      do j = 0, nthread - 1
        p_rhs(1 : nlrhs_val) = p_rhs(1 : nlrhs_val) + p_upd(1 : nlrhs_val, j)
      end do
      p_upd = zero

#if defined(SPLLT_PROFILING_FLOP)
      lflops = blkn * nrhs * nthread
#endif
    end if
#if defined(SPLLT_PROFILING_FLOP)
    flops = flops + lflops
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_ftac(2, threadID, timer)
#endif
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(1, threadID, timer, lflops)
#endif


    ! Perform triangular solve
    call slv_solve_ileave(blkn, blkn, rhsPtr(l) * nrhs + 1, &
      lcol(sa : sa + blkn * blkn - 1), 'Transpose    ', 'Non-unit', nrhs, rhs)

#if defined(SPLLT_PROFILING_FLOP)
    flops = blkn * blkn * nrhs
#endif
    offset = offset + blkn

    ! Deal with any left over trapezoidal part of diagonal block
    blkm = blkm - blkn
    if(blkm .gt. 0) then
      sa = sa + blkn * blkn
      call slv_fwd_update_ileave(blkm, blkn, col, offset, index,  &
        lcol(sa : sa + blkn * blkm - 1), blkn, nrhs,              &
        upd, tdu, rhs, n, rhsPtr, indir_rhs, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (blkn * nrhs * blkm)
#endif
    endif

#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_timer(threadID, timer, flops)
#endif
  end subroutine solve_fwd_block_work_ileave



  subroutine solve_fwd_update_work_ileave(blk, rhsPtr, indir_rhs, blkm,   &
      blkn, col, offset, index, lcol, blk_sa, nrhs, upd, tdu, rhs, n,     &
      xlocal, threadID, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(in)    :: blk
    integer,                  intent(in)    :: blkm
    integer,                  intent(in)    :: blkn
    integer,                  intent(in)    :: col
    integer,                  intent(in)    :: offset
    integer,                  intent(in)    :: blk_sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: tdu
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: threadID
    double precision,         intent(out)   :: flops
    integer,                  intent(in)    :: rhsPtr(:)
    integer,                  intent(in)    :: indir_rhs(:)
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)

    call slv_fwd_update_ileave(blkm, blkn, col, offset,             &
      index, lcol(blk_sa : blk_sa + blkn * blkm - 1), blkn, nrhs,   &
      upd, tdu, rhs, n, rhsPtr, indir_rhs, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (blkm * nrhs * blkn)
#endif

  end subroutine solve_fwd_update_work_ileave



  subroutine solve_fwd_node_ileave(nrhs, rhs, n, fkeep, node, &
      xlocal, rhs_local, tdu, task_manager, trace_id)
    use spllt_data_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n
    integer,                    intent(in)    :: node 
    integer,                    intent(in)    :: tdu
    real(wp),                   intent(inout) :: rhs(n*nrhs)
    real(wp),                   intent(inout) :: xlocal(:,:)
    real(wp),                   intent(inout) :: rhs_local(:)
    class(task_manager_base ),  intent(inout) :: task_manager
    integer, optional,          intent(in)    :: trace_id

    integer                   :: sa, en
    integer                   :: numcol, numrow ! #col/row in node 
    integer                   :: nc, nr         ! #block-col/block-row in node
    integer                   :: jj, ii
    integer                   :: dblk           ! Diagonal index 
    integer                   :: s_nb           ! Block size in node
    integer                   :: blk            ! Block index
    integer                   :: fwd_block_id
    integer                   :: fwd_update_id
    type(spllt_timer_t), save :: timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_fwd_node_ileave", timer)
#endif

    if(present(trace_id)) then
      fwd_block_id  = trace_id
      fwd_update_id = trace_id
    else
      fwd_block_id  = task_manager%trace_ids(trace_fwd_block_pos)
      fwd_update_id = task_manager%trace_ids(trace_fwd_update_pos)
    end if

    ! Get node info
    s_nb   = fkeep%nodes(node)%nb
    sa     = fkeep%nodes(node)%sa
    en     = fkeep%nodes(node)%en
    numcol = en - sa + 1
    numrow = size(fkeep%nodes(node)%index)
    nc     = (numcol-1) / s_nb + 1
    nr     = (numrow-1) / s_nb + 1
    
    ! Get first diag block in node
    dblk = fkeep%nodes(node)%blk_sa

    ! Loop over block columns
    do jj = 1, nc
       
      !
      ! Forward solve with block on diagoanl
      !
 !!   print *, "Submit fwd block ", dblk
      call spllt_tic("submit fwd block", 1, task_manager%workerID, timer)
      call task_manager%solve_fwd_block_il_task(dblk, nrhs, rhs_local,  &
        tdu, rhs, n, xlocal, fkeep, fwd_block_id)
      call spllt_tac(1, task_manager%workerID, timer)

      do ii = jj+1, nr

        blk = dblk+ii-jj

        !
        ! Forward update with off-diagonal
        !
 !!     print *, "Submit fwd update ", blk
        call spllt_tic("submit fwd update", 2, task_manager%workerID, timer)
        call task_manager%solve_fwd_update_il_task(blk, node, nrhs, rhs_local,&
          tdu, rhs, n, xlocal, fkeep, fwd_update_id)
        call spllt_tac(2, task_manager%workerID, timer)

      end do
      
      ! Update diag block in node          
      dblk = fkeep%bc(dblk)%last_blk + 1
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_fwd_node_ileave


  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Kernel of solve steps
  !
  subroutine solve_bwd_block_work_ileave(blk, rhsPtr, indir_rhs, blkm,    &
      blkn, col, offset, index, lcol, sa, nrhs, upd, tdu, rhs, n, xlocal, &
      threadID, nthread, flops)
    use spllt_data_mod
    use timer_mod
    implicit none
    integer,                  intent(in)    :: blk
    integer,                  intent(inout) :: blkm
    integer,                  intent(in)    :: blkn
    integer,                  intent(in)    :: col
    integer,                  intent(inout) :: offset
    integer,                  intent(inout) :: sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: tdu
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: nthread
    double precision,         intent(out)   :: flops
    integer,                  intent(in)    :: rhsPtr(:)
    integer,                  intent(in)    :: indir_rhs(:)
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),  target,        intent(inout) :: upd(:)
    real(wp),  target,        intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)

    integer           :: i, r, j, l
    integer           :: rowPtr_rhs, rowPtr_upd
    real(wp), pointer :: p_upd(:,:)
    real(wp), pointer :: p_rhs(:)
    double precision  :: lflops
    integer           :: nlrhs_val
#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer

    call spllt_open_timer(threadID, "solve_bwd_block_work_ileave", timer)
#endif

    flops = zero
    l     = indir_rhs(col)

    ! Perform retangular update from diagonal block
    if(blkm .gt. blkn) then
      call slv_bwd_update_ileave(blkm - blkn, blkn, col, offset + blkn,   &
        index, lcol(sa + blkn * blkn : sa + blkn * blkm - 1), blkn, nrhs, &
        upd, tdu, rhs, n, rhsPtr, indir_rhs, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (blkn * nrhs * (blkm - blkn))
#endif
    endif

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd block task reduction", 1, threadID, timer)
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_tic("bwd::reduction", 2, threadID, timer)
#endif
#if 0
    ! Sum contributions to rhs
    if(nthread .eq. 1) then
!     print *, "Nthread given to solve_bwd_block_work : ", nthread
!     print *, "Lighter reduction"
     !do r = 0, nrhs-1
     !  j = threadID + 1
     !  do i = col + r*ldr, col+n-1 + r*ldr
     !    rhs(i)    = rhs(i) + upd(i, j)
     !    upd(i,j)  = zero ! Reset in case of next fwd solve
     !  end do
     !end do
#if defined(SPLLT_PROFILING_FLOP)
    lflops = blkn * nrhs
#endif
    else
!     print *, "Nthread given to solve_bwd_block_work : ", nthread
     !do r = 0, nrhs-1
     !  do j = 1, nthread
     !    do i = col + r*ldr, col+n-1 + r*ldr
     !      rhs(i)    = rhs(i) + upd(i, j)
     !      upd(i,j)  = zero ! Reset in case of next fwd solve
     !    end do
     !  end do
     !end do
#if defined(SPLLT_PROFILING_FLOP)
    lflops = blkn * nrhs * nthread
#endif
#else
    ! Sum contributions to rhs
    if(nthread .eq. 1) then

      rowPtr_rhs  = rhsPtr(l) * nrhs
      rowPtr_upd  = rhsPtr(l) * nrhs * tdu + threadID * nrhs * blkn
      nlrhs_val   = blkn * nrhs
      p_rhs(1 : nlrhs_val) => rhs(rowPtr_rhs + 1 : rowPtr_rhs + nlrhs_val)

      j = threadID
 !!   print *, "Reduce upd(", rowPtr_upd + 1 + j * nlrhs_val, ",",  &
 !!     rowPtr_upd + (j + 1) * nlrhs_val, ") into rhs(",        &
 !!     rowPtr_rhs + 1, ',', rowPtr_rhs + nlrhs_val
      p_rhs = p_rhs + upd(rowPtr_upd + 1 : rowPtr_upd + blkn * nrhs)
      upd(rowPtr_upd + 1 : rowPtr_upd + blkn * nrhs) = zero

#if defined(SPLLT_PROFILING_FLOP)
    lflops = blkn * nrhs
#endif
    else
      rowPtr_rhs = rhsPtr(l) * nrhs
      rowPtr_upd = rhsPtr(l) * nrhs * tdu
      nlrhs_val = blkn * nrhs
      p_rhs(1 : nlrhs_val)                  => rhs(rowPtr_rhs + 1 : &
        rowPtr_rhs + nlrhs_val)
      p_upd(1 : nlrhs_val, 0 : nthread - 1) => upd(rowPtr_upd + 1 : &
        rowPtr_upd + nlrhs_val * tdu)

      do j = 0, nthread - 1
 !!     print *, "Reduce upd(", rowPtr_upd + 1 + j * nlrhs_val, ",",  &
 !!       rowPtr_upd + (j + 1) * nlrhs_val, ") into rhs(",        &
 !!       rowPtr_rhs + 1, ',', rowPtr_rhs + nlrhs_val
        p_rhs(1 : nlrhs_val) = p_rhs(1 : nlrhs_val) + p_upd(1 : nlrhs_val, j)
      end do
      p_upd = zero

#if defined(SPLLT_PROFILING_FLOP)
    lflops = blkn * nrhs * nthread
#endif
#endif
    end if
#if defined(SPLLT_PROFILING_FLOP)
    flops = flops + lflops
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_ftac(2, threadID, timer)
#endif
#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(1, threadID, timer, lflops)
#endif

    ! Perform triangular solve
 !! print *, "Call solve starting col =", rhsPtr(l) * nrhs + 1
    call slv_solve_ileave(blkn, blkn, rhsPtr(l) * nrhs + 1, &
      lcol(sa : sa + blkn * blkn - 1), 'Non-Transpose', 'Non-unit', nrhs, rhs)

#if defined(SPLLT_PROFILING_FLOP)
    flops = flops + blkn * blkn * nrhs
#endif

#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_timer(threadID, timer, flops)
#endif
  end subroutine solve_bwd_block_work_ileave



  subroutine solve_bwd_update_work_ileave(blk, rhsPtr, indir_rhs, blkm, &
      blkn, col, offset, index, lcol, blk_sa, nrhs, upd, tdu, rhs, n,   &
      xlocal, threadID, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(in)    :: blk
    integer,                  intent(in)    :: blkm
    integer,                  intent(in)    :: blkn
    integer,                  intent(in)    :: col
    integer,                  intent(in)    :: offset
    integer,                  intent(in)    :: blk_sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: tdu
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: threadID
    double precision,         intent(out)   :: flops
    integer,                  intent(in)    :: rhsPtr(:)
    integer,                  intent(in)    :: indir_rhs(:)
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)

    call  slv_bwd_update_ileave(blkm, blkn, col, offset,          &
      index, lcol(blk_sa : blk_sa + blkn * blkm - 1), blkn, nrhs, &
      upd, tdu, rhs, n, rhsPtr, indir_rhs, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (blkn * nrhs * blkm)
#endif

  end subroutine solve_bwd_update_work_ileave



  subroutine solve_bwd_node_ileave(nrhs, rhs, n, fkeep, node, xlocal, &
      rhs_local, tdu, task_manager, trace_id)
    use spllt_data_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n
    integer,                    intent(in)    :: tdu
    integer,                    intent(in)    :: node 
    real(wp),                   intent(inout) :: rhs(n*nrhs)
    real(wp),                   intent(inout) :: xlocal(:,:)
    real(wp),                   intent(inout) :: rhs_local(:)
    class(task_manager_base ),  intent(inout) :: task_manager
    integer, optional,          intent(in)    :: trace_id

    integer                 :: sa, en
    integer                 :: numcol, numrow ! #column/row in node 
    integer                 :: nc, nr         ! #block-column/block-row in node
    integer                 :: jj, ii
    integer                 :: dblk
    integer                 :: s_nb           ! Block size in node
    integer                 :: blk
    integer                 :: bwd_block_id
    integer                 :: bwd_update_id
    type(spllt_timer_t), save :: timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_bwd_node_ileave", timer)
#endif

    if(present(trace_id)) then
      bwd_block_id  = trace_id
      bwd_update_id = trace_id
    else
      bwd_block_id  = task_manager%trace_ids(trace_bwd_block_pos)
      bwd_update_id = task_manager%trace_ids(trace_bwd_update_pos)
    end if

    ! Get node info
    s_nb   = fkeep%nodes(node)%nb
    sa     = fkeep%nodes(node)%sa
    en     = fkeep%nodes(node)%en
    numcol = en - sa + 1
    numrow = size(fkeep%nodes(node)%index)
    nc     = (numcol-1) / s_nb + 1
    nr     = (numrow-1) / s_nb + 1 

    ! Get first diag block in node
    dblk = fkeep%bc(fkeep%nodes(node)%blk_en)%dblk

    ! Loop over block columns
    do jj = nc, 1, -1

      do ii = nr, jj+1, -1
        
        blk = dblk+ii-jj ! Block index

        !
        ! Backward update with block on diagoanl
        !
 !!     print *, "Submit bwd update ", blk
        call spllt_tic("submit bwd update", 1, task_manager%workerID, timer)
        call task_manager%solve_bwd_update_il_task(blk, node, nrhs, rhs_local, &
          tdu, rhs, n, xlocal, fkeep, bwd_update_id)
        call spllt_tac(1, task_manager%workerID, timer)
       !!$omp taskwait
      end do

      !
      ! Backward solve with block on diagoanl
      !
 !!   print *, "Submit bwd block ", dblk
      call spllt_tic("submit bwd block", 2, task_manager%workerID, timer)
      call task_manager%solve_bwd_block_il_task(dblk, nrhs, rhs_local, tdu, &
        rhs, n, xlocal, fkeep, bwd_block_id)
      call spllt_tac(2, task_manager%workerID, timer)
     !!$omp taskwait
     
      ! Update diag block in node       
      if (jj .gt. 1) dblk = fkeep%bc(dblk-1)%dblk
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_bwd_node_ileave


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! NEW
  !
  !
  subroutine slv_solve_ileave2(n, nelim, dest, trans, unit, nrhs, rhs, ldr)
    use spllt_data_mod
    implicit none

    integer,            intent(in)    :: n      ! leading dimension of 
                                                ! diag. block
    integer,            intent(in)    :: nelim  ! number eliminations 
                                                ! (immediate return if =0)
    real(wp),           intent(in)    :: dest(*)! holds destination block
    character(len=13),  intent(in)    :: trans  ! set to :
                                                ! 'Transpose    ' 
                                                ! for forward substitution 
                                                ! 'Non-Transpose' 
                                                ! for back substitution
    character(len=8),   intent(in)    :: unit   ! set to 
                                                ! 'Non-Unit' 
                                                ! for positive-definite case
    integer,            intent(in)    :: nrhs   ! number of right-hand sides
    real(wp),           intent(inout) :: rhs(*)
    integer,            intent(in)    :: ldr

    if(nelim.eq.0) return

    if(nrhs.eq.1) then
     !print *, "Rhs :", rhs(1 : nelim)
      call dtrsv('Upper', trans, unit, nelim, dest, n, rhs, 1)
     !print *, "dest", dest(1 : nelim * nelim)
     !print *, "Result of solving :", rhs(1 : nelim)
    else
     !print *, "Rhs :", rhs(1 : nelim * nrhs)
      call dtrsm('Left', 'Upper', trans, unit, nelim, nrhs, one, &
           dest, n, rhs, ldr)
     !print *, "dest", dest(1 : nelim * nelim)
     !print *, "Result of solving :", rhs(1 : nelim * nrhs)
    endif

  end subroutine slv_solve_ileave2



  subroutine slv_fwd_update_ileave2(m, nelim, dest, ldd, &
      n, nrhs, rhs, ldr, xlocal, ldx)
    use spllt_data_mod
    use timer_mod
 !$ use omp_lib, ONLY : omp_get_thread_num
    implicit none

    integer,  intent(in)    :: m              ! number of rows in block
    integer,  intent(in)    :: nelim          ! # eliminations 
                                              ! (immediate return if =0)
    real(wp), intent(in)    :: dest(m*ldd)    ! holds destination block
    integer,  intent(in)    :: ldd            ! leading dimension of block
    integer,  intent(in)    :: n
    integer,  intent(in)    :: nrhs
    real(wp), intent(in)    :: rhs(*)    ! rhs vector
    integer,  intent(in)    :: ldr
    real(wp), intent(out)   :: xlocal(*)
    integer,  intent(in)    :: ldx            ! leading dimension of block

    integer   :: i
    integer   :: j
    integer   :: k
    real(wp)  :: w ! temporary work value
    integer   :: threadID
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    type(spllt_timer_t), save :: timer
#endif

    threadID = 0
 !$ threadID = omp_get_thread_num()

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    call spllt_open_timer(threadID, "slv_fwd_update_ileave2", timer)
#endif

    if(nelim.eq.0) return

    ! forward substitution
    if(nrhs.eq.1) then
      !print *, "local RHS", rhs(1:nelim) 
      !print *, "Result before", xlocal(1 : m)
       if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single rhs, BLAS 2

          call dgemv('T', nelim, m, -one, dest, ldd, rhs, 1, one, xlocal, 1)

       else
!!! Single rhs, direct update
      !  print *, "Use homemade loop"
      !  print *, "local RHS", rhs(1:nelim) 
      !  print *, "Dest", dest(1 : m * nelim)
      !  print *, "Result before", xlocal(1 : m)
         j = 1
         do i = 1, m
           w = zero
           do k = 1, nelim
             w = w - dest(j)*rhs(k)
             j = j + 1
           end do
           j = j + (ldd-nelim)
           xlocal(i) = xlocal(i) + w
         end do
       endif
      !print *, "Result DGEMM", xlocal(1 : m)
    else
!!! Multiple rhs, BLAS 3
      call dgemm('T', 'N', m, nrhs, nelim, -one, dest, ldd, rhs, ldr, &
           one, xlocal, ldx)
    endif

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    call spllt_close_ftimer(threadID, timer)
#endif

  end subroutine slv_fwd_update_ileave2


  subroutine slv_bwd_update_ileave2(m, nelim, dest, ldd, &
      n, nrhs, xlocal, ldx, rhs, ldr)
    use spllt_data_mod
    use timer_mod
 !$ use omp_lib, ONLY : omp_get_thread_num
    implicit none

    integer,  intent(in)    :: m            ! number of rows in block
    integer,  intent(in)    :: nelim        ! # eliminations
                                            ! (immediate return if =0)
    integer,  intent(in)    :: ldd          ! leading dimension of block
    real(wp), intent(in)    :: dest(m*ldd)  ! holds block
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: n
    integer,  intent(in)    :: ldx
    integer,  intent(in)    :: ldr
    real(wp), intent(inout) :: rhs(*)       ! Update of the RHS
    real(wp), intent(in)    :: xlocal(*)    ! Data already gather in it

    integer   :: i
    integer   :: j
    integer   :: k
    real(wp)  :: w ! temporary work variable
    integer   :: threadID
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer
#endif

    threadID = 0
 !$ threadID = omp_get_thread_num()

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_open_timer(threadID, "slv_bwd_update_ileave2", timer)
#endif

    if(nelim.eq.0) return

    ! backward substitution
    if(nrhs.eq.1) then
     !print *, "local RHS", xlocal(1 : m) 
     !print *, "Result before", rhs(1 : nelim)
      if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single right-hand side, BLAS 2

        call dgemv('N', nelim, m, -one, dest, ldd, xlocal, 1, one, &
             rhs, 1)
      else
!!! Single right-hand side, direct update
        j = 1
        do i = 1, m
          w = xlocal(i)
          do k = 1, nelim
            rhs(k) = rhs(k) - dest(j)*w
            j = j + 1
          end do
          j = j + (ldd-nelim)
        end do
      endif
     !print *, "Result DGEMM", rhs(1 : nelim)
    else
!!! Multiple RHS, BLAS 3
      call dgemm('N', 'N', nelim, nrhs, m, -one, dest, ldd, xlocal, ldx, &
           one, rhs, ldr)
    endif

#if defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_ftimer(threadID, timer)
#endif

  end subroutine slv_bwd_update_ileave2


#if 0
  subroutine solve_fwd_block_work_ileave2(blkm, blkn, lcol, sa, nrhs, y, &
      ldy, xlocal, ldx, threadID, flops)
    use spllt_data_mod
    use timer_mod
    implicit none
    integer,                  intent(inout) :: blkm
    integer,                  intent(in)    :: blkn
    integer,                  intent(inout) :: sa
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: ldy
    integer,                  intent(in)    :: ldx
    integer,                  intent(in)    :: threadID
    double precision,         intent(out)   :: flops
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: y(:,:)
    real(wp),                 intent(inout) :: xlocal(:)

    integer           :: i
#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer

    call spllt_open_timer(threadID, "solve_fwd_block_work_ileave2", timer)
#endif

    ! Perform triangular solve
    call slv_solve_ileave2(blkn, blkn, lcol(sa : sa + blkn * blkn - 1), &
      'Transpose    ', 'Non-unit', nrhs, y, ldy)

#if defined(SPLLT_PROFILING_FLOP)
    flops = blkn * blkn * nrhs
#endif

    ! Deal with any left over trapezoidal part of diagonal block
    blkm = blkm - blkn
    if(blkm .gt. 0) then
      sa = sa + blkn * blkn
      call slv_fwd_update_ileave2(blkm, blkn, lcol(sa : sa + blkn * blkm - 1), &
        blkn, blkn, nrhs, y, ldy, xlocal, ldx)

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (blkn * nrhs * blkm)
#endif
    endif

#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_timer(threadID, timer, flops)
#endif
  end subroutine solve_fwd_block_work_ileave2



  subroutine solve_fwd_update_work_ileave2(blkm, blkn, lcol, blk_sa, &
      n, nrhs, rhs, ldr, xlocal, ldx, threadID, flops)
    use spllt_data_mod
    implicit none
    integer,                  intent(in)    :: blkm
    integer,                  intent(in)    :: blkn
    integer,                  intent(in)    :: blk_sa
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: nrhs
    integer,                  intent(in)    :: ldr
    integer,                  intent(in)    :: ldx
    integer,                  intent(in)    :: threadID ! useless ?
    double precision,         intent(out)   :: flops
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:)

    call slv_fwd_update_ileave2(blkm, blkn,                              &
      lcol(blk_sa : blk_sa + blkn * blkm - 1), blkn, n, nrhs, rhs, ldr, &
      xlocal, ldx)

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (blkm * nrhs * blkn)
#endif

  end subroutine solve_fwd_update_work_ileave2
#endif

  subroutine solve_fwd_node_ileave2(nrhs, rhs, n, fkeep, node, &
      task_manager, trace_id)
    use spllt_data_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n
    integer,                    intent(in)    :: node 
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    class(task_manager_base ),  intent(inout) :: task_manager
    integer, optional,          intent(in)    :: trace_id

    integer                   :: sa, en
    integer                   :: numcol, numrow ! #col/row in node 
    integer                   :: nc, nr         ! #block-col/block-row in node
    integer                   :: jj, ii
    integer                   :: dblk           ! Diagonal index 
    integer                   :: s_nb           ! Block size in node
    integer                   :: blk            ! Block index
    integer                   :: fwd_block_id
    integer                   :: fwd_update_id
    type(spllt_timer_t), save :: timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_fwd_node_ileave2", timer)
#endif

    if(present(trace_id)) then
      fwd_block_id  = trace_id
      fwd_update_id = trace_id
    else
      fwd_block_id  = task_manager%trace_ids(trace_fwd_block_pos)
      fwd_update_id = task_manager%trace_ids(trace_fwd_update_pos)
    end if

    ! Get node info
    s_nb   = fkeep%nodes(node)%nb
    sa     = fkeep%nodes(node)%sa
    en     = fkeep%nodes(node)%en
    numcol = en - sa + 1
    nc     = ceiling(numcol / real(s_nb))
    
    ! Get first diag block in node
    dblk = fkeep%nodes(node)%sblk_sa
   !call print_node_solve(fkeep, node)

    ! Loop over block columns
    do jj = 1, nc
       
      !
      ! Forward solve with block on diagoanl
      !
 !!   print *, "Submit fwd block ", dblk
      call spllt_tic("submit fwd block", 1, task_manager%workerID, timer)
      call task_manager%solve_fwd_block_il2_task(dblk, nrhs, n, rhs, &
        fkeep, trace_id)
     !call task_manager%solve_fwd_block_il_task(dblk, nrhs, rhs_local,  &
     !  tdu, rhs, n, xlocal, fkeep, fwd_block_id)
      call spllt_tac(1, task_manager%workerID, timer)

      do blk = dblk + 1, fkeep%sbc(dblk)%last_blk

       !blk = dblk+ii-jj

        !
        ! Forward update with off-diagonal
        !
 !!     print *, "Submit fwd update ", blk
        call spllt_tic("submit fwd update", 2, task_manager%workerID, timer)
        call task_manager%solve_fwd_update_il2_task(blk, node, nrhs, &
            n, rhs, fkeep, trace_id)
       !call task_manager%solve_fwd_update_il2_task(blk, node, nrhs, rhs_local,&
       !  tdu, rhs, n, xlocal, fkeep, fwd_update_id)
        call spllt_tac(2, task_manager%workerID, timer)
      end do
      
      ! Update diag block in node          
      dblk = fkeep%sbc(dblk)%last_blk + 1
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_fwd_node_ileave2


  subroutine solve_bwd_node_ileave2(nrhs, rhs, n, fkeep, node, &
      task_manager, trace_id)
    use spllt_data_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n
    integer,                    intent(in)    :: node 
    real(wp),                   intent(inout) :: rhs(n, nrhs)
    class(task_manager_base ),  intent(inout) :: task_manager
    integer, optional,          intent(in)    :: trace_id

    integer                 :: sa, en
    integer                 :: numcol, numrow ! #column/row in node 
    integer                 :: nc, nr         ! #block-column/block-row in node
    integer                 :: jj, ii
    integer                 :: dblk
    integer                 :: s_nb           ! Block size in node
    integer                 :: blk
    integer                 :: bwd_block_id
    integer                 :: bwd_update_id
    type(spllt_timer_t), save :: timer

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_open_timer(task_manager%workerID, "solve_bwd_node_ileave2", timer)
#endif

    if(present(trace_id)) then
      bwd_block_id  = trace_id
      bwd_update_id = trace_id
    else
      bwd_block_id  = task_manager%trace_ids(trace_bwd_block_pos)
      bwd_update_id = task_manager%trace_ids(trace_bwd_update_pos)
    end if

    ! Get node info
    s_nb   = fkeep%nodes(node)%nb
    sa     = fkeep%nodes(node)%sa
    en     = fkeep%nodes(node)%en
    numcol = en - sa + 1
    nc     = ceiling(numcol / real(s_nb))

    ! Get first diag block in node
    dblk = fkeep%sbc(fkeep%nodes(node)%sblk_en)%dblk

    ! Loop over block columns
    do jj = nc, 1, -1

 !!   print *, "update from", fkeep%sbc(dblk)%last_blk, "to", dblk + 1
      do blk = fkeep%sbc(dblk)%last_blk, dblk + 1, -1
        
        !
        ! Backward update with block on diagoanl
        !
 !!     print *, "Submit bwd update ", blk
       !print *, "blk", blk, "y INPUT", fkeep%sbc(blk)%p_upd
        call spllt_tic("submit bwd update", 1, task_manager%workerID, timer)
        call task_manager%solve_bwd_update_il2_task(blk, node, nrhs, & 
          n, rhs, fkeep, trace_id)
        call spllt_tac(1, task_manager%workerID, timer)
       !print *, "blk", blk, "x associated", fkeep%sbc(blk)%p_upd
      end do

      !
      ! Backward solve with block on diagoanl
      !
 !!   print *, "Submit bwd block ", dblk
     !print *, "dblk", dblk, "y INPUT", fkeep%sbc(dblk)%p_upd
      call spllt_tic("submit bwd block", 2, task_manager%workerID, timer)
      call task_manager%solve_bwd_block_il2_task(dblk, nrhs, &
        n, rhs, fkeep, trace_id)
      call spllt_tac(2, task_manager%workerID, timer)
     !print *, "dblk", dblk, "x associated", fkeep%sbc(dblk)%p_upd
     
      ! Update diag block in node       
      if (jj .gt. 1) dblk = fkeep%sbc(dblk-1)%dblk
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_bwd_node_ileave2

end module spllt_solve_kernels_mod
