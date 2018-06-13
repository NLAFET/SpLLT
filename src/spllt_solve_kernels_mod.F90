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
   !  do r = 0, nrhs - 1
   !    print *, "SOLVE org RHS", rhs(col + r * ldr : col + n - 1 + r * ldr)
   !  end do
      call dtrsm('Left', 'Upper', trans, unit, nelim, nrhs, one, &
           dest, n, rhs(col), ldr)
   !  do r = 0, nrhs - 1
   !    print *, "SOLVE upd RHS", rhs(col + r * ldr : col + n - 1 + r * ldr)
   !  end do
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
      !$omp taskwait

      do ii = jj+1, nr

        blk = dblk+ii-jj

        !
        ! Forward update with off-diagonal
        !
        call spllt_tic("submit fwd update", 2, task_manager%workerID, timer)
        call task_manager%solve_fwd_update_task(blk, node, nrhs, rhs_local,&
          rhs, ldr, xlocal, fkeep, fwd_update_id)
        call spllt_tac(2, task_manager%workerID, timer)
        !$omp taskwait

      end do
      
      ! Update diag block in node          
      dblk = fkeep%bc(dblk)%last_blk + 1
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_fwd_node


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! InterLeave subroutines



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
      call dtrsv('Upper', trans, unit, nelim, dest, n, rhs(col), 1)
    else
      call dtrsm('Left', 'Upper', trans, unit, nelim, nrhs, one, &
           dest, n, rhs(col), nelim)
    endif

  end subroutine slv_solve_ileave



  subroutine slv_fwd_update_ileave(m, nelim, col, offset, index, dest, &
      ldd, nrhs, upd, tdu, rhs, n, bdr, rhsPtr, xlocal)
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
    integer,  intent(in)    :: bdr
    integer,  intent(in)    :: rhsPtr(:)
    real(wp), intent(inout) :: upd(:)         ! vector to update
    real(wp), intent(in)    :: rhs(n*nrhs)    ! rhs vector
    real(wp), intent(out)   :: xlocal(*)

    integer   :: i, ind_i, ind_n
    integer   :: j
    integer   :: k
    integer   :: blk_id
    integer   :: r ! right hand side loop variable
    integer   :: l, l_sa, ll
    integer   :: nrow
    integer   :: lrow, rowPtr, lcol, ndblk
    real(wp)  :: w ! temporary work value
    !%%%  integer :: t_start, t_end, this_thread
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    type(spllt_timer_t), save :: timer
    integer                   :: threadID

    threadID = 0
 !$ threadID = omp_get_thread_num()

    call spllt_open_timer(threadID, "slv_fwd_update_ileave", timer)
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
             upd(index(i) * tdu) = upd(index(i) * tdu) + xlocal(j)
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
             do k = col, col + nelim - 1
                w = w - dest(j) * rhs(k)
                j = j + 1
             end do
             j = j + (ldd-nelim)
             upd(index(i) * tdu) = upd(index(i) * tdu) + w
          end do
       endif
    else
      lcol  = -1
      ndblk = size(rhsPtr) - 1
      do l = 1, ndblk
        if(col .eq. rhsPtr(l) + 1) then
          lcol = rhsPtr(l) * nrhs + 1
          exit
        end if
      end do
      if(lcol .eq. -1) print *, "/!\ Error ", col, "not in rhsPtr", rhsPtr(:) + 1

!!! Multiple rhs, BLAS 3
      call dgemm('T', 'N', m, nrhs, nelim, -one, dest, ldd, rhs(lcol), nelim, &
           zero, xlocal, m)

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
      call spllt_tic("fwd::scatter", 1, threadID, timer)
#endif
 !!   print *, "Compute scatter from xlocal to upd"
 !!   print *, "Local unknowns from", offset, 'to', offset + m - 1

      ! Place the pointer on the lowest block in RHS
      rowPtr  = -1
      ndblk   = size(rhsPtr) - 1

      do l = 1, ndblk
        if(index(offset) .le. rhsPtr(l + 1)) then
          rowPtr  = rhsPtr(l)
          ind_n   = rhsPtr(l + 1) - rowPtr
          l_sa    = l
          exit
        end if
      end do

      if(rowPtr .eq. -1) print *, "/!\ Error ", index(offset), "not in rhsPtr", rhsPtr(:) + 1

 !!   print *, "Found ", index(offset), "lower or equal", rhsPtr(l) + 1, "with l=", l
 !!   print *, "ThreadID", threadID

      ! Copy xlocal out
      do r = 0, nrhs-1
       !print *, "RHS ", r
        !Reset variable for next iteration
        i = offset
        l = l_sa
        rowPtr = rhsPtr(l)
        ind_n   = rhsPtr(l + 1) - rowPtr
 !!     print *, "********************"
        do j = 1 + r * m, (r + 1) * m
 !!       print *, "index(", i, ") = ", index(i), "starts in pos", rowPtr
          lrow  = index(i) - rowPtr - 1
 !!       print *, "lrow = ", lrow, "ind_n=", ind_n
          ind_i = rowPtr * nrhs * tdu + r * ind_n * tdu + lrow * tdu + threadID + 1
 !!       print *, "Update upd( ", ind_i, ') with xlocal(', j, ')'

          upd(ind_i) = upd(ind_i) + xlocal(j)

          if(j .ne. (r + 1) * m) then
            i = i + 1

 !!         print *, "TRY : index(", i, ") = ", index(i), "starts with", l
            rowPtr = -1
            do ll = l, ndblk
              if(index(i) .le. rhsPtr(ll + 1)) then
 !!             if(ll .gt. l) then
 !!               print *, "Shifted l to ", ll, "for", index(i), "into", rhsPtr(ll) + 1, "with l=", ll
 !!             end if
                l       = ll
                rowPtr  = rhsPtr(l)
                ind_n   = rhsPtr(l + 1) - rowPtr
                exit
              end if
            end do
            if(rowPtr .eq. -1) print *, "/!\ Error ", index(i), "not in rhsPtr", rhsPtr(:) + 1
          end if
        end do
      end do
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
      call spllt_ftac(1, threadID, timer)
#endif
    endif

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    call spllt_close_ftimer(threadID, timer)
#endif

  end subroutine slv_fwd_update_ileave
  


  subroutine solve_fwd_block_work_ileave(blk, rhsPtr, blkm, blkn, col,  &
      offset, index, lcol, sa, nrhs, upd, ldu, bdu, tdu, rhs, n, ldr,   &
      bdr, xlocal, threadID, nthread, flops)
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
    integer,                  intent(in)    :: ldu
    integer,                  intent(in)    :: bdu
    integer,                  intent(in)    :: tdu
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: ldr
    integer,                  intent(in)    :: bdr
    integer,                  intent(in)    :: threadID
    integer,                  intent(in)    :: nthread
    double precision,         intent(out)   :: flops
    integer,                  intent(in)    :: rhsPtr(:)
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)

    integer           :: i, r, j, blk_id, ndblk, l
    integer           :: rhs_shift, upd_shift
    integer           :: row_rhs, row_upd
    logical           :: found
    double precision  :: lflops, temp
#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    type(spllt_timer_t), save :: timer

    call spllt_open_timer(threadID, "solve_fwd_block_work_ileave", timer)
#endif

    flops   = zero
    found = .false.
    ndblk = size(rhsPtr) - 1
    do l = 1, ndblk
      if(col .eq. rhsPtr(l) + 1) then
        found = .true.
        exit
      end if
    end do
    if(.not. found) print *, "/!\ Error ", col, "not in rhsPtr", rhsPtr(:) + 1

 !! print *, "Consider block ", l, "from", rhsPtr(l)

#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("fwd block task reduction", 1, threadID, timer)
#endif
#if defined(SPLLT_SOLVE_KERNEL_GATHER)
call spllt_tic("fwd::reduction", 2, threadID, timer)
#endif
 !! print *, "compute reduction"
 !! print *, "___ rhsPtr(", l, ")=", rhsPtr(l)
 !! print *, "___ bdr ", bdr, "and bdu", bdu, "tdu", tdu
    ! Sum contributions to rhs
    if(nthread .eq. 1) then
      do r = 0, nrhs - 1
 !!     print *, "RHS ", r
       !row_rhs = rhs_shift + r * blkn
       !row_upd = upd_shift + r * blkn * nthread
        row_rhs = rhsPtr(l) * nrhs + r * bdr
        row_upd = rhsPtr(l) * nrhs  * tdu + r * bdu
 !!     print *, "Row Pointer in RHS : ", row_rhs, ", and in UPD", row_upd
        do i = 1, blkn
 !!       print *, "Row ", i, "Reduce ", row_rhs + i, 'with one thread', &
 !!         row_upd + i * tdu + threadID
          rhs(row_rhs + i) = rhs(row_rhs + i) + upd(row_upd + i * tdu + threadID)
        ! rhs(blk_id * ldr + i) = rhs(blk_id * ldr + i) + &
        !   upd(blk_id * ldu + i * bdu + threadID)
          ! Reset in case of bwd solve
          upd(row_upd + i + threadID) = zero
        ! upd(blk_id * ldu + i + threadID) = zero
        end do
      end do
#if defined(SPLLT_PROFILING_FLOP)
    lflops = blkn * nrhs
#endif
    else
      !Compute first index in rhs and upd to proceed
 !!  print *, "blkn : ", blkn
      do r = 0, nrhs - 1
 !!     print *, "RHS ", r
        !Compute the current row in the local memory spaces
        row_rhs = rhsPtr(l) * nrhs + r * blkn
        row_upd = rhsPtr(l) * nrhs * tdu + r * blkn * tdu
 !!     print *, "Row Pointer in RHS : ", row_rhs, ", and in UPD", row_upd
        !Reduce from nthread values to 1
        do i = 1, blkn
          temp = zero
          do j = 1, nthread
 !!         print *, "Reduce over threads", row_upd + (i - 1) * tdu + j
            temp = temp + upd(row_upd + (i - 1) * tdu + j)
          end do
 !!       print *, "To update row ", i, "in rhs(", row_rhs + i, ') = ', &
 !!         rhs(row_rhs + i)
          rhs(row_rhs + i) = rhs(row_rhs + i) + temp !&
        end do
      end do
      ! Reset in case of bwd solve
 !!   print *, "Reset from", rhsPtr(l) * nrhs * tdu + 1, 'to', &
 !!     row_upd + tdu * blkn
      upd(rhsPtr(l) * nrhs  * tdu + 1 : &
        row_upd + tdu * blkn) = zero !TODO fix it
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
 !! print *, "Call solve on col", col, "located in pos", &
 !!   1 + rhsPtr(l) * nrhs
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
      call slv_fwd_update_ileave(blkm, blkn, col, offset, index, &
        lcol(sa : sa + blkn * blkm - 1), blkn, nrhs,                      &
        upd, tdu, rhs, n, bdr, rhsPtr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (blkn * nrhs * blkm)
#endif
    endif

#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_timer(threadID, timer, flops)
#endif
  end subroutine solve_fwd_block_work_ileave



  subroutine solve_fwd_update_work_ileave(blk, rhsPtr, blkm, blkn, col,   &
      offset, index, lcol, blk_sa, nrhs, upd, ldu, bdu, tdu, rhs, n, bdr, &
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
    integer,                  intent(in)    :: ldu
    integer,                  intent(in)    :: bdu
    integer,                  intent(in)    :: tdu
    integer,                  intent(in)    :: n
    integer,                  intent(in)    :: bdr
    integer,                  intent(in)    :: threadID
    double precision,         intent(out)   :: flops
    integer,                  intent(in)    :: rhsPtr(:)
    integer,                  intent(inout) :: index(:)
    real(wp),                 intent(inout) :: lcol(:)
    real(wp),                 intent(inout) :: upd(:)
    real(wp),                 intent(inout) :: rhs(:)
    real(wp),                 intent(inout) :: xlocal(:,:)

    call slv_fwd_update_ileave(blkm, blkn, col, offset, &
      index, lcol(blk_sa : blk_sa + blkn * blkm - 1), blkn, nrhs,    &
      upd, tdu, rhs, n, bdr, rhsPtr, xlocal(:, threadID + 1))

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (blkm * nrhs * blkn)
#endif

  end subroutine solve_fwd_update_work_ileave



  subroutine solve_fwd_node_ileave(nrhs, rhs, n, ldr, bdr, fkeep, node, &
      xlocal, rhs_local, ldu, bdu, tdu, task_manager, trace_id)
    use spllt_data_mod
    use trace_mod
    use utils_mod
    use timer_mod
    use task_manager_mod
    implicit none

    type(spllt_fkeep), target,  intent(in)    :: fkeep
    integer,                    intent(in)    :: nrhs ! Number of RHS
    integer,                    intent(in)    :: n
    integer,                    intent(in)    :: ldr
    integer,                    intent(in)    :: bdr
    integer,                    intent(in)    :: node 
    integer,                    intent(in)    :: ldu
    integer,                    intent(in)    :: bdu
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
        ldu, fkeep%bc(dblk)%blkn, tdu, rhs, n, ldr, fkeep%bc(dblk)%blkn,&
        xlocal, fkeep, fwd_block_id)
      call spllt_tac(1, task_manager%workerID, timer)

      do ii = jj+1, nr

        blk = dblk+ii-jj

        !
        ! Forward update with off-diagonal
        !
 !!     print *, "Submit fwd update ", blk
        call spllt_tic("submit fwd update", 2, task_manager%workerID, timer)
        call task_manager%solve_fwd_update_il_task(blk, node, nrhs, rhs_local,&
          ldu, bdu, tdu, rhs, n, ldr, bdr, xlocal, fkeep, fwd_update_id)
        call spllt_tac(2, task_manager%workerID, timer)

      end do
      
      ! Update diag block in node          
      dblk = fkeep%bc(dblk)%last_blk + 1
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_fwd_node_ileave


end module spllt_solve_kernels_mod
