!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
!> \author    Sebastien Cayrols
module spllt_solve_kernels_mod
  
contains


  subroutine slv_solve(n, nelim, dest, trans, unit, nrhs, rhs, ldr)
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

  end subroutine slv_solve



  subroutine slv_fwd_update(m, nelim, dest, ldd, &
      n, nrhs, rhs, ldr, xlocal, ldx, reset)
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
    logical,  intent(in)    :: reset          ! When true, the DGEMM overwrites
                                              ! xlocal

    integer   :: i
    integer   :: j
    integer   :: k
    real(wp)  :: w ! temporary work value
    real(wp)  :: alpha, beta
    integer   :: threadID
#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    type(spllt_timer_t), save :: timer
#endif

    threadID = 0
 !$ threadID = omp_get_thread_num()


#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    call spllt_open_timer(threadID, "slv_fwd_update", timer)
#endif

    if(nelim.eq.0) return

    alpha = -one
    beta  = one
    if(reset) beta  = zero

    ! Fix NaN problem : it may have NaN in 
    ! the first use of the local x. To fix it, we have to reset it manually.
    ! In case, BLAS can handle it, this reset can be removed
    if(reset) xlocal(1 : m) = zero 

    ! forward substitution
    if(nrhs.eq.1) then
      !print *, "local RHS", rhs(1:nelim) 
      !print *, "Result before", xlocal(1 : m)
       if(m-nelim.gt.10 .and. nelim.gt.4) then
!!! Single rhs, BLAS 2

          call dgemv('T', nelim, m, alpha, dest, ldd, rhs, 1, beta, xlocal, 1)

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
           xlocal(i) = beta * xlocal(i) + w
         end do
       endif
      !print *, "Result DGEMM", xlocal(1 : m)
    else
!!! Multiple rhs, BLAS 3
      call dgemm('T', 'N', m, nrhs, nelim, alpha, dest, ldd, rhs, ldr, &
           beta, xlocal, ldx)
    endif

#if defined(SPLLT_SOLVE_KERNEL_SCATTER)
    call spllt_close_ftimer(threadID, timer)
#endif

  end subroutine slv_fwd_update


  subroutine slv_bwd_update(m, nelim, dest, ldd, &
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
    call spllt_open_timer(threadID, "slv_bwd_update", timer)
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

  end subroutine slv_bwd_update


#if 0
  subroutine solve_fwd_block_work(blkm, blkn, lcol, sa, nrhs, y, &
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

    call spllt_open_timer(threadID, "solve_fwd_block_work", timer)
#endif

    ! Perform triangular solve
    call slv_solve(blkn, blkn, lcol(sa : sa + blkn * blkn - 1), &
      'Transpose    ', 'Non-unit', nrhs, y, ldy)

#if defined(SPLLT_PROFILING_FLOP)
    flops = blkn * blkn * nrhs
#endif

    ! Deal with any left over trapezoidal part of diagonal block
    blkm = blkm - blkn
    if(blkm .gt. 0) then
      sa = sa + blkn * blkn
      call slv_fwd_update(blkm, blkn, lcol(sa : sa + blkn * blkm - 1), &
        blkn, blkn, nrhs, y, ldy, xlocal, ldx)

#if defined(SPLLT_PROFILING_FLOP)
      flops = flops + 2 * (blkn * nrhs * blkm)
#endif
    endif

#if defined(SPLLT_TIMER_TASKS) || defined(SPLLT_SOLVE_KERNEL_GATHER)
    call spllt_close_timer(threadID, timer, flops)
#endif
  end subroutine solve_fwd_block_work



  subroutine solve_fwd_update_work(blkm, blkn, lcol, blk_sa, &
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

    call slv_fwd_update(blkm, blkn,                              &
      lcol(blk_sa : blk_sa + blkn * blkm - 1), blkn, n, nrhs, rhs, ldr, &
      xlocal, ldx)

#if defined(SPLLT_PROFILING_FLOP)
    flops = 2 * (blkm * nrhs * blkn)
#endif

  end subroutine solve_fwd_update_work
#endif

  subroutine solve_fwd_node(nrhs, rhs, n, fkeep, node, &
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
    s_nb   = fkeep%nodes(node)%snb
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
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("submit fwd block", 1, task_manager%workerID, timer)
#endif
      call task_manager%solve_fwd_block_task(dblk, nrhs, n, rhs, &
        fkeep, trace_id)
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(1, task_manager%workerID, timer)
#endif

      do blk = dblk + 1, fkeep%sbc(dblk)%last_blk

       !blk = dblk+ii-jj

        !
        ! Forward update with off-diagonal
        !
 !!     print *, "Submit fwd update ", blk
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
        call spllt_tic("submit fwd update", 2, task_manager%workerID, timer)
#endif
        call task_manager%solve_fwd_update_task(blk, node, nrhs, &
            n, rhs, fkeep, trace_id)
       !call task_manager%solve_fwd_update_task(blk, node, nrhs, rhs_local,&
       !  tdu, rhs, n, xlocal, fkeep, fwd_update_id)
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
        call spllt_tac(2, task_manager%workerID, timer)
#endif
      end do
      
      ! Update diag block in node          
      dblk = fkeep%sbc(dblk)%last_blk + 1
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_fwd_node


  subroutine solve_bwd_node(nrhs, rhs, n, fkeep, node, &
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
    s_nb   = fkeep%nodes(node)%snb
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
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
        call spllt_tic("submit bwd update", 1, task_manager%workerID, timer)
#endif
        call task_manager%solve_bwd_update_task(blk, node, nrhs, & 
          n, rhs, fkeep, trace_id)
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
        call spllt_tac(1, task_manager%workerID, timer)
#endif
       !print *, "blk", blk, "x associated", fkeep%sbc(blk)%p_upd
      end do

      !
      ! Backward solve with block on diagoanl
      !
 !!   print *, "Submit bwd block ", dblk
     !print *, "dblk", dblk, "y INPUT", fkeep%sbc(dblk)%p_upd
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tic("submit bwd block", 2, task_manager%workerID, timer)
#endif
      call task_manager%solve_bwd_block_task(dblk, nrhs, &
        n, rhs, fkeep, trace_id)
#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
      call spllt_tac(2, task_manager%workerID, timer)
#endif
     !print *, "dblk", dblk, "x associated", fkeep%sbc(dblk)%p_upd
     
      ! Update diag block in node       
      if (jj .gt. 1) dblk = fkeep%sbc(dblk-1)%dblk
    end do

#if defined(SPLLT_TIMER_TASKS_SUBMISSION)
    call spllt_close_timer(task_manager%workerID, timer)
#endif
  end subroutine solve_bwd_node

end module spllt_solve_kernels_mod
