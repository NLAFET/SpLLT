 call task_manager%refresh_worker()
 threadID = task_manager%workerID
!threadID = omp_get_thread_num()
#if defined(SPLLT_TIMER_TASKS)
call spllt_tic("bwd block task", 13, threadID, timer)
#endif
call task_manager%incr_nrun()

if(ndep_lvl .le. chunk) then
! call solve_bwd_block_work(m, n, col, offset, p_index, p_lcol, sa, nrhs, &
!     p_upd, threadID, ldr, p_xlocal, trace_id, task_manager)
#if defined(SPLLT_OMP_TRACE)
   call trace_event_start(trace_id, threadID)
#endif
 
#if defined(SPLLT_VERBOSE)
   print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
   print *, p_dep
#endif
 
   ! Perform retangular update from diagonal block
   if(m .gt. n) then
     call slv_bwd_update(m - n, n, col, offset + n, p_index,     &
       p_lcol(sa + n * n : sa + n * m - 1), n, nrhs, p_rhs,      &
       p_upd(:, threadID + 1), ldr, p_xlocal(:, threadID + 1))
 
#if defined(SPLLT_PROFILING_FLOP)
     call task_manager%nflop_performed(2 * (n * nrhs * (m-n)) + zero)
#endif
   endif
 
   ! Sum contributions to rhs
   do r = 0, nrhs-1
     do j = 1, nthread
       do i = col + r*ldr, col+n-1 + r*ldr
         p_rhs(i)    = p_rhs(i) + p_upd(i, j)
         p_upd(i,j)  = zero ! Reset in case of next fwd solve
       end do
     end do
   end do
 
   ! Perform triangular solve
   call slv_solve(n, n, col, p_lcol(sa : sa + n * n - 1), &
        'Non-Transpose', 'Non-unit', nrhs, p_rhs, ldr)
 
#if defined(SPLLT_PROFILING_FLOP)
   call task_manager%nflop_performed(n * n * nrhs + zero)
#endif
 
#if defined(SPLLT_OMP_TRACE)
   call trace_event_stop(trace_id, threadID)
#endif
end if

#if defined(SPLLT_TIMER_TASKS)
call spllt_tac(13, threadID, timer)
#endif
