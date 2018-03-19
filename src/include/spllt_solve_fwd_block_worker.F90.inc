threadID  = omp_get_thread_num() + 1

if(ndep_lvl .le. chunk) then
#if defined(SPLLT_OMP_TRACE)
  call trace_event_start(trace_id, threadID - 1)
#endif

! print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
! print *, p_dep

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
#if defined(SPLLT_PROFILING_FLOP)
  scheduler%task_info(threadID)%nflop = scheduler%task_info(threadID)%nflop + &
    (n * n * nrhs)
#endif
  offset = offset + n

  ! Deal with any left over trapezoidal part of diagonal block
  m = m - n
  if(m .gt. 0) then
    sa = sa + n * n
    call slv_fwd_update(m, n, col, offset, p_index,             &
      p_lcol(sa : sa + n * m - 1), n, nrhs,                     &
      p_upd(:, threadID), ldr, p_rhs,                           &
      ldr, p_xlocal(:, threadID))
#if defined(SPLLT_PROFILING_FLOP)
    scheduler%task_info(threadID)%nflop = scheduler%task_info(threadID)%nflop &
      + 2 * (n * nrhs * m)
#endif
  endif

#if defined(SPLLT_OMP_TRACE)
  call trace_event_stop (trace_id, threadID - 1)
#endif
end if
