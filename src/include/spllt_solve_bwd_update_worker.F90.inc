threadID  = omp_get_thread_num()

if(ndep_lvl .le. chunk) then
  call trace_event_start(trace_id, threadID)

! print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
! print *, p_dep
  call  slv_bwd_update(m, n, col, offset, p_index,      &
        p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,   &
        p_rhs, p_upd(:, threadID + 1), ldr,             &
        xlocal(:, omp_get_thread_num() + 1))

  call trace_event_stop (trace_id, omp_get_thread_num())
end if
