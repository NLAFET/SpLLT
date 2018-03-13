threadID = omp_get_thread_num()

if(ndep_lvl .le. chunk) then
  call trace_event_start(trace_id, threadID)

! print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
! print *, p_dep

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
end if
