threadID  = omp_get_thread_num()

if(ndep_lvl .le. chunk) then
! print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : ]"
  call slv_fwd_update(m, n, col, offset, p_index,         &
    p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
    p_upd(:, threadID + 1), ldr, p_rhs,                   &
    ldr, p_xlocal(:, threadID + 1))
end if
