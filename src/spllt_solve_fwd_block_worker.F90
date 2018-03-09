threadID  = omp_get_thread_num()

if(ndep_lvl .le. chunk) then
! if(.true.) then
! if(.false.)then
!   print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
!   print *, p_dep
!    call spllt_solve_fwd_block_worker(m, n, nrhs, nthread, col, ldr, &
!      sa, offset, threadID, p_upd, p_rhs, p_lcol, p_index, p_xlocal)
! else
!   print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
!   print *, p_dep
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
      call slv_fwd_update(m, n, col, offset, p_index,             &
        p_lcol(sa : sa + n * m - 1), n, nrhs,                     &
        p_upd(:, threadID + 1), ldr, p_rhs,                       &
        ldr, p_xlocal(:, threadID + 1))
    endif
! end if
end if
