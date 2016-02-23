module spllt_kernels_mod
  implicit none

contains

  !*************************************************  
  !
  ! TASK_FACTORIZE_BLOCK (uses Lapack routine dpotrf and dtrsm)
  ! A_ii <- L_ii
  !
  subroutine spllt_factor_diag_block(m, n, dest)
    use spllt_mod
    implicit none
    
    integer, intent(in) :: m ! number of rows in dest
    integer, intent(in) :: n ! number of columns in dest
    real(wp), dimension(*), intent(inout) :: dest ! holds block
    ! on diagonal of factor L. It may not be square

    integer :: dpotrf_info ! error flag for dpotrf
    integer :: i, j ! Loop indices

    call dpotrf('Upper', n, dest, n, dpotrf_info)
    ! check for errors
    if(dpotrf_info.ne.0) return

    ! Do dtrsm with any remainder below diagonal block
    if(m.gt.n) then
       call dtrsm('Left', 'Upper', 'Transpose', 'Non-Unit', n, &
            m-n, one, dest, n, dest(1+n*n), n)
    endif

  end subroutine spllt_factor_diag_block  

end module spllt_kernels_mod
