program spllt_test
  use spllt_data_mod
  implicit none

  type(spllt_cntl) :: cntl

  write(*,'("SpLLT test")')

  call spllt_test_mat("matrix.rb", cntl)

  stop
  
contains

  subroutine spllt_test_mat(mf, cntl)
    use spllt_mod
    use spral_rutherford_boeing
    use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_PSDEF
    use spral_ssids
    use hsl_ma87_double
    use spllt_analyse_mod, only : spllt_analyse
    use spllt_stf_factorization_mod 
#if defined(SPLLT_USE_STARPU)
    use iso_c_binding
    use starpu_f_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_factorization_task_mod, only: spllt_data_unregister
#endif
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE) 
    use trace_mod
#endif
#elif defined(SPLLT_USE_PARSEC)
    use dague_f08_interfaces
    use spllt_parsec_mod
    use spllt_ptg_mod
#endif
    implicit none
    
    character(len=*), intent(in) :: mf ! matrix file
    type(spllt_cntl) :: cntl

    ! matrix reader options
    type(rb_read_options) :: rb_options
    integer :: rb_flag
    
    ! Matrix description
    character(len=200) :: matfile = ''
    integer :: m, n
    integer, dimension(:), allocatable :: ptr, row
    real(wp), dimension(:), allocatable :: val

    ! right-hand side and solution
    integer :: nrhs
    double precision, dimension(:, :), allocatable :: rhs, soln 
    double precision, dimension(:), allocatable :: res

    ! indexes
    integer :: r, i, j
    integer :: k

    ! ma87
    type(ma87_keep) :: keep
    type(ma87_control) :: control
    type(ma87_info) :: info

    ! SpLLT
    type(spllt_adata_type) :: adata
    type(spllt_data_type)  :: fdata
    integer, dimension(:), allocatable :: order

    ! timing
    integer :: start_t, stop_t, rate_t
    ! flags
    integer :: flag, more
    
    ! ssids options 
    type(ssids_options) :: ssids_opt
 
    ! ssids structures
    type(ssids_inform) :: inform ! stats
    type(ssids_akeep) :: akeep   ! analysis data
   
    ! test options
    type(spllt_options) :: options 

    ! stats
    real :: smanal, smfact, smaflop, smafact

    ! Set nrhs
    nrhs = 1
    
    call splllt_parse_args(options)

    cntl%nb   = options%nb
    cntl%ncpu = options%ncpu

    ! Set matrix file
    if (options%mat.ne.'') then
       matfile = options%mat
    else
       matfile = mf
    end if

    ! Set nemin
    if (options%nemin .gt. 0) then
       cntl%nemin = options%nemin
    end if

    ! Tree pruning
    if (options%prune_tree) cntl%prune_tree = .true.
    
    ! Read in a matrix
    write(*, "(a)") "Reading..."
    ! DEBUG ensure matrix is diag dominant
    rb_options%values = 3 ! Force diagonal dominance
    call rb_read(matfile, m, n, ptr, row, val, rb_options, rb_flag)
    if(rb_flag.ne.0) then
       print *, "Rutherford-Boeing read failed with error ", rb_flag
       stop
    endif
    write(*, "(a)") "ok"

    ! Make up a rhs associated with the solution x = 1.0
    allocate(rhs(n, nrhs), soln(n, nrhs))
    rhs = 0
    do r = 1, nrhs
       do i = 1, n
          do j = ptr(i), ptr(i+1)-1
             k = row(j)
             rhs(k, r) = rhs(k, r) + val(j)
             if(i.eq.k) cycle
             rhs(i, r) = rhs(i, r) + val(j)
          end do
       end do
    end do

    ! check matrix format is correct
    call cscl_verify(6, SPRAL_MATRIX_REAL_SYM_PSDEF, n, n, &
         ptr, row, flag, more)
    if(flag.ne.0) then
       print *, "CSCL_VERIFY failed: ", flag, more
       stop
    endif

    ! ordering
    allocate(order(n))

    ! Set options for analysis
    ssids_opt%ordering = 1 ! use Metis ordering
    ssids_opt%scaling = 0 ! no scaling

    ! Analyse SSIDS
    write(*, "(a)") "Analyse..."
    call system_clock(start_t, rate_t)
    call ssids_analyse(.false., n, ptr, row, akeep, &
         ssids_opt, inform, order, val=val)
    call system_clock(stop_t)
    print *, "Used order ", ssids_opt%ordering
    if (inform%flag < 0) then
       print *, "oops on analyse ", inform%flag
       stop
    endif
    write(*, "(a)") "ok"
    print *, "Analyse took ", (stop_t - start_t)/real(rate_t)
    !print *, "Used maximum memory of ", inform%maxmem
    smanal = (stop_t - start_t)/real(rate_t)
    print "(a,es10.2)", "Predict nfact = ", real(inform%num_factor)
    print "(a,es10.2)", "Predict nflop = ", real(inform%num_flops)
    print "(a6, i10)", "nparts", inform%nparts
    print "(a6, es10.2)", "cpu_flops", real(inform%cpu_flops)
    ! print "(a6, es10.2)", "gpu_flops", real(inform%gpu_flops)
    smaflop = real(inform%num_flops)
    smafact = real(inform%num_factor)

    ! Analyse SpLLT
    call spllt_analyse(adata, fdata, n, ptr, row, order, akeep, keep, cntl, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "error detected during analysis"
       stop
    endif

    ! Print elimination tree
    call spllt_print_atree(adata, keep, cntl)

#if defined(SPLLT_USE_OMP)
    !$omp parallel num_threads(cntl%ncpu)
    !$omp master
#endif
    ! Initialize solver
    ! goto 9999 ! DEBUG: jump init, factor, solve and finalize
    call spllt_init(cntl)

    ! factorize matrix
    ! goto 9998 ! DEBUG: jump factor and solve
    ! goto 9999 ! DEBUG: jump factor, solve and finalize
    write(*,'("Factorize...")')
    write(*,'("   nb: ", i6)') cntl%nb
    write(*,'("# cpu: ", i6)') cntl%ncpu
    call system_clock(start_t, rate_t)
    ! TODO create factorize method
#if defined(SPLLT_USE_STF) || defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)

#if defined(SPLLT_STF_LL)
    ! Unroll the DAG using a Left-Looking strategy
    call spllt_stf_ll_factorize(n, ptr, row, val, order, keep, control, info, fdata, cntl)

#else
    call spllt_stf_factorize(n, ptr, row, val, order, keep, control, info, adata, fdata, cntl)
    ! call MA87_factor(a%n, a%ptr, a%row, a%val, order, keep, control, info)
#endif

#elif defined(SPLLT_USE_PARSEC)
    call spllt_ptg_factorize(adata, val, keep, cntl, fdata, info)
#endif

#if defined(SPLLT_USE_STARPU)
    ! wait for task completion
    call starpu_f_task_wait_for_all()
#if defined(SPLLT_USE_GPU)
    ! put data back on home node e.g. from GPU to CPU
    call spllt_data_unregister(keep, fdata)
#endif
    call starpu_f_fxt_stop_profiling()
#elif defined(SPLLT_USE_OMP)
    !$omp taskwait
#elif (SPLLT_USE_PARSEC)
    ! write(*,'("[>] Parsec wait rank: ", i6)') rank
    call dague_context_wait(ctx)
#endif

    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "failed factorization"
    endif
    call system_clock(stop_t)
    write(*, "(a)") "ok"
    print *, "Factor took ", (stop_t - start_t)/real(rate_t)

    call spllt_finalize()
    ! #if defined(SPLLT_USE_STARPU)
    !     call system_clock(start_starpushutdown_t, rate_starpushutdown_t)
    !     call starpu_f_shutdown()
    !     call system_clock(stop_starpushutdown_t)
    !     write(*,'("[>] [spllt_test_mat] StarPU shutdown time: ", es10.3, " s")') &
    !          &(stop_starpushutdown_t - start_starpushutdown_t)/real(rate_starpushutdown_t)
#if defined(SPLLT_USE_OMP)
    !$omp end master
    !$omp end parallel
    ! #elif defined(SPLLT_USE_PARSEC)
    !     ! call dague_fini(ctx)
    !     call parsec_fini(ctx)
#endif
    
    write(*,'("Solve...")')
    soln = rhs ! init solution with RHS
    ! solve
    call MA87_solve(nrhs, n, soln, order, keep, control, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a,i4)") " fail on 1d solve", &
            info%flag
    endif

    ! call spllt_bwerr(a, x, b, resid)

    print *, "number bad cmp = ", count(abs(soln(1:n,1)-1.0).ge.1e-6)
    print *, "fwd error || ||_inf = ", maxval(abs(soln(1:n,1)-1.0))
    allocate(res(nrhs))
    call internal_calc_norm(n, ptr, row, val, soln, rhs, nrhs, res)
    print *, "bwd error scaled = ", res

    ! write(*,'("[>] [solve] bwderr ||Ax-b|| / (||A||||x|| + ||b||): ", es10.3)') resid    

9998 continue

    call MA87_finalise(keep, control)
    
  end subroutine spllt_test_mat

  subroutine internal_calc_norm(n, ptr, row, val, x_vec, b_vec, nrhs, res)
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in) :: nrhs
    real(wp), dimension(nrhs*n), intent(in) :: x_vec
    real(wp), dimension(nrhs*n), intent(in) :: b_vec
    real(wp), dimension(nrhs), intent(out) :: res

    integer :: i, j, k, r
    double precision, allocatable, dimension(:) :: x_norm
    real(wp), dimension(:), allocatable :: res_vec
    double precision :: temp
    double precision :: normA

    ! Find the residual
    allocate(res_vec(n*nrhs), x_norm(nrhs))
    res_vec = 0
    do i = 1, n
       do j = ptr(i), ptr(i+1)-1
          r = row(j)
          do k = 0, nrhs-1
             res_vec(i+k*n) = res_vec(i+k*n) + &
                  val(j) * x_vec(r+k*n)
          end do
          if(r.eq.i) cycle
          do k = 0, nrhs-1
             res_vec(r+k*n) = res_vec(r+k*n) + &
                  val(j) * x_vec(i+k*n)
          end do
       end do
    end do
    res_vec(:) = res_vec(:) - b_vec(:)

    ! Find matrix norm
    call matrix_inf_norm(n, ptr, row, val, normA)

    ! Find x norm
    do i = 1, nrhs
       x_norm(i) = 0
       do j =1, n
          x_norm(i) = max(x_norm(i), abs(x_vec((i-1)*n+j)))
          if(x_vec((i-1)*n+j).ne.x_vec((i-1)*n+j)) then ! Tests for NaN
             x_norm(i) = x_vec((i-1)*n+j)
             exit
          endif
       end do
    end do

    ! Scaled residual = ||r|| / ( ||A|| ||x|| + ||b|| )
    do i = 1, nrhs
       temp = normA * x_norm(i) + &
            maxval(abs(b_vec((i-1)*n+1:i*n)))
       if(temp .eq. 0) then
          res(i) = maxval(abs(res_vec((i-1)*n+1:i*n)))
       else
          res(i) = maxval(abs(res_vec((i-1)*n+1:i*n))) / temp
       endif
    end do
  end subroutine internal_calc_norm

  subroutine matrix_inf_norm(n, ptr, row, val, norm)
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    real(wp), intent(out) :: norm

    real(wp), allocatable, dimension(:) :: row_norm
    integer :: i

    allocate(row_norm(n))

    row_norm = 0
    do i = 1, ptr(n+1)-1
       row_norm(row(i)) = row_norm(row(i)) + abs(val(i))
    end do

    norm = maxval(row_norm) 
  end subroutine matrix_inf_norm

end program spllt_test
