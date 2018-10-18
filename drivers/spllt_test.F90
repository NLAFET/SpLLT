!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
!> \author    Sebastien Cayrols
program spllt_test
  use spllt_data_mod
  implicit none

  type(spllt_options) :: cntl

  write(*,'("SpLLT test")')

  call spllt_test_mat("matrix.rb", cntl)

  stop
  
contains

  subroutine spllt_test_mat(mf, cntl)
    use spllt_mod
    use spral_rutherford_boeing
    use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_PSDEF
    ! use spllt_stf_mod 
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
    use parsec_f08_interfaces
    use spllt_parsec_mod
    use spllt_parsec_factorization_mod
    use spllt_ptg_mod
#endif
    use spllt_solve_mod
    implicit none
    
    character(len=*), intent(in) :: mf ! matrix file
    type(spllt_options) :: cntl

    ! matrix reader options (Rutherford Boeing)
    type(rb_read_options) :: rb_options
    integer :: rb_flag

    ! Matrix description (Rutherford Boeing)
    character(len=200) :: matfile = ''
    integer :: m, n
    integer, dimension(:), allocatable :: ptr, row
    real(wp), dimension(:), allocatable :: val

    ! matrix reader options (Matrix Market)
    integer :: mm_flag
    integer :: nnz
    ! Matrix description (Rutherford Boeing)
    integer, dimension(:), allocatable :: indx, jndx
    real(wp), dimension(:), allocatable :: val_in


    ! right-hand side and solution
    integer :: nrhs
    double precision, dimension(:, :), allocatable :: rhs, soln 
    double precision, dimension(:), allocatable :: res

    ! indexes
    integer :: r, i, j
    integer :: k

    ! type(spllt_keep) :: keep
    type(spllt_inform) :: info

    ! SpLLT
    type(spllt_akeep) :: akeep
    type(spllt_fkeep), target  :: fkeep
    integer, dimension(:), allocatable :: order

    ! timing
    integer :: start_t, stop_t, rate_t
    ! flags
    integer :: flag, more
       
    ! test options
    type(spllt_options) :: options 

    ! stats
    real :: smanal, smfact, smaflop, smafact

#if defined(SPLLT_USE_PARSEC) && defined(SPLLT_USE_MPI)
    type(parsec_handle_t) :: gat_hdl
    type(c_ptr) :: base_desc
    integer(c_int) :: nbc
    type(c_ptr) :: bc_c
#endif

    ! Set nrhs
    nrhs = 1

   call spllt_parse_args(options, matfile, nrhs)

   cntl%nb = options%nb
   cntl%ncpu = options%ncpu

   if (matfile .eq. '') matfile = 'matrix.rb'

    if (options%nemin .gt. 0) then
       cntl%nemin = options%nemin
    end if

    if (options%prune_tree) cntl%prune_tree = .true. 

    write(*,*) '  mat: ', matfile
    write(*,*) '   nb: ', options%nb    
    write(*,*) ' ncpu: ', options%ncpu
    write(*,*) 'nemin: ', cntl%nemin
    
    ! Read in a matrix
    write(*, "(a)") "Reading..."
    if (options%fmt .eq. 'csc') then
       ! Rutherford boeing format

       ! DEBUG ensure matrix is diag dominant
       rb_options%values = 3 ! Force diagonal dominance
       call rb_read(matfile, m, n, ptr, row, val, rb_options, rb_flag)
       if(rb_flag.ne.0) then
          print *, "Rutherford-Boeing read failed with error ", rb_flag
          stop
       endif

    else if (options%fmt .eq. 'coo') then
       ! Matrix Market format

       ! read matrix
       call mm_read(matfile, m, n, nnz, indx, jndx, val_in, mm_flag)
       if(mm_flag.ne.0) then
          print *, "Matrix Market read failed with error ", mm_flag
          stop
       endif
       
       ! convert to csc format
       call coo_to_csc_double(m, n, nnz, indx, jndx, val_in, & 
            row, val, ptr, mm_flag)
       if(mm_flag.ne.0) then
          print *, "COO to CSC convertion failed with error ", mm_flag
          stop
       endif

       deallocate(indx, jndx, val_in)

       ! print *, "m,n,nnz: ", m,n,nnz
       ! print *, "size ptr: ", size(ptr)
       ! print *, "size ptr: ", size(val)
       ! print *, "mm_flag: ", mm_flag
    else
       write(*, "(a)") "Matrix format not suported"
       stop
    end if
    write(*, "(a)") "ok"

    print *, "n = ", n

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

    ! Analyse SpLLT
    write(*, "(a)") "Analyse..."
    call system_clock(start_t, rate_t)
    call spllt_analyse(akeep, fkeep, cntl, n, ptr, row,  info, order)
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "error detected during analysis"
       stop
    endif
    call system_clock(stop_t)
    write(*, "(a)") "ok"
    print *, "Analyse took ", (stop_t - start_t)/real(rate_t)
    print "(a,es10.2)", "Predict nfact = ", real(info%ssids_inform%num_factor)
    print "(a,es10.2)", "Predict nflop = ", real(info%ssids_inform%num_flops)
    smaflop = real(info%ssids_inform%num_flops)
    smafact = real(info%ssids_inform%num_factor)

    ! Print elimination tree
    call spllt_print_atree(akeep, fkeep, cntl)

#if defined(SPLLT_USE_OMP)
    !$omp parallel num_threads(cntl%ncpu)
    !$omp master
#endif
    ! Initialize solver
    ! goto 9999 ! DEBUG: jump init, factor, solve and finalize
    call spllt_init(cntl)

    ! Factor phase
    write(*, "(a)") "Factor..."
    write(*,'("   nb: ", i6)') cntl%nb
    write(*,'("# cpu: ", i6)') cntl%ncpu
    call system_clock(start_t, rate_t)
    call spllt_factor(akeep, fkeep, cntl, val, info)
    call spllt_wait() ! Wait for factorization to finish.
    if(info%flag .lt. 0) then
       write(*, "(a)") "failed factorization"
    endif
    call system_clock(stop_t)
    write(*, "(a)") "ok"
    print *, "Factor took ", (stop_t - start_t)/real(rate_t)

#if defined(SPLLT_USE_PARSEC) && defined(SPLLT_USE_MPI)

    write(*,*) "[spllt_test] gather blocks for solve"
    base_desc = alloc_desc()

    nbc = size(fkeep%bc,1)
    bc_c = c_loc(fkeep%bc(1))

    call data_init(base_desc, bc_c, nbc, nds, rank)

    gat_hdl = gather(fkeep%ddesc, base_desc, size(fkeep%bc,1), fkeep%maxmn)

    call parsec_enqueue(ctx, gat_hdl)
    call parsec_context_start(ctx)

    call parsec_context_wait(ctx)

#endif

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
    ! call MA87_solve(nrhs, n, soln, order, keep, control, info)
    call spllt_solve(fkeep, cntl, order, soln(:,1), info)
    ! if(info%flag .lt. spllt_success) then
    !    write(*, "(a,i4)") " fail on 1d solve", &
    !         info%flag
    ! endif

    ! call spllt_bwerr(a, x, b, resid)

    print *, "number bad cmp = ", count(abs(soln(1:n,1)-1.0).ge.1e-6)
    print *, "fwd error || ||_inf = ", maxval(abs(soln(1:n,1)-1.0))
    allocate(res(nrhs))
    call internal_calc_norm(n, ptr, row, val, soln, rhs, nrhs, res)
    print *, "bwd error scaled = ", res

9998 continue

    ! call MA87_finalise(keep, control)
    
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
    ! print *, "[internal_calc_norm] normA = ", normA

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
