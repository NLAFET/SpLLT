program spllt_test
  use hsl_ma87_double
  use spllt_mod
  implicit none

  type(spllt_cntl) :: cntl
  integer n, nnz
  
  write(*,'("[spllt test]")')
  
  ! n   = 384
  ! nnz = n*n

  ! n   = 512
  ! nnz = 10000

  ! call spllt_test_rand(n, nnz, cntl)

  call spllt_test_mat("matrix.rb", cntl)
  
  stop

contains

  subroutine spllt_test_mat(mf, cntl)
    use spllt_mod
    use spllt_stf_factorization_mod
    use spral_rutherford_boeing
    use hsl_mc68_integer
    use spllt_analyse_mod
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
    use hsl_mc56_double
    use hsl_mc69_double
    implicit none
    
    character(len=*), intent(in) :: mf
    type(spllt_cntl) :: cntl

    type(rb_reader_options) :: rb_options
    integer :: flag
    integer :: m, n

    type(spllt_options) :: options
    character(len=200) :: matfile    
    type(zd11_type) :: a
    integer :: nrhs
    integer, dimension(:), allocatable :: order
    real(wp) :: num_flops, num_factor, resid
    real(wp), dimension(:), allocatable :: b, x
    type(spllt_adata_type) :: adata
    type(spllt_data_type)  :: fdata

    ! mc68
    type(mc68_control) :: order_control
    type(mc68_info) :: order_info

    ! ma87
    type(ma87_keep) :: keep
    type(ma87_control) :: control
    type(ma87_info) :: info

    ! mc56
    type(mc56_control) :: control56
    integer :: info56

    ! timing
    integer :: start_t, stop_t, rate_t

    ! StarPU
#if defined(SPLLT_USE_STARPU) 
    ! when using StarPU
    ! integer(c_int) :: ret
    ! timing
    ! integer :: start_starpuinit_t, stop_starpuinit_t, rate_starpuinit_t
    ! integer :: start_starpushutdown_t, stop_starpushutdown_t, rate_starpushutdown_t
#endif

    call splllt_parse_args(options)

    cntl%nb   = options%nb
    cntl%ncpu = options%ncpu
    
    ! set matrix file
    if (options%mat.ne.'') then
       matfile = options%mat
    else
       matfile = mf
    end if

    ! set nemin
    if (options%nemin .gt. 0) then
       cntl%nemin = options%nemin
    end if

    ! tree pruning
    if (options%prune_tree) cntl%prune_tree = .true.

    ! write(*,*)"option mat: ", options%mat
    write(*,'("[spllt test mat] read matrix")')

    ! read matrix
    ! rb_options%values = 3 ! Make up values if necessary
    ! ! rb_options%values = 0 ! As per file
    ! rb_options%format = 1 ! Coordinate
    ! call rb_read(matfile, a%m, a%n, a%ptr, a%row, a%col, &
    !      a%val, rb_options, flag)

    control56%values = 3 ! make up values if necessary (posdef)
    call mc56_read(matfile, a, control56, info56)

    m = a%m
    n = a%n

    call mc69_cscl_clean(HSL_MATRIX_REAL_SYM_PSDEF, a%n, a%n, &
         a%ptr, a%row, flag, val=a%val)

    write(*,'("[>] generate ordering")')

    ! ordering
    allocate(order(a%n))

    ! amd ordering
    ! call amd_order(a, order)
    call mc68_order(3, a%n, a%ptr, a%row, order, &
         order_control, order_info)

    ! natural order
    ! do i = 1,n
    !    order(i) = i
    ! end do

    write(*,'("[>] perform analysis")')
    control%nb = cntl%nb
    control%nemin = cntl%nemin
    
    ! analysis
    call spllt_analyse(adata, fdata, a%n, a%ptr, a%row, order, keep, cntl, info)
    ! call MA87_analyse(a%n, a%ptr, a%row, order, keep, control, info)
    num_flops   = adata%num_flops
    num_factor  = real(adata%num_factor)
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "error detected during analysis"
       goto 9999
    endif

    write(*,'("[>] [analysis] matrix n  : ", es20.5)') real(a%n)    
    write(*,'("[>] [analysis] matrix nz : ", i40)') size(a%val)
    write(*,'("[>] [analysis] num factor: ", es20.5)') num_factor    
    write(*,'("[>] [analysis] num flops : ", es20.5)') num_flops    
    write(*,'("[>] [analysis] num nodes : ", i10)') adata%nnodes    

    call spllt_print_atree(adata, keep, cntl)

    ! generate rhs
    nrhs = 1
    allocate(b(m), x(n))
    call random_number(b)

#if defined(SPLLT_USE_OMP)
!$omp parallel num_threads(cntl%ncpu)
!$omp master
#endif
    ! initialize solver
    call spllt_init(cntl)

! #elif defined(SPLLT_USE_PARSEC)

!     ctx = parsec_init(cntl%ncpu, nds, rank)
!     write(*,'("[>] Parsec init    nodes: ", i6, ", rank: ", i6)') nds, rank
!     ! call dague_init(cntl%ncpu, ctx)
! #endif

    write(*,'("[>] factorize")')
    write(*,'("[>] [factorize]    nb: ", i6)') cntl%nb
    write(*,'("[>] [factorize] # cpu: ", i6)') cntl%ncpu
    ! factorize matrix
    call system_clock(start_t, rate_t)
    ! TODO create factorize method
#if defined(SPLLT_USE_STF) || defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)

#if defined(SPLLT_STF_LL)
    ! Unroll the DAG using a Left-Looking strategy
    call spllt_stf_ll_factorize(a%n, a%ptr, a%row, a%val, order, keep, control, info, fdata, cntl)
    
#else
    call spllt_stf_factorize(a%n, a%ptr, a%row, a%val, order, keep, control, info, adata, fdata, cntl)
    ! call MA87_factor(a%n, a%ptr, a%row, a%val, order, keep, control, info)
#endif

#elif defined(SPLLT_USE_PARSEC)
    call spllt_ptg_factorize(adata, a%val, keep, cntl, fdata, info)
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
    write(*,'("[>] [factorize] time: ", es10.3, " s")') (stop_t - start_t)/real(rate_t)

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

    write(*,'("[>] solve")')
    x = b
    ! solve
    call MA87_solve(x, order, keep, control, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a,i4)") " fail on 1d solve", &
            info%flag
    endif
    
    call spllt_bwerr(a, x, b, resid)
    
    write(*,'("[>] [solve] bwderr ||Ax-b|| / (||A||||x|| + ||b||): ", es10.3)') resid    

    call MA87_finalise(keep, control)

9999 continue

    return
  end subroutine spllt_test_mat

  subroutine spllt_test_rand(n, nnz, cntl)
    use hsl_mc68_integer
    use spllt_mod
    use spllt_analyse_mod
#if defined(SPLLT_USE_OMP)
!$ use omp_lib
#if defined(SPLLT_OMP_TRACE) 
    use trace_mod
#endif
#endif

#if defined(SPLLT_USE_STF) || defined(SPLLT_USE_OMP) || defined(SPLLT_USE_STARPU)
    use spllt_factorization_mod
#elif defined(SPLLT_USE_PARSEC)
    use dague_f08_interfaces
    use spllt_ptg_mod
#endif

    implicit none
    
    integer :: n, nnz
    type(spllt_cntl) :: cntl

    type(spllt_adata_type) :: adata
    type(spllt_data_type) :: fdata
    type(spllt_options) :: options

    ! mc68
    ! type(mc68_control) :: order_control
    ! type(mc68_info) :: order_info
    
    ! ma87
    type(ma87_keep) :: keep
    type(ma87_control) :: control
    type(ma87_info) :: info


    ! type(matrix_type) :: a
    integer :: iseed
    type(zd11_type) :: a
    real(wp), dimension(:), allocatable :: b, x
    integer :: nrhs
    integer, dimension(:), allocatable :: order
    real(wp) :: num_flops, resid

    call splllt_parse_args(options)

    cntl%nb   = options%nb
    cntl%ncpu = options%ncpu
    
    ! set nemin
    if (options%nemin .gt. 0) then
       cntl%nemin = options%nemin
    end if

    a%m = n
    a%n = n
    a%ne = nnz
    
    allocate(a%ptr(n+1))
    allocate(a%row(2*a%ne), a%col(2*a%ne), a%val(2*a%ne))

    write(*,'("[spllt test rand] generate random matrix")')

    call fa14id(iseed)

    ! generate matrix
    ! with spral
!!    call random_matrix_generate(state, SPRAL_MATRIX_REAL_SYM_PSDEF, & 
!!         & a%n, a%n, nnz, &
!!         & a%ptr, a%row, &
!!         & flag, &
!!         & val=a%val)

    call gen_random_posdef(a, a%ne, iseed)    

    write(*,'("[>] generate ordering")')

    ! ordering
    allocate(order(a%n))

    ! amd ordering
    call amd_order(a, order)
    ! call mc68_order(3, a%n, a%ptr, a%row, order, &
    !      order_control, order_info)

    ! natural order
    ! do i = 1,a%n
    !    order(i) = i
    ! end do

    write(*,'("[>] analyse")')

    ! analysis
    call spllt_analyse(adata, fdata, a%n, a%ptr, a%row, order, keep, cntl, info)
    ! call MA87_analyse(a%n, a%ptr, a%row, order, keep, control, info)
    num_flops = info%num_flops
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "error detected during analysis"
       goto 9999
    endif

    write(*,'("[>] [analysis] num flops: ", es10.3)') num_flops    
    write(*,'("[>] [analysis] num nodes: ", i10)') info%num_nodes    

    ! print atree in a dot format
    call spllt_print_atree(adata, keep, cntl)
    
    ! generate rhs
    nrhs = 1
    allocate(b(a%m), x(a%n))
    call random_number(b)

#if defined(SPLLT_USE_OMP)
!$omp parallel num_threads(cntl%ncpu)
!$omp master
#if defined(SPLLT_OMP_TRACE) 
    call trace_init(omp_get_num_threads())
    call trace_create_event('INIT_NODE', ini_nde_id)
    call trace_create_event('FACTO_BLK', fac_blk_id)
    call trace_create_event('SOLVE_BLK', slv_blk_id)
    call trace_create_event('UPDATE_BLK', upd_blk_id)
    call trace_create_event('UPDATE_BTW', upd_btw_id)
#endif

#elif defined(SPLLT_USE_PARSEC)
    call dague_init(cntl%ncpu, ctx)
#endif

    write(*,'("[>] factorize")')
    write(*,'("[>] [factorize]    nb: ", i6)') cntl%nb
    write(*,'("[>] [factorize] # cpu: ", i6)') cntl%ncpu

    ! factorize matrix
#if defined(SPLLT_USE_STF) || defined(SPLLT_USE_OMP) || defined(SPLLT_USE_STARPU) 
    call spllt_stf_factorize(a%n, a%ptr, a%row, a%val, order, keep, control, info, fdata, cntl)
#elif defined(SPLLT_USE_PARSEC)
    call spllt_ptg_factorize(adata, a%val, keep, cntl, fdata, info)
#endif
    ! call MA87_factor(a%n, a%ptr, a%row, a%val, order, keep, control, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a)") "failed factorization"
    endif

#if defined(SPLLT_USE_OMP)
!$omp taskwait
#elif (SPLLT_USE_PARSEC)
    call dague_context_wait(ctx)
#endif

#if defined(SPLLT_USE_OMP)
#if defined(SPLLT_OMP_TRACE) 
    call trace_log_dump_paje('trace')
#endif
!$omp end master
!$omp end parallel
#elif defined(SPLLT_USE_PARSEC)
    call dague_fini(ctx)
#endif

! goto 9999 ! no solve

    write(*,'("[>] solve")')

    x = b
    ! solve
    call MA87_solve(x, order, keep, control, info)
    if(info%flag .lt. spllt_success) then
       write(*, "(a,i4)") " fail on 1d solve", &
            info%flag
    endif
    
    call spllt_bwerr(a, x, b, resid)
    
    write(*,'("[>] [solve] bwderr ||Ax-b|| / (||A||||x|| + ||b||): ", es10.3)') resid    

    call MA87_finalise(keep, control)
    
9999 continue

    return
  end subroutine spllt_test_rand
  
end program spllt_test
