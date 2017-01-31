module spllt_mod
  use spllt_data_mod
  use hsl_ma87_double, only: block_type, node_type 
  ! use hsl_zd11_double
#if defined(SPLLT_USE_STARPU)
  use iso_c_binding
  use starpu_f_mod
#elif defined(SPLLT_USE_PARSEC)
  use dague_f08_interfaces
#endif
  implicit none
  
  ! interface gen_random_posdef
  !    module procedure gen_random_posdef
  ! end interface gen_random_posdef

  ! interface spllt_bwerr
  !    module procedure spllt_bwerr_1d
  ! end interface spllt_bwerr

contains

  ! initialize solver
  subroutine spllt_init(cntl)
#if defined(SPLLT_USE_STARPU)
    use iso_c_binding
    use starpu_f_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE) 
    use trace_mod
#endif
#elif defined(SPLLT_USE_PARSEC)
    use dague_f08_interfaces
    use spllt_parsec_mod
#endif

#if defined(SPLLT_USE_GPU)
    use magma
#endif
    implicit none

    type(spllt_cntl) :: cntl
#if defined(SPLLT_USE_STARPU)
    integer :: start_starpuinit_t, stop_starpuinit_t, rate_starpuinit_t
    integer(c_int) :: ret
#endif

#if defined(SPLLT_USE_STARPU)

    call system_clock(start_starpuinit_t, rate_starpuinit_t)
    ! initialize starpu
    ret = starpu_f_init(cntl%ncpu)
    call system_clock(stop_starpuinit_t)
    write(*,'("[>] [spllt_test_mat] StarPU init time: ", es10.3, " s")') &
         &(stop_starpuinit_t - start_starpuinit_t)/real(rate_starpuinit_t)
    call starpu_f_fxt_start_profiling()

#elif defined(SPLLT_USE_OMP)

#if defined(SPLLT_OMP_TRACE) 

    call trace_init(omp_get_num_threads())
    call trace_create_event('INIT_NODE', ini_nde_id)
    call trace_create_event('FACTO_BLK', fac_blk_id)
    call trace_create_event('SOLVE_BLK', slv_blk_id)
    call trace_create_event('UPDATE_BLK', upd_blk_id)
    call trace_create_event('UPDATE_BTW', upd_btw_id)

#endif

#elif defined(SPLLT_USE_PARSEC)

    ctx = parsec_init(cntl%ncpu, nds, rank)
    write(*,'("Parsec init    nodes: ", i6, ", rank: ", i6)') nds, rank
    ! call dague_init(cntl%ncpu, ctx)
#endif

#if defined(SPLLT_USE_GPU)
    call magmaf_init()
#endif

    return
  end subroutine spllt_init

  subroutine spllt_finalize()
#if defined(SPLLT_USE_STARPU)
    use iso_c_binding
    use starpu_f_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE) 
    use trace_mod
#endif
#elif defined(SPLLT_USE_PARSEC)
    use dague_f08_interfaces
    use spllt_parsec_mod
#endif
    ! use spllt_factorization_task_mod 
    implicit none

    integer :: start_starpushutdown_t, stop_starpushutdown_t, rate_starpushutdown_t

#if defined(SPLLT_USE_GPU)
    call magmaf_finalize()
#endif

#if defined(SPLLT_USE_STARPU)

    call system_clock(start_starpushutdown_t, rate_starpushutdown_t)
    call starpu_f_shutdown()
    call system_clock(stop_starpushutdown_t)
    write(*,'("[>] [spllt_test_mat] StarPU shutdown time: ", es10.3, " s")') &
         &(stop_starpushutdown_t - start_starpushutdown_t)/real(rate_starpushutdown_t)

#elif defined(SPLLT_USE_OMP)

#if defined(SPLLT_OMP_TRACE) 
    call trace_log_dump_paje('trace')
#endif

#elif defined(SPLLT_USE_PARSEC)

    ! call dague_fini(ctx)
    call parsec_fini(ctx)

#endif
    return
  end subroutine spllt_finalize

  ! subroutine spllt_dep_array_init(dep_arr, dep)
  !   implicit none

  !   type(spllt_dep_upd), dimension(:), pointer :: dep_arr
  !   type(spllt_dep_upd) :: dep

  !   if (.not. associated(dep_arr)) then
  !      allocate(dep_arr(1))
  !      dep_arr(1) = dep
  !   end if
    
  !   return
  ! end subroutine spllt_dep_array_init

  function spllt_dep_in_add(dep_in_arr, id_jk, id_ik, csrc, rsrc, sync)
    implicit none

    type(spllt_dep_upd_in), dimension(:), pointer :: dep_in_arr
    integer(long) :: id_jk, id_ik
    integer :: csrc, rsrc
    integer :: spllt_dep_in_add
    logical, optional :: sync

    integer :: sz
    logical :: s
    type(spllt_dep_upd_in), dimension(:), pointer :: new_dep_in_arr => null()

    if (present(sync)) then
       s = sync
    else
       s = .true.
    end if
    
    if (.not. associated(dep_in_arr)) then
       allocate(dep_in_arr(1))
       dep_in_arr(1)%id_jk = id_jk
       dep_in_arr(1)%id_ik = id_ik
       dep_in_arr(1)%csrc = csrc
       dep_in_arr(1)%rsrc = rsrc
       dep_in_arr(1)%sync = s
    else
       sz = size(dep_in_arr)
       allocate(new_dep_in_arr(sz+1))
       new_dep_in_arr(1:sz) = dep_in_arr(1:sz)       
       new_dep_in_arr(sz+1)%id_jk = id_jk
       new_dep_in_arr(sz+1)%id_ik = id_ik
       new_dep_in_arr(sz+1)%csrc = csrc
       new_dep_in_arr(sz+1)%rsrc = rsrc
       new_dep_in_arr(sz+1)%sync = s
       deallocate(dep_in_arr)
       dep_in_arr => new_dep_in_arr
    end if

    spllt_dep_in_add = size(dep_in_arr)

  end function spllt_dep_in_add

  function spllt_dep_out_add(dep_out_arr, id_ij, flow)
    implicit none

    type(spllt_dep_upd_out), dimension(:), pointer :: dep_out_arr
    integer(long) :: id_ij
    integer :: flow
    integer :: spllt_dep_out_add

    integer :: sz
    type(spllt_dep_upd_out), dimension(:), pointer :: new_dep_out_arr => null()
    
    if (.not. associated(dep_out_arr)) then
       allocate(dep_out_arr(1))
       dep_out_arr(1)%id_ij = id_ij
       dep_out_arr(1)%flow = flow
    else
       sz = size(dep_out_arr)
       allocate(new_dep_out_arr(sz+1))
       new_dep_out_arr(1:sz) = dep_out_arr(1:sz)       
       new_dep_out_arr(sz+1)%id_ij = id_ij
       new_dep_out_arr(sz+1)%flow = flow
       deallocate(dep_out_arr)
       dep_out_arr => new_dep_out_arr
    end if

    spllt_dep_out_add = size(dep_out_arr)

  end function spllt_dep_out_add

  ! subroutine spllt_dep_array_add(dep_arr, id_kk, id_jk, id_ik, id_ij)
  !   implicit none

  !   type(spllt_dep_upd), dimension(:), pointer :: dep_arr
  !   integer(long) :: id_kk, id_jk, id_ik, id_ij

  !   integer :: sz
  !   type(spllt_dep_upd), dimension(:), pointer :: new_dep_arr => null()

  !   sz = size(dep_arr)
    
  !   if (.not. associated(dep_arr)) then
  !      allocate(dep_arr(1))
  !      dep_arr(1)%id_kk = id_kk
  !      dep_arr(1)%id_jk = id_jk
  !      dep_arr(1)%id_ik = id_ik
  !      dep_arr(1)%id_ij = id_ij
  !   else
  !      allocate(new_dep_arr(sz+1))
  !      new_dep_arr(1:sz) = dep_arr(1:sz)
  !      new_dep_arr(sz+1)%id_kk = id_kk
  !      new_dep_arr(sz+1)%id_jk = id_jk
  !      new_dep_arr(sz+1)%id_ik = id_ik
  !      new_dep_arr(sz+1)%id_ij = id_ij       
  !      deallocate(dep_arr)
  !      dep_arr => new_dep_arr
  !   end if
    
  !   return
  ! end subroutine spllt_dep_array_add

  ! subroutine spllt_dep_list_init(dep_list, dep)
  !   implicit none

  !   type(spllt_dep_list), pointer :: dep_list
  !   type(spllt_dep_node), pointer :: dep

  !   if (.not. associated(dep_list)) then
  !      allocate(dep_list)
  !      dep_list%head => dep
  !      dep_list%tail => dep
  !   end if

  !   return
  ! end subroutine spllt_dep_list_init

  ! subroutine spllt_dep_list_push_back(dep_list, id_lik, id_ljk, id)
  !   implicit none

  !   type(spllt_dep_list), pointer :: dep_list
  !   integer(long) :: id_lik, id_ljk, id

  !   type(spllt_dep_node), pointer :: dep => null()

  !   allocate(dep)
    
  !   dep%upd%id_lik = id_lik
  !   dep%upd%id_ljk = id_ljk
  !   dep%upd%id     = id

  !   if (.not. associated(dep_list)) then
  !      call spllt_dep_list_init(dep_list, dep)
  !   else
  !      dep%prev           => dep_list%tail
  !      dep_list%tail%next => dep
  !      dep_list%tail      => dep
  !   end if

  !   return
  ! end subroutine spllt_dep_list_push_back

  ! subroutine spllt_realloc_1d(a, n)
  !   implicit none

  !   if (.not. allocated(a)) then
  !      allocate(a(n))
  !      return
  !   else
  !      if (size(a) .lt. n) then
  !         deallocate(a)
  !         allocate(a(n))
  !      end if
  !   end if
    
  !   return
  ! end subroutine spllt_realloc_1d

  subroutine spllt_print_atree(adata, keep, cntl)
    use spllt_data_mod
    use hsl_ma87_double
    implicit none    

    type(spllt_adata_type), intent(in)  :: adata
    type(MA87_keep), target, intent(in) :: keep
    type(spllt_cntl), intent(in) :: cntl

    integer :: snode, num_nodes
    type(node_type), pointer     :: node ! node in the atree
    integer :: m, n

    num_nodes = keep%info%num_nodes

    open(2, file="atree.dot")

    write(2, '("graph atree {")')
    write(2, '("node [")')
    write(2, '("style=filled")')
    write(2, '("]")')
    
    do snode=1,num_nodes

       node => keep%nodes(snode)
       m = size(node%index)
       n = node%en - node%sa + 1

       write(2, '(i10)', advance="no")snode
       write(2, '(" ")', advance="no")
       write(2, '("[")', advance="no")
       if (cntl%prune_tree) then
          if (adata%small(snode) .eq. 1) then
             write(2, '("fillcolor=grey ")', advance="no")
          else
             write(2, '("fillcolor=white ")', advance="no")
          endif
       else
          write(2, '("fillcolor=white ")', advance="no")
       endif
       write(2, '("label=""")', advance="no")  
       write(2, '("node:", i5,"\n")', advance="no")snode
       write(2, '("m:", i5,"\n")', advance="no")m
       write(2, '("n:", i5,"\n")', advance="no")n
       if (cntl%prune_tree) then
          write(2, '("small:", i5,"\n")', advance="no")adata%small(snode)
          write(2, '("weight:", f5.1,"\n")', advance="no") &
               & 100.0 * (real(adata%weight(snode), kind(1d0)) / real(adata%weight(num_nodes+1), kind(1d0)))
       end if
       write(2, '("""")', advance="no")         
       write(2, '("]")', advance="no")
       write(2, '(" ")')

       if(keep%nodes(snode)%parent .ne. -1) write(2, '(i10, "--", i10)')keep%nodes(snode)%parent, snode
    end do

    write(2, '("}")')

    close(2)
    
    return
  end subroutine spllt_print_atree
  
  ! subroutine gen_random_posdef(matrix, nza, iseed)
  !   implicit none
  !   type(zd11_type), intent(inout) :: matrix
  !   integer, intent(in) :: nza
  !   integer, intent(inout) :: iseed

  !   integer, dimension(10) :: icntl
  !   character(len=8) :: key
  !   integer, dimension(:), allocatable :: work
  !   integer :: i, j, k
  !   real(wp) :: tempv

  !   ! Setup stuff for ym11
  !   key = 'nothing '
  !   call ym11id(icntl, i)
  !   icntl(3) = 0 ! Symmetric
  !   allocate(work(2*matrix%n))

  !   ! Generate matrix
  !   call ym11ad(matrix%n, matrix%n, nza, matrix%ne, matrix%row, &
  !        matrix%val, matrix%ptr, work, icntl, key, iseed)

  !   ! Make matrix diagonally dominant, observing first entry in column
  !   ! is always the diagonal when matrix generated by ym11
  !   do k = 1, matrix%n
  !      tempv = zero
  !      do j = matrix%ptr(k)+1, matrix%ptr(k+1)-1
  !         tempv = tempv + abs(matrix%val(j))
  !         i = matrix%ptr(matrix%row(j))
  !         matrix%val(i) = matrix%val(i) + abs(matrix%val(j))
  !      end do
  !      i = matrix%ptr(k)
  !      matrix%val(i) = one + matrix%val(i) + tempv
  !   end do
  ! end subroutine gen_random_posdef

  ! subroutine spllt_bwerr_1d(a,x,b,res)

  !   type(zd11_type), intent(in) :: a
  !   real(wp), dimension(:), allocatable :: x, b
  !   real(wp) :: res

  !   integer :: i, j, k
  !   real(wp), dimension(:), allocatable :: work, r
  !   real(wp) :: anorm, atemp, bnorm, xnorm

  !   allocate(work(a%n), r(a%n))

  !   anorm = zero
  !   bnorm = zero
  !   xnorm = zero
  !   work = zero
    
  !   r = b
  !   do k = 1, a%n
  !      do j = a%ptr(k), a%ptr(k+1)-1
  !         i = a%row(j)
          
  !         atemp = a%val(j)
  !         r(i) = r(i) - x(k)*atemp 
  !         work(k) = work(k) + abs(atemp)
  !         if(i.eq.k) cycle
  !         r(k) = r(k) - x(i)*atemp 
  !         work(i) = work(i) + abs(atemp)
  !      end do
  !   end do

  !   anorm = maxval(abs(work))
  !   bnorm = maxval(abs(b))
  !   xnorm = maxval(abs(x))
  !   res   = maxval(abs(r))

  !   res = res/(anorm*xnorm + bnorm)

  ! end subroutine spllt_bwerr_1d
  
  ! subroutine compute_resid(nrhs,a,x,lx,rhs,lrhs,res)

  !   integer, intent(in) :: nrhs, lrhs, lx
  !   type(zd11_type), intent(in) :: a
  !   real(wp), intent(in) :: rhs(lrhs,nrhs)
  !   real(wp), intent(in) :: x(lx,nrhs)
  !   real(wp), dimension(:,:), allocatable, intent(inout) :: res
  !   real(wp), dimension(:), allocatable :: work

  !   integer :: i, j, k
  !   real(wp) :: anorm, atemp, bnorm(1:nrhs), xnorm(1:nrhs)

  !   if (allocated(res)) deallocate(res)
  !   allocate(res(a%n,nrhs),work(a%n))

  !   anorm = zero
  !   bnorm = zero
  !   xnorm = zero
  !   work = zero

  !   ! Check residual
  !   res(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
  !   do k = 1, a%n
  !      do j = a%ptr(k), a%ptr(k+1)-1
  !         i = a%row(j)
  !         atemp = a%val(j)
  !         res(i, 1:nrhs) = res(i, 1:nrhs) - atemp*x(k,1:nrhs)
  !         work(i) = work(i) + abs(atemp)
  !         if(i.eq.k) cycle
  !         res(k, 1:nrhs) = res(k, 1:nrhs) - atemp*x(i,1:nrhs)
  !         work(k) = work(k) + abs(atemp)
  !      end do
  !   end do

  !   do k = 1, a%n
  !      anorm = max(anorm,work(k))
  !      do i = 1,nrhs
  !         bnorm(i) = max(bnorm(i),abs(rhs(k,i)))
  !         xnorm(i) = max(xnorm(i),abs(x(k,i)))
  !      end do
  !   end do

  !   do k = 1,a%n
  !      do i = 1,nrhs
  !         res(k,i) = res(k,i)/(anorm*xnorm(i) + bnorm(i))
  !      end do
  !   end do

  ! end subroutine compute_resid

  ! subroutine amd_order(a,order)
  !   type(zd11_type), intent(in) :: a
  !   integer, dimension(:), allocatable :: order

  !   logical :: realloc_flag
  !   integer :: i, st
  !   integer, dimension(10) :: icntl, info
  !   real(wp), dimension(10) :: rinfo
  !   integer, dimension(:), allocatable :: work, ptr

  !   realloc_flag = .true.
  !   if(allocated(order)) realloc_flag = size(order).lt.a%n

  !   if(realloc_flag) then
  !      deallocate(order,stat=st)
  !      allocate(order(a%n))
  !   endif

  !   ! Initilise control
  !   call mc47id(icntl)
  !   icntl(1:2) = -1 ! Supress warnings and errors
  !   icntl(5) = huge(0) ! Largest integer

  !   ! Copy ptr data to work array
  !   allocate(ptr(a%n+1))
  !   ptr(:) = a%ptr(1:a%n+1)
  !   ! Copy row data to work array
  !   allocate(work(2*a%ptr(a%n+1) + 10*a%n))
  !   work(1:a%ptr(a%n+1)-1) = &
  !        a%row(1:a%ptr(a%n+1)-1)

  !   ! Perform AMD
  !   call mc47ad(a%n, a%ptr(a%n+1)-1, ptr, work, &
  !        size(work), icntl, info, rinfo)
  !   if(info(1).lt.0) then
  !      ! Failed for some reason
  !      do i = 1, a%n
  !         order(i) = i
  !      end do
  !      return
  !   endif

  !   ! Extract ordering
  !   order(1:a%n) = work(size(work)-a%n+1:size(work))

  ! end subroutine amd_order

  subroutine spllt_parse_args(options)
    implicit none

    type(spllt_options), intent(inout) :: options
    
    integer :: argnum, narg
    character(len=200) :: argval

    narg = command_argument_count()
    argnum = 1
    do while(argnum <= narg)
       call get_command_argument(argnum, argval)
       argnum = argnum + 1
       select case(argval)
       case("--nb")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%nb
       case("--ncpu")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%ncpu
       case("--mat")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%mat
       case("--nemin")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read( argval, * ) options%nemin
       case("--prune-tree")
          options%prune_tree = .true.
       case default
          write(*,'("Unrecognised command line argument: ", a20)'), argval
       end select
    end do

  end subroutine spllt_parse_args

  

end module spllt_mod
