module spllt_mod
  use hsl_ma87_double, only: block_type, node_type 
  use hsl_zd11_double
#if defined(SPLLT_USE_STARPU)
  use iso_c_binding
  use starpu_f_mod
#elif defined(SPLLT_USE_PARSEC)
  use dague_f08_interfaces
#endif
  implicit none

  ! type :: matrix_type
  !    integer :: n, ne
  !    integer, dimension(:), allocatable :: ptr, row, col
  !    real(wp), dimension(:), allocatable :: val
  ! end type matrix_type

  integer, parameter :: wp = kind(0d0)
#if defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_PARSEC) 
  integer, parameter :: long = c_long
#else
  integer, parameter :: long = selected_int_kind(18)
#endif

  real(wp), parameter :: one = 1.0_wp
  real(wp), parameter :: zero = 0.0_wp

  ! error flags
  integer, parameter :: spllt_success           = 0 
  integer, parameter :: spllt_error_allocation  = -1
  integer, parameter :: spllt_error_unknown     = -99 

  ! Default values
  ! node amalgamation parameter
  integer, parameter :: spllt_nemin_default = 32
  ! Block size with dense kernel
  integer, parameter :: spllt_nb_default = 256

#if defined(SPLLT_USE_OMP) && defined(SPLLT_OMP_TRACE) 
  integer, save :: ini_nde_id, fac_blk_id, slv_blk_id, upd_blk_id, upd_btw_id 
#endif

#if defined(SPLLT_USE_PARSEC)
  type(dague_context_t) :: ctx
  integer(c_int)        :: nds
  integer(c_int)        :: rank
#endif

  type spllt_cntl
     integer :: ncpu = 1 ! number of CPU workers
     integer :: nb   = 16 ! blocking size
     integer :: nemin = 32 ! node amalgamation parameter
     logical :: prune_tree = .false.
     integer :: min_width_blas  = 8      ! Minimum width of source block
         ! before we use an indirect update_between
  end type spllt_cntl

  ! node type
  type spllt_node_type
     type(node_type), pointer :: node
     integer :: num ! node id
#if defined(SPLLT_USE_STARPU)
     type(c_ptr)    :: hdl  ! StarPU handle
#endif     
  end type spllt_node_type

  ! type spllt_dep_upd
  !    integer(long) :: id_kk  = 0
  !    integer(long) :: id_jk = 0
  !    integer(long) :: id_ik = 0
  !    integer(long) :: id_ij  = 0
  ! end type spllt_dep_upd

  ! input dependency
  type spllt_dep_upd_in
     integer(long) :: id_jk = 0
     integer :: p1 = 0
     integer(long) :: id_ik = 0
     integer :: p2 = 0
     integer :: csrc = 0 ! first row in Ljk
     integer :: rsrc = 0 ! first row in Lik
     logical :: sync = .true.
  end type spllt_dep_upd_in

  ! input dependency
  type spllt_dep_upd_out
     integer :: flow = 0 ! 1 -> ljk and 2 -> ljk  
     integer(long) :: id_ij = 0
     integer :: p = 0
  end type spllt_dep_upd_out

  ! type spllt_dep_node
  !    type(spllt_dep_upd)           :: upd
  !    type(spllt_dep_node), pointer :: prev => null()     
  !    type(spllt_dep_node), pointer :: next => null()     
  ! end type spllt_dep_node

  ! type spllt_dep_list
  !    type(spllt_dep_node), pointer :: head => null()
  !    type(spllt_dep_node), pointer :: tail => null()
  ! end type spllt_dep_list

  ! block type  
  type spllt_bc_type
     type(block_type), pointer :: blk => null()! pointer to the block in keep
     real(wp), pointer :: c(:)
#if defined(SPLLT_USE_STARPU)
     type(c_ptr)    :: hdl  ! StarPU handle
#endif
     integer :: mem_node = 0 ! memory node where the block is allocated
     

#if defined(SPLLT_USE_PARSEC)
     ! store ids of blocks contributing to this block
     ! type(spllt_dep_list), pointer :: dep_in  => null()
     type(spllt_dep_upd_in), pointer :: dep_in(:)  => null()
     ! store ids of blocks for which this block contributes to 
     ! type(spllt_dep_list), pointer :: dep_out => null()
     type(spllt_dep_upd_out), pointer :: dep_out(:) => null()
#endif
  end type spllt_bc_type
  
  ! problem data (analyis)
  type spllt_adata_type
     integer :: nnodes
     integer :: n
     integer(long) :: num_factor = 0_long ! Number of entries in the factor.
     integer(long) :: num_flops = 0_long  ! Number of flops for factor.
     ! weight(i): weight of the subtree rooted at node i where weight
     ! corresponds to the number of flops
     integer(long), allocatable :: weight(:)
     integer, allocatable :: small(:)
  end type spllt_adata_type

  ! problem data (factorization)
  type spllt_data_type
     type(spllt_bc_type), allocatable :: bc(:) ! blocks
#if defined(SPLLT_USE_OMP)
     type(spllt_bc_type), allocatable :: workspace(:) ! workspaces
#else
     type(spllt_bc_type) :: workspace ! workspaces
#endif
     type(spllt_node_type), allocatable :: nodes(:) ! super nodes 
#if defined(SPLLT_USE_PARSEC)
     ! ids of diag blocks. size keep%nbcol
     integer(long), allocatable :: diags(:)
#endif
  end type spllt_data_type

  interface gen_random_posdef
     module procedure gen_random_posdef
  end interface gen_random_posdef

  interface spllt_bwerr
     module procedure spllt_bwerr_1d
  end interface spllt_bwerr

  type spllt_options
     integer :: ncpu = 1 ! number of CPU workers
     integer :: nb   = 16 ! blocking size
     character(len=100) :: mat = ''
     integer :: nemin = -1
  end type spllt_options

contains

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

  subroutine spllt_print_atree(keep)
    use hsl_ma87_double
    implicit none    

    type(MA87_keep), target, intent(in) :: keep

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
       write(2, '("fillcolor=white ")', advance="no")
       write(2, '("label=""")', advance="no")  
       write(2, '("node:", i5,"\n")', advance="no")snode
       write(2, '("m:", i5,"\n")', advance="no")m
       write(2, '("n:", i5,"\n")', advance="no")n
       write(2, '("""")', advance="no")         
       write(2, '("]")', advance="no")
       write(2, '(" ")')

       if(keep%nodes(snode)%parent .ne. -1) write(2, '(i10, "--", i10)')keep%nodes(snode)%parent, snode
    end do

    write(2, '("}")')

    close(2)
    
    return
  end subroutine spllt_print_atree
  
  subroutine gen_random_posdef(matrix, nza, iseed)
    implicit none
    type(zd11_type), intent(inout) :: matrix
    integer, intent(in) :: nza
    integer, intent(inout) :: iseed

    integer, dimension(10) :: icntl
    character(len=8) :: key
    integer, dimension(:), allocatable :: work
    integer :: i, j, k
    real(wp) :: tempv

    ! Setup stuff for ym11
    key = 'nothing '
    call ym11id(icntl, i)
    icntl(3) = 0 ! Symmetric
    allocate(work(2*matrix%n))

    ! Generate matrix
    call ym11ad(matrix%n, matrix%n, nza, matrix%ne, matrix%row, &
         matrix%val, matrix%ptr, work, icntl, key, iseed)

    ! Make matrix diagonally dominant, observing first entry in column
    ! is always the diagonal when matrix generated by ym11
    do k = 1, matrix%n
       tempv = zero
       do j = matrix%ptr(k)+1, matrix%ptr(k+1)-1
          tempv = tempv + abs(matrix%val(j))
          i = matrix%ptr(matrix%row(j))
          matrix%val(i) = matrix%val(i) + abs(matrix%val(j))
       end do
       i = matrix%ptr(k)
       matrix%val(i) = one + matrix%val(i) + tempv
    end do
  end subroutine gen_random_posdef

  subroutine spllt_bwerr_1d(a,x,b,res)

    type(zd11_type), intent(in) :: a
    real(wp), dimension(:), allocatable :: x, b
    real(wp) :: res

    integer :: i, j, k
    real(wp), dimension(:), allocatable :: work, r
    real(wp) :: anorm, atemp, bnorm, xnorm

    allocate(work(a%n), r(a%n))

    anorm = zero
    bnorm = zero
    xnorm = zero
    work = zero
    
    r = b
    do k = 1, a%n
       do j = a%ptr(k), a%ptr(k+1)-1
          i = a%row(j)
          
          atemp = a%val(j)
          r(i) = r(i) - x(k)*atemp 
          work(k) = work(k) + abs(atemp)
          if(i.eq.k) cycle
          r(k) = r(k) - x(i)*atemp 
          work(i) = work(i) + abs(atemp)
       end do
    end do

    anorm = maxval(abs(work))
    bnorm = maxval(abs(b))
    xnorm = maxval(abs(x))
    res   = maxval(abs(r))

    res = res/(anorm*xnorm + bnorm)

  end subroutine spllt_bwerr_1d
  
  subroutine compute_resid(nrhs,a,x,lx,rhs,lrhs,res)

    integer, intent(in) :: nrhs, lrhs, lx
    type(zd11_type), intent(in) :: a
    real(wp), intent(in) :: rhs(lrhs,nrhs)
    real(wp), intent(in) :: x(lx,nrhs)
    real(wp), dimension(:,:), allocatable, intent(inout) :: res
    real(wp), dimension(:), allocatable :: work

    integer :: i, j, k
    real(wp) :: anorm, atemp, bnorm(1:nrhs), xnorm(1:nrhs)

    if (allocated(res)) deallocate(res)
    allocate(res(a%n,nrhs),work(a%n))

    anorm = zero
    bnorm = zero
    xnorm = zero
    work = zero

    ! Check residual
    res(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
    do k = 1, a%n
       do j = a%ptr(k), a%ptr(k+1)-1
          i = a%row(j)
          atemp = a%val(j)
          res(i, 1:nrhs) = res(i, 1:nrhs) - atemp*x(k,1:nrhs)
          work(i) = work(i) + abs(atemp)
          if(i.eq.k) cycle
          res(k, 1:nrhs) = res(k, 1:nrhs) - atemp*x(i,1:nrhs)
          work(k) = work(k) + abs(atemp)
       end do
    end do

    do k = 1, a%n
       anorm = max(anorm,work(k))
       do i = 1,nrhs
          bnorm(i) = max(bnorm(i),abs(rhs(k,i)))
          xnorm(i) = max(xnorm(i),abs(x(k,i)))
       end do
    end do

    do k = 1,a%n
       do i = 1,nrhs
          res(k,i) = res(k,i)/(anorm*xnorm(i) + bnorm(i))
       end do
    end do

  end subroutine compute_resid

  subroutine spllt_print_err(iflag, context, st)
    use hsl_ma87_double
    implicit none
    
    integer, intent(in) :: iflag
    character (len=*), optional, intent(in) :: context
    integer, intent(in), optional :: st

    select case(iflag)
    case(spllt_error_allocation)
       write(*,*) 'allocation error'
    case default
       write(*,*) 'unknown error'
    end select

    return
  end subroutine spllt_print_err

  subroutine amd_order(a,order)
    type(zd11_type), intent(in) :: a
    integer, dimension(:), allocatable :: order

    logical :: realloc_flag
    integer :: i, st
    integer, dimension(10) :: icntl, info
    real(wp), dimension(10) :: rinfo
    integer, dimension(:), allocatable :: work, ptr

    realloc_flag = .true.
    if(allocated(order)) realloc_flag = size(order).lt.a%n

    if(realloc_flag) then
       deallocate(order,stat=st)
       allocate(order(a%n))
    endif

    ! Initilise control
    call mc47id(icntl)
    icntl(1:2) = -1 ! Supress warnings and errors
    icntl(5) = huge(0) ! Largest integer

    ! Copy ptr data to work array
    allocate(ptr(a%n+1))
    ptr(:) = a%ptr(1:a%n+1)
    ! Copy row data to work array
    allocate(work(2*a%ptr(a%n+1) + 10*a%n))
    work(1:a%ptr(a%n+1)-1) = &
         a%row(1:a%ptr(a%n+1)-1)

    ! Perform AMD
    call mc47ad(a%n, a%ptr(a%n+1)-1, ptr, work, &
         size(work), icntl, info, rinfo)
    if(info(1).lt.0) then
       ! Failed for some reason
       do i = 1, a%n
          order(i) = i
       end do
       return
    endif

    ! Extract ordering
    order(1:a%n) = work(size(work)-a%n+1:size(work))

  end subroutine amd_order

  subroutine splllt_parse_args(options)
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
       case default
          write(*,'("Unrecognised command line argument: ", a20)'), argval
       end select
    end do

  end subroutine splllt_parse_args

end module spllt_mod
