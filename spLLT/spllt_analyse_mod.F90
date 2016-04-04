module spllt_analyse_mod
  use spllt_mod
  implicit none

contains

  subroutine spllt_analyse(adata, n, ptr, row, order, cntl)
    use hsl_mc78_integer
    use hsl_mc34_double
    implicit none

    type(spllt_adata_type) :: adata
    integer, intent(in) :: n ! order of A
    integer, intent(in) :: row(:) ! row indices of lower triangular part
    integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
    integer, intent(inout), dimension(:) :: order
    ! order(i) must hold position of i in the pivot sequence. 
    ! On exit, holds the pivot order to be used by MA87_factor.
    type(spllt_cntl), intent(in) :: cntl

    integer, allocatable :: amap(:) ! map from user a to reordered a
    integer, allocatable :: aptr(:) ! column pointers of expanded matrix
    integer, allocatable :: arow(:) ! row pointers of expanded matrix
    integer, allocatable :: iw(:) ! work array
    
    integer :: num_nodes ! number of nodes in tree
    integer :: ne ! set to ptr(n+1) - 1
    integer :: nemin ! min. number of eliminations
    integer :: st ! stat parameter

    type(mc78_control) :: control78
    integer :: info78
    integer, dimension(:), allocatable :: sptr, sparent, rlist
    integer(long), dimension(:), allocatable :: rptr

    ne = ptr(n+1) - 1

    ! immediate return if n = 0
    if (n == 0) return

    ! expand the matrix
    
    ! allocate space for expanded matrix (held in aptr,arow)
    allocate (arow(2*ne-n),aptr(n+3),iw(n+1),amap(ptr(n+1)-1),stat=st)
    if (st /= 0) go to 9999

    arow(1:ne) = row(1:ne)
    aptr(1:n+1) = ptr(1:n+1)
    call mc34_expand(n, arow, aptr, iw)
    deallocate(iw,stat=st)

    nemin = cntl%nemin
    ! Check nemin (a node is merged with its parent if both involve
    ! fewer than nemin eliminations). If out of range, use the default
    if (nemin < 1) nemin = spllt_nemin_default
    
    control78%nemin = nemin
    control78%sort = .true.
    control78%lopt = .true.
    call mc78_analyse(n, aptr, arow, order, num_nodes, &
         sptr, sparent, rptr, rlist, control78, info78, nfact=adata%num_factor, &
         nflops=adata%num_flops)

    return

9999 continue

    ! TODO print error
    
    return
  end subroutine spllt_analyse

  ! tree pruning method. Inspired by the strategy employed in qr_mumps
  ! for pruning the atree
  subroutine spllt_prune_tree(num_nodes, sparent, nth)
    implicit none

    integer :: num_nodes ! number of nodes if the atree
    integer, dimension(:), intent(inout) :: sparent ! atree
    integer :: nth ! number of workers

    integer :: node, nlz
    integer, allocatable :: lzero(:)
    real(kind(1.d0)), allocatable   :: lzero_w(:), proc_w(:)

    ! initialize the l0 layer with the root nodes
    nlz = 0
    do node = 1,num_nodes
       if (sparent(node) .eq. 0) then
         nlz = nlz+1
         lzero(nlz) = node
       end if
    end do

    godown: do
       ! if (nth .eq. 1) exit ! serial execution process the whole tree as a subtree
       if(nlz .gt. nth*max(2.d0,(log(real(nth,kind(1.d0)))/log(2.d0))**2)) exit ! exit if already too many nodes in l0

       proc_w = 0.d0

    end do godown

    return
  end subroutine spllt_prune_tree

  ! subroutine spllt_node_flops()
  !   implicit none    

  !   return
  ! end subroutine spllt_node_flops

  subroutine spllt_symbolic(adata, nnodes, sptr, sparent, rptr)
    use spllt_mod
    implicit none

    type(spllt_adata_type), intent(inout) :: adata
    integer, intent(in) :: nnodes ! number of nodes if the atree
    integer, dimension(:), intent(in) :: sptr, sparent
    integer(long), dimension(:), intent(in) :: rptr

    integer :: node, parent
    integer(long) :: mm, m, n, j
    integer(long) :: nflops
    
    allocate(adata%weight(nnodes))

    do node = 1,nnodes

       parent = sparent(node)
       m = rptr(node+1) - rptr(node)
       n = sptr(node+1) - sptr(node)
       
       mm = m-n

       nflops = 0
       do j=1,n
          nflops = nflops + (mm+j)**2
       end do
       
       adata%weight(node) = adata%weight(node) + real(nflops, kind(1.d0)) 
       if (parent .gt. 0) then 
          adata%weight(parent) = adata%weight(parent) + adata%weight(node)
       end if

    end do

    return
  end subroutine spllt_symbolic

end module spllt_analyse_mod
