module spllt_data_mod
#if defined(SPLLT_USE_STARPU)
  use, intrinsic :: iso_c_binding
  use starpu_f_mod
#elif defined(SPLLT_USE_PARSEC)
  use parsec_f08_interfaces
#endif
  use spral_ssids_inform, only: ssids_inform
  implicit none

  integer, parameter :: k_dep = 2 ! Upper bound on the number of dep 
                                  !   of a normal block 
  logical, parameter :: use_omp_cases_method = .true.

  integer, parameter :: wp = kind(0d0)
#if defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_PARSEC) 
  integer, parameter :: long = c_long
#else
  integer, parameter :: long = selected_int_kind(18)
#endif

  real(wp), parameter :: one = 1.0_wp
  real(wp), parameter :: zero = 0.0_wp

  ! Error flags
  integer, parameter :: SPLLT_SUCCESS              = 0 
  integer, parameter :: SPLLT_ERROR_ALLOCATION     = -1
  integer, parameter :: SPLLT_WARNING_PARAM_VALUE  = -10
  integer, parameter :: SPLLT_ERROR_UNIMPLEMENTED  = -98
  integer, parameter :: SPLLT_ERROR_UNKNOWN        = -99 

  ! Default values
  integer, parameter :: nemin_default = 32 ! node amalgamation parameter
  integer, parameter :: nb_default = 256 ! block size with dense kernel

!#if defined(SPLLT_USE_OMP) && defined(SPLLT_OMP_TRACE) 
  integer, save :: ini_nde_id, fac_blk_id, slv_blk_id, upd_blk_id, upd_btw_id 
!#endif

#if defined(SPLLT_USE_PARSEC)
  type(parsec_context_t) :: ctx
  integer(c_int) :: nds
  integer(c_int) :: rank
#endif

  ! !*************************************************
  ! !  
  ! ! Data type for storing information for each block (BLK)
  ! ! The blocks are numbered 1,2,..., keep%final_blk
  ! type block_type
  !    ! Static info, which is set in ma87_analayse
  !    integer :: bcol            ! block column that blk belongs to
  !    integer :: blkm            ! height of block (number of rows in blk)
  !    integer :: blkn            ! width of block (number of columns in blk)
  !    integer(long) :: dblk      ! id of the block on the diagonal within the 
  !    ! block column to which blk belongs
  !    integer :: dep_initial     ! initial dependency count for block,
  !    integer(long) :: id        ! The block identitifier (ie, its number blk)
  !    integer(long) :: last_blk  ! id of the last block within the
  !    ! block column to which blk belongs
  !    integer :: node            ! node to which blk belongs
  !    integer :: sa              ! posn of the first entry of the
  !    ! block blk within the array that holds the block column of L
  !    ! that blk belongs to.

  !    ! Non-static info
  !    integer :: dep  ! dependency countdown/marker. Once factor or solve done,
  !    ! value is -2.
  !    ! for this block.
  !    ! Note: locks initialised in ma87_analyse and destroyed
  !    !       in ma87_finalise
  ! end type block_type


  ! useful type for representing dependencies betwewen blocks

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

  !*************************************************  
  !
  ! Data type that represents a single block column in L
  !
  type lfactor
     real(wp), dimension(:), allocatable :: lcol ! holds block column
  end type lfactor

  !*************************************************  
  ! Data type for storing mapping from user's matrix values into
  ! block columns of L
  type lmap_type
     integer(long) :: len_map ! length of map
     integer(long), allocatable :: map(:,:) ! holds map from user's val
     ! array into lfact(:)%lcol values as follows:
     ! lcol( map(1,i) ) += val( map(2,i) )     i = 1:lmap
     ! map is set at end of analyse phase using subroutines make_map
     ! and lcol_map, and is used in factorization phase by blk_col_add_a
  end type lmap_type

  ! block type  
  type spllt_block
     ! type(block_type), pointer :: blk => null()! pointer to the block in keep
     real(wp), pointer :: c(:)
#if defined(SPLLT_USE_STARPU)
     type(c_ptr)    :: hdl  ! StarPU handle
     type(c_ptr)    :: hdl2  ! DEBUG
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
     ! Static info, which is set in ma87_analyse
     integer :: bcol            ! block column that blk belongs to
     integer :: blkm            ! height of block (number of rows in blk)
     integer :: blkn            ! width of block (number of columns in blk)
     integer(long) :: dblk      ! id of the block on the diagonal within the 
     ! block column to which blk belongs
     integer :: dep_initial     ! initial dependency count for block,
     integer(long) :: id        ! The block identifier (ie, its number blk)
     integer(long) :: last_blk  ! id of the last block within the
     ! block column to which blk belongs
     integer :: node            ! node to which blk belongs
     integer :: sa              ! posn of the first entry of the
     ! block blk within the array that holds the block column of L
     ! that blk belongs to.

     ! Non-static info
     integer :: dep  ! dependency countdown/marker. Once factor or solve done,
     ! value is -2.
     ! for this block.
     ! Note: locks initialised in ma87_analyse and destroyed
     !       in ma87_finalise

     ! Additional components to handle the list of dependencies of the block
     ! List of blk indices dependencies in :
     !  - the forward step of the solve
     integer, allocatable :: fwd_dep(:)
     integer, allocatable :: fwd_update_dep(:)
     integer              :: fwd_solve_dep
     !  - the backward step of the solve
     integer              :: bwd_update_dep
     integer, allocatable :: bwd_solve_dep(:)
     integer, allocatable :: bwd_dep(:)

  end type spllt_block

  ! node type
  type spllt_node
     integer :: num ! node id
#if defined(SPLLT_USE_STARPU)
     type(c_ptr)    :: hdl  ! StarPU handle
     type(c_ptr)    :: hdl2  ! StarPU second handle
#endif
     type(spllt_block) :: buffer ! buffer for accumulating updates

     integer(long) :: blk_sa ! identifier of the first block in node
     integer(long) :: blk_en ! identifier of the last block in node

     integer :: nb ! Block size for nodal matrix
     ! If number of cols nc in nodal matrix is less than control%nb but 
     ! number of rows is large, the block size for the node is taken as 
     ! control%nb**2/nc, rounded up to a multiple of 8. The aim is for
     ! the numbers of entries in the blocks to be similar to those in the 
     ! normal case. 

     integer :: sa ! index (in pivotal order) of the first column of the node
     integer :: en ! index (in pivotal order) of the last column of the node

     integer, allocatable :: index(:) ! holds the permuted variable
     ! list for node. They are sorted into increasing order.
     ! index is set up by ma87_analyse

     integer :: nchild ! number of children node has in assembly tree
     integer, allocatable :: child(:) ! holds children of node
     integer :: parent ! Parent of node in assembly tree
     integer :: least_desc ! Least descendant in assembly tree

     ! List of rows of the node that are not present in the children
     integer, allocatable :: extra_row(:)
  end type spllt_node

  type spllt_workspace_i
     integer, pointer :: c(:)
#if defined(SPLLT_USE_STARPU)
     type(c_ptr)    :: hdl  ! StarPU handle
#endif
  end type spllt_workspace_i

  !*************************************************
  !
  ! Data type for control parameters
  type spllt_options
     integer :: print_level = 0 ! Controls diagnostic printing.
     ! Possible values are:
     !  < 0: no printing.
     !  0: error and warning messages only.
     !  1: as 0 plus basic diagnostic printing.
     !  > 1: as 1 plus some more detailed diagnostic messages.
     !  > 9999: debug (absolutely everything - really don't use this)
     integer :: ncpu = 1             ! Number of CPU workers
     integer :: nb   = 16            ! Blocking size
     character(len=100) :: mat = ''
     integer :: nemin = 32           ! nemin parameter for analysis
     logical :: prune_tree = .true.  ! Tree pruning
     character(len=3) :: fmt='csc'
     integer :: min_width_blas  = 8  ! Minimum width of source block
     ! before we use an indirect update_between
     integer :: nb_min = 32
     integer :: nb_max = 32
     integer :: nrhs_min = 1
     integer :: nrhs_max = 1
     logical :: nb_linear_comp = .false.
     logical :: nrhs_linear_comp = .false.
     logical :: ileave_solve = .false.
  end type spllt_options

  type spllt_tree_t
    integer :: num      ! Tree id
    integer :: nnode    ! # nodes in the tree
    integer :: node_sa  ! First node of the tree
    integer :: node_en  ! Last node of the tree
    integer :: nchild   ! #children tree
    integer, allocatable :: child(:) ! Holds children of the tree
    integer :: parent   ! Tree id the parent
  end type spllt_tree_t

  !*************************************************
  !
  ! data type for returning information to user.
  type spllt_inform
     type(ssids_inform) :: ssids_inform ! SSDIS stats
     integer :: flag = SPLLT_SUCCESS  ! Error return flag (0 on success)
     integer :: maxdepth = 0           ! Maximum depth of the tree.
     integer(long) :: num_factor = 0_long ! Number of entries in the factor.
     integer(long) :: num_flops = 0_long  ! Number of flops for factor.
     integer :: num_nodes = 0          ! Number of nodes
     integer :: stat = 0               ! STAT value on error return -1.
  end type spllt_inform

  !*************************************************
  !
  ! Data associated with input matrix know after analysis phase
  !
  type spllt_akeep
     !> Number of nodes in the assembly tree.
     integer :: nnodes
     !> Oder of the system.
     integer :: n
     !> Number of entries in the factor.
     integer(long) :: num_factor = 0_long
     integer(long) :: num_flops = 0_long  ! Number of flops for factor.
     ! weight(i): weight of the subtree rooted at node i where weight
     ! corresponds to the number of flops
     integer(long), allocatable :: weight(:)
     integer, allocatable :: small(:)
  end type spllt_akeep

  !*************************************************
  !
  ! Data type for data generated in factorise phase
  !
  type spllt_fkeep
     type(spllt_block), allocatable :: bc(:) ! blocks
#if defined(SPLLT_USE_OMP)
     type(spllt_block), allocatable :: workspace(:) ! workspaces
#else
     type(spllt_block) :: workspace ! workspaces
#endif
     type(spllt_node), allocatable :: nodes(:) ! supernodes 
#if defined(SPLLT_USE_PARSEC)
     ! ids of diag blocks. size keep%nbcol
     integer(long), allocatable :: diags(:)
     ! data descriptor
     type(c_ptr) :: ddesc
#endif

#if defined(SPLLT_USE_OMP)
     type(spllt_workspace_i), allocatable :: row_list(:), col_list(:) ! workspaces
#else
     type(spllt_workspace_i) :: row_list, col_list ! workspaces
#endif

#if defined(SPLLT_USE_OMP)
     type(spllt_workspace_i), allocatable :: map(:)
#else
     type(spllt_workspace_i) :: map ! workspaces
#endif
     !     private
     ! type(spllt_block), dimension(:), allocatable :: blocks ! block info
     integer, dimension(:), allocatable :: flag_array ! allocated to
     ! have size equal to the number of threads. For each thread, holds
     ! error flag
     integer(long) :: final_blk = 0 ! Number of blocks. Used for destroying
     ! locks in finalise
     type(spllt_inform) :: info ! Holds copy of info
     integer :: ndblk
     integer, allocatable :: rhsPtr(:)
     integer, allocatable :: indir_rhs(:)
     integer :: maxmn ! holds largest block dimension
     integer :: n  ! Order of the system.
     ! type(node_type), dimension(:), allocatable :: nodes ! nodal info
     integer :: nbcol = 0 ! number of block columns in L
     type(lfactor), dimension(:), allocatable :: lfact
     ! holds block cols of L
     type(lmap_type), dimension(:), allocatable :: lmap
     ! holds mapping from matrix values into lfact

!    logical, allocatable :: workspace_reset(:)
     type(spllt_tree_t), allocatable :: trees(:)
    !type(spllt_akeep), pointer :: p_akeep(:)
     integer, allocatable :: small(:)
     integer, allocatable :: assoc_tree(:)
  end type spllt_fkeep


  !*************************************************
  !  
!   ! Data type for user controls
!   type spllt_control

!      integer :: diagnostics_level = 0      ! Controls diagnostic printing.
!      ! Possible values are:
!      !  < 0: no printing.
!      !    0: error and warning messages only.
!      !    1: as 0 plus basic diagnostic printing.
!      !    2: as 1 plus some more detailed diagnostic messages.
!      !    3: as 2 plus all entries of user-supplied arrays.
!      integer :: nb    = spllt_nb_default ! Controls the size of the
!      ! blocks used within each node (used to set nb within node_type)
!      integer :: nemin = spllt_nemin_default    
!      ! Node amalgamation parameter. A child node is merged with its parent 
!      ! if they both involve fewer than nemin eliminations.
!      integer :: unit_diagnostics = 6    ! unit for diagnostic messages
!      ! Printing is suppressed if unit_diagnostics  <  0.
!      integer :: unit_error       = 6    ! unit for error messages
!      ! Printing is suppressed if unit_error  <  0.
!      integer :: unit_warning     = 6    ! unit for warning messages
!      ! Printing is suppressed if unit_warning  <  0.


! !!!! Undocumented
!      !**   integer :: time_out        = -1     ! If >= 0 some internal timings
!      !**      are printed on unit time_out. For HSL 2011 these are commented
!      !**      using comments !** so easy to search and uncomment
!      !%%%  integer :: unit_log        = -1     ! For profiling log output
!      !%%%     commented out for HSL 2011 using !%%%
!      !%%%  integer :: log_level       = 1      ! Level of profiling information
! !!! Note: commenting out use of time_out and unit_log means
!      !%%%     commented out for HSL 2011 using !%%%

!      integer :: min_width_blas  = 8      ! Minimum width of source block
!      ! before we use an indirect update_between

!   end type spllt_control

contains

  !!!!!!!!!!!!!!!!!!!!!
  ! TODO update for STARPU
  !
  subroutine spllt_deallocate_node(node, stat)
    type(spllt_node), intent(inout) :: node
    integer,          intent(out)   :: stat

    integer :: st

    stat = 0
    if(allocated(node%index)) then
      deallocate(node%index, stat=st)
      stat = stat + st
    end if
    if(allocated(node%child)) then
      deallocate(node%child, stat=st)
      stat = stat + st
    end if
    if(allocated(node%extra_row)) then
      deallocate(node%extra_row, stat=st)
      stat = stat + st
    end if
  end subroutine spllt_deallocate_node



  !!!!!!!!!!!!!!!!!!!!!
  ! TODO update for non OMP runtime
  !
  subroutine spllt_deallocate_block(blk, stat)
    type(spllt_block), intent(inout)  :: blk
    integer,           intent(out)    :: stat

    integer :: st

    stat = 0
    ! fwd data
    if(allocated(blk%fwd_dep)) then
      deallocate(blk%fwd_dep, stat=st)
      stat = stat + st
    end if
    if(allocated(blk%fwd_update_dep)) then
      deallocate(blk%fwd_update_dep, stat=st)
      stat = stat + st
    end if
    ! bwd data
    if(allocated(blk%bwd_dep)) then
      deallocate(blk%bwd_dep, stat=st)
      stat = stat + st
    end if
    if(allocated(blk%bwd_solve_dep)) then
      deallocate(blk%bwd_solve_dep, stat=st)
      stat = stat + st
    end if
  end subroutine spllt_deallocate_block



  !!!!!!!!!!!!!!!!!!!!!
  ! TODO update for STARPU
  !
  subroutine spllt_deallocate_workspace_i(workspace, stat)
    type(spllt_workspace_i), intent(inout)  :: workspace
    integer,                 intent(out)    :: stat

    if(associated(workspace%c)) then
      deallocate(workspace%c, stat=stat)
    end if
  end subroutine spllt_deallocate_workspace_i



  subroutine spllt_deallocate_lfactor(lfact, stat)
    type(lfactor), intent(inout)  :: lfact
    integer,       intent(out)    :: stat

    if(allocated(lfact%lcol)) then
      deallocate(lfact%lcol, stat=stat)
    end if
  end subroutine spllt_deallocate_lfactor



  subroutine spllt_deallocate_lmap_type(lmap, stat)
    type(lmap_type), intent(inout)  :: lmap
    integer,       intent(out)    :: stat

    if(allocated(lmap%map)) then
      deallocate(lmap%map, stat=stat)
    end if
  end subroutine spllt_deallocate_lmap_type



  !!!!!!!!!!!!!!!!!!!!!
  ! Deallocation of fkeep
  ! TODO update for non OMP runtime
  !
  subroutine spllt_deallocate_fkeep(fkeep, stat)
    type(spllt_fkeep), intent(inout)  :: fkeep
    integer,           intent(out)    :: stat

    integer :: i, st

    stat = 0
    if(allocated(fkeep%bc)) then
      do i = lbound(fkeep%bc,1), ubound(fkeep%bc,1)
        call spllt_deallocate_block(fkeep%bc(i), st)
        stat = stat + st
      end do
      deallocate(fkeep%bc, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%workspace)) then
      do i = lbound(fkeep%workspace,1), ubound(fkeep%workspace,1)
        call spllt_deallocate_block(fkeep%workspace(i), st)
        stat = stat + st
      end do
      deallocate(fkeep%workspace, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%row_list)) then
      do i = lbound(fkeep%row_list,1), ubound(fkeep%row_list,1)
        call spllt_deallocate_workspace_i(fkeep%row_list(i), st)
        stat = stat + st
      end do
      deallocate(fkeep%row_list, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%col_list)) then
      do i = lbound(fkeep%col_list,1), ubound(fkeep%col_list,1)
        call spllt_deallocate_workspace_i(fkeep%col_list(i), st)
        stat = stat + st
      end do
      deallocate(fkeep%col_list, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%map)) then
      do i = lbound(fkeep%map,1), ubound(fkeep%map,1)
        call spllt_deallocate_workspace_i(fkeep%map(i), st)
        stat = stat + st
      end do
      deallocate(fkeep%map, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%flag_array)) then
      deallocate(fkeep%flag_array, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%lfact)) then
      do i = lbound(fkeep%lfact,1), ubound(fkeep%lfact,1)
        call spllt_deallocate_lfactor(fkeep%lfact(i), st)
        stat = stat + st
      end do
      deallocate(fkeep%lfact, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%lmap)) then
      do i = lbound(fkeep%lmap,1), ubound(fkeep%lmap,1)
        call spllt_deallocate_lmap_type(fkeep%lmap(i), st)
        stat = stat + st
      end do
      deallocate(fkeep%lmap, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%small)) then
      deallocate(fkeep%small, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%assoc_tree)) then
      deallocate(fkeep%assoc_tree, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%rhsPtr)) then
      deallocate(fkeep%rhsPtr, stat=st)
      stat = stat + st
    end if
    if(allocated(fkeep%indir_rhs)) then
      deallocate(fkeep%indir_rhs, stat=st)
      stat = stat + st
    end if
  end subroutine spllt_deallocate_fkeep



  subroutine spllt_deallocate_akeep(akeep, stat)
    type(spllt_akeep), intent(inout)  :: akeep
    integer,           intent(out)    :: stat

    integer :: st

    stat = 0
    if(allocated(akeep%weight)) then
      deallocate(akeep%weight, stat=st)
      stat = stat + st
    end if
    if(allocated(akeep%small)) then
      deallocate(akeep%small, stat=st)
      stat = stat + st
    end if
  end subroutine spllt_deallocate_akeep 


  subroutine spllt_solve_workspace_size(fkeep, nworker, nrhs, size)
    implicit none
    type(spllt_fkeep),  intent(in)  :: fkeep
    integer,            intent(in)  :: nworker
    integer,            intent(in)  :: nrhs
    integer(long),      intent(out) :: size

    integer(long) :: n

    n = int(fkeep%n, long)

    size = n * nrhs + (fkeep%maxmn + n) * nrhs * nworker

  end subroutine spllt_solve_workspace_size

  !*************************************************  
  !
  ! Returns the destination block of an internal update task.
  ! Called by add_updates.  
  integer(long) function get_dest_block(src1, src2)

    type(spllt_block), intent(in) :: src1
    type(spllt_block), intent(in) :: src2

    integer(long) :: i
    integer :: sz

    ! Move to diagonal block of target column
    ! sz is the number of (row) blocks in src1
    sz = src1%last_blk - src1%dblk + 1 
    get_dest_block = src1%dblk
    do i = src1%dblk+1, src1%id
       get_dest_block = get_dest_block + sz
       sz = sz - 1
    end do

    ! Move to relevant row block in target col.
    get_dest_block = get_dest_block + src2%id - src1%id

  end function get_dest_block

  !*************************************************  
  !
  subroutine spllt_print_subtree(tree)
    implicit none
    type(spllt_tree_t), intent(in) :: tree

    print *, "Tree id : ", tree%num
    print *, "#node   : ", tree%nnode
    print *, "Node    : [", tree%node_sa, " , ", tree%node_en, "]"
    print *, "#child  : ", tree%nchild

  end subroutine spllt_print_subtree


  
  subroutine spllt_print_trees(fkeep)
    implicit none
    type(spllt_fkeep), target,  intent(inout) :: fkeep

    integer :: i

    do i = 0, size(fkeep%trees)
      call spllt_print_subtree(fkeep%trees(i))
    end do
  end subroutine spllt_print_trees



  subroutine spllt_init_tree(tree, num)
    implicit none
    type(spllt_tree_t), intent(inout) :: tree
    integer,            intent(in)    :: num

    tree%num      = num
    tree%nnode    = 0
    tree%node_sa  = 0
    tree%node_en  = 0
    tree%nchild   = 0
    tree%parent   = 0
  end subroutine spllt_init_tree



  subroutine spllt_create_subtree(akeep, fkeep)
  ! use utils_mod, ONLY : print_iarray
    implicit none
    type(spllt_akeep),          intent(in)    :: akeep
    type(spllt_fkeep), target,  intent(inout) :: fkeep

    integer :: ntree, num
    integer :: i, nnode
    integer :: st
    type(spllt_tree_t), pointer :: p_tree

    allocate(fkeep%small(akeep%nnodes))
    fkeep%small = akeep%small

  ! call print_iarray("Small", size(fkeep%small), fkeep%small, 1)

    ntree = count(akeep%small .eq. 1)

    allocate(fkeep%trees(ntree), stat=st)
    allocate(fkeep%assoc_tree(akeep%nnodes), stat=st)

    print *, "#Subtree : ", ntree

    if(ntree .gt. 0) then

      nnode   = akeep%nnodes
      num     = 1
      p_tree  => fkeep%trees(num)
      call spllt_init_tree(p_tree, num)

      do i = 1, nnode
        if( akeep%small(i) .lt. 0 .or. akeep%small(i) .eq. 1) then
          if( p_tree%nnode .eq. 0) then
            p_tree%node_en  = merge(-akeep%small(i), i, akeep%small(i) .ne. 1)
            p_tree%node_sa  = i
            p_tree%nnode    = p_tree%node_en - p_tree%node_sa + 1
            p_tree%nchild   = 0   ! none yet
            p_tree%parent   = num ! itself
          end if
          if(akeep%small(i) .eq. 1) then
            fkeep%assoc_tree(i) = num
    !       call spllt_print_subtree(p_tree)
            num     = num + 1
            if( num .le. ntree) then
              p_tree  => fkeep%trees(num)
              call spllt_init_tree(p_tree, num)
            end if
          end if
        end if
      end do
    end if

  end subroutine spllt_create_subtree
end module spllt_data_mod


