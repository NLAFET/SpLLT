module spllt_data_mod
#if defined(SPLLT_USE_STARPU)
  use, intrinsic :: iso_c_binding
  use starpu_f_mod
#elif defined(SPLLT_USE_PARSEC)
  use parsec_f08_interfaces
#endif
  use spral_ssids_inform, only: ssids_inform
  implicit none

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

#if defined(SPLLT_USE_OMP) && defined(SPLLT_OMP_TRACE) 
  integer, save :: ini_nde_id, fac_blk_id, slv_blk_id, upd_blk_id, upd_btw_id 
#endif

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
  end type spllt_options

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
     integer :: maxmn ! holds largest block dimension
     integer :: n  ! Order of the system.
     ! type(node_type), dimension(:), allocatable :: nodes ! nodal info
     integer :: nbcol = 0 ! number of block columns in L
     type(lfactor), dimension(:), allocatable :: lfact
     ! holds block cols of L
     type(lmap_type), dimension(:), allocatable :: lmap
     ! holds mapping from matrix values into lfact
  end type spllt_fkeep

  type spllt_omp_task
    integer :: ntask_run
  end type spllt_omp_task

  type spllt_omp_scheduler
    integer :: ntask_insert
    integer :: nfake_task_insert
    integer :: nthread_max
    integer :: nworker
    type(spllt_omp_task), pointer :: info_thread(:)
  end type spllt_omp_scheduler

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

  interface print_array
    module procedure print_darray
    module procedure print_iarray
  end interface print_array

contains

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

!!!!! TMP subroutine !!! TO move from it

  subroutine print_darray(array_name, n, val)
  ! use spllt_data_mod
    character(len=*)                      :: array_name
    integer,                intent(in)    :: n
    real(wp), dimension(n), intent(in)    :: val

    integer :: i
    print '(a)', array_name
    do i = 1, n
     !print '(f20.8)', val(i)
      print *, val(i)
    end do
  end subroutine print_darray

  subroutine print_iarray(array_name, n, val, display)
    character(len=*)                      :: array_name
    integer,                intent(in)    :: n
    integer, dimension(n),  intent(in)    :: val
    integer, optional,      intent(in)    :: display  ! 0 = Vertical,
                                                      ! 1 = Horizontal

    integer :: i
    integer :: disp

    if(present(display))then
      disp = display
    else
      disp = 0 ! Vertical
    end if

    print '(a)', array_name
    if(disp .eq. 0) then
      do i = 1, n
        print '(i9)', val(i)
      end do
    else
      do i = 1, n
        write(*, fmt="(i9)", advance="no") val(i)
      end do
      write(*,*) ""
    end if
  end subroutine print_iarray

  !Return the position in child_node of the rows of node that are 
  ! in child_node
  subroutine get_update_dep(fkeep, child_node, node, pos)
    
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(in)                   :: node
    integer, allocatable, intent(out)     :: pos(:)

    integer :: node_index
    integer :: child_node_index
    integer :: i, j, k
    integer, pointer, dimension(:) :: p_node_index, p_child_node_index

    p_node_index        => fkeep%nodes(node)%index
    p_child_node_index  => fkeep%nodes(child_node)%index


    i = 1
    j = 1
    k = 1
    do while(j .le. size(p_node_index) .and. k .le. size(p_child_node_index))
      if(p_node_index(j) .lt. p_child_node_index(k)) then
        j = j + 1
      else if(p_node_index(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        i = i + 1
        j = j + 1
        k = k + 1
      end if
    end do

    allocate(pos(i - 1))

    if( i .eq. 1 )then
      return
    end if

    !Reiterate to copy the position into pos
    i = 1
    j = 1
    k = 1
    do while(j .le. size(p_node_index) .and. k .le. size(p_child_node_index))
      if(p_node_index(j) .lt. p_child_node_index(k)) then
        j = j + 1
      else if(p_node_index(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        pos(i) = j
        i = i + 1
        j = j + 1
        k = k + 1
      end if
    end do

  end subroutine get_update_dep

  subroutine get_update_nblk(fkeep, child_node, blk_index, nblk)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(inout)                :: blk_index(:)
    integer, intent(out)                  :: nblk

    integer :: child_node_index
    integer :: i, j, k
!   integer :: nb
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp
!   integer :: lblk, nlblk, diff

    p_child_node_index  => fkeep%nodes(child_node)%index
!   nb          = fkeep%nodes(child_node)%nb
    cur_blk_dep = 0 !Impossible value but used as initialization

    nblk  = 0
    j     = 1
    k     = 1
    do while(j .le. size(blk_index) .and. k .le. size(p_child_node_index))
      if(blk_index(j) .lt. p_child_node_index(k)) then
        j = j + 1
      else if(blk_index(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
!       tmp = fkeep%nodes(child_node)%blk_en  -      &
!         ceiling(0.0 + size(p_child_node_index)/nb) +  &
!         ceiling(0.0 + (k - 1) / nb)
!       tmp   = fkeep%nodes(child_node)%blk_sa
!       lblk  = ceiling(0.0 + (k - 1) / nb)
!       nlblk = ceiling(0.0 + (size(p_child_node_index) - 1)/nb)
!       diff  = (fkeep%nodes(child_node)%blk_en - &
!         fkeep%bc(fkeep%nodes(child_node)%blk_en)%dblk)
!       if(lblk .lt. diff) then
!         tmp = tmp + lblk * (nlblk - 1) + 1
!       else
!         tmp = fkeep%bc(fkeep%nodes(child_node)%blk_en)%dblk + diff
!       end if

        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))

!       blk_index(j) = - blk_index(j)

        if(cur_blk_dep .lt. tmp) then
          cur_blk_dep = tmp
          nblk = nblk + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

  end subroutine get_update_nblk

  subroutine get_update_dep_blk(fkeep, child_node, blk_index, pos)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(inout)                :: blk_index(:)
    integer, intent(out)                  :: pos(:)

    integer :: child_node_index
    integer :: i, j, k
    integer :: nb
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp
    integer :: lblk, nlblk, diff

    p_child_node_index  => fkeep%nodes(child_node)%index
    nb = fkeep%nodes(child_node)%nb
    cur_blk_dep = 0 !Impossible value but used as initialization

!   print *, " in child node : ", child_node
    !Set variables
    cur_blk_dep = 0
    i           = 1
    j           = 1
    k           = 1
!   if(child_node .eq. 2) then
!     print *, "Index child : ", p_child_node_index
!     print *, "Ref index   : ", blk_index
!   end if
    do while(j .le. size(blk_index) .and. k .le. size(p_child_node_index))
      if(blk_index(j) .lt. p_child_node_index(k)) then
        j = j + 1
      else if(blk_index(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
!       tmp = fkeep%nodes(child_node)%blk_en  -      &
!         ceiling(0.0 + size(p_child_node_index)/nb) +  &
!         ceiling(0.0 + (k - 1) / nb)
!       pos(i) = fkeep%nodes(child_node)%blk_en  -      &
!         (size(p_child_node_index) - k)/nb
!       print *, "computed blk = ", pos(i)
        
        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))

!       if(tmp .eq. 131) then
!         print *, " j = ", j, " => blk_index(j) = ", blk_index(j), &
!           " , ", p_child_node_index(k)
!       end if
!       blk_index(j) = - blk_index(j)

   !    if(child_node .eq. 2) then
   !      print *, "Found row ", p_child_node_index(k), " in both"
   !      print *, "Size child_index ", size(p_child_node_index),             &
   !        "pos in child_index ", k
!  !      print *, "diff ", size(p_child_node_index) - k
!  !      print *, "offset local_blk ", (size(p_child_node_index) - k)/nb
   !      print *, "Child_node, last_blk ", fkeep%nodes(child_node)%blk_en,   &
   !        "max ncolBlk ", ceiling((size(p_child_node_index) + 0.0)/nb),       &
   !        "current colBlk", ceiling((k + 0.0) / nb)
   !      print *, "=> blk_dep : ", tmp
   !    end if

        if(cur_blk_dep .lt. tmp) then
          cur_blk_dep = tmp
          pos(i) = tmp
          i = i + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

  end subroutine get_update_dep_blk

  function get_child_dep_blk_id(fkeep, child_node, row, nrow)
    integer                               :: get_child_dep_blk_id
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(in)                   :: row
    integer, intent(in)                   :: nrow

    integer :: lblk, nlblk, diff, tmp, nb

    nb    = fkeep%nodes(child_node)%nb
    tmp   = fkeep%nodes(child_node)%blk_sa
    lblk  = ceiling((row + 0.0) / nb)
    nlblk = ceiling((nrow + 0.0)/nb)
    diff  = (fkeep%nodes(child_node)%blk_en - &
      fkeep%bc(fkeep%nodes(child_node)%blk_en)%dblk)
    if(lblk .le. (nlblk - diff)) then
!     print '(a, i2, a, i2, a, i2, a, i2)', "lblk ", lblk, " <= (", nlblk, &
!       " - ", diff, ") = ", (nlblk - diff)
     !tmp = tmp + (lblk - 1) * (nlblk - 1) + 1
      tmp = tmp + (lblk - 1) * ( nlblk + 1 - 0.5 * lblk )
    else
!     print '(a, i2, a, i2, a, i2, a, i2)', "lblk ", lblk, " > (", nlblk, &
!       " - ", diff, ") = ", (nlblk - diff)
      tmp = fkeep%bc(fkeep%nodes(child_node)%blk_en)%dblk + diff
    end if
    
    get_child_dep_blk_id = tmp

  end function get_child_dep_blk_id

  subroutine reduce_ind_and_get_ndep(fkeep, ind, nind, child_node, ndep)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(inout)                :: ind(:)
    integer, intent(inout)                :: nind       ! #elmt in ind
    integer, intent(in)                   :: child_node
    integer, intent(out)                  :: ndep       ! #block dep

    integer :: i, j, k
    integer :: nb, nval
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp
    integer :: lblk, nlblk, diff

    p_child_node_index  => fkeep%nodes(child_node)%index
    nb                  = fkeep%nodes(child_node)%nb
    cur_blk_dep         = 0 !Impossible value but used as initialization
    ndep                = 0

    !Set variables
    cur_blk_dep = 0
    i           = 1
    j           = 1
    k           = 1
    nval        = 0
    do while(j .le. nind .and. k .le. size(p_child_node_index))
      if(ind(j) .lt. p_child_node_index(k)) then
        nval = nval + 1
        ind(nval) = ind(j)
        j = j + 1
      else if(ind(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        
        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))
        
        if(cur_blk_dep .lt. tmp) then
          print *, "Dep with blk ", tmp
          cur_blk_dep = tmp
          ndep = ndep + 1
          i = i + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

    if(j .lt. nind) then
      nval = nval + 1
      ind(nval: nval + nind - j) = ind(j : nind)
      nind = nval + nind - j
    else 
      nind = nval
    end if


  end subroutine reduce_ind_and_get_ndep

  recursive subroutine getUpdateNDep(fkeep, node, ind, nind, ndep, lvl)
    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: node
    integer, intent(inout)        :: ind(:)
    integer, intent(inout)        :: nind
    integer, intent(inout)        :: ndep(:)
    integer, intent(in), optional :: lvl

    integer :: child, dep, nldep, offset, i, nsubind, nchild
    integer, allocatable :: subind(:)
    integer :: rlvl

    if(present(lvl)) then
      rlvl = lvl
    else
      rlvl = 1
    end if

!   print *, "[", rlvl, "] Get #dep of node ", node
!   call print_iarray("Starting Ind =   ", nind, ind, 1)
    
    nchild = fkeep%nodes(node)%nchild

    if(nind .eq. 0) then
      ndep(1) = 0
      return
    end if

    if(nchild .eq. 0) then
      return
    end if

    offset = 1
    allocate(subind(nind))

    do i = 1, nchild
      child = fkeep%nodes(node)%child(i)
      nldep = 0

      subind = ind
      nsubind = nind

!     print '(a, i3)', "Intersect with child node number : ", child
      call reduce_ind_and_get_ndep(fkeep, subind, nsubind, child, nldep)
!     print '(a, i2, a, i2)', "[", rlvl, "] nldep = ", nldep
!     call print_iarray("Updated ind", nsubind, subind, 1)
      ndep(offset) = nldep
      offset = offset + 1

      if(fkeep%nodes(child)%nchild .gt. 0) then
!       print '(a, i2, a, i2, a, i2)', "Call getUpdateNDep on child ",  &
!         child, " with a ndep space ",                                 &
!         offset, " to ", offset + fkeep%nodes(child)%nchild - 1

        call getUpdateNDep(fkeep, child, subind, nsubind,                 &
          ndep(offset : offset + fkeep%nodes(child)%nchild - 1),  &
          rlvl + 1)

!       print '(a, i2, a)', "[", rlvl, "] RAW ndep "
!       print *, ndep
!       print '(a, i2, a, i2, a)', "[", rlvl, "] #dep of child ", child, " are"
!       print *,  ndep(offset : offset + fkeep%nodes(child)%nchild - 1)

        offset = offset + fkeep%nodes(child)%nchild
      end if
    end do

!   print *, "[", rlvl, "] #Dep found ", ndep
    deallocate(subind)

  end subroutine getUpdateNDep

  !Get the dependencies after counting its number
  subroutine reduce_ind_and_get_dep(fkeep, ind, nind, child_node, dep)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(inout)                :: ind(:)
    integer, intent(inout)                :: nind       ! #elmt in ind
    integer, intent(in)                   :: child_node
    integer, intent(out)                  :: dep(:)     ! List of block dep

    integer :: i, j, k
    integer :: nb, nval
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp
    integer :: lblk, nlblk, diff, ndep

    p_child_node_index  => fkeep%nodes(child_node)%index
    nb                  = fkeep%nodes(child_node)%nb
    cur_blk_dep         = 0 !Impossible value but used as initialization
    ndep                = 1

    !Set variables
    cur_blk_dep = 0
    i           = 1
    j           = 1
    k           = 1
    nval        = 0
    do while(j .le. nind .and. k .le. size(p_child_node_index))
      if(ind(j) .lt. p_child_node_index(k)) then
        nval = nval + 1
        ind(nval) = ind(j)
        j = j + 1
      else if(ind(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        
        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))
        
        if(cur_blk_dep .lt. tmp) then
          print *, "Dep with blk ", tmp
          cur_blk_dep = tmp
          dep(ndep) = tmp
          ndep = ndep + 1
          i = i + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

    if(j .lt. nind) then
      nval = nval + 1
      ind(nval: nval + nind - j) = ind(j : nind)
      nind = nval + nind - j
    else 
      nind = nval
    end if


  end subroutine reduce_ind_and_get_dep

  recursive subroutine getUpdateDep(fkeep, node, ind, nind, dep, ndep, lvl)
    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: node
    integer, intent(inout)        :: ind(:)
    integer, intent(inout)        :: nind
    integer, intent(inout)        :: dep(:)
    integer, intent(in)           :: ndep(:)
    integer, intent(in), optional :: lvl

    integer :: child, nldep, offset, i, nsubind, nchild
    integer, allocatable :: subind(:)
    integer :: rlvl

    if(present(lvl)) then
      rlvl = lvl
    else
      rlvl = 1
    end if

    print *, "[", rlvl, "] Get #dep of node ", node
!   call print_iarray("Starting Ind =   ", nind, ind, 1)
    
    nchild = fkeep%nodes(node)%nchild

    if(nind .eq. 0) then
      return
    end if

    if(nchild .eq. 0) then
      return
    end if

    offset = 1
    allocate(subind(nind))

    do i = 1, nchild
      child = fkeep%nodes(node)%child(i)
      nldep = 0

      subind = ind
      nsubind = nind

      if(ndep(i+1) .gt. ndep(i)) then

        print '(a, i3)', "Intersect with child node number : ", child

        call reduce_ind_and_get_dep(fkeep, subind, nsubind, child, &
!         dep(offset : offset + ndep(i+1) -ndep(i) - 1))
          dep(ndep(i) : ndep(i + 1) - 1))
        print '(a, i2, a, i2)', "[", rlvl, "] ldep = "
!       print *, dep(offset : offset + ndep(i+1) - ndep(i) - 1)
        print *, dep(ndep(i) : ndep(i + 1) - 1)
  !     call print_iarray("Updated ind", nsubind, subind, 1)
  !     ndep(offset) = nldep
!       offset = offset + ndep(i)

        if(fkeep%nodes(child)%nchild .gt. 0) then
!         print '(a, i2, a, i2, a, i2)', "Call getUpdateNDep on child ",  &
!           child, " with a ndep space ",                                 &
!           offset, " to ", offset + fkeep%nodes(child)%nchild - 1
          
          nldep = ndep(i + 1 + fkeep%nodes(child)%nchild) - ndep(i + 1)

          if(nldep .gt. 0) then
!           call getUpdateDep(fkeep, child, subind, nsubind,      &
!             dep(offset : offset + nldep - 1),                   &
!             ndep(offset : offset + nldep - 1),                  &
!             rlvl + 1)
            call getUpdateDep(fkeep, child, subind, nsubind,      &
              dep(ndep(i+1) : ndep(i+1) + nldep - 1),             &
              ndep(ndep(i+1) : ndep(i+1) + nldep - 1),            &
              rlvl + 1)

            print '(a, i2, a)', "[", rlvl, "] RAW dep "
            print *, dep
            print '(a, i2, a, i2, a)', "[", rlvl, "] block dep of child ", child, " are"
!           print *,  dep(offset : offset + fkeep%nodes(child)%nchild - 1)
            print *,  dep(ndep(i+1) : ndep(i+1) + nldep - 1)
          end if
!         offset = offset + fkeep%nodes(child)%nchild
        end if
      end if
    end do

    print *, "[", rlvl, "] Blk Dep found ", dep

    deallocate(subind)

  end subroutine getUpdateDep

  subroutine reduce_ind_and_get_dep_blk(fkeep, ind, nind, child_node, blk_dep)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(inout)                :: ind(:)
    integer, intent(inout)                :: nind       ! #elmt in ind
    integer, intent(in)                   :: child_node
    integer, intent(out)                  :: blk_dep(:) ! Dep of the block

    integer :: i, j, k
    integer :: nb, nval
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp
    integer :: lblk, nlblk, diff

    p_child_node_index  => fkeep%nodes(child_node)%index
    nb = fkeep%nodes(child_node)%nb
    cur_blk_dep = 0 !Impossible value but used as initialization
!   print *, "Original index", ind
!   print *, "Filtered by ", p_child_node_index

    !Set variables
    cur_blk_dep = 0
    i           = 1
    j           = 1
    k           = 1
    nval        = 0
    do while(j .le. nind .and. k .le. size(p_child_node_index))
      if(ind(j) .lt. p_child_node_index(k)) then
        nval = nval + 1
        ind(nval) = ind(j)
!       print *, "Copy ", ind(j), "into position ", nval
        j = j + 1
      else if(ind(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        
        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))
        
        if(cur_blk_dep .lt. tmp) then
          cur_blk_dep = tmp
          blk_dep(i) = tmp
          i = i + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

    if(j .lt. nind) then
!     print *," Remain elements to copy ", ind(j : nind)
      nval = nval + 1
      ind(nval: nval + nind - j) = ind(j : nind)
      nind = nval + nind - j
    else 
      nind = nval
    end if

!   call print_iarray("... lead to ", nind, ind)

  end subroutine reduce_ind_and_get_dep_blk

! subroutine getUpdateDep(fkeep, node, ind, nind, dep)
!   type(spllt_fkeep), intent(in) :: fkeep
!   integer, intent(in)           :: node
!   integer, intent(inout)        :: ind
!   integer, intent(inout)        :: nind
!   integer, intent(inout)        :: dep(:)

!   integer :: child
!   integer, allocatable :: subind(:)

!   if(nind .eq. 0) then
!     return
!   end if

!   do i = 1, fkeep%nodes(node)%nchild
!     child = fkeep%nodes(node)%child(i)

!     allocate(subind(nind))
!     subind = ind
!     nsubind = nind

!     call reduce_ind_and_get_dep_blk(fkeep, ind, size(ind), child, dep)
!     call getUpdateDep(fkeep, child, subind, nsubind)

!   end do
! end subroutine getUpdateDep

end module spllt_data_mod


