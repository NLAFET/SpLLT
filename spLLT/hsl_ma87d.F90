! COPYRIGHT (c) 2009 Science and Technology Facilities Council
! Original date 27 April 2009, Version 1.0.0
!
! Written by: Jonathan Hogg, John Reid and Jennifer Scott

! Version 2.3.0
! See ChangeLog for version history

!
! To convert from double:
! * Change wp
! * Change _double
! * Change BLAS calls: dgemv, dgemm, dsyrk, dtrsm, dtrsv
! * Change LAPACK calls: dpotrf
!

module hsl_MA87_double
!$ use omp_lib
   use hsl_mc78_integer
   use hsl_mc34_double
   implicit none

   private
   public :: ma87_keep, ma87_control, ma87_info
   public :: ma87_analyse, ma87_factor, ma87_factor_solve, ma87_solve, &
      ma87_sparse_fwd_solve, ma87_finalise
   public :: ma87_get_n__

   ! factorization kernels
   public :: factorize_posdef, factor_diag_block, solv_col_block, update_block_block, &
        & update_between, update_direct, expand_buffer
   public :: copy_a_to_l, get_dest_block 
   ! derived types
   public :: node_type, block_type, lfactor, make_map, lcol_map

   ! Parameters
   ! Data kinds
   integer, parameter :: wp   = kind(0d0)
   integer, parameter :: long = selected_int_kind(18)

   ! Numerical constants
   real(wp), parameter :: one  = 1.0_wp
   real(wp), parameter :: zero = 0.0_wp

   ! Default values
   integer, parameter :: nemin_default = 32
      ! node amalgamation parameter
   integer, parameter :: nb_default = 256
      ! Block size with dense kernel
   integer, parameter :: pool_default = 25000
      ! size of task pool

   ! Symbolic constants
   ! These flag the different tasks within factor and solve
   integer, parameter :: TASK_DONE             = -1
   integer, parameter :: TASK_NONE             = 0
   integer, parameter :: TASK_FACTORIZE_BLOCK  = 1
   integer, parameter :: TASK_UPDATE_INTERNAL  = 3
   integer, parameter :: TASK_UPDATE_BETWEEN   = 4
   integer, parameter :: TASK_SOLVE_BLOCK      = 5
   integer, parameter :: TASK_SLV_FSLV         = 6
     ! Fwds solve on diag block
   integer, parameter :: TASK_SLV_FUPD         = 7
     ! Fwds update in solve
   integer, parameter :: TASK_SLV_BSLV         = 8
     ! Bwds solve on diag block
   integer, parameter :: TASK_SLV_BUPD         = 9
     ! Bwds update in solve

   ! Types of solve job
   integer, parameter :: SOLVE_JOB_ALL         = 0
   integer, parameter :: SOLVE_JOB_FWD         = 1
   integer, parameter :: SOLVE_JOB_BWD         = 2
   ! How processors share cache                    Example
   integer, parameter :: CACHE_COMPACT       = 1
      ! [0,1], [2,3], [4,5], [6,7]
   integer, parameter :: CACHE_SCATTER       = 2
      ! [0,4]. [1,5], [2,6], [3,7]
   integer, parameter :: CACHE_IDENTITY      = 3
      ! 0, 1, 2, 3, 4, 5, 6, 7

   ! Error flags
   integer, parameter :: MA87_SUCCESS               = 0
   integer, parameter :: MA87_ERROR_ALLOCATION      = -1
   integer, parameter :: MA87_ERROR_ORDER           = -2
   integer, parameter :: MA87_ERROR_NOT_POSDEF      = -3
   integer, parameter :: MA87_ERROR_X_SIZE          = -4
   integer, parameter :: MA87_ERROR_INFINITY        = -5
   integer, parameter :: MA87_ERROR_JOB_OOR         = -6
   integer, parameter :: MA87_ERROR_NBI_OOR         = -7
   integer, parameter :: MA87_ERROR_UNKNOWN         = -99

   ! warning flags
   integer, parameter :: MA87_WARNING_POOL_SMALL    = 1

   !*************************************************

   interface MA87_analyse
      module procedure MA87_analyse_double
   end interface

   interface MA87_factor
      module procedure MA87_factor_double
   end interface

   interface MA87_factor_solve
      module procedure MA87_factor_solve_one_double, &
                       MA87_factor_solve_mult_double
   end interface

   interface MA87_solve
      module procedure MA87_solve_one_double, MA87_solve_mult_double
   end interface

   interface MA87_sparse_fwd_solve
      module procedure MA87_sparse_fwd_solve_double
   end interface

   interface MA87_finalise 
      module procedure MA87_finalise_double
   end interface

   interface ma87_get_n__
      module procedure ma87_get_n_double
   end interface ma87_get_n__

   !*************************************************

   ! Data type for storing information for each block (BLK)
   ! The blocks are numbered 1,2,..., keep%final_blk
   type block_type
      ! Static info, which is set in ma87_analayse
      integer :: bcol            ! block column that blk belongs to
      integer :: blkm            ! height of block (number of rows in blk)
      integer :: blkn            ! width of block (number of columns in blk)
      integer(long) :: dblk      ! id of the block on the diagonal within the 
         ! block column to which blk belongs
      integer :: dep_initial     ! initial dependency count for block,
      integer(long) :: id        ! The block identitifier (ie, its number blk)
      integer(long) :: last_blk  ! id of the last block within the
         ! block column to which blk belongs
      integer :: node            ! node to which blk belongs
      integer :: sa              ! posn of the first entry of the
         ! block blk within the array that holds the block column of L
         ! that blk belongs to.

      ! Non-static info
      integer :: dep  ! dependency countdown/marker. Once factor or solve done,
                      ! value is -2.
!$    integer(omp_lock_kind) :: lock   ! Lock for altering dep
!$    integer(omp_lock_kind) :: alock  ! Lock for altering values in keep%lfact 
         ! for this block.
         ! Note: locks initialised in ma87_analyse and destroyed
         !       in ma87_finalise
   end type block_type

   !*************************************************

   ! Derived type for holding data for each node.
   ! This information is set up by ma87_analyse once the assembly tree
   ! has been constructed.
   type node_type
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

   end type node_type

   !*************************************************

   ! Data type that represents a single block column in L
   ! (allocated by ma87_analyse)
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

   !*************************************************

   ! Data type that contains counts and locks for the solve
   ! (one per a block column)
   type slv_count_type
      integer :: dep
      integer :: dblk
!$    integer(kind=omp_lock_kind) :: lock
   end type slv_count_type

   !*************************************************

   ! Data type for user controls
   type MA87_control

      integer :: diagnostics_level = 0      ! Controls diagnostic printing.
         ! Possible values are:
         !  < 0: no printing.
         !    0: error and warning messages only.
         !    1: as 0 plus basic diagnostic printing.
         !    2: as 1 plus some more detailed diagnostic messages.
         !    3: as 2 plus all entries of user-supplied arrays.
      integer :: nb    = nb_default ! Controls the size of the
         ! blocks used within each node (used to set nb within node_type)
      integer :: nemin = nemin_default    
         ! Node amalgamation parameter. A child node is merged with its parent 
         ! if they both involve fewer than nemin eliminations.
      integer :: pool_size       = pool_default ! Size of task pool arrays
      integer :: unit_diagnostics = 6    ! unit for diagnostic messages
         ! Printing is suppressed if unit_diagnostics  <  0.
      integer :: unit_error       = 6    ! unit for error messages
         ! Printing is suppressed if unit_error  <  0.
      integer :: unit_warning     = 6    ! unit for warning messages
         ! Printing is suppressed if unit_warning  <  0.


      !!!! Undocumented
   !**   integer :: time_out        = -1     ! If >= 0 some internal timings
   !**      are printed on unit time_out. For HSL 2011 these are commented
   !**      using comments !** so easy to search and uncomment
   !%%%  integer :: unit_log        = -1     ! For profiling log output
   !%%%     commented out for HSL 2011 using !%%%
   !%%%  integer :: log_level       = 1      ! Level of profiling information
   !!! Note: commenting out use of time_out and unit_log means
   !%%%     commented out for HSL 2011 using !%%%

   !!! that some subroutines have unused dummy arguments that
   !!! give warnings at compile time. We have not removed them
   !!! since uncommenting the above controls would then be more tedious.


      integer :: cache_tq_sz     = 100    ! Size of local task stack
      integer :: cache_layout    = CACHE_COMPACT ! Proc <=> cache mapping
      integer :: cache_cores     = 2      ! Number of cores per cache
      integer :: min_width_blas  = 8      ! Minimum width of source block
         ! before we use an indirect update_between
    
   end type MA87_control

   !*************************************************

   ! data type for returning information to user.
   type MA87_info 
      real(wp) :: detlog = 0            ! Holds logarithm of abs det A (or 0)
      integer :: flag = 0               ! Error return flag (0 on success)
      integer :: maxdepth = 0           ! Maximum depth of the tree.
      integer(long) :: num_factor = 0_long ! Number of entries in the factor.
      integer(long) :: num_flops = 0_long  ! Number of flops for factor.
      integer :: num_nodes = 0          ! Number of nodes
      integer :: pool_size = pool_default  ! Maximum size of task pool used
      integer :: stat = 0               ! STAT value on error return -1.
   end type MA87_info

   !*************************************************

   ! Data type for communication between threads and routines
   type ma87_keep
 !     private
      type(block_type), dimension(:), allocatable :: blocks ! block info
      integer, dimension(:), allocatable :: flag_array ! allocated to
         ! have size equal to the number of threads. For each thread, holds
         ! error flag
      integer(long) :: final_blk = 0 ! Number of blocks. Used for destroying
         ! locks in finalise
      type(ma87_info) :: info ! Holds copy of info
      integer :: maxmn ! holds largest block dimension
      integer :: n  ! Order of the system.
      type(node_type), dimension(:), allocatable :: nodes ! nodal info
      integer :: nbcol = 0 ! number of block columns in L
      type(lfactor), dimension(:), allocatable :: lfact
         ! holds block cols of L
      type(lmap_type), dimension(:), allocatable :: lmap
         ! holds mapping from matrix values into lfact
   end type ma87_keep

   !*************************************************

   ! Data type for a task
   type dagtask
      integer :: task_type    ! One of TASK_FACTORIZE_BLOCK, ...
      integer(long) :: dest   ! id of the target (destination) block
      integer(long) :: src1   ! 
         ! if task_type = TASK_UPDATE_INTERNAL, src1 holds the id of the first 
         ! source block
         ! if task_type = TASK_UPDATE_BETWEEN, src1 holds the id of a block 
         ! in the block column of the source node that is used
         ! in updating dest.
      integer(long) :: src2   
         ! if task_type = TASK_UPDATE_INTERNAL, src2 holds the id of the second 
         ! source block
         ! (src1 and src2 are blocks belonging to the same block column
         ! of the source node with src1 .le. src2)
         ! src2 is not used by the other tasks
      integer :: csrc(2)
      integer :: rsrc(2)
         ! for an UPDATE_BETWEEN task, we need to hold some additional
         ! information, which locates the source blocks rsrc and csrc
         ! within the source block col.
         ! This info. is set up the subroutine add_between_updates
   end type dagtask

   !*************************************************

   ! Data type for storing tasks we need to do.
   ! Task pool is held as a collection of 4 stacks for different priorities 
   ! using the same space. 
   ! Each stack is a linked list with its head given by an element of prihead.
   ! There are also local stacks for each cache. 

   type taskstack
      integer :: max_pool_size = 0 ! max. number of tasks that are in
         ! the task pool at any one time during the factorization.
      logical :: abort = .false.   ! true if we have aborted
      integer :: active ! Number of active threads 
         ! (number of tasks in execution)
      type(dagtask), dimension(:,:), allocatable :: ctasks ! local task stacks.
         ! allocated to have size (control%cache_tq_sz, ncache), where
         ! ncache is number of caches
      integer, dimension(:), allocatable :: cheads   ! Heads for local stacks.
         ! allocated to have size equal to number of caches
!$    integer(omp_lock_kind), dimension(:), allocatable :: clocks 
         ! Locks for local stacks.
      integer :: freehead  ! Holds the head of linked list of
         ! entries in the task pool that are free
!$    integer(omp_lock_kind) :: lock   ! lock so only one thread at a time
         ! can read/alter the task pool 
      integer :: lowest_priority_value = huge(0) ! 
         ! lowest priority value of the tasks in the pool.
         ! The priority value for each of the different types of task is
         !  1. factor             Highest priority 
         !  2. solve 
         !  3. update_internal 
         !  4. update_between     Lowest priority
      integer, dimension(:), allocatable :: next  ! next task in linked list.
         ! allocated to have size pool_size. Reallocated if initial setting
         ! for pool_size found to be too small.
      integer :: pool_size   ! sizes of task pool arrays next and tasks. 
         ! Initialised to control%pool_size
      integer :: prihead(4)  ! Holds the heads of the linked lists for tasks
         ! with priority values 1,2,3,4. 
      type(dagtask), dimension(:), allocatable :: tasks ! Holds tasks.
         ! allocated to have size pool_size. Reallocated if initial setting
         ! for pool_size found to be too small.
      integer :: total       ! Total number of tasks in pool
   !**   real, dimension(:), allocatable :: waiting  ! Allocated to have size
   !**  ! equal to the number of threads. Used to hold times the threads spent 
   !**  ! waiting if control%time_out >= 0

   end type taskstack
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Analyse phase.
! The user inputs the pivot order and lower
! triangular parts of A. Structure is expanded.
! Supervariables are computed
! and then the assembly tree is constructed and the data structures
! required by the factorization are set up.
! There is no checking of the user's data.
!
subroutine MA87_analyse_double(n, ptr, row, order, keep, control, info)
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   integer, intent(inout), dimension(:) :: order
      ! order(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by MA87_factor.
   ! For details of keep, control, info : see derived type descriptions
   type(MA87_keep), intent(out) :: keep
   type(MA87_control), intent(in) :: control
   type(MA87_info), intent(inout) :: info

   ! Local arrays
   integer, allocatable :: amap(:) ! map from user a to reordered a
   integer, allocatable :: aptr(:) ! column pointers of expanded matrix
   integer, allocatable :: arow(:) ! row pointers of expanded matrix
   integer, allocatable :: count(:) ! used for depth
      ! first search of tree to keep track of level we are at
   integer, allocatable :: cnode(:) ! used for depth
      ! first search of tree
   integer, allocatable :: iw(:) ! work array
   integer, allocatable :: lsize(:) ! lsize(ie) holds the 
      ! length of the list of variables associated with node ie. 
   integer, allocatable :: map(:) ! Allocated to have size n.
      ! used in computing dependency counts. For each row k in 
      ! j-th block column of a node, map(k1) is set to j
   integer, allocatable :: perm(:) ! inverse permutation.
      ! perm(i) holds the variable that is i-th in the pivot sequence.
      ! Also used for checking user-supplied permutation.
   integer, allocatable :: roots(:) ! Holds roots of the forest.
   integer, allocatable :: node_map(:) ! used to number nodes 
      ! contiguously after construction of tree is complete

   ! Local scalars.
   integer :: a_nb ! block size of anode
   integer :: anode ! ancestoral node of snode in tree
   integer(long) :: blk ! temporary variable for holding block number
   integer :: blkn ! number of columns within block column
   integer :: cptr ! pointer into rows of snode that correspond
      ! to columns of an ancestor anode
   integer :: cb ! index of block column within node
   integer :: col_used ! used in computing number of cols in block col.
   integer :: ci ! do loop variable. current block col.
   integer(long) :: dblk ! diagonal block within block column
   integer :: en ! holds keep%nodes(snode)%en
   integer :: i ! temporary variable
   integer :: j ! temporary variable
   integer :: jj ! temporary variable
   integer :: jb ! block index in anode
   integer :: k
   integer :: k1 ! temporary variable
   integer :: l_nb ! set to block size of snode (keep%nodes(snode)%nb)
   integer :: mp  ! copy of control%unit_diagnostics
   integer :: ne ! set to ptr(n+1) - 1
   integer :: nemin ! min. number of eliminations (see control%nemin)
   integer :: node ! a node in tree
   integer :: num_nodes ! number of nodes in tree
   integer :: numcol ! number of cols in node (en-sa+1)
   integer :: numrow ! number of rows in a block column
   integer :: row_used ! used in computing number of rows in a block.
   integer :: sa ! holds keep%nodes(snode)%sa
   integer :: size_anode ! size(keep%nodes(anode)%index)
   integer :: st ! stat parameter
   integer :: sz ! number of blocks in a block column of node
   integer :: swidth ! number of block columns in node
  !** integer :: t_start, t_end, t_rate

   type(mc78_control) :: control78
   integer :: par
   integer :: info78
   integer, dimension(:), allocatable :: sptr, sparent, rlist
   integer(long), dimension(:), allocatable :: rptr

   ! Possible error returns:
   !  MA87_ERROR_ALLOCATION   Allocation error
   !  MA87_ERROR_ORDER        Error in order


   ! initialise
   info%flag = 0
   info%num_factor = 0_long
   info%num_flops = 0_long
   info%num_nodes = 0
   info%maxdepth = 0
   info%stat = 0

   keep%n = n
   ne = ptr(n+1) - 1

   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
      write (mp,'(/a)') ' On entry to MA87_analyse:'
      write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
      write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
      write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
      write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
      write (mp,'(a,i15)') ' control%nemin             =  ',control%nemin
      write (mp,'(a,i15)') ' control%nb                =  ',control%nb
      write (mp,'(a,i15)') ' n                         =  ',n
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then
      write (mp,'(a)') ' ptr = '
      write (mp,'(5i15)') ptr(1:n+1)

      write (mp,'(a)') ' row = '
      write (mp,'(5i15)') row(1:ne)

      ! Print out pivot order.
      write (mp,'(a)') ' User-supplied elimination order :'
      i = min(size(order),n)
      write (mp,'(5i15)') order(1:i)

   else if (control%diagnostics_level==2 .and. mp>=0) then

      write (mp,'(a)') ' ptr(1:min(5,n+1)) = '
      write (mp,'(5i15)') ptr(1:min(5,n+1))

      write (mp,'(a)') ' row(1:min(5,ne)) =  '
      write (mp,'(5i15)') row(1:min(5,ne))

      i = min(10,n)
      i = min(size(order),i)
      write (mp,'(a,2(/5i12))')  &
         ' User-supplied elimination order :', order(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'

   end if

   ! immediate return if n = 0
   if (n == 0) return

   ! expand the matrix

   ! allocate space for expanded matrix (held in aptr,arow)
   allocate (arow(2*ne-n),aptr(n+3),iw(n+1),amap(ptr(n+1)-1),stat=st)
   if (st /= 0) go to 490

   arow(1:ne) = row(1:ne)
   aptr(1:n+1) = ptr(1:n+1)
   call mc34_expand(n, arow, aptr, iw)

   deallocate(iw,stat=st)

   nemin = control%nemin
   ! Check nemin (a node is merged with its parent if both involve
   ! fewer than nemin eliminations). If out of range, use the default
   if (nemin < 1) nemin = nemin_default

   ! Check the user-supplied array order and set the inverse in perm.
   if (size(order).lt.n) then
      info%flag = MA87_ERROR_ORDER
      call MA87_print_flag(info%flag, control, context='MA87_analyse')
      go to 500
   end if

   deallocate (perm,stat=st)
   allocate (perm(n),stat=st)
   if (st /= 0) go to 490
   perm(:) = 0
   k1 = 0
   do i = 1, n
      jj = order(i)
      if (jj < 1 .or. jj > n) exit
      if (perm(jj) /= 0) exit ! Duplicate found
      perm(jj) = i
   end do
   if (i-1 /= n) then
      info%flag = MA87_ERROR_ORDER
      call MA87_print_flag(info%flag, control, context='MA87_analyse')
      go to 500
   end if

   control78%nemin = nemin
   control78%sort = .true.
   control78%lopt = .true.
   call mc78_analyse(n, aptr, arow, order, num_nodes, &
      sptr, sparent, rptr, rlist, control78, info78, nfact=info%num_factor, &
      nflops=info%num_flops)

   info%num_nodes = num_nodes
   !**************************************
   ! Set up nodal data structures
   ! For each node, hold data in keep%nodes(node) 
   deallocate(keep%nodes, stat=st)

   allocate(keep%nodes(-1:num_nodes),stat=st)
   if (st /= 0) go to 490

   keep%nodes(0)%blk_en = 0
   keep%nodes(1)%blk_sa = 1
   keep%nodes(1)%sa = 1

   ! loop over root nodes
   keep%nodes(:)%nchild = 0
   do node = 1, num_nodes
      keep%nodes(node)%sa = sptr(node)
      keep%nodes(node)%en = sptr(node+1)-1

      par = sparent(node)
      keep%nodes(node)%parent = par
      if(par .le. num_nodes) then
         keep%nodes(par)%nchild = keep%nodes(par)%nchild + 1
      else
         keep%nodes(node)%parent = -1
      endif

      ! determine and record the block size for node
      ! note we are careful in case l_nb**2 overflows (in fact 1+l_nb must
      ! not overflow at the end), and limit the answer to huge(l_nb)/2
      l_nb = control%nb
      if (l_nb < 1) l_nb = nb_default
      ! DEBUG deactivated the following oprimisation to avoid small blocks
      ! l_nb = min(huge(l_nb)/2_long, &
      !    (l_nb**2_long) / min(sptr(node+1)-sptr(node), l_nb) )
      ! l_nb = (l_nb-1) / 8 + 1
      ! l_nb = 8 * l_nb
      keep%nodes(node)%nb = l_nb !

      ! Copy row list into keep%nodes
      allocate(keep%nodes(node)%index(rptr(node+1)-rptr(node)),stat=st)
      if (st /= 0) go to 490
      j = 1
      do i = rptr(node), rptr(node+1)-1
         keep%nodes(node)%index(j) = rlist(i)
         j = j + 1
      end do

      ! Allocate space to store child nodes
      allocate(keep%nodes(node)%child(keep%nodes(node)%nchild), stat=st)
      if(st.ne.0) goto 490

      ! Calculate number j of blocks in node and set
      ! keep%nodes(node)%blk_en
      sz = (rptr(node+1)-rptr(node) - 1) / l_nb + 1
      j = 0
      do i = keep%nodes(node)%sa, keep%nodes(node)%en, l_nb
         j = j + sz
         sz = sz - 1
      end do
      keep%nodes(node)%blk_en = keep%nodes(node-1)%blk_en + j

      ! if node is not the final node, record first block
      ! for the next node (node+1)
      if (node < num_nodes)  &
         keep%nodes(node+1)%blk_sa = keep%nodes(node)%blk_en + 1
   end do

   ! set keep%final_blk to hold total number of blocks.
   keep%final_blk = keep%nodes(num_nodes)%blk_en

   ! Add children to nodes, use sptr as a counter as it has fufilled its purpose
   sptr(:) = 0
   do node = 1, num_nodes
      par = sparent(node)
      if(par.gt.num_nodes) cycle
      sptr(par) = sptr(par) + 1
      keep%nodes(par)%child(sptr(par)) = node
   end do

   ! Setup least descendants, to allow easy walk of subtrees
   do node = -1, num_nodes
      ! initialise least descendat to self
      keep%nodes(node)%least_desc = node
   end do
   do node = 1, num_nodes
      ! walk up tree from leaves. A parent's least descendant is either this
      ! nodes least descendant (itself in case of a leaf), or its existing
      ! one if that is smaller.
      anode = keep%nodes(node)%parent
      keep%nodes(anode)%least_desc = &
         min(keep%nodes(node)%least_desc, keep%nodes(anode)%least_desc)
   end do

   !**************************************   
   ! Fill out block information. 

  !** call system_clock(t_start)

   deallocate(keep%blocks,stat=st)
   allocate(keep%blocks(keep%final_blk),stat=st)
   if(st.ne.0) go to 490

   ! Loop over the nodes. Number the blocks in the first node
   ! contiguously, then those in the second node, and so on.
   ! Each node has a number of block columns; the blocks within
   ! each block column are numbered contiguously.
   ! Also add up the number of block columns and store largest block dimension.
   blk = 1
   keep%nbcol = 0
   keep%maxmn = 0
   do node = 1, num_nodes

      sa = keep%nodes(node)%sa
      en = keep%nodes(node)%en
      numcol = en - sa + 1
      numrow = size(keep%nodes(node)%index)

      ! l_nb is the size of the blocks
      l_nb = keep%nodes(node)%nb

      ! sz is number of blocks in the current block column
      sz = (numrow - 1) / l_nb + 1

      ! cb is the index of the block col. within node
      cb = 0
      col_used = 0

      ! Loop over the block columns in node. 
      do ci = sa, en, l_nb
         k = 1 ! use k to hold position of block within block column
         ! increment count of block columns
         keep%nbcol = keep%nbcol + 1

         cb = cb + 1

         ! blkn is the number of columns in the block column.
         ! For all but the last block col. blkn = l_nb.
         blkn = min(l_nb, numcol-col_used)
         col_used = col_used + blkn

         dblk = blk

         ! loop over the row blocks (that is, loop over blocks in block col)
         row_used = 0 
         do blk = dblk, dblk+sz-1
            ! store identity of block
            keep%blocks(blk)%id       = blk

            ! store number of rows in the block.
            ! For all but the last block, the number of rows is l_nb
            keep%blocks(blk)%blkm     = min(l_nb, numrow-row_used)
            row_used = row_used + keep%blocks(blk)%blkm

            ! store number of columns in the block.
            keep%blocks(blk)%blkn     = blkn

            keep%maxmn = max(keep%maxmn, &
                             keep%blocks(blk)%blkm,  keep%blocks(blk)%blkn)

            ! store position of the first entry of the block within the
            ! block column of L
            keep%blocks(blk)%sa       = k

            ! store identity of diagonal block within current block column
            keep%blocks(blk)%dblk     = dblk

            ! store identity of last block within current block column
            keep%blocks(blk)%last_blk = dblk + sz - 1

            ! store node the blk belongs to
            keep%blocks(blk)%node     = node

            ! initialise dependency count
            keep%blocks(blk)%dep_initial = cb

            ! store identity of block column that blk belongs to
            keep%blocks(blk)%bcol     = keep%nbcol

!$          call omp_init_lock(keep%blocks(blk)%lock)
!$          call omp_init_lock(keep%blocks(blk)%alock)

            ! increment k by number of entries in block
            k = k + keep%blocks(blk)%blkm * keep%blocks(blk)%blkn

         end do

         ! Diagonal block has no dependency for factor(dblk)
         keep%blocks(dblk)%dep_initial = cb - 1 

         ! decrement number of row blocks and rows in next block column
         sz = sz - 1
         numrow = numrow - l_nb
      end do
   end do

 !**  if(control%time_out.ge.0) then
 !**     call system_clock(t_end, t_rate)
 !**     write(control%time_out,"(a,es12.4)") "fill block took ", &
 !**     (t_end - t_start) / real(t_rate)
 !**  end if

   !
   ! Compute dependency counts
   ! Note: This might be more efficient if implemented in left-looking
   ! (all descendants) rather than right-looking (all ancestors) fashion.
   !
 !**  call system_clock(t_start)
   allocate (map(n),stat=st)
   if(st.ne.0) go to 490

   ! loop over nodes
   ! FIXME: this loop can be particuarly expensive, and should perhaps be
   ! reordered so the map biulding loop only occurs once for each node.
   ! (eg PARSEC/Ga41As41H72 is particularly bad)
   do node = 1, num_nodes
      ! anode is initially the parent of node. Later it will be the
      ! grandparent, then great grandparent as we move up the tree to the root
      anode = keep%nodes(node)%parent
      ! numcol is number of columns in node
      numcol = keep%nodes(node)%en - keep%nodes(node)%sa + 1

      ! initialise cptr 
      cptr = 1 + numcol

      ! set swidth to number of block columns in node
      l_nb = keep%nodes(node)%nb
      swidth = (numcol-1)/l_nb + 1

      ! loop over ancestors of node
      do while(anode.gt.0)
         ! if we have finished with node, move to next node
         if(cptr.gt.size(keep%nodes(node)%index)) exit

         ! If we have skipped an anode (eg if its only a parent because of 
         ! other nodes in the subtree) we skip the current anode
         if(keep%nodes(node)%index(cptr).gt.keep%nodes(anode)%en) then
            anode = keep%nodes(anode)%parent
            cycle
         endif

         ! Build a map of anode's blocks. 
         ! Within the matrix for anode, the block columns are numbered
         ! 1,2,3... For each row k1 in jb-th block column,
         ! map(k1) is set to jb.

         a_nb = keep%nodes(anode)%nb
         jb = 1 ! Block
         ! loop over the block columns in anode
         size_anode = size(keep%nodes(anode)%index)
         do i = 1, size_anode, a_nb
            ! loop over the rows in the block column
            do k = i, min(i+a_nb-1, size_anode)
               k1 = keep%nodes(anode)%index(k)
               map(k1) = jb
            end do
            jb = jb + 1
         end do
         !print *, "   map = ", map

         ! Loop over affected block columns
         call calc_dep(cptr, node, anode, keep%nodes, keep%blocks, &
            swidth, map)

         ! Move up the tree to the parent of anode
         anode = keep%nodes(anode)%parent 
      end do
   end do

   ! Note: following two functions could probably be combined with mc34 call
   ! above, but have been left as is to ease maintenance
   allocate(keep%lmap(keep%nbcol),stat=st)
   if(st.ne.0) go to 490
   call make_map(n, order, ptr, row, aptr, arow, amap)
   call lcol_map(aptr, arow, num_nodes, keep%nodes, keep%blocks, &
      keep%lmap, map, amap, st)
   if(st.ne.0) goto 490

 !**  if(control%time_out.ge.0) then
 !**     call system_clock(t_end)
 !**     write(control%time_out,"(a,es12.4)") &
 !**        "calculating initial dependencies took ", &
 !**         (t_end - t_start) / real(t_rate)
 !**  end if

   if (mp < 0) go to 500

   if (control%diagnostics_level >= 1) then
      write (mp,'(/a)') ' Leaving MA87_analyse with:'
      write (mp,'(a,i15)')    ' flag              = ',info%flag
      write (mp,'(a,i15)')    ' maxdepth          = ',info%maxdepth
      write (mp,'(a,es15.5)') ' num_factor        = ',real(info%num_factor)
      write (mp,'(a,i15)')    ' num_nodes         = ',info%num_nodes
      write (mp,'(a,i15)')    ' stat              = ',info%stat
   end if
   if (control%diagnostics_level>2) then
      ! Print out pivot order.
      write (mp,'(a)') ' On exit, elimination order :'
      write (mp,'(5i15)') order(1:n)
   else if (control%diagnostics_level==2) then
      i = min(10,n)
      write (mp,'(a,2(/5i12))')  &
         ' On exit, elimination order :', order(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'
   end if

   go to 500

   490 info%flag = MA87_ERROR_ALLOCATION
       info%stat = st
       call MA87_print_flag(info%flag, control, context='MA87_analyse',st=st)

   500 continue
   ! before returning take copy of components of info set by MA87_analyse
   keep%info%flag         = info%flag
   keep%info%num_factor   = info%num_factor
   keep%info%num_flops    = info%num_flops
   keep%info%num_nodes    = info%num_nodes
   keep%info%maxdepth     = info%maxdepth
   keep%info%stat         = info%stat

   deallocate (arow,stat=st)
   deallocate (aptr,stat=st)
   deallocate (amap,stat=st)
   deallocate (count,stat=st)
   deallocate (cnode,stat=st)
   deallocate (iw,stat=st)
   deallocate (lsize,stat=st)
   deallocate (map,stat=st)
   deallocate (perm,stat=st)
   deallocate (roots,stat=st)
   deallocate (node_map,stat=st)

end subroutine MA87_analyse_double

! Make a map from original A to reordered half matrix A
! The reordered half matrix's pattern is returned in nptr and nrow
subroutine make_map(n, perm, optr, orow, nptr, nrow, map)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n+1), intent(in) :: optr
   integer, dimension(optr(n+1)-1), intent(in) :: orow
   integer, dimension(n+3), intent(out) :: nptr ! extra space used for tricks
   integer, dimension(optr(n+1)-1), intent(out) :: nrow
   integer, dimension(optr(n+1)-1), intent(out) :: map

   integer :: i, k, l
   integer(long) :: j

   nptr(:) = 0

   ! Count number of entries in each column of new matrix (at offset 2)
   do i = 1, n
      l = perm(i)
      do j = optr(i), optr(i+1)-1
         k = perm(orow(j))
         if(k<l) then
            nptr(k+2) = nptr(k+2) + 1
         else
            nptr(l+2) = nptr(l+2) + 1
         endif
      end do
   end do

   ! Calculate column starts (at offset 1)
   nptr(1:2) = 1
   do i = 2, n
      nptr(i+1) = nptr(i) + nptr(i+1)
   end do

   ! Now build map
   do i = 1, n
      l = perm(i)
      do j = optr(i), optr(i+1)-1
         k = perm(orow(j))
         if(k<l) then
            map(nptr(k+1)) = j
            nrow(nptr(k+1)) = l
            nptr(k+1) = nptr(k+1) + 1
         else
            map(nptr(l+1)) = j
            nrow(nptr(l+1)) = k
            nptr(l+1) = nptr(l+1) + 1
         endif
      end do
   end do
end subroutine make_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build mapping on per block column basis from user's val to block col's lcol
! This routine uses the reordered half matrix and map from the make_map routine
subroutine lcol_map(aptr, arow, num_nodes, nodes, blocks, lmap, map, amap, st)
   ! Reordered lower triangle of reordered matrix held using aptr and arow
   integer, dimension(:), intent(in) :: aptr
   integer, dimension(:), intent(in) :: arow
   integer, intent(in) :: num_nodes
   type(node_type), dimension(-1:), intent(in) :: nodes ! Node info
   type(block_type), dimension(:), intent(in) :: blocks ! block info
   type(lmap_type), dimension(:), intent(out) :: lmap ! output lcol map
   integer, dimension(:), intent(out) :: map ! work array
   integer, dimension(:), intent(in) :: amap ! map set up by make_map
   integer, intent(out) :: st

   ! Local scalars
   integer :: bcol ! block column
   integer :: cb ! Temporary variable
   integer :: col ! do loop index
   integer(long) :: dblk ! set to keep%nodes(snode)%blk_sa (first block
     ! in snode which is, of course, a diagonal block)
   integer :: en ! set keep%nodes(snode)%en
   integer(long) :: i ! Temporary variable. global row index
   integer(long) :: j ! Temporary variable
   integer :: l_nb ! set to keep%nodes(snode)%nb
   integer(long) :: offset
   integer :: sa ! set keep%nodes(snode)%sa
   integer :: snode ! node
   integer :: swidth ! set to keep%blocks(dblk)%blkn (number of columns
     ! in block column to which dblk belongs)

   integer :: k

   st = 0

   do snode = 1, num_nodes ! loop over nodes
      ! Build a map from global to local indices
      do j = 1, size(nodes(snode)%index)
        i = nodes(snode)%index(j)
        map(i) = j - 1
      end do

      ! Fill in lfact by block columns
      dblk = nodes(snode)%blk_sa

      l_nb = nodes(snode)%nb
      sa = nodes(snode)%sa
      en = nodes(snode)%en

      do cb = sa, en, l_nb
         bcol = blocks(dblk)%bcol

         offset = blocks(dblk)%sa - (cb-sa)*blocks(dblk)%blkn
         swidth = blocks(dblk)%blkn

         k = 1
         lmap(bcol)%len_map = aptr(min(cb+l_nb-1,en)+1) - aptr(cb)
         allocate(lmap(bcol)%map(2,lmap(bcol)%len_map), stat=st)
         if(st.ne.0) return

         ! Loop over columns in the block column bcol
         do col = cb, min(cb+l_nb-1, en)
            ! loop over rows in column col
            do j = aptr(col), aptr(col+1)-1
               i = arow(j)
               i = map(i)
               i = offset + i*swidth
               lmap(bcol)%map(1,k) = i ! destination in lfact(:)%lcol
               lmap(bcol)%map(2,k) = amap(j) ! source
               k = k + 1
            end do
            offset = offset + 1
         end do
         ! move to next block column in snode
         dblk = blocks(dblk)%last_blk + 1
      end do
   end do
end subroutine lcol_map

!****************************************************************************

!
! Factorisation phase.
!
subroutine MA87_factor_double(n, ptr, row, val, order, keep, control, info)

   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   real(wp), intent(in) :: val(:) ! matrix values
   integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
      ! since the analyse phase)
   type(MA87_keep), intent(inout) :: keep ! see description of derived type
   type(MA87_control), intent(in) :: control ! see description of derived type
   type(MA87_info), intent(out) :: info ! see description of derived type

   integer :: i ! temporary variable
   integer :: mp ! copy of control%unit_diagnostics
   integer :: ne ! set to ptr(n+1) - 1
  !** integer :: t_start, t_end, t_rate
   real(wp) :: soln(0)


   ! Reset components of info and return if an
   ! unrecoverable error was encountered earlier
   info%flag = keep%info%flag
   select case (info%flag)
   case(MA87_ERROR_ALLOCATION)
      ! Allocation error on previous subroutine call
      return ! Unrecoverable
   case(MA87_ERROR_ORDER)
      ! Incorrect order fed to previous subroutine call
      return ! Unrecoverable
   case(MA87_ERROR_UNKNOWN)
      ! Unknown error on previous subroutine call
      return ! Unrecoverable
   case default
      ! Otherwise assume recoverable error, reset flag to 0 and continue
      info%flag = 0
   end select
   info%num_factor   = keep%info%num_factor
   info%num_flops    = keep%info%num_flops
   info%num_nodes    = keep%info%num_nodes
   info%maxdepth     = keep%info%maxdepth
   info%stat         = keep%info%stat

   ne = ptr(n+1) - 1
   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
      write (mp,'(/a)') ' On entry to MA87_factor:'
      write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
      write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
      write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
      write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
      write (mp,'(a,i15)') ' control%pool_size         =  ', &
        control%pool_size
      write (mp,'(a,i15)') ' n                         =  ',n
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then

      write (mp,'(a)') ' ptr = '
      write (mp,'(5i15)') ptr(1:n+1)

      write (mp,'(a)') ' row = '
      write (mp,'(5i15)') row(1:ne)

      write (mp,'(a)') ' val = '
      write (mp,'(5es14.6)') val(1:ne)

      write (mp,'(a)') ' Elimination order :'
      write (mp,'(5i15)') order(1:n)

   else if (control%diagnostics_level==2 .and. mp>=0) then

      write (mp,'(a)') ' ptr(1:min(5,n+1)) = '
      write (mp,'(5i15)') ptr(1:min(5,n+1))

      write (mp,'(a)') ' row(1:min(5,ne)) =  '
      write (mp,'(5i15)') row(1:min(5,ne))

      write (mp,'(a)') ' val(1:min(5,ne)) =  '
      write (mp,'(5es14.6)') val(1:min(5,ne))

      i = min(10,n)
      write (mp,'(a,2(/5i12))')  &
         ' Elimination order :', order(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'

   end if

   ! immediate return if n = 0
   if (n == 0) return

   ! Ready to perform the sparse factorization
 !**  call system_clock(t_start)

   call factorize_posdef(n, val, order, keep, control, info, 0, 0, soln)

   if (info%flag < 0) then
      keep%info%flag  = info%flag
      return
   end if

 !**  if(control%time_out.ge.0) then
 !**     call system_clock(t_end, t_rate)
 !**     write(control%time_out,"(a,es12.4)") "factorization took ", &
 !**     (t_end - t_start) / real(t_rate)
 !**  end if

   if (control%diagnostics_level >= 1 .and. mp >= 0) then
      write (mp,'(/a)')       ' Leaving MA87_factor with:'
      write (mp,'(a,i15)')    ' flag              = ',info%flag
      write (mp,'(a,i15)')    ' num_nodes         = ',info%num_nodes
      write (mp,'(a,i15)')    ' num_factor        = ',info%num_factor
      write (mp,'(a,i15)')    ' num_flops         = ',info%num_flops
      write (mp,'(a,i15)')    ' pool_size         = ',info%pool_size
      write (mp,'(a,i15)')    ' stat              = ',info%stat
   end if

   ! Take a copy of any components of info that may have changed
   keep%info%flag          = info%flag
   keep%info%num_nodes     = info%num_nodes
   keep%info%num_factor    = info%num_factor
   keep%info%num_flops     = info%num_flops
   keep%info%pool_size     = info%pool_size

end subroutine MA87_factor_double

!****************************************************************************

!
! Combined factor+solve (simplified interface for single rhs).
!
subroutine MA87_factor_solve_one_double(n, ptr, row, val, order, keep, control,&
      info, x)
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   real(wp), intent(in) :: val(:) ! matrix values
   integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
      ! since the analyse phase)
   type(MA87_keep), intent(inout) :: keep ! see description of derived type
   type(MA87_control), intent(in) :: control ! see description of derived type
   type(MA87_info), intent(out) :: info ! see description of derived type
   real(wp), intent(inout) :: x(keep%n) ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i) is the corresponding component of the right-hand side.
      ! On exit, if i has been used to index a variable,
      ! x(i) holds solution for variable i.

   call MA87_factor_solve_mult_double(n, ptr, row, val, order, keep, &
      control, info, 1, keep%n, x)

end subroutine MA87_factor_solve_one_double

!****************************************************************************

!
! Combined factor+solve.
!
subroutine MA87_factor_solve_mult_double(n, ptr, row, val, order, keep, &
      control, info, nrhs, lx, x)

   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   real(wp), intent(in) :: val(:) ! matrix values
   integer, intent(in) :: order(:) ! holds pivot order (must be unchanged
      ! since the analyse phase)
   type(MA87_keep), intent(inout) :: keep ! see description of derived type
   type(MA87_control), intent(in) :: control ! see description of derived type
   type(MA87_info), intent(out) :: info ! see description of derived type
   integer, intent(in) :: nrhs ! number of right-hand sides to solver for.
   integer, intent(in) :: lx ! first dimension of x
   real(wp), intent(inout) :: x(lx,nrhs) ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i,j) is the corresponding component of the
      ! right-hand side for the jth system (j = 1,2,..., nrhs).
      ! On exit, if i has been used to index a variable,
      ! x(i,j) holds solution for variable i to system j

   integer :: i ! temporary variable
   integer :: j ! temporary variable
   integer :: mp ! copy of control%unit_diagnostics
   integer :: ne ! set to ptr(n+1) - 1
 !**  integer :: t_start, t_end, t_rate
   integer :: st ! stat parameter

   real(wp), dimension(:), allocatable :: soln ! allocated to have 
     ! size n*nrhs.
     ! used to hold reordered rhs and then overwritten by reordered solution.

   ! Reset components of info and return if an error was encountered earlier
   info%flag = keep%info%flag
   select case (info%flag)
   case(MA87_ERROR_ALLOCATION)
      ! Allocation error on previous subroutine call
      return ! Unrecoverable
   case(MA87_ERROR_ORDER)
      ! Incorrect order fed to previous subroutine call
      return ! Unrecoverable
   case(MA87_ERROR_UNKNOWN)
      ! Unknown error on previous subroutine call
      return ! Unrecoverable
   case default
      ! Otherwise assume recoverable error, reset flag to 0 and continue
      info%flag = 0
   end select
   info%num_factor   = keep%info%num_factor
   info%num_flops    = keep%info%num_flops
   info%num_nodes    = keep%info%num_nodes
   info%maxdepth     = keep%info%maxdepth
   info%stat         = keep%info%stat

   ne = ptr(n+1) - 1
   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
      write (mp,'(/a)') ' On entry to MA87_factor:'
      write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
      write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
      write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
      write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
      write (mp,'(a,i15)') ' control%pool_size         =  ', &
        control%pool_size
      write (mp,'(a,i15)') ' n                         =  ',n
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then

      write (mp,'(a)') ' ptr = '
      write (mp,'(5i15)') ptr(1:n+1)

      write (mp,'(a)') ' row = '
      write (mp,'(5i15)') row(1:ne)

      write (mp,'(a)') ' val = '
      write (mp,'(5es14.6)') val(1:ne)

      write (mp,'(a)') ' Elimination order :'
      write (mp,'(5i15)') order(1:n)

   else if (control%diagnostics_level==2 .and. mp>=0) then

      write (mp,'(a)') ' ptr(1:min(5,n+1)) = '
      write (mp,'(5i15)') ptr(1:min(5,n+1))

      write (mp,'(a)') ' row(1:min(5,ne)) =  '
      write (mp,'(5i15)') row(1:min(5,ne))

      write (mp,'(a)') ' val(1:min(5,ne)) =  '
      write (mp,'(5es14.6)') val(1:min(5,ne))

      i = min(10,n)
      write (mp,'(a,2(/5i12))')  &
         ' Elimination order :', order(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'

   end if

   ! immediate return if n = 0
   if (n == 0) return

   if (lx < n .or. nrhs < 1) then
      info%flag = MA87_ERROR_X_SIZE
      call MA87_print_flag(info%flag, control, context='MA87_factor_solve')
      return
   end if

   allocate (soln(n*nrhs),stat=st)

   if (st.ne.0) then
      info%flag = MA87_ERROR_ALLOCATION
      info%stat = st
      call MA87_print_flag(info%flag, control, context='MA87_factor_solve', &
        st=st)
      return
   end if

   ! Reorder rhs

   do i = 1, nrhs
      do j = 1, n
         soln((i-1)*n + order(j)) = x(j, i)
      end do
   end do

   ! Ready to perform the sparse factorization
 !**  call system_clock(t_start)

   call factorize_posdef(n, val, order, keep, control, info, nrhs, n, soln)

   if (info%flag < 0) then
      keep%info%flag  = info%flag
      go to 10
   end if

 !**  if(control%time_out.ge.0) then
 !**     call system_clock(t_end, t_rate)
 !**     write(control%time_out,"(a,es12.4)") "factorization took ", &
 !**     (t_end - t_start) / real(t_rate)
 !**  end if

   ! Perform back substitution
   call solve_posdef(SOLVE_JOB_BWD, nrhs, soln, n, keep, control, info)
   if (info%flag.lt.0) go to 10

   !
   ! Reorder soln
   !
   do i = 1, nrhs
      do j = 1, n
         x(j, i) = soln((i-1)*n + order(j))
      end do
   end do

   if (control%diagnostics_level >= 1 .and. mp >= 0) then
      write (mp,'(/a)')       ' Leaving MA87_factor with:'
      write (mp,'(a,i15)')    ' flag              = ',info%flag
      write (mp,'(a,i15)')    ' num_nodes         = ',info%num_nodes
      write (mp,'(a,i15)')    ' num_factor        = ',info%num_factor
      write (mp,'(a,i15)')    ' num_flops         = ',info%num_flops
      write (mp,'(a,i15)')    ' pool_size         = ',info%pool_size
      write (mp,'(a,i15)')    ' stat              = ',info%stat
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then
      write (mp,'(a)') ' Computed solution for first right-hand side :'
      write (mp,'(5es14.6)') x(1:n,1)
   else if (control%diagnostics_level==2 .and. mp>=0) then
      i = min(10,n)
      write (mp,'(a)') ' Computed solution for first right-hand side :'
      write (mp,'(5es14.6)') x(1:i,1)
      if (i < n) write (mp,'(a)') '  . . . . . .'
   end if

 10 deallocate(soln,stat=st)

   ! Take a copy of any components of info that may have changed
   keep%info%flag          = info%flag
   keep%info%num_nodes     = info%num_nodes
   keep%info%num_factor    = info%num_factor
   keep%info%num_flops     = info%num_flops
   keep%info%pool_size     = info%pool_size
   keep%info%stat          = info%stat

end subroutine MA87_factor_solve_mult_double

!*************************************************

!
! Solve phase. simplified interface for a single rhs
!
subroutine MA87_solve_one_double(x,order,keep,control,info,job)
   type(MA87_keep), intent(inout) :: keep
   real(wp), intent(inout) :: x(keep%n) ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i) is the corresponding component of the right-hand side.
      ! On exit, if i has been used to index a variable,
      ! x(i) holds solution for variable i.
   integer, intent(in) :: order(:) ! pivot order. must be unchanged
   ! For details of keep, control, info : see derived type description
   type(MA87_control), intent(in) :: control
   type(MA87_info), intent(out) :: info
   integer, optional, intent(in) :: job  ! used to indicate whether
      ! partial solution required
      ! job = 0 or absent: complete solve performed
      ! job = 1 : forward eliminations only (PLx = b)
      ! job = 2 : backsubs only ((PL)^Tx = b)

   call MA87_solve_mult_double(1, keep%n, x, order, keep, &
      control, info, job)

end subroutine MA87_solve_one_double

!*************************************************

! Solve phase. Optionally performs only the forward, diagonal or backward sub.

subroutine MA87_solve_mult_double(nrhs,lx,x,order,keep, &
      control,info,job)
   integer, intent(in) :: nrhs ! number of right-hand sides to solver for.
   integer, intent(in) :: lx ! first dimension of x
   real(wp), intent(inout) :: x(lx,nrhs) ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i,j) is the corresponding component of the
      ! right-hand side for the jth system (j = 1,2,..., nrhs).
      ! On exit, if i has been used to index a variable,
      ! x(i,j) holds solution for variable i to system j
   integer, intent(in) :: order(:) ! pivot order. must be unchanged
   ! For details of keep, control, info : see derived type description
   type(MA87_keep), intent(inout) :: keep
   type(MA87_control), intent(in) :: control
   type(MA87_info), intent(out) :: info
   integer, optional, intent(in) :: job  ! used to indicate whether
      ! partial solution required
      ! job = 0 or absent: complete solve performed
      ! job = 1 : forward eliminations only (PLX = B). 
      ! job = 2 : backsubs only ((PL)^TX = B)

   integer :: i
   integer :: j
   integer :: local_job ! set to job or 0 if job not present
   integer :: mp ! set to control%unit_diagnostics
   integer :: n ! order of system
   integer :: st ! stat parameter
   real(wp), dimension(:), allocatable :: soln ! allocated to have size n*nrhs.
     ! used to hold reordered rhs and then overwritten by reorder solution. 

   ! Reset components of info and return if an error was encountered earlier
   info%flag = keep%info%flag
   if (info%flag < 0) return
   info%maxdepth        = keep%info%maxdepth
   info%num_factor      = keep%info%num_factor
   info%num_flops       = keep%info%num_flops
   info%num_nodes       = keep%info%num_nodes
   info%pool_size       = keep%info%pool_size
   info%stat            = keep%info%stat

   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
     write (mp,'(/a)') ' On entry to MA87_solve:'
     write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
     write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
     write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
     write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
     write (mp,'(a,i15)') ' control%pool_size         =  ', &
        control%pool_size
     write (mp,'(a,i15)') ' nrhs                      =  ',nrhs
     write (mp,'(a,i15)') ' lx                        =  ',lx
     if (present(job)) &
     write (mp,'(a,i15)') ' job                       =  ',job
   end if

   local_job = 0
   if (present(job)) then
      select case (job)
      case (SOLVE_JOB_ALL,SOLVE_JOB_FWD,SOLVE_JOB_BWD)
         ! everything OK
      case default
         info%flag = MA87_ERROR_JOB_OOR
         ! note: do not set keep%info%flag so can recall after correction
         call MA87_print_flag(info%flag, control, context='MA87_solve')
         return
      end select
      local_job = job
   end if

   n = keep%n

   ! immediate return if n = 0
   if (n == 0) return

   if (lx < n .or. nrhs < 1) then
      info%flag = MA87_ERROR_X_SIZE
      ! note: do not set keep%info%flag so can recall after correction
      call MA87_print_flag(info%flag, control, context='MA87_solve')
      return
   end if

   !
   ! Reorder rhs
   !
   deallocate(soln,stat=st)
   allocate(soln(n*nrhs),stat=st)
   if (st.ne.0) then
      info%flag = MA87_ERROR_ALLOCATION
      info%stat = st
      call MA87_print_flag(info%flag, control, context='MA87_solve',st=st)
      return
   end if

   do i = 1, nrhs
      do j = 1, n
         soln((i-1)*n + order(j)) = x(j, i)
      end do
   end do

   call solve_posdef(local_job, nrhs, soln, n, keep, control, info)
   if (info%flag.lt.0) go to 10

   !
   ! Reorder soln
   !
   do i = 1, nrhs
      do j = 1, n
         x(j, i) = soln((i-1)*n + order(j))
      end do
   end do

   if (control%diagnostics_level >= 1 .and. mp >= 0) then
      write (mp,'(/a)') ' Leaving MA87_solve with:'
      write (mp,'(a,i15)') ' flag              = ',info%flag
      write (mp,'(a,i15)') ' stat              = ',info%stat
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then
      write (mp,'(a)') ' Computed solution for first right-hand side :'
      write (mp,'(5es14.6)') x(1:n,1)
   else if (control%diagnostics_level==2 .and. mp>=0) then
      i = min(10,n)
      write (mp,'(a)') ' Computed solution for first right-hand side :'
      write (mp,'(5es14.6)') x(1:i,1)
      if (i < n) write (mp,'(a)') '  . . . . . .'
   end if

 10 deallocate(soln,stat=st)

   ! Take a copy of any components of info that may have changed
   keep%info%flag       = info%flag
   keep%info%stat       = info%stat

end subroutine MA87_solve_mult_double

!*************************************************
!
! Sparse forward solve (sparse rhs).
!
subroutine MA87_sparse_fwd_solve_double(nbi,bindex,b,order,invp,nxi,index,x,&
      w,keep,control,info)
   integer, intent(in) :: nbi ! number of nonzero entries in right-hand side
   integer, intent(inout) :: bindex(:) ! On entry, first nbi entries must 
      !  hold indices of  nonzero entries in the right-hand side. 
      !  Only first nbi entries are accessed. Used within code so that 
      !  bindex holds indices of entries in permuted right-hand side 
      !  but reset before return.
   real(wp), intent(in) :: b(:) ! If bindex(i)=k, b(k) must hold the k-th
      ! nonzero component of right-hand side; other entries of b not accessed.
   integer, intent(in) :: order(:) ! pivot order. must be unchanged so that  
      ! order(i)  must hold the  position of variable i in the pivot sequence
   integer, intent(in) :: invp(:) ! must hold inverse pivot order so that 
      ! invp(j) holds the j-th pivot. 
   integer, intent(out) :: nxi ! number of nonzero entries in the solution.
   integer, intent(out) :: index(:) ! First nxi entries holds indices of
      ! nonzero entries in solution. 
   real(wp), intent(inout) :: x(:) ! Must be set by the user to zero on entry.
      ! If index(i)=k, on exit x(k) holds the k-th
      ! nonzero component of solution; all other entries of x are zero.
   real(wp), intent(out) :: w(:) ! work array of size n at least n.
   ! For details of keep, control, info : see derived type description
   type(MA87_keep), intent(inout) :: keep
   type(MA87_control), intent(in) :: control
   type(MA87_info), intent(out) :: info

   integer :: i
   integer :: j
   integer :: k
   integer :: mp ! set to control%unit_diagnostics
   integer :: n ! order of system
   integer :: st ! stat parameter

   ! Reset components of info and return if an error was encountered earlier
   info%flag = keep%info%flag
   if (info%flag < 0) return
   info%maxdepth        = keep%info%maxdepth
   info%num_factor      = keep%info%num_factor
   info%num_flops       = keep%info%num_flops
   info%num_nodes       = keep%info%num_nodes
   info%pool_size       = keep%info%pool_size
   info%stat            = keep%info%stat

   ! Perform appropriate printing
   mp = control%unit_diagnostics
   if (control%diagnostics_level>=1 .and. mp>=0) then
     write (mp,'(/a)') ' On entry to MA87_sparse_fwd_solve:'
     write (mp,'(a,i15)') ' control%diagnostics_level =  ', &
         control%diagnostics_level
     write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
     write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
     write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
     write (mp,'(a,i15)') ' control%pool_size         =  ', &
        control%pool_size
     write (mp,'(a,i15)') ' nbi                       =  ',nbi
   end if

   n = keep%n

   ! immediate return if n = 0
   if (n == 0) return

   ! immediate return if nbi <= 0 .or. nbi > n
   if (nbi <= 0 .or. nbi > n) then
      info%flag = MA87_ERROR_NBI_OOR
      call MA87_print_flag(info%flag, control, context='MA87_sparse_fwd_solve')
      return
   end if

   ! Put non-zero entries of permuted right-hand side into x and
   ! overwrite bindex with permuted bindex 
   do i = 1, nbi
      j = bindex(i)
      x(order(j)) = b(j)
      bindex(i) = order(j)
   end do

   call sparse_fwd_solve_posdef(x, n, nbi, bindex, nxi, index, &
      keep, control, info)
   if (info%flag.lt.0) go to 10

   !
   ! Reorder the solution (but only touch the nxi nonzero entries).
   ! This is where we need a work array w (note: could overwrite b
   ! but then non zero entries of b would be given by index)
   !
   k = 0
   do i = 1, nxi
      j = index(i)
      if (x(j).eq.zero) cycle
      k = k + 1
      index(k) = invp(j)
      w(index(k)) = x(j)
      ! reset to zero
      x(j) = zero
   end do
   nxi = k

   ! copy back into x
   do i = 1,nxi
      j = index(i)
      x(j) = w(j)
   end do

   ! reset bindex
   do i = 1, nbi
      j = bindex(i)
      bindex(i) = invp(j)
   end do

   if (control%diagnostics_level >= 1 .and. mp >= 0) then
      write (mp,'(/a)') ' Leaving MA87_sparse_fwd_solve with:'
      write (mp,'(a,i15)') ' flag              = ',info%flag
      write (mp,'(a,i15)') ' stat              = ',info%stat
      write (mp,'(a,i15)') ' nxi               = ',nxi
   end if

   if (control%diagnostics_level > 2 .and. mp>=0) then
      write (mp,'(a)') ' Index array :'
      write (mp,'(8i8)') index(1:nxi)
      write (mp,'(a)') ' Computed solution :'
      write (mp,'(5es14.6)') x(1:n)
   else if (control%diagnostics_level==2 .and. mp>=0) then
      i = min(10,nxi)
      write (mp,'(a)') ' Index array :'
      write (mp,'(8i8)') index(1:nxi)
      if (i < nxi) write (mp,'(a)') '  . . . . . .'
      i = min(10,n)
      write (mp,'(a)') ' Computed solution :'
      write (mp,'(5es14.6)') x(1:i)
      if (i < n) write (mp,'(a)') '  . . . . . .'
   end if

   10 continue

   ! Take a copy of any components of info that may have changed
   keep%info%flag       = info%flag
   keep%info%stat       = info%stat


end subroutine MA87_sparse_fwd_solve_double

!*************************************************
! This routine must be called after all other calls to routines
! in the package.

subroutine MA87_finalise_double(keep, control)

   type(MA87_keep), intent(inout) :: keep    ! See derived-type declaration
   type(MA87_control), intent(in) :: control ! See derived-type declaration

   integer :: i
   integer :: st     ! stat parameter

   if (control%diagnostics_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(/a)') &
         ' Entering MA87_finalise'
   end if

   if(allocated(keep%lfact)) then
      do i = 1, keep%nbcol
         deallocate(keep%lfact(i)%lcol,stat=st)
      end do
      deallocate(keep%lfact,stat=st)
      keep%nbcol=0    
   endif
   if(allocated(keep%blocks)) then
!$       do i = 1, keep%final_blk
!$          call omp_destroy_lock(keep%blocks(i)%lock)
!$          call omp_destroy_lock(keep%blocks(i)%alock)
!$       end do
      keep%final_blk = 0
      deallocate(keep%blocks)
   endif

   deallocate(keep%nodes,stat=st)
   deallocate(keep%flag_array,stat=st)

end subroutine MA87_finalise_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine factorize_posdef(n, val, order, keep, control, info, nrhs, ldr, rhs)

   integer, intent(in) :: n ! dimension of system
   real(wp), dimension(*), intent(in) :: val
   integer, intent(in) :: order(:)  ! holds pivot order
   type(MA87_keep), intent(inout) :: keep ! see description of derived type
   type(MA87_control), intent(in) :: control ! see description of derived type
   type(MA87_info), intent(inout) :: info ! see description of derived type
   integer, intent(in) :: nrhs  ! number of right-hand sides (maybe = 0)
   integer, intent(in) :: ldr  ! leading extent of rhs
   real(wp), intent(inout) :: rhs(ldr*nrhs)  ! On entry holds rhs data. 
      ! Overwritten by partial solution (forward substitution performed).
   
   ! local derived types
   type(dagtask) :: task ! see description of derived type
   type(taskstack) :: stack ! see description of derived type
   
   ! local arrays
   real(wp), dimension(:), allocatable :: detlog ! per thread sum of log pivot
   integer, dimension(:), allocatable ::  invp ! used to hold inverse ordering
   integer, dimension(:), allocatable ::  map ! allocated to have size n.
     ! used in copying entries of user's matrix a into factor storage 
     ! (keep%fact).
   real(wp), dimension(:,:), allocatable ::  rhs_local ! Local right-hand 
     ! side arrays. allocated to have size (nrhs*ldr,0:total_threads)

   ! local scalars
   integer(long) :: blk ! block identity
   integer :: blkm,blkn ! number of rows/cols in block
   integer(long) :: dblk ! diagonal block within block column
   integer :: en ! holds keep%nodes(snode)%en
   integer :: flag ! Error flag
   integer(long) :: i
   integer :: j
   integer :: l_nb ! set to block size of snode (keep%nodes(snode)%nb)
   integer :: nbcol ! number of block column
   integer :: num_nodes ! number of nodes
   integer :: pool_size ! Normally set to control%pool_size
   integer :: sa ! holds keep%nodes(snode)%sa
   integer :: size_bcol ! number of entries in the block column (sum over the
     ! row blocks in the column)
   integer :: snode 
   integer :: st ! stat parameter
   integer :: sz ! number of blocks in a block column of snode
 !**  integer :: t_start, t_end, t_rate
   integer :: this_thread
   integer :: total_threads ! number of threads being used


   ! Initialise
   flag = 0

   total_threads = 1
!$ total_threads = omp_get_max_threads()

   call zero_task(task)

   num_nodes = keep%info%num_nodes

 !**  call system_clock(t_start)

   ! reset initial dependencies (this may be a second or subsequent
   ! call to factorize phase after single analyse)
   do i = 1,keep%final_blk
      keep%blocks(i)%dep = keep%blocks(i)%dep_initial 
   end do

   ! Initialize task pool
   pool_size = control%pool_size
   if (pool_size < 1) pool_size = pool_default

   call init_stack(stack, pool_size, control, flag, st)

   ! check for allocation error
   if (flag == MA87_ERROR_ALLOCATION) go to 10

   ! Add initial tasks (i.e. those with dependency count equal to 0)
   ! to global task pool to avoid necessity of workstealing
   task%task_type = TASK_FACTORIZE_BLOCK
   do i = 1, keep%final_blk
      if(keep%blocks(i)%dep.ne.0) cycle
      task%dest = i
      call add_task_g(stack, task, control, flag, st)
      ! check for allocation error
      if(flag < 0) go to 10
   end do

 !**  if(control%time_out.ge.0) then
 !**     call system_clock(t_end)
 !**     write(control%time_out,"(a,es12.4)") "init took ", &
 !**        (t_end - t_start) / real(t_rate)
 !**  end if

   ! Set up inverse permutation
   deallocate (invp,stat=st)
   allocate (invp(n),stat=st)
   if(st.ne.0) go to 10

   do j = 1, n
      invp(order(j)) = j
   end do

   ! Allocate factor storage (in keep%lfact)
   deallocate (keep%lfact,stat=st)
   allocate (keep%lfact(keep%nbcol),stat=st)
   if(st.ne.0) go to 10

   blk = 1
   nbcol = 0
   ! loop over the nodes
   do snode = 1, num_nodes
      ! Loop over the block columns in snode, allocating space 
      ! l_nb is the size of the blocks and sz is number of
      ! blocks in the current block column
      l_nb = keep%nodes(snode)%nb
      sz = (size(keep%nodes(snode)%index) - 1) / l_nb + 1
      sa = keep%nodes(snode)%sa
      en = keep%nodes(snode)%en

      size_bcol = 0
      do i = sa, en, l_nb
         nbcol = nbcol + 1
         size_bcol = 0
         dblk = blk
         ! loop over the row blocks
         do blk = dblk, dblk+sz-1
            blkm = keep%blocks(blk)%blkm
            blkn = keep%blocks(blk)%blkn
            size_bcol = size_bcol + blkm*blkn
         end do
         sz = sz - 1
         allocate (keep%lfact(nbcol)%lcol(size_bcol),stat=st)
         if(st.ne.0) go to 10
      end do
   end do

 !**  call system_clock(t_start, t_rate)

   ! Allocate parallel error array
   deallocate(keep%flag_array,stat=st)
   allocate(keep%flag_array(0:total_threads-1),stat=st)
   if(st.ne.0) go to 10

   ! Allocate local right-hand side arrays
   deallocate(rhs_local,stat=st)
   allocate(rhs_local(nrhs*ldr,0:total_threads-1),stat=st)
   if(st.ne.0) go to 10

   ! Allocate information arrays
   deallocate(detlog,stat=st)
   allocate(detlog(0:total_threads-1),stat=st)

10 if(st.ne.0) then
      info%flag = MA87_ERROR_ALLOCATION
      info%stat = st
      call cleanup_stack(stack)
      call MA87_print_flag(info%flag, control, context='MA87_factor',st=st)
      return
   endif

   ! initialise local error flags
   keep%flag_array = 0

   ! initialise rhs_local
   rhs_local(1:nrhs*ldr,0) = rhs(1:nrhs*ldr)
   rhs_local(1:nrhs*ldr,1:total_threads-1) = zero

   !
   ! Copy matrix values across from a into keep%lfact
   !
 !**  call system_clock(t_start)
   st = 0
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(map, flag, this_thread, st) &
!$OMP SHARED(control, detlog, info, invp, keep, ldr, n, nrhs, num_nodes, &
!$OMP    order, rhs, rhs_local, total_threads, stack, val &
!** !$OMP        ,t_start,t_end,t_rate &
!$OMP        )


      this_thread = 0
!$    this_thread = omp_get_thread_num()

      allocate(map(n),stat=st)
      if (st.ne.0) then
         keep%flag_array(this_thread) = MA87_ERROR_ALLOCATION
         info%stat = st
         if(control%diagnostics_level.ge.0 .and. control%unit_error.ge.0) &
            write(control%unit_error, "(/a,i4)") &
            " MA87_factor: Allocation error, stat = ", st
         go to 15
      end if
!$OMP BARRIER
      if(any(keep%flag_array(:).lt.0)) go to 15
      call copy_a_to_l(n,num_nodes,val,map,keep)

!$OMP BARRIER
!$OMP SINGLE
      deallocate(invp,stat=st)

   !**   if(control%time_out.ge.0) then
   !**      call system_clock(t_end,t_rate)
   !**      write(control%time_out,"(a,es12.4)") "copy matrix took ", &
   !**         (t_end - t_start) / real(t_rate)
   !**   end if

   !**   call system_clock(t_start)
!$OMP END SINGLE NOWAIT

      ! Perform actual factorization (and possibly forward subs).
 
      keep%flag_array(this_thread) = info%flag

      call task_dispatch(keep%nbcol, keep%lfact, map, stack, keep%blocks, &
         keep%nodes, control, keep%flag_array(this_thread), st, &
         nrhs, ldr, rhs, total_threads, rhs_local, keep%maxmn, &
         detlog(this_thread))

      ! check for errors
      if(keep%flag_array(this_thread).lt.0) call set_abort(stack)

      if(keep%flag_array(this_thread).eq.MA87_ERROR_ALLOCATION) &
         info%stat = st

      ! Reductions
   15 continue
!$OMP BARRIER
!$OMP SINGLE
         flag = info%flag ! Preserve any warnings from before PARALLEL
         info%flag = minval(keep%flag_array(:))
         if(info%flag.ge.0) &
            info%flag = max(flag, maxval(keep%flag_array(:)))
         info%pool_size = stack%max_pool_size
         info%detlog = sum(detlog(:))
!$OMP END SINGLE
!$OMP END PARALLEL

 !**  if(control%time_out.ge.0) then
 !**     call system_clock(t_end)
 !**     write(control%time_out,"(a,es12.4)") "task_dispatch took ", &
 !**        (t_end - t_start) / real(t_rate)
 !**     write(control%time_out,"(a,50es12.4)")  &
 !**        "Waiting time = ", stack%waiting(:)
 !**  end if

   ! write (6,*) 'max task pool size ', stack%max_pool_size

   call cleanup_stack(stack)

end subroutine factorize_posdef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! this subroutine copies entries in the lower triangular
   ! part of the user-supplied matrix A
   ! into the factor storage, which is keep%lfact.
   ! This makes use of the mapping set up in the analyse phase.
   ! Within each block column, entries are held by rows
   ! (equivalent to holding upper triangular part by columns).

subroutine copy_a_to_l(n,num_nodes,val,map,keep)

   integer, intent(in) :: n      ! order of matrix 
   integer, intent(in) :: num_nodes ! number of nodes in assembly tree
   real(wp), dimension(*), intent(in) :: val ! user's matrix values
   integer, intent(out) :: map(n)  ! mapping array. Reset for each node
     ! so that, if variable (row) i is involved in node,
     ! map(i) is set to its local row index
   type(MA87_keep), intent(inout) :: keep ! on exit, matrix a copied
     ! into relevant part of keep%lfact

   ! Local scalars
   integer :: bcol ! block column
   integer :: cb ! Temporary variable
   integer(long) :: dblk ! set to keep%nodes(snode)%blk_sa (first block
     ! in snode which is, of course, a diagonal block)
   integer :: en ! set keep%nodes(snode)%en
   integer(long) :: i ! Temporary variable. global row index
   integer(long) :: j ! Temporary variable
   integer :: l_nb ! set to keep%nodes(snode)%nb
   integer(long) :: offset
   integer :: sa ! set keep%nodes(snode)%sa
   integer :: snode ! node
   integer :: swidth ! set to keep%blocks(dblk)%blkn (number of columns
     ! in block column to which dblk belongs)
   integer :: sz ! set to size of keep%lfact(nbcol)%lcol

!$OMP DO
   do snode = 1, num_nodes ! loop over nodes
      ! Build a map from global to local indices
      do j = 1, size(keep%nodes(snode)%index)
        i = keep%nodes(snode)%index(j)
        map(i) = j - 1
      end do

      ! Fill in keep%lfact by block columns
      dblk = keep%nodes(snode)%blk_sa

      l_nb = keep%nodes(snode)%nb
      sa = keep%nodes(snode)%sa
      en = keep%nodes(snode)%en

      do cb = sa, en, l_nb
         bcol = keep%blocks(dblk)%bcol
         sz = size(keep%lfact(bcol)%lcol)
         ! Zero the block column. 
         keep%lfact(bcol)%lcol(1:sz) = zero

         offset = keep%blocks(dblk)%sa - (cb-sa)*keep%blocks(dblk)%blkn
         swidth = keep%blocks(dblk)%blkn

         do i = 1, keep%lmap(bcol)%len_map
            keep%lfact(bcol)%lcol(keep%lmap(bcol)%map(1,i)) = &
               val(keep%lmap(bcol)%map(2,i))
         end do

         ! move to next block column in snode
         dblk = keep%blocks(dblk)%last_blk + 1
      end do
   end do
!$OMP END DO

end subroutine copy_a_to_l

!*************************************************

   ! Performs solve. Called by driver MA87_solve.

subroutine solve_posdef(job,  nrhs, rhs, ldr, keep, control, info)

   integer, intent(in) :: job                ! controls full or partial solve
   integer, intent(in) :: nrhs               ! number of rhs
   integer, intent(in) :: ldr                ! leading extent of rhs
   real(wp), intent(inout) :: rhs(ldr*nrhs)  ! On entry holds rhs data. 
      ! Overwritten by solution.
   type(MA87_keep), intent(inout) :: keep
   type(MA87_control), intent(in) :: control
   type(MA87_info), intent(inout) :: info

   integer :: bcol ! blocks(blk)%bcol
   integer(long) :: blk ! block identifier
   type(slv_count_type), dimension(:), allocatable :: counts
   integer :: flag ! temporary store for info%flag
   integer :: i
   integer :: idx ! index into row list
   integer :: j
   integer :: last ! last block column seen
   integer, dimension(:), allocatable :: map ! map of block columns
   integer :: maxmn ! holds largest block dimension
   integer :: node ! temporary variable (for looping over nodes)
   integer, dimension(:), allocatable :: npi ! node previous index
   integer :: num_nodes ! number of nodes
   integer :: offset ! offset into index for current block
   integer :: pool_size ! inital size of task pool
   real(wp), dimension(:,:), allocatable :: rhs_local
   integer :: sa  ! keep%blocks(blk)%sa
   integer :: st
   type(taskstack) :: stack
   type(dagtask) :: task
   integer :: this_thread
   integer :: total_threads
   real(wp), dimension(:), allocatable :: xlocal ! allocated to have
      ! size (maxmn*nrhs).

   total_threads = 1
!$ total_threads = omp_get_max_threads()

   num_nodes = keep%info%num_nodes
   ! Allocate workspace
   maxmn  = keep%maxmn
   deallocate(xlocal,stat=st)
   allocate(xlocal(maxmn*nrhs),stat=st)
   if(st.ne.0) go to 10

   ! Initialize task pool
   pool_size = control%pool_size
   if (pool_size < 1) pool_size = pool_default
   call init_stack(stack, pool_size, control, flag, st)
   ! Note: only possible error is allocation, discard flag as irrelvant.
   if(st.ne.0) go to 10

   ! Allocate parallel error array
   deallocate(keep%flag_array,stat=st)
   allocate(keep%flag_array(0:total_threads-1),stat=st)
   if(st.ne.0) go to 10
   keep%flag_array(:) = 0

   ! Set up map for which blocks variables belong to
   allocate(map(keep%n),stat=st)
   if (st.ne.0) go to 10

   bcol = 1
   do node = 1, num_nodes
      do sa = keep%nodes(node)%sa, keep%nodes(node)%en, keep%nodes(node)%nb
         j = min(sa + keep%nodes(node)%nb - 1, keep%nodes(node)%en)
         map(sa:j) = bcol
         bcol = bcol + 1
      end do
   end do

   call zero_task(task)
   if (job == SOLVE_JOB_ALL .or. job == SOLVE_JOB_FWD) then
      
      ! Setup for Forwards Solve
      ! initialise counts
      allocate(counts(keep%nbcol),stat=st)
      if (st.ne.0) go to 10

      do bcol = 1, keep%nbcol
         counts(bcol)%dep = 0
!$       call omp_init_lock(counts(bcol)%lock)
      end do

      ! Determine dependency counts for each block column
      ! by counting number of blocks that update it.
      bcol = 1
      do node = 1, num_nodes
         blk = keep%nodes(node)%blk_sa
         offset = 1
         do while(blk.le.keep%nodes(node)%blk_en)
            counts(bcol)%dblk = blk
            idx = offset
            do blk = blk, keep%blocks(blk)%last_blk
               last = bcol ! stores last block column updated. Set to bcol
                  ! initially so bcol is ignored for adding dependencies (does
                  ! not depend on itself)
               do idx = idx, idx+keep%blocks(blk)%blkm-1
                  j = map(keep%nodes(node)%index(idx))
                  if(last.ne.j) then
                     ! this block hasn't updated bcol j yet, add a dependency
                     counts(j)%dep = counts(j)%dep + 1
                     last = j
                  endif
               end do
            end do
            bcol = bcol + 1
            offset = offset + keep%nodes(node)%nb
         end do
      end do

      ! Add leaf nodes ready for fwd solve
      task%task_type = TASK_SLV_FSLV
      do node = 1, num_nodes
         blk = keep%nodes(node)%blk_sa

         do while(blk.le.keep%nodes(node)%blk_en)
            bcol = keep%blocks(blk)%bcol
            if(counts(bcol)%dep.eq.0) then
               task%dest = blk
               call add_task(stack, task, control, info%flag, info%stat)
               if(info%flag.lt.0) return
            endif

            blk = keep%blocks(blk)%last_blk+1
         end do
      end do

   else
      ! ensure counts is allocated
      allocate(counts(1),stat=st)
      if (st.ne.0) go to 10

   end if ! job == all or fwd

   if (job == SOLVE_JOB_ALL .or. job == SOLVE_JOB_BWD) then
 
      ! Setup for Backwards Solve

      ! Determine dependency counts for blocks in backwards solve
      ! this being the number of block columns it needs solved before
      ! all its data dependencies are met.
      allocate(npi(num_nodes),stat=st)
      if (st.ne.0) goto 10

      bcol = 1
      call zero_task(task)
      do node = 1, num_nodes
         npi(node) = size(keep%nodes(node)%index)
         blk = keep%nodes(node)%blk_sa
         offset = 1
         do while(blk.le.keep%nodes(node)%blk_en)
            idx = offset
            do blk = blk, keep%blocks(blk)%last_blk
               last = bcol ! stores last block column updated. Set to bcol
                  ! initially so bcol is ignored for adding dependencies (does
                  ! not depend on itself)
               keep%blocks(blk)%dep = 0
               if(blk .eq. keep%blocks(blk)%dblk) then
                  ! diagonal block, account for rest of column updates to
                  ! complete before adding solve task
                  keep%blocks(blk)%dep = keep%blocks(blk)%last_blk - blk
               endif
               do idx = idx, idx+keep%blocks(blk)%blkm-1
                  j = map(keep%nodes(node)%index(idx))
                  if(last.ne.j) then
                     ! we are simply counting the number of different block
                     ! columns that we reach
                     keep%blocks(blk)%dep = keep%blocks(blk)%dep + 1
                     last = j
                  endif
               end do
               if(keep%blocks(blk)%dep.eq.0 .and. job.eq.SOLVE_JOB_BWD) then
                  ! Note: only add task if we are doing a backward solve only
                  task%task_type = TASK_SLV_BSLV
                  task%dest = blk
                  call add_task(stack, task, control, info%flag, info%stat)
                  if(info%flag.lt.0) return
               endif
            end do
            bcol = bcol + 1
            offset = offset + keep%nodes(node)%nb
         end do
      end do

   else
      ! ensure npi is allocated
      allocate(npi(1),stat=st)
      if (st.ne.0) go to 10

   end if ! job == all or bwd

   !
   ! Now execute the establish task graph in parallel
   !
   allocate(rhs_local(ldr*nrhs, 0:total_threads-1),stat=st)

10 if(st.ne.0) then
      info%flag = MA87_ERROR_ALLOCATION
      info%stat = st
      call cleanup_stack(stack)
      call MA87_print_flag(info%flag, control, context='MA87_solve',st=st)
      return
   endif
   rhs_local(:,:) = zero

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(this_thread) &
!$OMP SHARED(control, counts, info, job, keep, ldr, map, maxmn, npi, nrhs, &
!$OMP    rhs, rhs_local, st, stack, total_threads)
   this_thread = 0
!$ this_thread = omp_get_thread_num()

   call solve_task_dispatch(keep%nbcol, keep%lfact, map, npi, stack, &
      keep%blocks, keep%nodes, counts, control, keep%flag_array(this_thread), &
      st, nrhs, ldr, rhs, total_threads, rhs_local, maxmn, job)


   if(keep%flag_array(this_thread).lt.0) call set_abort(stack)

   if(keep%flag_array(this_thread).eq.MA87_ERROR_ALLOCATION) &
      info%stat = st
!$OMP END PARALLEL

   ! Reduce flag_array nicely
   flag = info%flag
   info%flag = minval(keep%flag_array(:))
   if(info%flag.ge.0) &
      info%flag = max(flag, maxval(keep%flag_array(:)))

   !
   ! Dependency error checking - reneable by commenting if statement if
   ! needed for debugging
   !
   !if(.false.) then
   !if(.true.) then
      ! Check for errors in dep counting, fwd solve
   !   if(job.eq.SOLVE_JOB_ALL .or. job.eq.SOLVE_JOB_FWD) then
   !      do i = 1, size(counts)
   !         if(counts(i)%dep.ne.0) then
   !            print *, "fdep(", i, ") = ", counts(i)%dep
   !         endif
   !      end do
   !   endif

      ! Check for any mistakes in dep counting, bwd solve
   !   if(job.eq.SOLVE_JOB_ALL .or. job.eq.SOLVE_JOB_BWD) then
   !      do i = 1, size(keep%blocks)
   !         if(keep%blocks(i)%dep.ne.0) then
   !            print *, "bdep(", i, ") = ", keep%blocks(i)%dep
   !         endif
   !      end do
   !   endif
   !endif

   !
   ! Cleanup things that require explicit stuff
   !
   call cleanup_stack(stack)

   ! Destroy locks
!$ if(job.eq.SOLVE_JOB_ALL .or. job.eq.SOLVE_JOB_FWD) then
!$    do i = 1, size(counts)
!$       call omp_destroy_lock(counts(i)%lock)
!$    end do
!$ endif

end subroutine solve_posdef

!*************************************************
!
! Performs solve. Called by driver MA87_solve.
!
subroutine sparse_fwd_solve_posdef(rhs, n, nb, bindex, nxi, index, &
      keep, control, info)
   integer, intent(in) :: n          
   integer, intent(in) :: nb ! number of nonzer entries in rhs on input 
   real(wp), intent(inout) :: rhs(n)  ! On entry holds rhs data. 
      ! Overwritten by solution.
   integer, intent(in) :: bindex(nb) ! holds indices of nonzero entries 
      ! in rhs
   integer, intent(out) :: nxi ! holds number of number of nonzero entries
      ! in solution.
   integer, intent(out) :: index(n)
   type(MA87_keep), intent(inout) :: keep
   type(MA87_control), intent(in) :: control
   type(MA87_info), intent(inout) :: info

   integer :: bcol ! blocks(blk)%bcol
   integer(long) :: blk ! block identifier
   type(slv_count_type), dimension(:), allocatable :: counts
   integer :: flag ! temporary store for info%flag
   integer :: i
   integer :: idx ! index into row list
   integer :: j
   integer :: last ! last block column seen
   integer :: ldr ! set to n
   integer, dimension(:), allocatable :: map ! map of block columns
   integer :: maxmn ! holds largest block dimension
   integer, dimension(:), allocatable :: npi ! not used 
   integer :: nrhs ! set to 1 (single rhs only)
   integer :: node ! temporary variable (for looping over nodes)
   integer :: num_nodes ! number of nodes
   integer :: offset ! offset into index for current block
   integer :: pool_size ! inital size of task pool
   real(wp), dimension(:,:), allocatable :: rhs_local
   integer :: row
   integer :: sa  ! keep%blocks(blk)%sa
   integer :: st
   type(taskstack) :: stack
   type(dagtask) :: task
   integer :: this_thread
   integer :: total_threads
   real(wp), dimension(:), allocatable :: xlocal ! allocated to have
      ! size (maxmn*nrhs).
   logical, dimension(:), allocatable :: zflag ! allocated to have size nbcol.
      ! initialise to false. Set zflag(j) to true if block col. j of
      ! solution has one or more nonzero entries. 

   total_threads = 1
!$ total_threads = omp_get_max_threads()

   num_nodes = keep%info%num_nodes
   ! Allocate workspace
   maxmn  = keep%maxmn
   deallocate(xlocal,stat=st)
   allocate(xlocal(maxmn),stat=st)
   if(st.ne.0) go to 10

   ! Initialize task pool
   pool_size = control%pool_size
   if (pool_size < 1) pool_size = pool_default
   call init_stack(stack, pool_size, control, flag, st)
   ! Note: only possible error is allocation, discard flag as irrelvant.
   if(st.ne.0) go to 10

   ! Allocate parallel error array
   deallocate(keep%flag_array,stat=st)
   allocate(keep%flag_array(0:total_threads-1),stat=st)
   if(st.ne.0) go to 10
   keep%flag_array(:) = 0

   ! Set up map for which blocks variables belong to
   ! Note: if we are going to do lots of calls to ma87_sparse_fwd_solve, 
   ! we could allocate and set this on first call and then retain.
   ! would need another keep-type structure to keep map unchanged between
   ! calls (or have it as an argument and tell user not to change 
   ! and would have to set
   ! a flag so we could identify if the call is the first after the
   ! factorize phase.
   allocate(map(keep%n),stat=st)
   if (st.ne.0) go to 10

   bcol = 1
   do node = 1, num_nodes
      do sa = keep%nodes(node)%sa, keep%nodes(node)%en, keep%nodes(node)%nb
         j = min(sa + keep%nodes(node)%nb - 1, keep%nodes(node)%en)
         map(sa:j) = bcol
         bcol = bcol + 1
      end do
   end do

   call zero_task(task)
      
      ! initialise counts
      allocate(counts(keep%nbcol),zflag(keep%nbcol),stat=st)
      if (st.ne.0) go to 10

      do bcol = 1, keep%nbcol
         counts(bcol)%dep = 0
!$       call omp_init_lock(counts(bcol)%lock)
      end do

      ! initialise zflag for nonzeros in right-hand side vector
      zflag(1:keep%nbcol) = .false.
      do i = 1,nb
         row = bindex(i)
         bcol = map(row)
         zflag(bcol) = .true.
      end do

      ! Determine dependency counts for each block column
      ! by counting number of blocks that update it.
      bcol = 1
      do node = 1, num_nodes
         blk = keep%nodes(node)%blk_sa
         offset = 1
         do while(blk.le.keep%nodes(node)%blk_en)
            counts(bcol)%dblk = blk
            idx = offset
            if (zflag(bcol)) then
               do blk = blk, keep%blocks(blk)%last_blk
                  last = bcol ! stores last block column updated. Set to bcol
                     ! initially so bcol is ignored for adding dependencies
                     ! (does not depend on itself)
                  do idx = idx, idx+keep%blocks(blk)%blkm-1
                     j = map(keep%nodes(node)%index(idx))
                     if(last.ne.j) then
                        ! this block hasn't updated bcol j yet, add a dependency
                        counts(j)%dep = counts(j)%dep + 1
                        zflag(j) = .true.
                        last = j
                     end if
                  end do
               end do
            else
               blk = keep%blocks(blk)%last_blk + 1
            end if
            bcol = bcol + 1
            offset = offset + keep%nodes(node)%nb
         end do
      end do

      ! Add leaf nodes ready for fwd solve
      task%task_type = TASK_SLV_FSLV
      do node = 1, num_nodes
         blk = keep%nodes(node)%blk_sa

         do while(blk.le.keep%nodes(node)%blk_en)
            bcol = keep%blocks(blk)%bcol
            if(counts(bcol)%dep.eq.0 .and. zflag(bcol)) then
               task%dest = blk
               call add_task(stack, task, control, info%flag, info%stat)
               if(info%flag.lt.0) return
            endif

            blk = keep%blocks(blk)%last_blk+1
         end do
      end do
   !
   ! Now execute the establish task graph in parallel
   !
   allocate(rhs_local(n, 0:total_threads-1),npi(1),stat=st)

10 if(st.ne.0) then
      info%flag = MA87_ERROR_ALLOCATION
      info%stat = st
      call cleanup_stack(stack)
      call MA87_print_flag(info%flag, control, &
          context='MA87_sparse_fwd_solve',st=st)
      return
   endif
   rhs_local(:,:) = zero
   ldr = n
   nrhs = 1

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(this_thread) &
!$OMP SHARED(control, counts, info, keep, ldr, map, maxmn, npi, nrhs, &
!$OMP    rhs, rhs_local, st, stack, total_threads)
   this_thread = 0
!$ this_thread = omp_get_thread_num()

   call solve_task_dispatch(keep%nbcol, keep%lfact, map, npi, stack, &
      keep%blocks, keep%nodes, counts, control, keep%flag_array(this_thread), &
      st, nrhs, ldr, rhs, total_threads, rhs_local, maxmn, SOLVE_JOB_FWD)


   if(keep%flag_array(this_thread).lt.0) call set_abort(stack)

   if(keep%flag_array(this_thread).eq.MA87_ERROR_ALLOCATION) &
      info%stat = st
!$OMP END PARALLEL

   ! Reduce flag_array nicely
   flag = info%flag
   info%flag = minval(keep%flag_array(:))
   if(info%flag.ge.0) &
      info%flag = max(flag, maxval(keep%flag_array(:)))

   ! Cleanup things that require explicit stuff
   !
   call cleanup_stack(stack)

   ! Set nxi and index
   nxi = 0
   bcol = 1
   do node = 1, num_nodes
      ! loop over block columns in node
      do sa = keep%nodes(node)%sa, keep%nodes(node)%en, keep%nodes(node)%nb
         if (zflag(bcol)) then
            ! block is nonzero so nonzero in solution.
            j = min(sa + keep%nodes(node)%nb - 1, keep%nodes(node)%en)
            do i = sa,j
               nxi = nxi + 1
               index(nxi) = i
            end do
         end if
         bcol = bcol + 1
      end do
   end do

end subroutine sparse_fwd_solve_posdef

!*************************************************************
!
! This routine exploits the fact that in any given index list, the nodes
! corresponding to its indices must be solved for in reverse order. We can
! thus store, for each node, the index that we expect to be solved for next.
! This is what npi contains.
!
subroutine bwd_reduce_dep(fnode, bcol, nodes, blocks, map, npi, stack, &
      control, info, st)

   integer, intent(in) :: fnode ! final node to search. i.e. node of bcol
   integer, intent(in) :: bcol ! block column which has now finished
   type(node_type), dimension(-1:), intent(in) :: nodes
   type(block_type), dimension(:), intent(inout) :: blocks
   integer, dimension(:), intent(in) :: map
   integer, dimension(:), intent(inout) :: npi ! node previous index
   type(taskstack), intent(inout) :: stack
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(inout) :: st

   integer(long) :: blk
   integer :: ben
   integer :: bsa
   integer(long) :: dblk
   integer :: j
   integer :: k
   integer :: node
   type(dagtask) :: task

   call zero_task(task)

   !
   ! We have a post order on the tree. Thus we can enumerate a subtree as a
   ! consecutive set of nodes, from the least descendant to the current node.
   ! Thus we avoid checking nodes that can't possibly be relevant
   !
   do node = nodes(fnode)%least_desc, fnode
      ! Check if the next block to check is relevant
      j = npi(node)
      k = nodes(node)%index(j)
      if(map(k).ne.bcol) cycle

      ben = (j-1) / nodes(node)%nb + 1 ! last block row affected

      ! Loop backwards over indices to find the next bcol's position
      ! then store it for later
      do while(j.gt.1 .and. map(k).eq.bcol)
         j = j - 1
         k = nodes(node)%index(j)
      end do
      npi(node) = j

      ! first row affected by bcol is actually j+1
      bsa = (j+1 - 1) / nodes(node)%nb + 1 ! first block row affected

      ! Loop over affected block rows
      do blk = bsa, ben
         ! This block row is a hit, loop over all blocks in it reducing
         ! dependencies
         k = blk - 1
         dblk = nodes(node)%blk_sa
         do j = nodes(node)%sa, nodes(node)%en, nodes(node)%nb
            dblk = dblk + k
            if(k.lt.0) exit ! quit if we try to pass the diagonal
            if(blocks(dblk)%bcol.ge.bcol) exit ! quit if we hit src block
!$          call omp_set_lock(blocks(dblk)%lock)
            blocks(dblk)%dep = blocks(dblk)%dep - 1
            if(blocks(dblk)%dep.eq.0) then
               ! Add relevant task
               task%task_type = TASK_SLV_BUPD
               if(dblk.eq.blocks(dblk)%dblk) &
                  task%task_type = TASK_SLV_BSLV
               task%dest = dblk
               call add_task(stack, task, control, info, st)
            endif
!$          call omp_unset_lock(blocks(dblk)%lock)
            dblk = blocks(dblk)%last_blk + 1
            k = k - 1
         end do ! blk row
      end do ! blks
   end do ! nodes
end subroutine bwd_reduce_dep

!*************************************************************

subroutine fwd_reduce_dep(stack, counts, map, offset, m, index, &
      control, info, st)
   type(taskstack), intent(inout) :: stack
   type(slv_count_type), dimension(:), intent(inout) :: counts
   integer, dimension(:), intent(in) :: map
   integer, intent(in) :: offset
   integer, intent(in) :: m
   integer, dimension(:), intent(in) :: index
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(inout) :: st

   integer :: i
   integer :: last
   type(dagtask) :: task

   call zero_task(task)

   ! Loop over indices, updating dependencies for affected blocks
   last = -1
   do i = offset, offset+m-1
      if(last.ne.map(index(i))) then
         last = map(index(i))
!$       call omp_set_lock(counts(last)%lock)
         counts(last)%dep = counts(last)%dep - 1
         if(counts(last)%dep .eq. 0) then
            task%task_type = TASK_SLV_FSLV
            task%dest = counts(last)%dblk
            call add_task(stack, task, control, info, st)
         endif
!$       call omp_unset_lock(counts(last)%lock)
      endif
   end do
end subroutine fwd_reduce_dep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main task dispatch routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Main worker routine
!
subroutine task_dispatch(nbcol, lfact, map, stack, blocks, nodes, &
      control, info, st, nrhs, ldr, rhs, total_threads, rhs_local, maxmn, &
      detlog)
   integer, intent(in) :: nbcol  ! Size of lfact (no. block cols)
   type(lfactor), dimension(nbcol), intent(inout) :: lfact ! Entries in 
      ! block columns of L
   integer, dimension(:), intent(inout) :: map     ! private work array
   type(taskstack), intent(inout) :: stack         ! task pool
   type(block_type), dimension(:), intent(inout) :: blocks ! block info
   type(node_type), dimension(-1:), target, intent(in) :: nodes ! Node info
   type(MA87_control), intent(in) :: control
   integer, intent(out) :: info ! error flag
   integer, intent(out) :: st ! stat parameter
   integer, intent(in) :: nrhs  ! number of right-hand sides (maybe = 0)
   integer, intent(in) :: ldr  ! leading extent of rhs
   real(wp), intent(inout) :: rhs(ldr*nrhs)  ! On entry holds rhs data. 
      ! Overwritten by partial solution (forward substitution performed).
   integer, intent(in) :: total_threads  ! number of threads
   real(wp), intent(inout) :: rhs_local(ldr*nrhs,0:total_threads-1)  ! Local
      ! right-hand side arrays (one for each thread).
   integer, intent(in) :: maxmn ! max block dimension (=keep%maxmn)
   real(wp), intent(out) :: detlog ! sum of abs values of pivots on this thread 

   integer :: bcol ! block col that task%dest belongs to
   integer :: bcol_src ! block col that task%src1 belongs to
   integer(long) :: blk ! task%dest
   integer :: col ! global index of fist col. in dcol
   real(wp), dimension(:), allocatable :: buffer ! update_buffer workspace
   integer, dimension(:), allocatable :: col_list ! update_buffer workspace
   integer(long) :: dblk ! block on diagonal
   integer :: dcol ! local block column index
   integer :: den ! end column of current block column
   type(node_type), pointer :: dnode
   integer :: dsa ! start column of current block column
   integer :: flag
   type(dagtask) :: task
   integer :: i
   integer(long) :: ii ! do loop variable
   integer(long) :: id ! set to blocks(blk)%id
   integer(long) :: id1 ! set to blocks(dblk)%id
   integer :: j
   integer :: m ! set to blocks(blk)%blkm (number of rows in blk)
   integer :: m1 ! set to blocks(dblk)%blkm (number of rows in dblk)
   integer :: m2 ! set to blocks(task%src2)%blkm 
   integer :: n ! set to blocks(blk)%blkn (number of columns in blk and dblk)
   integer :: n1 ! set to blocks(task%src1)%blkn (=blocks(task%src2)%blkn)
   integer :: offset ! offset into index for current block
   integer, dimension(:), allocatable :: row_list ! update_buffer workspace
   integer :: sa ! set to blocks(blk)%sa 
   integer :: sa1 ! set to blocks(dblk)%sa
   integer :: sa2 ! set to blocks(task%src2)%sa
   integer(long) :: src ! task%src1
  !%%% integer :: t_start, t_end
   integer :: this_thread
   integer :: csrc(2), rsrc(2) ! used for update_between tasks to
     ! locate source blocks


   this_thread = 0
!$ this_thread = omp_get_thread_num()

  ! if(control%diagnostics_level.gt.2 .and. control%unit_diagnostics.ge.0) &
  !    write(control%unit_diagnostics, "(a,i4,a)") &
  !       "Thread ", this_thread, " joined worker pool."

   ! Initialize things
   info = 0; st = 0             ! By default everything went OK
   call zero_task(task)         ! Set everything to zero
   task%task_type = TASK_NONE   ! Needs to be set to prevent 
                                ! dispatched incrementing
   detlog = zero         ! log(abs value of det A) = sum abs values of pivots
                         ! accumlated on a per thread basis

   allocate(col_list(1), row_list(1), buffer(maxmn*nrhs), stat=st)
   if (st.ne.0) then
      info = MA87_ERROR_ALLOCATION
      call MA87_print_flag(info, control, context='MA87_factor',st=st)
      return
   end if

   ! Mark thread as active
!$ call omp_set_lock(stack%lock)
   stack%active = stack%active + 1
!$ call omp_unset_lock(stack%lock)

   !
   ! Main loop
   !
   do
      !
      ! Retrieve next task to perform and wait till all is ready
      !
   !%%%   if(control%unit_log.gt.0) call system_clock(t_start)

      call get_task(stack, task, control, info, st)
      if(info.lt.0) return

  !%%%    if(control%unit_log.gt.0) then
  !%%%       call system_clock(t_end)
  !%%%       call log_task(control, this_thread, t_start, t_end, "GT", &
  !%%%          int(task%task_type,long))
  !%%%    endif

!$OMP FLUSH(lfact)

      !
      !if(control%diagnostics_level.gt.2 .and. control%unit_diagnostics.ge.0 &
      !      .and. task%task_type.ne.TASK_NONE .and. task%type.ne.TASK_DONE)  &
      !   then
      !   write(control%unit_diagnostics, "(i3,2a)") &
      !      this_thread, " got task ",print_job(task) 
      !endif

      !
      ! Perform task and update dependencies
      !
      select case(task%task_type)
      case(TASK_DONE)
         exit

      case(TASK_NONE) ! Job not finished but no tasks available, spin.
         cycle

      case(TASK_FACTORIZE_BLOCK) ! Factorize diagonal block

         blk = task%dest

         ! use n, m, and sa to point to where blk is stored in lfact
         n  = blocks(blk)%blkn
         m  = blocks(blk)%blkm
         sa = blocks(blk)%sa
         id = blocks(blk)%id

         ! bcol is block column that blk belongs to
         bcol = blocks(blk)%bcol
 
         call factor_diag_block(n, m, id, lfact(bcol)%lcol(sa:sa+n*m-1),   &
            control, flag, detlog)

         if(flag.ne.0) then
            if(flag.gt.0) flag = nodes(blocks(blk)%node)%sa + flag - 1
            info = MA87_ERROR_NOT_POSDEF
            if(control%diagnostics_level.ge.0 .and. &
               control%unit_error.ge.0) then
               write (control%unit_error,'(/a,i3)') &
              ' Error return from MA87_factor. Error flag = ', info
               write (control%unit_error,'(a)') &
              ' Matrix is not positive definite.'
               write(control%unit_error, "(/a,i10)") &
                  " non-positive pivot at column ", flag
            end if
            return
         endif

         ! dnode is node to which blk belongs
         dnode => nodes(blocks(blk)%node)

         ! Find block column dcol of dnode that blk belongs to. The block
         ! cols are numbered locally within dnode as 1,2,3,...

         dcol = bcol - blocks(dnode%blk_sa)%bcol + 1
         dsa  = dnode%sa + (dcol-1)*dnode%nb
         den  = min(dsa+dnode%nb-1, dnode%en)

         if(nrhs.gt.0) then
            ! forward substitution
            ! sum contributions from rhs_local
            do j = 1,nrhs
               i = (j-1)*ldr
               do i = i+dsa, i+den
                  rhs(i) = sum(rhs_local(i,:))
               end do
            end do

            call slv_solve(n, n, dsa, lfact(bcol)%lcol(sa:sa+n*n-1), &
               'Transpose    ', 'Non-unit', nrhs, rhs, ldr, control, id)

            ! remains to deal with trapezoidal part (if any)
            m = m - n
            if(m.gt.0) then
               ! hold start location
               offset = 1 + (dcol-1)*dnode%nb + n
               sa = sa + n*n
               call slv_fwd_update(m, n, dsa, offset, dnode%index,    &
                  lfact(bcol)%lcol(sa:sa+n*m-1), n, nrhs,             &
                  rhs_local(:,this_thread), ldr, rhs, ldr, buffer,    &
                  control, id)
            endif
         endif

!$OMP FLUSH

         ! Mark blk as done
         blocks(blk)%dep = -2

         ! Add relevant updates

         call add_between_updates(nodes, blocks, blk, &
            blocks(blk)%node, dcol, stack, map, control, info, st)
         if(info.lt.0) return
         
     !%%%    if(control%unit_log.gt.0 .and. control%log_level.gt.1) &
     !%%%       call system_clock(t_start)

         ! Reduce dependency counts
         do ii = blk+1, blocks(blk)%last_blk
            call reduce_dep(stack, blocks(ii), control, info, st)
            if(info.lt.0) return
         end do

     !%%%    if(control%unit_log.gt.0 .and. control%log_level.gt.1) then
     !%%%       call system_clock(t_end)
     !%%%       call log_task(control, this_thread, t_start, t_end, "FP")
     !%%%    endif

      case(TASK_SOLVE_BLOCK) ! Solve column block with diagonal (_trsm)

         blk = task%dest

         ! use n, m, and sa to point to where blk is stored in lfact
         n  = blocks(blk)%blkn
         m  = blocks(blk)%blkm
         sa = blocks(blk)%sa
         id = blocks(blk)%id
 
         ! dblk is identifier of block on diagonal
         dblk = blocks(blk)%dblk
         m1  = blocks(dblk)%blkm
         sa1 = blocks(dblk)%sa
         id1 = blocks(dblk)%id

         ! bcol is block column that blk and dblk belong to
         bcol = blocks(blk)%bcol

         call solv_col_block(m, n, id, lfact(bcol)%lcol(sa:sa+n*m-1), &
              id1, lfact(bcol)%lcol(sa1:sa1+n*m1), control)

         if(nrhs.gt.0) then
           ! dnode is node to which blk belongs
           dnode => nodes(blocks(blk)%node)

           ! dcol is the local index of the
           ! block col. of dnode that blk belongs to and
           ! col is the global index of the first column of dcol

           dcol = bcol - blocks(dnode%blk_sa)%bcol + 1
           col = dnode%sa + (dcol-1)*dnode%nb
      
           offset = 1 + (dcol-1)*dnode%nb + (blk-dblk)*dnode%nb

           call slv_fwd_update(m, n, col, offset, dnode%index,                &
             lfact(bcol)%lcol(sa:sa+n*m-1), n, nrhs, rhs_local(:,this_thread),&
             ldr, rhs, ldr, buffer, control, id)

         end if

!$OMP FLUSH(lfact, blocks)

     !%%%    if(control%unit_log.gt.0 .and. control%log_level.gt.1)    &
     !%%%       call system_clock(t_start)

         call add_updates(stack, map, nodes, blocks, &
            blk, blocks(blk)%node, control, info, st)
         if(info.lt.0) return

     !%%%    if(control%unit_log.gt.0 .and. control%log_level.gt.1) then
     !%%%       call system_clock(t_end)
     !%%%       call log_task(control, this_thread, t_start, t_end, "SP")
     !%%%    endif


      case(TASK_UPDATE_INTERNAL) ! update one block by two others
         ! in the same node

         blk = task%dest
         n   = blocks(blk)%blkn
         m   = blocks(blk)%blkm
         sa  = blocks(blk)%sa

         ! src1 and src2 belong to same block column
         sa1 = blocks(task%src1)%sa
         n1  = blocks(task%src1)%blkn
         m1  = blocks(task%src1)%blkm

         sa2 = blocks(task%src2)%sa
         m2  = blocks(task%src2)%blkm

         bcol = blocks(blk)%bcol
         bcol_src = blocks(task%src1)%bcol

         call update_block_block(m, n, lfact(bcol)%lcol(sa:sa+n*m-1),  &
            blocks(blk), n1, lfact(bcol_src)%lcol(sa1:sa1+n1*m1-1),    &
            lfact(bcol_src)%lcol(sa2:sa2+n1*m2-1), control)

         ! reduce dependency count for task%dest
         call reduce_dep(stack, blocks(blk), control, info, st)
         if(info.lt.0) return

      case(TASK_UPDATE_BETWEEN) ! update a block with information from 
                                ! another node

         blk = task%dest
         n   = blocks(blk)%blkn
         m   = blocks(blk)%blkm
         sa  = blocks(blk)%sa

         bcol = blocks(blk)%bcol

         ! task%src1 is a block in the source node
         src = task%src1
         n1  = blocks(src)%blkn
         bcol_src = blocks(src)%bcol
         csrc(1:2) = task%csrc(1:2)
         rsrc(1:2) = task%rsrc(1:2)

         call update_between(m, n, blk, nodes(blocks(blk)%node),  &
            n1, task%src1, nodes(blocks(task%src1)%node),         &
            lfact(bcol)%lcol(sa:sa+m*n-1),                        &
            lfact(bcol_src)%lcol(csrc(1):csrc(1)+csrc(2)-1),      &
            lfact(bcol_src)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1),      &
            blocks, row_list, col_list, buffer, control, info, st)
         if(info.lt.0) return

         ! reduce dependency count for task%dest
         call reduce_dep(stack, blocks(blk), control, info, st)
         if(info.lt.0) return

      case default
         info = MA87_ERROR_UNKNOWN
         call MA87_print_flag(info, control, context='MA87_factor')
      end select

   end do

   !if(control%diagnostics_level.gt.2 .and. control%unit_diagnostics.ge.0) &
   !   write(control%unit_diagnostics, "(a,i4,a)") &
   !      "Thread ", this_thread, " left worker pool."

end subroutine task_dispatch

!*************************************************
!
! This routine adds updates and is called after a solve has completed. It
! checks if any other solves have completed, and if it is the final one of
! a pairwise dependency to complete it adds an update task.
! Inter-nodal updates are also added.
!
subroutine add_updates(stack, map, nodes, blocks, blk, node, &
      control, info, st)

   type(taskstack), intent(inout) :: stack ! holds tasks that have to be done.
   integer, dimension(:), intent(inout) :: map ! used for inter-nodal updates.
   type(node_type), dimension(-1:), target, intent(in) :: nodes
   type(block_type), dimension(:), intent(inout) :: blocks
   integer(long), intent(in) :: blk ! block that solve just completed on
   integer, intent(in) :: node ! node to which blk belongs
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag
   integer, intent(inout) :: st ! stat parameter

   integer(long) :: d ! id of diagonal block in block col. that blk belongs to
   integer(long) :: i
   integer(long) :: last ! id of final block that generates an update 
     ! within node
   integer :: j
   integer :: nb ! block size nodes(node)%nb
   integer :: nbcol ! number of block cols in node
   integer :: numcol ! number of cols in node
   integer :: numrow ! number of rows in node
   integer :: sz ! number of blocks in first block col of node
   type(dagtask) :: task
   type(node_type), pointer :: dnode   

   ! Work out useful information on blocks in node
   nb = nodes(node)%nb
   numrow = size(nodes(node)%index)
   numcol = nodes(node)%en - nodes(node)%sa + 1

   !
   ! Determine final block that generates an update within node
   ! (it is the last block in a block row that has diagonal elements)
   !
   ! sz is number of (row) blocks in the first block column of node
   ! nbcol is the number of block cols in node

   sz = (numrow - 1) / nb + 1
   nbcol = (numcol - 1) / nb + 1 
   last = blocks(blk)%last_blk - (sz - nbcol)

   !
   ! Acquire lock
   !
   d = blocks(blk)%dblk
!$ call omp_set_lock(blocks(d)%lock)
   blocks(blk)%dep = -2

   call zero_task(task)
   task%task_type = TASK_UPDATE_INTERNAL

   ! Add updates as they are pairwise valid if reordering has been done
   task%src2 = blocks(blk)%id
   do i = blocks(blk)%dblk+1, min(blocks(blk)%id, last)
      ! Has solve for block i been done?
      if(blocks(i)%dep.ne.-2) cycle 
      task%src1 = i
      task%dest = get_dest_block(blocks(task%src1), blocks(task%src2))
      call add_task(stack, task, control, info, st)
      if(info.lt.0) go to 10
   end do
   if(blocks(blk)%id.le.last) then
      task%src1 = blocks(blk)%id
      do i = blocks(blk)%id+1, blocks(blk)%last_blk
         ! Has solve for block i been done?
         if(blocks(i)%dep.ne.-2) cycle 
         task%src2 = i
         task%dest = get_dest_block(blocks(task%src1), blocks(task%src2))
         call add_task(stack, task, control, info, st)
         if(info.lt.0) go to 10
      end do
   end if

   ! add UPDATE_BETWEEN tasks.

   dnode => nodes(node)

   ! Find block column j of dnode that blk belongs to. The block
   ! cols are numbered locally within dnode as 1,2,3,...

   j = blocks(blk)%bcol - blocks(dnode%blk_sa)%bcol + 1

   call add_between_updates(nodes, blocks, blk, node, &
      j, stack, map, control, info, st)

   ! Release lock
   !
   10 continue
   d = blocks(blk)%dblk
!$ call omp_unset_lock(blocks(d)%lock)

end subroutine add_updates

!*************************************************

! This function calculates column of a node we are on

integer function calc_col(node, block)
   type(node_type), intent(in) :: node
   type(block_type), intent(in) :: block

   calc_col = (size(node%index)-1)/node%nb + 1 ! no. row blks for node

   calc_col = calc_col - (block%last_blk - block%dblk + 1) + 1 ! column of node

end function calc_col

!*************************************************

! Following the completion of the factor or
! solve task for block src in block column scol
! of node snode, add to the task pool those inter-nodal updates to  
! ancestors of snode that are now ready.

subroutine add_between_updates(nodes, blocks, src, snode, scol, &
      stack, map, control, info, st)

   type(node_type), dimension(-1:), intent(in) :: nodes
   type(block_type), dimension(:), intent(inout) :: blocks
   integer(long), intent(in) :: src ! Id of block for which solve just done.
   integer, intent(in) :: snode ! Node containing src.
   integer, intent(in) :: scol  ! Index within snode of block column of src. 
   type(taskstack), intent(inout) :: stack ! holds tasks that have to be done.
   integer, dimension(:), intent(inout) :: map ! Workarray to hold map from row 
     ! indices to block indices in ancestor node. 
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(inout) :: st

   type(dagtask) :: task ! Used to send a task to the pool
   integer :: a_nb  ! Block size of anode
   integer :: anode ! Ancestor of snode
   integer(long) :: blk   ! Block id
   integer :: cb    ! Local index of column block in anode
   integer :: cptr  ! Position in snode of the first row 
     ! matching a column of the current block column of anode.
   integer :: cptr2  ! Position in snode of the last row 
     ! matching a column of the current block column of anode.
   integer(long) :: dblk ! id of diagonal block of anode
   integer :: i
   integer :: ilast
   integer :: jb ! Block index in anode
   integer :: jlast ! Last column in the cb-th block column of anode
   integer :: k
   integer :: k1
   integer :: n ! set to blocks(src)%blkn
   integer :: numcol ! number of cols in snode
   integer :: rb ! Index of block row in snode
   integer(long) :: rblk ! id of block in scol containing row 
      ! nodes(snode)%index(cptr).
   integer :: size_anode ! size(nodes(anode)%index)
   integer :: size_snode ! size(nodes(snode)%index)
   integer :: s_nb   ! Block size of snode
 
   logical :: map_done ! True if map has been built for anode.
   logical :: new ! True if the rows of snode corresponding to the block 
     ! column or block row of anode overlap src.
   logical :: newc ! True if the rows of snode corresponding to the block 
     ! column of anode overlap src. 
   logical :: ready ! True if all blocks needed to update current block in 
     ! anode are ready

   call zero_task(task)
   task%task_type = TASK_UPDATE_BETWEEN
   task%src1 = src

   ! cache some values in variables
   size_snode = size(nodes(snode)%index)
   s_nb = nodes(snode)%nb
   n = blocks(src)%blkn

   !  Loop over ancestors of snode
   anode = nodes(snode)%parent
   !  Initialize cptr to correspond to the first row of the rectangular part of 
   !  the snode matrix.
   numcol = nodes(snode)%en - nodes(snode)%sa + 1
   cptr = 1 + numcol

   do while(anode.gt.0)

      ! Skip columns that come from other children
      do cptr = cptr, size_snode
         if(nodes(snode)%index(cptr).ge.nodes(anode)%sa) exit
      end do
      if(cptr.gt.size_snode) exit ! finished with snode

      map_done = .false. ! We will only build a map when we need it
      a_nb = nodes(anode)%nb

      ! Loop over affected block columns of anode
bcols: do
         if(cptr.gt.size_snode) exit
         if(nodes(snode)%index(cptr).gt.nodes(anode)%en) exit

         ! find id of the block in scol that contains the cptr-th row
         rblk = (cptr-1)/s_nb - (scol-1) + blocks(src)%dblk
  
         ! If we've gone past the src block then we are done.
         if(rblk.gt.src) exit

         ! compute local index of block column in anode and find the id of 
         ! its diagonal block
         cb = (nodes(snode)%index(cptr) - nodes(anode)%sa)/a_nb + 1
         dblk = nodes(anode)%blk_sa
         do jb = 2, cb
            dblk = blocks(dblk)%last_blk + 1
         end do

         ! Find cptr2
         jlast = min(nodes(anode)%sa + cb*a_nb - 1, nodes(anode)%en)
         do cptr2 = cptr,size_snode
            if(nodes(snode)%index(cptr2) > jlast) exit
         end do
         cptr2 = cptr2 - 1 

        ! Check blocks of snode,scol containing rows cptr to cptr2 are complete
        ! Be careful as cptr need not be at the start of a block
         newc = .false.
         do rb = (cptr-1)/s_nb, (cptr2-1)/s_nb
            blk = rb - (scol-1) + blocks(src)%dblk ! block id
            newc = newc .or. (blk.eq.src)
            if(blocks(blk)%dep.ne.-2) then
               cptr = cptr2 + 1
               cycle bcols
            endif
         end do

         ! Set info for source block csrc (hold start location
         ! and number of entries)
         k = (scol-1)*s_nb + 1
         task%csrc(1) = 1 + (cptr-k)*n
         task%csrc(2) = (cptr2 - cptr + 1)*n

         ! Build a map of anode's blocks if this is first for anode
         if(.not.map_done) then
            ! The indices for each row block in anode are mapped to a local row
            ! block index.
            size_anode = size(nodes(anode)%index)
            map_done = .true.
            jb = 1
            do i = 1, size_anode, a_nb
               do k = i, size_anode
                  k1 = nodes(anode)%index(k)
                  map(k1) = jb
               end do
               jb = jb + 1 
            end do
         endif

         ! Loop over the blocks of snode
         ready = .true.
         new = newc
         jb = map(nodes(snode)%index(cptr))
         i = cptr
         ilast = i ! Set start of current block
         do blk = rblk, blocks(rblk)%last_blk
            ready = ready .and. (blocks(blk)%dep.eq.-2)
            new = new .or. (blk.eq.src)
            rb = blk - blocks(blk)%dblk + scol ! block index 
            do i = i, min(size_snode, 1+rb*s_nb-1)
               k = map(nodes(snode)%index(i))
               if(k.ne.jb) then
                  ! Moved to a new block in anode
                  if(new .and. ready) then
                     task%dest = dblk + jb - cb
                     ! Set info for source block rsrc (hold start location
                     ! and number of entries)
                     task%rsrc(1) = 1 + (ilast-(scol-1)*s_nb-1)*n
                     task%rsrc(2) = (i - ilast)*n
                     call add_task(stack, task, control, info, st)
                     if(info.lt.0) return
                  endif
                  ready = blocks(blk)%dep.eq.-2
                  new = newc .or. (blk.eq.src)
                  jb = k
                  ilast = i ! Update start of current block
               endif
            end do
         end do
         if(new .and. ready) then
            task%dest = dblk+jb-cb
            task%rsrc(1) = 1 + (ilast-(scol-1)*s_nb-1)*n
            task%rsrc(2) = (i - ilast)*n
            call add_task(stack, task, control, info, st)
            if(info.lt.0) return
         endif

         ! Move cptr down, ready for next block column of anode
         cptr = cptr2 + 1

      end do bcols

      ! Move up the tree
      anode = nodes(anode)%parent
   end do

end subroutine add_between_updates


!*************************************************
!
! Reduce the given dependency of a block by one; 
! if it is then zero add the resultant factor or solve task
! to the local task stack (or, if full, to task pool)
!
subroutine reduce_dep(stack, block, control, info, st)

   type(taskstack), intent(inout) :: stack ! holds tasks that have to be done.
      ! it is altered if the dependency for block is reduced to 0.
   type(block_type), intent(inout) :: block ! block for which dependency
      ! count is to be reduced
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag
   integer, intent(inout) :: st ! stat parameter

   ! Local
   type(dagtask) :: task ! used to hold the task to be added

!$ call omp_set_lock(block%lock)

   ! Decrement dep count
   block%dep = block%dep - 1

   !
   ! If dep count is zero then add the appropriate task
   !
   if(block%dep.ne.0) then
!$    call omp_unset_lock(block%lock)
      return
   endif

   call zero_task(task)

   if(block%id.eq.block%dblk) then
      ! diagonal block so task is a factorization
      task%task_type = TASK_FACTORIZE_BLOCK
   else
      task%task_type = TASK_SOLVE_BLOCK
   endif
   task%dest = block%id

   block%dep = -1 ! we need this so that we don't end up adding
      ! internal update tasks twice.

   call add_task(stack, task, control, info, st)
   ! If info is non-zero we're about to return anyway

!$ call omp_unset_lock(block%lock)

end subroutine reduce_dep

!**********************************************

subroutine solve_task_dispatch(nbcol, lfact, map, npi, stack, blocks, nodes, &
      counts, control, info, st, nrhs, ldr, rhs, total_threads, rhs_local,   &
      maxmn, job)

   integer, intent(in) :: nbcol  ! Size of lfact (no. block cols)
   type(lfactor), dimension(nbcol), intent(inout) :: lfact ! Entries in 
      ! block columns of L
   integer, dimension(:), intent(inout) :: map     ! private work array
   integer, dimension(:), intent(inout) :: npi     ! cache for bwd_reduce_dep
   type(taskstack), intent(inout) :: stack         ! task pool
   type(block_type), dimension(:), intent(inout) :: blocks ! block info
   type(node_type), dimension(-1:), intent(in) :: nodes ! Node info
   type(slv_count_type), dimension(:), intent(inout) :: counts
   type(MA87_control), intent(in) :: control
   integer, intent(out) :: info ! error flag
   integer, intent(out) :: st ! stat parameter
   integer, intent(in) :: nrhs  ! number of right-hand sides (maybe = 0)
   integer, intent(in) :: ldr  ! leading extent of rhs
   real(wp), intent(inout) :: rhs(ldr*nrhs)  ! On entry holds rhs data. 
      ! Overwritten by partial solution (forward substitution performed).
   integer, intent(in) :: total_threads  ! number of threads
   real(wp), intent(inout) :: rhs_local(ldr*nrhs,0:total_threads-1)  ! Local
      ! right-hand side arrays (one for each thread).
   integer, intent(in) :: maxmn ! max block dimension (=keep%maxmn)
   integer, intent(in) :: job

   integer :: bcol ! block col that task%dest belongs to
   integer(long) :: blk ! task%dest
   integer :: col ! global index of fist col. in dcol
   integer :: dblk ! diagonal block of current column
   integer :: dcol ! local block column index
   integer :: i ! loop variable
   integer :: j ! loop variable
   integer :: m ! set to blocks(blk)%blkm (number of rows in blk)
   integer :: n ! set to blocks(blk)%blkn (number of columns in blk and dblk)
   integer :: node
   integer :: r ! current right hand side, 0 indexed
   integer :: offset ! offset into index for current block
   integer :: sa ! set to blocks(blk)%sa 
   type(dagtask) :: task
 !%%%  integer :: t_start, t_end
   integer :: this_thread
   real(wp), dimension(:), allocatable :: xlocal ! update_buffer workspace

   this_thread = 0
!$ this_thread = omp_get_thread_num()

   ! Initialize things
   info = 0; st = 0        ! By default everything went OK
   call zero_task(task)    ! Set everything to zero
   task%task_type = TASK_NONE   ! Needs to be set to prevent 
                                ! dispatched incrementing

   allocate(xlocal(maxmn*nrhs), stat=st)
   if (st.ne.0) then
      info = MA87_ERROR_ALLOCATION
      call MA87_print_flag(info, control, context='MA87_solve',st=st)
      return
   end if

   ! Mark thread as active
!$ call omp_set_lock(stack%lock)
   stack%active = stack%active + 1
!$ call omp_unset_lock(stack%lock)

   !
   ! Main loop
   !
   do
      !
      ! Retrieve next task to perform and wait till all is ready
      !
    !%%%  if(control%unit_log.gt.0) call system_clock(t_start)

      call get_task(stack, task, control, info, st)
      if(info.lt.0) return

      !write(6, "(i3,2a)") &
      !   this_thread, " got task ",print_job(task) 

    !%%%  if(control%unit_log.gt.0) then
    !%%%     call system_clock(t_end)
    !%%%     call log_task(control, this_thread, t_start, t_end, "GT", &
    !%%%        int(task%task_type,long))
    !%%%  endif

!$OMP FLUSH(stack,rhs)

      !
      ! Perform task and update dependencies
      !
      select case(task%task_type)
      case(TASK_DONE)
         exit

      case(TASK_NONE) ! Job not finished but no tasks available, spin.
         cycle

      case(TASK_SLV_FSLV) ! Perform forward solve with diagonal block

         blk = task%dest
         node = blocks(blk)%node

         ! Establish variables describing block
         n        = blocks(blk)%blkn
         m        = blocks(blk)%blkm
         sa       = blocks(blk)%sa
         bcol     = blocks(blk)%bcol
         dcol     = bcol - blocks(nodes(node)%blk_sa)%bcol + 1
         col      = nodes(node)%sa + (dcol-1)*nodes(node)%nb
         offset   = col - nodes(node)%sa + 1

         ! Sum contributions to rhs
         do r = 0, nrhs-1
            do j = 0, total_threads-1
               do i = col + r*ldr, col+n-1 + r*ldr
                  rhs(i) = rhs(i) + rhs_local(i, j)
                  rhs_local(i,j) = zero ! Reset in case of bwd solve
               end do
            end do
         end do

         ! Perform triangular solve
         call slv_solve(n, n, col, lfact(bcol)%lcol(sa:sa+n*n-1), &
            'Transpose    ', 'Non-unit', nrhs, rhs, ldr, control, &
             blocks(blk)%id)
         offset = offset + n

!$OMP FLUSH

         ! Add a bwd solve if appropriate
         if(blocks(blk)%dep.eq.0 .and. job.eq.SOLVE_JOB_ALL) then
            task%task_type = TASK_SLV_BSLV
            task%dest = blk
            call add_task(stack, task, control, info, st)
            if(info.lt.0) return
         endif

         ! Loop over sub-diagonal blocks adding tasks
         task%task_type = TASK_SLV_FUPD
         do blk = blocks(blk)%dblk+1, blocks(blk)%last_blk
            task%dest = blk
            call add_task(stack, task, control, info, st)
            if(info.lt.0) return
         end do

         ! Deal with any left over trapezoidal part of diagonal block
         m = m - n
         if(m.gt.0) then
            sa = sa + n*n
            call slv_fwd_update(m, n, col, offset, nodes(node)%index,    &
               lfact(bcol)%lcol(sa:sa+n*m-1), n, nrhs,                   &
               rhs_local(:,this_thread), ldr, rhs, ldr, xlocal, control, &
               blocks(blk)%id)
!$OMP FLUSH
            call fwd_reduce_dep(stack, counts, map, offset, m, &
               nodes(node)%index, control, info, st)
         endif

      case(TASK_SLV_FUPD) ! Forward update with off-diagonal
         blk = task%dest

         node = blocks(blk)%node

         ! Establish variables describing block
         n        = blocks(blk)%blkn
         m        = blocks(blk)%blkm
         sa       = blocks(blk)%sa
         bcol     = blocks(blk)%bcol
         dcol     = bcol - blocks(nodes(node)%blk_sa)%bcol + 1
         col      = nodes(node)%sa + (dcol-1)*nodes(node)%nb

         offset   = col - nodes(node)%sa + 1 ! diagonal blk
         offset   = offset + (blk-blocks(blk)%dblk) * nodes(node)%nb ! this blk

         call slv_fwd_update(m, n, col, offset, nodes(node)%index,            &
            lfact(bcol)%lcol(sa:sa+n*m-1), n, nrhs, rhs_local(:,this_thread), &
            ldr, rhs, ldr, xlocal, control, blocks(blk)%id)
!$OMP FLUSH
         call fwd_reduce_dep(stack, counts, map, offset, m, &
            nodes(node)%index, control, info, st)

      case(TASK_SLV_BSLV) ! Backward solve with diagonal block
         blk = task%dest
         node = blocks(blk)%node

         ! Establish variables describing block column
         n      = blocks(blk)%blkn
         m      = blocks(blk)%blkm
         sa     = blocks(blk)%sa
         bcol   = blocks(blk)%bcol
         col    = calc_col(nodes(node), blocks(blk)) ! current bcol
         col    = nodes(node)%sa + (col-1)*nodes(node)%nb
         offset = col - nodes(node)%sa + 1

         ! Perform and retangular update from diagonal block
         if(m.gt.n) then
            call slv_bwd_update(m-n, n, col, offset+n, nodes(node)%index, &
               lfact(bcol)%lcol(sa+n*n:sa+n*m-1), n, nrhs, rhs,           &
               rhs_local(:,this_thread), ldr, xlocal, control, blocks(blk)%id)
         endif

         ! Sum contributions to rhs
         do r = 0, nrhs-1
            do j = 0, total_threads-1
               do i = col + r*ldr, col+n-1 + r*ldr
                  rhs(i) = rhs(i) + rhs_local(i, j)
               end do
            end do
         end do

         ! Perform triangular solve
         call slv_solve(n, n, col, lfact(bcol)%lcol(sa:sa+n*n-1), &
            'Non-Transpose', 'Non-unit', nrhs, rhs, ldr, control, &
             blocks(blk)%id)

!$OMP FLUSH

         ! Reduce dependencies of blocks associated with this bcol's vars
         call bwd_reduce_dep(node, bcol, nodes, blocks, map, npi, stack, &
            control, info, st)
         if(info.lt.0) return

      case(TASK_SLV_BUPD) ! Backward update with off-diagonal
         blk      = task%dest

         ! Establish variables describing block
         n        = blocks(blk)%blkn
         m        = blocks(blk)%blkm
         sa       = blocks(blk)%sa
         bcol     = blocks(blk)%bcol
         node     = blocks(blk)%node
         dcol     = bcol - blocks(nodes(node)%blk_sa)%bcol + 1
         col      = nodes(node)%sa + (dcol-1)*nodes(node)%nb

         offset   = col - nodes(node)%sa + 1 ! diagonal blk
         offset   = offset + (blk-blocks(blk)%dblk) * nodes(node)%nb ! this blk

         call slv_bwd_update(m, n, col, offset, nodes(node)%index,  &
            lfact(bcol)%lcol(sa:sa+n*m-1), n, nrhs, rhs,            &
            rhs_local(:,this_thread), ldr, xlocal, control, blocks(blk)%id)

!$OMP FLUSH

         ! Reduce diagonal block's dependencies, add bslv task if appropriate
         dblk = blocks(blk)%dblk
!$       call omp_set_lock(blocks(dblk)%lock)
         blocks(dblk)%dep = blocks(dblk)%dep - 1
         if(blocks(dblk)%dep.eq.0) then
            task%task_type = TASK_SLV_BSLV
            task%dest = dblk
            call add_task(stack, task, control, info, st)
            if(info.lt.0) return
         endif
!$       call omp_unset_lock(blocks(dblk)%lock)

      case default
         info = MA87_ERROR_UNKNOWN
         call MA87_print_flag(info, control, context='MA87_solve')
       !  if(control%diagnostics_level.ge.0 .and. control%unit_error.ge.0) &
       !     write(control%unit_error, "(/a,i3,a,i8)") &
       !     " MA87_factor: Internal Error ", info, &
       !     " Unknown task type encountered type = ", task%task_type
         exit
      end select

   end do

end subroutine solve_task_dispatch


!*************************************************

! Returns the destination block of an internal update task.
! Called by add_updates.

integer(long) function get_dest_block(src1, src2)

   type(block_type), intent(in) :: src1
   type(block_type), intent(in) :: src2

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Numerical block operation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! TASK_SOLVE_BLOCK
! Solve using factorization of diag. block (uses dtrsm)
! A_ij <- A_ij A_ii^-1
! dest <- dest diag^-1
!
subroutine solv_col_block(m, n, id, dest, id1, diag, control)

   integer, intent(in) :: m ! number of rows in dest
   integer, intent(in) :: n ! number of columns in dest
   integer(long), intent(in) :: id ! block identifier for dest
   integer(long), intent(in) :: id1 ! block identifier for diag
   real(wp), dimension(*), intent(inout) :: dest ! holds destination block
   real(wp), dimension(*), intent(inout) :: diag ! block
     ! on diagonal of factor L in same block col as dest.
   type(MA87_control), intent(in) :: control
   

 !%%%  integer :: t_start, t_end, this_thread

  !%%% if(control%unit_log.gt.0) call system_clock(t_start)

   call dtrsm('Left', 'Upper', 'Transpose', 'Non-Unit', n, m, &
      one, diag, n, dest, n)

!%%%   if(control%unit_log.gt.0) then
!%%%      call system_clock(t_end)
!%%%      this_thread = 0
!%%% !$    this_thread = omp_get_thread_num()
!%%%      call log_task(control, this_thread, t_start, t_end, "TS", &
!%%%         id, id1)
!%%%   endif
end subroutine solv_col_block

!*************************************************

!
! TASK_FACTORIZE_BLOCK (uses Lapack routine dpotrf and dtrsm)
! A_ii <- L_ii
!
subroutine factor_diag_block(n, m, id, dest, control, dpotrf_info, detlog)

   integer, intent(in) :: m ! number of rows in dest
   integer, intent(in) :: n ! number of columns in dest
   integer(long), intent(in) :: id ! block identifier for dest
   real(wp), dimension(*), intent(inout) :: dest ! holds block
      ! on diagonal of factor L. It may not be square
   type(MA87_control), intent(in) :: control
   integer, intent(out) :: dpotrf_info ! error flag for dpotrf
   real(wp), intent(inout) :: detlog ! accumlated sum of absolute values
      ! of pivots on this thread

   integer :: i, j ! Loop indices
 !%%%  integer :: t_start, t_end, this_thread ! used for timings

 !%%%  if(control%unit_log.gt.0) call system_clock(t_start)


   call dpotrf('Upper', n, dest, n, dpotrf_info)
   ! check for errors
   if(dpotrf_info.ne.0) return

   ! Work out detlog
   j = 1
   do i = 1, n
      detlog = detlog + log(dest(j))
      j = j + n + 1
   end do

   ! Do dtrsm with any remainder below diagonal block
   if(m.gt.n) then
      call dtrsm('Left', 'Upper', 'Transpose', 'Non-Unit', n, &
                  m-n, one, dest, n, dest(1+n*n), n)
   endif

!%%%   if(control%unit_log.gt.0) then
!%%%      call system_clock(t_end)
!%%%      this_thread = 0
!%%%!$    this_thread = omp_get_thread_num()
!%%%      call log_task(control, this_thread, t_start, t_end, "DF", id)
!%%%   endif

end subroutine factor_diag_block

!*************************************************

!
! TASK_UPDATE_INTERNAL
! A_ik <- A_ik - A_ij A_kj^T
! dest <- dest - src2 src1^T
! Remember that the blocks are stored by rows.
! dest, src1 and src2 all belong to the same node.
!
subroutine update_block_block(m, n, dest, blk, n1, src1, src2, control)

   integer, intent(in) :: m ! number of rows in dest
   integer, intent(in) :: n ! number of columns in dest
   real(wp), dimension(*), intent(inout) :: dest ! holds block in L
     ! that is to be updated. 
   type(block_type), intent(inout) :: blk ! destination block  
   integer, intent(in) :: n1 ! number of columns in src1 and src2
   real(wp), dimension(*), intent(in) :: src1
   real(wp), dimension(*), intent(in) :: src2
   type(MA87_control), intent(in) :: control

   logical :: diag ! set to true if dest is the diagonal block
!%%%   integer :: t_start, t_end, this_thread

   diag = (blk%dblk.eq.blk%id)

!%%%   if(control%unit_log.gt.0) call system_clock(t_start)

!$ call omp_set_lock(blk%alock)

   if(diag) then
      call dsyrk('U', 'T', n, n1, -one, src1, n1, one, dest, n)

      if(m.gt.n) then
         ! Do a dgemm on any remainder
         call dgemm('T', 'N', n, m-n, n1, -one,                   &
            src1, n1, src2(1+n*n1), n1, one, dest(1+n*n), n)
      endif
   else
      ! dest is an off-diagonal block
      call dgemm('T', 'N', n, m, n1, -one, src1, n1, src2, n1, one, dest, n)
   endif

!$ call omp_unset_lock(blk%alock)

!%%%   if(control%unit_log.gt.0) then
!%%%      call system_clock(t_end)
!%%%      this_thread = 0
!%%% !$    this_thread = omp_get_thread_num()
!%%%      call log_task(control, this_thread, t_start, t_end, "UW")
!%%%   endif

end subroutine update_block_block

!*************************************************

!   Given a destination block dest, update_between performs the update
!                     L_dest <-- L_dest - L_rsrc (L_csrc)^T
!   where L_dest is a submatrix of the block dest of an ancestor
!   of the node snode and L_rsrc and L_csrc are submatrices of contiguous
!   rows that belong to the same block column of snode as the block src
!   (this block col. has local index scol).
!   The first row of L_rsrc is the first row
!   of the block column that corresponds to a row in the block dest and the 
!   last row of L_rsrc is the last row of the block column that corresponds to 
!   a row in the block dest. Similarly, the first/last row of L_csrc is the
!   first/last row of the block column that corresponds to a column in the 
!   block dest. The set of rows and columns of dest thus
!   determine which two sets of contiguous rows in scol are involved.
!   Unless the number of entries updated is very small, use BLAS 3 kernel 
!   gemm or syrk by placing its result in a buffer from which we add the 
!   update into the appropriate entries of the
!   destination block dest.

subroutine update_between(m, n, blk, dnode, n1, src, snode, dest, csrc,    &
      rsrc, blocks, row_list, col_list, buffer, control, info, st)

   integer, intent(in) :: m  ! number of rows in destination block
   integer, intent(in) :: n  ! number of columns in destination block
   integer(long), intent(in) :: blk ! identifier of destination block
   type(node_type), intent(in) :: dnode ! Node to which blk belongs
   integer :: n1 ! number of columns in source block column
   integer(long), intent(in) :: src  ! identifier of block in source block col
   type(node_type), intent(in) :: snode ! Node to which src belongs
   real(wp), dimension(*), intent(inout) :: dest ! holds block in L
     ! that is to be updated.
   real(wp), dimension(*), intent(in) :: csrc ! holds csrc block
   real(wp), dimension(*), intent(in) :: rsrc ! holds rsrc block
   type(block_type), dimension(:), intent(inout) :: blocks
   real(wp), dimension(:), allocatable :: buffer
   integer, dimension(:), allocatable :: row_list ! reallocated to min size m
   integer, dimension(:), allocatable :: col_list ! reallocated to min size n
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info   
   integer, intent(inout) :: st


   ! Local scalars
   integer :: cptr ! used in determining the rows that belong to the
     ! source block csrc
   integer :: col_list_sz ! initialise to 0. then increment while
     ! rows involed in csrc are recorded, until
     ! holds the number of columns in blk (= number
     ! of rows in csrc)
   integer :: dcen ! index of end column in dcol
   integer :: dcol ! index of block column that blk belongs to in dnode
   integer :: dcsa ! index of first column in dcol
   logical :: diag ! set to true if blk is the diagonal block
   integer :: dptr
   integer :: dptr_sa
   integer :: drsa, dren ! point to first and last rows of destination
     ! block blk
   integer :: i

   integer :: ndiag ! set to int(s1en-s1sa+1) if blk is a
     ! block on diagonal and 0 ow. so is number of triangular rows of update
   integer :: row_list_sz ! initialise to 0. then increment while
     ! rows involed in rsrc are recorded, until
     ! holds the number of rows in blk (= number of rows in rsrc)
   integer :: rptr ! used in determining the rows that belong to the
     ! source block rsrc
   integer :: scol ! index of block column that src belongs to in snode
   integer :: s1sa, s1en ! point to the first and last rows of
     ! the block csrc within scol
   integer :: s2sa, s2en ! point to the first and last rows of
     ! the block rsrc within scol
   integer :: size_dnode ! size(dnode%index)
   integer :: size_snode ! size(snode%index)

 !%%%  integer :: t_start, t_end, this_thread

 !%%%  if(control%unit_log.gt.0) call system_clock(t_start)

   ! set diag to true if blk is block on diagonal
   diag = (blocks(blk)%dblk.eq.blocks(blk)%id)

   ! Make a list of incident csrc rows (ie. columns of blk)
   !
   ! Initialize lists
   if(size(col_list).lt.n) then
      deallocate(col_list, stat=st)
      allocate(col_list(n), stat=st)
      if (st.ne.0) go to 10
   endif
   if(size(row_list).lt.m) then
      deallocate(row_list, stat=st)
      allocate(row_list(m), stat=st)
      if (st.ne.0) go to 10
   endif

 10   if (st.ne.0) then
         info = MA87_ERROR_ALLOCATION
         call MA87_print_flag(info, control, context='MA87_factor',st=st)
         return
      end if

   col_list_sz = 0
   row_list_sz = 0

   size_dnode = size(dnode%index)
   size_snode = size(snode%index)

   ! Find block column dcol of dnode that blk belongs to. The block
   ! cols are numbered locally within dnode as 1,2,3,...

   dcol = blocks(blk)%bcol - blocks(dnode%blk_sa)%bcol + 1

   ! Set dcsa and dcen to hold indices
   ! of start and end columns in dcol (global column indices)
   dcsa = dnode%sa + (dcol-1)*dnode%nb                
   dcen = min(dnode%sa + dcol*dnode%nb-1, dnode%en)

   ! Find block column scol of snode that src belongs to. 
   scol = blocks(src)%bcol - blocks(snode%blk_sa)%bcol + 1

   ! Set cptr to point to the first row in csrc
   cptr = 1 + min(snode%en-snode%sa+1, (scol-1)*snode%nb)
 
   ! loop while row index within scol is less the index
   ! of the first column in blk

   do while(snode%index(cptr).lt.dcsa)
      cptr = cptr + 1
      if(cptr.gt.size_snode) return ! No incident columns
   end do

   ! Set s1sa to point to first row in csrc
   s1sa = cptr 

   ! Now record the rows in csrc. Local row numbers
   ! are held in col_list(1:slen-slsa+1)
   do while(snode%index(cptr).le.dcen)
      col_list_sz = col_list_sz + 1
      col_list(col_list_sz) = snode%index(cptr) - dcsa + 1
      cptr = cptr + 1
      if(cptr.gt.size_snode) exit ! No more rows
   end do

   ! Set slen to point to last row in csrc
   s1en = cptr - 1 
   
   ! Loop over rsrc rows, building row list. Identify required data, form
   ! outer product of it into buffer.

   ! Find first and last rows of destination block
   i = dcol + blocks(blk)%id - blocks(blk)%dblk ! block in snode
   drsa = dnode%index(1 + (i-1)*dnode%nb)
   dren = dnode%index(min(1 + i*dnode%nb - 1, size_dnode))

   ! Find first row in rsrc
   rptr = s1sa
   do while(snode%index(rptr).lt.drsa)
      rptr = rptr + 1
      if(rptr.gt.size_snode) return ! No incident row! Shouldn't happen.
   end do
   s2sa = rptr ! Points to first row in rsrc

   ! Find the first row of destination block
   i = blk - blocks(blk)%dblk + 1 ! row block of blk column
   dptr_sa = 1 + (dcol-1 + i-1)*dnode%nb

   ! Now record the rows in rsrc. Local row numbers
   ! are held in row_list(1:s2en-s2sa+1)
   dptr = dptr_sa ! Pointer for destination block

   do rptr = s2sa, size_snode
      if(snode%index(rptr).gt.dren) exit
      do while(dnode%index(dptr).lt.snode%index(rptr))
         dptr = dptr + 1
      end do
      row_list_sz = row_list_sz + 1
      row_list(row_list_sz) = dptr - dptr_sa + 1
   end do
   s2en = rptr - 1 ! Points to last row in rsrc

   if(n1.ge.control%min_width_blas) then
      ! High flop/buffer sz ratio => perform operations into buffer with BLAS
      if(size(buffer).lt.row_list_sz*col_list_sz) then
         deallocate(buffer, stat=st)
         allocate(buffer(row_list_sz*col_list_sz), stat=st)
         if (st.ne.0) then
            info = MA87_ERROR_ALLOCATION
            call MA87_print_flag(info, control, context='MA87_factor',st=st)
            return
         end if
      endif

      if(diag) then
         ! blk is a block on diagonal
         ndiag = int(s1en-s1sa+1)
         call dsyrk('U', 'T', ndiag, n1, -one, csrc,                  &
            n1, zero, buffer, col_list_sz)

         if(s2en-s2sa+1-ndiag.gt.0) then
            call dgemm('T', 'N', ndiag, int(s2en-s2sa+1-ndiag),       &
               n1, -one, csrc, n1, rsrc(1+n1*ndiag), n1, zero,        &
               buffer(1+col_list_sz*ndiag), col_list_sz)
         endif
      else
         ! Off diagonal block
         ndiag = 0
         call dgemm('T', 'N', int(s1en-s1sa+1), int(s2en-s2sa+1), n1, &
            -one, csrc, n1, rsrc, n1, zero, buffer, col_list_sz)
      endif

      !
      ! Apply update
      !

      ! Acquire lock on destination block
!$    call omp_set_lock(blocks(blk)%alock)

      call expand_buffer(dest, n, row_list, row_list_sz, &
         col_list, col_list_sz, ndiag, buffer)

      ! Release lock
!$    call omp_unset_lock(blocks(blk)%alock)

   else
      ! Low flop/buffer ratio => perform update operation directly
      ! set ndiag if blk is a diagonal block
      ndiag = 0
      if(diag) ndiag = int(s1en-s1sa+1)

!$    call omp_set_lock(blocks(blk)%alock)

      call update_direct(n, dest, n1, csrc, rsrc, row_list, row_list_sz, &
         col_list, col_list_sz, ndiag)

!$    call omp_unset_lock(blocks(blk)%alock)

   endif

!%%%   if(control%unit_log.gt.0) then
!%%%      call system_clock(t_end)
!%%%      this_thread = 0
!%%% !$    this_thread = omp_get_thread_num()
!%%%      call log_task(control, this_thread, t_start, t_end, "UB", blk,  &
!%%%         int(snode%sa,long), int(scol,long))
!%%%   endif

end subroutine update_between

!*************************************************

! Performs an update_between task directly rather than via a buffer.

subroutine update_direct(n, dest, n1, csrc, rsrc,  &
      row_list, rls, col_list, cls, ndiag)

   integer, intent(in) :: n ! number of columns in dest
   real(wp), dimension(*), intent(inout) :: dest ! holds destination block
   integer, intent(in) :: n1 ! number of columns in source blocks
   real(wp), dimension(*), intent(in) :: csrc, rsrc ! source blocks
   integer, intent(in) :: rls  ! size of row_list
   integer, intent(in) :: row_list(rls) ! local row indices for
     ! rows in rsrc (= local row indices for dest)
   integer, intent(in) :: cls  ! size of col_list
   integer, intent(in) :: col_list(cls) ! local row indices for
     ! rows in csrc (= local column indices for dest)
   integer, intent(in) :: ndiag ! Number of triangular rows of update

   integer :: cptr
   integer :: i
   integer :: j
   integer :: k1, k2, k3, k4, l
   integer :: rptr
   real(wp) :: work1, work2, work3, work4

   ! Note: ndiag is either equal to 0 or, if dest is a diagonal block, it
   ! is equal to the number of cols cls in dest (and this is
   ! equal to the number of rows rls in dest).
   ! first treat the case when dest is on the diagonal.
   ! loop over the rows of dest
   do j = 1, ndiag
      cptr = (row_list(j)-1)*n
      rptr = 1 + (j-1)*n1
      do i = 1, j
         k1 = 1 + (i-1)*n1
         work1 = zero
         do l = rptr, rptr + n1 - 1
            work1 = work1 + csrc(k1)*rsrc(l)
            k1 = k1 + 1
         end do
         k1 = cptr + col_list(i)
         dest(k1) = dest(k1) - work1
      end do
   end do

   ! Now consider the case when dest is not a diagonal block
   do j = ndiag+1, rls
      cptr = (row_list(j)-1)*n
      rptr = 1 + (j-1)*n1
      i = 4*int(cls/4)
      do i = 1, i, 4
         k1 = 1 + (i+0-1) * n1
         k2 = 1 + (i+1-1) * n1
         k3 = 1 + (i+2-1) * n1
         k4 = 1 + (i+3-1) * n1
         work1 = zero; work2 = zero; work3 = zero; work4 = zero
         do l = rptr, rptr + n1 - 1
            work1 = work1 + csrc(k1)*rsrc(l); k1 = k1 + 1
            work2 = work2 + csrc(k2)*rsrc(l); k2 = k2 + 1
            work3 = work3 + csrc(k3)*rsrc(l); k3 = k3 + 1
            work4 = work4 + csrc(k4)*rsrc(l); k4 = k4 + 1
         end do
         k1 = cptr + col_list(i+0); dest(k1) = dest(k1) - work1
         k2 = cptr + col_list(i+1); dest(k2) = dest(k2) - work2
         k3 = cptr + col_list(i+2); dest(k3) = dest(k3) - work3
         k4 = cptr + col_list(i+3); dest(k4) = dest(k4) - work4
      end do
      i = 4*int(cls/4) + 1
      do i = i, cls
         k1 = 1 + (i-1)*n1
         work1 = zero
         do l = rptr, rptr + n1 - 1
            work1 = work1 + csrc(k1)*rsrc(l)
            k1 = k1 + 1
         end do
         k1 = cptr + col_list(i)
         dest(k1) = dest(k1) - work1
      end do
   end do

end subroutine update_direct

!*************************************************

! Optimised sparse expansion by lists
! Entry buffer(i,j) gets added to a(row_list(i),col_list(j))
! Note: by uncommenting the i2 and mm_prefetch lines we can speed this up
!       for the Intel compiler by explicitly prefetching the next row.

subroutine expand_buffer(a, blkn, row_list, rls, col_list, cls, ndiag, buffer)

   real(wp), dimension(*), intent(inout) :: a ! holds L
   integer, intent(in) :: blkn ! number of cols in destination block
   integer, intent(in) :: rls ! size of row_list
   integer, intent(in) :: row_list(rls)
   integer, intent(in) :: cls ! size of col_list
   integer, intent(in) :: col_list(cls)
   integer, intent(in) :: ndiag ! Number of triangular rows of update
   real(wp), intent(in) :: buffer(rls*cls)

   integer :: i, j, k, rptr, cptr
   !integer :: i2

   rptr = 1
   do j = 1, ndiag
      cptr = 1 + (row_list(j)-1)*blkn - 1
      do i = 1, j
         k = cptr + col_list(i)
         a(k) = a(k) + buffer(rptr)
         rptr = rptr + 1
      end do
      rptr = rptr + (cls-j)
   end do
   do j = ndiag+1, rls-1
      cptr = (row_list(j)-1) * blkn
      !i2 = (row_list(j+1)-row_list(j)) * blkn
      do i = 1, cls
         k = cptr + col_list(i)
         a(k) = a(k) + buffer(rptr)
         !call mm_prefetch(a(k+i2),1)
         rptr = rptr + 1
      end do
   end do

   if(ndiag.lt.rls) then
      cptr = (row_list(rls)-1) * blkn
      do i = 1, cls
         k = cptr + col_list(i)
         a(k) = a(k) + buffer(rptr)
         rptr = rptr + 1
      end do
   endif

end subroutine expand_buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine updates the dependency count for blocks in the ancestor
! anode of snode.

subroutine calc_dep(cptr, snode, anode, nodes, blocks, swidth, map)

   integer, intent(inout) :: cptr 
     ! pointer into the row indices for snode. On entry points
     ! to the first row corresponding to the column set for anode.
     ! On exit, updated ready for the next ancestor.
     ! It could feasibly point beyond the column set for anode if
     ! there are no entries for anode's columns in the row set of snode.
   integer, intent(in) :: snode ! node
   integer, intent(in) :: anode ! ancestor of snode
   type(node_type), dimension(-1:), intent(in) :: nodes ! node info
   type(block_type), dimension(:), intent(inout) :: blocks ! block info. On 
     ! exit, dependency count updated
   integer, intent(in) :: swidth ! number of block cols in snode
   integer, dimension(:), intent(in) :: map ! For each row k in j-th 
     ! block column of anode, map(k) is set to j

   integer :: cb ! index of block column in anode
   integer(long) :: dblk ! id of diagonal block in anode
   integer :: jlast ! last column of block in anode
   integer :: nb ! set to nodes(anode)%nb (block size for anode)
   integer :: i, jb, k, k1
   integer :: size_snode ! number of rows in snode

   nb = nodes(anode)%nb
   size_snode = size(nodes(snode)%index)
   do 
      if(cptr.gt.size_snode) exit
      if (nodes(snode)%index(cptr).gt.nodes(anode)%en) exit

      ! Find block column of anode
      cb = 1 + (nodes(snode)%index(cptr) - nodes(anode)%sa) / nb

      ! Find diagonal block in block column cb

      ! loop over block columns. blocks(dblk)%last_blk is the last
      ! block in the block column to which dblk belongs and
      ! so blocks(dblk)%last_blk + 1 is first block in the next
      ! block column, which is the diagonal block in that block column

      dblk = nodes(anode)%blk_sa ! first block in anode
      do i = 1, cb-1
         dblk = blocks(dblk)%last_blk + 1
      end do

      ! Increment dep count for each block involved.
      ! loop over rows in snode
      jb = -1 ! Last block
      do i = cptr, size_snode
         k1 = nodes(snode)%index(i)
         k = map(k1)
         ! k is local block number. When we reach a new block,
         ! we increment dep by the number of block columns in snode
         ! (since each block col. in snode will have to update this block)
         if(k.ne.jb) blocks(dblk+k-cb)%dep_initial = &
            blocks(dblk+k-cb)%dep_initial + swidth
         jb = k
      end do

      ! Move cptr to first row in another block of anode
      jlast = min(nodes(anode)%sa + cb*nb - 1, nodes(anode)%en)
      do cptr = cptr, size_snode
         if(nodes(snode)%index(cptr) > jlast) exit
      end do

   end do

end subroutine calc_dep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Initialize the task pool, which is held using a variable of derived
! type taskstack
!
subroutine init_stack(stack, pool_size, control, info, st)

   type(taskstack), intent(out) :: stack ! see description of derived type
   integer, intent(in) :: pool_size
   type(MA87_control), intent(in) :: control  ! see description of derived type
   integer, intent(out) :: info ! error flag. Only possible error is
     ! an allocation error
   integer, intent(out) :: st ! stat parameter

   ! local variables
   integer :: i
   integer :: ncache
   integer :: total_threads


   ! Initialise
   info = 0; st = 0
   total_threads = 1
!$ total_threads = omp_get_max_threads()

   stack%pool_size = pool_size
   stack%total = 0
   stack%active = 0
   stack%freehead = 1
   stack%abort = .false.

   ncache = calc_cache(total_threads-1, control)

   deallocate(stack%ctasks,stat=st)
   deallocate(stack%cheads,stat=st)
!$ deallocate(stack%clocks,stat=st)
   deallocate(stack%tasks,stat=st)
!**deallocate(stack%waiting,stat=st)
   deallocate(stack%next,stat=st)

   allocate(stack%ctasks(control%cache_tq_sz, ncache),  &
            stack%cheads(ncache),   &
!$          stack%clocks(ncache),   &
            stack%tasks(pool_size), &
            stack%next(pool_size),  &
 !**        stack%waiting(0:total_threads-1), &
            stat=st)

!$ call omp_init_lock(stack%lock)

   if (st.ne.0) then 
     info = MA87_ERROR_ALLOCATION
     return
   end if

   ! Initialise stack
   do i = 1, stack%pool_size-1
      stack%next(i) = i + 1
   end do
   stack%next(stack%pool_size) = -1
 !**  stack%waiting(:) = 0.0
   stack%cheads(:) = 0
   stack%prihead(1:4) = -1 ! empty

!$ do i = 1, ncache
!$    call omp_init_lock(stack%clocks(i))
!$ end do

end subroutine init_stack

!*************************************************
!
! Free any resources involved in task pool
!
subroutine cleanup_stack(stack)

   type(taskstack), intent(inout) :: stack ! see description of derived type

   ! local variables
   integer :: i ! temporary variable
   integer :: st ! stat parameter

!$ if(allocated(stack%tasks)) call omp_destroy_lock(stack%lock)

   deallocate(stack%tasks,stat=st)
   deallocate(stack%ctasks,stat=st)
   deallocate(stack%cheads,stat=st)
!**deallocate(stack%waiting,stat=st)
   deallocate(stack%next,stat=st)

!$  if(allocated(stack%clocks)) then
!$    do i = 1, size(stack%clocks)
!$       call omp_destroy_lock(stack%clocks(i))
!$    end do
!$    deallocate(stack%clocks,stat=st)
!$  endif

end subroutine cleanup_stack

!*************************************************
!
! This subroutine ensures all components of task are defined
!
subroutine zero_task(task)

   type(dagtask), intent(out) :: task

   task%task_type = 0
   task%dest = 0
   task%src1 = 0
   task%src2 = 0
   task%csrc(:) = 0
   task%rsrc(:) = 0
end subroutine zero_task

!*************************************************

!
! Add a task to the local stack or the task pool.
! In fact we add it to the cache's local stack, if this is full then it
! gets thrown to the task pool. 
!
subroutine add_task(stack, task, control, info, st)

   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(dagtask), intent(in) :: task ! task to be added
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(inout) :: st

   integer :: cache, i, j, this_thread

   this_thread = 0
!$ this_thread = omp_get_thread_num()

   cache = calc_cache(this_thread, control)

!$ call omp_set_lock(stack%clocks(cache))

   if(stack%cheads(cache).eq.control%cache_tq_sz) then
      ! Local stack is full.
      ! Clear lower half of stack to task pool.
!$    call omp_set_lock(stack%lock)
      do i = 1, control%cache_tq_sz / 2
         call add_task_g(stack, stack%ctasks(i,cache), &
            control, info, st, locked=.true.)
         if(info.lt.0) then
!$          call omp_unset_lock(stack%lock)
!$          call omp_unset_lock(stack%clocks(cache))
            return
         endif
      end do
!$    call omp_unset_lock(stack%lock)
      j = 1
      do i = control%cache_tq_sz / 2 + 1, control%cache_tq_sz
         stack%ctasks(j,cache) = stack%ctasks(i,cache)
         j = j + 1
      end do
      stack%cheads(cache) = stack%cheads(cache) - control%cache_tq_sz / 2
   endif

   ! Add to top of local stack
   stack%cheads(cache) = stack%cheads(cache) + 1
   stack%ctasks(stack%cheads(cache), cache) = task

!$ call omp_unset_lock(stack%clocks(cache))

end subroutine add_task

!*************************************************

!
! Add a task to the task pool at a given priority.
! The tasks have different priorities and those with the same
! priority are held within the pool using a linked list
!
subroutine add_task_g(stack, task, control, info, st, locked)

   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(dagtask), intent(in) :: task ! task to be added to task pool.
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag (only possible error is
     ! an allocation error)
   integer, intent(inout) :: st ! stat parameter
   logical, optional, intent(in) :: locked ! indicates whether stack%lock
     ! is set. If not present, stack%lock will be set and then unset before
     ! return

   integer :: i
   integer :: insert ! holds free location in task pool
   integer :: priority ! priority of task
   type(dagtask), dimension(:), allocatable :: temp_tasks ! only allocated if
     ! task pool is full and has to be reallocated, in which case
     ! temp_task is used to hold a temporary copy.
   integer, dimension(:), allocatable :: temp_next ! used only if
     ! task pool has to be increased in size.

   ! set priority according to type of task
   select case(task%task_type)
   case(TASK_FACTORIZE_BLOCK)
      priority = 1
   case(TASK_SOLVE_BLOCK)
      priority = 2
   case(TASK_UPDATE_INTERNAL)
      priority = 3
   case(TASK_UPDATE_BETWEEN)
      priority = 4
   case(TASK_SLV_FSLV, TASK_SLV_BSLV)
      priority = 1
   case(TASK_SLV_FUPD, TASK_SLV_BUPD)
      priority = 2
   case default
      priority = 5
   end select

!$ if(.not.present(locked)) call omp_set_lock(stack%lock)

   !
   ! Find a free spot
   !
   insert = stack%freehead
   if(insert.eq.-1) then
      ! We have overflowed the task pool, we need to reallocate it

      ! Copy tasks
      allocate(temp_tasks(stack%pool_size),stat=st)
      if(st.ne.0) go to 10
      temp_tasks(:) = stack%tasks(:)

      deallocate(stack%tasks,stat=st)
      allocate(stack%tasks(stack%pool_size*2),stat=st)
      if(st.ne.0) go to 10

      stack%tasks(1:stack%pool_size) = temp_tasks(:)
      deallocate(temp_tasks,stat=st)

      ! Copy next
      allocate(temp_next(stack%pool_size),stat=st)
      if(st.ne.0) go to 10
      temp_next(:) = stack%next(:)

      deallocate(stack%next,stat=st)
      allocate(stack%next(stack%pool_size*2),stat=st)
      if(st.ne.0) go to 10

      stack%next(1:stack%pool_size) = temp_next(:)
      deallocate(temp_next,stat=st)

      ! Extend linked list
      stack%freehead = stack%pool_size + 1
      do i = stack%pool_size+1, stack%pool_size*2-1
         stack%next(i) = i + 1
      end do
      stack%next(stack%pool_size*2) = -1

      ! Increase stored pool_size
      stack%pool_size = stack%pool_size*2
      insert = stack%freehead

      info = MA87_WARNING_POOL_SMALL

      if(control%diagnostics_level.ge.0 .and. control%unit_warning.ge.0) &
         write(control%unit_warning, "(/a,i3,2a,i8)")&
         " MA87_factor: Warning:", info, " Task pool overflow, ", &
         "size increased to = ", stack%pool_size
   end if

   stack%freehead = stack%next(insert)

   !
   ! Place task in pool in position insert and add into the linked list of 
   ! tasks with same priority
   !
   stack%tasks(insert) = task
   stack%next(insert) = stack%prihead(priority)
   stack%prihead(priority) = insert
   stack%total = stack%total + 1

   !
   ! Update highest priority (the task with the highest priority
   ! is the one with the lowest priority value)
   !
   stack%lowest_priority_value = min(stack%lowest_priority_value, priority)

!$ if(.not.present(locked)) call omp_unset_lock(stack%lock)
   return ! Successful exit

   ! Error handling in case of allocation failure
   10 continue
   info = MA87_ERROR_ALLOCATION
   call MA87_print_flag(info, control, context='MA87_factor',st=st)
   stack%abort = .true.
!$ if(.not.present(locked)) call omp_unset_lock(stack%lock)
   return

end subroutine add_task_g

!*************************************************
!
! Figure out which cache this thread belongs to - allows easy changes between
! machines by requiring only one function be altered.
!
integer function calc_cache(thread, control)

   integer, intent(in) :: thread
   type(MA87_control), intent(in) :: control

   integer :: total_threads

   select case (control%cache_layout)
   case(CACHE_COMPACT)
      calc_cache = thread / control%cache_cores + 1
   case(CACHE_SCATTER)
      total_threads = 1
!$    total_threads = omp_get_max_threads()
      total_threads = max(1, total_threads/control%cache_cores)
         ! (Avoid div by zero)
      calc_cache = mod(thread, total_threads) + 1
   case default ! CACHE_IDENITY
      calc_cache = thread + 1
   end select

end function calc_cache

!*************************************************
!
! Get a task; if none remain then end.
! If we can't find any work in the local task stack then we first try the
! task pool, and if this doesn't contain any, we steal tasks from
! other caches.
!
subroutine get_task(stack, task, control, info, st)

   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(dagtask), intent(inout) :: task
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag
   integer, intent(inout) :: st ! stat parameter

   integer :: this_thread, cache

   this_thread = 0
!$ this_thread = omp_get_thread_num()
   cache = calc_cache(this_thread, control)

!$OMP FLUSH(stack)
   !
   ! Check for error termination (no lock required)
   !
   if(stack%abort) then
      ! Error abort
      task%task_type = TASK_DONE
      if(task%task_type.ne.TASK_NONE) then
!$       call omp_set_lock(stack%lock)
         ! decrease count of acive threads
         stack%active = stack%active - 1
!$       call omp_unset_lock(stack%lock)
      endif
!$OMP FLUSH
      return
   endif

!$ call omp_set_lock(stack%clocks(cache))
   if(stack%cheads(cache).ne.0) then
      ! local task stack contains tasks ... take the one off the top and return
      task = stack%ctasks(stack%cheads(cache),cache)
      stack%cheads(cache) = stack%cheads(cache) - 1
!$    call omp_unset_lock(stack%clocks(cache))
      return
   endif
!$ call omp_unset_lock(stack%clocks(cache))

!$ call omp_set_lock(stack%lock)

   ! reduce count of number of active threads
   stack%active = stack%active - 1

   ! If no task in local task stack then we /must/ get one.
   ! First attempt to take a task from the task pool.
   ! If this pool is empty, search for largest local stack
   ! belonging to another cache. If found, move tasks in bottom half 
   ! from this local stack to the task pool (workstealing).
   ! Then take task of highest priority from the pool as next task.
 
   task%task_type = TASK_NONE
   do while(task%task_type.eq.TASK_NONE)
      call get_task_from_pool(stack, task)

      if(info.lt.0) then
!$       call omp_unset_lock(stack%lock)
         return
      endif
      if(task%task_type.eq.TASK_NONE) then
         ! Check if we've finished
         if(stack%active.eq.0) then
!$          call omp_unset_lock(stack%lock)
            task%task_type = TASK_DONE
            return
         endif

!$       call omp_unset_lock(stack%lock)

         ! Spin until a task become available
         call await_task(stack,control)
         ! immediate return if we have to abort
         if(stack%abort) return

         if(stack%cheads(cache).ne.0) then
            ! tasks available in local task stack. Take one from the top.
!$          call omp_set_lock(stack%clocks(cache))
            if(stack%cheads(cache).ne.0) then
               task = stack%ctasks(stack%cheads(cache),cache)
               stack%cheads(cache) = stack%cheads(cache) - 1
!$             call omp_unset_lock(stack%clocks(cache))
!$             call omp_set_lock(stack%lock)
               exit
            endif
!$          call omp_unset_lock(stack%clocks(cache))
         endif

         if(stack%total.le.0) then
            ! nothing left in task pool so look to steal tasks
            ! from another thread. 
            call worksteal(stack, control, info, st)
            if(info.lt.0) return
         endif
!$       call omp_set_lock(stack%lock)
         cycle
      endif
   enddo

   if(task%task_type.ne.TASK_DONE) stack%active = stack%active + 1

!$ call omp_unset_lock(stack%lock)

end subroutine get_task

!*************************************************
!
! Look to other caches to steal work from.
! Steals from the largest local stack. Moves bottom half
! of this stack into the task pool and then moves the
! remaining tasks down the stack
!
subroutine worksteal(stack, control, info, st)

   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(MA87_control), intent(in) :: control
   integer, intent(inout) :: info ! error flag
   integer, intent(inout) :: st ! stat parameter

   integer :: i, j
   integer :: mi ! index oflargest local stack (-1 if all are empty)
   integer :: msz ! size of largest locak stack
   integer :: cache, total_threads, this_thread

   total_threads = 1
!$ total_threads = omp_get_num_threads()
   this_thread = 0
!$ this_thread = omp_get_thread_num()
   cache = calc_cache(this_thread-1, control)

   ! Find biggest stack to steal from
   msz = -1
   mi = -1
   do i = 1, calc_cache(total_threads-1, control)
      if(cache.eq.i) cycle
      if(stack%cheads(i).gt.msz) then
!$       call omp_set_lock(stack%clocks(i))
         if(stack%cheads(i).gt.msz) then
!$          if(mi.ne.-1) call omp_unset_lock(stack%clocks(mi))
            mi = i
            msz = stack%cheads(i)
         else ! Its no longer bigger, release lock
!$          call omp_unset_lock(stack%clocks(i))
         endif
      endif
   end do
   if(mi.eq.-1) return ! No other caches

   msz = stack%cheads(mi)

   ! Steal half from bottom of the stack mi
!$ call omp_set_lock(stack%lock)
   do i = 1, msz / 2
      call add_task_g(stack, stack%ctasks(i,mi), &
         control, info, st, locked=.true.)
      if(info.lt.0) then
!$       call omp_unset_lock(stack%lock)
         return
      endif
   end do
!$ call omp_unset_lock(stack%lock)
   ! move the remaining tasks down to occupied freed up space
   j = 1
   do i = msz / 2 + 1, msz
      stack%ctasks(j,mi) = stack%ctasks(i,mi)
      j = j + 1
   end do
   stack%cheads(mi) = stack%cheads(mi) - msz / 2
!$ call omp_unset_lock(stack%clocks(mi))
end subroutine worksteal

!*************************************************
!
! Get a task from the task pool. Want a task with as low a priority
! value as possible.
!
subroutine get_task_from_pool(stack, task)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(dagtask), intent(inout) :: task

   integer :: priority ! priority of task
   integer :: ptr ! pointer into stack of tasks
   !
   ! Check there exists a task for us to get
   !
   if(stack%total.le.0) then

      task%task_type = TASK_NONE
      return
   else

      ! Find a task with the lowest priority value (that is, the
      ! task which has the highest priority)
      ! start with the current lowest priority value and increase priority
      ! until a non empty linked list of tasks is found

      priority = stack%lowest_priority_value
      do while(stack%prihead(priority).eq.-1)
         priority = priority + 1
      end do
      stack%lowest_priority_value = priority

      stack%max_pool_size = max(stack%max_pool_size,stack%total)

      ! Grab a task from the stack of tasks
      !
      ptr = stack%prihead(priority)
      task = stack%tasks(ptr)
      stack%prihead(priority) = stack%next(ptr)
      stack%next(ptr) = stack%freehead
      stack%freehead = ptr

      ! reduce count of tasks in stack
      stack%total = stack%total - 1

   endif
end subroutine get_task_from_pool
 
!*************************************************
!
! Spin until either some work crops up or notification of abort
! received or all work has been performed so nothing to wait for
!
subroutine await_task(stack, control)

   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool
   type(MA87_control), intent(in) :: control

   integer :: cache, this_thread
  !** integer :: t_start, t_end, t_rate

   this_thread = 0
!$ this_thread = omp_get_thread_num()
   cache = calc_cache(this_thread, control)

 !**  call system_clock(t_start, t_rate)

!$OMP CRITICAL (await_task_idle)
   do
!$OMP FLUSH(stack)
      if(stack%abort .or. stack%total.gt.0 .or. &
         (stack%total.eq.0 .and. stack%active.eq.0) .or. &
          any(stack%cheads(:) .gt. 4) .or. stack%cheads(cache).ne.0) exit
   end do
!$OMP END CRITICAL (await_task_idle)

 !**  if (control%time_out.ge.0) then
 !**    call system_clock(t_end)
 !**    stack%waiting(this_thread) = &
 !**       stack%waiting(this_thread) + (t_end-t_start)/real(t_rate)
 !**  end if

end subroutine await_task

!*************************************************
!
! Sets the abort flag, called on error
!
subroutine set_abort(stack)
   type(taskstack), intent(inout) :: stack ! holds local stacks and task pool

   ! Set the state
!$ call omp_set_lock(stack%lock)
   stack%abort = .true.
!$ call omp_unset_lock(stack%lock)
end subroutine set_abort

!*************************************************
!
! TASK_BSOLV, TASK_FSOLV
! B_j <- L_jj^-1 B_j
! B_j <- L_jj^-T B_j
!
! Note: While diagonal blocks may be trapezoidal, this is handled at the
! level calling this subroutine through a call to slv_fwd_update or
! slv_bwd_update

subroutine slv_solve(n, nelim, col, dest, trans, unit,  &
     nrhs, rhs, ldr, control, id)

   integer, intent(in) :: n ! leading dimension of diag. block
   integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
   integer, intent(in) :: col ! start of block column variables in rhs
   real(wp), dimension(*), intent(in) :: dest ! holds destination block
   character(len=13), intent(in) :: trans ! set to 
      ! 'Transpose    ' for forward substitution and to 
      ! 'Non-Transpose' for back substitution
   character(len=8), intent(in) :: unit ! set to 
      ! 'Non-Unit' for positive-definite case
   integer, intent(in) :: nrhs ! number of right-hand sides
   integer, intent(in) :: ldr ! leading extent of rhs
   real(wp), intent(inout) :: rhs(ldr*nrhs)
   type(MA87_control), intent(in) :: control
   integer(long), intent(in) :: id

!%%%   integer :: this_thread
!%%%   integer :: t_start, t_end

   if(nelim.eq.0) return

!%%%   if(control%unit_log.gt.0) call system_clock(t_start)

   if(nrhs.eq.1) then
      call dtrsv('Upper', trans, unit, nelim, dest, n, rhs(col), 1)
   else
      call dtrsm('Left', 'Upper', trans, unit, nelim, nrhs, one, &
         dest, n, rhs(col), ldr)
   endif

!%%%   if(control%unit_log.gt.0) then
!%%%      call system_clock(t_end)
!%%%      this_thread = 0
!%%% !$    this_thread = omp_get_thread_num()
!%%%      if(trans(1:1).eq.'T') then
!%%%         call log_task(control, this_thread, t_start, t_end, "FS", id)
!%%%      else
!%%%         call log_task(control, this_thread, t_start, t_end, "BS", id)
!%%%      endif
!%%%   endif

end subroutine slv_solve

!*************************************************
!
! TASK_FUPD
! B_j <- B_j - L_ij B_i
!
subroutine slv_fwd_update(m, nelim, col, offset, index, dest, ldd, nrhs, &
      upd, ldu, rhs, ldr, xlocal, control, id)

   integer, intent(in) :: m ! number of rows in block
   integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
   integer, intent(in) :: col ! start of block column variables in rhs
   integer, intent(in) :: offset ! offset into index we start at
   integer, dimension(*), intent(in) :: index
   integer, intent(in) :: ldd ! leading dimension of block
   real(wp), dimension(m*ldd), intent(in) :: dest ! holds destination block
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldu  ! leading extent of upd
   real(wp), intent(inout) :: upd(ldu*nrhs) ! vector to update
   integer, intent(in) :: ldr  ! leading extent of rhs
   real(wp), intent(in) :: rhs(ldr*nrhs) ! rhs vector
   real(wp), dimension(*), intent(out) :: xlocal
   type(MA87_control), intent(in) :: control
   integer(long), intent(in) :: id

   integer :: i
   integer :: j
   integer :: k
   integer :: r ! right hand side loop variable
   real(wp) :: w ! temporary work value
 !%%%  integer :: t_start, t_end, this_thread

   if(nelim.eq.0) return


  !%%% if(control%unit_log.gt.0) call system_clock(t_start)

   ! forward substitution
   if(nrhs.eq.1) then
      if(m-nelim.gt.10 .and. nelim.gt.4) then
         !!! Single rhs, BLAS 2

         call dgemv('T', nelim, m, -one, dest, ldd, rhs(col), 1, zero, &
            xlocal, 1)
   
         ! Copy xlocal out
         j = 1
         do i = offset, offset+m-1
            upd(index(i)) = upd(index(i)) + xlocal(j)
            j = j + 1
         end do
      else
         !!! Single rhs, direct update
         j = 1
         do i = offset, offset+m-1
            w = zero
            do k = col, col+nelim-1
               w = w - dest(j)*rhs(k)
               j = j + 1
            end do
            j = j + (ldd-nelim)
            upd(index(i)) = upd(index(i)) + w
         end do   
      endif
   else
      !!! Multiple rhs, BLAS 3
      call dgemm('T', 'N', m, nrhs, nelim, -one, dest, ldd, rhs(col), ldr, &
         zero, xlocal, m)
   
      ! Copy xlocal out
      j = 1
      do i = offset, offset+m-1
         do r = 0, nrhs-1
            upd(index(i)+r*ldu) = upd(index(i)+r*ldu) + xlocal(j+r*m)
         end do
         j = j + 1
      end do
   endif

!%%%   if(control%unit_log.gt.0) then
!%%%      this_thread = 0
!%%% !$    this_thread = omp_get_thread_num()
!%%%      call system_clock(t_end)
!%%%      call log_task(control, this_thread, t_start, t_end, "FU", id)
!%%%   endif
end subroutine slv_fwd_update

!*************************************************
!
! TASK_BUPD
! B_i <- B_i - L_ij^-T B_j
!
subroutine slv_bwd_update(m, nelim, col, offset, index, dest, ldd, nrhs, rhs, &
      upd, ldr, xlocal, control, id)
   integer, intent(in) :: m ! number of rows in block
   integer, intent(in) :: nelim ! number eliminations (immediate return if =0)
   integer, intent(in) :: col ! start of block column variables in rhs
   integer, intent(in) :: offset ! offset into index we start at
   integer, dimension(*), intent(in) :: index
   integer, intent(in) :: ldd ! leading dimension of block
   real(wp), dimension(m*ldd), intent(in) :: dest ! holds block
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldr  ! leading extent of rhs
   real(wp), intent(inout) :: rhs(ldr*nrhs)
   real(wp), intent(inout) :: upd(ldr*nrhs)
   real(wp), dimension(*), intent(out) :: xlocal
   type(MA87_control), intent(in) :: control
   integer(long), intent(in) :: id

   integer :: i
   integer :: j
   integer :: k
   integer :: r ! right hand side loop variable
   real(wp) :: w ! temporary work variable
 !%%%  integer :: t_start, t_end, this_thread

   if(nelim.eq.0) return

 !%%%  if(control%unit_log.gt.0) call system_clock(t_start)

   ! backward substitution
   if(nrhs.eq.1) then
      if(m-nelim.gt.10 .and. nelim.gt.4) then
         !!! Single right-hand side, BLAS 2

         ! Copy xlocal in
         j = 1
         do i = offset, offset+m-1
            xlocal(j) = rhs(index(i))
            j = j + 1
         end do

         call dgemv('N', nelim, m, -one, dest, ldd, xlocal, 1, one, &
            upd(col), 1)
      else
         !!! Single right-hand side, direct update
         j = 1
         do i = offset, offset+m-1
            w = rhs(index(i))
            do k = col, col + nelim - 1
               upd(k) = upd(k) - dest(j)*w
               j = j + 1
            end do
            j = j + (ldd-nelim)
         end do
      endif
   else
      !!! Multiple RHS, BLAS 3

      ! Copy xlocal in
      j = 1
      do i = offset, offset+m-1
         do r = 0, nrhs-1
            xlocal(j+r*m) = rhs(index(i)+r*ldr)
         end do
         j = j + 1
      end do

      call dgemm('N', 'N', nelim, nrhs, m, -one, dest, ldd, xlocal, m, &
         one, upd(col), ldr)
   endif

!%%%   if(control%unit_log.gt.0) then
!%%%      this_thread = 0
!%%% !$    this_thread = omp_get_thread_num()
!%%%      call system_clock(t_end)
!%%%      call log_task(control, this_thread, t_start, t_end, "BU", id)
!%%%   endif

end subroutine slv_bwd_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debugging, logging and error handling routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Log a task
!
!%%%subroutine log_task(control, thread, start, finish, task, opt1, opt2, opt3)
!%%%    type(ma87_control), intent(in) :: control
!%%%   integer, intent(in) :: thread
!%%%   integer, intent(in) :: start
!%%%   integer, intent(in) :: finish
!%%%   character(len=2), intent(in) :: task
!%%%   integer(long), optional, intent(in) :: opt1
!%%%   integer(long), optional, intent(in) :: opt2
!%%%   integer(long), optional, intent(in) :: opt3

!%%%   integer :: arg

!%%%   arg = 0
!%%%   if(present(opt1)) arg = arg + 1
!%%%   if(present(opt2)) arg = arg + 1
!%%%   if(present(opt3)) arg = arg + 1

!%%%   select case(arg)
!%%%   case(0)
!%%%      write(control%unit_log, "(i4, 2i12, 1x, a2)") &
!%%%         thread, start, finish, task
!%%%   case(1)
!%%%      write(control%unit_log, "(i4, 2i12, 1x, a2, 1x, 1i12)") &
!%%%         thread, start, finish, task, opt1
!%%%   case(2)
!%%%      write(control%unit_log, "(i4, 2i12, 1x, a2, 1x, 2i12)") &
!%%%         thread, start, finish, task, opt1, opt2
!%%%   case(3)
!%%%      write(control%unit_log, "(i4, 2i12, 1x, a2, 1x, 3i12)") &
!%%%         thread, start, finish, task, opt1, opt2, opt3
!%%%   end select
!%%%end subroutine log_task

!*************************************************
!
! Converts a task into a character representation thereof
! Beware not to cause nested i/O !!!
!
!function print_job(task)
!   character(len=26) :: print_job
!   type(dagtask), intent(in) :: task
!
!   select case(task%task_type)
!   case(TASK_FACTORIZE_BLOCK)
!      write(print_job, "(a,'(',i5,')')") "factorize_block", task%dest
!   case(TASK_SOLVE_BLOCK)
!      write(print_job, "(a,'(',i5,',',i5,')')") "solve_block", &
!      task%dest, task%src1
!   case(TASK_UPDATE_INTERNAL)
!      write(print_job, "(a,'(',i5,2(',',i5),')')") "within", task%dest, &
!         task%src1, task%src2
!   case(TASK_UPDATE_BETWEEN)
!      write(print_job, "(a,'(',i5,',',i5,')')") "between", &
!         task%dest, task%src1
!   case(TASK_SLV_FSLV)
!      write(print_job, "(a,'(',i5,')')") "slv_fslv", &
!         task%dest
!   case(TASK_SLV_FUPD)
!      write(print_job, "(a,'(',i5,')')") "slv_fupd", &
!         task%dest
!   case(TASK_SLV_BSLV)
!      write(print_job, "(a,'(',i5,')')") "slv_bslv", &
!         task%dest
!   case(TASK_SLV_BUPD)
!      write(print_job, "(a,'(',i5,')')") "slv_bupd", &
!         task%dest
!   case(TASK_DONE)
!      write(print_job, "(a)") "FINISH"
!   case(TASK_NONE)
!      write(print_job, "(a)") "WAIT"
!   case default
!      write(print_job, "(a)") "UNKNOWN TASK"
!   end select
!end function print_job

!*************************************************

! printing of error messages

subroutine ma87_print_flag(iflag, control, context, st)

   integer, intent(in) :: iflag
   type(ma87_control), intent(in) :: control
   integer, intent(in), optional :: st
   ! context: is an optional assumed size character array of intent(in).
   ! It describes the context under which the error occured
   character (len=*), optional, intent(in) :: context

   integer :: nout

      nout = control%unit_error
      if (control%diagnostics_level < 0) nout = -1
      if (nout < 0) return
      write (nout,'(/3a,i3)') ' Error return from ',trim(context),&
         '. Error flag = ', iflag


   select case(iflag)
   case(MA87_ERROR_ALLOCATION)
      if (present(st)) write (nout,'(a,i8)') &
         ' Allocation error. stat parameter = ', st
   case(MA87_ERROR_JOB_OOR)
      write (nout,'(a)') ' job out of range.'
   case(MA87_ERROR_NBI_OOR)
      write (nout,'(a)') ' nbi out of range.'
   case(MA87_ERROR_ORDER)
      write (nout,'(a)') ' Error in user-supplied elimination order'
   case(MA87_ERROR_X_SIZE)
      write (nout,'(a)') ' Error in size of x. lx or nrhs too small'
   case(MA87_ERROR_INFINITY)
      write (nout,'(a)') ' IEEE infinities found in factorization'

   ! Unknown error
   case default
      write (nout,'(a)') ' Unexpected Error. Please report.'
   end select

end subroutine MA87_print_flag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routines to help interfaces get at private data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure integer function ma87_get_n_double(keep)
   type(ma87_keep), intent(in) :: keep
   ma87_get_n_double = keep%n
end function ma87_get_n_double

end module hsl_MA87_double
