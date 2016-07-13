module spllt_data_mod
#if defined(SPLLT_USE_STARPU)
  use iso_c_binding
  use starpu_f_mod
#elif defined(SPLLT_USE_PARSEC)
  use dague_f08_interfaces
#endif
  use hsl_ma87_double, only: block_type, node_type 
  use hsl_zd11_double
  implicit none

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

  ! user control
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
     ! data descriptor
     type(c_ptr) :: ddesc
#endif

#if defined(SPLLT_USE_OMP)
     type(spllt_bc_type), allocatable :: row_list(:), col_list(:) ! workspaces
#else
     type(spllt_bc_type) :: row_list, col_list ! workspaces
#endif
  end type spllt_data_type

  type spllt_options
     integer :: ncpu = 1 ! number of CPU workers
     integer :: nb   = 16 ! blocking size
     character(len=100) :: mat = ''
     integer :: nemin = -1
  end type spllt_options

end module spllt_data_mod
