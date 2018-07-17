module task_manager_mod
  use spllt_data_mod
  use worker_info_mod
  implicit none

  integer, parameter :: ntrace_id = 16
  character(len=ntrace_id), parameter :: task_manager_trace_names(ntrace_id) = &
    [character(len=20) :: &
      "INIT_NODE"  ,      &
      "FACTO_BLK"  ,      &
      "SOLVE_BLK"  ,      &
      "UPDATE_BLK" ,      &
      "UPDATE_BTW" ,      &
      "fwd_update" ,      &
      "fwd_block"  ,      &
      "bwd_update" ,      &
      "bwd_block"  ,      &
      "fwd_submit" ,      &
      "bwd_submit" ,      &
      "chk_err"    ,      &
      "fwd_submitTree",   &
      "bwd_submitTree",   &
      "fwd_subtree",      &
      "bwd_subtree"       &
      ]
  integer, parameter :: no_trace                  =  0
  integer, parameter :: trace_init_node_pos       =  1
  integer, parameter :: trace_facto_blk_pos       =  2
  integer, parameter :: trace_solve_blk_pos       =  3
  integer, parameter :: trace_update_blk_pos      =  4
  integer, parameter :: trace_update_btw_pos      =  5
  integer, parameter :: trace_fwd_update_pos      =  6
  integer, parameter :: trace_fwd_block_pos       =  7
  integer, parameter :: trace_bwd_update_pos      =  8
  integer, parameter :: trace_bwd_block_pos       =  9
  integer, parameter :: trace_fwd_submit_pos      = 10
  integer, parameter :: trace_bwd_submit_pos      = 11
  integer, parameter :: trace_chk_err_pos         = 12
  integer, parameter :: trace_fwd_submit_tree_pos = 13
  integer, parameter :: trace_bwd_submit_tree_pos = 14
  integer, parameter :: trace_fwd_subtree_pos     = 15
  integer, parameter :: trace_bwd_subtree_pos     = 16

  type, abstract :: task_manager_base
    integer                       :: nworker
    integer                       :: masterWorker
    integer                       :: workerID
    integer, pointer              :: trace_ids(:)
    type(worker_info_t), pointer  :: worker_info(:)
  contains
    procedure(task_manager_init_iface),           deferred :: init
    procedure(task_manager_reset_iface),          deferred :: reset
    procedure(task_manager_print_iface),          deferred :: print
    procedure(task_manager_incr_alloc_iface),     deferred :: incr_alloc
    procedure(task_manager_refresh_master_iface), deferred :: refresh_master
    procedure(task_manager_refresh_worker_iface), deferred :: refresh_worker
    procedure(task_manager_deallocate_iface),     deferred :: deallocate


    procedure(solve_fwd_block_task_iface),  deferred :: solve_fwd_block_task
    procedure(solve_fwd_update_task_iface), deferred :: solve_fwd_update_task
    procedure(solve_bwd_block_task_iface),  deferred :: solve_bwd_block_task
    procedure(solve_bwd_update_task_iface), deferred :: solve_bwd_update_task
    procedure(solve_fwd_subtree_task_iface),deferred :: solve_fwd_subtree_task
    procedure(solve_bwd_subtree_task_iface),deferred :: solve_bwd_subtree_task

    procedure(task_manager_get_nflop_iface),deferred :: get_nflop_performed

    procedure(solve_fwd_block_task_il_iface),  deferred :: solve_fwd_block_il_task
    procedure(solve_fwd_update_task_il_iface), deferred :: solve_fwd_update_il_task
    procedure(solve_bwd_block_task_il_iface),  deferred :: solve_bwd_block_il_task
    procedure(solve_bwd_update_task_il_iface), deferred :: solve_bwd_update_il_task
    procedure(solve_fwd_subtree_task_il_iface),deferred :: solve_fwd_subtree_il_task
    procedure(solve_bwd_subtree_task_il_iface),deferred :: solve_bwd_subtree_il_task

    procedure(solve_fwd_block_task_il2_iface),  deferred :: solve_fwd_block_il2_task
    procedure(solve_fwd_update_task_il2_iface), deferred :: solve_fwd_update_il2_task
    procedure(solve_bwd_block_task_il2_iface),  deferred :: solve_bwd_block_il2_task
    procedure(solve_bwd_update_task_il2_iface), deferred :: solve_bwd_update_il2_task
    procedure(solve_fwd_subtree_task_il2_iface),deferred :: solve_fwd_subtree_il2_task
    procedure(solve_bwd_subtree_task_il2_iface),deferred :: solve_bwd_subtree_il2_task
  end type task_manager_base

  abstract interface 
    !!!!!!!!!!!!!!!!!!!!!
    ! Initialize the task_manager
    !
    subroutine task_manager_init_iface(self, trace_names, stat)
      import task_manager_base
      class(task_manager_base), target,     intent(inout) :: self
      character(len=*), optional,         intent(in)    :: trace_names(:)
      integer,          optional,         intent(out)   :: stat
    end subroutine task_manager_init_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Reset variables
    !
    subroutine task_manager_reset_iface(self)
      import task_manager_base
      class(task_manager_base), intent(inout) :: self
    end subroutine task_manager_reset_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Refresh master id
    !
    subroutine task_manager_refresh_master_iface(self)
      import task_manager_base
      class(task_manager_base), intent(inout) :: self
    end subroutine task_manager_refresh_master_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Refresh current worker id
    !
    subroutine task_manager_refresh_worker_iface(self)
      import task_manager_base
      class(task_manager_base), intent(inout) :: self
    end subroutine task_manager_refresh_worker_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Print the state of the task manager
    !
    subroutine task_manager_print_iface(self, msg, option)
      import task_manager_base
      class(task_manager_base), intent(in)  :: self
      character(len=*), optional            :: msg
      integer, optional                     :: option
    end subroutine task_manager_print_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Increment the number of allocation counted by the task manager
    !
    subroutine task_manager_incr_alloc_iface(self, stat)
      import task_manager_base
      class(task_manager_base), intent(inout) :: self
      integer,                  intent(in)    :: stat
    end subroutine task_manager_incr_alloc_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Deallocate the task manager
    !
    subroutine task_manager_deallocate_iface(self)
      import task_manager_base
      class(task_manager_base), intent(inout) :: self
    end subroutine task_manager_deallocate_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Get number of flop performed at that point
    !
    subroutine task_manager_get_nflop_iface(self, nflop, thn)
      import task_manager_base
      class(task_manager_base), intent(inout) :: self
      double precision,         intent(out)   :: nflop
      integer, optional,        intent(in)    :: thn
    end subroutine task_manager_get_nflop_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a forward block task by the task manager
    !
    subroutine solve_fwd_block_task_iface(task_manager, dblk, nrhs, upd, &
      rhs, ldr, xlocal, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: dblk 
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: ldr
      real(wp), target,           intent(inout) :: upd(:, :)
      real(wp), target,           intent(inout) :: rhs(ldr * nrhs)
      real(wp), target,           intent(inout) :: xlocal(:, :)
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      integer, optional,          intent(in)    :: trace_id

    end subroutine solve_fwd_block_task_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a forward update task by the task manager
    !
    subroutine solve_fwd_update_task_iface(task_manager, blk, node, nrhs, upd, &
        rhs, ldr, xlocal, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: blk
      integer,                    intent(in)    :: node
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: ldr
      real(wp), target,           intent(inout) :: upd(:,:)        
      real(wp), target,           intent(in)    :: rhs(ldr*nrhs)
      real(wp), target,           intent(out)   :: xlocal(:,:)
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      integer, optional,          intent(in)    :: trace_id
    end subroutine solve_fwd_update_task_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a backward block task by the task manager
    !
    subroutine solve_bwd_block_task_iface(task_manager, dblk, nrhs, upd, rhs, &
        ldr, xlocal, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: dblk
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: ldr
      real(wp), target,           intent(inout) :: upd(:, :)
      real(wp), target,           intent(inout) :: rhs(ldr * nrhs)
      real(wp), target,           intent(inout) :: xlocal(:, :)
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      integer, optional,          intent(in)    :: trace_id
    end subroutine solve_bwd_block_task_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a backward update task by the task manager
    !
    subroutine solve_bwd_update_task_iface(task_manager, blk, node, nrhs, upd, &
        rhs, ldr, xlocal, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: blk
      integer,                    intent(in)    :: node 
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: ldr
      real(wp), target,           intent(inout) :: upd(:,:)
      real(wp), target,           intent(inout) :: rhs(ldr * nrhs)
      real(wp), target,           intent(inout) :: xlocal(:,:)
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      integer, optional,          intent(in)    :: trace_id
    end subroutine solve_bwd_update_task_iface

    subroutine solve_fwd_subtree_task_iface(task_manager, nrhs, rhs, ldr, &
        fkeep, tree, xlocal, rhs_local)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base ),  intent(inout) :: task_manager
      integer,                    intent(in)    :: nrhs ! Number of RHS
      real(wp),                   intent(inout) :: rhs(ldr*nrhs)
      integer,                    intent(in)    :: ldr  ! Leading dimension 
                                                        ! of RHS
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      type(spllt_tree_t),         intent(in)    :: tree
      real(wp),                   intent(inout) :: xlocal(:,:)
      real(wp),                   intent(inout) :: rhs_local(:,:)
    end subroutine solve_fwd_subtree_task_iface

    subroutine solve_bwd_subtree_task_iface(task_manager, nrhs, rhs, ldr, &
        fkeep, tree, xlocal, rhs_local)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base ),  intent(inout) :: task_manager
      integer,                    intent(in)    :: nrhs ! Number of RHS
      real(wp),                   intent(inout) :: rhs(ldr*nrhs)
      integer,                    intent(in)    :: ldr  ! Leading dimension 
                                                        ! of RHS
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      type(spllt_tree_t),         intent(in)    :: tree
      real(wp),                   intent(inout) :: xlocal(:,:)
      real(wp),                   intent(inout) :: rhs_local(:,:)
    end subroutine solve_bwd_subtree_task_iface

    !!!!!!!!!!!!!!!!!!!!!
    !   Interleave interfaces
    !!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a forward block task by the task manager
    !
    subroutine solve_fwd_block_task_il_iface(task_manager, dblk, nrhs, upd, &
      tdu, rhs, n, xlocal, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: dblk 
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: tdu
      integer,                    intent(in)    :: n
      real(wp), target,           intent(inout) :: upd(:)
      real(wp), target,           intent(inout) :: rhs(n * nrhs)
      real(wp), target,           intent(inout) :: xlocal(:, :)
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      integer, optional,          intent(in)    :: trace_id

    end subroutine solve_fwd_block_task_il_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a forward update task by the task manager
    !
    subroutine solve_fwd_update_task_il_iface(task_manager, blk, node, nrhs, &
        upd, tdu, rhs, n, xlocal, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: blk
      integer,                    intent(in)    :: node
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: tdu
      integer,                    intent(in)    :: n
      real(wp), target,           intent(inout) :: upd(:)
      real(wp), target,           intent(in)    :: rhs(n*nrhs)
      real(wp), target,           intent(out)   :: xlocal(:,:)
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      integer, optional,          intent(in)    :: trace_id
    end subroutine solve_fwd_update_task_il_iface

    subroutine solve_fwd_subtree_task_il_iface(task_manager, nrhs, rhs, n, &
        fkeep, tree, xlocal, rhs_local, tdu)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base ),  intent(inout) :: task_manager
      integer,                    intent(in)    :: nrhs ! Number of RHS
      integer,                    intent(in)    :: tdu
      integer,                    intent(in)    :: n
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      type(spllt_tree_t),         intent(in)    :: tree
      real(wp),                   intent(inout) :: rhs(n*nrhs)
      real(wp),                   intent(inout) :: xlocal(:,:)
      real(wp),                   intent(inout) :: rhs_local(:)
    end subroutine solve_fwd_subtree_task_il_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a backward block task by the task manager
    !
    subroutine solve_bwd_block_task_il_iface(task_manager, dblk, nrhs, upd, &
      tdu, rhs, n, xlocal, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: dblk 
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: tdu
      integer,                    intent(in)    :: n
      real(wp), target,           intent(inout) :: upd(:)
      real(wp), target,           intent(inout) :: rhs(n*nrhs)
      real(wp), target,           intent(inout) :: xlocal(:, :)
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      integer, optional,          intent(in)    :: trace_id

    end subroutine solve_bwd_block_task_il_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a forward update task by the task manager
    !
    subroutine solve_bwd_update_task_il_iface(task_manager, blk, node, nrhs, &
        upd, tdu, rhs, n, xlocal, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: blk
      integer,                    intent(in)    :: node
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: tdu
      integer,                    intent(in)    :: n
      real(wp), target,           intent(inout) :: upd(:)
      real(wp), target,           intent(inout) :: rhs(n*nrhs)
      real(wp), target,           intent(inout) :: xlocal(:,:)
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      integer, optional,          intent(in)    :: trace_id
    end subroutine solve_bwd_update_task_il_iface

    subroutine solve_bwd_subtree_task_il_iface(task_manager, nrhs, rhs, n, &
        fkeep, tree, xlocal, rhs_local, tdu)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base ),  intent(inout) :: task_manager
      integer,                    intent(in)    :: nrhs ! Number of RHS
      integer,                    intent(in)    :: tdu
      integer,                    intent(in)    :: n
      type(spllt_fkeep), target,  intent(in)    :: fkeep
      type(spllt_tree_t),         intent(in)    :: tree
      real(wp),                   intent(inout) :: rhs(n*nrhs)
      real(wp),                   intent(inout) :: xlocal(:,:)
      real(wp),                   intent(inout) :: rhs_local(:)
    end subroutine solve_bwd_subtree_task_il_iface

    subroutine solve_fwd_block_task_il2_iface(task_manager, dblk, nrhs, &
      n, rhs, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: dblk 
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: n
      real(wp), target,           intent(inout) :: rhs(n, nrhs)
      type(spllt_fkeep), target,  intent(inout) :: fkeep
      integer, optional,          intent(in)    :: trace_id

    end subroutine solve_fwd_block_task_il2_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a forward update task by the task manager
    !
    subroutine solve_fwd_update_task_il2_iface(task_manager, blk, node, nrhs, &
        n, rhs, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: blk
      integer,                    intent(in)    :: node
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: n
      real(wp), target,           intent(in)    :: rhs(n, nrhs)
      type(spllt_fkeep), target,  intent(inout) :: fkeep
      integer, optional,          intent(in)    :: trace_id
    end subroutine solve_fwd_update_task_il2_iface

    subroutine solve_bwd_block_task_il2_iface(task_manager, dblk, nrhs, &
      n, rhs, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: dblk 
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: n
      real(wp), target,           intent(inout) :: rhs(n, nrhs)
      type(spllt_fkeep), target,  intent(inout) :: fkeep
      integer, optional,          intent(in)    :: trace_id

    end subroutine solve_bwd_block_task_il2_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a forward update task by the task manager
    !
    subroutine solve_bwd_update_task_il2_iface(task_manager, blk, node, nrhs, &
        n, rhs, fkeep, trace_id)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base),   intent(inout) :: task_manager
      integer,                    intent(in)    :: blk
      integer,                    intent(in)    :: node
      integer,                    intent(in)    :: nrhs
      integer,                    intent(in)    :: n
      real(wp), target,           intent(in)    :: rhs(n, nrhs)
      type(spllt_fkeep), target,  intent(inout) :: fkeep
      integer, optional,          intent(in)    :: trace_id
    end subroutine solve_bwd_update_task_il2_iface

    subroutine solve_fwd_subtree_task_il2_iface(task_manager, nrhs, rhs, n, &
        fkeep, tree)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base ),  intent(inout) :: task_manager
      integer,                    intent(in)    :: nrhs ! Number of RHS
      integer,                    intent(in)    :: n
      type(spllt_fkeep), target,  intent(inout) :: fkeep
      type(spllt_tree_t),         intent(in)    :: tree
      real(wp),                   intent(inout) :: rhs(n, nrhs)
    end subroutine solve_fwd_subtree_task_il2_iface

    subroutine solve_bwd_subtree_task_il2_iface(task_manager, nrhs, rhs, n, &
        fkeep, tree)
      use spllt_data_mod
      import task_manager_base

      class(task_manager_base ),  intent(inout) :: task_manager
      integer,                    intent(in)    :: nrhs ! Number of RHS
      integer,                    intent(in)    :: n
      type(spllt_fkeep), target,  intent(inout) :: fkeep
      type(spllt_tree_t),         intent(in)    :: tree
      real(wp),                   intent(inout) :: rhs(n, nrhs)
    end subroutine solve_bwd_subtree_task_il2_iface

  end interface

contains

end module task_manager_mod
