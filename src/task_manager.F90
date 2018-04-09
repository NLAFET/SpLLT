module task_manager_mod
  use spllt_data_mod
  implicit none

  integer, parameter :: ntrace_id = 12
  character(len=12), parameter :: task_manager_trace_names(ntrace_id) = &
    [character(len=10) :: &
      "INIT_NODE",        &
      "FACTO_BLK",        &
      "SOLVE_BLK",        &
      "UPDATE_BLK",       &
      "UPDATE_BTW",       &
      "fwd_update",       &
      "fwd_block",        &
      "bwd_update",       &
      "bwd_block",        &
      "fwd_submit",       &
      "bwd_submit",       &
      "chk_err" ]
  integer, parameter :: trace_init_node_pos   =  1
  integer, parameter :: trace_facto_blk_pos   =  2
  integer, parameter :: trace_solve_blk_pos   =  3
  integer, parameter :: trace_update_blk_pos  =  4
  integer, parameter :: trace_update_btw_pos  =  5
  integer, parameter :: trace_fwd_update_pos  =  6
  integer, parameter :: trace_fwd_block_pos   =  7
  integer, parameter :: trace_bwd_update_pos  =  8
  integer, parameter :: trace_bwd_block_pos   =  9
  integer, parameter :: trace_fwd_submit_pos  = 10
  integer, parameter :: trace_bwd_submit_pos  = 11
  integer, parameter :: trace_chk_err_pos     = 12

  type, abstract :: task_manager_base
    integer           :: nworker
    integer           :: masterWorker
    integer           :: workerID
    integer, pointer  :: trace_ids(:)
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
    ! Submission of a forward block task by the task manager
    !
    subroutine solve_fwd_block_task_iface(task_manager, dblk, nrhs, upd, &
      rhs, ldr, xlocal, fkeep)
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

    end subroutine solve_fwd_block_task_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a forward update task by the task manager
    !
    subroutine solve_fwd_update_task_iface(task_manager, blk, node, nrhs, upd, &
        rhs, ldr, xlocal, fkeep)
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
    end subroutine solve_fwd_update_task_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a backward block task by the task manager
    !
    subroutine solve_bwd_block_task_iface(task_manager, dblk, nrhs, upd, rhs, &
        ldr, xlocal, fkeep)
      use spllt_data_mod
      import task_manager_base
      class(task_manager_base),   intent(inout) :: task_manager
      integer, intent(in)                       :: dblk
      integer, intent(in)                       :: nrhs
      integer, intent(in)                       :: ldr
      real(wp), target, intent(inout)           :: upd(:, :)
      real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
      real(wp), target, intent(inout)           :: xlocal(:, :)
      type(spllt_fkeep), target, intent(in)     :: fkeep
    end subroutine solve_bwd_block_task_iface

    !!!!!!!!!!!!!!!!!!!!!
    ! Submission of a backward update task by the task manager
    !
    subroutine solve_bwd_update_task_iface(task_manager, blk, node, nrhs, upd, &
        rhs, ldr, xlocal, fkeep)
      use spllt_data_mod
      import task_manager_base
      class(task_manager_base),   intent(inout) :: task_manager
      integer, intent(in)                       :: blk
      integer, intent(in)                       :: node 
      integer, intent(in)                       :: nrhs
      integer, intent(in)                       :: ldr
      real(wp), target, intent(inout)           :: upd(:,:)
      real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
      real(wp), target, intent(inout)           :: xlocal(:,:)
      type(spllt_fkeep), target, intent(in)     :: fkeep
    end subroutine solve_bwd_update_task_iface

  end interface

contains

end module task_manager_mod