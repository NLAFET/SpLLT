module worker_info_mod
  
  type worker_info_t
    double precision  :: nflop              = 0.0 ! # flop performed
    double precision  :: nflop_performed    = 0.0 ! # flop performed locally
    integer           :: ntask_run          = 0   ! #task run by this worker
    integer           :: ntask_insert       = 0   ! #task insert to the runtim
    integer           :: nfake_task_insert  = 0   ! #fake task insert
    integer           :: narray_allocated   = 0   ! #allocation
  contains
    procedure :: print  => print_info
    procedure :: init   => init_worker
    procedure :: reset  => reset_worker
  end type worker_info_t

contains

  subroutine print_info(self, msg, worker_id)
    implicit none
    class(worker_info_t), intent(in)  :: self
    character(len=*),     intent(in)  :: msg
    integer, optional,    intent(in)  :: worker_id

    integer :: id

    id = -1

    if(present(worker_id)) id = worker_id

    print *, msg, " : ", id
    print '(a, i9)',      "#task insert        : ", self%ntask_insert
    print '(a, i9)',      "#fake task insert   : ", self%nfake_task_insert
    print '(a, i9)',      "#task run           : ", self%ntask_run
    print '(a, i9)',      "#array allocate     : ", self%narray_allocated
    print '(a, es10.2)',  "#flop               : ", self%nflop
    print '(a, es10.2)',  "#flop_performed     : ", self%nflop_performed

  end subroutine print_info



  subroutine init_worker(self)
    implicit none
    class(worker_info_t), intent(out) :: self

    call self%reset()
    
  end subroutine init_worker



  subroutine reset_worker(self)
    implicit none
    class(worker_info_t), intent(out) :: self

    self%nflop               = 0.0
    self%nflop_performed     = 0.0
    self%ntask_run           = 0
    self%ntask_insert        = 0
    self%nfake_task_insert   = 0
    self%narray_allocated    = 0
    
  end subroutine reset_worker



end module worker_info_mod
