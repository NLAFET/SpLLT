module spllt_data_ciface                                                       
  use iso_c_binding
  implicit none

  type, bind(C) :: spllt_options_t
     integer(C_INT)     :: print_level
     integer(C_INT)     :: nrhs
     integer(C_INT)     :: ncpu
     integer(C_INT)     :: nb
!    character(len=100) :: mat
     integer(C_INT)     :: nemin
     integer(C_INT)     :: prune_tree
!    character(len=3)   :: fmt
     integer(C_INT)     :: min_width_blas
     integer(C_INT)     :: nb_min
     integer(C_INT)     :: nb_max
     integer(C_INT)     :: nrhs_min
     integer(C_INT)     :: nrhs_max
     integer(C_INT)     :: nb_linear_comp
     integer(C_INT)     :: nrhs_linear_comp
     integer(C_INT)     :: chunk
  end type spllt_options_t

  type, bind(C) :: spllt_inform_t
!   !type(ssids_inform) :: ssids_inform
     integer(C_INT) :: flag
     integer(C_INT) :: maxdepth
     integer(C_INT) :: num_factor
     integer(C_INT) :: num_flops
     integer(C_INT) :: num_nodes
     integer(C_INT) :: stat
  end type spllt_inform_t

contains

  subroutine copy_options_in(coptions, foptions)
! subroutine copy_options_in(coptions, foptions, cindexed)
    use spllt_data_mod
    implicit none
    type(spllt_options_t),  intent(in)  :: coptions
    type(spllt_options),    intent(out) :: foptions
!   logical,                intent(out) :: cindexed

!   cindexed                  = (coptions%array_base .eq. 0)
    foptions%print_level      = coptions%print_level
    foptions%ncpu             = coptions%ncpu
    foptions%nb               = coptions%nb
   !foptions%mat              = coptions%mat
    foptions%nemin            = coptions%nemin
    foptions%prune_tree       = (coptions%prune_tree .ne. 0)
   !foptions%fmt              = coptions%fmt
    foptions%min_width_blas   = coptions%min_width_blas
    foptions%nb_min           = coptions%nb_min
    foptions%nb_max           = coptions%nb_max
    foptions%nrhs_min         = coptions%nrhs_min
    foptions%nrhs_max         = coptions%nrhs_max
    foptions%nb_linear_comp   = (coptions%nb_linear_comp .ne. 0)
    foptions%nrhs_linear_comp = (coptions%nrhs_linear_comp .ne. 0)
    foptions%chunk            = coptions%chunk

  end subroutine copy_options_in



  subroutine copy_inform_out(finform, cinform)
    use spllt_data_mod
    implicit none
    type(spllt_inform),  intent(in)  :: finform
    type(spllt_inform_t),    intent(out) :: cinform

    cinform%flag       = finform%flag
    cinform%maxdepth   = finform%maxdepth
    cinform%num_factor = finform%num_factor
    cinform%num_flops  = finform%num_flops
    cinform%num_nodes  = finform%num_nodes
    cinform%stat       = finform%stat

  end subroutine copy_inform_out



end module spllt_data_ciface


  subroutine spllt_c_task_manager_init(ctask_manager) &
      bind(C, name="spllt_task_manager_init")
    use task_manager_omp_mod
    use spllt_data_ciface
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),  intent(inout) :: ctask_manager

    type(task_manager_omp_t),  pointer :: ftask_manager

    if(C_ASSOCIATED(ctask_manager)) then
      call C_F_POINTER(ctask_manager, ftask_manager)
    else
      allocate(ftask_manager)
      ctask_manager = C_LOC(ftask_manager)
    end if

    call ftask_manager%init()

  end subroutine spllt_c_task_manager_init



  subroutine spllt_c_wait() &
      bind(C, name="spllt_wait")
    use spllt_mod
    implicit none

    call spllt_wait()

  end subroutine spllt_c_wait



  subroutine spllt_c_analyse(cakeep, cfkeep, coptions, n, cptr, crow, &
      cinfo, corder) bind(C, name="spllt_analyse")
    use spllt_analyse_mod
    use spllt_data_mod
    use spllt_data_ciface
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),            intent(inout) :: cakeep
    type(C_PTR),            intent(inout) :: cfkeep
    type(spllt_options_t),  intent(in)    :: coptions
    integer(C_INT),         value         :: n
    type(C_PTR),            value         :: cptr
    type(C_PTR),            value         :: crow
    type(spllt_inform_t),   intent(out)   :: cinfo
    type(C_PTR),            value         :: corder

    type(spllt_akeep),  pointer :: fakeep
    type(spllt_fkeep),  pointer :: ffkeep
    type(spllt_options)         :: foptions
    integer(C_INT),     pointer :: fptr(:)
    integer(C_INT),     pointer :: frow(:)
    type(spllt_inform)          :: finfo
    integer(C_INT),     pointer :: forder(:)

    call copy_options_in(coptions, foptions)

    if(C_ASSOCIATED(cakeep)) then
      call C_F_POINTER(cakeep, fakeep)
    else
      allocate(fakeep)
      cakeep = C_LOC(fakeep)
    end if

    if(C_ASSOCIATED(cfkeep)) then
      call C_F_POINTER(cfkeep, ffkeep)
    else
      allocate(ffkeep)
      cfkeep = C_LOC(ffkeep)
    end if

    if(C_ASSOCIATED(cptr)) then
      call C_F_POINTER(cptr, fptr, shape=(/ n + 1 /))
    else
      nullify(fptr)
    end if

    if(C_ASSOCIATED(crow)) then
      call C_F_POINTER(crow, frow, shape=(/ fptr(n + 1) - 1 /))
    else
      nullify(frow)
    end if

    if(C_ASSOCIATED(corder)) then
      call C_F_POINTER(corder, forder, shape=(/ n /))
    else
      nullify(forder)
    end if

    call spllt_analyse(fakeep, ffkeep, foptions, n, fptr, frow, finfo, forder)

    call copy_inform_out(finfo, cinfo)
  end subroutine spllt_c_analyse



  subroutine spllt_c_factor(cakeep, cfkeep, coptions, nnz, cval, cinfo)  &
      bind(C, name="spllt_factor")
    use spllt_mod
    use spllt_data_mod
    use spllt_data_ciface
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),            value         :: cakeep
    type(C_PTR),            value         :: cfkeep
    type(spllt_options_t),  intent(in)    :: coptions
    integer(C_INT),         value         :: nnz
    type(C_PTR),            value         :: cval
    type(spllt_inform_t),   intent(out)   :: cinfo

    type(spllt_akeep),  pointer :: fakeep
    type(spllt_fkeep),  pointer :: ffkeep
    type(spllt_options)         :: foptions
    real(wp),           pointer :: fval(:)
    type(spllt_inform)          :: finfo

    call copy_options_in(coptions, foptions)

    if(.not. C_ASSOCIATED(cakeep)) then
      write (stderr,*) "Error, akeep provided by the user is empty"
    end if
    call C_F_POINTER(cakeep, fakeep)

    if(.not. C_ASSOCIATED(cfkeep)) then
      write (stderr,*) "Error, fkeep provided by the user is empty"
    end if
    call C_F_POINTER(cfkeep, ffkeep)

    if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape=(/ nnz /))
    else
      nullify(fval)
    end if

    call spllt_factor(fakeep, ffkeep, foptions, fval, finfo)
   !call spllt_wait()

    call copy_inform_out(finfo, cinfo)

  end subroutine spllt_c_factor



  subroutine spllt_c_prepare_solve(cakeep, cfkeep, nb, nrhs, cworksize, cinfo) &
      bind(C, name="spllt_prepare_solve")
    use spllt_data_mod
    use spllt_solve_dep_mod
    use spllt_data_ciface
    use timer_mod
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),            value         :: cakeep
    type(C_PTR),            value         :: cfkeep
    integer(C_INT),         value         :: nb
    integer(C_INT),         value         :: nrhs
    integer(C_LONG),        intent(out)   :: cworksize
    type(spllt_inform_t),   intent(out)   :: cinfo

    type(spllt_akeep),  pointer :: fakeep
    type(spllt_fkeep),  pointer :: ffkeep
    type(spllt_inform)          :: finfo
    integer                     :: st
    integer(long)               :: fworksize

    if(.not. C_ASSOCIATED(cakeep)) then
      write (stderr,*) "Error, akeep provided by the user is empty"
    end if
    call C_F_POINTER(cakeep, fakeep)

    if(.not. C_ASSOCIATED(cfkeep)) then
      write (stderr,*) "Error, fkeep provided by the user is empty"
    end if
    call C_F_POINTER(cfkeep, ffkeep)

    call spllt_init_timer(st)
    if(st .ne. 0) then
      write (stderr,*) "Warning, timer environment is not set"
    end if

   !print *, "Create subtrees"
    call spllt_create_subtree(fakeep, ffkeep)

   !print *, "Compute Sblocks with nb = ", nb, "nrhs = ", nrhs
    call get_solve_blocks(ffkeep, nb, nrhs, fworksize, ffkeep%sbc)

   !print *, "Compute dep"
    call spllt_compute_solve_dep(ffkeep, stat = st)

    cworksize = fworksize
    call copy_inform_out(finfo, cinfo)

  end subroutine spllt_c_prepare_solve



  subroutine spllt_c_set_mem_solve(cakeep, cfkeep, nb, nrhs, worksize, &
      cy, cworkspace, cinfo) bind(C, name="spllt_set_mem_solve")
    use spllt_data_mod
    use spllt_solve_dep_mod
    use spllt_data_ciface
    use timer_mod
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),            value         :: cakeep
    type(C_PTR),            value         :: cfkeep
    integer(C_INT),         value         :: nb
    integer(C_INT),         value         :: nrhs
    integer(C_LONG),        value         :: worksize
    type(C_PTR),            value         :: cy
    type(C_PTR),            value         :: cworkspace
    type(spllt_inform_t),   intent(out)   :: cinfo

    type(spllt_akeep),  pointer :: fakeep
    type(spllt_fkeep),  pointer :: ffkeep
    type(spllt_inform)          :: finfo
    real(wp),           pointer :: fy(:)   
    real(wp),           pointer :: fw(:)   
    integer                     :: st

    if(.not. C_ASSOCIATED(cakeep)) then
      write (stderr,*) "Error, akeep provided by the user is empty"
    end if
    call C_F_POINTER(cakeep, fakeep)

    if(.not. C_ASSOCIATED(cfkeep)) then
      write (stderr,*) "Error, fkeep provided by the user is empty"
    end if
    call C_F_POINTER(cfkeep, ffkeep)

    if(.not. C_ASSOCIATED(cy)) then
      write (stderr,*) "Error, y provided by the user is empty"
    end if
    call C_F_POINTER(cy, fy, shape=(/ ffkeep%n * nrhs /))

    if(.not. C_ASSOCIATED(cworkspace)) then
      write (stderr,*) "Error, workspace provided by the user is empty"
    end if
    call C_F_POINTER(cworkspace, fw, shape=(/ worksize /))

    call sblock_assoc_mem(ffkeep, nb, nrhs, fy, fw, ffkeep%sbc)

    call copy_inform_out(finfo, cinfo)
  end subroutine spllt_c_set_mem_solve



 !subroutine spllt_c_solve_workspace_size(cfkeep, nworker, nrhs, csize) &
 !    bind(C, name="spllt_solve_workspace_size")
 !  use spllt_data_mod
 !  use spllt_solve_mod
 !  use spllt_data_ciface
 !  use task_manager_omp_mod
 !  use ISO_Fortran_env, only: stderr => ERROR_UNIT
 !  implicit none

 !  type(C_PTR),            value         :: cfkeep
 !  integer(C_INT),         value         :: nworker
 !  integer(C_INT),         value         :: nrhs
 ! !type(C_PTR),            intent(out)   :: csize
 !  integer(C_LONG),        intent(out)   :: csize

 !  type(spllt_fkeep),  pointer :: ffkeep
 !  integer(long)               :: fsize

 !  if(.not. C_ASSOCIATED(cfkeep)) then
 !    write (stderr,*) "Error, fkeep provided by the user is empty"
 !  end if
 !  call C_F_POINTER(cfkeep, ffkeep)

 !  call spllt_solve_workspace_size(ffkeep, nworker, nrhs, fsize)
 !  csize = fsize

 !end subroutine spllt_c_solve_workspace_size


  subroutine spllt_c_solve(cfkeep, coptions, corder, nrhs, cx, cinfo, job) &
      bind(C, name="spllt_solve")
    use spllt_mod
    use spllt_data_mod
    use spllt_solve_mod
    use spllt_data_ciface
    use task_manager_omp_mod
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),            value         :: cfkeep
    type(spllt_options_t),  intent(in)    :: coptions
    type(C_PTR),            value         :: corder
    integer(C_INT),         value         :: nrhs
    type(C_PTR),            value         :: cx
    type(spllt_inform_t),   intent(out)   :: cinfo
    integer(C_INT),         value         :: job

    type(spllt_fkeep),  pointer :: ffkeep
    type(spllt_options)         :: foptions
    integer,            pointer :: forder(:)
    real(wp),           pointer :: fx(:)   
    type(spllt_inform)          :: finfo
    type(task_manager_omp_t)    :: task_manager

    if(.not. C_ASSOCIATED(cfkeep)) then
      write (stderr,*) "Error, fkeep provided by the user is empty"
    end if
    call C_F_POINTER(cfkeep, ffkeep)

    call copy_options_in(coptions, foptions)

    if(C_ASSOCIATED(corder)) then
      call C_F_POINTER(corder, forder, shape=(/ ffkeep%n /))
    else
      nullify(forder)
    end if

    if(.not. C_ASSOCIATED(cx)) then
      write (stderr,*) "Error, x/rhs provided by the user is empty"
    end if

   !call C_F_POINTER(cx, fx, shape=(/ ffkeep%n, nrhs/))
    call C_F_POINTER(cx, fx, shape=(/ ffkeep%n /))

    call task_manager%init()

    call spllt_solve_mult_double_worker(ffkeep, foptions, nrhs, fx, &
      job, task_manager, finfo)

    call spllt_wait()

    call task_manager%deallocate()

    call copy_inform_out(finfo, cinfo)

  end subroutine spllt_c_solve



  subroutine spllt_c_solve_worker(cfkeep, coptions, corder, nrhs, cx, &
      cinfo, job, cwork, cworksize, ctask_manager) &
      bind(C, name="spllt_solve_worker")
    use spllt_mod
    use spllt_data_mod
    use spllt_solve_mod
    use spllt_data_ciface
    use task_manager_omp_mod
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),            value         :: cfkeep
    type(spllt_options_t),  intent(in)    :: coptions
    type(C_PTR),            value         :: corder
    integer(C_INT),         value         :: nrhs
    type(C_PTR),            value         :: cx
    type(spllt_inform_t),   intent(out)   :: cinfo
    integer(C_INT),         value         :: job
    type(C_PTR),            value         :: cwork
    integer(C_LONG),        value         :: cworksize
    type(C_PTR),            value         :: ctask_manager

    type(spllt_fkeep),        pointer :: ffkeep
    type(spllt_options)               :: foptions
    integer,                  pointer :: forder(:)
    real(wp),                 pointer :: fx(:)   
    type(spllt_inform)                :: finfo
    real(wp),                 pointer :: fwork(:)
    type(task_manager_omp_t), pointer :: ftask_manager

    if(.not. C_ASSOCIATED(cfkeep)) then
      write (stderr,*) "Error, fkeep provided by the user is empty"
    end if
    call C_F_POINTER(cfkeep, ffkeep)

    call copy_options_in(coptions, foptions)

    if(C_ASSOCIATED(corder)) then
      call C_F_POINTER(corder, forder, shape=(/ ffkeep%n /))
    else
      nullify(forder)
    end if

    if(.not. C_ASSOCIATED(cx)) then
      write (stderr,*) "Error, x/rhs provided by the user is empty"
    end if

   !if(.not. C_ASSOCIATED(cwork)) then
   !  write (stderr,*) "Error, the workspace provided by the user is empty"
   !end if

    if(.not. C_ASSOCIATED(ctask_manager)) then
      write (stderr,*) "Error, the task_manager provided by the user is empty"
    end if

    call C_F_POINTER(cx, fx, shape=(/ ffkeep%n /))
    call C_F_POINTER(cwork, fwork, shape=(/ cworksize /))
    call C_F_POINTER(ctask_manager, ftask_manager)

    call ftask_manager%refresh_master()

    call spllt_solve_mult_double_worker(ffkeep, foptions, nrhs, fx, &
      job, ftask_manager, finfo)

    call copy_inform_out(finfo, cinfo)

  end subroutine spllt_c_solve_worker



  subroutine spllt_c_check_backward_error(n, cptr, crow, cval, nrhs, cx, crhs)&
      bind(C, name="spllt_chkerr")
    use utils_mod
    use spllt_data_mod
    use spllt_data_ciface
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    integer(C_INT), value  :: n
    type(C_PTR),    value  :: cptr
    type(C_PTR),    value  :: crow
    type(C_PTR),    value  :: cval
    integer(C_INT), value  :: nrhs
    type(C_PTR),    value  :: cx
    type(C_PTR),    value  :: crhs

    integer(C_INT),     pointer :: fptr(:)
    integer(C_INT),     pointer :: frow(:)
    real(wp),           pointer :: fval(:)
    real(wp),           pointer :: fx(:,:)
    real(wp),           pointer :: frhs(:,:)

    if(.not. C_ASSOCIATED(cptr)) then
      write (stderr,*) "Error, ptr provided by the user is empty"
    end if

    if(.not. C_ASSOCIATED(crow)) then
      write (stderr,*) "Error, row provided by the user is empty"
    end if

    if(.not. C_ASSOCIATED(cval)) then
      write (stderr,*) "Error, val provided by the user is empty"
    end if

    if(.not. C_ASSOCIATED(cx)) then
      write (stderr,*) "Error, x provided by the user is empty"
    end if

    if(.not. C_ASSOCIATED(crhs)) then
      write (stderr,*) "Error, rhs provided by the user is empty"
    end if

    call C_F_POINTER(cptr, fptr, shape=(/ n + 1 /))
    call C_F_POINTER(crow, frow, shape=(/ fptr(n + 1) - 1 /))
    call C_F_POINTER(cval, fval, shape=(/ fptr(n + 1) - 1 /))
    call C_F_POINTER(cx,   fx,   shape=(/ n, nrhs /))
    call C_F_POINTER(crhs, frhs, shape=(/ n, nrhs /))

    call check_backward_error(n, fptr, frow, fval, nrhs, fx, frhs)

  end subroutine spllt_c_check_backward_error


  subroutine spllt_c_deallocate_akeep(cakeep, cstat)  &
      bind(C, name="spllt_deallocate_akeep")
    use spllt_data_mod
    use spllt_data_ciface
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none
    type(C_PTR),    intent(inout) :: cakeep
    integer(C_INT), intent(out)   :: cstat

    integer :: stat
    type(spllt_akeep),  pointer :: fakeep

    if(.not. C_ASSOCIATED(cakeep)) then
      write (stderr,*) "Error, akeep provided by the user is empty"
      cstat = 1
      return
    end if
    call C_F_POINTER(cakeep, fakeep)

    call spllt_deallocate_akeep(fakeep, stat)
    
    if(stat .ne. 0) then
      write (stderr,*) "Error in deallocation of akeep : ", stat
    else
      deallocate(fakeep)
      cakeep  = C_NULL_PTR
    endif
    cstat   = stat

  end subroutine spllt_c_deallocate_akeep



  subroutine spllt_c_deallocate_fkeep(cfkeep, cstat)  &
      bind(C, name="spllt_deallocate_fkeep")
    use spllt_data_mod
    use spllt_data_ciface
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none
    type(C_PTR),    intent(inout) :: cfkeep
    integer(C_INT), intent(out)   :: cstat

    integer                     :: fstat
    type(spllt_fkeep),  pointer :: ffkeep

    if(.not. C_ASSOCIATED(cfkeep)) then
      write (stderr,*) "Error, fkeep provided by the user is empty"
      cstat = 1
      return
    end if
    call C_F_POINTER(cfkeep, ffkeep)

    call spllt_deallocate_fkeep(ffkeep, fstat)
    
    if(fstat .ne. 0) then
      write (stderr,*) "Error in deallocation of fkeep : ", fstat
    else
      deallocate(ffkeep)
      cfkeep  = C_NULL_PTR
    endif
    cstat   = fstat

  end subroutine spllt_c_deallocate_fkeep



  subroutine spllt_c_task_manager_deallocate(ctask_manager, cstat) &
      bind(C, name="spllt_task_manager_deallocate")
    use task_manager_omp_mod
    use spllt_data_ciface
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),    intent(inout) :: ctask_manager
    integer(C_INT), intent(out)   :: cstat

    type(task_manager_omp_t),  pointer  :: ftask_manager
    integer                             :: fstat

    fstat = 0
    if(.not. C_ASSOCIATED(ctask_manager)) then
      write (stderr,*) "Error, task_manager provided by the user is empty"
      cstat = 1
      return
    end if

    call C_F_POINTER(ctask_manager, ftask_manager)
    
    call ftask_manager%deallocate()

    if(fstat .ne. 0) then
      write (stderr,*) "Error in deallocation of fkeep : ", fstat
    else
      deallocate(ftask_manager)
      ctask_manager  = C_NULL_PTR
    endif
    cstat   = fstat
  end subroutine spllt_c_task_manager_deallocate



  subroutine spllt_c_all(cakeep, cfkeep, coptions, n, nnz, nrhs, nb, cptr, &
    crow, cval, cx, crhs, cinfo) bind(C, name="spllt_all")
    use spllt_analyse_mod
    use spllt_solve_dep_mod
    use spllt_solve_mod
    use spllt_data_mod
    use spllt_data_ciface
    use task_manager_omp_mod
    use timer_mod
    use utils_mod
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),            intent(inout) :: cakeep
    type(C_PTR),            intent(inout) :: cfkeep
    type(spllt_options_t),  intent(in)    :: coptions
    integer(C_INT),         value         :: n
    integer(C_INT),         value         :: nnz
    integer(C_INT),         value         :: nrhs
    integer(C_INT),         value         :: nb
    type(C_PTR),            value         :: cptr
    type(C_PTR),            value         :: crow
    type(C_PTR),            value         :: cval
    type(C_PTR),            value         :: cx
    type(C_PTR),            value         :: crhs
    type(spllt_inform_t),   intent(out)   :: cinfo

    type(spllt_akeep),  pointer :: fakeep
    type(spllt_fkeep),  pointer :: ffkeep
    type(spllt_options)         :: foptions
    integer,       allocatable  :: forder(:)
    integer(C_INT),     pointer :: fptr(:)
    integer(C_INT),     pointer :: frow(:)
    real(wp),           pointer :: fval(:)
    type(spllt_inform)          :: finfo
    integer                     :: st
    integer(long)               :: fworksize
    real(wp),     allocatable   :: fy(:)   
    real(wp),     allocatable   :: fw(:)   
    real(wp),           pointer :: fx(:,:)   
    real(wp),           pointer :: frhs(:,:)
    type(task_manager_omp_t)    :: task_manager

    call copy_options_in(coptions, foptions)

    if(C_ASSOCIATED(cakeep)) then
      call C_F_POINTER(cakeep, fakeep)
    else
      allocate(fakeep)
      cakeep = C_LOC(fakeep)
    end if

    if(C_ASSOCIATED(cfkeep)) then
      call C_F_POINTER(cfkeep, ffkeep)
    else
      allocate(ffkeep)
      cfkeep = C_LOC(ffkeep)
    end if

    if(C_ASSOCIATED(cptr)) then
      call C_F_POINTER(cptr, fptr, shape=(/ n + 1 /))
    else
      nullify(fptr)
    end if

    if(C_ASSOCIATED(crow)) then
      call C_F_POINTER(crow, frow, shape=(/ fptr(n + 1) - 1 /))
    else
      nullify(frow)
    end if

    if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape=(/ nnz /))
    else
      nullify(fval)
    end if

    if(.not. C_ASSOCIATED(crhs)) then
      write (stderr,*) "Error, rhs provided by the user is empty"
    end if

    call C_F_POINTER(crhs, frhs, shape=(/ n, nrhs /))

    if(.not. C_ASSOCIATED(cx)) then
      write (stderr,*) "Error, x/rhs provided by the user is empty"
    end if

   !call C_F_POINTER(cx, fx, shape=(/ ffkeep%n, nrhs/))
    call C_F_POINTER(cx, fx, shape=(/ ffkeep%n, nrhs /))

    print *, "Call of spllt_analyse"
    call spllt_analyse(fakeep, ffkeep, foptions, n, fptr, frow, finfo, forder)

    print *, "Call of spllt_factor"
    call spllt_factor(fakeep, ffkeep, foptions, fval, finfo)
    call spllt_wait()

    call spllt_init_timer(st)
    if(st .ne. 0) then
      write (stderr,*) "Warning, timer environment is not set"
    end if

    print *, "Create subtrees"
    call spllt_create_subtree(fakeep, ffkeep)

    print *, "Compute Sblocks with nb = ", nb, "nrhs = ", nrhs
    call get_solve_blocks(ffkeep, nb, nrhs, fworksize, ffkeep%sbc)

    print *, "Compute dep"
    call spllt_compute_solve_dep(ffkeep, stat = st)

    allocate(fy(n * nrhs), stat=st)
    allocate(fw(fworksize), stat=st)
    call sblock_assoc_mem(ffkeep, nb, nrhs, fy, fw, ffkeep%sbc)

    call task_manager%init()

    call spllt_solve_mult_double_worker(ffkeep, foptions, nrhs, fx, &
      0, task_manager, finfo)
    call spllt_wait()

    call check_backward_error(n, fptr, frow, fval, nrhs, fx, frhs)

    call copy_inform_out(finfo, cinfo)
  end subroutine spllt_c_all
