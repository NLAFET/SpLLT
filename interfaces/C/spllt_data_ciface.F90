module spllt_data_ciface                                                       
  use iso_c_binding
  implicit none

! type, bind(C) :: spllt_lfactor_t
!    real(C_DOUBLE), allocatable :: lcol(:)
! end type spllt_lfactor_t

! type, bind(C) :: spllt_lmap_type_t
!    integer(C_INT)               :: len_map
!    integer(C_INT), allocatable  :: map(:,:)
! end type spllt_lmap_type_t

! type, bind(C) :: spllt_block_t
!    real(C_PTR)                  :: c(:)
!    integer(C_INT)               :: mem_node
!    integer(C_INT)               :: bcol
!    integer(C_INT)               :: blkm
!    integer(C_INT)               :: blkn
!    integer(C_INT)               :: dblk
!    integer(C_INT)               :: dep_initial
!    integer(C_INT)               :: id
!    integer(C_INT)               :: last_blk
!    integer(C_INT)               :: node
!    integer(C_INT)               :: sa
!    integer(C_INT)               :: dep
!    integer(C_INT), allocatable  :: fwd_dep(:)
!    integer(C_INT), allocatable  :: fwd_update_dep(:)
!    integer(C_INT)               :: fwd_solve_dep
!    integer(C_INT)               :: bwd_update_dep
!    integer(C_INT), allocatable  :: bwd_solve_dep(:)
!    integer(C_INT), allocatable  :: bwd_dep(:)
! end type spllt_block_t

! type, bind(C) :: spllt_node_t
!    integer(C_INT)               :: num
!    type(spllt_block_t)          :: buffer
!    integer(C_INT)               :: blk_sa
!    integer(C_INT)               :: blk_en
!    integer(C_INT)               :: nb
!    integer(C_INT)               :: sa
!    integer(C_INT)               :: en
!    integer(C_INT), allocatable  :: index(:)
!    integer(C_INT)               :: nchild
!    integer(C_INT), allocatable  :: child(:)
!    integer(C_INT)               :: parent
!    integer(C_INT)               :: least_desc
!    integer(C_INT), allocatable  :: extra_row(:)
! end type spllt_node_t

! type, bind(C) :: spllt_workspace_i_t
!    integer(C_PTR) :: c(:)
! end type spllt_workspace_i_t

  type, bind(C) :: spllt_options_t
     integer(C_INT)     :: print_level
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
  end type spllt_options_t

! type spllt_tree_t
!   integer(C_INT)              :: num
!   integer(C_INT)              :: nnode
!   integer(C_INT)              :: node_sa
!   integer(C_INT)              :: node_en
!   integer(C_INT)              :: nchild
!   integer(C_INT), allocatable :: child(:)
!   integer(C_INT)              :: parent
! end type spllt_tree_t

  type, bind(C) :: spllt_inform_t
!   !type(ssids_inform) :: ssids_inform
     integer(C_INT) :: flag
     integer(C_INT) :: maxdepth
     integer(C_INT) :: num_factor
     integer(C_INT) :: num_flops
     integer(C_INT) :: num_nodes
     integer(C_INT) :: stat
  end type spllt_inform_t

  type, bind(C) :: spllt_akeep_t
    integer(C_INT)  :: nnodes
    integer(C_INT)  :: n
    integer(C_INT)  :: num_factor
    integer(C_INT)  :: num_flops
    type(C_PTR)     :: weight
    type(C_PTR)     :: small
  end type spllt_akeep_t

  type, bind(C) :: spllt_fkeep_t
     type(C_PTR)          :: bc
     type(C_PTR)          :: workspace
     type(C_PTR)          :: nodes
     type(C_PTR)          :: row_list
     type(C_PTR)          :: col_list
     type(C_PTR)          :: map
     type(C_PTR)          :: flag_array
     integer(C_INT)       :: final_blk
     type(spllt_inform_t) :: info
     integer(C_INT)       :: maxmn
     integer(C_INT)       :: n
     integer(C_INT)       :: nbcol
     type(C_PTR)          :: lfact
     type(C_PTR)          :: lmap
     type(C_PTR)          :: trees
     type(C_PTR)          :: small
     type(C_PTR)          :: assoc_tree
  end type spllt_fkeep_t

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

    !$omp single
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
    !$omp end single
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

    !$omp single
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
    call spllt_wait()

    call copy_inform_out(finfo, cinfo)
    !$omp end single

  end subroutine spllt_c_factor



  subroutine spllt_c_prepare_solve(cakeep, cfkeep, cinfo) &
      bind(C, name="spllt_prepare_solve")
    use spllt_data_mod
    use spllt_solve_dep_mod
    use spllt_data_ciface
    use timer_mod
    use ISO_Fortran_env, only: stderr => ERROR_UNIT
    implicit none

    type(C_PTR),            value         :: cakeep
    type(C_PTR),            value         :: cfkeep
    type(spllt_inform_t),   intent(out)   :: cinfo

    type(spllt_akeep),  pointer :: fakeep
    type(spllt_fkeep),  pointer :: ffkeep
    type(spllt_inform)          :: finfo
    integer                     :: st

    !$omp single
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

    call spllt_create_subtree(fakeep, ffkeep)
    call spllt_compute_solve_dep(ffkeep)

    call copy_inform_out(finfo, cinfo)
    !$omp end single

  end subroutine spllt_c_prepare_solve



  subroutine spllt_c_solve(cfkeep, coptions, corder, nrhs, cx, cinfo, job) &
      bind(C, name="spllt_solve")
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

    !$omp single
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

    call C_F_POINTER(cx, fx, shape=(/ ffkeep%n /))

    call task_manager%init()

    call spllt_solve_mult_double(ffkeep, foptions, forder, nrhs, fx, finfo, &
      job, task_manager)

    call task_manager%deallocate()

    call copy_inform_out(finfo, cinfo)
    !$omp end single

  end subroutine spllt_c_solve



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
    real(wp),           pointer :: fx(:)
    real(wp),           pointer :: frhs(:)

    !$omp single
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
    call C_F_POINTER(cx,   fx,   shape=(/ n /))
    call C_F_POINTER(crhs, frhs, shape=(/ n /))

    call check_backward_error(n, fptr, frow, fval, nrhs, fx, frhs)
    !$omp end single

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

    integer :: stat
    type(spllt_fkeep),  pointer :: ffkeep

    if(.not. C_ASSOCIATED(cfkeep)) then
      write (stderr,*) "Error, fkeep provided by the user is empty"
      cstat = 1
      return
    end if
    call C_F_POINTER(cfkeep, ffkeep)

    call spllt_deallocate_fkeep(ffkeep, stat)
    
    if(stat .ne. 0) then
      write (stderr,*) "Error in deallocation of fkeep : ", stat
    else
      deallocate(ffkeep)
      cfkeep  = C_NULL_PTR
    endif
    cstat   = stat

  end subroutine spllt_c_deallocate_fkeep

! subroutine fkeep_default_control(cfkeep) bind(C, name="spllt_fkeep_init")
!   use spllt_data_mod
!   use spllt_data_ciface
!   implicit none
!   type(spllt_fkeep_t), intent(out) :: cfkeep

!   type(spllt_fkeep) :: ffkeep

!    cfkeep%bc          = C_NULL_PTR
!    cfkeep%workspace   = C_NULL_PTR 
!    cfkeep%nodes       = C_NULL_PTR 
!    cfkeep%row_list    = C_NULL_PTR 
!    cfkeep%col_list    = C_NULL_PTR 
!    cfkeep%map         = C_NULL_PTR 
!    cfkeep%flag_array  = C_NULL_PTR 
!    cfkeep%final_blk   = ffkeep%final_blk
!    cfkeep%maxmn       = ffkeep%maxmn
!    cfkeep%n           = ffkeep%n
!    cfkeep%nbcol       = ffkeep%nbcol
!    cfkeep%lfact       = C_NULL_PTR 
!    cfkeep%lmap        = C_NULL_PTR 
!    cfkeep%trees       = C_NULL_PTR 
!    cfkeep%small       = C_NULL_PTR 
!    cfkeep%assoc_tree  = C_NULL_PTR 
!    call copy_inform_out(ffkeep%info, cfkeep%info)

! end subroutine fkeep_default_control



! subroutine akeep_default_control(cakeep) bind(C, name="spllt_akeep_init")
!   use spllt_data_mod
!   use spllt_data_ciface
!   implicit none
!   type(spllt_akeep_t), intent(out) :: cakeep

!   type(spllt_akeep) :: fakeep

!   cakeep%nnodes     = fakeep%nnodes
!   cakeep%n          = fakeep%n
!   cakeep%num_factor = fakeep%num_factor
!   cakeep%num_flops  = fakeep%num_flops
!   cakeep%weight     = C_NULL_PTR
!   cakeep%small      = C_NULL_PTR

! end subroutine akeep_default_control

