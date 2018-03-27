module utils_mod


  interface print_array
    module procedure print_darray
    module procedure print_iarray
  end interface print_array

  interface timer_log_dump
    module procedure timer_log_dump_one
    module procedure timer_log_dump_mult
  end interface timer_log_dump

  interface flop_log_dump
    module procedure flop_log_dump_one
    module procedure flop_log_dump_mult
  end interface flop_log_dump
contains

  subroutine print_darray(array_name, n, val)
    use spllt_data_mod
    character(len=*),       intent(in)    :: array_name
    integer,                intent(in)    :: n
    real(wp), dimension(n), intent(in)    :: val

    integer :: i

    print '(a)', array_name
    do i = 1, n
      print *, val(i)
    end do
  end subroutine print_darray

  subroutine print_iarray(array_name, n, val, display)
    character(len=*),       intent(in)    :: array_name
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

  subroutine print_blk_index(array_name, n, val, display)
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
      write(*, fmt="(i9)", advance="no") val(1)
      do i = 1, n - 1
        if(val(i+1) - val(i) .gt. 1) then
          print '(a,i9)', " : ", val(i)
          if(i + 1 .lt. n) then
            write(*, fmt="(i9)", advance="no") val(i+1)
          end if
        end if
      end do
      print '(a,i9)', ":", val(n)
    else
      write(*, fmt="(i9)", advance="no") val(1)
      do i = 1, n - 1
        if(val(i+1) - val(i) .gt. 1) then
          write(*, fmt='(a,i9,a)', advance="no") " : ", val(i), ","
          if(i + 1 .lt. n) then
            write(*, fmt="(i9)", advance="no") val(i+1)
          end if
        end if
      end do
      print '(a,i9)', ":", val(n)
    end if
  end subroutine print_blk_index

  subroutine print_node(fkeep, node_num)
    use spllt_data_mod
    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: node_num

    integer :: ncol, last_blk, first_blk, i, j, nrow

    first_blk = fkeep%nodes(node_num)%blk_sa
    last_blk  = fkeep%nodes(node_num)%blk_en
    ncol      = fkeep%bc(last_blk)%bcol - fkeep%bc(first_blk)%bcol + 1
    nrow      = fkeep%bc(first_blk)%last_blk - first_blk + 1
    do i = 1, nrow
      do j = 1, min(i, ncol)
        write(*, fmt="(i9)", advance="no") &
          first_blk + int((j - 1) * ( nrow + 1 - 0.5 * j )) + (i - j)
      end do
      write (*,*) ""
    end do
  end subroutine print_node

  !Compute res = b - Ax 
  subroutine compute_residual(n, ptr, row, val, nrhs, x, b, res)
    use spllt_data_mod
    integer, intent(in)                         :: n
    integer, dimension(n+1), intent(in)         :: ptr
    integer, dimension(ptr(n+1)-1), intent(in)  :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in)                         :: nrhs
    real(wp), dimension(n,nrhs), intent(in)     :: x
    real(wp), dimension(n,nrhs), intent(in)     :: b
    real(wp), dimension(n,nrhs), intent(inout)  :: res

    ! Find the residual
    res = 0

    call compute_Ax(n, ptr, row, val, nrhs, x, res)

    res = b - res
  end subroutine compute_residual
  
  !Compute Ax
  subroutine compute_Ax(n, ptr, row, val, nrhs, x, res)
    use spllt_data_mod
    integer, intent(in)                         :: n
    integer, dimension(n+1), intent(in)         :: ptr
    integer, dimension(ptr(n+1)-1), intent(in)  :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in)                         :: nrhs
    real(wp), dimension(n,nrhs), intent(in)     :: x
    real(wp), dimension(n,nrhs), intent(inout)  :: res

    integer :: i, j, k, r
    res = 0
    do i = 1, n
      do j = ptr(i), ptr(i+1)-1
        r = row(j)
        do k = 1, nrhs
          res(r, k) = res(r, k) + val(j) * x(i, k)
          if(r .eq. i) cycle
          res(i, k) = res(i, k) + val(j) * x(r, k)
        end do
      end do
    end do
  end subroutine compute_Ax

  subroutine matrix_norm_1(n, ptr, row, val, norm)
    use spllt_data_mod
    integer, intent(in)                         :: n
    integer, dimension(n+1), intent(in)         :: ptr
    integer, dimension(ptr(n+1)-1), intent(in)  :: row    ! Unused
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    real(wp), intent(out)                       :: norm

    integer   :: i, j
    real(wp)  :: sum_col

    norm = 0
    do i = 1, n
      sum_col = 0
      do j =  ptr(i), ptr(i + 1) - 1
        sum_col = sum_col + abs(val(j))
      end do
      norm = merge(norm, sum_col, norm > sum_col)
    end do
  end subroutine matrix_norm_1

  subroutine vector_norm_2(n, val, norm)
    use spllt_data_mod
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: val(:,:)
    real(wp), intent(out) :: norm(:)
    
    integer :: i
    norm = 0
  
    do i = 1, n
      norm(:) = norm(:) + val(i,:) * val(i, :)
    end do
    norm = sqrt(norm)

  end subroutine vector_norm_2

  subroutine matrix_norm_max(n, ptr, row, val, norm)
    use spllt_data_mod
    integer, intent(in)                         :: n
    integer, dimension(n+1), intent(in)         :: ptr
    integer, dimension(ptr(n+1)-1), intent(in)  :: row    ! Unused
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    real(wp), intent(out)                       :: norm

    integer :: i

    norm = 0
    do i = 1, ptr(n+1)-1
      norm = merge(norm, abs(val(i)), norm > abs(val(i)))
    end do
  end subroutine matrix_norm_max

  subroutine print_task_stat(task, msg)
    use spllt_data_mod
    type(spllt_omp_task_stat), intent(in)   :: task
    character(len=*), optional, intent(in)  :: msg

    if(present(msg)) then
      print *, msg
    end if
    print '(a, i9)', "max #dep of a blk   : ", task%max_dep
    print '(a, i1, a, i9)', "#blk with #dep>", k_dep,"    : ", task%nblk_kdep

  end subroutine print_task_stat

  subroutine print_omp_task_stat(msg, task_id, task)
    use spllt_data_mod
    character(len=*), intent(in)          :: msg
    integer, intent(in)                   :: task_id
    type(spllt_omp_task_stat), intent(in) :: task

    print *, msg, " : ", task_id
    print '(a, i9)',    "#task insert        : ", task%ntask_insert
    print '(a, i9)',    "#fake task insert   : ", task%nfake_task_insert
    print '(a, i9)',    "#task run           : ", task%ntask_run
    print '(a, i9)',    "#array allocate     : ", task%narray_allocated
    print '(a, es10.2)', "#flop               : ", task%nflop

  end subroutine print_omp_task_stat

  subroutine print_scheduler(sched)
    use spllt_data_mod
    type(spllt_omp_scheduler), intent(in) :: sched

    integer :: i

    print *, "Scheduler state"
    print '(a, i3)', "workerID      :", sched%workerID
    print '(a, i3)', "masterWorker  :", sched%masterWorker
    print '(a, i3)', "nworker       :", sched%nworker
    print '(a, i3)', "nthread_max   :", sched%nthread_max
    do i = 1, sched%nworker
      call print_omp_task_stat("Init omp info task", i, sched%task_info(i))
      call print_task_stat(sched%task_info(i))
    end do

  end subroutine print_scheduler
  
  subroutine spllt_omp_init_task_info(task_stat)
    use spllt_data_mod
    type(spllt_omp_task_stat), intent(out) :: task_stat

    task_stat%nflop               = 0.0
    task_stat%ntask_run           = 0
    task_stat%ntask_insert        = 0
    task_stat%nfake_task_insert   = 0
    task_stat%max_dep             = 0
    task_stat%narray_allocated    = 0
    task_stat%nblk_kdep           = 0
    
  end subroutine spllt_omp_init_task_info

  subroutine spllt_omp_reset_scheduler(scheduler)
    use spllt_data_mod
 !$ use omp_lib, only : omp_get_num_threads, omp_get_thread_num
    type(spllt_omp_scheduler), intent(inout) :: scheduler

    scheduler%workerID                = 0
 !$ scheduler%workerID                = omp_get_thread_num()
    scheduler%nworker                 = 1
 !$ scheduler%nworker                 = omp_get_num_threads()
    scheduler%masterWorker            = scheduler%workerID
    scheduler%nthread_max             = scheduler%nworker

  end subroutine spllt_omp_reset_scheduler



  subroutine spllt_omp_init_scheduler(scheduler, trace_names, stat)
    use spllt_data_mod
    use trace_mod, only : trace_create_event
 !$ use omp_lib, only : omp_get_num_threads, omp_get_thread_num
    type(spllt_omp_scheduler), target,  intent(inout) :: scheduler
    character(len=*), optional,         intent(in)    :: trace_names(:)
    integer, optional,                  intent(out)   :: stat

    integer                             :: st1, st2, i, ntrace_id
    type(spllt_omp_task_stat), pointer  :: p_task_info

    st1 = 0
    st2 = 0

    call spllt_omp_reset_scheduler(scheduler)

    allocate(scheduler%task_info(0:scheduler%nthread_max-1), stat=st1)
    if(st1 .eq. 0) then
      do i = lbound(scheduler%task_info, 1), ubound(scheduler%task_info, 1)
        p_task_info => scheduler%task_info(i)
        call spllt_omp_init_task_info(p_task_info)
      end do
    end if

    if(present(trace_names)) then
      ntrace_id = size(trace_names)
      allocate(scheduler%trace_ids(ntrace_id), stat=st2)
      if(st2 .eq. 0) then
        do i = 1, ntrace_id
          call trace_create_event(trace_names(i), scheduler%trace_ids(i))
!         print *, "Create id ", scheduler%trace_ids(i), " for step ", &
!           trace_names(i)
        end do
      end if
    else
      scheduler%trace_ids => null()
    end if

    if(present(stat)) then
      stat = st1 + st2
    end if

!   call print_scheduler(scheduler)

  end subroutine spllt_omp_init_scheduler



  subroutine spllt_scheduler_alloc(scheduler, stat)
    use spllt_data_mod
    type(spllt_omp_scheduler), intent(inout)  :: scheduler
    integer,                   intent(in)     :: stat

    if(stat .ne. 0) then
      write(0,'(a)') "Error of alloc"
    end if
    scheduler%task_info(scheduler%workerID)%narray_allocated =  &
      scheduler%task_info(scheduler%workerID)%narray_allocated  &
      + merge(1, 0, stat .eq. 0)

  end subroutine spllt_scheduler_alloc



  subroutine spllt_update_omp_task_info(task_info, ntask, nftask)
    use spllt_data_mod
    type(spllt_omp_task_stat), intent(inout)  :: task_info
    integer, intent(in)                       :: ntask  ! #task insert
    integer, intent(in)                       :: nftask ! #fake task

    task_info%nfake_task_insert = task_info%nfake_task_insert + nftask
    task_info%ntask_insert      = task_info%ntask_insert + ntask

  end subroutine spllt_update_omp_task_info



  subroutine spllt_update_task_info(task_info, ndep)
    use spllt_data_mod
    type(spllt_omp_task_stat), intent(inout)  :: task_info
    integer, intent(in)                       :: ndep   ! #dep of the block

    task_info%nblk_kdep         = task_info%nblk_kdep + &
      merge(1, 0, ndep .gt. k_dep)
    task_info%max_dep           = merge(ndep, task_info%max_dep, &
      ndep .gt. task_info%max_dep)

  end subroutine spllt_update_task_info



  subroutine timer_log_dump_mult(header, timer, ofile)

    character, intent(in)         :: header*(*)
    double precision, intent(in)  :: timer(:,:)
    character, intent(in)         :: ofile*(*)

    integer :: n, i

    n = size(timer,1)

    open(4, file="timer_"//ofile//".data", action='write')
    
    write(4,'(a, a)') "#", trim(header)
    do i = 1, n
      write(4, *) timer(i,:)
    end do

    close(4)
  end subroutine timer_log_dump_mult

  subroutine timer_log_dump_one(header, timer, ofile)

    character, intent(in)         :: header*(*)
    double precision, intent(in)  :: timer(:)
    character, intent(in)         :: ofile*(*)

    integer :: n, i

    n = size(timer,1)

    open(4, file="timer_"//ofile//".data", action='write')
    
    write(4,'(a, a)') "#", trim(header)
    do i = 1, n
      write(4, *) timer(i)
    end do

    close(4)
  end subroutine timer_log_dump_one



  subroutine flop_log_dump_one(header, flop, ofile)

    character, intent(in)         :: header*(*)
    double precision, intent(in)  :: flop(:)
    character, intent(in)         :: ofile*(*)

    integer :: n, i

    n = size(flop,1)

    open(4, file="flop_"//ofile//".data", action='write')
    
    write(4,'(a, a)') "#", trim(header)
    do i = 1, n
      write(4, *) flop(i)
    end do

    close(4)
  end subroutine flop_log_dump_one



  subroutine flop_log_dump_mult(header, flop, ofile)

    character, intent(in)         :: header*(*)
    double precision, intent(in)  :: flop(:,:)
    character, intent(in)         :: ofile*(*)

    integer :: n, i

    n = size(flop,1)

    open(4, file="flop_"//ofile//".data", action='write')
    
    write(4,'(a, a)') "#", trim(header)
    do i = 1, n
      write(4, *) flop(i,:)
    end do

    close(4)
  end subroutine flop_log_dump_mult



  subroutine compute_range(vmin, vmax, linear_mode, val)
    integer, intent(in)               :: vmin
    integer, intent(in)               :: vmax
    logical, intent(in)               :: linear_mode
    integer, allocatable, intent(out) :: val(:)

    integer :: nval, offset

    if(linear_mode) then
      nval = int((vmax + 0.0) / vmin)
      allocate(val(nval))
      do i=1, nval
        val(i) = vmin * i
      end do
    else
      offset = int(log(real(vmin))/log(2.0))
      nval  = int(log(real(vmax))/log(2.0)) - offset + 1
      allocate(val(nval))
      do i=1, nval
        val(i) = 2 **(i - 1 + offset)
      end do
    end if

  end subroutine compute_range
! subroutine permute_darray(n, val, perm, val_perm, trans)
!   integer,                intent(in)      :: n
!   real(wp), dimension(n), intent(in)      :: val
!   integer,  dimension(n), intent(in)      :: perm
!   real(wp), dimension(n), intent(out)     :: val_perm
!   character(len=1), optional, intent(in)  :: trans

!   integer           :: i
!   character(len=1)  :: permute_type

!   if(.not. present(trans)) then
!     permute_type = 'N'
!   else
!     permute_type = trans
!   end if

!   if(permute_type == 'T') then
!     do i = 1, n
!       val_perm(perm(i)) = val(i)
!     end do
!   else
!     do i = 1, n
!       val_perm(i) = val(perm(i))
!     end do
!   end if
! end subroutine permute_darray

! subroutine permute_array(n, val, perm, val_perm)
!   integer,                intent(in)    :: n
!   real(wp), dimension(n), intent(in)    :: val
!   integer,  dimension(n), intent(in)    :: perm
!   real(wp), dimension(n), intent(out)   :: val_perm

!   integer :: i

!   do i = 1, n
!     val_perm(i) = val(perm(i))
!   end do
! end subroutine permute_array
end module utils_mod
