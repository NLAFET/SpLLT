module timer_mod
  use get_wtime_mod
  use iso_c_binding
  implicit none

  integer, parameter  :: max_steps = 15

  type spllt_steps
   !character(len=:), pointer :: p
    character(len=200)        :: name(0 : max_steps)
  end type spllt_steps

  ! Contains timers of each step of a subroutine
  type spllt_timer
   !character(len=200)        :: step_name(0 : max_steps)
    type(spllt_steps), pointer  :: steps
    double precision,  pointer  :: start(:,:)
    double precision,  pointer  :: stop(:,:)
    double precision,  pointer  :: swap(:)
    integer,           pointer  :: ncall(:,:)
    integer                     :: state
    double precision,  pointer  :: time(:,:)
    double precision,  pointer  :: time_min(:,:)
    double precision,  pointer  :: time_max(:,:)
  end type spllt_timer

  ! Record the timer of all functions used
  type spllt_timers
    integer                         :: max_ntimers  = 10
    integer                         :: ntimers      = 0
    type(spllt_timer), allocatable  :: timers(:)
  end type spllt_timers

  type(spllt_timers), target, save :: all_timers

contains

  subroutine spllt_init_timer(stat)

    integer, intent(out) :: stat

    allocate(all_timers%timers(all_timers%max_ntimers), stat = stat)
    if(stat .ne. 0) then
      print *, "Error, can not allocate the timers"
    end if

  end subroutine spllt_init_timer



  subroutine save_timer(local_timer)
    type(spllt_timer), intent(in) :: local_timer

!   type(spllt_timer), allocatable :: buf(:)

    ! Increase size if necessary
    if(all_timers%ntimers .eq. all_timers%max_ntimers) then
      print *, "Error, can not get a new timer, no more space"
      stop
!     allocate(buf(2 * all_timers%max_ntimers))
!     buf(1 : all_timers%max_ntimers) = all_timers%timers
!     deallocate(all_timers%timers)
!     call move_alloc(buf, all_timers%timers)
!     all_timers%max_ntimers = all_timers%max_ntimers * 2
    end if

    all_timers%ntimers = all_timers%ntimers + 1
    all_timers%timers(all_timers%ntimers) = local_timer

  end subroutine save_timer



  subroutine spllt_open_timer(nthread, thn, fun_name, timer)
    integer,                    intent(in)    :: nthread
    integer,                    intent(in)    :: thn
    character(len=*),           intent(in)    :: fun_name
    type(spllt_timer),          intent(inout) :: timer
    
    integer :: step_id

    step_id = 0

    if(timer%state .eq. 0) then
      allocate(timer%start    (0 : max_steps, 0 : nthread - 1))
      allocate(timer%time     (0 : max_steps, 0 : nthread - 1))
      allocate(timer%time_min (0 : max_steps, 0 : nthread - 1))
      allocate(timer%time_max (0 : max_steps, 0 : nthread - 1))
      allocate(timer%ncall    (0 : max_steps, 0 : nthread - 1))
      allocate(timer%swap     (0 : nthread - 1))
      allocate(timer%steps)
      timer%state = 1
      timer%ncall(:, :) = 0
      timer%steps%name(step_id) = fun_name
      call save_timer(timer)
    end if

    timer%start(step_id, thn) = omp_get_wtime()

  end subroutine spllt_open_timer



  subroutine spllt_tic(step_name, step_id, thn, timer)
    character(len=*),           intent(in)    :: step_name
    integer,                    intent(in)    :: step_id
    integer,                    intent(in)    :: thn
    type(spllt_timer),          intent(inout) :: timer

    timer%start(step_id, thn) = omp_get_wtime()
    timer%steps%name(step_id)  = step_name

  end subroutine spllt_tic



  subroutine spllt_tac(step_id, thn, timer)
    integer,                    intent(in)    :: step_id
    integer,                    intent(in)    :: thn
    type(spllt_timer),          intent(inout) :: timer


    timer%swap(thn) = omp_get_wtime()
    timer%start(step_id, thn)  = timer%swap(thn) - timer%start(step_id, thn)

    ! Compute elapsed time
    ! Compute min/max
    if(timer%ncall(step_id, thn) .eq. 0) then
      timer%time(step_id, thn)     = timer%start(step_id, thn)
      timer%time_min(step_id, thn) = timer%start(step_id, thn)
      timer%time_max(step_id, thn) = timer%start(step_id, thn)
    else
      timer%time(step_id, thn)   = timer%time(step_id, thn) + &
        timer%start(step_id, thn)
      if(timer%time_min(step_id, thn) .gt. timer%start(step_id, thn)) then
        timer%time_min(step_id, thn) = timer%start(step_id, thn)
      end if
      if(timer%time_max(step_id, thn) .lt. timer%start(step_id, thn)) then
        timer%time_max(step_id, thn) = timer%start(step_id, thn)
      end if
    end if

    ! Counter
    timer%ncall(step_id, thn) = timer%ncall(step_id, thn) + 1
    ! Prepare in case of another call to tac
    timer%start(step_id, thn) = timer%swap(thn)

  end subroutine spllt_tac



  subroutine spllt_close_timer(thn, timer)
    integer,                    intent(in)    :: thn
    type(spllt_timer),          intent(inout) :: timer
    
    integer :: step_id

    step_id = 0

    timer%swap(thn) = omp_get_wtime()
    timer%start(step_id, thn)  = timer%swap(thn) - timer%start(step_id, thn)

    ! Compute elapsed time
    ! Compute min/max
    if(timer%ncall(step_id, thn) .eq. 0) then
      timer%time(step_id, thn)     = timer%start(step_id, thn)
      timer%time_min(step_id, thn) = timer%start(step_id, thn)
      timer%time_max(step_id, thn) = timer%start(step_id, thn)
    else
      timer%time(step_id, thn)   = timer%time(step_id, thn) + &
        timer%start(step_id, thn)
      if(timer%time_min(step_id, thn) .gt. timer%start(step_id, thn)) then
        timer%time_min(step_id, thn) = timer%start(step_id, thn)
      end if
      if(timer%time_max(step_id, thn) .lt. timer%start(step_id, thn)) then
        timer%time_max(step_id, thn) = timer%start(step_id, thn)
      end if
    end if
    ! Counter
    timer%ncall(step_id, thn) = timer%ncall(step_id, thn) + 1
    
   !print *, "[th ", thn, "] Timer of ", trim(timer%step_name(0)), &
   !  " with ncall = ", timer%ncall(0,thn), " is ", timer%time(0, thn)

  end subroutine spllt_close_timer



  subroutine spllt_print_timer(thn, timer)
    integer,            intent(in) :: thn
    type(spllt_timer),  intent(in) :: timer

    integer :: i != 0

    i = 0
!   print *, "i = ", i

    print '(a, i3, a, a20, a, es10.2, a, i7, a, es10.2, a, es10.2, a)', &
      "[Th: ", thn, "] ", trim(timer%steps%name(i)), " : ",        &
      timer%time(i, thn), " s (ncall ", timer%ncall(i, thn), ") [",     &
      timer%time_min(i, thn), ' , ', timer%time_max(i, thn), ' ]' 

    do i = 1, ubound(timer%ncall, 1)
      if(timer%ncall(i, thn) .gt. 0) then
        if(timer%ncall(i, thn) .eq. 1) then
          print '(20x, a20, a, es10.2)', trim(timer%steps%name(i)), " : ", timer%time(i, thn)
        else
          print '(20x, a20, a, es10.2, a, i7, a, es10.2, a, es10.2, a)',  &
            trim(timer%steps%name(i)), " : ", timer%time(i, thn), " s (", &
            timer%ncall(i, thn), ") [ ",                                  &
            timer%time_min(i, thn), ' , ', timer%time_max(i, thn), ' ]'
        end if
      end if
    end do

  end subroutine spllt_print_timer

  subroutine spllt_print_timers(thn)
    integer, intent(in) :: thn

    integer :: i, th

    print *, "=========================================="
    print *, "                  TIMER"
    print *, "=========================================="
    do th = 0, thn - 1
      do i = 1, all_timers%ntimers
        if(all_timers%timers(i)%ncall(0,th) .gt. 0) then
          call spllt_print_timer(th, all_timers%timers(i))
        end if
      end do
    end do
  end subroutine spllt_print_timers

end module timer_mod
