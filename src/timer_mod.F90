module timer_mod
  use get_wtime_mod
  use iso_c_binding
  implicit none

  integer, parameter  :: max_steps      = 15
  integer, parameter  :: max_thread     = 64
  integer, parameter  :: max_ntimer_th  = 20

  type spllt_steps
    character(len=200)        :: name(0 : max_steps)
  end type spllt_steps

  ! Contains timers of each step of a subroutine
 !type spllt_timer
 !  type(spllt_steps), pointer  :: steps          => null()
 !  double precision,  pointer  :: start(:,:)     => null()
 !  double precision,  pointer  :: stop(:,:)      => null()
 !  double precision,  pointer  :: swap(:)        => null()
 !  integer,           pointer  :: ncall(:,:)     => null()
 !  integer                     :: state          = 0
 !  double precision,  pointer  :: time(:,:)      => null()
 !  double precision,  pointer  :: time_min(:,:)  => null()
 !  double precision,  pointer  :: time_max(:,:)  => null()
 !end type spllt_timer

  type spllt_p_timer_t
    type(spllt_timer_t), pointer :: p => null()
  end type spllt_p_timer_t

  type spllt_timer_t
    type(spllt_steps) :: steps
    double precision  :: start(   0 : max_steps, 0 : max_thread - 1)
    double precision  :: stop(    0 : max_steps, 0 : max_thread - 1)
    double precision  :: swap(    0 : max_thread - 1)
    integer           :: ncall(   0 : max_steps, 0 : max_thread - 1)  = 0
    double precision  :: time(    0 : max_steps, 0 : max_thread - 1)  = 0.0
    double precision  :: time_min(0 : max_steps, 0 : max_thread - 1)  = 1e30
    double precision  :: time_max(0 : max_steps, 0 : max_thread - 1)  = 0.0
    integer           :: status(  0 : max_thread - 1)
    double precision  :: flop(    0 : max_steps, 0 : max_thread - 1)  = 0.0
    double precision  :: flop_min(0 : max_steps, 0 : max_thread - 1)  = 1e30
    double precision  :: flop_max(0 : max_steps, 0 : max_thread - 1)  = 0.0
  end type spllt_timer_t

  ! Record the timer of all functions used
  type spllt_timers_t
    integer                             :: ntimer_th
    integer                             :: nthread
    integer              , allocatable  :: ntimer(:) 
    type(spllt_p_timer_t), allocatable  :: timers(:,:)
  end type spllt_timers_t

  type(spllt_timers_t), target, save :: all_timers

 !type spllt_flop_t
 !  double precision :: flop_rate(0 : max_steps, 0 : max_thread - 1)
 !end type spllt_flop_t
   interface spllt_tac
      module procedure spllt_tac
      module procedure spllt_tac_flop
   end interface
   interface spllt_close_timer
      module procedure spllt_close_timer
      module procedure spllt_close_timer_flop
   end interface
contains

  subroutine spllt_init_timer(stat, nthread, ntimer_th)

    integer,            intent(out) :: stat
    integer, optional,  intent(in)  :: nthread
    integer, optional,  intent(in)  :: ntimer_th

    if(present(ntimer_th)) then
      all_timers%ntimer_th = ntimer_th
    else
      all_timers%ntimer_th = max_ntimer_th
    end if

    if(present(nthread)) then
      all_timers%nthread = nthread
    else
      all_timers%nthread = max_thread
    end if

    allocate(all_timers%timers(all_timers%ntimer_th, &
      0 : all_timers%nthread - 1), &
      stat = stat)

    if(stat .ne. 0) then
      write (0,*) "Error, can not allocate the timers"
    end if
  
    allocate(all_timers%ntimer(0 : all_timers%nthread - 1), stat = stat)

    if(stat .ne. 0) then
      write (0,*) "Error, can not allocate the timers"
    end if

    all_timers%ntimer(:) = 0.0

  end subroutine spllt_init_timer



  subroutine save_timer(local_timer, thn)
    implicit none
    type(spllt_timer_t),  target, intent(in)  :: local_timer
    integer,                      intent(in)  :: thn

    ! Can not increase size if necessary
    if(all_timers%ntimer(thn) .eq. all_timers%ntimer_th) then
      write (0,*) "Error, can not add a new timer, no more space"
      stop
    end if

    all_timers%ntimer(thn) = all_timers%ntimer(thn) + 1
    all_timers%timers(all_timers%ntimer(thn), thn)%p => local_timer

  end subroutine save_timer



  subroutine spllt_open_timer(thn, fun_name, timer)
    integer,                    intent(in)    :: thn
    character(len=*),           intent(in)    :: fun_name
    type(spllt_timer_t),        intent(inout) :: timer
    
    integer :: step_id

    step_id = 0

    if(timer%status(thn) .eq. 0) then
      timer%status(thn) = 1
      timer%steps%name(step_id) = fun_name
      call save_timer(timer, thn)
    end if

    timer%start(step_id, thn) = omp_get_wtime()

  end subroutine spllt_open_timer



  subroutine spllt_tic(step_name, step_id, thn, timer)
    character(len=*),           intent(in)    :: step_name
    integer,                    intent(in)    :: step_id
    integer,                    intent(in)    :: thn
    type(spllt_timer_t),        intent(inout) :: timer

    timer%start(step_id, thn) = omp_get_wtime()
    timer%steps%name(step_id)  = step_name

  end subroutine spllt_tic



  subroutine spllt_tac(step_id, thn, timer)
    integer,                    intent(in)    :: step_id
    integer,                    intent(in)    :: thn
    type(spllt_timer_t),        intent(inout) :: timer


    timer%swap(thn) = omp_get_wtime()
    timer%start(step_id, thn)  = timer%swap(thn) - timer%start(step_id, thn)

    ! Compute elapsed time
    ! Compute min/max
    timer%time(step_id, thn)   = timer%time(step_id, thn) + &
      timer%start(step_id, thn)
    if(timer%time_min(step_id, thn) .gt. timer%start(step_id, thn)) then
      timer%time_min(step_id, thn) = timer%start(step_id, thn)
    end if
    if(timer%time_max(step_id, thn) .lt. timer%start(step_id, thn)) then
      timer%time_max(step_id, thn) = timer%start(step_id, thn)
    end if

    ! Counter
    timer%ncall(step_id, thn) = timer%ncall(step_id, thn) + 1
    ! Prepare in case of another call to tac
    timer%start(step_id, thn) = timer%swap(thn)

  end subroutine spllt_tac



  subroutine spllt_tac_flop(step_id, thn, timer, flop)
    integer,                    intent(in)    :: step_id
    integer,                    intent(in)    :: thn
    type(spllt_timer_t),        intent(inout) :: timer
    double precision,           intent(in)    :: flop

    double precision :: rate

    timer%swap(thn) = omp_get_wtime()
    timer%start(step_id, thn)  = timer%swap(thn) - timer%start(step_id, thn)


    ! Compute elapsed time
    ! Compute min/max
    timer%time(step_id, thn)   = timer%time(step_id, thn) + &
      timer%start(step_id, thn)
    if(timer%time_min(step_id, thn) .gt. timer%start(step_id, thn)) then
      timer%time_min(step_id, thn) = timer%start(step_id, thn)
    end if
    if(timer%time_max(step_id, thn) .lt. timer%start(step_id, thn)) then
      timer%time_max(step_id, thn) = timer%start(step_id, thn)
    end if

    ! Flop rate count
    rate = flop / timer%start(step_id, thn)
    timer%flop(step_id, thn)   = timer%flop(step_id, thn) + flop
   !print '(a, a, a, es10.2, a, es10.2)', "Flop/s of ", &
   !  trim(timer%steps%name(step_id)), " : ", rate,     &
   !  " sum : ", timer%flop(step_id, thn)
    if(timer%flop_min(step_id, thn) .gt. rate) then
      timer%flop_min(step_id, thn) = rate
    end if
    if(timer%flop_max(step_id, thn) .lt. rate) then
      timer%flop_max(step_id, thn) = rate
    end if

    ! Counter
    timer%ncall(step_id, thn) = timer%ncall(step_id, thn) + 1
    ! Prepare in case of another call to tac
    timer%start(step_id, thn) = timer%swap(thn)

  end subroutine spllt_tac_flop



  subroutine spllt_ftac(step_id, thn, timer)
    integer,                    intent(in)    :: step_id
    integer,                    intent(in)    :: thn
    type(spllt_timer_t),        intent(inout) :: timer


    timer%swap(thn) = omp_get_wtime()
    timer%start(step_id, thn)  = timer%swap(thn) - timer%start(step_id, thn)

    ! Compute elapsed time
    timer%time(step_id, thn)   = timer%time(step_id, thn) + &
      timer%start(step_id, thn)
    ! Counter
    timer%ncall(step_id, thn) = timer%ncall(step_id, thn) + 1
    ! Prepare in case of another call to tac
    timer%start(step_id, thn) = timer%swap(thn)

  end subroutine spllt_ftac



  subroutine spllt_close_ftimer(thn, timer)
    integer,                    intent(in)    :: thn
    type(spllt_timer_t),        intent(inout) :: timer
    
    integer :: step_id

    step_id = 0

    timer%swap(thn)           = omp_get_wtime()
    timer%start(step_id, thn) = timer%swap(thn) - timer%start(step_id, thn)

    ! Compute elapsed time
    ! Compute min/max
    timer%time(step_id, thn)   = timer%time(step_id, thn) + &
      timer%start(step_id, thn)
    ! Counter
    timer%ncall(step_id, thn) = timer%ncall(step_id, thn) + 1
    
  end subroutine spllt_close_ftimer



  subroutine spllt_close_timer(thn, timer)
    integer,                    intent(in)    :: thn
    type(spllt_timer_t),        intent(inout) :: timer
    
    integer :: step_id

    step_id = 0

    timer%swap(thn)           = omp_get_wtime()
    timer%start(step_id, thn) = timer%swap(thn) - timer%start(step_id, thn)

    ! Compute elapsed time
    ! Compute min/max
    timer%time(step_id, thn)   = timer%time(step_id, thn) + &
      timer%start(step_id, thn)
    if(timer%time_min(step_id, thn) .gt. timer%start(step_id, thn)) then
      timer%time_min(step_id, thn) = timer%start(step_id, thn)
    end if
    if(timer%time_max(step_id, thn) .lt. timer%start(step_id, thn)) then
      timer%time_max(step_id, thn) = timer%start(step_id, thn)
    end if
    ! Counter
    timer%ncall(step_id, thn) = timer%ncall(step_id, thn) + 1
    
  end subroutine spllt_close_timer



  subroutine spllt_close_timer_flop(thn, timer, flop)
    integer,                    intent(in)    :: thn
    type(spllt_timer_t),        intent(inout) :: timer
    double precision,           intent(in)    :: flop
    
    integer :: step_id
    double precision :: rate

    step_id = 0

    timer%swap(thn)           = omp_get_wtime()
    timer%start(step_id, thn) = timer%swap(thn) - timer%start(step_id, thn)

    ! Compute elapsed time
    ! Compute min/max
    timer%time(step_id, thn)   = timer%time(step_id, thn) + &
      timer%start(step_id, thn)
    if(timer%time_min(step_id, thn) .gt. timer%start(step_id, thn)) then
      timer%time_min(step_id, thn) = timer%start(step_id, thn)
    end if
    if(timer%time_max(step_id, thn) .lt. timer%start(step_id, thn)) then
      timer%time_max(step_id, thn) = timer%start(step_id, thn)
    end if

    ! Flop rate count
    rate = flop / timer%start(step_id, thn)
   !print '(a, es10.2)', "Close timer => GigaFlop/s : ", rate
    timer%flop(step_id, thn)   = timer%flop(step_id, thn) + flop
    if(timer%flop_min(step_id, thn) .gt. rate) then
      timer%flop_min(step_id, thn) = rate
    end if
    if(timer%flop_max(step_id, thn) .lt. rate) then
      timer%flop_max(step_id, thn) = rate
    end if

    ! Counter
    timer%ncall(step_id, thn) = timer%ncall(step_id, thn) + 1
    
  end subroutine spllt_close_timer_flop



  subroutine spllt_print_timer_color(thn, timer)
    integer,              intent(in) :: thn
    type(spllt_timer_t),  intent(in) :: timer

    integer :: i

    i = 0

    if(timer%ncall(i, thn) .gt. 0) then
      print '(a, i3, a, a30, a, es10.2, a, i7, a, es10.2, a, es10.2, a)', &
        achar(27) // '[82m ' //                                           &
        "[Th: ", thn, "] ", trim(timer%steps%name(i)), " : ",             &
        timer%time(i, thn), " s (ncall ", timer%ncall(i, thn), ") [",     &
        timer%time_min(i, thn), ' , ', timer%time_max(i, thn), " ]"       &
        // achar(27) // '[0m' 
      if(timer%flop(i, thn) .gt. 0.0) then
        print '(20x, a33, es10.2, a, es10.2, a, es10.2, a)',          &
          achar(27) // '[33m ' //                                     &
          "Flop_rate : ", timer%flop(i, thn) / timer%time(i, thn),    &
          " [", timer%flop_min(i, thn), ' , ', timer%flop_max(i, thn),&
          " ]"                                                        &
          // achar(27) // '[0m'
      end if
    else
      print '(a, i3, a, a20, a)',                   &
        achar(27) // '[31m ' //                     &
        "[Th: ", thn, "] ", trim('sub steps'), ":"  &
        // achar(27) // '[0m'
    endif

    do i = 1, ubound(timer%ncall, 1)
      if(timer%ncall(i, thn) .gt. 0) then
        if(timer%ncall(i, thn) .eq. 1) then
          print '(20x, a30, a, es10.2, a)',                       &
            achar(27) // '[32m ' //                               &
            trim(timer%steps%name(i)), " : ", timer%time(i, thn), &
            achar(27) // '[0m'
          if(timer%flop(i, thn) .gt. 0.0) then
            print '(20x, a33, es10.2, a)',                              &
              achar(27) // "[33m " //                                   &
              "Flop_rate : ", timer%flop(i, thn) / timer%time(i, thn),  &
              achar(27) //'[0m'
          end if
        else
          print '(20x, a30, a, es10.2, a, i7, a, es10.2, a, es10.2, a)',  &
            achar(27) // '[32m ' //                                       &
            trim(timer%steps%name(i)), " : ", timer%time(i, thn), " s (", &
            timer%ncall(i, thn), ") [ ",                                  &
            timer%time_min(i, thn), ' , ', timer%time_max(i, thn), ' ]'   &
            // achar(27) //'[0m'

          if(timer%flop(i, thn) .gt. 0.0) then
            print '(20x, a33, es10.2, a, es10.2, a, es10.2, a)',              &
              achar(27) // "[33m " //                                         &
              "Flop_rate : ", timer%flop(i, thn) / timer%time(i, thn), " [",  &
              timer%flop_min(i, thn), ' , ', timer%flop_max(i, thn), ' ]'     &
              // achar(27) // '[0m'
          end if
        end if
      end if
    end do
    print *, ""

  end subroutine spllt_print_timer_color



  subroutine spllt_print_timer(thn, timer)
    integer,              intent(in) :: thn
    type(spllt_timer_t),  intent(in) :: timer

    integer :: i

    i = 0

    if(timer%ncall(i, thn) .gt. 0) then
      print '(a, i3, a, a30, a, es10.2, a, i7, a, es10.2, a, es10.2, a)', &
        "[Th: ", thn, "] ", trim(timer%steps%name(i)), " : ",             &
        timer%time(i, thn), " s (ncall ", timer%ncall(i, thn), ") [",     &
        timer%time_min(i, thn), ' , ', timer%time_max(i, thn), " ]"
      if(timer%flop(i, thn) .gt. 0.0) then
        print '(20x, a33, es10.2, a, es10.2, a, es10.2, a)',          &
          "Flop_rate : ", timer%flop(i, thn) / timer%time(i, thn),    &
          " [", timer%flop_min(i, thn), ' , ', timer%flop_max(i, thn),&
          " ]"
      end if
    else
      print '(a, i3, a, a20, a)',                   &
        "[Th: ", thn, "] ", trim('sub steps'), ":"
    endif

    do i = 1, ubound(timer%ncall, 1)
      if(timer%ncall(i, thn) .gt. 0) then
        if(timer%ncall(i, thn) .eq. 1) then
          print '(20x, a30, a, es10.2, a)',                       &
            trim(timer%steps%name(i)), " : ", timer%time(i, thn)
          if(timer%flop(i, thn) .gt. 0.0) then
            print '(20x, a33, es10.2, a)',                              &
              "Flop_rate : ", timer%flop(i, thn) / timer%time(i, thn)
          end if
        else
          print '(20x, a30, a, es10.2, a, i7, a, es10.2, a, es10.2, a)',  &
            trim(timer%steps%name(i)), " : ", timer%time(i, thn), " s (", &
            timer%ncall(i, thn), ") [ ",                                  &
            timer%time_min(i, thn), ' , ', timer%time_max(i, thn), ' ]'

          if(timer%flop(i, thn) .gt. 0.0) then
            print '(20x, a33, es10.2, a, es10.2, a, es10.2, a)',              &
              "Flop_rate : ", timer%flop(i, thn) / timer%time(i, thn), " [",  &
              timer%flop_min(i, thn), ' , ', timer%flop_max(i, thn), ' ]'
          end if
        end if
      end if
    end do
    print *, ""

  end subroutine spllt_print_timer



  subroutine spllt_print_timers(thn)
    integer, optional, intent(in) :: thn

    integer :: t
    integer :: th, th_min, th_max
    type(spllt_p_timer_t), pointer :: p_timer

    if(present(thn)) then
      th_min = thn
      th_max = thn
    else
      th_min = 0
      th_max = all_timers%nthread - 1
    end if

    print *, "=========================================="
    print *, "                  TIMER"
    print *, "=========================================="
    do th = th_min, th_max
      do t = 1, all_timers%ntimer_th
        if(associated(all_timers%timers(t, th)%p)) then
          p_timer => all_timers%timers(t, th)
          if(p_timer%p%status(th) .eq. 1) then
            if(sum(p_timer%p%ncall(:, th)) .gt. 0) then
              call spllt_print_timer(th, p_timer%p)
            end if
          end if
        end if
      end do
    end do
  end subroutine spllt_print_timers

end module timer_mod
