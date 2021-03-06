!> \file
!> \copyright 2018 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Florent Lopez
!> \author    Sebastien Cayrols
module trace_mod
  use get_wtime_mod
  use iso_c_binding
  implicit none

  type event_type
     integer :: id, thn
     real(kind(1.d0)) :: start, stop
  end type event_type

  integer, save               :: trace_nth
  integer, parameter          :: trace_nstep=2
  logical, allocatable, save  :: pendings(:,:)
  logical, allocatable, save  :: pauses(:,:)
  real(c_double), save        :: timezero, start, stop
  integer, parameter          :: maxevents=15000, maxtypes=20
  real(c_double), save        :: ttimes(1:maxtypes)
  integer                     :: trace_cpt
  integer, save               :: ids(1:maxtypes)=(/(trace_cpt, &
    trace_cpt=1,maxtypes)/)
  type(event_type), allocatable, save :: events(:,:)
  real(c_double), allocatable, save   :: stops(:)
  real(c_double), allocatable, save   :: starts(:,:)

  character(len=20), save     :: labels(maxtypes)
  character(len=7)            :: colors(maxtypes)
  integer, save               :: nevtype
  integer, allocatable , save :: nevents(:)

contains

  subroutine trace_init(nth)
    use get_wtime_mod
    implicit none
    integer :: nth ! number of threads/workers involved in the
                   ! execution
    
    trace_nth = nth
    allocate(pendings(maxtypes, -trace_nstep : trace_nth - 1))
    allocate(pauses  (maxtypes, -trace_nstep : trace_nth - 1))
    allocate(nevents(           -trace_nstep : trace_nth - 1))
    allocate(stops(             -trace_nstep : trace_nth - 1))
    allocate(starts(maxtypes,   -trace_nstep : trace_nth - 1))
    pendings(:,:) = .false.
    pauses(:,:)   = .false.
    nevtype     = 0
    nevents(:)  = 0
    ttimes(:)   = 0
    allocate(events(            -trace_nstep : trace_nth - 1, maxevents))
    colors(1:(trace_nth+trace_nstep))= (/&
!     '#d38d5f', '#ffdd55', &
!     '#8dd35f', '#80b3ff', &
!     '#e580ff', '#ac9d93', &
!     '#bcd35f', '#a61e22', &
!     '#5542d7', '#2F4F4F', &
!     '#6600ff', '#ff0066', &
!     '#ff00ff'             &
      '#e6194b', '#3cb44b', &
      '#ffe119', '#0082c8', &
      '#f58231', '#911eb4', &
      '#46f0f0', '#f032e6', &
      '#d2f53c', '#fabebe', &
      '#008080', '#e6beff', &
      '#aa6e28', '#fffac8', &
      '#800000', '#aaffc3', &
      '#808000', '#ffd8b1', &
      '#000080', '#808080', &
      '#ff00ff', '#6600ff'  &
      /)
    timezero = omp_get_wtime()
    return

  end subroutine trace_init

  integer function trace_create_event_id(label)
    implicit none
    character, intent(in) :: label*(*)

    nevtype                       = nevtype + 1
    trace_create_event_id         = nevtype
    labels(trace_create_event_id) = label

  end function trace_create_event_id

  subroutine trace_create_event(label, id)
    implicit none
    character :: label*(*)
    integer   :: id

    nevtype = nevtype + 1
    id      = nevtype
    labels(id) = label

    return
    
  end subroutine trace_create_event

  subroutine trace_event_start(id, thn)
    use get_wtime_mod
    implicit none

    integer :: id, thn, old_id

    if(id .eq. 0) return

    if(pendings(id, thn)) then
      write(*,'("Tracing error!!! events nesting not supported")')
      write(*,*) 'Th ', thn, ' already recorded as doing trace ', id
      print *, "PENDINGS : ", pendings(:, thn)
      return
    end if

    if(any(pendings(:, thn))) then
      if(count(pendings(:, thn)) .gt. 1) then
        write(*,'("Tracing error!!! count events not supported")')
        write(*,*) 'Th ', thn, ' already recorded as doing trace ', &
          sum(ids, mask=pendings(:,thn))
        return
      end if
          
      old_id = sum(ids, mask=pendings(:,thn))
      call trace_event_stop(old_id, thn)
      pauses(old_id,thn) = .true.
    end if

    pendings(id, thn) = .true.
    starts(id, thn) = omp_get_wtime()

    return
  end subroutine trace_event_start



  subroutine trace_event_stop(id, thn)
    use get_wtime_mod
    implicit none

    integer :: id, thn, old_id

    if(id .eq. 0) return
    
    stops(thn)                                = omp_get_wtime()
    nevents(thn)                              = nevents(thn) + 1
    events(thn, min(nevents(thn),maxevents))  = event_type(id, thn, starts(id, thn)-timezero, stops(thn)-timezero)
    ttimes(id)                                = ttimes(id)+stops(thn)-starts(id, thn)
    pendings(id, thn)                         = .false.

    if(any(pauses(:, thn))) then
      if(count(pauses(:, thn)) .gt. 1) then
        write(*,'("Tracing error!!! events nesting not supported")')
        write(*,*) 'Th ', thn, ' already recorded as doing trace ', &
          sum(ids, mask=pauses(:,thn))
        print *, "PAUSES : ", pauses(:, thn)
        return
      end if
          
      old_id = sum(ids, mask=pauses(:,thn))
      starts(old_id, thn) = stops(thn)
      pauses(old_id,thn) = .false.
    end if

    return
  end subroutine trace_event_stop

  subroutine trace_log_dump_paje(ofile)

    character :: ofile*(*)

    integer :: i, j
    real(kind(1.d0)) :: tottime, r, g, b
    real(kind(1.d0)), parameter :: h=20.d0, scale=10000

    
! #ifndef paje
    
!     open(4, file=ofile//".svg", action='write')

!     write(4,'("<svg>")')
!     do i=0, maxth-1
!        do j=1, min(nevents(i),maxevents)
!           write(4,'("<rect style=""fill:",a7,";stroke:#000000;stroke-width:0;fill-opacity:1"" '//&
!                & 'height=""",f5.1,""" width=""",f10.2,""" y=""",f10.2,""" x=""",f10.2,""" />")')&
!                & colors(events(i,j)%id), h, &
!                & (events(i,j)%stop-events(i,j)%start)*scale, &
! !                & h*(i-1), events(i,j)%start*scale
!                & real(h)*real(maxth-i), events(i,j)%start*scale
!        end do
!     end do
    
!     tottime = sum(ttimes)

!     do i=1, nevtype
!        write(4,'("<text x=""0"" y=""",f10.2,""" font-size=""",i3,""" fill=""",a7,""">",a20," -- ",f4.1,"%</text>")')&
!             & real(h)*real(maxth+i+2),floor(h),colors(i),labels(i),(ttimes(i)/tottime)*100
!     end do

!     write(4,'("</svg>")')

! #else

    open(4, file=ofile//".paje", action='write')

    write(4,'("%EventDef PajeDefineContainerType 1")')
    write(4,'("% Alias string ")')
    write(4,'("% ContainerType string ")')
    write(4,'("% Name string ")')
    write(4,'("%EndEventDef")')
    write(4,'("%EventDef PajeDefineStateType 3")')
    write(4,'("% Alias string ")')
    write(4,'("% ContainerType string ")')
    write(4,'("% Name string ")')
    write(4,'("%EndEventDef ")')
    write(4,'("%EventDef PajeDefineEntityValue 6")')
    write(4,'("% Alias string  ")')
    write(4,'("% EntityType string  ")')
    write(4,'("% Name string  ")')
    write(4,'("% Color color ")')
    write(4,'("%EndEventDef  ")')
    write(4,'("%EventDef PajeCreateContainer 7")')
    write(4,'("% Time date  ")')
    write(4,'("% Alias string  ")')
    write(4,'("% Type string  ")')
    write(4,'("% Container string  ")')
    write(4,'("% Name string  ")')
    write(4,'("%EndEventDef  ")')
    write(4,'("%EventDef PajeDestroyContainer 8")')
    write(4,'("% Time date  ")')
    write(4,'("% Name string  ")')
    write(4,'("% Type string  ")')
    write(4,'("%EndEventDef  ")')
    write(4,'("%EventDef PajeSetState 10")')
    write(4,'("% Time date  ")')
    write(4,'("% Type string  ")')
    write(4,'("% Container string  ")')
    write(4,'("% Value string  ")')
    write(4,'("%EndEventDef ")')
    write(4,'("1 CT_Prog   0       ''Program''")')
    write(4,'("1 CT_Thread CT_Prog ''Thread''")')
    write(4,'("3 ST_ThreadState CT_Thread ''Thread State''")')
    do i=1, nevtype
       read(colors(i)(2:3),'(z2)')j
       r =  real(j,kind(1.d0))/255.d0
       read(colors(i)(4:5),'(z2)')j
       g =  real(j,kind(1.d0))/255.d0
       read(colors(i)(6:7),'(z2)')j
       b =  real(j,kind(1.d0))/255.d0
       write(4,'("6 ",a20," ST_ThreadState ''",a20,"''  ''",3(f8.6,x),"''")')labels(i),labels(i),r,g,b
    end do

    write(4,'("6 idle ST_ThreadState ''Idle''  ''255 255 255''")')

    write(4,'("7 0.000000 C_Prog CT_Prog 0 ''Programme''")')
    
    do i = -trace_nstep, -1
      if(nevents(i) .gt. 0) then
        write(4,'("7  0.000000 MASTER",i3.3," CT_Thread C_Prog ''Master ",i3,"''")')-i,-i
      end if
    end do
    do i=0, trace_nth-1
      if(nevents(i) .gt. 0) then
        write(4,'("7  0.000000 C_Thread",i3.3," CT_Thread C_Prog ''Thread ",i3,"''")')i,i
      end if
    end do
            
 
    do i = -trace_nstep, -1
      do j=1, min(nevents(i),maxevents)
        write(4,'("10 ",f15.6," ST_ThreadState MASTER",i3.3,x,a20)')events(i,j)%start, -i, labels(events(i,j)%id)
        write(4,'("10 ",f15.6," ST_ThreadState MASTER",i3.3," idle")')events(i,j)%stop, -i
      end do
    end do

    do i=0, trace_nth-1
      do j=1, min(nevents(i),maxevents)
        write(4,'("10 ",f15.6," ST_ThreadState C_Thread",i3.3,x,a20)')events(i,j)%start, i, labels(events(i,j)%id)
        write(4,'("10 ",f15.6," ST_ThreadState C_Thread",i3.3," idle")')events(i,j)%stop, i
      end do
    end do


    do i = -trace_nstep, -1
      if(nevents(i) .gt. 0) then
        write(4,'("8 ",f15.6," MASTER",i3.3," CT_Thread")') maxval(events(:,:)%stop), -i
      end if
    end do
    do i=0, trace_nth-1
      if(nevents(i) .gt. 0) then
        write(4,'("8 ",f15.6," C_Thread",i3.3," CT_Thread")')maxval(events(:,:)%stop),i
      end if
    end do
    write(4,'("8 ",f15.6," C_Prog CT_Prog")')maxval(events(:,:)%stop)
    

! #endif

    close(4)

    deallocate(events)

    return
    
  end subroutine trace_log_dump_paje

end module trace_mod
