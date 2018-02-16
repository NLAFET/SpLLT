module get_wtime_mod
  implicit none

  interface
     function omp_get_wtime() bind(C)
       use iso_c_binding
       real(c_double) :: omp_get_wtime 
     end function omp_get_wtime
  end interface

contains

! function get_fwtime()
!   real :: get_fwtime
!   character(len=10) :: time
!   call date_and_time(TIME=time)
!   read(time, *) get_fwtime
! end function get_fwtime

end module get_wtime_mod
