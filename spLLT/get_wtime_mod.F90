module get_wtime_mod
  implicit none

  interface
     function get_wtime() bind(C)
       use iso_c_binding
       real(c_double) :: get_wtime 
     end function get_wtime
  end interface

end module get_wtime_mod
