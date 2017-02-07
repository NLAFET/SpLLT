!!! papif_omp - a simple wrapper for using PAPI with Fortran/OpenMP.
!!! Assumes that all OpenMP threads are bound and that the binding does not change in run-time.

module papif_omp

  implicit none

contains

  subroutine papif_omp_init(info)

    use, intrinsic :: iso_c_binding
    include 'f90papi.h'
    implicit none

    integer(c_int), intent(out) :: info

    info = int(PAPI_VER_CURRENT, c_int)

    call PAPIF_library_init(info)
    if (info .ne. int(PAPI_VER_CURRENT, c_int)) return

    call PAPIF_thread_init(papif_omp_gettid, info)
    if (info .ne. int(PAPI_OK, c_int)) return

  end subroutine papif_omp_init

  function papif_omp_gettid()

    use omp_lib
    use, intrinsic :: iso_c_binding

    implicit none

    integer(c_long_long) :: papif_omp_gettid

    papif_omp_gettid = int(omp_get_thread_num(), c_long_long)

  end function papif_omp_gettid

  subroutine papif_omp_register(nt, info)

    use omp_lib
    use, intrinsic :: iso_c_binding
    ! not really needed
    ! include 'f90papi.h'
    implicit none

    integer, intent(in) :: nt
    integer(c_int), intent(out) :: info

    integer :: thr
    integer(c_int) :: check

    if (nt .eq. 0) then
       thr = max(omp_get_max_threads(), 1)
    else if (nt .gt. 0) then
       thr = min(omp_get_max_threads(), nt)
    else
       ! positive info, to avoid confusion with PAPI error codes
       info = 1_c_int
       return
    end if

    check = 0_c_int
    !$omp parallel num_threads(thr), reduction(min:check)
    call PAPIF_register_thread(check)
    !$omp end parallel
    info = check

  end subroutine papif_omp_register

  subroutine papif_omp_unregister(nt, info)

    use omp_lib
    use, intrinsic :: iso_c_binding
    ! not really needed
    ! include 'f90papi.h'
    implicit none

    integer, intent(in) :: nt
    integer(c_int), intent(out) :: info

    integer :: thr
    integer(c_int) :: check

    if (nt .eq. 0) then
       thr = max(omp_get_max_threads(), 1)
    else if (nt .gt. 0) then
       thr = min(omp_get_max_threads(), nt)
    else
       ! positive info, to avoid confusion with PAPI error codes
       info = 1_c_int
       return
    end if

    check = 0_c_int
    !$omp parallel num_threads(thr), reduction(min:check)

!!! PAPIF_unregster_thread
!!! NOT unregister (typo)!

    call PAPIF_unregster_thread(check)
    !$omp end parallel
    info = check

  end subroutine papif_omp_unregister

  subroutine papif_omp_free(info)

    use, intrinsic :: iso_c_binding
    include 'f90papi.h'
    implicit none

    integer(c_int), intent(out) :: info

    call PAPIF_shutdown
    info = int(PAPI_OK, c_int)

  end subroutine papif_omp_free

end module papif_omp
