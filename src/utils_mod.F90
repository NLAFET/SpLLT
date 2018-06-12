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

  interface pack_rhs
    module procedure linear_pack_rhs
    module procedure vector_pack_rhs
  end interface pack_rhs

  interface unpack_rhs
    module procedure linear_unpack_rhs
    module procedure vector_unpack_rhs
  end interface unpack_rhs
contains

  subroutine print_darray(array_name, n, val, display)
    use spllt_data_mod
    character(len=*),       intent(in)    :: array_name
    integer,                intent(in)    :: n
    real(wp), dimension(n), intent(in)    :: val
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
        print *, val(i)
      end do
    else
      do i = 1, n
        write(*, fmt="(es10.2)", advance="no") val(i)
      end do
      write(*,*) ""
    end if
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
  

  subroutine check_backward_error(n, ptr, row, val, nrhs, x, b) 
    use spllt_data_mod
    implicit none

    integer, intent(in)                         :: n
    integer, dimension(n+1), intent(in)         :: ptr
    integer, dimension(ptr(n+1)-1), intent(in)  :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in)                         :: nrhs
    real(wp), dimension(n,nrhs), intent(in)     :: x
    real(wp), dimension(n,nrhs), intent(in)     :: b

    logical           :: bwd_error_ok
    integer           :: i
    real(wp)          :: norm_max
    real(wp)          :: res(n, nrhs)
    double precision  :: normRes(nrhs)
    double precision  :: normRHS(nrhs)
    double precision  :: solNorm(nrhs)
    double precision  :: err

    call compute_residual(n, ptr, row, val, nrhs, x, b, res)

    call vector_norm_2(n, res, normRes)
    call vector_norm_2(n,   b, normRHS)
    call vector_norm_2(n,   x, solNorm)
    call matrix_norm_max(n, ptr, row, val, norm_max)

    bwd_error_ok = .true.
    do i = 1, nrhs
      err = normRes(i) / (normRHS(i) + norm_max * solNorm(i))
     !if(normRes(i) / normRHS(i) .gt. 1e-14) then
      if(err .gt. 1e-14) then
        write(0, "(a, i4, a, i4, a, es10.2)") "Wrong Bwd error for ", i, &
          "/", nrhs, " : ", err
        bwd_error_ok = .false.
      end if
    end do
    if(bwd_error_ok) then
      write(0, "(a)") "Backward error... ok"
    end if
  end subroutine check_backward_error


  !*************************************************
  !
  ! This function calculates column of a node we are on  
  integer function calc_col(node, bc)
    use spllt_data_mod    
    implicit none

    type(spllt_node), intent(in) :: node
    type(spllt_block), intent(in) :: bc

    calc_col = (size(node%index)-1)/node%nb + 1 ! no. row blks for node

    calc_col = calc_col - (bc%last_blk - bc%dblk + 1) + 1 ! column of node

  end function calc_col

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

! subroutine print_csc_matrix(m, n, nnz, colPtr, rowInd, val)
!   use spllt_data_mod
!   implicit none
!   integer,  intent(in)  :: m         ! Number of rows
!   integer,  intent(in)  :: n         ! Number of columns
!   integer,  intent(in)  :: nnz       ! Number of entries
!   integer,  intent(in)  :: colPtr(:) ! Elements that index rowInd array
!   integer,  intent(in)  :: rowInd(:) ! Stores the row index
!   real(wp), intent(in)  :: val(:)    ! Entry array in CSC format
!   
!   integer :: i, j, k
!   real(wp), allocatable :: work(:, :)

!   call print_array("colPtr", 15, colPtr, 1)
!   call print_array("rowInd", 30  , rowInd, 1)
!   call print_array("val   ", 30  ,    val, 1)

!   return

!   allocate(work(m, n))
!   work = -9999.8888

!   do i = 1, n
!     do j = colPtr(i), colPtr(i + 1) - 1
!       work(rowInd(j), i) = val(j)
!     end do
!   end do
!  !do j = 1, m
!  !  do i = 1, n
!  !    if(rowInd(colPtr(i) + k) .eq. j) then
!  !      write(*, fmt="(es10.2)", advance="no") val(colPtr(i) + k)
!  !      k = k + 1
!  !    else
!  !      write(*, fmt="(a5)", advance="no") "."
!  !    end if
!  !  end do
!  !end do

!   do i = 1, m
!     do j = 1, n
!       if(work(i,j) .eq. -9999.8888) then
!         write(*, fmt="(a10)", advance="no") "."
!       else
!         write(*, fmt="(es10.2)", advance="no") work(i, j)
!       end if
!     end do
!     print *, ""
!   end do

! end subroutine print_csc_matrix
! 


! subroutine print_csr_matrix(m, n, nnz, rowPtr, colInd, val)
!   use spllt_data_mod
!   implicit none
!   integer,  intent(in)  :: m         ! Number of rows
!   integer,  intent(in)  :: n         ! Number of columns
!   integer,  intent(in)  :: nnz       ! Number of entries
!   integer,  intent(in)  :: rowPtr(:) ! Elements that index rowInd array
!   integer,  intent(in)  :: colInd(:) ! Stores the row index
!   real(wp), intent(in)  :: val(:)    ! Entry array in CSR format
!   
!   integer :: i, j, k

!  !call print_array("rowPtr", m + 1, rowPtr, 1)
!  !call print_array("colInd", nnz  , colInd, 1)
!  !call print_array("val   ", nnz  ,    val, 1)
!   call print_array("rowPtr", 15, rowPtr, 1)
!   call print_array("colInd", 30  , colInd, 1)
!   call print_array("val   ", 30  ,    val, 1)
!  !do i = 1, m
!  !  k = 0
!  !  do j = 1, n
!  !    if(j .eq. colInd(rowPtr(i) + k)) then
!  !      write(*, fmt="(es10.2)", advance="no") val(rowPtr(i) + k)
!  !      k = k + 1
!  !    else
!  !      write(*, fmt="(a10)", advance="no") "."
!  !    end if
!  !  end do
!  !  print *, ""
!  !end do

! end subroutine print_csr_matrix

  subroutine linear_pack_rhs(nrhs, rhs, n, ldr, bdr, rhs_pack, st)
    use spllt_data_mod
    implicit none
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: n
    integer,  intent(in)    :: ldr
    integer,  intent(in)    :: bdr
    real(wp), intent(in)    :: rhs(n,nrhs)
    real(wp), intent(inout) :: rhs_pack(n*nrhs)
    integer,  intent(out)   :: st

    integer(long) :: size_nrhs
    integer       :: nblk, blk, r
    integer       :: last_entry

    size_nrhs   = int(n, long) * nrhs
    last_entry  = mod(n, bdr)
    nblk        = n / bdr

!   print *, "Nblk : ", nblk, 'with a rest of ', last_entry

    do blk = 0, nblk - 1
      do r = 1, nrhs
!       print *, 'Fill rhs_pack from ', 1 + ldr * blk + bdr * (r - 1),   &
!         ' to ', ldr * blk + bdr * r, 'with', 1 + bdr * blk, 'to', &
!         bdr * (blk + 1)
        rhs_pack(1 + ldr * blk + bdr * (r - 1) : ldr * blk + bdr * r ) =&
          rhs(1 + bdr * blk : bdr * (blk + 1), r)
      end do
    end do
    if (last_entry .gt. 0) then
      do r = 1, nrhs
        if( r .lt. nrhs) then
!         print *, 'Fill rhs_pack from ', 1 + ldr * blk + bdr * (r - 1),   &
!           ' to ', ldr * blk + last_entry * r, ' with ', 1 + bdr * blk,&
!           ' to ', n
          rhs_pack(1 + ldr * blk + bdr * (r - 1) : ldr * blk + last_entry * r) = &
            rhs(1 + bdr * blk : n, r)
        else
!         print *, 'Fill xil from ', 1 + ldr * blk + last_entry * (r - 1), &
!           ' to ', size_nrhs, ' with ', 1 + bdr * blk, ' to ', n
          rhs_pack(1 + ldr * blk + last_entry * (r - 1) : size_nrhs) =&
            rhs(1 + bdr * blk : n, r)
        end if
      end do
    end if

  end subroutine linear_pack_rhs



  subroutine linear_unpack_rhs(nrhs, rhs_pack, n, ldr, bdr, rhs, st)
    use spllt_data_mod
    implicit none
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: n
    integer,  intent(in)    :: ldr
    integer,  intent(in)    :: bdr
    real(wp), intent(in)    :: rhs_pack(n*nrhs)
    real(wp), intent(inout) :: rhs(n,nrhs)
    integer,  intent(out)   :: st

    integer(long) :: size_nrhs
    integer       :: nblk, blk, r
    integer       :: last_entry

    size_nrhs = int(n, long) * nrhs
    last_entry  = mod(n, bdr)
    nblk      = n / bdr

    do blk = 0, nblk - 1!- merge(1, 0, last_blk .gt. 0)
      do r = 1, nrhs
!       print *, 'RHS ', r, 'Fill x from ', 1 + bdr * blk, ' to ', &
!         bdr * (blk + 1)
        rhs(1 + bdr * blk : bdr * (blk + 1), r) = &
          rhs_pack(1 + ldr * blk + bdr * (r - 1) : ldr * blk + bdr * r )
      end do
    end do
    if (last_entry .gt. 0) then
      do r = 1, nrhs
        if( r .lt. nrhs) then
!         print *, 'RHS ', r, 'Fill x from ', 1 + bdr * blk, ' to ', n, &
!           'with', 1 + ldr * blk + last_entry * (r - 1), 'to', &
!           ldr * blk + last_entry * r
          rhs(1 + bdr * blk : n, r) = &
            rhs_pack(1 + ldr * blk + last_entry * (r - 1) : ldr * blk + last_entry * r )
        else
!         print *, 'RHS ', r, 'Fill x from ', 1 + bdr * blk, ' to ', n, &
!           'with', 1 + ldr * blk + last_entry * (r - 1), 'to', size_nrhs
          rhs(1 + bdr * blk : n, r) = &
            rhs_pack(1 + ldr * blk + last_entry * (r - 1) :  size_nrhs )
        end if
      end do
    end if

  end subroutine linear_unpack_rhs



  subroutine vector_pack_rhs(nrhs, rhs, n, v, rhs_pack, st)
    use spllt_data_mod
    implicit none
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: n
    integer,  intent(in)    :: v(:)
    real(wp), intent(in)    :: rhs(n,nrhs)
    real(wp), intent(inout) :: rhs_pack(n*nrhs)
    integer,  intent(out)   :: st

    integer(long) :: size_nrhs
    integer       :: nblk, blk, r
    integer       :: rowPtr, bdr

    size_nrhs   = int(n, long) * nrhs
    nblk        = size(v) - 1

    print *, "Nblk : ", nblk

    do blk = 1, nblk
      rowPtr  = v(blk) * nrhs
      bdr     = v(blk + 1) - v(blk)
      do r = 0, nrhs - 1
        print *, 'Fill rhs_pack from ', 1 + rowPtr + bdr * r,       &
          ' to ', rowPtr + bdr * (r + 1), 'with', v(blk) + 1, 'to', &
          v(blk) + bdr
        rhs_pack(1 + rowPtr + bdr * r : rowPtr + bdr * (r + 1)) =&
          rhs(v(blk) + 1: v(blk) + bdr, r + 1)
      end do
    end do

  end subroutine vector_pack_rhs



  subroutine vector_unpack_rhs(nrhs, rhs_pack, n, v, rhs, st)
    use spllt_data_mod
    implicit none
    integer,  intent(in)    :: nrhs
    integer,  intent(in)    :: n
    integer,  intent(in)    :: v(:)
    real(wp), intent(in)    :: rhs_pack(n*nrhs)
    real(wp), intent(inout) :: rhs(n,nrhs)
    integer,  intent(out)   :: st

    integer(long) :: size_nrhs
    integer       :: nblk, blk, r
    integer       :: bdr, rowPtr

    size_nrhs   = int(n, long) * nrhs
    nblk        = size(v) - 1

    print *, "Nblk : ", nblk

    do blk = 1, nblk
      rowPtr  = v(blk) * nrhs
      bdr     = v(blk + 1) - v(blk)
      do r = 0, nrhs - 1
        print *, 'Fill x from ', 1 + rowPtr + bdr * r,       &
          ' to ', rowPtr + bdr * (r + 1), 'with', v(blk) + 1, 'to', &
          v(blk) + bdr
        rhs(v(blk) + 1: v(blk) + bdr, r + 1)  = &
          rhs_pack(1 + rowPtr + bdr * r : rowPtr + bdr * (r + 1))
      end do
    end do

  end subroutine vector_unpack_rhs



end module utils_mod
