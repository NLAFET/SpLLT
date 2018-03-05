module utils_mod

  use spllt_data_mod

  interface print_array
    module procedure print_darray
    module procedure print_iarray
  end interface print_array

contains

  subroutine print_darray(array_name, n, val)
  ! use spllt_data_mod
    character(len=*)                      :: array_name
    integer,                intent(in)    :: n
    real(wp), dimension(n), intent(in)    :: val

    integer :: i
    print '(a)', array_name
    do i = 1, n
     !print '(f20.8)', val(i)
      print *, val(i)
    end do
  end subroutine print_darray

  subroutine print_iarray(array_name, n, val, display)
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

  function contain(fkeep, blk1, blk2) result(isIn)
    use spllt_solve_dep_mod, only : getPointerBlkIndex
    
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: blk1
    integer, intent(in)                   :: blk2

    integer :: j, k
    integer, pointer :: p_blk1_index(:), p_blk2_index(:)
    logical :: isIn

    j = 1
    k = 1
    isIn = .false.

    call getPointerBlkIndex(fkeep, blk1, p_blk1_index)
    call getPointerBlkIndex(fkeep, blk2, p_blk2_index)

    do while(j .le. size(p_blk1_index) .and. k .le. size(p_blk2_index))
      if(p_blk1_index(j) .lt. p_blk2_index(k)) then
        j = j + 1
      else if(p_blk1_index(j) .gt. p_blk2_index(k)) then
        k = k + 1
      else
        isIn = .true.
        return
      end if
    end do
  end function contain

  !Compute res = b - Ax 
  subroutine compute_residual(n, ptr, row, val, nrhs, x, b, res)
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in) :: nrhs
    real(wp), dimension(n,nrhs), intent(in) :: x
    real(wp), dimension(n,nrhs), intent(in) :: b
    real(wp), dimension(n,nrhs), intent(inout) :: res

    integer :: i, j, k, r
    ! Find the residual
    !allocate(res(n,nrhs))
    res(:,:) = 0
    call compute_Ax(n, ptr, row, val, nrhs, x, res)
    res = b - res
  end subroutine compute_residual
  
  !Compute Ax
  subroutine compute_Ax(n, ptr, row, val, nrhs, x, res)
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in) :: nrhs
    real(wp), dimension(n,nrhs), intent(in) :: x
    real(wp), dimension(n,nrhs), intent(inout) :: res

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
!   use spllt_data_mod
    integer, intent(in)                         :: n
    integer, dimension(n+1), intent(in)         :: ptr
    integer, dimension(ptr(n+1)-1), intent(in)  :: row    ! Unused
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    real(wp), intent(out)                       :: norm

    real(wp) :: i, j, sum_col

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
!   use spllt_data_mod
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
