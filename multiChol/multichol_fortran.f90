subroutine multichol_fort_il(num_vec, m, n, a, lda, info) bind(C)
   use iso_c_binding
   implicit none

   integer(C_INT), value :: num_vec
   integer(C_INT), value :: m
   integer(C_INT), value :: n
   real(C_DOUBLE), dimension(*), intent(inout) :: a
   integer(C_INT), value :: lda
   integer(C_INT), dimension(*), intent(inout) :: info

   integer, parameter :: veclen = 4

   integer :: i
   integer :: sublda, set, col, row, upd
   integer :: aoffset

   real(C_DOUBLE), dimension(veclen) :: diag, rdiag, rval, a_ij, a_ik, a_jk
   !dir$ assume_aligned a:64

   ! AVX version
   sublda = veclen*lda
   do set = 1, num_vec
      aoffset = (set-1)*veclen*lda*n
      do col = 1, n
         ! Take square root of diagonal (or error if non-positive)
         diag(:) = a(aoffset+(col-1)*(sublda+veclen)+1:aoffset+(col-1)*(sublda+veclen)+veclen)
         if(any(diag(:) < 0.0)) then
            ! A diagonal entry is less than zero: flag in info, but apply ops
            ! regardless
            do i = 1, veclen
               if(diag(i) < 0.0) then
                  if(info((set-1)*veclen+i).eq.0) &
                     info((set-1)*veclen+i) = col
               end if
            end do
         end if
         diag(:) = sqrt(diag(:))
         rdiag(:) = 1.0 / diag(:)
         a(aoffset+(col-1)*(sublda+veclen)+1:aoffset+(col-1)*(sublda+veclen)+veclen) = diag(:)
         ! Divide column by a
         do row = col+1, m
            rval(:) = a(aoffset+(col-1)*sublda+(row-1)*veclen+1:aoffset+(col-1)*sublda+(row-1)*veclen+veclen)
            rval(:) = rval(:) * rdiag(:)
            a(aoffset+(col-1)*sublda+(row-1)*veclen+1:aoffset+(col-1)*sublda+(row-1)*veclen+veclen) = rval(:)
         end do
         ! Apply outer product update to remaining columns
         do upd = col+1, n
            ! NB: Aim to do a_ij -= a_ik * a_jk
            !     i = row
            !     j = upd
            !     k = col
            a_jk(:) = a(aoffset+(col-1)*sublda+(upd-1)*veclen+1:aoffset+(col-1)*sublda+(upd-1)*veclen+veclen)
            do row = upd, n
               a_ij(:) = a(aoffset+(upd-1)*sublda+(row-1)*veclen+1:aoffset+(upd-1)*sublda+(row-1)*veclen+veclen)
               a_ik(:) = a(aoffset+(col-1)*sublda+(row-1)*veclen+1:aoffset+(col-1)*sublda+(row-1)*veclen+veclen)
               a_ij(:) = a_ij(:) - a_ik(:) * a_jk(:)
               a(aoffset+(upd-1)*sublda+(row-1)*veclen+1:aoffset+(upd-1)*sublda+(row-1)*veclen+veclen) = a_ij(:)
            end do ! row
         end do ! upd
      end do ! col
   end do ! set
end subroutine multichol_fort_il
