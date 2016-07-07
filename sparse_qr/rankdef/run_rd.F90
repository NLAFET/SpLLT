#include "qrm_common.h"

program run_rd

  use dqrm_mod
  implicit none

  type(dqrm_spmat_type)          :: qrm_mat, qrm_mat_tilde
  character(len=100) :: matfile = ''
  character(len=200) :: argval
  integer                        :: info, ounit
  integer                        :: nth
  integer                        :: nrhs
  character                      :: transp
  ! regularization parameters
  real(kind(1.d0)) :: sigma, esigma
  !solve
  integer :: i
  real(kind(1.d0)), allocatable, target :: b(:,:), x(:,:), r(:,:)
  real(kind(1.d0)), allocatable, target :: b_tilde(:,:)
  real(kind(1.d0))                      :: anrm
  real(kind(1.d0)), allocatable         :: rnrm(:), onrm(:), bnrm(:), xnrm(:)
  integer                        :: iseed(4)=(/1,1,1,1/)
  ! timing
  real(kind(1.d0))               :: t1, ta, tf, ts
  
  ounit = 6 ! standard output
  nth = 2
  
  call qrm_init(nth)

  ! initialize the matrix data structure. 
  call qrm_spmat_init(qrm_mat)

  ! set parameters
  ! get matrix
  call get_command_argument(1, argval)
  read( argval, *)matfile

  ! sigma parameter
  call get_command_argument(2, argval)
  read( argval, *)esigma
  sigma = 10**esigma
  
  ! output
  call qrm_set(qrm_mat, 'qrm_ounit', ounit)
  
  ! data partitioning parameters
  call qrm_set(qrm_mat, 'qrm_nb', 128)
  call qrm_set(qrm_mat, 'qrm_mb', 256)
  call qrm_set(qrm_mat, 'qrm_ib', 64)
  
  call qrm_set(qrm_mat, 'qrm_keeph', 2) ! keep factors for solve

  call qrm_readmat(matfile, qrm_mat, .true., info)
  __QRM_INFO_CHECK(info,'qrm_test','qrm_readmat',10)

  if(ounit.gt.0) write(ounit,'("Matrix ready. M:",i7,"  N:",i7,"  NZ:",i10)')qrm_mat%m, qrm_mat%n, qrm_mat%nz
  if(ounit.gt.0) write(ounit,'(" ")')

  ! set up regularized matrix
  call qrm_spmat_copy(qrm_mat, qrm_mat_tilde, .true.)
  
  ! regularize
  ! sigma = 1e-10
  call dqrm_tikhonov(qrm_mat_tilde, sigma)
  if(ounit.gt.0) write(ounit,'("Regularized Matrix. M:",i7,"  N:",i7,"  NZ:",i10)')qrm_mat%m, qrm_mat%n, qrm_mat%nz
  if(ounit.gt.0) write(ounit,'(" ")')

  if(qrm_mat%m .ge. qrm_mat%n) then
     transp='n'
  else
     if(ounit.gt.0) write(ounit,'("Transpose")')
     transp='t'
  end if

  ! analyse
  call qrm_analyse(qrm_mat, transp, info)
  if(ounit.gt.0) write(ounit,'("  Estimated total flops at facto               : ",i20)') &
       & qrm_mat%gstats(qrm_e_facto_flops_)

  call qrm_analyse(qrm_mat_tilde, transp, info)

  if(ounit.gt.0) write(ounit,'("  Estimated total flops at facto (regularized) : ",i20)') &
       & qrm_mat_tilde%gstats(qrm_e_facto_flops_)

  ! factorize
  t1 = qrm_swtime()
  call qrm_factorize(qrm_mat_tilde, transp, info)
  tf = qrm_swtime()-t1

  ! solve
  if(ounit.gt.0) write(ounit,'("Starting Solve")')
  nrhs = 1
  call qrm_alloc(b, qrm_mat%m, nrhs)
  ! modified rhs corresponding to regularized problem 
  !
  ! b_tilde = |b|
  !           |0|
  !
  call qrm_alloc(b_tilde, qrm_mat_tilde%m, nrhs)
  call qrm_alloc(r, qrm_mat_tilde%m, nrhs)
  call qrm_alloc(x, qrm_mat%n, nrhs)
     
  call qrm_alloc(bnrm, nrhs)
  call qrm_alloc(rnrm, nrhs)
  call qrm_alloc(xnrm, nrhs)
  call qrm_alloc(onrm, nrhs)
  
  call dlarnv(2, iseed, size(b), b(1,1))
  r = b
  b_tilde = 0
  b_tilde(1:qrm_mat%m,:) = b(1:qrm_mat%m,:)
  
  if(transp .eq. 'n') then
     call qrm_apply(qrm_mat_tilde, 't', b_tilde)
     call qrm_solve(qrm_mat_tilde, 'n', b_tilde, x)
  else if(transp .eq. 't') then
     call qrm_solve(qrm_mat_tilde, 't', b_tilde, x)
     call qrm_apply(qrm_mat_tilde, 'n', x)
  end if
  
  ! compute the residual
  call qrm_residual_norm(qrm_mat, r, x, rnrm)
  call qrm_vecnrm(x, size(x,1), '2', xnrm)
  call qrm_vecnrm(b, size(b,1), '2', bnrm)
  call qrm_matnrm(qrm_mat, 'f', anrm)
  call qrm_residual_orth(qrm_mat, r, onrm)   
  if(ounit.gt.0) then
     write(ounit,'(" ")')

     write(ounit,'("||A||          = ",e10.2)')anrm
     write(ounit,'(" ")')
     write(ounit,'("             ||b||         ||x||         ||r||/||A||   ||A^tr||/||r||")')
     write(ounit,'("---------------------------------------------------------------------")')
     do i=1, nrhs
        write(ounit,'("RHS ",i3,"  : ",4(e10.2,4x))')i,bnrm(i),xnrm(i),rnrm(i),onrm(i)
     end do
  end if

10 continue

  call qrm_spmat_destroy(qrm_mat, all=.true.)

  if(ounit.gt.0) write(ounit,'(" ")')
  if(ounit.gt.0) write(ounit,'("Done.")')
  if(ounit.gt.0) write(ounit,'(" ")')
  if(ounit.gt.0) write(ounit,'(" ")')

  if(ounit.gt.0) write(ounit,'("  Time to do the facto     : ",es10.3)')tf

contains

  subroutine dqrm_tikhonov(qrm_mat, gamma)
    implicit none

    type(dqrm_spmat_type)  :: qrm_mat
    real(kind(1.d0))              :: gamma

    integer                :: ounit
    integer                :: i
    integer, parameter     :: ione=1
    ! real(kind(1.d0))              :: fnorma, dnrm2

    integer, pointer, dimension(:) :: irn_tmp, jcn_tmp 
    real(kind(1.d0)), pointer, dimension(:)   :: val_tmp
    
    ! fnorma = dnrm2(qrm_mat%nz, qrm_mat%val(1), ione)
    call qrm_get('qrm_ounit', ounit)

    if(ounit.gt.0) write(ounit,'("Tikhonov regularization with gamma=",es10.2)')gamma
    
    allocate(irn_tmp(qrm_mat%nz+min(qrm_mat%m,qrm_mat%n)))
    allocate(jcn_tmp(qrm_mat%nz+min(qrm_mat%m,qrm_mat%n)))
    allocate(val_tmp(qrm_mat%nz+min(qrm_mat%m,qrm_mat%n)))

    irn_tmp(1:qrm_mat%nz) = qrm_mat%irn(1:qrm_mat%nz)
    jcn_tmp(1:qrm_mat%nz) = qrm_mat%jcn(1:qrm_mat%nz)
    val_tmp(1:qrm_mat%nz) = qrm_mat%val(1:qrm_mat%nz)
    
    ! call qrm_realloc(qrm_mat%irn, qrm_mat%nz+min(qrm_mat%m,qrm_mat%n), copy=.true.)
    ! call qrm_realloc(qrm_mat%jcn, qrm_mat%nz+min(qrm_mat%m,qrm_mat%n), copy=.true.)
    ! call qrm_realloc(qrm_mat%val, qrm_mat%nz+min(qrm_mat%m,qrm_mat%n), copy=.true.)

    deallocate(qrm_mat%irn, qrm_mat%jcn, qrm_mat%val)
    qrm_mat%irn => irn_tmp 
    qrm_mat%jcn => jcn_tmp 
    qrm_mat%val => val_tmp 
    
    do i=1, min(qrm_mat%m,qrm_mat%n)
       qrm_mat%val(qrm_mat%nz+i) = gamma
       if(qrm_mat%m .ge. qrm_mat%n) then
          qrm_mat%irn(qrm_mat%nz+i) = qrm_mat%m+i
          qrm_mat%jcn(qrm_mat%nz+i) = i
       else
          qrm_mat%jcn(qrm_mat%nz+i) = qrm_mat%n+i
          qrm_mat%irn(qrm_mat%nz+i) = i
       end if
    end do

    qrm_mat%nz = qrm_mat%nz+min(qrm_mat%m,qrm_mat%n)
    if(qrm_mat%m .ge. qrm_mat%n) then
       qrm_mat%m = qrm_mat%m + qrm_mat%n
    else
       qrm_mat%n = qrm_mat%n + qrm_mat%m
    end if
    return
  end subroutine dqrm_tikhonov

  subroutine parse_args(regularize)
    implicit none
    
    logical, intent(inout) :: regularize

    integer :: argnum, narg
    character(len=200) :: argval
    
    narg = command_argument_count()
    argnum = 1
    do while(argnum <= narg)
       call get_command_argument(argnum, argval)
       argnum = argnum + 1
       select case(argval)
       case("-r")
          regularize = .true.
       end select
    end do
    
  end subroutine parse_args
  
end program run_rd
