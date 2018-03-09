#define L_EVALN(_N,_X,_DEP_LIST,_COMPONENT,_A,_B) L/**/_N(_X,_DEP_LIST,_COMPONENT,_A,_B)
#define L(_N,_X,_DEP_LIST,_COMPONENT,_A,_B) L_EVALN(_N,_X,_DEP_LIST,_COMPONENT,_A,_B)
                                                                                
#define L1(X,_DEP_LIST,_COMPONENT,_A,_B) X(_DEP_LIST,_COMPONENT,_A*1+_B)           
#define L2(X,_DEP_LIST,_COMPONENT,_A,_B) L1(X,_DEP_LIST,_COMPONENT,_A,_B),X(_DEP_LIST,_COMPONENT,_A*2+_B)
#define L3(X,_DEP_LIST,_COMPONENT,_A,_B) L2(X,_DEP_LIST,_COMPONENT,_A,_B),X(_DEP_LIST,_COMPONENT,_A*3+_B)
#define L4(X,_DEP_LIST,_COMPONENT,_A,_B) L3(X,_DEP_LIST,_COMPONENT,_A,_B),X(_DEP_LIST,_COMPONENT,_A*4+_B)
                                                                                   
#define EVAL_PRAGMA(_X) !$_X                                                    
#define DO_PRAGMA(_X) EVAL_PRAGMA(_X)                                           
                                                                                
#define OMP_DEPS(_VAR,_COMPONENT,_INC) _VAR(_COMPONENT(_INC))                   
                                                                                   
#define N_OMP_DEPS(_N,_DEP_LIST,_COMPONENT,_A,_B) L(_N,OMP_DEPS,_DEP_LIST,_COMPONENT,_A,_B)
                                                                                   
#define OMP_TASK(_INIT_DEP, _EXTRA) DO_PRAGMA(omp task depend(out: _INIT_DEP) depend(in: _EXTRA))
#define OMP_CONT_TASK(_EXTRA) DO_PRAGMA(omp depend(in: _EXTRA))                    
                                                                                   
#define VARIABLE_OMP_DEP(_N,_DEP_LIST,_COMPONENT,_A,_B)\
OMP_TASK(OMP_DEPS(_DEP_LIST,_COMPONENT,_A*1+_B),N_OMP_DEPS(_N,_DEP_LIST,_COMPONENT,_A,_B))
                                                                                   
#define VAR_OMP_DEP_CONTD(_N,_DEP_LIST,_COMPONENT,_A,_B) OMP_CONT_TASK(N_OMP_DEPS(_N,_DEP_LIST,_COMPONENT,_A,_B))
                                                                                   
#define FPRIV_VAR_DECL(_SA) !$omp firstprivate(m, n, nrhs, col, ldr, _SA, offset, threadID)             
#define FPRIV_ADD_VAR_DECL(_VAR) !$omp firstprivate(_VAR)
#define FPRIV_PTR_DECL !$omp firstprivate(p_upd, p_rhs, p_lcol, p_index, p_xlocal)
#define FPRIV_VAR_DEP_DECL !$omp firstprivate(p_bc, p_dep)
#define FPRIV_VAR_LOOP_DECL !$omp firstprivate(chunk, ndep_lvl, alpha, beta)
#define RELEASE_BLK(_BLK) !$omp firstprivate(_BLK) depend(inout: p_bc(_BLK))

module spllt_solve_task_mod
contains
  
  !*************************************************  
  !
  ! Forward solve with block on diagoanl
  !
  subroutine spllt_solve_fwd_block_task(dblk, nrhs, upd, rhs, ldr, xlocal, &
      fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    use spllt_solve_dep_mod
    implicit none
    
    integer, intent(in)                       :: dblk !Index of diagonal block
    integer, intent(in)                       :: nrhs !Number of RHS
    integer, intent(in)                       :: ldr  !Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:, :)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target,  intent(inout)          :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(spllt_omp_scheduler), intent(inout)  :: scheduler
    
    ! Node info
    integer :: sa
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa
    integer :: bcol, dcol, col
    integer :: offset
    integer :: node
    integer :: dep
    integer :: i, j, r, chunkth
    integer :: threadID, nthread
    integer :: ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    type(spllt_block), pointer  :: p_blk_dep_update

    real(wp), pointer :: p_upd(:,:)
    real(wp), pointer :: p_xlocal(:,:)
    real(wp), pointer :: p_rhs(:)
    integer,  pointer :: p_dep(:)
    integer,  pointer :: p_child_blk_index(:)
    integer :: new_method
    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted

    type(spllt_block), pointer :: p_bc(:)
        
    nthread   = omp_get_num_threads()

    ! Get block info
    node      = fkeep%bc(dblk)%node
    m         = fkeep%bc(dblk)%blkm
    n         = fkeep%bc(dblk)%blkn
    sa        = fkeep%bc(dblk)%sa
    bcol      = fkeep%bc(dblk)%bcol ! Current block column
    dcol      = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col       = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
    p_bc      => fkeep%bc
    
    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs
    
    new_method = 1
!   call fwd_update_dependency(fkeep, dblk, p_dep)
!   p_dep => fkeep%bc(dblk)%fwd_update_dep
    p_dep => fkeep%bc(dblk)%fwd_dep

    ndep = size(p_dep)

    if(new_method .eq. 1) then

!     threadID  = omp_get_thread_num()
      chunk     = 0
      ndep_lvl  = 0

      if(ndep .eq. 0) then
        !$omp task                                &
        include 'spllt_solve_fwd_block_omp_decl.F90'

        include 'spllt_solve_fwd_block_worker.F90'

        !$omp end task
      else
        chunk = 10 ! Do not use chunk = 1 ; a non-sence
        lvl   = 1
        alpha = 1
        ndep_lvl = ndep ! #dep local to the lvl
        all_task_submitted = .false.

        do while(.not. all_task_submitted)
        
          nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)

          beta = 1 - alpha

          do chunkth = 1, nchunk
            chunk_size = merge(ndep_lvl - (chunkth - 1) * chunk, chunk, &
              chunkth * chunk .gt. ndep_lvl)

            select case(chunk_size)

              case(0)
                print *, "No dep ??"

              !
              !This file contains the remaining cases that are generated through a script
              !
              include 'spllt_fwd_block_cases.F90'

            end select

            beta = beta + chunk_size
            if(ndep_lvl .le. chunk) then
              all_task_submitted = .true.
            else
              scheduler%task_info(scheduler%workerID)%nfake_task_insert =   &
                scheduler%task_info(scheduler%workerID)%nfake_task_insert   &
                + 1
            end if

          end do
          if(ndep_lvl .gt. chunk) then
            ndep_lvl = nchunk
            lvl = lvl + 1
            alpha = alpha * chunk
          end if
        end do
      end if

    else
      dep  = ndep


   !do dep = 1, ndep - 1

!  !  p_lcol_update     => fkeep%lfact(fkeep%bc(p_dep(dep))%bcol)%lcol
!  !  p_blk_dep_update  => fkeep%bc(p_dep(dep))

 ! !  !$omp task firstprivate(p_lcol_update, p_blk_dep_update)      &
 ! !  !$omp firstprivate(p_dep)                              &
 ! !  !$omp firstprivate(dblk, dep)                                 &
 ! !  !$omp private(threadID)                                       &
 ! !  !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
 ! !  !$omp depend(inout: p_dep(1))

   !  !$omp task firstprivate(dep, p_lcol_sa, dep_lcol_sa)      &
   !  !$omp firstprivate(p_bc, p_dep)                    &
   !  !$omp private(threadID)                                   &
   ! !!$omp depend(in: p_lcol_sa(dep)%p(dep_lcol_sa(dep)))      &
   !  !$omp depend(in: p_bc(p_dep(dep)))                 &
   !  !$omp depend(inout: p_dep(1))
   !  threadID = omp_get_thread_num()
!  !  print '(a, i3, a, i3)', "SLV Fake Task dep ", dblk, " [in : ", &
!  !    p_dep(dep)

   !  !$omp end task

   !end do

      p_lcol_update     => fkeep%lfact(fkeep%bc(p_dep(dep))%bcol)%lcol
      p_blk_dep_update  => fkeep%bc(p_dep(dep))
    !$omp task firstprivate(m, n, col, offset)                      &
    !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)              &
!   !$omp firstprivate(dep_update, dblk)                            &
    !$omp firstprivate(dblk)                                        &
    !$omp firstprivate(nthread)                                     &
    !$omp firstprivate(p_blk_dep_update, p_lcol_update)             &
    !$omp firstprivate(p_dep, dep, p_bc)                     &
    !$omp firstprivate(p_upd, p_rhs, p_xlocal)                      &
    !$omp private(threadID, r, j, i)                                &
   !!$omp depend(in: p_lcol_update(p_blk_dep_update%sa))            &
    !$omp depend(in: p_bc(p_dep(dep)))                       &
    !$omp depend(in: p_dep(1))                               &
   !!$omp depend(inout: p_lcol(sa))                                  
    !$omp depend(in   : p_bc(dblk))                                 &
    !$omp depend(  out: p_bc(dblk))


    call trace_event_start(trace_id, omp_get_thread_num())

    threadID  = omp_get_thread_num()
!   print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
!   print *, p_dep

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          p_rhs(i)    = p_rhs(i) + p_upd(i, j)
          p_upd(i,j)  = zero ! Reset in case of bwd solve
        end do
      end do
    end do


    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa:sa+n*n-1),    &
         'Transpose    ', 'Non-unit', nrhs, p_rhs, ldr)
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m .gt. 0) then
       sa = sa + n * n
       call slv_fwd_update(m, n, col, offset, p_index,                &
            p_lcol(sa : sa + n * m - 1), n, nrhs,                     &
            p_upd(:, threadID + 1), ldr, p_rhs,                       &
            ldr, p_xlocal(:, threadID + 1))
    endif

!   deallocate(p_dep)
    call trace_event_stop (trace_id, omp_get_thread_num())

    !$omp end task
    end if

    scheduler%task_info(scheduler%workerID)%max_dep = merge(ndep, &
      scheduler%task_info(scheduler%workerID)%max_dep,            &
      ndep .gt. scheduler%task_info(scheduler%workerID)%max_dep)
    scheduler%task_info(scheduler%workerID)%nblk_require_fake_task =    &
      scheduler%task_info(scheduler%workerID)%nblk_require_fake_task    &
      + merge(1, 0, ndep .gt. 1)
    if(new_method .eq. 0) then
      scheduler%task_info(scheduler%workerID)%nfake_task_insert =   &
        scheduler%task_info(scheduler%workerID)%nfake_task_insert   &
        + merge(ndep - 1, 0, ndep - 1 .gt. 0)
    end if
    scheduler%task_info(scheduler%workerID)%ntask_insert =   &
      scheduler%task_info(scheduler%workerID)%ntask_insert + 1
  end subroutine spllt_solve_fwd_block_task


  subroutine spllt_solve_fwd_update_task(blk, node, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num
    use trace_mod
    use spllt_solve_dep_mod
    implicit none
    
    integer, intent(in)                       :: blk  ! Index of block
    integer, intent(in)                       :: node
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:,:)        
    real(wp), target, intent(in)              :: rhs(ldr*nrhs)
    real(wp), target, intent(out)             :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(spllt_omp_scheduler), intent(inout)  :: scheduler

    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID
    integer                     :: dep, ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    real(wp), pointer           :: p_lcol_solve(:)
    integer, pointer            :: p_dep(:)
    integer                     :: blk_dep_solve
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)
    type(spllt_block), pointer  :: p_bc(:)
    integer :: j
    integer :: new_method
    integer :: chunk, chunk_size, ndep_lvl, lvl, nchunk, beta, alpha
    logical :: all_task_submitted

    ! Establish variables describing block
    n         = fkeep%bc(blk)%blkn
    m         = fkeep%bc(blk)%blkm
    blk_sa    = fkeep%bc(blk)%sa
    bcol      = fkeep%bc(blk)%bcol
    dcol      = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col       = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset    = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
    offset    = offset + (blk-fkeep%bc(blk)%dblk) &
      * fkeep%nodes(node)%nb ! this blk
    p_index   => fkeep%nodes(node)%index
    p_lcol    => fkeep%lfact(bcol)%lcol
    p_bc      => fkeep%bc

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

!   call fwd_update_dependency(fkeep, blk, p_dep)
!   blk_dep_solve  = fwd_solve_dependency(fkeep, blk)
!   p_dep  => fkeep%bc(blk)%fwd_update_dep
    blk_dep_solve = fkeep%bc(blk)%fwd_solve_dep

    p_dep  => fkeep%bc(blk)%fwd_dep
    ndep = size(p_dep)

    new_method = 1

    if(new_method .eq. 1) then

!     threadID  = omp_get_thread_num()
      chunk     = 0
      ndep_lvl  = 0

      if(ndep .eq. 0) then
        !$omp task                                &
        include 'spllt_solve_fwd_update_omp_decl.F90'

        include 'spllt_solve_fwd_update_worker.F90'
   !    threadID  = omp_get_thread_num()
!  !    print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : ]"
   !    call slv_fwd_update(m, n, col, offset, p_index,         &
   !      p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
   !      p_upd(:, threadID + 1), ldr, p_rhs,                   &
   !      ldr, p_xlocal(:, threadID + 1))

        !$omp end task
      else

        chunk = 10 ! Do not use chunk = 1 ; a non-sence
        lvl   = 1
        alpha = 1
        ndep_lvl = ndep ! #dep local to the lvl
        all_task_submitted = .false.
  !     print *, "================"
  !     print *, "Chunk ", chunk

        do while(.not. all_task_submitted)
        
          nchunk = ceiling( (ndep_lvl  + 0.0 ) / chunk)
  !       print *, "LVL : ", lvl, " #dep ", ndep_lvl, " => #chunk ", nchunk

          beta = 1 - alpha

          do j = 1, nchunk
            chunk_size = merge(ndep_lvl - (j - 1) * chunk, chunk, &
              j * chunk .gt. ndep_lvl)
  !         print *, "chunk id ", j, " with a chunk_size : ", chunk_size

  !         print *, "Select case of lvl ", lvl, " with alpha ", alpha, " and beta ", beta
            select case(chunk_size)

              case(0)
                print *, "No dep ?? "

              !
              !This file contains the remaining cases that are generated through a script
              !
              include 'spllt_fwd_update_cases.F90'
          !   case(1)
          !   VARIABLE_OMP_DEP(1,p_bc,p_dep,alpha,beta) &
          !   include 'spllt_solve_fwd_update_omp_decl.F90'
          !  !FPRIV_VAR_DECL(blk_sa)                    &
          !  !FPRIV_PTR_DECL                            &
          !  !FPRIV_VAR_DEP_DECL                        &
          !  !FPRIV_VAR_LOOP_DECL                       &
          !  !RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(2)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta) &
          !   FPRIV_VAR_DECL(blk_sa)                    &
          !   FPRIV_PTR_DECL                            &
          !   FPRIV_VAR_DEP_DECL                        &
          !   FPRIV_VAR_LOOP_DECL                       &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(3)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta)   &
          !   VAR_OMP_DEP_CONTD(1,p_bc,p_dep,alpha,beta+2)&
          !   FPRIV_VAR_DECL(blk_sa)                      &
          !   FPRIV_PTR_DECL                              &
          !   FPRIV_VAR_DEP_DECL                          &
          !   FPRIV_VAR_LOOP_DECL                         &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(4)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta)   &
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+2)&
          !   FPRIV_VAR_DECL(blk_sa)                      &
          !   FPRIV_PTR_DECL                              &
          !   FPRIV_VAR_DEP_DECL                          &
          !   FPRIV_VAR_LOOP_DECL                         &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(5)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta)   &
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+2)&
          !   VAR_OMP_DEP_CONTD(1,p_bc,p_dep,alpha,beta+4)&
          !   FPRIV_VAR_DECL(blk_sa)                      &
          !   FPRIV_PTR_DECL                              &
          !   FPRIV_VAR_DEP_DECL                          &
          !   FPRIV_VAR_LOOP_DECL                         &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(6)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta)   &
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+2)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+4)&
          !   FPRIV_VAR_DECL(blk_sa)                      &
          !   FPRIV_PTR_DECL                              &
          !   FPRIV_VAR_DEP_DECL                          &
          !   FPRIV_VAR_LOOP_DECL                         &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(7)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta)   &
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+2)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+4)&
          !   VAR_OMP_DEP_CONTD(1,p_bc,p_dep,alpha,beta+6)&
          !   FPRIV_VAR_DECL(blk_sa)                      &
          !   FPRIV_PTR_DECL                              &
          !   FPRIV_VAR_DEP_DECL                          &
          !   FPRIV_VAR_LOOP_DECL                         &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(8)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta)   &
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+2)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+4)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+6)&
          !   FPRIV_VAR_DECL(blk_sa)                      &
          !   FPRIV_PTR_DECL                              &
          !   FPRIV_VAR_DEP_DECL                          &
          !   FPRIV_VAR_LOOP_DECL                         &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(9)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta)   &
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+2)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+4)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+6)&
          !   VAR_OMP_DEP_CONTD(1,p_bc,p_dep,alpha,beta+8)&
          !   FPRIV_VAR_DECL(blk_sa)                      &
          !   FPRIV_PTR_DECL                              &
          !   FPRIV_VAR_DEP_DECL                          &
          !   FPRIV_VAR_LOOP_DECL                         &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

          !   case(10)
          !   VARIABLE_OMP_DEP(2,p_bc,p_dep,alpha,beta)   &
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+2)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+4)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+6)&
          !   VAR_OMP_DEP_CONTD(2,p_bc,p_dep,alpha,beta+8)&
          !   FPRIV_VAR_DECL(blk_sa)                      &
          !   FPRIV_PTR_DECL                              &
          !   FPRIV_VAR_DEP_DECL                          &
          !   FPRIV_VAR_LOOP_DECL                         &
          !   RELEASE_BLK(blk)

          !   threadID  = omp_get_thread_num()
          !   if(ndep_lvl .le. chunk) then
!         !     print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!         !     print *, p_dep
          !     call slv_fwd_update(m, n, col, offset, p_index,         &
          !       p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
          !       p_upd(:, threadID + 1), ldr, p_rhs,                   &
          !       ldr, p_xlocal(:, threadID + 1))
          !   end if

          !   !$omp end task

            end select
            beta = beta + chunk_size
            if(ndep_lvl .le. chunk) then
              all_task_submitted = .true.
            else
              scheduler%task_info(scheduler%workerID)%nfake_task_insert =   &
                scheduler%task_info(scheduler%workerID)%nfake_task_insert   &
                + 1
            end if

          end do
          if(ndep_lvl .gt. chunk) then
            ndep_lvl = nchunk
            lvl = lvl + 1
            alpha = alpha * chunk
          end if
        end do
      end if
!     !$omp taskwait
    else
    
    dep = ndep

  ! do dep = 1, ndep - 1

  !!  p_lcol_update => fkeep%lfact(fkeep%bc(p_dep(dep))%bcol)%lcol
  !!  p_blk_dep_update  => fkeep%bc(p_dep(dep))

  !!  !$omp task firstprivate(p_lcol_update, p_blk_dep_update)      &
  !!  !$omp firstprivate(p_dep)                              &
! !!  !$omp firstprivate(blk, dep)                                  &
  !!  !$omp private(threadID)                                       &
  !!  !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
  !!  !$omp depend(inout: p_dep(1))
  !   
  !   !$omp task firstprivate(p_lcol_update, p_blk_dep_update)      &
  !   !$omp firstprivate(p_dep, p_bc, dep)                   &
! !   !$omp firstprivate(blk, dep)                                  &
  !   !$omp private(threadID)                                       &
  !   !$omp depend(in: p_bc(p_dep(dep)))                     &
  !   !$omp depend(inout: p_dep(1))

  !   threadID = omp_get_thread_num()
! !   print '(a, i3, a, i3)', "UPD Fake Task dep of ", blk, &
! !     " [in : ", p_dep(dep)

  !   !$omp end task

  ! end do

    p_lcol_update => fkeep%lfact(fkeep%bc(p_dep(dep))%bcol)%lcol
    p_lcol_solve  => fkeep%lfact(fkeep%bc(blk_dep_solve )%bcol)%lcol

    p_blk_dep_update  => fkeep%bc(p_dep(dep))
    p_blk_dep_solve   => fkeep%bc(blk_dep_solve)

    !$omp task firstprivate(m, n, col, offset)                    &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(p_blk_dep_update, p_blk_dep_solve)         &
    !$omp firstprivate(p_lcol_update, p_lcol_solve)               &
    !$omp firstprivate(blk)                                       &
    !$omp firstprivate(p_bc, p_dep, dep)                   &
    !$omp private(threadID)                                       &
   !!$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
   !!$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(in: p_bc(p_dep(dep)))                     &
    !$omp depend(in: p_bc(blk_dep_solve))                         &
    !$omp depend(in: p_dep(1))                             &
   !!$omp depend(inout: p_lcol(blk_sa))                            
    !$omp depend(inout: p_bc(blk))

    call trace_event_start(trace_id, omp_get_thread_num())
    
    threadID  = omp_get_thread_num()
    print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
    print *, p_dep
    print *, blk_dep_solve

    call slv_fwd_update(m, n, col, offset, p_index,         &
      p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,         &
      p_upd(:, threadID + 1), ldr, p_rhs,                   &
      ldr, p_xlocal(:, threadID + 1))

    deallocate(p_dep)
    call trace_event_stop (trace_id, omp_get_thread_num())

    !$omp end task

  end if
    scheduler%task_info(scheduler%workerID)%max_dep = merge(ndep, &
      scheduler%task_info(scheduler%workerID)%max_dep,            &
      ndep .gt. scheduler%task_info(scheduler%workerID)%max_dep)
    scheduler%task_info(scheduler%workerID)%nblk_require_fake_task =    &
      scheduler%task_info(scheduler%workerID)%nblk_require_fake_task +  &
      merge(1, 0, ndep .gt. 1)
    if(new_method .eq. 0) then
      scheduler%task_info(scheduler%workerID)%nfake_task_insert =   &
        scheduler%task_info(scheduler%workerID)%nfake_task_insert   &
        + merge(ndep - 1, 0, ndep - 1 .gt. 0)
    end if
    scheduler%task_info(scheduler%workerID)%ntask_insert =   &
      scheduler%task_info(scheduler%workerID)%ntask_insert + 1
  end subroutine spllt_solve_fwd_update_task

  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         
  subroutine spllt_solve_bwd_block_task(dblk, nrhs, upd, rhs, ldr, xlocal, &
      fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    use spllt_solve_dep_mod
    implicit none

    integer, intent(in)                       :: dblk ! Index of diagonal block
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:, :)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:, :)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(spllt_omp_scheduler), intent(inout)  :: scheduler
    
    ! Node info
    integer                     :: sa
    ! Block info
    integer                     :: m, n ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: node
    integer                     :: dep
    integer                     :: i, j, r
    integer                     :: threadID, nthread
    integer                     :: ndep
    integer                     :: blk_dep_update
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    real(wp), pointer           :: p_lcol_solve(:)
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve

    real(wp), pointer :: p_upd(:,:)
    real(wp), pointer :: p_xlocal(:,:)
    real(wp), pointer :: p_rhs(:)
    integer,  pointer :: p_dep_solve(:)

    nthread = omp_get_num_threads()

    ! Get block info
    node    = fkeep%bc(dblk)%node
    n       = fkeep%bc(dblk)%blkn
    m       = fkeep%bc(dblk)%blkm
    sa      = fkeep%bc(dblk)%sa
    bcol    = fkeep%bc(dblk)%bcol ! Current block column
    col     = calc_col(fkeep%nodes(node), fkeep%bc(dblk)) ! current bcol
    col     = fkeep%nodes(node)%sa + (col-1)*fkeep%nodes(node)%nb
    offset  = col - fkeep%nodes(node)%sa + 1
    p_index => fkeep%nodes(node)%index
    p_lcol  => fkeep%lfact(bcol)%lcol

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

!   call bwd_solve_dependency(fkeep, dblk, p_dep_solve)
    p_dep_solve => fkeep%bc(dblk)%bwd_solve_dep
    ndep = size(p_dep_solve)

    do dep = 1, ndep - 1

      p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol
      p_blk_dep_solve => fkeep%bc(p_dep_solve(dep))

      !$omp task firstprivate(p_lcol_solve, p_blk_dep_solve)        &
      !$omp firstprivate(p_dep_solve)                               &
      !$omp firstprivate(dblk, dep)                                 &
      !$omp private(threadID)                                       &
      !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
      !$omp depend(inout: p_dep_solve(1))

!     threadID = omp_get_thread_num()
!     print '(a, i3, a, i3)', "SLV Fake Task dep ", dblk, " [in : ",&
!       p_dep_solve(dep)

      !$omp end task

    end do

!   blk_dep_update  = bwd_update_dependency(fkeep, dblk)
    blk_dep_update  = fkeep%bc(dblk)%bwd_update_dep
    p_lcol_update   => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
    p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol

    p_blk_dep_update  => fkeep%bc(blk_dep_update)   
    p_blk_dep_solve   => fkeep%bc(p_dep_solve(dep))

    !$omp task firstprivate(m, n, col, offset, node, bcol)        &
    !$omp firstprivate(sa, nrhs, ldr, p_index, p_lcol)            &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(p_dep_solve, dblk)                         &
    !$omp firstprivate(p_lcol_update, p_blk_dep_update)           &
    !$omp firstprivate(p_lcol_solve, p_blk_dep_solve)             &
    !$omp firstprivate(blk_dep_update)                            &
    !$omp firstprivate(nthread)                                   &
    !$omp private(threadID, r, j, i)                              &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
    !$omp depend(in: p_dep_solve(1))                              &
    !$omp depend(inout: p_lcol(sa))                                  

    call trace_event_start(trace_id, omp_get_thread_num())

    threadID = omp_get_thread_num()
!   print '(a, i3, a, i3)', "SLV      Task dep of ", dblk, " [in : "
!   print *, p_dep_solve
!   print *, blk_dep_update

    ! Perform retangular update from diagonal block
    if(m .gt. n) then
       call slv_bwd_update(m - n, n, col, offset + n, p_index,      &
            p_lcol(sa + n * n : sa + n * m - 1), n, nrhs, p_rhs,    &
            p_upd(:, threadID + 1), ldr, p_xlocal(:, threadID + 1))
    endif

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          p_rhs(i)    = p_rhs(i) + p_upd(i, j)
        end do
      end do
    end do

    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa : sa + n * n - 1), &
         'Non-Transpose', 'Non-unit', nrhs, p_rhs, ldr)
    
    call trace_event_stop (trace_id, omp_get_thread_num())
    !$omp end task

    scheduler%task_info(scheduler%workerID)%max_dep = merge(ndep, &
      scheduler%task_info(scheduler%workerID)%max_dep,            &
      ndep .gt. scheduler%task_info(scheduler%workerID)%max_dep)
    scheduler%task_info(scheduler%workerID)%nblk_require_fake_task =    &
      scheduler%task_info(scheduler%workerID)%nblk_require_fake_task    &
      + merge(1, 0, ndep .gt. 1)
    scheduler%task_info(scheduler%workerID)%nfake_task_insert =   &
      scheduler%task_info(scheduler%workerID)%nfake_task_insert   &
      + merge(ndep - 1, 0, ndep - 1 .gt. 0)
    scheduler%task_info(scheduler%workerID)%ntask_insert =   &
      scheduler%task_info(scheduler%workerID)%ntask_insert + 1
  end subroutine spllt_solve_bwd_block_task

  !*************************************************  
  !
  ! Backward solve with block on diagoanl
  !         
  subroutine spllt_solve_bwd_udpate_task(blk, node, nrhs, upd, rhs, ldr, &
      xlocal, fkeep, trace_id, scheduler)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads
    use trace_mod
    use spllt_solve_dep_mod
    implicit none

    integer, intent(in)                       :: blk  ! Index of block 
    integer, intent(in)                       :: node 
    integer, intent(in)                       :: nrhs ! Number of RHS
    integer, intent(in)                       :: ldr  ! Leading dimension of RHS
    real(wp), target, intent(inout)           :: upd(:,:)
    real(wp), target, intent(inout)           :: rhs(ldr * nrhs)
    real(wp), target, intent(inout)           :: xlocal(:,:)
    type(spllt_fkeep), target, intent(in)     :: fkeep
    integer, intent(in)                       :: trace_id
    type(spllt_omp_scheduler), intent(inout)  :: scheduler
    
    ! Block info
    integer                     :: m, n         ! Block dimension
    integer                     :: blk_sa
    integer                     :: bcol, dcol, col
    integer                     :: offset
    integer                     :: threadID
    integer                     :: dep, ndep
    integer, pointer            :: p_index(:)
    real(wp), pointer           :: p_lcol(:)
    real(wp), pointer           :: p_lcol_update(:)
    real(wp), pointer           :: p_lcol_solve(:)
    integer, pointer            :: p_dep_solve(:)
    integer                     :: blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_update
    type(spllt_block), pointer  :: p_blk_dep_solve
    real(wp)         , pointer  :: p_upd(:,:)
    real(wp)         , pointer  :: p_xlocal(:,:)
    real(wp)         , pointer  :: p_rhs(:)

    ! Establish variables describing block
    n       = fkeep%bc(blk)%blkn
    m       = fkeep%bc(blk)%blkm
    blk_sa  = fkeep%bc(blk)%sa
    bcol    = fkeep%bc(blk)%bcol
    dcol    = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
    col     = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb
    offset  = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
    offset  = offset + (blk-fkeep%bc(blk)%dblk) * fkeep%nodes(node)%nb !this blk
    p_index => fkeep%nodes(node)%index
    p_lcol  => fkeep%lfact(bcol)%lcol

    p_upd     => upd
    p_xlocal  => xlocal
    p_rhs     => rhs

!   blk_dep_update = bwd_update_dependency(fkeep, blk)
!   call bwd_solve_dependency(fkeep, blk, p_dep_solve)
    blk_dep_update = fkeep%bc(blk)%bwd_update_dep
    p_dep_solve => fkeep%bc(blk)%bwd_solve_dep
    ndep = size(p_dep_solve)

    do dep = 1, ndep - 1

      p_lcol_solve    => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol
      p_blk_dep_solve => fkeep%bc(p_dep_solve(dep))

      !$omp task firstprivate(p_lcol_solve, p_blk_dep_solve)        &
      !$omp firstprivate(p_dep_solve)                               &
      !$omp firstprivate(blk, dep)                                  &
      !$omp private(threadID)                                       &
      !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
      !$omp depend(inout: p_dep_solve(1))

!     threadID = omp_get_thread_num()
!     print '(a, i3, a, i3)', "UPD Fake Task dep ", blk, " [in : ", &
!       p_dep_solve(dep)

      !$omp end task

    end do

    p_lcol_update => fkeep%lfact(fkeep%bc(blk_dep_update)%bcol)%lcol
    p_lcol_solve  => fkeep%lfact(fkeep%bc(p_dep_solve(dep))%bcol)%lcol

    p_blk_dep_update  => fkeep%bc(blk_dep_update)
    p_blk_dep_solve   => fkeep%bc(p_dep_solve(dep))

    !$omp task firstprivate(m, n, col, offset, node, bcol)        &
    !$omp firstprivate(blk_sa, nrhs, ldr, p_index, p_lcol)        &
    !$omp firstprivate(p_xlocal, p_rhs, p_upd)                    &
    !$omp firstprivate(p_blk_dep_update, p_blk_dep_solve)         &
    !$omp firstprivate(p_lcol_update, p_lcol_solve)               &
    !$omp firstprivate(blk)                                       &
    !$omp private(threadID)                                       &
    !$omp depend(in: p_lcol_solve(p_blk_dep_solve%sa))            &
    !$omp depend(in: p_lcol_update(p_blk_dep_update%sa))          &
    !$omp depend(in: p_dep_solve(1))                              &
    !$omp depend(inout: p_lcol(blk_sa))                            

    call trace_event_start(trace_id, omp_get_thread_num())

    threadID  = omp_get_thread_num()
!   print '(a, i3, a)', "UPD      Task dep of ", blk, " [in : "
!   print *, blk_dep_update
!   print *, p_dep_solve

    call  slv_bwd_update(m, n, col, offset, p_index,      &
          p_lcol(blk_sa : blk_sa + n * m - 1), n, nrhs,   &
          p_rhs, p_upd(:, threadID + 1), ldr,             &
          xlocal(:, omp_get_thread_num() + 1))

    deallocate(p_dep_solve)
    call trace_event_stop (trace_id, omp_get_thread_num())

    !$omp end task

    scheduler%task_info(scheduler%workerID)%max_dep = merge(ndep, &
      scheduler%task_info(scheduler%workerID)%max_dep,            &
      ndep .gt. scheduler%task_info(scheduler%workerID)%max_dep)
    scheduler%task_info(scheduler%workerID)%nblk_require_fake_task =    &
      scheduler%task_info(scheduler%workerID)%nblk_require_fake_task    &
      + merge(1, 0, ndep .gt. 1)
    scheduler%task_info(scheduler%workerID)%nfake_task_insert =   &
      scheduler%task_info(scheduler%workerID)%nfake_task_insert   &
      + merge(ndep - 1, 0, ndep - 1 .gt. 0)
    scheduler%task_info(scheduler%workerID)%ntask_insert =   &
      scheduler%task_info(scheduler%workerID)%ntask_insert + 1
  end subroutine spllt_solve_bwd_udpate_task

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

  subroutine spllt_solve_fwd_block_worker(m, n, nrhs, nthread, col, ldr, sa, &
      offset, threadID, p_upd, p_rhs, p_lcol, p_index, p_xlocal)
    use spllt_data_mod
    use spllt_solve_kernels_mod
    integer, intent(inout)        :: m
    integer, intent(in)           :: n
    integer, intent(in)           :: nrhs
    integer, intent(in)           :: nthread
    integer, intent(in)           :: col
    integer, intent(in)           :: ldr
    integer, intent(inout)        :: sa
    integer, intent(inout)        :: offset
    integer, intent(in)           :: threadID
    real(wp), pointer, intent(in) :: p_upd(:,:)
    real(wp), pointer, intent(in) :: p_rhs(:)
    real(wp), pointer, intent(in) :: p_lcol(:)
    integer, pointer, intent(in)  :: p_index(:)
    real(wp), pointer, intent(in) :: p_xlocal(:,:)


    integer :: i,j,k

    ! Sum contributions to rhs
    do r = 0, nrhs-1
      do j = 1, nthread
        do i = col + r*ldr, col+n-1 + r*ldr
          p_rhs(i)    = p_rhs(i) + p_upd(i, j)
          p_upd(i,j)  = zero ! Reset in case of bwd solve
        end do
      end do
    end do


    ! Perform triangular solve
    call slv_solve(n, n, col, p_lcol(sa:sa+n*n-1),    &
      'Transpose    ', 'Non-unit', nrhs, p_rhs, ldr)
    offset = offset + n

    ! Deal with any left over trapezoidal part of diagonal block
    m = m - n
    if(m .gt. 0) then
      sa = sa + n * n
      call slv_fwd_update(m, n, col, offset, p_index,             &
        p_lcol(sa : sa + n * m - 1), n, nrhs,                     &
        p_upd(:, threadID + 1), ldr, p_rhs,                       &
        ldr, p_xlocal(:, threadID + 1))
    endif

  end subroutine spllt_solve_fwd_block_worker

end module spllt_solve_task_mod
