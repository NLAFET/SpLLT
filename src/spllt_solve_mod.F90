module spllt_solve_mod
  ! use hsl_ma87_double
! integer, external :: omp_get_thread_num, omp_get_num_threads
   use omp_lib, ONLY : omp_get_thread_num, omp_get_num_threads

   interface spllt_solve
      module procedure spllt_solve_one_double
      module procedure spllt_solve_mult_double
   end interface

contains
  
  !*************************************************
  !
  !
  ! Solve phase. simplified interface for a single rhs
  !
  subroutine spllt_solve_one_double(fkeep, options, order, x, info, job)
    use spllt_data_mod
    implicit none

    type(spllt_fkeep),    intent(in)  :: fkeep          ! Factorization data
    type(spllt_options),  intent(in)  :: options        ! User-supplied options
    integer,              intent(in)  :: order(fkeep%n) ! pivot order
    ! For details of fkeep, control, info : see derived type description        
    real(wp),           intent(inout) :: x(fkeep%n)     ! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i) is the corresponding component of the right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i) holds solution for variable i.
    type(spllt_inform),   intent(out) :: info
    integer, optional,    intent(in)  :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)

    integer :: step ! Step selector in solve phase

    ! immediate return if n = 0
    if (fkeep%n == 0) return

    ! Get the step of the solve to perform
    if(present(job)) then
      step = job
    else
      step = 0
    end if
    
    call spllt_solve_mult_double(fkeep, options, order, 1, x, info, step)

  end subroutine spllt_solve_one_double

  subroutine spllt_solve_mult_double(fkeep, options, order, nrhs, x, info, job, work)
    use spllt_data_mod
    implicit none

    type(spllt_fkeep),    intent(in)  :: fkeep          ! Factorization data
    type(spllt_options),  intent(in)  :: options        ! User-supplied options
    integer,              intent(in)  :: order(fkeep%n) ! pivot order
    integer,              intent(in)  :: nrhs           ! Number of RHS
    ! For details of fkeep, control, info : see derived type description        
    real(wp),           intent(inout) :: x(fkeep%n,nrhs)! On entry, x must
    ! be set so that if i has been used to index a variable,
    ! x(i, j) is the corresponding component of the j-th right-hand side.
    ! On exit, if i has been used to index a variable,
    ! x(i, j) holds solution for variable i of rhs-th right-hand side
    type(spllt_inform),   intent(out) :: info
    integer, optional,    intent(in)  :: job  ! used to indicate whether
    ! partial solution required
    ! job = 0 or absent: complete solve performed
    ! job = 1 : forward eliminations only (PLx = b)
    ! job = 2 : backsubs only ((PL)^Tx = b)
    real(wp), optional, target, dimension(fkeep%n, nrhs) :: work

    integer           :: j        ! Iterator
    integer           :: n        ! Order of the system
    integer           :: st       ! stat parameter
    integer           :: step     ! Step selector in solve phase
    character(len=30) :: context  ! Name of the subroutine
    real(wp), dimension(:, :), pointer :: soln
   !real(wp), dimension(:, :), allocatable :: soln

    n = fkeep%n
    context = 'spllt_solve_mult_double'

    ! immediate return if n = 0
    if (n == 0) return

    ! ma87
    ! type(ma87_control) :: ma_control
    ! type(ma87_keep) :: ma_keep
    ! type(ma87_info) :: ma_info
    
    ! Set control for HSL_MA87
    ! ma_control%nb = options%nb

    ! Set keep for HSL_MA87
    ! ma_keep%n = keep%n
    ! allocate(keep%nodes(-1:info%num_nodes+1))

    ! ma_keep = keep

    ! call spllt_solve_mult_double(1, keep%n, x, order, keep, &
    !      control, info, job)
    ! TODO Use ssids solve ?
    ! call MA87_solve(x, order, ma_keep, ma_control, ma_info, job)
    ! call MA87_solve(nrhs, n, soln, order, ma_keep, ma_control, ma_info)

    !
    ! Reorder rhs
    !
   !deallocate(soln,stat=st)

    if(.not.present(work)) then
      allocate(soln(n, nrhs), stat = st)
    else
      soln => work
    end if
    
    ! do i = 1, nrhs
    !    do j = 1, n
    !       soln((i-1)*n + order(j)) = x(j, i)
    !    end do
    ! end do

    ! Get the step of the solve to perform
    if(present(job)) then
      step = job
    else
      step = 0
    end if
    
    select case(step)
      case(0)
       !
       ! Reorder x
       !
        do j = 1, n
           soln(order(j),:) = x(j,:)
        end do

        ! Forward solve
        call solve_fwd(nrhs, soln, n, fkeep)

        ! Backward solve
        call solve_bwd(nrhs, soln, n, fkeep)

       !
       ! Reorder soln
       !
        do j = 1, n
           x(j,:) = soln(order(j),:)
        end do

      case(1)
        do j = 1, n
           soln(order(j),:) = x(j,:)
        end do

        call solve_fwd(nrhs, soln, n, fkeep)

        x = soln
      
      case(2)
        soln = x

        call solve_bwd(nrhs, soln, n, fkeep)

        do j = 1, n
           x(j,:) = soln(order(j),:)
        end do

      case default
        info%flag = SPLLT_WARNING_PARAM_VALUE
        write (0, '(a, i2, a, i4)') "Unknown requested job = ", step, &
          " returned code : ", info%flag
        return
    end select

    if(.not.present(work)) then
      deallocate(soln)
    end if

  end subroutine spllt_solve_mult_double

  !*************************************************
  !
  ! Forward solve routine
  subroutine solve_fwd(nrhs, rhs, ldr, fkeep)
    use spllt_data_mod
    use spllt_solve_task_mod
    use spllt_solve_kernels_mod
    use trace_mod
    implicit none

    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: nrhs ! Number of RHS
    integer, intent(in)           :: ldr  ! Leading dimension of RHS
    real(wp), intent(inout)       :: rhs(ldr*nrhs)

    ! real(wp) :: xlocal(keep%n)
    integer :: num_node
    integer :: node
    integer :: i, j
    integer :: sa, en, blk_sa
    integer :: numcol, numrow ! Number of column/row in node 
    integer :: nc, nr         ! Number of block-column/block-row in node
    integer :: jj, ii
    integer :: bcol, dcol, col, offset
    integer :: dblk           ! Diagonal index 
    integer :: s_nb           ! Block size in node
    integer :: m, n           ! Block dimension 
    integer :: blk            ! Block index
    integer :: st             ! Stat parameter
    integer :: fwd_update_id, fwd_block_id
    integer :: nthread, threadID
    real(wp), allocatable :: xlocal(:,:)    ! update_buffer workspace
    real(wp), allocatable :: rhs_local(:,:) ! update_buffer workspace

    nthread   = omp_get_num_threads()
    threadID  = omp_get_thread_num()
    print *, "nthreads = ", nthread

   !call trace_init(nthread)

   !call trace_create_event("fwd_update", fwd_update_id)
   !call trace_create_event("fwd_block", fwd_block_id)

    print *, "[spllt_solve_mod] solve_fwd"

    ! Allocate workspace
    allocate(xlocal(fkeep%maxmn*nrhs, nthread), &
      stat=st) !May reduce 
    allocate(rhs_local(ldr*nrhs, nthread), stat=st)

    ! initialise rhs_local
    xlocal    = zero 
    rhs_local(:,:) = zero

    num_node = fkeep%info%num_nodes
    
    do node = 1, num_node

       ! Get node info
       s_nb = fkeep%nodes(node)%nb
       sa = fkeep%nodes(node)%sa
       en = fkeep%nodes(node)%en
       numcol = en - sa + 1
       numrow = size(fkeep%nodes(node)%index)
       nc = (numcol-1) / s_nb + 1
       nr = (numrow-1) / s_nb + 1 
       
       ! Get first diag block in node
       dblk = fkeep%nodes(node)%blk_sa

       ! Loop over block columns
       do jj = 1, nc
          
          !
          ! Forward solve with block on diagoanl
          !
!         call print_darray("rhs_local before fwd block task",  &
!           ldr * nrhs, rhs_local(:, threadID + 1))

          call spllt_solve_fwd_block_task(dblk, nrhs, rhs_local, rhs, ldr, &
            xlocal, fkeep, fwd_block_id)

!         call print_darray("rhs_local after fwd block task",   &
!           ldr * nrhs, rhs_local(:, threadID + 1))
          
          do ii = jj+1, nr

             blk = dblk+ii-jj

             !
             ! Forward update with off-diagonal
             !
!            call print_darray("rhs_local before fwd update task",  &
!              ldr * nrhs, rhs_local(:, threadID + 1))

             call spllt_solve_fwd_update_task(blk, node, nrhs, rhs_local, &
               rhs, ldr, xlocal, fkeep, fwd_update_id)

!            call print_darray("rhs_local after fwd update task",   &
!              ldr * nrhs, rhs_local(:, threadID + 1))

          end do
          
          ! Update diag block in node          
          dblk = fkeep%bc(dblk)%last_blk + 1
       end do
       !$omp taskwait
              
    end do

    ! Deallocate workspace
    !$omp taskwait
    deallocate(xlocal, rhs_local)

   !call trace_log_dump_paje('trace_fwd.out')

  end subroutine solve_fwd

  subroutine solve_bwd(nrhs, rhs, ldr, fkeep)
    use spllt_data_mod
    use spllt_solve_task_mod
    use spllt_solve_kernels_mod
    use trace_mod
    implicit none

    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: nrhs   ! Number of RHS
    integer, intent(in)           :: ldr    ! Leading dimension of RHS
    real(wp), intent(inout)       :: rhs(ldr, nrhs)

    integer :: num_node
    ! Node info
    integer :: node
    integer :: sa, en
    integer :: numcol, numrow ! Number of column/row in node 
    integer :: nc, nr ! Number of block-column/block-row in node
    integer :: s_nb ! Block size in node
    integer :: jj, ii
    integer :: dblk
    ! Block info
    integer :: m, n ! Block dimension
    integer :: blk_sa 
    integer :: bcol ! Global block-column index
    integer :: dcol, col, offset
    integer :: blk ! Block index
    real(wp), dimension(:,:), allocatable :: xlocal ! update_buffer workspace
    integer :: st ! Stat parameter
    integer :: bwd_update_id, bwd_block_id
!   integer :: nthreads

!   nthreads = omp_get_num_threads()

!   call trace_init(nthreads)

!   call trace_create_event("bwd_update", bwd_update_id)
!   call trace_create_event("bwd_block", bwd_block_id)
    print *, "[spllt_solve_mod] solve_bwd"

    ! Allocate workspace
    allocate(xlocal(fkeep%maxmn*nrhs, omp_get_num_threads()), stat=st)

    num_node = fkeep%info%num_nodes
    ! print *, "num_node: ", num_node
    
    do node = num_node, 1, -1

       ! Get node info
       s_nb = fkeep%nodes(node)%nb
       sa = fkeep%nodes(node)%sa
       en = fkeep%nodes(node)%en
       numcol = en - sa + 1
       numrow = size(fkeep%nodes(node)%index)
       nc = (numcol-1) / s_nb + 1
       nr = (numrow-1) / s_nb + 1 

       ! Get first diag block in node
       dblk = fkeep%bc(fkeep%nodes(node)%blk_en)%dblk

       ! print *, "[solve_bwd] node = ", node, ", dblk = ", dblk

       ! Loop over block columns
       do jj = nc, 1, -1

          do ii = nr, jj+1, -1
             
             blk = dblk+ii-jj ! Block index

             call spllt_solve_bwd_udpate_task(blk, node, nrhs, rhs, &
               ldr, xlocal, fkeep, bwd_update_id)
             !$omp taskwait

             !
             ! Backward update with block on diagoanl
             !

             ! Establish variables describing block
 !           n        = fkeep%bc(blk)%blkn
 !           m        = fkeep%bc(blk)%blkm
 !           blk_sa       = fkeep%bc(blk)%sa
 !           bcol     = fkeep%bc(blk)%bcol
 !           ! node     = fkeep%bc(blk)%node
 !           dcol     = bcol - fkeep%bc(fkeep%nodes(node)%blk_sa)%bcol + 1
 !           col      = fkeep%nodes(node)%sa + (dcol-1)*fkeep%nodes(node)%nb

 !           offset   = col - fkeep%nodes(node)%sa + 1 ! diagonal blk
 !           offset   = offset + (blk-fkeep%bc(blk)%dblk) * fkeep%nodes(node)%nb ! this blk

 !           !$omp task depend(out: xlocal(:, omp_get_thread_num() + 1))    &
 !           !$omp depend(inout: rhs(offset : offset + m - 1, :))           &
 !           !$omp firstprivate(m, n, col, offset, node, bcol)              &
 !           !$omp firstprivate(blk_sa, nrhs, ldr)                          &
 !           !$omp shared(fkeep, xlocal, rhs)
!!           print *, "[BW UPDATE] Thread id ", omp_get_thread_num()
 !           call slv_bwd_update(m, n, col, offset, fkeep%nodes(node)%index, &
 !                fkeep%lfact(bcol)%lcol(blk_sa:blk_sa+n*m-1), n, nrhs, rhs, &
 !                rhs, ldr, xlocal(:, omp_get_thread_num() + 1))
 !           !$omp end task
             
          end do

          !
          ! Backward solve with block on diagoanl
          !
          call spllt_solve_bwd_block_task(dblk, nrhs, rhs, ldr, &
            xlocal, fkeep, bwd_block_id)
          !$omp taskwait
          
          ! Update diag block in node       
          if (jj .gt. 1) dblk = fkeep%bc(dblk-1)%dblk
       end do
       
    end do

    ! Deallocate workspace
    !$omp taskwait
    deallocate(xlocal)
!   call trace_log_dump_paje('trace.out')

  end subroutine solve_bwd
  
end module spllt_solve_mod
