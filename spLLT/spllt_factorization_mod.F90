module spllt_factorization_mod
  use spllt_mod
  implicit none
  
contains

  ! init node
  ! subroutine spllt_factorize_block_task(node)
  !   implicit none

  !   type(node_type), intent(inout) :: node ! node in the atree

  !   return
  ! end subroutine spllt_factorize_block_task

  ! factorize block 
  ! _potrf
  subroutine spllt_factorize_block_task(node, bc, lfact, prio)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none
    
    type(spllt_node_type), intent(in)          :: node
    type(spllt_bc_type), target, intent(inout) :: bc ! block to be factorized    
#if defined(SPLLT_USE_OMP)
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
#else
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    integer, optional :: prio

    integer :: m, n, bcol, sa
    integer(long) :: id
    type(block_type), pointer :: blk ! block to be factorized
    integer :: p
#if defined(SPLLT_USE_STARPU)
    type(c_ptr) :: node_hdl
#elif defined(SPLLT_USE_OMP)
    real(wp), dimension(:), pointer :: lcol
    real(wp), dimension(:), pointer :: bc_c
#if defined(SPLLT_OMP_TRACE)
    integer :: th_id
#endif
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)
    node_hdl = c_null_ptr
    if (node%node%blk_sa .eq. bc%blk%id) then
       node_hdl = node%hdl
       ! write(*,*)'Test'
    end if
    call spllt_starpu_insert_factorize_block_c(bc%hdl, node_hdl, p)

    ! call spllt_starpu_insert_factorize_block(bc, p)
#elif defined(SPLLT_USE_OMP)

    bc_c => bc%c
    blk => bc%blk
    
!$omp task private(m, n, sa, bcol, bc, lcol) &
!$omp    & firstprivate(blk, bc_c) &
!$omp    & shared(lfact) &
!$omp    & depend(inout:bc_c(1:1))

! !$omp    & depend(inout:lcol(sa:n*m-1))

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(fac_blk_id, th_id)
#endif

    bcol = blk%bcol
    lcol => lfact(bcol)%lcol

    m    = blk%blkm
    n    = blk%blkn
    sa   = blk%sa

    ! write(*,*)"bcol: ", bcol
    call spllt_factor_diag_block(m, n, lcol(sa:sa+n*m-1))

#if defined(SPLLT_OMP_TRACE)
    ! write(*,*)"fac_blk_id: ", fac_blk_id
    call trace_event_stop(fac_blk_id, th_id)
#endif
!$omp end task

#else    

    blk => bc%blk

    m    = blk%blkm
    n    = blk%blkn
    id   = blk%id 
    bcol = blk%bcol

    sa   = blk%sa

    ! factorize_block
    ! call factor_diag_block(n, m, id, &
    !      & lfact(bcol)%lcol(sa:sa+n*m-1),   &
    !      & control, flag, detlog(0))

    call spllt_factor_diag_block(m, n, lfact(bcol)%lcol(sa:sa+n*m-1))
#endif
    
    return
  end subroutine spllt_factorize_block_task

  ! _trsm
  subroutine spllt_solve_block_task(bc_kk, bc_ik, lfact, prio)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none

    type(spllt_bc_type), intent(inout) :: bc_kk, bc_ik ! block to be factorized    
#if defined(SPLLT_USE_OMP)
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
#else
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    integer, optional :: prio 
    
    integer :: m, n, bcol, sa
    integer :: d_m, d_n, d_sa
    integer(long) :: id, d_id
    integer :: p
#if defined(SPLLT_USE_OMP)
    real(wp), dimension(:), pointer :: lcol
    real(wp), dimension(:), pointer :: bc_kk_c, bc_ik_c
#if defined(SPLLT_OMP_TRACE)
    integer :: th_id
#endif
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)    
    call spllt_starpu_insert_solve_block_c(bc_kk%hdl, bc_ik%hdl, p)
    ! call spllt_starpu_insert_solve_block(bc_kk, bc_ik, 0)

#elif defined(SPLLT_USE_OMP)
    
    bc_kk_c => bc_kk%c
    bc_ik_c => bc_ik%c
    
!$omp task private(m, n, sa, d_m, d_n, d_sa) &
!$omp    & firstprivate(bc_kk, bc_ik) & 
!$omp    & firstprivate(bc_kk_c, bc_ik_c) & 
!$omp    & private(lcol, bcol) &
!$omp    & shared(lfact) &
!$omp    & depend(in:bc_kk_c(1:1)) &
!$omp    & depend(inout:bc_ik_c(1:1))

! !$omp    & depend(in:lcol(d_sa:d_n*d_m)) &
! !$omp    & depend(inout:lcol(sa:n*m-1))

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(slv_blk_id, th_id)
#endif

    ! bcol is block column that blk and dblk belong to
    bcol = bc_kk%blk%bcol    
    lcol => lfact(bcol)%lcol

    ! bc_kk
    d_m  = bc_kk%blk%blkm
    d_n  = bc_kk%blk%blkn
    d_sa = bc_kk%blk%sa

    ! bc_ik
    n  = bc_ik%blk%blkn
    m  = bc_ik%blk%blkm
    sa = bc_ik%blk%sa

    call spllt_solve_block(m, n, &
         & lcol(sa:sa+n*m-1), & 
         & lcol(d_sa:d_sa+d_n*d_m))

#if defined(SPLLT_OMP_TRACE)
    ! write(*,*)"fac_blk_id: ", fac_blk_id
    call trace_event_stop(slv_blk_id, th_id)
#endif

!$omp end task

#else    

    ! bc_kk
    d_m  = bc_kk%blk%blkm
    d_n  = bc_kk%blk%blkn
    d_sa = bc_kk%blk%sa
    d_id = bc_kk%blk%id

    ! bc_ik
    n  = bc_ik%blk%blkn
    m  = bc_ik%blk%blkm
    sa = bc_ik%blk%sa
    id = bc_ik%blk%id
    
    ! bcol is block column that blk and dblk belong to
    bcol = bc_kk%blk%bcol    

    ! solve_block task
    ! call solv_col_block(m, n, id, & 
    !      & lfact(bcol)%lcol(sa:sa+n*m-1), &
    !      & d_id, lfact(bcol)%lcol(d_sa:d_sa+d_n*d_m), &
    !      & control)

    call spllt_solve_block(m, n, &
         & lfact(bcol)%lcol(sa:sa+n*m-1), & 
         & lfact(bcol)%lcol(d_sa:d_sa+d_n*d_m))

#endif

    return
  end subroutine spllt_solve_block_task

  ! syrk/gemm (same node)
  ! A_ij <- A_ij - A_ik A_jk^T
  subroutine spllt_update_block_task(bc_ik, bc_jk, bc_ij, lfact, prio)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
    !$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none
    
    ! type(block_type), intent(inout) :: bc_ik, bc_jk, bc_ij ! block to be updated    
    type(spllt_bc_type), intent(inout) :: bc_ik, bc_jk, bc_ij
#if defined(SPLLT_USE_OMP)
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
#else
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    integer, optional :: prio 

    integer :: n1, m1, sa1, n2, m2, sa2, n, m, sa
    integer :: bcol1, bcol, bcol2
    integer :: p, d
#if defined(SPLLT_USE_OMP)
    logical :: is_diag
    real(wp), dimension(:), pointer :: lcol1, lcol2, lcol
    real(wp), dimension(:), pointer :: bc_ik_c, bc_jk_c, bc_ij_c
#if defined(SPLLT_OMP_TRACE)
    integer :: th_id
#endif
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)
    if (bc_ij%blk%dblk .eq. bc_ij%blk%id) then
       d = 1
    else
       d = 0
    end if

    call spllt_starpu_insert_update_block_c(bc_ik%hdl, bc_jk%hdl, bc_ij%hdl, d, p)
#elif defined(SPLLT_USE_OMP)

    bc_ik_c => bc_ik%c
    bc_jk_c => bc_jk%c
    bc_ij_c => bc_ij%c
    
!$omp task private(n1, m1, sa1, bcol1, n2, m2, sa2, bcol2, &
!$omp    & n, m, sa, bcol, lcol, lcol1) shared(lfact) &
!$omp    & private(is_diag) &
!$omp    & firstprivate(bc_ik, bc_jk, bc_ij) &
!$omp    & firstprivate(bc_ik_c, bc_jk_c, bc_ij_c) &
!$omp    & depend(in:bc_ik_c(1:1), bc_jk_c(1:1)) &
!$omp    & depend(inout: bc_ij_c(1:1))

! !$omp    & depend(in:lcol1(sa1:n1*m1-1)) &
! !$omp    & depend(in:lcol2(sa2:n1*m2-1)) &
! !$omp    & depend(inout: lcol(sa:n*m-1))

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(upd_blk_id, th_id)
#endif

        ! bc_ik
    n1  = bc_ik%blk%blkn
    m1  = bc_ik%blk%blkm
    sa1 = bc_ik%blk%sa
    bcol1 = bc_ik%blk%bcol
    lcol1 => lfact(bcol1)%lcol

    ! bc_jk
    n2  = bc_jk%blk%blkn
    m2  = bc_jk%blk%blkm
    sa2 = bc_jk%blk%sa
    bcol2 = bc_jk%blk%bcol
    lcol2 => lfact(bcol2)%lcol

    ! bc_ij
    n  = bc_ij%blk%blkn
    m  = bc_ij%blk%blkm
    sa = bc_ij%blk%sa
    bcol = bc_ij%blk%bcol
    lcol => lfact(bcol)%lcol

    is_diag = bc_ij%blk%dblk.eq.bc_ij%blk%id 

    ! write(*,*)"bc_ik id: ", bc_ik%id, ", m: ", m1, ", n: ", n1
    ! write(*,*)"bc_jk id: ", bc_jk%id, ", m: ", m2, ", n: ", n2
    ! write(*,*)"bc_ij id: ", bc_ij%id, ", m: ", m, ", n: ", n

    ! write(*,*) "size lfact(bcol1)%lcol: ", size(lfact(bcol1)%lcol)
    ! write(*,*) "sa2+n1*m2-1: ", sa2+n1*m2-1

    call spllt_update_block(m, n, &
         & lcol(sa:sa+n*m-1), &
         & is_diag, n1, &
         & lcol1(sa2:sa2+n1*m2-1), &
         & lcol1(sa1:sa1+n1*m1-1))

#if defined(SPLLT_OMP_TRACE)
    ! write(*,*)"fac_blk_id: ", fac_blk_id
    call trace_event_stop(upd_blk_id, th_id)
#endif

!$omp end task

#else

    ! bc_ik
    n1  = bc_ik%blk%blkn
    m1  = bc_ik%blk%blkm
    sa1 = bc_ik%blk%sa
    bcol1 = bc_ik%blk%bcol
    
    ! bc_jk
    n2  = bc_jk%blk%blkn
    m2  = bc_jk%blk%blkm
    sa2 = bc_jk%blk%sa
    bcol2 = bc_jk%blk%bcol
    
    ! bc_ij
    n  = bc_ij%blk%blkn
    m  = bc_ij%blk%blkm
    sa = bc_ij%blk%sa
    
    bcol = bc_ij%blk%bcol

    ! write(*,*)"bc_ik id: ", bc_ik%id, ", m: ", m1, ", n: ", n1
    ! write(*,*)"bc_jk id: ", bc_jk%id, ", m: ", m2, ", n: ", n2
    ! write(*,*)"bc_ij id: ", bc_ij%id, ", m: ", m, ", n: ", n
    
    ! write(*,*) "size lfact(bcol1)%lcol: ", size(lfact(bcol1)%lcol)
    ! write(*,*) "sa2+n1*m2-1: ", sa2+n1*m2-1

    call spllt_update_block(m, n, &
         & lfact(bcol)%lcol(sa:sa+n*m-1), &
         & bc_ij%blk%dblk.eq.bc_ij%blk%id, n1, &
         & lfact(bcol1)%lcol(sa2:sa2+n1*m2-1), &
         & lfact(bcol1)%lcol(sa1:sa1+n1*m1-1))

#endif

    return
  end subroutine spllt_update_block_task

  ! syrk/gemm (inter-node)
  subroutine spllt_update_between_task(dbc, snode, a_bc, anode, &
       ! & csrc, rsrc, &
       & cptr, cptr2, rptr, rptr2, &
       & row_list, col_list, workspace, &
       & lfact, blocks, bcs, &
       & control, prio)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
!$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none

    ! type(block_type), intent(inout)     :: a_bc ! dest block
    ! type(block_type), intent(in)        :: bc ! src block
    type(spllt_bc_type), intent(inout)     :: a_bc ! dest block
    ! type(spllt_bc_type), intent(in)        :: bc ! src block
    type(spllt_bc_type), intent(in)        :: dbc ! diag block in source node
#if defined(SPLLT_USE_STARPU)
    type(node_type), target                :: snode ! src node
    type(spllt_node_type), target          :: anode ! dest node
#else
    type(node_type)                        :: snode ! src node
    type(spllt_node_type)                  :: anode ! dest node
#endif
    integer :: cptr, cptr2, rptr, rptr2 
!    integer :: csrc(2), rsrc(2) ! used for update_between tasks to
    integer, dimension(:), allocatable  :: row_list, col_list
#if defined(SPLLT_USE_OMP)
    type(spllt_bc_type), allocatable :: workspace(:)
    type(spllt_bc_type), target      :: bcs(:) ! block info.
#else
    type(spllt_bc_type)              :: workspace
    type(spllt_bc_type)              :: bcs(:) ! block info.
#endif
    ! real(wp), dimension(:), allocatable :: buffer ! update_buffer workspace
    type(block_type), dimension(:)      :: blocks ! block info. 

#if defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)
    type(lfactor), allocatable, target, intent(inout) :: lfact(:)
#else
    type(lfactor), allocatable, intent(inout) :: lfact(:)
#endif
    type(MA87_control), intent(in) :: control     
    ! integer, intent(out) :: st ! TODO error managment
    integer, optional :: prio 

    ! integer :: info ! TODO error managment
    
    integer :: p
    integer :: m, n, n1, sa
    integer :: bcol, bcol1
    integer(long) :: id, id1
    integer :: dcol, scol
    integer :: csrc, csrc2, rsrc, rsrc2
    integer :: s_nb

#if defined(SPLLT_USE_STARPU) || defined(SPLLT_USE_OMP)
    integer(long) :: blk_sa, blk_en, nb_blk, dblk 
#endif  

#if defined(SPLLT_USE_STARPU)
    integer :: blkn, ljk_sa
    integer :: nhljk, nhlik
    integer(c_int) :: ljk_m, ljk_n
    integer(long) :: blk
    type(c_ptr), dimension(:), allocatable :: lik_handles, ljk_handles
    ! type(c_ptr) :: ljk_hdl
    type(c_ptr) :: snode_c, anode_c, a_bc_c
#endif

#if defined(SPLLT_USE_OMP)
    real(wp), dimension(:), pointer :: lcol1, lcol
    type(block_type), pointer :: a_blk ! block to be updated
    integer, dimension(:), allocatable  :: rlst, clst
    integer :: th_id
    integer :: csrc_sa, rsrc_sa, csrc_en, rsrc_en
    ! type(spllt_bc_type), pointer :: bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en
    real(wp), dimension(:), pointer :: a_bc_c, bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en, dbc_c
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

    s_nb = snode%nb    ! block size in source node
    n1  = dbc%blk%blkn ! width of column

    bcol1 = dbc%blk%bcol
    scol = bcol1 - blocks(snode%blk_sa)%bcol + 1

    bcol = a_bc%blk%bcol
    dcol = bcol - blocks(anode%node%blk_sa)%bcol + 1
    ! dcol = a_bc%blk%bcol ! blocks(anode%node%blk_sa)%bcol + 1

#if defined(SPLLT_USE_STARPU)
! #if 0
    ! write(*,*)"associated blk: ", associated(a_bc%blk)
    dblk = dbc%blk%id

    ! ljk
    blk_sa = (cptr -1)/s_nb - (scol-1) + dblk
    blk_en = (cptr2-1)/s_nb - (scol-1) + dblk
    nb_blk = blk_en-blk_sa+1

    allocate(ljk_handles(nb_blk))
    ! ljk_m = 0 ! compute height of ljk factor
    nhljk=0
    do blk=blk_sa,blk_en
       nhljk=nhljk+1
       ljk_handles(nhljk) = bcs(blk)%hdl
       ! ljk_m = ljk_m + bcs(blk)%blk%blkm
    end do
    ! ljk_n = bc%blk%blkn

    ! lik
    blk_sa = (rptr -1)/s_nb - (scol-1) + dblk
    blk_en = (rptr2-1)/s_nb - (scol-1) + dblk
    nb_blk = blk_en-blk_sa+1

    allocate(lik_handles(nb_blk))
    nhlik=0
    do blk=blk_sa,blk_en
       nhlik=nhlik+1
       lik_handles(nhlik) = bcs(blk)%hdl
    end do
    
    csrc  = 1 + (mod(cptr-1, s_nb))*n1
    csrc2 = (cptr2 - cptr + 1)*n1

    rsrc  = 1 + (mod(rptr-1, s_nb))*n1
    rsrc2 = (rptr2 - rptr + 1)*n1

    snode_c = c_loc(snode)
    anode_c = c_loc(anode%node)
    a_bc_c  = c_loc(a_bc%blk)
    
    ! write(*,*)"dcol: ", dcol

    call spllt_starpu_insert_update_between_c(&
         & lik_handles, nhlik, &
         & ljk_handles, nhljk, &
         & a_bc%hdl, &
         & snode_c, scol, &
         & anode_c, a_bc_c, dcol, &
         & csrc, csrc2, rsrc, rsrc2, &
         & control%min_width_blas, &
         & workspace%hdl, &
         & anode%hdl, &
         & p &
         &)

    ! call test_insert_c(lik_handles, nhlik, ljk_handles, nhljk)

    ! call starpu_f_task_wait_for_all()

    deallocate(ljk_handles)
    deallocate(lik_handles)
    
#elif defined(SPLLT_USE_OMP)


    ! ! lik
    ! blk_sa = (rptr -1)/s_nb - (scol-1) + dblk
    ! blk_en = (rptr2-1)/s_nb - (scol-1) + dblk
    ! nb_blk = blk_en-blk_sa+1
    
    ! write(*,*) "update_between_task"
    ! call update_between(m, n, id, anode, &
    !      & n1, id1, snode, &
    !      & lfact(bcol)%lcol(sa:sa+m*n-1), &
    !      & lfact(bcol1)%lcol(csrc(1):csrc(1)+csrc(2)-1), &
    !      & lfact(bcol1)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1), &
    !      & blocks, row_list, col_list, buffer, &
    !      & control, info, st)

    dblk = dbc%blk%id
    ! ljk
    blk_sa = (cptr -1)/s_nb - (scol-1) + dblk
    blk_en = (cptr2-1)/s_nb - (scol-1) + dblk
    nb_blk = blk_en-blk_sa+1
    bc_jk_sa => bcs(blk_sa)%c
    bc_jk_en => bcs(blk_en)%c
    ! write(*,*) "nb_blk: ", nb_blk
    ! write(*,*) "nb_blk:", nb_blk, "blk_sa: ", blk_sa, ", blk_en: ", blk_en
    csrc_sa = blocks(blk_sa)%sa
    csrc_en = blocks(blk_en)%sa

    ! csrc_sa = 1 + (cptr_sa-(scol-1)*s_nb-1)*n1
    ! csrc_le = (cptr_en - cptr_sa + 1)*n1

    ! lik
    blk_sa = (rptr -1)/s_nb - (scol-1) + dblk
    blk_en = (rptr2-1)/s_nb - (scol-1) + dblk
    nb_blk = blk_en-blk_sa+1
    bc_ik_sa => bcs(blk_sa)%c
    bc_ik_en => bcs(blk_en)%c

    rsrc_sa = blocks(blk_sa)%sa
    rsrc_en = blocks(blk_en)%sa

    a_bc_c => a_bc%c

    dbc_c => dbc%c
! !$omp taskwait
    
!$omp task firstprivate(dbc, rlst, clst, cptr, cptr2, rptr, rptr2, &
!$omp    & dcol, scol, snode, anode) &    
!$omp    & firstprivate(a_bc, bc_jk_sa, bc_jk_en, bc_ik_sa, bc_ik_en) &
!$omp    & private(lcol, lcol1, m, n, sa, th_id, bcol, bcol1, &
!$omp    & a_blk, s_nb, n1) &    
!$omp    & shared(workspace, control, lfact) &
!$omp    & depend(in:bc_jk_sa(1:1), bc_jk_en(1:1)) &
!$omp    & depend(in:bc_ik_sa(1:1), bc_ik_en(1:1)) &
! !$omp    & depend(inout:dbc_c(1:1)) &
!$omp    & depend(inout:a_bc_c(1:1))

! !$omp    & depend(inout:lcol(sa:m*n-1)) &
! !$omp    & depend(in:lcol1(csrc_sa:1), lcol1(csrc_en:1)) &
! !$omp    & depend(in:lcol1(rsrc_sa:1), lcol1(rsrc_en:1))

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"upd_btw_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(upd_btw_id, th_id)
#endif

    s_nb = snode%nb    ! block size in source node
    n1  = dbc%blk%blkn ! width of column
    
    bcol1 = dbc%blk%bcol
    bcol  = a_bc%blk%bcol

    lcol1 => lfact(bcol1)%lcol
    lcol  => lfact(bcol)%lcol

    ! write(*,*)"thread id: ", th_id
    m  = a_bc%blk%blkm
    n  = a_bc%blk%blkn
    sa = a_bc%blk%sa

    th_id = omp_get_thread_num()
    ! write(*,*)"thread id: ", th_id

    a_blk => a_bc%blk

    csrc  = 1 + (cptr-(scol-1)*s_nb-1)*n1
    csrc2 = (cptr2 - cptr + 1)*n1

    rsrc  = 1 + (rptr-(scol-1)*s_nb-1)*n1
    rsrc2 = (rptr2 - rptr + 1)*n1

    allocate(rlst(1), clst(1))

    call spllt_update_between(m, n, a_blk, dcol, anode%node, &
         & n1, scol, snode, &
         & lcol(sa:sa+m*n-1), &
         & lcol1(csrc:csrc+csrc2-1), &
         & lcol1(rsrc:rsrc+rsrc2-1), &
         & rlst, clst, workspace(th_id)%c, &
         & control%min_width_blas)

    deallocate(rlst, clst)

#if defined(SPLLT_OMP_TRACE)
    ! write(*,*)"fac_blk_id: ", fac_blk_id
    call trace_event_stop(upd_btw_id, th_id)
#endif

!$omp end task
    
#else
    
    m  = a_bc%blk%blkm
    n  = a_bc%blk%blkn
    id = a_bc%blk%id
    sa = a_bc%blk%sa
    
    ! write(*,*) "update_between_task"
    ! call update_between(m, n, id, anode, &
    !      & n1, id1, snode, &
    !      & lfact(bcol)%lcol(sa:sa+m*n-1), &
    !      & lfact(bcol1)%lcol(csrc(1):csrc(1)+csrc(2)-1), &
    !      & lfact(bcol1)%lcol(rsrc(1):rsrc(1)+rsrc(2)-1), &
    !      & blocks, row_list, col_list, buffer, &
    !      & control, info, st)

    csrc  = 1 + (cptr-(scol-1)*s_nb-1)*n1
    csrc2 = (cptr2 - cptr + 1)*n1

    rsrc  = 1 + (rptr-(scol-1)*s_nb-1)*n1
    rsrc2 = (rptr2 - rptr + 1)*n1

    call spllt_update_between(m, n, a_bc%blk, dcol, anode%node, &
         & n1, scol, snode, &
         & lfact(bcol)%lcol(sa:sa+m*n-1), &
         & lfact(bcol1)%lcol(csrc:csrc+csrc2-1), &
         & lfact(bcol1)%lcol(rsrc:rsrc+rsrc2-1), &
         & row_list, col_list, workspace%c, &
         & control%min_width_blas)

#endif

    return
  end subroutine spllt_update_between_task

  ! init node
  subroutine spllt_init_node_task(node, n, val, map, keep, prio)
    use spllt_mod
    use hsl_ma87_double
    use spllt_kernels_mod
#if defined(SPLLT_USE_STARPU)
    use spllt_starpu_factorization_mod
#elif defined(SPLLT_USE_OMP)
!$ use omp_lib
#if defined(SPLLT_OMP_TRACE)
    use trace_mod
#endif
#endif
    implicit none
    
    type(spllt_node_type), intent(in) :: node
    integer, intent(in) :: n      ! order of matrix 
    real(wp), dimension(*), target, intent(in) :: val ! user's matrix values
    integer, target, intent(out) :: map(n)  ! mapping array. Reset for each node
     ! so that, if variable (row) i is involved in node,
     ! map(i) is set to its local row index
    type(MA87_keep), target, intent(inout) :: keep ! on exit, matrix a copied
     ! into relevant part of keep%lfact
    integer, optional :: prio 

    integer :: p
#if defined(SPLLT_USE_STARPU)
    type(c_ptr) :: val_c, map_c, keep_c
#endif
#if defined(SPLLT_USE_OMP)
    integer :: snum
    integer :: th_id
#endif

    if (present(prio)) then
       p = prio
    else
       p = 0
    end if

#if defined(SPLLT_USE_STARPU)
    val_c  = c_loc(val(1)) 
    map_c  = c_loc(map(1)) 
    keep_c = c_loc(keep)

    call spllt_insert_init_node_task_c(node%hdl, &
         & node%num, n, val_c, map_c, keep_c, p)
#elif defined(SPLLT_USE_OMP)

! !$omp task default(shared) private(snum) shared(n, val, map, keep)
!$omp task private(snum) shared(n, val, map, keep) 

#if defined(SPLLT_OMP_TRACE)
    th_id = omp_get_thread_num()
    ! write(*,*)"fac_blk_id: ", fac_blk_id, ", th_id: ", th_id
    call trace_event_start(ini_nde_id, th_id)
#endif

    snum = node%num

    call spllt_init_node(snum, n, val, map, keep)

#if defined(SPLLT_OMP_TRACE)
    call trace_event_stop(ini_nde_id, th_id)
#endif

!$omp end task
! !$omp taskwait
#else
    call spllt_init_node(node%num, n, val, map, keep)
#endif
    
    return
  end subroutine spllt_init_node_task

end module spllt_factorization_mod
