module spllt_solve_dep_mod

! subroutine spllt_compute_blk_solve_dep(fkeep, blk)
! subroutine spllt_compute_solve_dep(fkeep)
! subroutine get_update_dep(fkeep, child_node, node, pos)
! subroutine get_update_nblk(fkeep, child_node, blk_index, nblk)
! subroutine get_update_dep_blk(fkeep, child_node, blk_index, pos)
! subroutine reduce_ind_and_get_ndep(fkeep, ind, nind, child_node, ndep)
! recursive subroutine getUpdateNDep(fkeep, node, ind, nind, ndep, lvl)
! subroutine reduce_ind_and_get_dep(fkeep, ind, nind, child_node, dep)
! recursive subroutine getUpdateDep(fkeep, node, ind, nind, dep, ndep, lvl)
! subroutine reduce_ind_and_get_dep_blk(fkeep, ind, nind, child_node, blk_dep)
! subroutine getPointerBlkIndex(fkeep, blk, p)
! subroutine getSolveNDep(fkeep, node, ind, nind, ndep)
! subroutine getSolveDep(fkeep, node, ind, nind, dep, ndep)
! subroutine fwd_update_dependency(fkeep, blk, dep)
! subroutine bwd_solve_dependency(fkeep, blk, dep)

  use spllt_data_mod

contains

  subroutine spllt_compute_blk_solve_dep(fkeep, blk)
   !use utils_mod
    type(spllt_fkeep), target, intent(inout)  :: fkeep
    integer, intent(in)                       :: blk

    integer :: ndep, nldep, repr, old_node, cpt
    integer :: dep_blk, dep_node
    integer, pointer :: p_small(:)

    p_small => fkeep%small

    !!!!!!!!!!!!!!!!!!!!!
    ! FWD dep
    !
    call fwd_update_dependency(fkeep, blk, fkeep%bc(blk)%fwd_update_dep)
    fkeep%bc(blk)%fwd_solve_dep   = fwd_solve_dependency(fkeep, blk)

    ndep  = 0

    if(p_small(fkeep%bc(blk)%node) .eq. 0) then

      nldep = 0
      old_node  = 0
      if(fkeep%bc(blk)%fwd_update_dep(1) .ne. blk) then
        do i = 1, size(fkeep%bc(blk)%fwd_update_dep)
          dep_blk   = fkeep%bc(blk)%fwd_update_dep(i)
          dep_node  = fkeep%bc(dep_blk)%node
          repr      = p_small(dep_node)
          ! If the blk index is part of a subtree, count
          if(repr .lt. 0) then 
            repr = -repr
            if(repr .gt. old_node) then
              nldep = nldep + 1
              old_node = repr
            end if
          else if(repr .eq. 1) then
            if(dep_node .gt. old_node) then
              nldep = nldep + 1
              old_node = dep_node
            end if
          else
            nldep = nldep + 1
          end if
        end do
        ndep = ndep + nldep
      end if
      if(fkeep%bc(blk)%fwd_solve_dep .ne. blk) then
        ndep = ndep + 1
      end if
    else
      if(fkeep%bc(blk)%fwd_update_dep(1) .ne. blk) then
        ndep = ndep + size(fkeep%bc(blk)%fwd_update_dep)
      end if
      if(fkeep%bc(blk)%fwd_solve_dep .ne. blk) then
        ndep = ndep + 1
      end if
    end if
    
#if defined(SOLVE_TASK_LOCKED)
    if(ndep .gt. 0) then
#endif
      allocate(fkeep%bc(blk)%fwd_dep(ndep))
      if(p_small(fkeep%bc(blk)%node) .eq. 0) then
        cpt = 1
        old_node  = 0
        if(fkeep%bc(blk)%fwd_update_dep(1) .ne. blk) then
          do i = 1, size(fkeep%bc(blk)%fwd_update_dep)
            dep_blk   = fkeep%bc(blk)%fwd_update_dep(i)
            dep_node  = fkeep%bc(dep_blk)%node
            repr      = p_small(dep_node)
            ! If the blk index is part of a subtree, count
            if(repr .lt. 0) then
              repr = -repr
              if(repr .gt. old_node) then
                fkeep%bc(blk)%fwd_dep(cpt) = fkeep%nodes(repr)%blk_en
                old_node = repr
                cpt = cpt + 1
              end if
            else if(repr .eq. 1) then
              if(dep_node .gt. old_node) then
                fkeep%bc(blk)%fwd_dep(cpt) = fkeep%nodes(dep_node)%blk_en
                old_node = dep_node
                cpt = cpt + 1
              end if
            else
              fkeep%bc(blk)%fwd_dep(cpt) = fkeep%bc(blk)%fwd_update_dep(i)
              cpt = cpt + 1
            end if
          end do
        end if
        if(fkeep%bc(blk)%fwd_solve_dep .ne. blk) then
          dep_blk   = fkeep%bc(blk)%fwd_solve_dep
          dep_node  = fkeep%bc(dep_blk)%node
          repr      = p_small(dep_node)
          ! If the blk index is part of a subtree, count
          if(repr .lt. 0) then
            repr = -repr
            if(repr .gt. old_node) then
              fkeep%bc(blk)%fwd_dep(ndep) = fkeep%nodes(repr)%blk_en
            end if
          else if(repr .eq. 1) then
            if(dep_node .gt. old_node) then
              fkeep%bc(blk)%fwd_dep(cpt) = fkeep%nodes(dep_node)%blk_en
              old_node = dep_node
            end if
          else
            fkeep%bc(blk)%fwd_dep(ndep) = fkeep%bc(blk)%fwd_solve_dep
          end if
        end if
      else
        if(fkeep%bc(blk)%fwd_update_dep(1) .ne. blk) then
          fkeep%bc(blk)%fwd_dep(1:size(fkeep%bc(blk)%fwd_update_dep)) = &
            fkeep%bc(blk)%fwd_update_dep
        end if
        if(fkeep%bc(blk)%fwd_solve_dep .ne. blk) then
          fkeep%bc(blk)%fwd_dep(ndep) = fkeep%bc(blk)%fwd_solve_dep
        end if
      end if
#if defined(SOLVE_TASK_LOCKED)
    else
      if(blk .eq. fkeep%bc(blk)%dblk) then
        allocate(fkeep%bc(blk)%fwd_dep(1))
        fkeep%bc(blk)%fwd_dep(1) = 1
      end if
    endif
#endif

    !!!!!!!!!!!!!!!!!!!!!
    ! BWD dep
    !
    !Compute basic dependencies
    fkeep%bc(blk)%bwd_update_dep  = bwd_update_dependency(fkeep, blk)
    call bwd_solve_dependency(fkeep, blk, fkeep%bc(blk)%bwd_solve_dep)

    ! Count the number of dep different to the block blk in order to
    ! allocate the right size and copy them into it
    ndep = 0
    if(fkeep%bc(blk)%bwd_solve_dep(1) .ne. blk) then
      ndep = ndep + size(fkeep%bc(blk)%bwd_solve_dep)
    end if
    if(fkeep%bc(blk)%bwd_update_dep .ne. blk) then
      ndep = ndep + 1
    end if
    
#if defined(SOLVE_TASK_LOCKED)
    if(ndep .gt. 0) then
#endif
      allocate(fkeep%bc(blk)%bwd_dep(ndep))
      
      ! Copy array if the first element is not blk
      ! This assumes that the content is either only blk or a list that
      ! does not contain blk
      if(fkeep%bc(blk)%bwd_solve_dep(1) .ne. blk) then
        fkeep%bc(blk)%bwd_dep(1:size(fkeep%bc(blk)%bwd_solve_dep)) = &
          fkeep%bc(blk)%bwd_solve_dep
      end if
      ! Copy the dependency only if it is not blk
      if(fkeep%bc(blk)%bwd_update_dep .ne. blk) then
        fkeep%bc(blk)%bwd_dep(ndep) = fkeep%bc(blk)%bwd_update_dep
      end if
#if defined(SOLVE_TASK_LOCKED)
    else
      if(blk .eq. fkeep%bc(blk)%dblk) then
        allocate(fkeep%bc(blk)%bwd_dep(1))
        fkeep%bc(blk)%bwd_dep(1) = fkeep%info%num_nodes
      end if
    end if
#endif
  end subroutine spllt_compute_blk_solve_dep



  ! This routine computes the dependencies of the blocks in fkeep
  subroutine spllt_compute_solve_dep(fkeep)
    use spllt_tree_stat_mod 
    type(spllt_fkeep), target, intent(inout)  :: fkeep

    integer                 :: i
    type(spllt_tree_stat_t) :: task_info_fwd
    type(spllt_tree_stat_t) :: task_info_bwd

    call spllt_init_tree_stat(task_info_fwd)
    call spllt_init_tree_stat(task_info_bwd)

    ! For each block i, compute its dependencies for forward and 
    ! backward steps,
    do i = 1, fkeep%nodes(fkeep%info%num_nodes)%blk_en
      call spllt_compute_blk_solve_dep(fkeep, i)
    end do
    ! then, update the list of dependencies for the backward case
    ! and for each tree.
    do i = 1, size(fkeep%trees)
      call spllt_update_blk_dep(fkeep, fkeep%trees(i))
    end do

    ! Statistics
    do i = 1, fkeep%nodes(fkeep%info%num_nodes)%blk_en
      call spllt_update_tree_stat(task_info_fwd, size(fkeep%bc(i)%fwd_dep))
      call spllt_update_tree_stat(task_info_bwd, size(fkeep%bc(i)%bwd_dep))
    end do

    call print_tree_stat(task_info_fwd, "FWD STAT")
    call print_tree_stat(task_info_bwd, "BWD STAT")
  end subroutine spllt_compute_solve_dep



  subroutine spllt_compute_solve_extra_row(fkeep)
    type(spllt_fkeep), intent(inout)  :: fkeep

    integer :: i

    do i = 1, fkeep%info%num_nodes
      call spllt_compute_extra_row(fkeep, i)
    end do
  end subroutine spllt_compute_solve_extra_row



  !Return the position in child_node of the rows of node that are 
  ! in child_node
  subroutine get_update_dep(fkeep, child_node, node, pos)
    
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(in)                   :: node
    integer, allocatable, intent(out)     :: pos(:)

    integer :: i, j, k
    integer, pointer, dimension(:) :: p_node_index, p_child_node_index

    p_node_index        => fkeep%nodes(node)%index
    p_child_node_index  => fkeep%nodes(child_node)%index


    i = 1
    j = 1
    k = 1
    do while(j .le. size(p_node_index) .and. k .le. size(p_child_node_index))
      if(p_node_index(j) .lt. p_child_node_index(k)) then
        j = j + 1
      else if(p_node_index(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        i = i + 1
        j = j + 1
        k = k + 1
      end if
    end do

    allocate(pos(i - 1))

    if( i .eq. 1 )then
      return
    end if

    !Reiterate to copy the position into pos
    i = 1
    j = 1
    k = 1
    do while(j .le. size(p_node_index) .and. k .le. size(p_child_node_index))
      if(p_node_index(j) .lt. p_child_node_index(k)) then
        j = j + 1
      else if(p_node_index(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        pos(i) = j
        i = i + 1
        j = j + 1
        k = k + 1
      end if
    end do

  end subroutine get_update_dep



  subroutine get_update_nblk(fkeep, child_node, blk_index, nblk)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(inout)                :: blk_index(:)
    integer, intent(out)                  :: nblk

    integer :: j, k
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp

    p_child_node_index  => fkeep%nodes(child_node)%index
    cur_blk_dep = 0 !Impossible value but used as initialization

    nblk  = 0
    j     = 1
    k     = 1
    do while(j .le. size(blk_index) .and. k .le. size(p_child_node_index))
      if(blk_index(j) .lt. p_child_node_index(k)) then
        j = j + 1
      else if(blk_index(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else

        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))

        if(cur_blk_dep .lt. tmp) then
          cur_blk_dep = tmp
          nblk = nblk + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

  end subroutine get_update_nblk



  subroutine get_update_dep_blk(fkeep, child_node, blk_index, pos)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(inout)                :: blk_index(:)
    integer, intent(out)                  :: pos(:)

    integer :: i, j, k
    integer :: nb
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp

    p_child_node_index  => fkeep%nodes(child_node)%index
    nb = fkeep%nodes(child_node)%nb
    cur_blk_dep = 0 !Impossible value but used as initialization

    !Set variables
    cur_blk_dep = 0
    i           = 1
    j           = 1
    k           = 1
    do while(j .le. size(blk_index) .and. k .le. size(p_child_node_index))
      if(blk_index(j) .lt. p_child_node_index(k)) then
        j = j + 1
      else if(blk_index(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        
        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))

        if(cur_blk_dep .lt. tmp) then
          cur_blk_dep = tmp
          pos(i) = tmp
          i = i + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

  end subroutine get_update_dep_blk



  function get_child_dep_blk_id(fkeep, child_node, row, nrow)
    integer                               :: get_child_dep_blk_id
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(in)                   :: row
    integer, intent(in)                   :: nrow

    integer :: lblk, nlblk, diff, tmp, nb

    nb    = fkeep%nodes(child_node)%nb
    tmp   = fkeep%nodes(child_node)%blk_sa
    lblk  = ceiling((row + 0.0) / nb)
    nlblk = ceiling((nrow + 0.0)/nb)
    diff  = (fkeep%nodes(child_node)%blk_en - &
      fkeep%bc(fkeep%nodes(child_node)%blk_en)%dblk + 0) ! #blk on the last nbol
    if(lblk .le. (nlblk - diff)) then
      tmp = tmp + (lblk - 1) * ( nlblk + 1 - 0.5 * lblk )
    else
      tmp = fkeep%bc(fkeep%nodes(child_node)%blk_en)%dblk + lblk - &
        (nlblk - diff) - 0
    end if

    get_child_dep_blk_id = tmp

  end function get_child_dep_blk_id



  subroutine reduce_ind_and_get_ndep(fkeep, ind, nind, child_node, ndep)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(inout)                :: ind(:)
    integer, intent(inout)                :: nind       ! #elmt in ind
    integer, intent(in)                   :: child_node
    integer, intent(out)                  :: ndep       ! #block dep

    integer :: i, j, k
    integer :: nb, nval
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp

    p_child_node_index  => fkeep%nodes(child_node)%index
    nb                  = fkeep%nodes(child_node)%nb
    cur_blk_dep         = 0 !Impossible value but used as initialization
    ndep                = 0

    !Set variables
    cur_blk_dep = 0
    i           = 1
    j           = 1
    k           = 1
    nval        = 0
    do while(j .le. nind .and. k .le. size(p_child_node_index))
      if(ind(j) .lt. p_child_node_index(k)) then
        nval = nval + 1
        ind(nval) = ind(j)
        j = j + 1
      else if(ind(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        
        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))

        
        if(cur_blk_dep .lt. tmp) then
!         print *, "Dep with blk ", tmp
          cur_blk_dep = tmp
          ndep = ndep + 1
          i = i + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

    if(j .lt. nind) then
      nval = nval + 1
      ind(nval: nval + nind - j) = ind(j : nind)
      nind = nval + nind - j
    else 
      nind = nval
    end if


  end subroutine reduce_ind_and_get_ndep



  recursive subroutine getUpdateNDep(fkeep, node, ind, nind, ndep, lvl)
    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: node
    integer, intent(inout)        :: ind(:)
    integer, intent(inout)        :: nind
    integer, intent(inout)        :: ndep(:)
    integer, intent(in), optional :: lvl

    integer :: child, nldep, offset, i, nsubind, nchild
    integer, allocatable :: subind(:)
    integer :: rlvl

    if(present(lvl)) then
      rlvl = lvl
    else
      rlvl = 1
    end if

    nchild = fkeep%nodes(node)%nchild

    if(nind .eq. 0) then
      ndep(1) = 0
      return
    end if

    if(nchild .eq. 0) then
      return
    end if

    offset = nchild + 1
    allocate(subind(nind))

    do i = 1, nchild
      child = fkeep%nodes(node)%child(i)
      nldep = 0

      subind = ind
      nsubind = nind

      call reduce_ind_and_get_ndep(fkeep, subind, nsubind, child, nldep)
      ndep(i) = nldep

      if(fkeep%nodes(child)%nchild .gt. 0) then

        call getUpdateNDep(fkeep, child, subind, nsubind, ndep(offset :       &
          offset + child - fkeep%nodes(child)%least_desc - 1), &
          rlvl + 1)

        offset = offset + fkeep%nodes(child)%nchild
      end if
    end do

    deallocate(subind)

  end subroutine getUpdateNDep



  !Get the dependencies after counting its number
  subroutine reduce_ind_and_get_dep(fkeep, ind, nind, child_node, dep)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(inout)                :: ind(:)
    integer, intent(inout)                :: nind       ! #elmt in ind
    integer, intent(in)                   :: child_node
    integer, intent(out)                  :: dep(:)     ! List of block dep

    integer :: i, j, k
    integer :: nb, nval
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp
    integer :: ndep

    p_child_node_index  => fkeep%nodes(child_node)%index
    nb                  = fkeep%nodes(child_node)%nb
    cur_blk_dep         = 0 !Impossible value but used as initialization
    ndep                = 1

    !Set variables
    cur_blk_dep = 0
    i           = 1
    j           = 1
    k           = 1
    nval        = 0
    do while(j .le. nind .and. k .le. size(p_child_node_index))
      if(ind(j) .lt. p_child_node_index(k)) then
        nval = nval + 1
        ind(nval) = ind(j)
        j = j + 1
      else if(ind(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        
        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))
        
        if(cur_blk_dep .lt. tmp) then
!         print *, "Dep with blk ", tmp
          cur_blk_dep = tmp
          dep(ndep) = tmp
          ndep = ndep + 1
          i = i + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

    if(j .lt. nind) then
      nval = nval + 1
      ind(nval: nval + nind - j) = ind(j : nind)
      nind = nval + nind - j
    else 
      nind = nval
    end if


  end subroutine reduce_ind_and_get_dep



  recursive subroutine getUpdateDep(fkeep, node, ind, nind, dep, ndep, lvl)
    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: node
    integer, intent(inout)        :: ind(:)
    integer, intent(inout)        :: nind
    integer, intent(inout)        :: dep(:)
    integer, intent(in)           :: ndep(:)
    integer, intent(in), optional :: lvl

    integer :: child, nldep, offset, i, nsubind, nchild
    integer, allocatable :: subind(:)
    integer :: rlvl

    if(present(lvl)) then
      rlvl = lvl
    else
      rlvl = 1
    end if

    nchild = fkeep%nodes(node)%nchild

    if(nind .eq. 0) then
      return
    end if

    if(nchild .eq. 0) then
      return
    end if

    offset = nchild + 1
    allocate(subind(nind))

    do i = 1, nchild
      child = fkeep%nodes(node)%child(i)
      nldep = 0

      subind = ind
      nsubind = nind

      if(ndep(i+1) .gt. ndep(i)) then

        call reduce_ind_and_get_dep(fkeep, subind, nsubind, child, &
          dep(ndep(i) : ndep(i + 1) - 1))

        if(fkeep%nodes(child)%nchild .gt. 0) then
          
          nldep = ndep(offset + child - fkeep%nodes(child)%least_desc ) - &
            ndep(offset)

          if(nldep .gt. 0) then
            call getUpdateDep(fkeep, child, subind, nsubind,      &
              dep(ndep(offset)  : ndep(offset) + nldep - 1),      &
              ndep(ndep(offset) : ndep(offset) + nldep - 1),      &
              rlvl + 1)

            offset = offset + nldep
          end if
        end if
      end if
    end do

    deallocate(subind)

  end subroutine getUpdateDep



  subroutine reduce_ind_and_get_dep_blk(fkeep, ind, nind, child_node, blk_dep)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(inout)                :: ind(:)
    integer, intent(inout)                :: nind       ! #elmt in ind
    integer, intent(in)                   :: child_node
    integer, intent(out)                  :: blk_dep(:) ! Dep of the block

    integer :: i, j, k
    integer :: nb, nval
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp

    p_child_node_index  => fkeep%nodes(child_node)%index
    nb = fkeep%nodes(child_node)%nb
    cur_blk_dep = 0 !Impossible value but used as initialization

    !Set variables
    cur_blk_dep = 0
    i           = 1
    j           = 1
    k           = 1
    nval        = 0
    do while(j .le. nind .and. k .le. size(p_child_node_index))
      if(ind(j) .lt. p_child_node_index(k)) then
        nval = nval + 1
        ind(nval) = ind(j)
        j = j + 1
      else if(ind(j) .gt. p_child_node_index(k)) then
        k = k + 1
      else
        
        tmp = get_child_dep_blk_id(fkeep, child_node, k, &
          size(p_child_node_index))
        
        if(cur_blk_dep .lt. tmp) then
          cur_blk_dep = tmp
          blk_dep(i) = tmp
          i = i + 1
        end if
        j = j + 1
        k = k + 1
      end if
    end do

    if(j .lt. nind) then
      nval = nval + 1
      ind(nval: nval + nind - j) = ind(j : nind)
      nind = nval + nind - j
    else 
      nind = nval
    end if

  end subroutine reduce_ind_and_get_dep_blk



  subroutine getPointerBlkIndex(fkeep, blk, p)
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: blk
    integer, pointer, intent(out)         :: p(:)

    integer :: node        
    integer :: nb          
    integer :: blkm        
    integer :: blk_sa      
    integer :: bcol_blk_sa 
    integer :: bcol        
    integer :: local_blk   
    integer :: blk_ind_sa  

    node        = fkeep%bc(blk)%node
    nb          = fkeep%nodes(node)%nb    
    blkm        = fkeep%bc(blk)%blkm
    blk_sa      = fkeep%nodes(node)%blk_sa
    bcol_blk_sa = fkeep%bc(blk_sa)%bcol
    bcol        = fkeep%bc(blk)%bcol
    local_blk   = blk - fkeep%bc(blk)%dblk ! In the bcol
    blk_ind_sa  = nb * (local_blk + (bcol - bcol_blk_sa)) + 1
    p           => fkeep%nodes(node)%index(blk_ind_sa : blk_ind_sa + blkm - 1)

  end subroutine getPointerBlkIndex



  subroutine getSolveNDep(fkeep, node, ind, nind, ndep)
    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: node
    integer, intent(inout)        :: ind(:)
    integer, intent(inout)        :: nind
    integer, intent(inout)        :: ndep(:)

    integer :: nldep, i
    integer :: parent, nparent

    parent = node
    nparent = 0
    do while(parent .lt. fkeep%info%num_nodes)
      nparent = nparent + 1
      parent = fkeep%nodes(parent)%parent
    end do

    if(nparent .eq. 0) then
      ndep(1) = 0
      return
    end if

    parent = node

    do i = 1, nparent
      parent = fkeep%nodes(parent)%parent
      nldep = 0

      if(parent .le. fkeep%info%num_nodes) then
        call reduce_ind_and_get_ndep(fkeep, ind, nind, parent, nldep)
        ndep(i) = nldep
      end if
    end do

  end subroutine getSolveNDep



  subroutine getSolveDep(fkeep, node, ind, nind, dep, ndep)
    type(spllt_fkeep), intent(in) :: fkeep
    integer, intent(in)           :: node
    integer, intent(inout)        :: ind(:)
    integer, intent(inout)        :: nind
    integer, intent(inout)        :: dep(:)
    integer, intent(in)           :: ndep(:)

    integer :: i
    integer :: parent, nparent

    parent = node
    nparent = 0

    do while(parent .lt. fkeep%info%num_nodes)
      nparent = nparent + 1
      parent = fkeep%nodes(parent)%parent
    end do

    parent = node

    do i = 1, nparent
      parent = fkeep%nodes(parent)%parent

      if(parent .le. fkeep%info%num_nodes) then
        call reduce_ind_and_get_dep(fkeep, ind, nind, parent, &
          dep(ndep(i) : ndep(i + 1) - 1))
      end if

    end do

  end subroutine getSolveDep



  subroutine fwd_update_dependency(fkeep, blk, dep)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)     :: fkeep
    integer, intent(in)               :: blk    ! Index of block 
!   integer, pointer,     intent(out) :: dep(:) !List of dependencies
    integer, allocatable, intent(out) :: dep(:) !List of dependencies
    
    integer :: previous_dblk
    integer :: last_previous_dblk
    integer :: diff_bcol
    integer :: diff_previous_bcol
    integer :: node
    integer :: blk_sa, blk_en
    integer :: bcol_blk_sa, bcol
    integer :: blkm, nb, local_blk, offset
    integer :: child_node, i
    integer :: nchild, nnode_child, nind
    integer :: cur_node
    integer :: new_method
    integer, allocatable :: nblk_child(:) ! Remove by using a workspace
    integer, allocatable :: node_child(:) ! Remove by using a workspace
    integer, allocatable :: node_index(:) ! Remove by using a workspace
    integer, allocatable :: node_child_bis(:) ! Remove by using a workspace

    node    = fkeep%bc(blk)%node
    blk_sa  = fkeep%nodes(node)%blk_sa
    blk_en  = fkeep%nodes(node)%blk_en
    bcol_blk_sa = fkeep%bc(blk_sa)%bcol
    bcol = fkeep%bc(blk)%bcol

    nb      = fkeep%nodes(node)%nb    
    blkm    = fkeep%bc(blk)%blkm
    local_blk = blk - fkeep%bc(blk)%dblk
    offset  = nb * ( local_blk + (bcol - bcol_blk_sa)) + 1
    
    new_method = 1

!   print '(a, i3, a, i3)', "Initial dep of blk ", blk,        &
!     " is ", fkeep%bc(blk)%dep_initial
!   print *, "Trying to print from ", offset, " to ", offset + blkm - 1
!   call print_iarray("block row", blkm,                &
!     fkeep%nodes(node)%index(offset: offset + blkm - 1))

!   if(blk .eq. blk_sa) then
!     print '(a, i3, a, i3, a, i3)', "In node ", node, " blk_sa = ",    &
!       blk_sa, " and blk_en ", blk_en
!     print *, "Node index ", fkeep%nodes(node)%index
!     print *, "Children nodes are ", fkeep%nodes(node)%child
!   end if
!       call get_update_dep(fkeep, child_node, node, pos)
!       print *, "This correponds in index of the child node as "
!       call print_iarray("rows in child_node", size(pos),    &
!         fkeep%nodes(node)%index(pos))
!       deallocate(pos)

!   print *, "Get dependencies of blk ", blk

    !If blk belongs to the first block column of the node,
    ! then we consider the node dependencies of the block
    if(fkeep%nodes(node)%blk_sa .eq. fkeep%bc(blk)%dblk) then

      nchild  = node - fkeep%nodes(node)%least_desc
      allocate(nblk_child(nchild + 1), node_child(nchild))
      nblk_child = zero

      nnode_child = 1
      cur_node    = node

      do i = 1, nchild
        node_child(nnode_child : nnode_child + fkeep%nodes(cur_node)%nchild - 1) = fkeep%nodes(cur_node)%child
        nnode_child = nnode_child + fkeep%nodes(cur_node)%nchild
        cur_node = node_child(i)
      end do

      allocate(node_index(size(fkeep%nodes(node)%index)))
      node_index = fkeep%nodes(node)%index

      node_index  = fkeep%nodes(node)%index
      if(new_method .eq. 1) then
        if(nchild .gt. 0) then

          allocate(node_child_bis(nchild + 1))
          node_child_bis    = zero
          node_child_bis(1) = one

          nind  = size(node_index(offset : offset + blkm - 1))
!         print *, "=============== GET number of update for block ", blk

          call getUpdateNDep(fkeep, node, node_index(offset: offset + blkm - 1),&
            nind, node_child_bis(2 : nchild+1))
!         call print_iarray("=====> #dep in NODE_CHILD ", nchild, &
!           node_child_bis(2 : nchild + 1), 1)

          !Restore the array
          node_index = fkeep%nodes(node)%index

          !PostTreatment of ndep to compute the accsum
          do i = 1, nchild
            node_child_bis(i + 1) = node_child_bis(i) + node_child_bis(i + 1)
          end do
!         call print_iarray("PostTreated node_child_bis ", nchild + 1, &
!           node_child_bis, 1)

          if(node_child_bis(nchild + 1) .gt. 1) then
            allocate(dep(node_child_bis(nchild + 1) - 1))

            call getUpdateDep(fkeep, node, node_index(offset : offset + blkm - 1),&
              nind, dep, node_child_bis)
!           print *, "=====> dep of ", blk, " in NODE_CHILD are ", dep

!           deallocate(dep)
          else
            allocate(dep(1))
            dep(1) = blk
          end if
        else
          allocate(dep(1))
          dep(1) = blk
        end if
      else
!     node_index = fkeep%nodes(node)%index

        do i = 1, nchild
  !       child_node = fkeep%nodes(fkeep%nodes(node)%least_desc + i - 1)%num
          child_node = node_child(i)

          call get_update_nblk(fkeep, child_node,                     &
  !         fkeep%nodes(node)%index(offset : offset + blkm - 1),       &
            node_index(offset : offset + blkm - 1),                    &
            nblk_child(i + 1))

!         print *, nblk_child(i + 1), " found in ", child_node

          if(nblk_child(i+1) .gt. 0) then
            nblk_child(i+1) = nblk_child(i+1) + nblk_child(i)
          else
            nblk_child(i+1) = nblk_child(i)
          end if
        end do

!       call print_iarray("nblk_child", nchild + 1, nblk_child)
        node_index = fkeep%nodes(node)%index

        if(nblk_child(nchild+1) .gt. 0) then
          allocate(dep(nblk_child(nchild + 1)))

          do i = 1, nchild
            if(nblk_child(i + 1) .ne. nblk_child(i)) then
  !           child_node = fkeep%nodes(fkeep%nodes(node)%least_desc + i - 1)%num
              child_node = node_child(i)
              call get_update_dep_blk(fkeep, child_node,            &
  !             fkeep%nodes(node)%index(offset : offset + blkm - 1), &
                node_index(offset : offset + blkm - 1), &
                dep(nblk_child(i) + 1 : nblk_child(i+1)))
            end if
          end do

!         print *, "For blk ", blk
!         call print_iarray("Node dependencies found ", size(dep), dep)
        else

          allocate(dep(1))
          dep(1) = blk

  !       call print_iarray("Simulated dependencies ",  &
  !         size(dep), dep)

        end if

        deallocate(nblk_child)
      end if
    else

      allocate(dep(1))

!     if(bcol .gt. bcol_blk_sa) then
      previous_dblk       = fkeep%bc(fkeep%bc(blk)%dblk - 1)%dblk
      last_previous_dblk  = fkeep%bc(previous_dblk)%last_blk
      diff_bcol           = blk - fkeep%bc(blk)%dblk
      diff_previous_bcol  = last_previous_dblk - previous_dblk

      if(last_previous_dblk .gt. (previous_dblk + diff_bcol)) then
        dep(1) = previous_dblk + diff_bcol + 1
      end if
!     end if
!     call print_iarray("Local dependencies ", size(dep), dep)
    end if

!   call print_iarray("Returned dep from fwd_update_dependency",  &
!     size(dep), dep)

  end subroutine fwd_update_dependency



  integer function bwd_update_dependency(fkeep, blk)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)   :: fkeep
    integer, intent(in)             :: blk  ! Index of block 
    
    bwd_update_dependency = blk

    if(fkeep%bc(blk)%last_blk .ne. blk) then
      bwd_update_dependency = blk + 1
    end if

  end function bwd_update_dependency



  integer function fwd_solve_dependency(fkeep, blk)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)   :: fkeep
    integer, intent(in)             :: blk  ! Index of block 

    if(fkeep%bc(blk)%dblk .ne. blk) then
      fwd_solve_dependency = fkeep%bc(blk)%dblk
    else
      fwd_solve_dependency = blk
    end if

  end function fwd_solve_dependency



  subroutine bwd_solve_dependency(fkeep, blk, dep)
    use spllt_data_mod

    type(spllt_fkeep), intent(in)     :: fkeep
    integer, intent(in)               :: blk    ! Index of block 
!   integer, pointer, intent(out)     :: dep(:) !List of dependencies
    integer, allocatable, intent(out) :: dep(:) !List of dependencies

    integer :: offset
    integer :: node, parent, max_node, nparent
    integer :: nb, blkm, blk_sa, blk_en
    integer :: bcol_blk_en, bcol_blk_sa, bcol, nbcol
    integer :: lblk, nlblk, dblk
    integer, allocatable :: node_parent(:) ! Remove by using a workspace
    integer, allocatable :: node_index(:) ! Remove by using a workspace

    node      = fkeep%bc(blk)%node
    parent    = fkeep%nodes(node)%parent
    max_node  = fkeep%info%num_nodes
    nparent   = max_node - node

    blk_sa      = fkeep%nodes(node)%blk_sa
    blk_en      = fkeep%nodes(node)%blk_en
    bcol_blk_sa = fkeep%bc(blk_sa)%bcol
    bcol_blk_en = fkeep%bc(blk_en)%bcol
    bcol        = fkeep%bc(blk)%bcol
    nbcol       = bcol_blk_en - bcol_blk_sa + 1

    nb      = fkeep%nodes(node)%nb    
    blkm    = fkeep%bc(blk)%blkm
    dblk    = fkeep%bc(blk)%dblk
    local_blk = blk - dblk
    offset  = nb * ( local_blk + (bcol - bcol_blk_sa)) + 1

    lblk = local_blk + bcol - bcol_blk_sa + 1
    nlblk = fkeep%bc(blk_sa)%last_blk - blk_sa + 1

    allocate(node_index(size(fkeep%nodes(node)%index)))
    node_index = fkeep%nodes(node)%index

    !Check if blk belongs to the last bcol of the node
    ! or blk is a diagonal block
    ! or blk belongs to L_{21}, 
    ! where the L_{21} is the off diagonal block of L such that 
    ! L = [L_{11} ; L_{21}], with L_{11} is the triangular diagonal block
    if(bcol .eq. bcol_blk_en  &
      .or. blk .eq. dblk      &
      .or. lblk .gt. nbcol) then
      if(nparent .gt. 0) then

        allocate(node_parent(nparent + 1))
        node_parent    = zero
        node_parent(1) = one

        nind  = size(node_index(offset : offset + blkm - 1))
!       print *, "=============== GET number of update for block ", blk

        call getSolveNDep(fkeep, node, node_index(offset: offset + blkm - 1),&
          nind, node_parent(2 : nparent+1))
       !call print_iarray("=====> #dep in NODE_parent ", nparent, &
       !  node_parent(2 : nparent + 1), 1)

        !Restore the array
        node_index  = fkeep%nodes(node)%index
        nind        = size(node_index(offset : offset + blkm - 1))

        !PostTreatment of ndep to compute the accsum
        do i = 1, nparent
          node_parent(i + 1) = node_parent(i) + node_parent(i + 1)
        end do
       !call print_iarray("PostTreated node_parent_bis ", nparent + 1, &
       !  node_parent, 1)

        if(node_parent(nparent + 1) .gt. 1) then
          allocate(dep(node_parent(nparent + 1) - 1))

          call getSolveDep(fkeep, node, node_index(offset : offset + blkm - 1),&
            nind, dep, node_parent)
!         print *, "=====> dep of ", blk, " in NODE_parent are ", dep

          !       deallocate(dep)
        else
          allocate(dep(1))
          dep(1) = blk
        end if
      else

        allocate(dep(1))

        if(bcol .lt. bcol_blk_en) then
!         print *, "lblk = ", lblk, "/", nlblk
          dep(1) = blk_sa + (lblk - 1) * ( nlblk + 1 - 0.5 * lblk )
!         print *, "dep =====> ", dep(1)
        else
          dep(1) = blk
        end if
      end if
    else

      allocate(dep(1))

      dep(1) = blk_sa + (lblk - 1) * ( nlblk + 1 - 0.5 * lblk )
      
    end if

  end subroutine bwd_solve_dependency



  function contain(fkeep, blk1, blk2) result(isIn)
    use spllt_data_mod
    
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

  

  subroutine spllt_compute_extra_row(fkeep, node)
    type(spllt_fkeep),  target, intent(inout) :: fkeep
    integer,            intent(in)            :: node

    integer               :: i, nchild, child
    integer               :: blk_sa, last_blk
    integer               :: nblk
    integer               :: nrow, nval
    logical, allocatable  :: is_present(:)
    integer, pointer      :: p_index(:)
    integer, pointer      :: p_child_index(:)
    integer, pointer      :: p_extra_row(:)

    ! Number of nodes to visit
    nchild    = fkeep%nodes(node)%nchild
    blk_sa    = fkeep%nodes(node)%blk_sa
    last_blk  = fkeep%bc(blk_sa)%last_blk
    nblk      = last_blk - blk_sa + 1
    nrow      = nblk * fkeep%nodes(node)%nb + fkeep%bc(last_blk)%blkm
    p_index   => fkeep%nodes(node)%index
    nrow      = size(p_index, 1)

    if(nchild .eq. 0) then
      allocate(fkeep%nodes(node)%extra_row(nrow))
      fkeep%nodes(node)%extra_row(:) = p_index(:)
      return
    end if

    allocate(is_present(nrow))
    is_present(:) = .true.

    nval = nrow
    do i = 1, nchild
      child = fkeep%nodes(node)%child(i)
      j = 1
      k = 1
      p_child_index => fkeep%nodes(child)%index
      nchild_row = size(p_child_index, 1)
      do while(j .le. nrow .and. k .le. nchild_row .and. nval .gt. 0)
        if(p_index(j) .lt. p_child_index(k)) then
          j = j + 1
        else if(p_index(j) .gt. p_child_index(k)) then
          k = k + 1
        else
          if(is_present(j)) then
            nval = nval - 1
            is_present(j) = .false.
          end if
          j = j + 1
          k = k + 1
        end if
      end do
    end do

    allocate(fkeep%nodes(node)%extra_row(nval))
    j = 1
    i = 1
    p_extra_row => fkeep%nodes(node)%extra_row
    do while(i .le. nval .and. j .le. nrow)
      if(is_present(j)) then
        p_extra_row(i) = p_index(j)
        i = i + 1
      end if
      j = j + 1
    end do
    deallocate(is_present)

  end subroutine spllt_compute_extra_row


  ! This routine changes the list of dependencies of the last block of the tree
  ! given in parameter
  subroutine spllt_update_blk_dep(fkeep, tree)
    use utils_mod
    type(spllt_fkeep),  intent(inout) :: fkeep
    type(spllt_tree_t), intent(in)    :: tree

    integer :: blk, i, nblk
    integer :: blk_sa, blk_en
    integer :: st, cpt
    integer :: dep_blk, ndep
    logical, allocatable :: is_dep(:)

    ! Considering the subtree given in parameter, the algorithm
    ! computes the dependencies of all its blocks,
    ! and replaces the dependencies of the last block, denoted blk_en,
    ! with the dependencies of all the blocks in the tree that are 
    ! greater than blk_en.


    nblk = fkeep%nodes(fkeep%info%num_nodes)%blk_en
    ndep = 0
    blk_sa = fkeep%nodes(tree%node_sa)%blk_sa
    blk_en = fkeep%nodes(tree%node_en)%blk_en

    ! Set blocks greater that blk_en as to not be a dependency
    allocate(is_dep(blk_en + 1 : nblk), stat = st)
    is_dep(:) = .false.

   !! Count the number of actual dependency blocks of blk_en
   !if(fkeep%bc(blk_en)%bwd_update_dep .ne. blk_en) then
   !  ndep = ndep + 1
   !end if

    ! For each block blk in the tree, for each dep_blk in the list of dep of 
    ! the current blk, if the dep_blk is greater than blk_en, i.e. it does not
    ! belongs to the tree, then this block is recorded and counted only once.
    do blk = blk_sa, blk_en
      do i = 1, size(fkeep%bc(blk)%bwd_solve_dep)
        dep_blk = fkeep%bc(blk)%bwd_solve_dep(i)
        if(dep_blk .gt. blk_en) then
          ndep = ndep + merge(0, 1, is_dep(dep_blk)) ! reverse to count it once
          is_dep(dep_blk) = .true.
        end if
      end do
    end do

    ! Rebuild the bwd_dep list of the last block of the tree
    deallocate(fkeep%bc(blk_en)%bwd_dep)
    allocate(fkeep%bc(blk_en)%bwd_dep(ndep))

    ! Copy the dependency blocks in ascending order
    cpt = 1
    do blk = blk_en + 1, nblk
      if(is_dep(blk)) then
        fkeep%bc(blk_en)%bwd_dep(cpt) = blk
        cpt = cpt + 1
      end if
    end do

   !if(fkeep%bc(blk_en)%bwd_update_dep .ne. blk_en) then
   !  fkeep%bc(blk_en)%bwd_dep(cpt) = fkeep%bc(blk_en)%bwd_update_dep
   !end if

    deallocate(is_dep)

  end subroutine spllt_update_blk_dep



  subroutine spllt_compute_rhs_block(fkeep, st)
    use spllt_data_mod
    implicit none
    type(spllt_fkeep), target,  intent(inout) :: fkeep
    integer,                    intent(out)   :: st

    integer :: nnode, node
    integer :: ndblk
    integer :: sa, en
    integer :: ncol
    integer :: i
    integer :: dblk, blk
    integer :: nb
    integer :: nlblk
    integer, pointer :: rhsPtr(:)

    nnode = fkeep%info%num_nodes
    ndblk = 0
    ! Count the number of diagonal block
    do node = 1, nnode
      sa    = fkeep%nodes(node)%sa
      en    = fkeep%nodes(node)%en
      ncol  = en - sa + 1
      nb    = fkeep%nodes(node)%nb
      nlblk = ncol/nb + merge(1, 0, mod(ncol, nb) .gt. 0) ! #local diag block
!!    print *, "Node ", node, "has", nlblk, "dblk"
      ndblk  = ndblk + nlblk
    end do
!!  print *, "Total number of dblk", ndblk

    fkeep%ndblk = ndblk
    allocate(fkeep%rhsPtr(ndblk + 1), stat=st)

    rhsPtr => fkeep%rhsPtr
    rhsPtr(1) = 0
    blk       = 0 ! #blk already treated

    ! Fill in rhsPtr array
    do node = 1, nnode
      sa    = fkeep%nodes(node)%sa
      en    = fkeep%nodes(node)%en
      ncol  = en - sa + 1
      nb    = fkeep%nodes(node)%nb
      nlblk = ncol/nb + merge(1, 0, mod(ncol, nb) .gt. 0) ! #local diag block
      dblk  = fkeep%nodes(node)%blk_sa
      do i = 1, nlblk
!!      print *, fkeep%bc(dblk)%blkn, "ncol of blk", dblk
        rhsPtr(blk + i + 1 ) = rhsPtr(blk + i) + fkeep%bc(dblk)%blkn
        dblk = fkeep%bc(dblk)%last_blk + 1
      end do
      blk = blk + nlblk
    end do

!!  print *, "Result of rhsPtr", rhsPtr

  end subroutine spllt_compute_rhs_block

! integer function bwd_solve_dependency(fkeep, blk)
!   use spllt_data_mod

!   type(spllt_fkeep), intent(in)   :: fkeep
!   integer, intent(in)             :: blk  ! Index of block 

!   integer :: next_dblk    
!   integer :: dist

!   dist = blk - fkeep%bc(blk)%dblk
!   next_dblk = fkeep%bc(blk)%last_blk + 1

!   if(dist .gt. 0) then
!     bwd_solve_dependency = next_dblk + dist - 1
!   else
!     bwd_solve_dependency = blk
!   end if

! end function bwd_solve_dependency
end module spllt_solve_dep_mod
