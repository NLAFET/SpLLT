module spllt_solve_dep_mod

  use spllt_data_mod

contains

  !Return the position in child_node of the rows of node that are 
  ! in child_node
  subroutine get_update_dep(fkeep, child_node, node, pos)
    
    type(spllt_fkeep), target, intent(in) :: fkeep
    integer, intent(in)                   :: child_node
    integer, intent(in)                   :: node
    integer, allocatable, intent(out)     :: pos(:)

    integer :: node_index
    integer :: child_node_index
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

    integer :: child_node_index
    integer :: i, j, k
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

    integer :: child_node_index
    integer :: i, j, k
    integer :: nb
    integer, pointer :: p_child_node_index(:)
    integer :: cur_blk_dep, tmp
    integer :: lblk, nlblk, diff

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
    integer :: lblk, nlblk, diff
    integer, pointer :: p_child_blk_index(:)

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

    integer :: child, dep, nldep, offset, i, nsubind, nchild
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
    integer :: lblk, nlblk, diff, ndep

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
    integer :: lblk, nlblk, diff

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

end module spllt_solve_dep_mod
