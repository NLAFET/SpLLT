forall nodes snode in post-order
  call alloc(snode) ! allocate data structures
   
  call init(snode) ! initianlize node structure
end do

forall nodes snode in post-order
  ! factorize node
  do k=1..n in snode
    call factorize(blk(k,k)) ! factorize block
    ! solve block
    do i=k+1..m in snode
        call solve(blk(k,k), blk(i,k))
    end do
    ! update block
    do j=k+1..n in snode
      do i=k+1..m in snode
        call update(blk(j,k), blk(i,k), blk(i,j))
      end do
    end do

    ! update ancestor nodes
    forall ancestors(snode) anode 
      do j=k+1..p(anode) in snode
        do i=k+1..m in snode
           call update_btw(blk(j,k), blk(i,k),
                           a_blk(rmap(i), cmap(j)))
        end do
      end do
    end do

  end do
end do
