forall nodes snode in post-order
   ! allocate data structures
   call alloc(snode)
   ! initianlize node structure
   call init(snode) 
end do

forall nodes snode in post-order

  ! factorize node
  do k=1..n in snode
    call factorize(blk(k,k))
    
    do i=k+1..m in snode
        call solve(blk(k,k), blk(i,k))
    end do

    do j=k+1..n in snode
       do i=k+1..m in snode
          call update(blk(j,k), blk(i,k), blk(i,j))
       end do
    end do
    
    forall ancestors(snode) anode 
      do j=k+1..p in snode
         do i=j..q in snode
            call update_btw(blk(j,k), blk(i,k), a_blk(rmap(i), cmap(j)))
         end do
      end do
    end do

  end do
end do
