forall nodes snode in post-order
   ! allocate data structures
   call alloc(snode)
   ! initianlize node structure
   call submit(init, snode:W) 
end do

forall nodes snode in post-order

  ! factorize node
  do k=1..n in snode
    call submit(factorize, snode:R, blk(k,k):RW)
    
    do i=k+1..m in snode
        call submit(solve, blk(k,k):R, blk(i,k):RW)
    end do

    do j=k+1..n in snode
       do i=k+1..m in snode
          call submit(update, blk(j,k):R, blk(i,k):R, blk(i,j):RW)
       end do
    end do

    forall ancestors(snode) anode 
      do j=k+1..p(anode) in snode
         do i=k+1..q(anode) in snode
            call submit(update_between, blk(j,k):R, blk(i,k):R, a_blk(rmap(i), cmap(j)):RW)
         end do
      end do
    end do

  end do
end do
