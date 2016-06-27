forall nodes snode in post-order
   call alloc(snode) ! allocate data structures
   
   call submit(init, snode:W) ! initianlize node structure
end do

forall nodes snode in post-order
  ! factorize node
  do k=1..n in snode
    call submit(factorize, snode:R, blk(k,k):RW) ! factorize block
    
    do i=k+1..m in snode
        call submit(solve, blk(k,k):R, blk(i,k):RW) ! perform solve
    end do

    do j=k+1..n in snode
       do i=k+1..m in snode
          call submit(update, blk(j,k):R, blk(i,k):R, blk(i,j):RW)
       end do
    end do

    ! update ancestor nodes
    forall ancestors(snode) anode 
      do j=k+1..p(anode) in snode
         do i=k+1..m in snode
            call submit(update_btw, blk(j,k):R, blk(i,k):R, a_blk(rmap(i), cmap(j)):RW)
         end do
      end do
    end do

  end do
end do
