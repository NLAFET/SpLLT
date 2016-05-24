forall nodes snode in post-order
   ! allocate data structures
   call alloc(snode)
   ! initianlize node structure
   call submit(init, snode:W) 
end do

forall nodes snode in post-order

  ! factorize node
  do k=1..n in snode
    ! factorize diagonal block
    call submit(factorize, snode:R, blk(k,k):RW)
    
    do i=k+1..m in snode
       ! perform triangular solve w.r.t diag block
        call submit(solve, blk(k,k):R, blk(i,k):RW)
    end do

    do j=k+1..n
       do i=k+1..m
          call submit(update, blk(j,k):R, blk(i,k):R, blk(i,j):RW)
       end do
    end do

    forall ancestors(snode) anode
      ! udpate ancestor nodes
      call submit(update_between, snode:R, anode:RW)
    end do

  end do
end do
