forall nodes snode in post-order
   ! allocate data structures
   call submit(alloc, snode:RW)
   ! initianlize node structure
   call submit(init, snode:RW) 
end do

forall nodes snode in post-order

  ! factorize node
  ! call submit(factorize, snode:RW)
  do k=1..n in snode
    ! factorize diagonal block
    call submit(factorize, blk(k,k):RW)
    
    do i=k+1..m in snode
       ! perform triangular solve w.r.t diag block
        call submit(solve, blk(k,k):R, blk(i,k):Rw)
    end do

    do j=k+1..n
       do i=k+1..m
          call submit(update, blk(j,k):R, blk(i,k):R, blk(i,j):RW)
       end do
    end do

    forall ancestors(snode) anode
      ! udpate ancestor nodes
      call submit(update, snode:R, anode:RW)
    end do

  end do
end do
