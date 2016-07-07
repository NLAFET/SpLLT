forall nodes snode in post-order
   call alloc(snode) ! allocate data structures
   
   call init(snode) ! initianlize node structure
end do

forall nodes snode in post-order
  ! factorize node
   call factorize(snode)

   ! update ancestor nodes
   forall ancestors(snode) anode 
     call update_btw(snode, anode)
   end do

  end do
end do















!
