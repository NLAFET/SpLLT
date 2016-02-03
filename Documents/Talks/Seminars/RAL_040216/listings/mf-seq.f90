forall fronts f in topological order

   ! compute front structure 
   call activate(f)
   ! allocate and initialize front
   call init(f)

   ! front assembly
   forall children c of f
      call assemble(c, f)
      ! Deactivate child
      call deactivate(c)
   end do

   ! front factorization
   call factorize(f)
end do
