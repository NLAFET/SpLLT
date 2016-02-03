do f=1, nfronts                               ! in postorder
  call activate(f)                            ! activate front
  call submit(init, f:RW)                     ! init front

  do c=1, f%nchildren                         ! for all the children of f
    do i=1,c%m
      do j=1,c%n
        call submit(assemble, c(i,j):R, f:RW) ! assemble block(i,j) of c
      end do
    end do
    call submit(deactivate, c:RW)             ! Deactivate child
  end do

  ca_facto: do k=1, min(f%m,f%n)
    do s=0, log2(f%m-k+1)
      do i = k, f%n, 2**s
        if(s.eq.0) then
          call submit(_geqrt, f(i,k):RW)
          do j=k+1, f%n
            call submit(_gemqrt, f(i,k):R, f(i,j):RW)
          end do
        else
          l = i+2**(s-1)
          call submit(_tpqrt, f(i,k):RW, f(l,k):RW)
          do j=k+1, front%n
            call submit(_tpmqrt, f(l,k):R, f(i,j):RW, f(l,j):RW)
          end do
        end if
      end do
    end do
  end do ca_facto
end do
call wait_tasks_completion()                  ! wait for the tasks to be executed
