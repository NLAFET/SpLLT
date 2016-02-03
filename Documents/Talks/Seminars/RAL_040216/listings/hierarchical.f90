forall outer panels o_p=1...o_n in f
   ! partition (outer) block column f(o_p) into
   ! i_n inner block columns f(o_p,1) .. f(o_p,i_n)
   call submit(partition, f(o_p):R, f(o_p,1):W, .., f(o_p,i_n):W)     |\alert{$\longleftarrow$ partition block}|

   forall inner panels i_p=1..i_n
      ! panel reduction of inner block column i_p
      call submit(_geqrt, f(i_p):RW)

      ! update (inner) column in_u with panel p
      forall inner blockcolumns i_u=i_p+1..i_n in f(o_p)
        call submit(_gemqrt, f(i_p):R, f(i_u):RW)
      end do

      ! update outer block column o_u with panel i_p
      forall outer blockcolumns o_u=o_p+1..o_n
        call submit(_gemqrt, f(i_p):R,f(o_u):RW)
      end do
   end do

   ! unpartition (outer) block column
   call submit(unpartition, f(o_p,1):R..f(o_p,i_n):R, f(o_p):W)       |\alert{$\longleftarrow$ unpartition block}|
end do
