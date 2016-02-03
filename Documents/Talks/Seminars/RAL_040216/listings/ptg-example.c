N     [type = int]

task_f(i) /* Task name */

i = 1..N-1 /* Execution space declaration for parameter i */

: x(i) /* Task must be executed on the node where x(i) is stored */

  /* Task reads x(i) from memory ... */
  RW X <- x(i)     
  /* ... and sends it to task_g(i) */
  -> X task_g(i)
BODY
  X = f(X)   /* Code executed by the task */
END

task_g(i) /* Task name */

i = 1..N-1 /* Execution space declaration for parameter i */

: y(i) /* Task must be executed on the node where x(i) is stored */


  /* Task reads x(i) fron task_f(i)... */
  R X  <- X task_f(i)                        
  /* ... y(i-1) from task_g(i-1)... */
  R Y1 <- (i > 1) ? Y2 task_g(i-1) : y(i-1)
  /* ... and sends y(i) to task_g(i+1) */
  W Y2 -> (i < N-1) ? Y2 task_g(i+1)       

BODY
  Y2 = h(X, Y1) /* Code executed by the task */
END
