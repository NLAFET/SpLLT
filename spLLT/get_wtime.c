#include <stdio.h>
#include <sys/time.h>

double get_wtime() {
   struct timeval tv;
   gettimeofday(&tv,NULL);
   return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
