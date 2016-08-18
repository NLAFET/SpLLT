#define EXPAND_BUFFER_NTX 8  // Number of threads x
#define EXPAND_BUFFER_NTY 8   // Number of threads y

void __global__ cu_expand_buffer1(int n, double *a, double *buffer, 
                                  int *row_list, int *col_list, int rls, int cls,
                                  int ndiag) {

   
   int i, j, imax;

   for (i = 0; i < rls; i++) {
      int row = row_list[i]-1;
      
      imax = cls-1;
      if (i < ndiag) imax = i;
      for (j = 0; j <= imax; j++) {
         int col = col_list[j]-1;
         a[row*n + col] += buffer[i*cls + j];         
      }
   }
}

void __global__ cu_expand_buffer(int n, double *a, double *buffer, 
                                 int *row_list, int rls, int *col_list, int cls,
                                 int ndiag) {

   int i = (blockIdx.x * blockDim.x) + threadIdx.x;   
   int j = (blockIdx.y * blockDim.y) + threadIdx.y;

   int imax; 

   if (i < rls) {
      int row = row_list[i]-1;
      imax = cls-1;
      if (i < ndiag) imax = i;
      if (j <= imax) {
         int col = col_list[j]-1;
         a[row*n + col] += buffer[i*cls + j];
      }
   }
}

extern "C" {
void spllt_cu_expand_buffer(int n, double *a, double *buffer,
                            int *row_list, int rls, int *col_list, int cls,
                            int ndiag,
                            const cudaStream_t stream) {
   
   // dim3 threads(EXPAND_BUFFER_NTX, EXPAND_BUFFER_NTY);
   // dim3 blocks(rls/EXPAND_BUFFER_NTX, cls/EXPAND_BUFFER_NTY);

   dim3 threads(EXPAND_BUFFER_NTX, EXPAND_BUFFER_NTY);
   dim3 blocks((rls-1)/EXPAND_BUFFER_NTX+1, (cls-1)/EXPAND_BUFFER_NTY+1);

   cu_expand_buffer<<<blocks, threads, 0, stream>>>(n, a, buffer, row_list, rls, col_list, cls, ndiag);

   // cu_expand_buffer1<<<1, 1, 0, stream>>>(n, a, buffer, row_list, col_list, rls, cls, ndiag);
}
}
