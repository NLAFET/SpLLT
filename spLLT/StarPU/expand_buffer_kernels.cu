#define EXPAND_BUFFER_NTX 8  // Number of threads x
#define EXPAND_BUFFER_NTY 8   // Number of threads y

void __global__ cu_expand_buffer(double *buffer, 
                                 int row_list, int col_list, int row_list_sz, int col_list_sz) {

   int j = (blockIdx.y * blockDim.y) + threadIdx.y;
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;

   if (j < col_list_sz) {
      int col = col_list[j]-1
      if (i < row_list_sz) {
         int row = row_list[i]-1
         
      }
   }
      
}

void spllt_cu_expand_buffer(double *buffer,
                            int row_list, int col_list, int row_list_sz, int col_list_sz, 
                            const cudaStream_t *stream) {
   
   dim3 threads(EXPAND_BUFFER_NTX, EXPAND_BUFFER_NTY);

   dim3 blocks(row_list_sz/threads.x, col_list_sz/threads.y);

   cu_expand_buffer<<<blocks, threads, 0, *stream>>>(row_list, col_list, row_list_sz, col_list_sz);
}
