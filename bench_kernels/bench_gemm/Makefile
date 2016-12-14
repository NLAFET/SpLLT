include ../make.inc

DEPS = Makefile ../make.inc ../common.hpp

all: bench_cuBLAS_gemm.exe bench_MKL_gemm.exe

bench_cuBLAS_gemm.exe: bench_cuBLAS_gemm.cu $(DEPS)
	$(NVCC) -I$(CUDADIR_INC) $(CPPFLAGS) $(CXXFLAGS) bench_cuBLAS_gemm.cu -o $@ -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" $(LDFLAGS) $(CUDA_LIBS) $(LIBS)

bench_MKL_gemm.exe: bench_MKL_gemm.cu $(DEPS)
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) bench_MKL_gemm.cu -o $@ $(LDFLAGS) $(LIBS)

clean:
	-rm -f bench_MKL_gemm.exe
	-rm -f bench_cuBLAS_gemm.exe
