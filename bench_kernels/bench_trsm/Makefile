include ../make.inc

DEPS = Makefile ../make.inc ../common.hpp

all: bench_cuBLAS_trsm.exe bench_MKL_trsm.exe

bench_cuBLAS_trsm.exe: bench_cuBLAS_trsm.cu $(DEPS)
	$(NVCC) -I$(CUDADIR_INC) $(CPPFLAGS) $(CXXFLAGS) bench_cuBLAS_trsm.cu -o $@ -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" $(LDFLAGS) $(CUDA_LIBS) $(LIBS)

bench_MKL_trsm.exe: bench_MKL_trsm.cu $(DEPS)
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) bench_MKL_trsm.cu -o $@ $(LDFLAGS) $(LIBS)

clean:
	-rm -f bench_MKL_trsm.exe
	-rm -f bench_cuBLAS_trsm.exe
