include ../make.inc

DEPS = Makefile ../make.inc ../common.hpp

all: bench_cuBLAS_syrk.exe bench_MKL_syrk.exe

bench_cuBLAS_syrk.exe: bench_cuBLAS_syrk.cu $(DEPS)
	$(NVCC) -I$(CUDADIR_INC) $(CPPFLAGS) $(CXXFLAGS) bench_cuBLAS_syrk.cu -o $@ -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" $(LDFLAGS) $(CUDA_LIBS) $(LIBS)

bench_MKL_syrk.exe: bench_MKL_syrk.cu $(DEPS)
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) bench_MKL_syrk.cu -o $@ $(LDFLAGS) $(LIBS)

clean:
	-rm -f bench_MKL_syrk.exe
	-rm -f bench_cuBLAS_syrk.exe
