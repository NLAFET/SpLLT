include ../make.inc

DEPS = Makefile ../make.inc ../common.hpp

all: bench_cuSOLVER_potrf.exe bench_MAGMA_potrf_gpu.exe bench_MAGMA_potf2_gpu.exe bench_MKL_potrf.exe

bench_cuSOLVER_potrf.exe: bench_cuSOLVER_potrf.cu $(DEPS)
	$(NVCC) -I$(CUDADIR_INC) $(CPPFLAGS) $(CXXFLAGS) bench_cuSOLVER_potrf.cu -o $@ -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" $(LDFLAGS) -lcusolver $(CUDA_LIBS) $(LIBS)

bench_MAGMA_potrf_gpu.exe: bench_MAGMA_potrf_gpu.cu $(DEPS)
	$(NVCC) -DADD_ -I$(MAGMADIR_INC) -I$(CUDADIR_INC) $(CPPFLAGS) $(CXXFLAGS) bench_MAGMA_potrf_gpu.cu -o $@ -L$(MAGMADIR_LIB) -Xlinker="-rpath=$(MAGMADIR_LIB)" -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" $(LDFLAGS) $(MAGMA_LIBS) -lcusparse $(CUDA_LIBS) $(LIBS)

bench_MAGMA_potf2_gpu.exe: bench_MAGMA_potf2_gpu.cu $(DEPS)
	$(NVCC) -DADD_ -I$(MAGMADIR_INC) -I$(CUDADIR_INC) $(CPPFLAGS) $(CXXFLAGS) bench_MAGMA_potf2_gpu.cu -o $@ -L$(MAGMADIR_LIB) -Xlinker="-rpath=$(MAGMADIR_LIB)" -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" $(LDFLAGS) $(MAGMA_LIBS) -lcusparse $(CUDA_LIBS) $(LIBS)

bench_MKL_potrf.exe: bench_MKL_potrf.cu $(DEPS)
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) bench_MKL_potrf.cu -o $@ $(LDFLAGS) $(LIBS)

clean:
	-rm -f bench_MKL_potrf.exe
	-rm -f bench_MAGMA_potf2_gpu.exe
	-rm -f bench_MAGMA_potrf_gpu.exe
	-rm -f bench_cuSOLVER_potrf.exe
