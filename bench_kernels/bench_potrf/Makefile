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

bench_ESSL_potrf.exe: bench_ESSL_potrf.cpp $(DEPS)
	xlC_r -DNDEBUG -I.. -std=gnu++11 -O -qfloat=nans:spnans:subnormals -qmaxmem=-1 -qsmp=omp -qarch=pwr8 -qtune=pwr8:smt8 bench_ESSL_potrf.cpp -o $@ -L/usr/lib64 -lesslsmp -lessl $(HOME)/lapack/libtmglib.a -lesslsmp -lessl $(HOME)/lapack/liblapack.a -lesslsmp -lessl $(HOME)/lapack/librefblas.a -L/opt/ibm/xlf/15.1.3/lib -lxlf90_r -lxlfmath

bench_OpenBLAS_potrf.exe: bench_OpenBLAS_potrf.cpp $(DEPS)
	g++ -DNDEBUG -I.. -I$(NLAFETPREFIX)/OpenBLAS/include -std=gnu++11 -O3 -fopenmp -march=native bench_OpenBLAS_potrf.cpp -o $@ -L$(NLAFETPREFIX)/OpenBLAS/lib -lopenblas_omp -lm

clean:
	-rm -f bench_OpenBLAS_potrf.exe
	-rm -f bench_ESSL_potrf.exe
	-rm -f bench_MKL_potrf.exe
	-rm -f bench_MAGMA_potf2_gpu.exe
	-rm -f bench_MAGMA_potrf_gpu.exe
	-rm -f bench_cuSOLVER_potrf.exe
