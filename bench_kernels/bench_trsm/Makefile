include ../make.inc

DEPS = Makefile ../make.inc ../common.hpp

all: bench_MKL_trsm.exe bench_cuBLAS_trsm.exe

bench_MKL_trsm.exe: bench_MKL_trsm.cu $(DEPS)
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) bench_MKL_trsm.cu -o $@ $(LDFLAGS) $(LIBS)

bench_cuBLAS_trsm.exe: bench_cuBLAS_trsm.cu $(DEPS)
	$(NVCC) -I$(CUDADIR_INC) $(CPPFLAGS) $(CXXFLAGS) bench_cuBLAS_trsm.cu -o $@ -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" $(LDFLAGS) $(CUDA_LIBS) $(LIBS)

bench_ESSL_trsm.exe: bench_ESSL_trsm.cpp $(DEPS)
	xlC_r -DNDEBUG -I.. -std=gnu++11 -O -qfloat=nans:spnans:subnormals -qmaxmem=-1 -qsmp=omp -qarch=pwr8 -qtune=pwr8:smt8 bench_ESSL_trsm.cpp -o $@ -L/usr/lib64 -lesslsmp -lessl $(HOME)/lapack/libtmglib.a -lesslsmp -lessl $(HOME)/lapack/liblapack.a -lesslsmp -lessl $(HOME)/lapack/librefblas.a -L/opt/ibm/xlf/15.1.3/lib -lxlf90_r -lxlfmath

bench_OpenBLAS_trsm.exe: bench_OpenBLAS_trsm.cpp $(DEPS)
	g++ -DNDEBUG -I.. -I$(NLAFETPREFIX)/OpenBLAS/include -std=gnu++11 -O3 -fopenmp -march=native bench_OpenBLAS_trsm.cpp -o $@ -L$(NLAFETPREFIX)/OpenBLAS/lib -lopenblas_omp -lm

clean:
	-rm -f bench_OpenBLAS_trsm.exe
	-rm -f bench_ESSL_trsm.exe
	-rm -f bench_cuBLAS_trsm.exe
	-rm -f bench_MKL_trsm.exe
