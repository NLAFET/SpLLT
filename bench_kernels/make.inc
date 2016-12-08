MAGMADIR ?= /usr/local/magma
MAGMADIR_LIB = $(MAGMADIR)/lib
MAGMADIR_INC = $(MAGMADIR)/include
MAGMA_LIBS = -lmagma

CUDADIR ?= /usr/local/cuda
CUDADIR_LIB = $(CUDADIR)/lib64
CUDADIR_INC = $(CUDADIR)/include
CUDA_LIBS = -lcusolver -lcublas

LAPACKDIR ?= $(MKLROOT)
LAPACKDIR_LIB = $(LAPACKDIR)/lib/intel64_lin
LAPACKDIR_INC = $(LAPACKDIR)/include
LAPACK_LIBS = -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core

NVCC ?= nvcc
ifdef NDEBUG
CPPFLAGS = -DNDEBUG -I..
CXXFLAGS = -O$(NDEBUG)
else
CPPFLAGS = -I..
CXXFLAGS = -O0 -g -G -src-in-ptx
endif
CXXFLAGS += -std=c++11 -arch=sm_35 -res-usage -Xcompiler="-fopenmp"
CPPFLAGS += -DADD_ -I$(MAGMADIR_INC) -I$(CUDADIR_INC) -I$(LAPACKDIR_INC)
LIBS = $(MAGMA_LIBS) $(CUDA_LIBS) $(LAPACK_LIBS) -lgomp -lrt -lpthread -lm
LDFLAGS = -L$(MAGMADIR_LIB) -Xlinker="-rpath=$(MAGMADIR_LIB)" -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" -L$(LAPACKDIR_LIB) -Xlinker="-rpath=$(LAPACKDIR_LIB)"
