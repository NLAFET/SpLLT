MAGMADIR ?= /usr/local/magma
MAGMADIR_INC = $(MAGMADIR)/include
MAGMADIR_LIB = $(MAGMADIR)/lib
MAGMA_LIBS = -lmagma

CUDADIR ?= /usr/local/cuda
CUDADIR_INC = $(CUDADIR)/include
CUDADIR_LIB = $(CUDADIR)/lib64
CUDA_LIBS = -lcusolver -lcublas

LAPACKDIR ?= $(MKLROOT)
LAPACKDIR_INC = $(LAPACKDIR)/include
LAPACKDIR_LIB = $(LAPACKDIR)/lib/intel64_lin
LAPACK_LIBS = -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core

NVCC ?= nvcc
ifdef NDEBUG
CPPFLAGS = -DNDEBUG -I..
CXXFLAGS = -O$(NDEBUG)
else # debug
CPPFLAGS = -I..
CXXFLAGS = -O0 -g -G -src-in-ptx
endif
CPPFLAGS += -DADD_ -I$(MAGMADIR_INC) -I$(CUDADIR_INC) -I$(LAPACKDIR_INC)
CXXFLAGS += -std=c++11 -arch=sm_35 -res-usage -Xcompiler="-fopenmp"
LDFLAGS = -L$(MAGMADIR_LIB) -Xlinker="-rpath=$(MAGMADIR_LIB)" -L$(CUDADIR_LIB) -Xlinker="-rpath=$(CUDADIR_LIB)" -L$(LAPACKDIR_LIB) -Xlinker="-rpath=$(LAPACKDIR_LIB)"
LIBS = $(MAGMA_LIBS) $(CUDA_LIBS) $(LAPACK_LIBS) -lgomp -lrt -lpthread -lm

ifdef USE_COMPLEX
CPPFLAGS += -DUSE_COMPLEX=$(USE_COMPLEX)
endif

ifdef USE_FLOAT
CPPFLAGS += -DUSE_FLOAT=$(USE_FLOAT)
endif
