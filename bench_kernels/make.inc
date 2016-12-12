MAGMADIR ?= /usr/local/magma
MAGMADIR_INC = $(MAGMADIR)/include
MAGMADIR_LIB = $(MAGMADIR)/lib
MAGMA_LIBS = -lmagma

CUDADIR ?= /usr/local/cuda
CUDADIR_INC = $(CUDADIR)/include
CUDADIR_LIB = $(CUDADIR)/lib64
CUDA_LIBS = -lcublas

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
CPPFLAGS += -I$(LAPACKDIR_INC)
CXXFLAGS += -std=c++11 -arch=sm_35 -res-usage -Xcompiler="-fopenmp"
LDFLAGS = -L$(LAPACKDIR_LIB) -Xlinker="-rpath=$(LAPACKDIR_LIB)"
LIBS = $(LAPACK_LIBS) -lgomp -lrt -lpthread -lm

ifdef USE_COMPLEX
CPPFLAGS += -DUSE_COMPLEX=$(USE_COMPLEX)
endif

ifdef USE_FLOAT
CPPFLAGS += -DUSE_FLOAT=$(USE_FLOAT)
endif
