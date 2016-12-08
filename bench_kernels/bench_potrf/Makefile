include ../make.inc

DEPS = Makefile ../make.inc ../common.hpp

all: bench_potrf.exe

bench_potrf.exe: bench_potrf.cu $(DEPS)
	$(NVCC) $(CPPFLAGS) $(CXXFLAGS) bench_potrf.cu -o $@ $(LDFLAGS) $(LIBS)

clean:
	-rm -f bench_potrf.exe
