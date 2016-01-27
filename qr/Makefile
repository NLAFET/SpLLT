CXX = g++
#CXXFLAGS = -std=c++11 -Wall -g -O2 -fsanitize=address -march=native -mavx2
CXXFLAGS = -std=c++11 -Wall -g -O2 -march=native# -mavx2 -mfma
#FC = gfortran-5
#FCFLAGS = -Wall -g -O2 -march=native

#CXX = icc
#CXXFLAGS = -std=c++11 -g -ip -Wall -march=native
FC = ifort
FCFLAGS = -g -O3 -ip -march=native -opt-report

LIBS = -llapack -lblas

# VTUNE_DIR = /opt/intel/vtune_amplifier_xe
# CXXFLAGS += -I $(VTUNE_DIR)/include
# LIBS += -L$(VTUNE_DIR)/lib64 -littnotify -lpthread -ldl

all: bench

bench: qr.o bench.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	-rm bench
	-rm *.o

qr.o: qr.hh
bench.o: qr.hh problem.hh cpu_wrappers.hh
