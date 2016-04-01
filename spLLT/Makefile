include Make.inc

#Select compiler and options from below
############################
CC = gcc
FC = gfortran

CFLAGS = -Wall -pedantic -g -O2
#FFLAGS = -Wall -pedantic -fbounds-check -march=native -g -fbacktrace -fopenmp -ffpe-trap=invalid,zero,overflow
#FFLAGS = -Wall -pedantic -fbounds-check -g -fbacktrace -fopenmp -ffpe-trap=invalid,zero,overflow
#FFLAGS = -Wall -pedantic -fbounds-check -march=native -g -Wunused -fbacktrace -fopenmp
#FFLAGS += -fprofile-arcs -ftest-coverage
FFLAGS = -Wall -pedantic -g -O2 -fopenmp
############################
#FC = ifort
#FFLAGS = -g -C -static -openmp
############################

# Things you might add to DEFS:

# SPLLT_STARPU_NOSUB: wait for all the tasks to be submitted to StarPU
# before starting the factorization

#DEFS =
DEFS = -DSPLLT_USE_STARPU # -DSPLLT_STARPU_NOSUB

CFLAGS += $(DEFS)
FFLAGS += $(DEFS)

CFLAGS += $(ISTARPU) $(ICUDA) $(IHWLOC)
FFLAGS += $(ISPRAL) 

LIBS = $(LBLAS) $(LLAPACK) $(LSPRAL) $(LSTARPU) $(LMETIS)

all: run_ma87 spllt_starpu_test

# HSL deps
DTDEPS := $(HSLPACKDIR)/fa14/fa14d.f
DTDEPS += $(HSLPACKDIR)/mc34/mc34d.f
DTDEPS += $(HSLPACKDIR)/ym11/ym11d.f
DTDEPS += $(HSLPACKDIR)/mc47/mc47d.f
DTDEPS += $(HSLPACKDIR)/mc59/mc59d.f
DTDEPS += $(HSLPACKDIR)/mc54/mc54d.f
DTDEPS += $(HSLPACKDIR)/mc56/mc56d.f
DTDEPS += $(HSLPACKDIR)/mc56/mc56i.f

# HSL f90 deps
COMMON90 := $(HSLPACKDIR)/hsl_mc78/hsl_mc78i.f90
DDEPS90  := $(HSLPACKDIR)/hsl_mc34/hsl_mc34d.f90
DTDEPS90 := $(HSLPACKDIR)/hsl_mc34/hsl_mc34d.f90
DTDEPS90 += $(HSLPACKDIR)/hsl_fa14/hsl_fa14d.f90
DTDEPS90 += $(HSLPACKDIR)/hsl_zd11/hsl_zd11d.f90
DTDEPS90 += $(HSLPACKDIR)/archive/hsl_zd13d.f90
DTDEPS90 += $(HSLPACKDIR)/hsl_mc56/hsl_mc56d.f90
DTDEPS90 += $(HSLPACKDIR)/hsl_zb01/hsl_zb01i.f90
DTDEPS90 += $(HSLPACKDIR)/hsl_mc68/hsl_mc68i.f90
DTDEPS90 += $(HSLPACKDIR)/hsl_mc69/hsl_mc69d.f90

# object files
objs := hsl_ma87d.o
objs +=	spllt_kernels_mod.o
objs +=	spllt_stf_mod.o spllt_factorization_mod.o spllt_mod.o

include Makefile.rules.inc

# hsl_ma87dt: hsl_ma87dt.o
# 	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)
# spllt_starpu_test:  DEFS += -DSPLLT_USE_STARPU
spllt_starpu_test: spllt_test.F90 $(objs) $(starpuobjs) common90.o dtdeps.o dtdeps90.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

spllt_test: spllt_test.F90 $(objs) common90.o dtdeps.o dtdeps90.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

spllt_testings: hsl_ma87dt.o dtdeps.o dtdeps90.o hsl_ma87d.o ddeps90.o common90.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

run_ma87: run_ma87.f90 dtdeps.o dtdeps90.o hsl_ma87d.o common90.o spllt_mod.o # ddeps90.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

## Generic compilation rules

# The following rules build a source file from a variable
dtdeps.f: $(DTDEPS)
	cat $^ > $@
common90.f90: $(COMMON90)
	cat $^ > $@
ddeps90.f90: $(DDEPS90)
	cat $^ > $@
dtdeps90.f90: $(DTDEPS90)
	cat $^ > $@

clean:
	-rm -f spllt_starpu_test spllt_test spllt_testings
	-rm -f run_ma87
	-rm -f *.o *.mod
	-rm -f metis.f *deps.f *deps90.f90 common.f common90.f90
