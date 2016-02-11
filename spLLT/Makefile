#Select compiler and options from below
############################
FC = gfortran
FFLAGS = -Wall -pedantic -fbounds-check -march=native -g -fbacktrace -fopenmp
#FFLAGS = -Wall -pedantic -fbounds-check -march=native -g -Wunused -fbacktrace -fopenmp
#FFLAGS += -fprofile-arcs -ftest-coverage
############################
#FC = ifort
#FFLAGS = -g -C -static -openmp
############################

FFLAGS += $(ISPRAL)

LIBS = $(LBLAS) $(LLAPACK) $(LSPRAL)

all: spllt_test

# HSL deps
DTDEPS = $(HSLPACKDIR)/fa14/fa14d.f \
	$(HSLPACKDIR)/mc34/mc34d.f \
	$(HSLPACKDIR)/ym11/ym11d.f \
	$(HSLPACKDIR)/mc47/mc47d.f \
	$(HSLPACKDIR)/mc59/mc59d.f \
        $(HSLPACKDIR)/mc54/mc54d.f

# HSL f90 deps
COMMON90 = $(HSLPACKDIR)/hsl_mc78/hsl_mc78i.f90
DDEPS90  = $(HSLPACKDIR)/hsl_mc34/hsl_mc34d.f90
DTDEPS90 = $(HSLPACKDIR)/hsl_fa14/hsl_fa14d.f90 \
	$(HSLPACKDIR)/hsl_zd11/hsl_zd11d.f90

# Module build dependencies
spllt_test.f90: hsl_ma87d.o spllt_stf_mod.o spllt_mod.o dtdeps.o dtdeps90.o
spllt_stf_mod.o: hsl_ma87d.o spllt_factorization_mod.o spllt_mod.o
spllt_factorization_mod.o: spllt_mod.o 
spllt_mod.o: dtdeps90.o
hsl_ma87dt.o: hsl_ma87d.o dtdeps90.o
hsl_ma87d.o: ddeps90.o common90.o 
ddeps90.o: common90.o

# hsl_ma87dt: hsl_ma87dt.o
# 	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

spllt_test: spllt_test.f90 hsl_ma87d.o spllt_stf_mod.o spllt_factorization_mod.o spllt_mod.o ddeps90.o common90.o dtdeps.o dtdeps90.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

spllt_testings: hsl_ma87dt.o dtdeps.o dtdeps90.o hsl_ma87d.o ddeps90.o common90.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

## Generic compilation rules

# These two rules tell us how to compile module files from source files
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

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
	-rm -f *.o *.mod
	-rm -f metis.f *deps.f *deps90.f90 common.f common90.f90
