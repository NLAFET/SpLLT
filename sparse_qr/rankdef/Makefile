# compilers
CC      = gcc
FC      = gfortran

# compiler options
FCFLAGS = -g -O2
CFLAGS  = -g -O2

# link options
LDFLAGS  = $(FCFLAGS) 

OBJS = run_rd.o

INCLUDES = $(IQRM) $(ISTARPU)
LIBS = $(LQRM) $(LSTARPU) $(LMETIS) $(LLAPACK) $(LBLAS)

run_rd: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean: 
	(rm -f *.o *.mod)

%.mod: %.o
	@true

%.o: %.c
	$(CC)  $(CFLAGS)  $(CINCLUDES) $(INCLUDES) $(CDEFS) $(PREC) -c $<
%.o: %.f
	$(FC)  $(FCFLAGS) $(FINCLUDES) $(INCLUDES) $(FDEFS) $(DEFINE_PREPEND)$(PREC) -c $<
%.o: %.f90				       
	$(FC)  $(FCFLAGS) $(FINCLUDES) $(INCLUDES) $(FDEFS) $(DEFINE_PREPEND)$(PREC) -c $<
%.o: %.F
	$(FC)  $(FCFLAGS) $(FINCLUDES) $(INCLUDES) $(FDEFS) $(DEFINE_PREPEND)$(PREC) -c $<
%.o: %.F90
	$(FC)  $(FCFLAGS) $(FINCLUDES) $(INCLUDES) $(FDEFS) $(DEFINE_PREPEND)$(PREC) -c $<
