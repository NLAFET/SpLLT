# These two rules tell us how to compile module files from source files
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<
