# These two rules tell us how to compile module files from source files
%.c: %.jdf
	$(DAGUEPP) --noline -i $< -o $*
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<
%.o: %.f
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<
%.o: %.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<
%.o: %.cu
	$(NVCC) $(NVCC_FLAGS) -dw $<
