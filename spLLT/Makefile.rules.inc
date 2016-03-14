usestarpu=$(filter %USE_STARPU%, $(DEFS))
#usestarpu=

# starpusrc := StarPU/spllt_starpu_factorization_mod.o
# starpusrc += StarPU/spllt_starpu_kernels.o

# starpusrc := StarPU/starpu_f_mod.F90 StarPU/starpu_f_wrappers.c
# starpuobjs += $(patsubst %.F90, %.o, $(patsubst %.F90, %.o, $(notdir $(starpusrc))))

#ifdef usestarpu
# VPATH = StarPU
# starpuobjs += starpu_f_mod.o starpu_f_wrappers.o
# starpuobjs += StarPU/starpu_f_mod.o StarPU/starpu_f_wrappers.o
#endif

# Module build dependencies
spllt_test.F90: hsl_ma87d.o spllt_stf_mod.o spllt_mod.o dtdeps.o dtdeps90.o
spllt_stf_mod.o: hsl_ma87d.o spllt_factorization_mod.o spllt_mod.o
spllt_factorization_mod.o: spllt_mod.o spllt_kernels_mod.o spllt_starpu_factorization_mod.o 
spllt_kernels_mod.o: spllt_mod.o hsl_ma87d.o 
spllt_mod.o: dtdeps90.o  
# spllt_starpu_kernels.o: spllt_starpu_factorization_mod.o

ifdef usestarpu
spllt_stf_mod.o: starpu_f_mod.o
spllt_mod.o: starpu_f_mod.o
spllt_starpu_factorization_mod.o: spllt_mod.o spllt_starpu_kernels.o
starpu_f_mod.o: starpu_f_wrappers.o
endif

hsl_ma87dt.o: hsl_ma87d.o dtdeps90.o
hsl_ma87d.o: ddeps90.o common90.o 
ddeps90.o: common90.o
