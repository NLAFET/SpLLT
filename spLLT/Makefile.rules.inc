usestarpu=$(filter %USE_STARPU%, $(DEFS))
#usestarpu=

# Module build dependencies
spllt_test.F90: hsl_ma87d.o spllt_stf_mod.o spllt_mod.o dtdeps.o dtdeps90.o
spllt_stf_mod.o: hsl_ma87d.o spllt_factorization_mod.o spllt_mod.o starpu_f_mod.o
spllt_starpu_factorization_mod.o: spllt_mod.o spllt_starpu_kernels.o
spllt_factorization_mod.o: spllt_mod.o spllt_kernels_mod.o spllt_starpu_factorization_mod.o 
spllt_kernels_mod.o: spllt_mod.o hsl_ma87d.o 
spllt_mod.o: dtdeps90.o starpu_f_mod.o  
# spllt_starpu_kernels.o: spllt_starpu_factorization_mod.o

ifneq ($(strip $(usestarpu)),)
starpu_f_mod.o: starpu_f_wrappers.o
endif

hsl_ma87dt.o: hsl_ma87d.o dtdeps90.o
hsl_ma87d.o: ddeps90.o common90.o 
ddeps90.o: common90.o
