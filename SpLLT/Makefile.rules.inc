#usestarpu=$(filter %USE_STARPU, $(DEFS))
usestarpu=$(findstring USE_STARPU, $(DEFS))
useomp=$(findstring USE_OMP, $(DEFS))
useomptrace=$(findstring OMP_TRACE, $(DEFS))
useparsec=$(findstring USE_PARSEC, $(DEFS))
usegpu=$(findstring USE_GPU, $(DEFS))
#usestarpu=

# starpusrc := StarPU/starpu_f_mod.F90 StarPU/starpu_f_wrappers.c
# starpuobjs += $(patsubst %.F90, %.o, $(patsubst %.F90, %.o, $(notdir $(starpusrc))))

ifneq ($(usegpu),)
CFLAGS += $(IMAGMA)
FFLAGS += $(IMAGMA)
LIBS += $(LMAGMA)
LIBS += -L$(CUDADIR)/lib64 -lcusolver

gpuobjs := expand_buffer_kernels.o
endif

ifneq ($(usestarpu),)
VPATH += StarPU
starpuobjs := spllt_starpu_factorization_mod.o spllt_starpu_kernels.o starpu_f_mod.o starpu_f_wrappers.o
#starpuobjs :=
endif

ifneq ($(useomptrace),)
ompobjs += trace_mod.o get_wtime.o get_wtime_mod.o
endif

ifneq ($(useparsec),)
VPATH += PaRSEC

INCLUDES += -I./PaRSEC
CFLAGS += -m64 -mcx16
LIBS += $(LPARSEC) $(LHWLOC) $(LMPI)

parsecobjs := spllt_ptg_mod.o
parsecobjs += factorize.o
parsecobjs += spllt_parsec_factorization_mod.o
parsecobjs += spllt_parsec_blk_data.o
parsecobjs += spllt_parsec_wrap.o spllt_parsec_mod.o
endif

# Module build dependencies
spllt_test.o: hsl_ma87d.o spllt_data_mod.o spllt_analyse_mod.o spllt_factorization_mod.o spllt_mod.o dtdeps.o dtdeps90.o
spllt_factorization_mod.o: hsl_ma87d.o spllt_factorization_task_mod.o spllt_data_mod.o spllt_error_mod.o
spllt_ptg_mod.o: spllt_data_mod.o
spllt_factorization_task_mod.o: spllt_data_mod.o spllt_kernels_mod.o 
spllt_kernels_mod.o: spllt_data_mod.o hsl_ma87d.o 
spllt_analyse_mod.o: spllt_utils_mod.o spllt_data_mod.o hsl_ma87d.o
spllt_utils_mod.o: spllt_data_mod.o
spllt_mod.o: hsl_ma87d.o dtdeps90.o
spllt_data_mod.o: hsl_ma87d.o dtdeps90.o
spllt_c_interface.o: spllt_data_mod.o
# spllt_starpu_kernels.o: spllt_starpu_factorization_mod.o

ifneq ($(usestarpu),)
spllt_test.o: starpu_f_mod.o
spllt_factorization_mod.o: starpu_f_mod.o
spllt_mod.o: starpu_f_mod.o
spllt_data_mod.o: starpu_f_mod.o
spllt_factorization_task_mod.o: spllt_starpu_factorization_mod.o 
spllt_starpu_factorization_mod.o: spllt_data_mod.o spllt_starpu_kernels.o
starpu_f_mod.o: starpu_f_wrappers.o
endif

ifneq ($(useomptrace),)
spllt_factorization_task_mod.o: trace_mod.o
trace_mod.o: get_wtime_mod.o
endif

ifneq ($(useparsec),)
spllt_ptg_mod.o: spllt_parsec_factorization_mod.o spllt_parsec_mod.o
factorize.o: spllt_parsec_blk_data.o
endif

hsl_ma87dt.o: hsl_ma87d.o dtdeps90.o
hsl_ma87d.o: ddeps90.o common90.o 
ddeps90.o: common90.o

ifneq ($(useparsec),)
spllt_mod.o: spllt_factorization_task_mod.o 
endif