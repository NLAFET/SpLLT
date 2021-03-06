#BSUB -q numanlys-cpu
# # BSUB -q numanlys-gpu
# #BSUB -n 56
# #BUSB -W 02:00
# #BSUB -o %J.log
# #BSUB -e %J.err
#BSUB -o run.log
#BSUB -e run.err
#BSUB -x
# #BSUB -app no_turbo

case $HOSTNAME in
    gauss)
        module purge
        module load use.own
        module load gnu/comp/default
        module load gnu/mkl/seq/11.2.0
        module load hwloc/1.11.2
        module load fxt/0.3.1
        module load starpu/trunk
        module load metis/4.0.3
        module load hsl/latest
        module load spral/trunk
        ;;
    cn202.scarf.rl.ac.uk | cn255.scarf.rl.ac.uk)
        module purge
        module load use.own
        module load gcc/5.3.0
        module load intel/mkl/11.3.1.150
        module load hwloc/1.11.2
        module load starpu/trunk-nogpu-nofxt
        module load parsec-icldistcomp/trunk
        module load metis/4.0.3
        module load hsl/latest
        module load spral/trunk
        ;;
    cn1g01.gpu.rl.ac.uk)
        module purge
        module load use.own
        module load gcc/5.3.0
        module load cuda/8.0.44
        module load intel/mkl/11.2.0.090
        module load hwloc/gpu-1.11.4
        module load metis/4.0.3
        module load starpu/trunk-gpu
        module load hsl/latest
        module load spral/master-gnu-5.3.0
esac

# Choose build
# - starpu: StarPU version of SpLLT
# - starpu_prune: StarPU version of SpLLT + tree prunnig strategy
# - gnu_omp: SpLLT OpenMP version. GNU compiler 
# - gnu_intel: SpLLT OpenMP version. Intel compiler 
# - parsec: PaRSEC version of SpLLT

build="gnu_omp_prune"

build_dir=`pwd`
id=`whoami`
outdir=../data
#outsuffix="_NOSUB"
outsuffix=

echo "[run_tests] build dir: $build_dir"
#matrices=(JGD_Trefethen/Trefethen_20000)

mkdir -p $outdir
mkdir -p $outdir/ma87

case $build in
    stf)
        mkdir -p $outdir/stf
        ;;
    parsec)
        mkdir -p $outdir/parsec
        mkdir -p $outdir/parsec/traces
        # make clean
        # make parsec
        ;;
    starpu|starpu_prune)
        mkdir -p $outdir/starpu
        mkdir -p $outdir/starpu/traces
        ;;
    starpu_nested_stf)
        mkdir -p $outdir/starpu_nested_stf
        mkdir -p $outdir/starpu_nested_stf/traces
        ;;
    gnu_omp|gnu_omp_prune)
        mkdir -p $outdir/omp
        mkdir -p $outdir/omp/gnu
        mkdir -p $outdir/omp/gnu/traces
        ;;
    intel_omp|intel_omp_prune)
        mkdir -p $outdir/omp
        mkdir -p $outdir/omp/intel
        mkdir -p $outdir/omp/intel/traces
        ;;
    ma87)
        mkdir -p $outdir/ma87
        ;;
esac

ncpu_list=(27)
nb_list=(256 384 512 768 1024 1536)
nemin_list=(32)

. ./run_tests_loop.sh
