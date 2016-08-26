#BSUB -q numanlys-cpu
#BSUB -n 56
#BUSB -W 02:00
#BSUB -o %J.log
#BSUB -e %J.err
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
        module load starpu/trunk-nogpu
        # module load starpu/trunk-nogpu-nofxt
        module load metis/4.0.3
        module load hsl/latest
        module load spral/trunk
        ;;
esac

build_dir=`pwd`
id=`whoami`
outdir=data
#outsuffix="_NOSUB"
outsuffix=

echo "[run_tests] build dir: $build_dir"
#matrices=(JGD_Trefethen/Trefethen_20000)

ncpu_list=(27)
nb_list=(768)
# nb_list=(384)
nemin_list=(32)
trace_dir=/tmp
prof_file=prof_file_scarf462_0
# for matrix in ${matrices[@]}

mkdir -p $outdir
mkdir -p $outdir/spllt_starpu
mkdir -p $outdir/spllt_starpu/traces

#./run_tests.sh
. ./run_tests_trace.sh
