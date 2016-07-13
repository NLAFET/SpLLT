#BSUB -q numanlys-cpu
#BSUB -n 56
#BUSB -W 02:00
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
esac

build="ma87"
build_dir=`pwd`
id=`whoami`
outdir=data
#outsuffix="_NOSUB"
outsuffix=

echo "[run_tests] build dir: $build_dir"
#matrices=(JGD_Trefethen/Trefethen_20000)

prof_file=prof_file_scarf462_0
# for matrix in ${matrices[@]}

declare -a problems=("matrix:Schmid/thermal2 ncpu:27 nb:1024 ne:32"
    "matrix:AMD/G3_circuit ncpu:27 nb:512 ne:32")

mkdir -p $outdir
mkdir -p $outdir/ma87

case $build in
    stf)
        mkdir -p $outdir/spllt_stf
        ;;
    parsec)
        mkdir -p $outdir/spllt_parsec
        mkdir -p $outdir/spllt_parsec/traces
        # make clean
        # make parsec
        ;;
    starpu)
        mkdir -p $outdir/spllt_starpu
        mkdir -p $outdir/spllt_starpu/traces
        ;;
    gnu_omp)
        mkdir -p $outdir/spllt_omp
        mkdir -p $outdir/spllt_omp/gnu
        mkdir -p $outdir/spllt_omp/gnu/traces
        ;;
    intel_omp)
        mkdir -p $outdir/spllt_omp
        mkdir -p $outdir/spllt_omp/intel
        mkdir -p $outdir/spllt_omp/intel/traces
        ;;
    ma87)
        mkdir -p $outdir/ma87
        ;;
esac

#./run_tests.sh
. ./run_test_list.sh
