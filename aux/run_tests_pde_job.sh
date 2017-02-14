#BSUB -q numanlys-cpu
#BSUB -n 56
#BUSB -W 02:00
# #BSUB -o %J.log
# #BSUB -e %J.err
#BSUB -o run.log
#BSUB -e run.err
#BSUB -x
#BSUB -app no_turbo

module purge
module load use.own
module load gcc/6.1.0
module load intel/mkl/11.3.1.150
module load hwloc/1.11.4
module load starpu/trunk
module load parsec-icldistcomp/trunk
module load metis/4.0.3
module load hsl/latest
module load spral/master-gnu-6.1.0

build="starpu_prune"

build_dir=`pwd`
id=`whoami`
outdir=
#outsuffix="_NOSUB"
outdir=../data

echo "[run_tests] build dir: $build_dir"

# Create output directory if necessary
mkdir -p $outdir

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
    starpu|starpu_prune)
        mkdir -p $outdir/spllt_starpu
        mkdir -p $outdir/spllt_starpu/traces
        ;;
    starpu_nested_stf)
        mkdir -p $outdir/spllt_starpu_nested_stf
        mkdir -p $outdir/spllt_starpu_nested_stf/traces
        ;;
    gnu_omp|gnu_omp_prune)
        mkdir -p $outdir/spllt_omp
        mkdir -p $outdir/spllt_omp/gnu
        mkdir -p $outdir/spllt_omp/gnu/traces
        ;;
    intel_omp|intel_omp_prune)
        mkdir -p $outdir/spllt_omp
        mkdir -p $outdir/spllt_omp/intel
        mkdir -p $outdir/spllt_omp/intel/traces
        ;;
    ma87)
        mkdir -p $outdir/ma87
        ;;
esac

# Poisson 3D
# msz_list=(32 64 96 128 160)

# Poisson 2D
msz_list=(32 64 96 128 256 512 1024)
nb_list=(256 384 512 768 1024)
nemin_list=(32)

ncpu_list=(27)

# . ./run_tests_poisson3d.sh
. ./run_tests_poisson2d.sh
