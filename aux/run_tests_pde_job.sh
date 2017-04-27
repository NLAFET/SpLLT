#BSUB -q numanlys-cpu
#BSUB -n 56
#BUSB -W 02:00
# #BSUB -o %J.log
# #BSUB -e %J.err
#BSUB -o run.log
#BSUB -e run.err
#BSUB -x
# #BSUB -app no_turbo

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

build="ma87"

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
        ;;
    starpu_trace)
        mkdir -p $outdir/starpu/traces
        module unload starpu/trunk
        module load starpu/trunk-fxt

        trace_dir=/tmp
        prof_file=prof_file_scarf462_0

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

# Poisson 3D
msz_list=(20 60 80 100 120 140 160)
# msz_list=(64)

# Poisson 2D
# msz_list=(32 64 96 128 256 512 1024)

# Biharmonic 2D
# msz_list=(32 64 128 256 512)

nb_list=(256 384 512 768 1024)
nemin_list=(32)

# ncpu_list=(7 14 21)
ncpu_list=(28)

. ./run_tests_helmholtz3d.sh
# . ./run_tests_poisson3d.sh
# . ./run_tests_poisson2d.sh
# . ./run_tests_biharmonic2d.sh
