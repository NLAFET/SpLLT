#BSUB -q numanlys-gpu
#BSUB -n 32
#BUSB -W 02:00
# #BSUB -o %J.log
# #BSUB -e %J.err
#BSUB -o run.log
#BSUB -e run.err
#BSUB -x
#BSUB -app no_turbo

case $HOSTNAME in
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

build_dir=`pwd`
id=`whoami`
outdir=../data

mkdir -p $outdir
mkdir -p $outdir/ssids-gpuonly

for matrix in `cat list.matrix`
do

    ~/builds/spral-gpuonly/spral_ssids --pos
done
