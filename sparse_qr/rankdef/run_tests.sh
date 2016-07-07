outdir=data
mkdir -p $outdir

sigma_list=(-6 -12)

. ./run_tests_loop.sh
