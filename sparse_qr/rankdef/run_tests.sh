outdir=data
mkdir -p $outdir

sigma_list=(-2 -6 -12)

. ./run_tests_loop.sh
