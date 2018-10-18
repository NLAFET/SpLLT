#!/bin/bash
#This script checks the consistency of the solve phasis


#
#            ENV
#

BUILD="."
EXEC="test_solve_phasis"
nb=32
nrhs=2
ncpu=4
nworker=4

nmat=0 # 0 means all matrices in the subset directory

if [ $# -ge 1 ]; then
  nmat=$1
  echo "nmat $nmat used"
fi

#Size of the tile in the nodes of the tree
l_nb=""
for i in `seq -s ' ' 5 7`; do l_nb+=" $((2**i))";done
#Number of rhs
l_nrhs=`seq -s ' ' 1 10`
for i in `seq -s ' ' 4 7`; do l_nrhs+=" $((2**i))";done
#Number of worker i.e. threads used during the solve
l_nworker=`seq -s ' ' 4`
#Number of cpu used to generate the subtrees during the analysis
l_ncpu=$l_nworker

l_matrices=`find ${BUILD}/matrices/subset -iname "*.rb"`

if [ $nmat -ge 1 ];then
  l_matrices=`echo $l_matrices | cut -d ' ' -f 1-$nmat`
  echo "Run on subset of matrices : $l_matrices"
fi

#FUNCTIONS
function checkRes()
{
  checkRes_output_file=$1

  checkRes_NTEST=`grep "NTEST" $checkRes_output_file | awk '{print $NF}'`
  checkRes_res=`grep "Backward" $checkRes_output_file | tr '/' ' ' | awk '{if($5 == $6) sum++} END {if(sum>0) print sum; else print 0}'`
 #echo "#RES : $checkRes_res"
  if [ $checkRes_res -lt $checkRes_NTEST ]; then 
    grep "Backward" $checkRes_output_file | tr '/' ' '  | awk '{if($5 != $6) printf("%d != %d\n", $5, $6)}'
   #echo "ERROR occured : $success / $cpt STOP"
   #exit 1
  fi
}

function run()
{
  run_exec=$1
  run_param=$2
  run_output=$3
  run_errput=$4

  $run_exec $run_param >$run_output 2>$run_errput

  checkRes $output_file
  run_NTEST=$checkRes_NTEST
  if [ $checkRes_res -ne $run_NTEST ]; then
    echo "ERROR with : $run_exec $run_param >$run_output 2>$run_errput"
  fi
  success=$((success + checkRes_res))
  cpt=$((cpt + 1 * run_NTEST))
}

#STAT
cpt=0
success=0
output_file=/tmp/res_spllt

for mat in $l_matrices; do
  echo "Check ${mat}"
  param="--mat $mat"

  #Strong Scaling on the number of cpu
  echo -e "\t--ss_ncpu"
  for worker in $l_nworker; do
    export OMP_NUM_THREADS=$worker
    lparam="$param --nworker $worker --nb $nb --nrhs $nrhs --ncpu $worker"
    run ./$EXEC "$lparam" /dev/null $output_file
  done

  #Strong Scaling on the number of worker
  echo -e "\t--ss_nworker"
  export OMP_NUM_THREADS=$nworker
  for worker in $l_nworker; do
    export OMP_NUM_THREADS=$
    lparam="$param --nworker $worker --nb $nb --nrhs $nrhs --ncpu $ncpu"
    run ./$EXEC "$lparam" /dev/null $output_file
  done

  #Strong Scaling on the number of rhs
  echo -e "\t--ss_nrhs"
  export OMP_NUM_THREADS=$nworker
  for rhs in $l_nrhs; do
    lparam="$param --nworker $nworker --nb $nb --nrhs $rhs --ncpu $ncpu"
    run ./$EXEC "$lparam" /dev/null $output_file
  done

  #Strong Scaling on the number of rhs
  echo -e "\t--ss_nrhs without pruning"
  export OMP_NUM_THREADS=$nworker
  for rhs in $l_nrhs; do
    lparam="$param --nworker $nworker --nb $nb --nrhs $rhs --ncpu $ncpu"
    lparam+=" --no-prune-tree"
    run ./$EXEC "$lparam" /dev/null $output_file
  done

  #Strong Scaling on the block size
  echo -e "\t--ss_nb"
  for b in $l_nb; do
    lparam="$param --nworker $nworker --nb $b --nrhs $nrhs --ncpu $ncpu"
    run ./$EXEC "$lparam" /dev/null $output_file
  done
done


echo '======================='
echo "# Success : $success / $cpt"

