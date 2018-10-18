#!/bin/bash

#=======================
# GENERATED body of a case
#

#input : 
#       $1 is the number of tasks to declare
function declareTask(){
  caseOutput+="VARIABLE_OMP_DEP($1,${p_memory_space},"
  caseOutput+="${p_dep_array},${alpha},${beta})   &\n"
}

#input : 
#       $1 is the number of tasks to declare
#       $2 is the offset, i.e., the number of tasks already declared
function additionalDep(){
  caseOutput+="VAR_OMP_DEP_CONTD($1,${p_memory_space},"
  caseOutput+="${p_dep_array},${alpha},${beta}+$2)&\n"
}

#input : 
#       $1 is the file that contains the remaining omp declaration
function addRemainingOMPDecl(){
# file="../src/spllt_solve_fwd_block_omp_decl.F90"
  file=$1
  caseOutput+="`cat ${file}`\n"
}

#input : 
#       $1 is the file that contains the work to assign to the task
function addWork(){
# file="spllt_solve_fwd_block_worker.F90"
  file=$1
  caseOutput+="\nMACRO_SYMBOL include \"${file}\"\n"
  caseOutput+='\n!$omp end task\n'
}

#block="\nthreadID  = omp_get_thread_num()\n"
#block+="if(ndep_lvl .le. chunk) then\n"
#block+="! print '(a, i3, a, i3)', 'SLV      Task dep of ', dblk, ' [in : '\n"
#block+="! print *, p_dep\n"
#block+="  call spllt_solve_fwd_block_worker(m, n, nrhs, nthread, col, ldr, &\n"
#block+="    sa, offset, threadID, p_upd, p_rhs, p_lcol, p_index, p_xlocal)\n"
#block+="end if\n"
#block+='\n!$omp end task\n'

#=======================
# PARAMETERS
#
chunk=2
offset=0
NCASE_MIN=1
NCASE=1

p_memory_space="p_bc"
p_dep_array="p_dep"
alpha="alpha"
beta="beta"
caseOutput=""

# files
output="tmp_cases.F90"
workFile="spllt_solve_fwd_block_worker.F90"
ompDeclFile="../src/spllt_solve_fwd_block_omp_decl.F90"

#print in bold text
function bold(){

  echo -e "\033[1m$1\033[0m"

}

# @param $1 a string
function write() {
  echo -e "$1" >>${output}
}

function print_chelp(){
  echo -e "\t\t\t\tGeneral Commands Manual\t\t\t\t"
  bold    "NAME"
  echo -e "\t$0 - Generates list of fortran cases"
  bold    "SYNOPSIS"
  echo -e "\t$0 [OPTIONS]"
  bold    "DESCRIPTION"
  echo -e "\t${DESC}"
  bold    "OPTIONS"
  echo -e "\t-n         < #case     : default ${NCASE}>"
  echo -e "\t\t\tthe number of cases to generate"
  echo -e "\t--n_min    < integer   : default ${NCASE_MIN}>"
  echo -e "\t\t\tthe first generated case"
  echo -e "\t-o         < filename  : default ${output}>"
  echo -e "\t\t\tthe output filename where the cases are generated"
  echo -e "\t-w         < filename  : default ${workFile}>"
  echo -e "\t\t\tthe filename that contains the work of the task."
  echo -e "\t\t\tIt has to be located to the same folder as the"
  echo -e "\t\t\t generated included file will be placed"
  echo -e "\t--ompFile  < filename  : default ${ompDeclFile}>"
  echo -e "\t\t\tthe filename that contains the remaining omp instruction."
  echo -e "\t\t\tThis variable has to contain the path to it"
}

DESC="This script generates cases that create omp task with variable list "
DESC+="of dependencies in spllt code."

if [ $# -gt 0 ];then
  while [ $# -gt 0 ]; do
    case $1 in
      "-h" | "--help")
        print_chelp
        exit 0
        ;;
      "-n")
        shift
        NCASE=$1
        shift
        ;;
      "--n_min")
        shift
        NCASE_MIN=$1
        shift
        ;;
      "-o")
        shift
        if [ -e $1 ]; then
          echo -e "File $1 already exists. Overwrite (Y/n) ? [n] "
          read input
          if [ "${input}" == 'Y' -o "${input}" == 'y' ]; then
            output=$1
          else
            echo "Default output file named ${output} is used"
          fi
        else
          output=$1
        fi
        shift
        ;;
      "-w")
        shift
        workFile=$1
        shift
        ;;
      "--ompFile")
        shift
        ompDeclFile=$1
        if [ ! -e ${ompDeclFile} ]; then
          echo "Error, File ${ompDeclFile} does not exists."
          echo "Aborted script"
          exit 6
        fi
        shift
        ;;
      *)
        echo -e "Error in parameter $1\n\tPlease use -h option."
        exit 1
    esac
  done
fi

NCASE=$((NCASE_MIN + NCASE - 1))

echo "Generated ${NCASE} cases written in ${output}"

echo "" > ${output}

write '#include "macro_omp.def"\n'

if [ ${NCASE_MIN} -eq 0 ]; then
  ndep=0
  timer_id=$((ndep + 1))
  caseOutput+="MACRO_SYMBOL if defined(SPLLT_TIMER_TASKS_SUBMISSION)\n"
  caseOutput+="call spllt_tic('CASE(${ndep})',"
  caseOutput+="${timer_id}, task_manager%workerID, timer)\n"
  caseOutput+="MACRO_SYMBOL endif\n"
  caseOutput+='!$omp task &\n'
  addRemainingOMPDecl ${ompDeclFile}
  addWork ${workFile}  
  caseOutput+="\nMACRO_SYMBOL if defined(SPLLT_TIMER_TASKS_SUBMISSION)\n"
  caseOutput+="call spllt_tac(${timer_id}, task_manager%workerID, timer)\n"
  caseOutput+="MACRO_SYMBOL endif\n"
  write "${caseOutput}"
else
  for ndep in `seq ${NCASE_MIN} ${NCASE}`; do
    offset=0
    nldep=${ndep}
    timer_id=$((ndep + 1))
    
    caseOutput="case(${ndep})\n\n"
    caseOutput+="MACRO_SYMBOL if defined(SPLLT_TIMER_TASKS_SUBMISSION)\n"
    caseOutput+="call spllt_tic('CASE(${ndep})',"
    caseOutput+="${timer_id}, task_manager%workerID, timer)\n"
    caseOutput+="MACRO_SYMBOL endif\n"
    if [ ${ndep} -eq 1 ]; then
      declareTask ${ndep}
    else
      declareTask ${chunk}
      if [ ${ndep} -gt ${chunk} ]; then
        ndep=$((ndep - chunk))
        offset=$((offset + chunk))
        for nldep in `seq ${chunk} ${chunk} ${ndep}`; do
          additionalDep ${chunk} ${offset}
          offset=$((offset + chunk))
        done
        if [ $((ndep % chunk)) -ne 0 ]; then
          additionalDep $((ndep % chunk)) ${offset}
        fi
      fi
    fi
    addRemainingOMPDecl ${ompDeclFile}
    addWork ${workFile}  
    caseOutput+="\nMACRO_SYMBOL if defined(SPLLT_TIMER_TASKS_SUBMISSION)\n"
    caseOutput+="call spllt_tac(${timer_id}, task_manager%workerID, timer)\n"
    caseOutput+="MACRO_SYMBOL endif\n"
    write "${caseOutput}"
  done
fi
