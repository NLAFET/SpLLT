#!/bin/bash

NCASE=10
preprocess=0
generate=0
err=0
OUTPUTPATH="../src/include"
OMPDECLPATH="../src/include"
file_suffixe=''

#print in bold text
function bold(){

  echo -e "\033[1m$1\033[0m"

}

function print_chelp(){
  echo -e "\t\t\t\tGeneral Commands Manual\t\t\t\t"
  bold    "NAME"
  echo -e "\t$0 - Creates Fortran include files"
  bold    "SYNOPSIS"
  echo -e "\t$0 [OPTIONS]"
  bold    "DESCRIPTION"
  echo -e "\t${DESC}"
  bold    "OPTIONS"
  echo -e "\t-n <integer : default ${NCASE}>"
  echo -e "\t\t\tthe number of cases to generate"
  echo -e "\t-c"
  echo -e "\t\t\tgenerates raw cases files"
  echo -e "\t-E"
  echo -e "\t\t\tpreprocesses the raw cases files"
  echo -e "\t--ileave"
  echo -e "\t\t\tconsider interleave files"
  echo -e "\t--ileave2"
  echo -e "\t\t\tconsider new interleave files"
}

DESC="This script generates the Fortran file that should be included into "
DESC+="a source code.\n\tThe general usage is to generate the raw files and "
DESC+="then preprocess them."
DESC+="\n\tThis corresponds to a standard usage of : $0 -c -E"

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
      "-c")
        generate=1
        shift
        ;;
      "-E")
        preprocess=1
        shift
        ;;
      "--ileave")
        interleave=1
        file_suffixe='_il'
        shift
        ;;
      "--ileave2")
        interleave=1
        file_suffixe='_il2'
        shift
        ;;
      *)
        echo -e "Error in parameter $1\n\tPlease use -h option."
        exit 1
    esac
  done
fi

if [ ${preprocess} -eq 0 -a ${generate} -eq 0 ]; then
  echo "Nothing to do. Please rerun with -h option"
  exit 0
fi

raw_name=("fwd_block${file_suffixe}" "fwd_update${file_suffixe}" 
  "bwd_block${file_suffixe}" "bwd_update${file_suffixe}" 
  "bwd_node${file_suffixe}")
for step in ${raw_name[*]}; do

  raw_file="raw_${step}_cases.F90"
  output="${OUTPUTPATH}/spllt_${step}_cases.F90.inc"
 #workFile="${OUTPUTPATH}/spllt_solve_${step}_worker.F90.inc"
  workFile="spllt_solve_${step}_worker.F90.inc"
  ompDeclFile="${OMPDECLPATH}/spllt_solve_${step}_omp_decl.F90.inc"

  if [ ${generate} -eq 1 ]; then
    bash generator_case_omp.sh -n ${NCASE} -o ${raw_file} -w ${workFile} --ompFile ${ompDeclFile}
    err=$?
  fi
  if [ ${err} -eq 0 ]; then
    if [ ${preprocess} -eq 1 ]; then
      if [ -e ${raw_file} ]; then
        echo "Preprocess ${raw_file} with gfortran : output ${output}"
        gfortran -E -P ${raw_file} -o ${output}
      else
        echo "Unable to preprocess ${raw_file} : file does not exist."
      fi
    fi
  else
    echo "An error occurs during the generation of ${raw_file}."
  fi
done
exit 0
