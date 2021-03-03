#!/bin/bash

PATHS_TO_RUN=("." "run_sga" "graph")

NUMBER_OF_CORES=2

ARGUMENTS="-p -j ${NUMBER_OF_CORES} --dryrun --use-conda --debug-dag" #"--rulegraph"
#--report report.html"

########################################################################

BASE_PATH=$(pwd)

for path in "${PATHS_TO_RUN[@]}"
do
  cd ${path} || exit
  echo "In path: ${path}"
  snakemake ${ARGUMENTS}
  cd ${BASE_PATH} || exit
done
