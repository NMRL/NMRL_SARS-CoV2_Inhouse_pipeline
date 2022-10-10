#!/bin/bash

eval "$(conda shell.bash hook)" 
conda activate /mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/tools/rbase_env

PATH=${1}
CUR_BATCH_PATH=${2}
MUT_STAT_PATH=${3}
OUTPUT_PATH=${4}
FOUND_PATH=${5}
REPORT_NOTEBOOK=${6}

/mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/tools/rbase_env/bin/papermill /mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/generate_weekly_report_with_plots.ipynb -p path ${PATH} -p curr_batch_path ${CUR_BATCH_PATH} -p mut_stat_path ${MUT_STAT_PATH} -p output_path ${OUTPUT_PATH} -p found_path ${FOUND_PATH} ${REPORT_NOTEBOOK}