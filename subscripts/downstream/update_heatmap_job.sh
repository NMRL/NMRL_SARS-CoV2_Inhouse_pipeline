#!/bin/bash
#PBS -N update_mut_matrix
#PBS -l walltime=03:00:00
#PBS -l procs=20
#PBS -l pmem=2g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

# eval "$(conda shell.bash hook)" 
source /opt/exp_soft/conda/anaconda3/bin/activate /opt/exp_soft/conda/anaconda3
conda activate /mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/tools/rbase_env

MUTATION_FILES='/mnt/home/groups/nmrl/cov_analysis/mutation_heatmap/mutation_files'
SCRIPTS_DIR='/mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/subscripts/downstream'
SUMMARY_PATH='/mnt/home/groups/nmrl/cov_analysis/analysis_history/summary_file*'

cd ${SCRIPTS_DIR}
python update_heatmap_data.py ${MUTATION_FILES} ${SUMMARY_PATH}