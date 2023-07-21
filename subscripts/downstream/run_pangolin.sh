#!/bin/bash
#PBS -N run_pangolin
#PBS -l walltime=00:40:00
#PBS -l procs=20
#PBS -l pmem=2g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

#THE SCRIPT IS USED TO RUN PANGOLIN FROM SINGULARITY CONTAINER WITH DOWNSTREAM/MUTATION_REPORT.PY SCRIPT
source /opt/exp_soft/conda/anaconda3/bin/activate /opt/exp_soft/conda/anaconda3
ctrl_arg1=${1:-1} #IF NO PARAMETER VALUE PROVIDED, SET 1ST PARAMETER VALUE TO 1 (PERFORM PANGOLIN TYPING)
pango_sif_path=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "pangolin.sif") #PATH TO PANGOLIN CONTAINER
ctrl_arg2=${2:-"~/"} #IF NO PARAMETER SUPPLIED, SET TO 0, ELSE SET TO PROVIDED VALUE - PATH WHERE TO RUN PANGOLIN
module load singularity
now=$(date +"%m_%d_%Y")
if [ "$ctrl_arg2" == "~/" ]
then
    cd /mnt/home/groups/nmrl/cov_analysis/reports/report_${now}
else
    cd ${ctrl_arg2}
fi

awk '{print}' *.fasta > ${now}_combined.fasta
echo singularity run $pango_sif_path pangolin $(pwd)/${now}_combined.fasta -t 20 --outfile $(pwd)/${now}_lineage_report.csv # RUNNING PANGOLIN AND PROVIDING TIME-DEPENDENT OUTPUT FILE NAMING
singularity --silent run --bind $(pwd):$(pwd) $pango_sif_path pangolin $(pwd)/${now}_combined.fasta -t 20 --outfile $(pwd)/${now}_lineage_report.csv # RUNNING PANGOLIN AND PROVIDING TIME-DEPENDENT OUTPUT FILE NAMING

