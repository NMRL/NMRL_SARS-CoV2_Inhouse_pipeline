#!/bin/bash

#THE SCRIPT IS USED TO RUN PANGOLIN FROM SINGULARITY CONTAINER WITH DOWNSTREAM/MUTATION_REPORT.PY SCRIPT
ctrl_arg1=${1:-1} #IF NO PARAMETER VALUE PROVIDED, SET 1ST PARAMETER VALUE TO 1 (PERFORM PANGOLIN TYPING)
pango_sif_path=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "pangolin.sif") #PATH TO PANGOLIN CONTAINER
ctrl_arg2=${2:-"~/"} #IF NO PARAMETER SUPPLIED, SET TO 0, ELSE SET TO PROVIDED VALUE - PATH WHERE TO RUN PANGOLIN
module load singularity
now=$(date +"%m_%d_%Y")
if [ "$ctrl_arg2" == "~/" ]
then
    cd ../reports/report_${now}
else
    cd ${ctrl_arg2}
fi

if [ $ctrl_arg1 -eq 1 ] #IF THE FIRST PASSED PARAMETER IS EQUAL TO 1, PERFORM PANGOLIN TYPING
then
    awk '{print}' *.fasta > ${now}_combined.fasta
    singularity run $pango_sif_path pangolin ${now}_combined.fasta --outfile ${now}_lineage_report.csv # RUNNING PANGOLIN AND PROVIDING TIME-DEPENDENT OUTPUT FILE NAMING
elif [ $ctrl_arg1 -eq 0 ] #IF THE FIRST PASSED PARAMETER IS EQUAL TO 0, SKIP PANGOLIN TYPING & REPORT TO THE TERMINAL THAT PANGOLIN TYPING IS SKIPPED
then
    echo "Pangolin typing skipped."
fi
