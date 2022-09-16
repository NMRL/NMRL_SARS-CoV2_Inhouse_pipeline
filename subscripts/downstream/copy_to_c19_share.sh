#!/bin/bash
#PBS -N update_c19_data_share
#PBS -l walltime=24:00:00
#PBS -l procs=48
#PBS -q long
#PBS -j oe
#PBS -A rakus


###This script can be used to copy bam, fasta & vcf files from covid_output/ to c19data_share based on sample ids provided as list in c19_share_update.txt
#Usage:
#1.Add sample ids to the c19_share_update.txt file (replace old ids)
#2.add n(dry-run) to the RSYNC_OPTS string
#3.run the script to see if it works correctly
#4.Fix the script or the file list if there are some problems
#5.Remove n from RSYNC_OPTS and run the script when you are sure that everything is fine

cd "/home/groups/nmrl/cov_analysis"
C19_SHARE_PATH="/home/groups/c19data_share/NMRL/"
RSYNC_OPTS='rsync -ghpv --ignore-existing'
SAMPLE_COUNT=$(cat c19_share_update.txt | grep -v 'processing' | wc -l)
RUN_PATH=${1} #IF NO FOLDER NAME IS PROVIDED, THE SCRIPT WILL LOOK IN ALL COVID_OUTPUT, WHICH TAKES LONGER
#g - keep groups
#p - keep permissions
#n - dry run (for testing)
#v - more information to standard output
#h - human-readable


#BUILDING LIST OF FILES
cat c19_share_update.txt | xargs -n1 -P48 -I% find ./covid_output/ -type f -name %"_sorted.bam" > bam_list.txt
cat c19_share_update.txt | xargs -n1 -P48 -I% find ./covid_output/ -type f -name %".vcf" > vcf_list.txt 
cat c19_share_update.txt | xargs -n1 -P48 -I% find ./covid_output/ -type f -name %"_consensus.fasta" > fasta_list.txt


#COPY FILES USING MULTIPROCESSING
cat bam_list.txt | xargs -n1 -P48 -I% ${RSYNC_OPTS} % ${C19_SHARE_PATH}bam/
cat vcf_list.txt | xargs -n1 -P48 -I% ${RSYNC_OPTS} % ${C19_SHARE_PATH}vcf/
cat fasta_list.txt | xargs -n1 -P48 -I% ${RSYNC_OPTS} % ${C19_SHARE_PATH}fasta/


#REMOVE LIST FILES
rm bam_list.txt vcf_list.txt fasta_list.txt

# if [ -d "./covid_output/$RUN_PATH" ]
# then
#   cat c19_share_update.txt | while IFS="" read a; do $RSYNC_OPTS $(find ./covid_output/${RUN_PATH} -type f -name "${a}_sorted.bam") ${C19_SHARE_PATH}bam/; $RSYNC_OPTS $(find ./covid_output/${RUN_PATH} -type f -name "${a}.vcf") ${C19_SHARE_PATH}vcf/; $RSYNC_OPTS $(find ./covid_output/${RUN_PATH} -type f -name "${a}_consensus.fasta") ${C19_SHARE_PATH}fasta/; done
# else
#   echo 'Error: ${RUN_PATH} not found.'
#   exit 1
# fi

BAM_CONTENT=$(ls ${C19_SHARE_PATH}bam/)
VCF_CONTENT=$(ls ${C19_SHARE_PATH}vcf/)
FASTA_CONTENT=$(ls ${C19_SHARE_PATH}fasta/)
BAM_COUNT=$(cat c19_share_update.txt | while IFS="" read a; do echo "$BAM_CONTENT" | grep ${a}_sorted.bam ; done | wc -l)
VCF_COUNT=$(cat c19_share_update.txt | while IFS="" read a; do echo "$VCF_CONTENT" | grep ${a}.vcf ; done | wc -l)
CONSENSUS_COUNT=$(cat c19_share_update.txt | while IFS="" read a; do echo "$FASTA_CONTENT" | grep ${a}_consensus.fasta ; done | wc -l)

if [ $BAM_COUNT -eq $SAMPLE_COUNT ]
then
    echo "BAM:OK ${BAM_COUNT} matches ${SAMPLE_COUNT}"
else
    echo "BAM:WARNING! ${BAM_COUNT} does not match ${SAMPLE_COUNT}"
fi

if [ $VCF_COUNT -eq $SAMPLE_COUNT ]
then
    echo "VCF:OK ${VCF_COUNT} matches ${SAMPLE_COUNT}"
else
    echo "VCF:WARNING! ${VCF_COUNT} does not match ${SAMPLE_COUNT}"
fi 

if [ $CONSENSUS_COUNT -eq $SAMPLE_COUNT ]
then
    echo "FASTA:OK ${CONSENSUS_COUNT} matches ${SAMPLE_COUNT}"
else
    echo "FASTA:WARNING! ${CONSENSUS_COUNT} does not match ${SAMPLE_COUNT}"
fi

dry_run=$(python -c "print('n' in '${RSYNC_OPTS}'.split(' ')[1])") #CHECK IF IT IS DRY RUN, IF SO, THIS STEP SHOULD BE SKIPPED

if [ $dry_run == 'False' ]; then
    echo 'Updating metadata file and reaplying permissions'
    rm ${C19_SHARE_PATH}summary_file*
    cp downstream/summary_file* ${C19_SHARE_PATH}
    chmod -R 775 ${C19_SHARE_PATH}bam/
    chmod -R 775 ${C19_SHARE_PATH}vcf/
    chmod -R 775 ${C19_SHARE_PATH}fasta/
    chmod 775 ${C19_SHARE_PATH}summary_file*

fi
