#### Scope: Installation, configuration & usage 

# General note:
- Pipeline is designed to be run on [HPC cluster](https://hpc.rtu.lv/?lang=en) with a PBS job scheduler.
- It is required that [singularity](https://docs.sylabs.io) is available both on login node and computational nodes of the HPC.

# Prerequisites
* Linux machine with *sudo* access and <br>
    - [git](https://github.com/git-guides/install-git)<br>
    - [conda](https://docs.anaconda.com/anaconda/install/linux/)<br>
    - [python >= 3.7](https://www.python.org/downloads/)<br>
    - [singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps)


# Pipeline installation
1. Clone the repository to a machine sudo access:
```
git clone https://github.com/NMRL/SARS-CoV2-pipe/
```
2. Build the containers (sudo credentials will be required)
``` 
cd SARS-CoV2-pipe/config_files/s_recipes/
./build_singularity.sh 
```
* The images should appear under SARS-CoV2_assembly/tools/images/ folder.<br>
* Copy them to the HPC system where you plan to run the pipeline.

3. Log into the HPC cluster and clone the pipeline to a desired location:
```
git clone https://github.com/NMRL/SARS-CoV2-pipe/
```


# Configuration
* All configuration should be performed on the login node of HPC cluster. 

## Assembly mode configuration
* Assembly mode is used to call variants and generate consensus sequence for one or many samples.
* The only aggregation here is done using [multiqc](https://multiqc.info/) to assess the quality of given sample batch.
    1. Open **SARS-CoV2_assembly/config_files/yaml/config_modular.yaml**
        - **home_dir** - full path to the pipeline folder
        - **work_dir** - full path to the folder where intermediate processing results should be saved
        - **fastq_sif** - full path to fastq_processing.sif copied from local machine (see Installation section for details)
        - **multiqc_sif** - full path to fqc.sif copied from local machine (see Installation section for details)
        - **assembly_subscripts** - full path to SARS-CoV2_assembly/subscripts/assembly/ folder
        - **latest_id_file** - full path to SARS-CoV2_assembly/resources/assembly/current_id file (used to store latest used processing ids)
        - Replace all occurrences /mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/ with full path to SARS-CoV2_assembly/ in your system.
    2. Open SARS-CoV2_assembly/config_files/json/module_data.json and <br> 
    replace all occurrences of /mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/ <br> 
    with full path to SARS-CoV2_assembly/ in your system.
    3. Configuring [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
        - In **SARS-CoV2_assembly/config_files/yaml/config_modular.yaml**
            - set **fq_screen_config** to **SARS-CoV2_assembly/config_files/fastq_screen.conf**
        - Open SARS-CoV2_assembly/snakefiles/cov_assembly Snakefile
            - Under the rule fastq_screening/shell section:<br>
            replace **/home/groups/nmrl/** with **path to the folder where SARS-CoV2_assembly was cloned to**
        - By-default FastQ Screen will produce results only for SarsCoV2, please check the [documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/details.html?highlight=database)<br>
        for details on how to add more sequences into the database.
        - Configuration file can be modified according to local needs: SARS-CoV2_assembly/config_files/fastq_screen.conf
    4. Configuring [Qualimap](http://qualimap.conesalab.org/)
        - Copy the latest image from docker hub into a desired location on HPC:
        ```
        singularity pull docker://pegi3s/qualimap:latest
        ```
        - In **SARS-CoV2_assembly/config_files/yaml/config_modular.yaml**
            - set **qualimap_sif_path** to **full path to the resulting image**
    5. Open **SARS-CoV2_assembly/config_files/yaml/cluster.yaml**
        - Set all occurrences of **account** to the your HPC cluster account (this account will be charged according to number of CPU-hours used)
        - Set **outdir** to the full path to the folder where log files for individual jobs should be saved. The folder must exist and wit will grow in size very fast.

## Downstream mode configuration
* Downstream mode is used to batch-wise lineage detection and to combine the output for many samples into several report files. 
* The reports are designed to serve specific needs, but the scripts and templates are available and can be customized according to the requirements.

1. Create the environment used for downstream analysis:
    ```
    conda env create --file SARS-CoV2_assembly/config_files/conda_defs/rbase_env.yml
    ```
    - The environment can be placed in a shared directory using [--prefix](https://stackoverflow.com/questions/37926940/how-to-specify-new-environment-location-for-conda-create) argument
2. Open **SARS-CoV2_assembly/subscripts/covipipe_classes.py**:
    - Under #static paths section, replace all occurrences of /mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/ with full path to SARS-CoV2_assembly/ in your system
    - Set the static paths based on pathing you chose
        - self.report_folder_path = "full path to the folder where reports for different sample batches should be saved"
        - self.mutation_file_folder = "full path to the folder where variant files for all samples should be stored for aggregation into heatmap"
        - self.log_folder_path = "path to the folder where logs for all jobs are stored (see Assembly configuration/5)"
        - self.papermill_path = "full path to rbase_env(see Downstream configuration/1)/bin/papermill"
        - self.share_update_file_path - used to specify temp file that is used to copy files to a different location for sharing (probably useless function).
3. Further adjustment should be done to the scripts found in SARS-CoV2_assembly/subscripts/downstream/ based on your requirements.

# Usage
## First time use
```
python SARS-CoV2_assembly/covipipe.py -s -m assembly -i {INPUT_PATH} -o {OUTPUT_PATH} -t
```
- INPUT_PATH - path to a folder containing fastq.gz files named according to illumina (id_R{1,2}_001.fastq.gz) or (id_{1,2}.fastq.gz) convention.
- OUTPUT_PATH - path to a folder where output files should be stored after processing
- t - dry run flag - used to verify that all outputs can be generated from corresponding inputs and not syntax errors present in snakefiles.
- s - flag to install snakemake on HPC login node. 
    - By-default it is installed under /mnt/home/{username}/.conda/envs/mamba_env/envs/snakemake environment/
    - If the installation process fails, it might help to adjust the pathing in SARS-CoV2_assembly/subscripts/src/utilities.py - install_snakemake() function.
- for more options use:
```
python SARS-CoV2_assembly/covipipe.py --help
```
## Running the analysis
```
python SARS-CoV2_assembly/covipipe.py -s -m assembly -i {INPUT_PATH} -o {OUTPUT_PATH}
```