import re, os, pandas as pd
from Bio import SeqIO
from subscripts.covipipe_utilities import covipipe_housekeeper as hk
from subscripts.src.modules import (
    Module, #base pipeline wrapper class
    module_data, #module metadata 
    pipeline_path #path to pipeline home folder
)



###############################################################
# Defining module classes to extend basic module functionality
###############################################################


class Covid_assembly(Module):
    '''Class extends Module and implements pipeline-specific file processing methods'''
    latest_processing_id = None #latest processing id read from file
    processing_id_dict = {} #to map processing ids to original ids for the current batch of samples
    renamed_result_file_map = {} #to map file paths annotated with processing ids and with run id removed to the raw file paths
    backward_result_file_map = {} #to map file restored annotated paths to specific annotated path
    renamed_fastq_file_map = {} #to store original and illumina-formatted fastq file name
    backward_fastq_file_map = {} #storing reversed renamed_fastq_file_map to restore original fastq file names after processing

    def fill_input_dict(self): #extends Module class
        '''Extends fill_input_dict method of the Module class to run integrity check on fastq files before adding them to the input dictionary.'''
        super(Covid_assembly, self).fill_input_dict() #running fill_input_dict of Module class - input_dict is filled by all fastq files found in folders and subfolders of input dir
        with open(f'{self.output_path}input_integrity_check.log', "w+") as logfile: #to log the results of integrity check
            for fastq in self.input_dict["001.fastq.gz"]: #for all detected fastq files
                if not hk.check_fastq_integrity(fastq): #if integrity check fails
                    logfile.write(f'{fastq} FAILED\n') #write to log as failed
                    if 'R1' in fastq: #if read 1 is corrupted
                        self.input_dict["001.fastq.gz"].remove(fastq) #remove from further processing
                        self.input_dict["001.fastq.gz"].remove(fastq.replace("R1", "R2"))
                    else: #if read 2 is corrupted
                        self.input_dict["001.fastq.gz"].remove(fastq) #remove from further processing
                        self.input_dict["001.fastq.gz"].remove(fastq.replace("R2", "R1"))
                # else: #if integrity intact
                #     logfile.write(f'{fastq} OK\n') #log as ok and proceed to next file


    def map_fastq_to_illumina(self):
        '''Generates illumina-format fastq file names for every fastq.gz file found in self.input_path and stores in self.renamed_fastq_file_map and self.backward_fastq_file_map'''
        fastq_path_list = hk.parse_folder(folder_pth_str=self.input_path, file_fmt_str='_[1,2].fastq.gz')
        prefix_patterns = [r'CO-[0-9]{5}_LVA[0-9]{3}_'] #Eurofins prefixes
        suffix_patterns = [r'_lib[0-9]{6}'] #Eurofins suffixes
        for path in fastq_path_list:
            new_path = path
            #replace prefix
            for prefix in prefix_patterns: 
                if re.search(prefix,new_path).group(0) is not None: #if pattern is detected in file name
                    new_path = re.sub(prefix,"",new_path) #replace
                    break #done with prefixes
            #replace suffix
            for suffix in suffix_patterns:
                if re.search(suffix,new_path).group(0) is not None: #if pattern is detected in file name
                    new_path = re.sub(suffix,"",new_path) #replace
                    break #done with suffixes
            #convert fastq part to illumina format
            if "_1.fastq.gz" in new_path: new_path = new_path.replace("_1.fastq.gz", "_R1_001.fastq.gz") #read_1
            if "_2.fastq.gz" in new_path: new_path = new_path.replace("_2.fastq.gz", "_R2_001.fastq.gz") #read_2
            self.renamed_fastq_file_map[path] = new_path #save to forward map
            self.backward_fastq_file_map[new_path] = path #save to backward map


    def convert_fastq_names(self, to_illumina:bool=True):
        '''
        Changes fastq file names for every fastq.gz file in self.input_path, using self.renamed_fastq_file_map. 
        If to_illumina set to False, performs renaming using self.backward_fastq_file_map to restore original file names.
        Creates renaming log files in the self.output_path, indicating succesful and failed renaming attempts.
        '''
        if to_illumina: #preprocessing fastq files
            with open(f'{self.output_path}fastq_forward_renaming.log', "w+") as rename_log: #log renaming process
                for path in self.renamed_fastq_file_map: 
                    try: #attempt renaming
                        os.rename(path, self.renamed_fastq_file_map[path])
                        rename_log.write(f'{path} {self.renamed_fastq_file_map[path]} OK\n') #log successful renaming
                    except OSError: #if failed
                        rename_log.write(f'{path} {self.renamed_fastq_file_map[path]} FAILED:{error}\n') #log renaming error message
        else: #restoring names of fastq files to original
            with open(f'{self.output_path}fastq_backward_renaming.log', "w+") as rename_log:
                for path in self.backward_fastq_file_map: 
                    try:
                        os.rename(path, self.backward_fastq_file_map[path])
                        rename_log.write(f'{path} {self.backward_fastq_file_map[path]} OK\n')
                    except OSError as error:
                        rename_log.write(f'{path} {self.backward_fastq_file_map[path]} FAILED:{error}\n')
                

    def restore_annotated_results(self):
        '''Checks if output directory contains annotations log file. If it is there, reads the contents and restores original name of all files according to the log.'''
        if os.path.isfile(f'{self.output_path}result_annotation.log') and os.stat(f'{self.output_path}result_annotation.log').st_size > 0:
            df = pd.read_table(f'{self.output_path}result_annotation.log', sep=" ", header=None)
            for i, path in enumerate(df[0]):
                self.backward_result_file_map[str(path)] = str(df[1][i])
                os.system(f'mv {str(df[1][i])} {str(path)} 2> /dev/null')


    def compute_processing_ids(self):
        '''Calculates a range of processing ids for current batch of samples and stores it in self.processing_id_dict'''
        latest_processing_id = hk.read_plaintext_file(self.config_file['latest_id_file'])[0] #reading latest id
        processing_id_list = [latest_processing_id[:3]+str(int(latest_processing_id[3:])+i+1).zfill(6) for i in range(len(self.sample_sheet['sample_id']))] #generating processing ids
        self.processing_id_dict = {ids[0]:f'{ids[1]}_{ids[0]}' for ids in zip(list(self.sample_sheet['sample_id']),processing_id_list)} #generating numbered ids
        if not self.backward_result_file_map:
            hk.overwrite_plaintext_file(self.config_file['latest_id_file'],processing_id_list[-1]) #updating processing id after range of ids was used 


    def fill_output_file_maps(self):
        '''Saves original-formatted file path pair in self.renamed_result_file_map'''
        for file_path in self.config_file["assembly_target_files"]: #generating new file path
            for id in self.processing_id_dict:
                if id in file_path:
                    new_path = re.sub(id, self.processing_id_dict[id], file_path) #adding numbered id
                    new_path = re.sub(r"(_S[0-9]{1}|_S[0-9]{2}|_S[0-9]{3})", "", new_path) #removing illumina run sample number
                    if 'qualimap' in file_path: #sample id in directory name and directory should be named
                        old_path = os.path.dirname(file_path)
                        new_path = os.path.dirname(new_path)
                        self.renamed_result_file_map[old_path] = new_path
                    else:
                        self.renamed_result_file_map[file_path] = new_path #storing for forward renaming
                    break


    def annotate_processed_files(self):
        '''Renames processed files using self.renamed_result_file_map and creates a renaming log file in the output directory'''
        with open(f'{self.output_path}result_annotation.log', 'w+') as ann_log: #opening log file for writing
            if self.backward_result_file_map:
                for path in self.backward_result_file_map: 
                    if os.path.isfile(path) or os.path.isdir(path): #if output file or directory exists
                        ann_log.write(f'{path} {self.backward_result_file_map[path]}\n') #log change
                        os.system(f'mv {path} {self.backward_result_file_map[path]} 2> /dev/null')
            else:
                for path in self.renamed_result_file_map: 
                    if os.path.isfile(path) or os.path.isdir(path): #if output file or directory exists
                        ann_log.write(f'{path} {self.renamed_result_file_map[path]}\n') #log change
                        os.system(f'mv {path} {self.renamed_result_file_map[path]} 2> /dev/null')


    def switch_sample_ids(self, new_to_old:bool = False):
        '''
        Replaces sample_id column in self.sample_sheet with ids that contain processing number and lack run number.
        If new_to_old is set to True, performs the reverse operation.
        '''
        if new_to_old:
            temp_dict = {re.sub(r"(_S[0-9]{1}|_S[0-9]{2}|_S[0-9]{3})", "", self.processing_id_dict[id]):id for id in self.processing_id_dict} #removing run ids
            self.sample_sheet = hk.map_replace_column(self.sample_sheet, temp_dict, 'sample_id', 'sample_id')
        else:
            temp_dict = {id:re.sub(r"(_S[0-9]{1}|_S[0-9]{2}|_S[0-9]{3})", "", self.processing_id_dict[id]) for id in self.processing_id_dict} #removing run ids
            self.sample_sheet = hk.map_replace_column(self.sample_sheet, temp_dict, 'sample_id', 'sample_id')


class Covid_downstream(Module):
    '''Class extends Module and implements pipeline-specific downstream processing methods'''
    metadata_table = None #To store SPKC metadata
    sequence_stats_table = None #To store statistics computed from consensus sequences


    def read_metadata_table(self):
        '''Reading metadata table to memeory'''


    def get_sequence_stats(self, path_to_multifasta:str): 
        '''
        Given path to multifasta, uses calculates GC content, N content and sequence length for each sequence.
        Fills sequence.sequence_stats_table with this data.
        '''


    def extract_specific_mutations(self, path_to_filter_file:str):
        '''Extracts mutation data from an annotated and csv-converted vcf file generated by the pipeline. Returns a dictionary where each mutation is maped to 0(not found) or 1 (found)'''


    def extract_coverage_depth(self, path_to_flagstat_report:str):
        '''Extracts median coverage and average coverage from samtools flagstat report. Returns tuple (median_cov:float, average_cov:float)'''




###############################################
# Defining wrapper functions to call from main
###############################################



def run_all(args, num_jobs):
    '''Wrapper function to run all modules sequentially.'''
    assembly = Covid_assembly(
            module_name='assembly',
            input_path=args.input,
            module_config=args.config,
            output_path=args.output_dir,
            run_mode=args.submit_modules,
            dry_run=args.dry_run,
            force_all=args.force_all,
            rule_graph=args.rule_graph,
            pack_output=args.pack_output,
            unpack_output=args.unpack_output,
            job_name=module_data['assembly']['job_name'],
            patterns=module_data['assembly']['patterns'],
            targets=module_data['assembly']['targets'],
            requests=module_data['assembly']['requests'],
            snakefile_path=module_data['snakefiles']['assembly'],
            cluster_config_path=module_data['cluster_config']
            )

    downstream = Covid_downstream(
        module_name='downstream', 
        input_path=assembly.output_path, 
        module_config=assembly.config_file, 
        output_path=args.output_dir, 
        run_mode=args.submit_modules,
        dry_run=args.dry_run,
        force_all=args.force_all,
        rule_graph=args.rule_graph,
        pack_output=args.pack_output,
        unpack_output=args.unpack_output,
        job_name=module_data['downstream']['job_name'],
        patterns=module_data['downstream']['patterns'],
        targets=module_data['downstream']['targets'],
        requests=module_data['downstream']['requests'],
        snakefile_path=module_data['snakefiles']['downstream'],
        cluster_config_path=module_data['cluster_config']
        )

    #Running assembly
    assembly.fill_input_dict()
    assembly.fill_sample_sheet()
    if assembly.unfold_output: assembly.unfold_output()
    assembly.make_output_dir()
    assembly.write_sample_sheet()
    assembly.fill_target_list()
    assembly.add_module_targets()
    assembly.add_output_dir()
    assembly.write_module_config()
    assembly.files_to_wd()
    try:
        assembly.run_module(job_count=num_jobs)
    except Exception as e:
        assembly.clear_working_directory() #to avoid manually moving files back to input
        raise e
    assembly.check_module_output()
    assembly.write_sample_sheet()
    assembly.clear_working_directory()

    #Connecting assembly to downstream
    downstream.receive_sample_sheet(assembly.supply_sample_sheet())
    if downstream.dry_run == "" and downstream.rule_graph == "": 
        samples_cleared = downstream.remove_invalid_samples(connect_from_module_name='assembly') #in dry run mode none of the rules are executed, hence all samples will be removed, causing error
        downstream.save_removed()
        if samples_cleared == 1: 
            if assembly.pack_output: assembly.fold_output()
            raise Exception('Missing files requested by downstream')

    #Running downstream
    downstream.fill_input_dict()
    downstream.add_fasta_samples()
    downstream.write_sample_sheet()
    downstream.fill_target_list()
    downstream.add_module_targets()
    downstream.write_module_config()
    downstream.files_to_wd()
    try:
        downstream.run_module(job_count=num_jobs)
    except Exception as e:
        downstream.clear_working_directory()
        raise e
    downstream.check_module_output()
    downstream.write_sample_sheet()
    downstream.clear_working_directory()


def run_assembly(args, num_jobs):
    '''Wrapper function to run only assembly module.'''
    assembly = Covid_assembly(
        module_name='assembly',
        input_path=args.input,
        module_config=args.config,
        output_path=args.output_dir,
        run_mode=None,
        dry_run=args.dry_run,
        force_all=args.force_all,
        rule_graph=args.rule_graph,
        pack_output=None,
        unpack_output=None,
        job_name=module_data['assembly']['job_name'],
        patterns=module_data['assembly']['patterns'],
        targets=module_data['assembly']['targets'],
        requests=module_data['assembly']['requests'],
        snakefile_path=module_data['snakefiles']['assembly'],
        cluster_config_path=module_data['cluster_config']
    )
    assembly.make_output_dir()
    assembly.map_fastq_to_illumina()
    assembly.convert_fastq_names()
    assembly.restore_annotated_results()
    assembly.fill_input_dict()
    assembly.fill_sample_sheet()
    assembly.write_sample_sheet()
    assembly.fill_target_list()
    assembly.add_module_targets()
    assembly.add_output_dir()
    assembly.write_module_config()
    assembly.files_to_wd()
    try:
        assembly.run_module(job_count=num_jobs)
    except Exception as e:
        assembly.clear_working_directory()
        raise e
    assembly.check_module_output()
    assembly.compute_processing_ids()
    assembly.fill_output_file_maps()
    assembly.annotate_processed_files()
    assembly.switch_sample_ids()
    assembly.write_sample_sheet()
    assembly.clear_working_directory()
    assembly.convert_fastq_names(to_illumina=False)
    
    #Housekeeping
    hk.asign_perm_rec(path_to_folder=assembly.output_path)
    hk.asign_perm_rec(f"{pipeline_path}/covipipe_job_logs/")
    hk.name_job_logs(pipeline_name='covipipe')
    if args.clean_job_logs:
        hk.remove_old_files(f"{pipeline_path}/covipipe_job_logs/")


def run_downstream(args, num_jobs):
    '''Wrapper function to run only downstream module.'''