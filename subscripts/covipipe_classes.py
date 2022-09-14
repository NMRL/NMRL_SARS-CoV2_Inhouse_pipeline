import re, os, pandas as pd, concurrent.futures
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

    def fill_input_dict(self, check_integrity:bool=False): #extends Module class
        '''
        Extends fill_input_dict method of the Module class to run integrity check on fastq files before adding them to the input dictionary.
        The integrity check runs by-default. Set check_integrity = False to avoid it.
        '''
        super(Covid_assembly, self).fill_input_dict() #running fill_input_dict of Module class - input_dict is filled by all fastq files found in folders and subfolders of input dir
        
        if check_integrity:
            fastq_total_count = len(self.input_dict["001.fastq.gz"]) #count total fastq found
            print(f"Performing fastq file integrity check:\n{fastq_total_count} total fastq files.") #announcing to the terminal
            failed_log_list = [] #to store data about fastq that failed integrity test
            with concurrent.futures.ProcessPoolExecutor() as executor: #using multiprocessing to speed-up (one check takes about 2.3 sec which translates to about 30 min for 384 nextseq batch using single thread)
                    results = [executor.submit(hk.check_fastq_integrity, fastq) for fastq in self.input_dict["001.fastq.gz"]] #submitting function calls to different processes
                    processed_count = 0 #counter for processed samples
                    for f in concurrent.futures.as_completed(results): #collecting results from different processes
                        result = f.result() #accessing results
                        if not result[0]:
                            failed_log_list.append(f'{result[1]} FAILED\n') #save result as failed to write to log later
                            if 'R1' in result[1]: #if read 1 is corrupted
                                self.input_dict["001.fastq.gz"].remove(result[1]) #remove from further processing
                                self.input_dict["001.fastq.gz"].remove(result[1].replace("R1", "R2"))
                            else: #if read 2 is corrupted
                                self.input_dict["001.fastq.gz"].remove(result[1]) #remove from further processing
                                self.input_dict["001.fastq.gz"].remove(result[1].replace("R2", "R1"))
                        processed_count += 1 #counting processed files
                        hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            with open(f'{self.output_path}input_integrity_check.log', "w+") as logfile: #loggin failed results
                for record in failed_log_list: logfile.write(f"{record}\n")
            print(f"Fastq integrity check complete!\n") #report completion
            if len(failed_log_list) > 0: #if some files have failed
                print(f"Please check {self.output_path}input_integrity_check.log for details on files that failed.")
            else: #if no files failed
                print(f"All fastq files are intact!")
                            


    def map_fastq_to_illumina(self):
        '''Generates illumina-format fastq file names for every fastq.gz file found in self.input_path and stores in self.renamed_fastq_file_map and self.backward_fastq_file_map'''
        fastq_path_list = hk.parse_folder(folder_pth_str=self.input_path, file_fmt_str='_[1,2].fastq.gz')
        fastq_total_count = len(fastq_path_list)
        if fastq_total_count > 0: #if non-illumina formatted samples present
            with concurrent.futures.ProcessPoolExecutor() as executor:
                print(f"Generating illumina format fastq file names:\n{fastq_total_count} total fastq files.") #announcing to the terminal
                results = [executor.submit(hk.name_formatter, path) for path in fastq_path_list] #submitting function calls to different processes
                processed_count = 0 #counter for processed samples
                for f in concurrent.futures.as_completed(results):
                    result = f.result()
                    self.renamed_fastq_file_map[result[0]] = result[1]
                    self.backward_fastq_file_map[result[1]] = result[0]
                    processed_count += 1
                    hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            print(f"Illumina name generation complete!\n") #report completion

    def convert_fastq_names(self, to_illumina:bool=True):
        '''
        Changes fastq file names for every fastq.gz file in self.input_path, using self.renamed_fastq_file_map. 
        If to_illumina set to False, performs renaming using self.backward_fastq_file_map to restore original file names.
        Creates renaming log files in the self.output_path, indicating succesful and failed renaming attempts.
        '''

        if to_illumina: #preprocessing fastq files
            renaming_logs = []
            fastq_total_count = len(self.renamed_fastq_file_map)
            if fastq_total_count > 0:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    print(f"Renaming fastq files to illumina format names:\n{fastq_total_count} total fastq files.") #announcing to the terminal
                    results = [executor.submit(hk.renamer, path, self.renamed_fastq_file_map[path]) for path in self.renamed_fastq_file_map] #submitting function calls to different processes
                    processed_count = 0

                    for f in concurrent.futures.as_completed(results):
                        result = f.result()
                        renaming_logs.append(result)
                        processed_count += 1
                        hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
                print(f"Fastq renaming to Illumina format finished!\n") #report completion
                with open(f'{self.output_path}fastq_forward_renaming.log', "w+") as rename_log: #log renaming process
                    for record in renaming_logs: rename_log.write(record)
                
        else: #restoring names of fastq files to original
            renaming_logs = []
            fastq_total_count = len(self.backward_fastq_file_map)
            if fastq_total_count > 0:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    print(f"Restoring original (non-illumina) fastq file names:\n{fastq_total_count} total fastq files.") #announcing to the terminal
                    results = [executor.submit(hk.renamer, path, self.backward_fastq_file_map[path]) for path in self.backward_fastq_file_map] #submitting function calls to different processes
                    processed_count = 0

                    for f in concurrent.futures.as_completed(results):
                        result = f.result()
                        renaming_logs.append(result)
                        processed_count += 1
                        hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
                print(f"Fastq renaming to original format finished!\n") #report completion
                with open(f'{self.output_path}fastq_backward_renaming.log', "w+") as rename_log: #log renaming process
                    for record in renaming_logs: rename_log.write(record)
                

    def restore_annotated_results(self):
        '''Checks if output directory contains annotations log file. If it is there, reads the contents and restores original name of all files according to the log.'''
        if os.path.isfile(f'{self.output_path}result_annotation.log') and os.stat(f'{self.output_path}result_annotation.log').st_size > 0:
            df = pd.read_table(f'{self.output_path}result_annotation.log', sep=" ", header=None)
            fastq_total_count = len(df)
            processed_count = 0
            with concurrent.futures.ProcessPoolExecutor() as executor:
                print(f"Removing processing id annotation from output file names for reprocessing:\n{fastq_total_count} total output files.") #announcing to the terminal
                results = [] 
                for i, path in enumerate(df[0]):#submitting function calls to different processes
                    self.backward_result_file_map[str(path)] = str(df[1][i])
                    results.append(executor.submit(hk.renamer, str(df[1][i]), str(path)))

                for _ in concurrent.futures.as_completed(results):
                    processed_count += 1
                    hk.printProgressBar(processed_count, fastq_total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            print(f"Finished removing processing ids for reprocessing!\n") #report completion


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
                    new_path = re.sub(r"_S[0-9]*\/", "", new_path) #removing illumina run sample number
                    if 'qualimap' in file_path: #sample id in directory name and directory should be named
                        old_path = os.path.dirname(file_path)
                        new_path = os.path.dirname(new_path)
                        self.renamed_result_file_map[old_path] = new_path
                    else:
                        self.renamed_result_file_map[file_path] = new_path #storing for forward renaming
                    break


    def annotate_processed_files(self):
        '''Renames processed files using self.renamed_result_file_map and creates a renaming log file in the output directory'''
        if self.backward_result_file_map:
            total_count = len(self.backward_result_file_map)
            processed_count = 0
            annotation_logs = []
            with concurrent.futures.ProcessPoolExecutor() as executor:
                print(f"Reannotating output files with processing ids:\n{total_count} total output files.") #announcing to the terminal
                results = [executor.submit(hk.renamer, path, self.backward_result_file_map[path]) for path in self.backward_result_file_map] #submitting function calls to different processes
                for f in concurrent.futures.as_completed(results):
                    result = f.result()
                    annotation_logs.append(result)
                    processed_count += 1
                    hk.printProgressBar(processed_count, total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)
            print(f"Result annotation finished!\n") #report completion
            with open(f'{self.output_path}result_annotation.log', 'w+') as ann_log: #opening log file for writing  
                for log in annotation_logs: ann_log.write(log)         
        else:
            total_count = len(self.renamed_result_file_map)
            processed_count = 0
            annotation_logs = []
            with concurrent.futures.ProcessPoolExecutor() as executor:
                print(f"Annotating output files with processing ids:\n{total_count} total output files.") #announcing to the terminal
                results = [executor.submit(hk.renamer, path, self.renamed_result_file_map[path]) for path in self.renamed_result_file_map] #submitting function calls to different processes
                for f in concurrent.futures.as_completed(results):
                    result = f.result()
                    annotation_logs.append(result)
                    processed_count += 1
                    hk.printProgressBar(processed_count, total_count, prefix = 'Progress:', suffix = 'Complete', length = 50)

            with open(f'{self.output_path}result_annotation.log', 'w+') as ann_log: #opening log file for writing  
                for log in annotation_logs: ann_log.write(log)

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
        '''Reading metadata table to memory'''


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
    assembly.fill_input_dict(check_integrity=args.fq_integrity_check)
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