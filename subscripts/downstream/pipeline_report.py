#!/mnt/beegfs2/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/tools/rbase_env/bin/python
##########
#IMPORTS
##########

import sys, os, pandas as pd, re, time, concurrent.futures, subprocess, argparse, shutil, pathlib, filecmp, numpy as np, warnings
from datetime import datetime


#################
#PATHS & SETUPS
#################

#SUPPRESSING PANDAS WARNINGS
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

#PATHS
subprocess_path = f'/mnt/beegfs2/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/'
covid_output_path = "/mnt/beegfs2/home/groups/nmrl/cov_analysis/covid_output/"
report_folder_path = f'/mnt/beegfs2/home/groups/nmrl/cov_analysis/reports'
raw_folder_path = f'/mnt/beegfs2/home/groups/nmrl/cov_analysis/raw'
filter_path = f'/mnt/beegfs2/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/resources/downstream/report_filters.txt'
default_metadata_path = "/mnt/beegfs2/home/groups/nmrl/cov_analysis/metadata/spkc_latest_3_month.csv"
log_folder_path = "/mnt/beegfs2/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/covipipe_job_logs/"

#TIMESTAMPS & STATIC STRINGS
pipeline_string = 'bwa 0.7.17-r1198-dirty mem'
timestr_fasta = time.strftime("%m_%d_%Y")
analysis_institution = 'NMRL'

#FILTERING SETUPS
f_value = 65 #SET MAX ACCEPTABLE FREQUENCY VALUE TO CONFIRM THAT MUTATION IS DETECTED (% of reads)
p_value = 0.05 #SET MAX ACCEPTABLE P VALUE TO CONFIRM THAT MUTATION IS DETECTED (ALTERNATIVE FOR f_value)

#CONTROLLERS
skip_excel_mining = False
skip_sequence_stats = False
rename_fasta_header = True
use_p_filter = False
use_f_filter = True

#FILE FORMAT STRINGS TO LOOK FOR IN PIPELINE OUTPUT FILES
context_path_map = {
    "ann.csv":[],
    "seq_depth.txt":[],
    "mapped_report.txt":[],
    "consensus.fasta":[],
    "fastp_report.json":[],
    "cutadapt_log.txt":[],
    "ivar_log.txt":[]
}




##################
#HELPER FUNCIONS
##################

def parse_arguments(parser:argparse.ArgumentParser):
    '''
    Verifies that arguments were passed to the script.
    If no arguments provided, prints help message and stops the script.
    Else returns argparse.Namespace object providing access the supplied argument values.
    '''
    if len(sys.argv)==1: #NO COMMAND-LINE ARGUMENTS PROVIDED
        parser.print_help(sys.stderr) #DISPLAY HELP
        sys.exit(1) # STOP THE SCRIPT
    return parser.parse_args() #NAMESPACE OBJECT WITH REFERENCES TO ARGUMENT VALUES


def check_date_format(date_string:str, valid_format:str="[0-9]{4}-[0-9]{2}-[0-9]{2}"):
    '''Used to verify that provided dates are in correct format.'''
    return re.match(valid_format, date_string)


def validate_arguments(args:argparse.Namespace):
    '''
    Given a namespace object containing values of parsed arguments, returns a dictionary of validated values.
    If validation fails, prints the error-indicating message to stdout and returns None.
    '''
    validated_arguments = {'skip_pango':args.skip_pango}
    if (args.start_date == None or args.end_date == None) and args.name_list == None: 
        #IF NO DATE AND NO SAMPLE ID LIST PROVIDED
        print(
            """
            If using dates, both d1 and d2 should be provided. 
            If one or both dates are missing, a csv file with sample id column should be provided.\n
            """
            )
        parser.print_help(sys.stderr) #PRINT HELP & STOP THE SCRIPT
        sys.exit(1)
    elif (args.start_date != None and args.end_date != None) and args.name_list != None: 
        #IF SAMPLE ID LIST AND BOTH SEQUENCING DATES ARE PROVIDED
        print(
            """
            WARNING: Both sample id list path and sequencing date range were provided. 
            Using sample id list to calculate mutation statistics.\n
            """
            )
        validated_arguments["id_list_path"] = args.name_list
        validated_arguments["metadata_file_path"] = args.metadata
        validated_arguments["date_1"] = None
        validated_arguments["date_2"] = None
    else: 
        #DATES OR LIST OF IDS ARE PROVIDED, NOT BOTH
        validated_arguments["id_list_path"] = args.name_list
        validated_arguments["metadata_file_path"] = args.metadata
        validated_arguments["date_1"] = args.start_date
        validated_arguments["date_2"] = args.end_date
    return validated_arguments


def validate_dates(valid_args:dict, date_format:str="%Y-%m-%d"):
    '''
    Given dictionary containing "date_1" and "date_2" as keys and date strings as values, checks if date format is corrects.
    If format is invalid, prints error message and stops the script.
    Else checks if dates represent valid time interval (e.g. date_1 <= date_2).
    If not, prints error message and stops the script.
    Else sets "date_1" and "date_2" values in dictionary to corresponding datetime objects.
    Returns None.
    '''

    checks = [check_date_format(valid_args["date_1"]), check_date_format(valid_args["date_2"])]
    if all(checks): #DATE FORMATS ARE VALID
        range_check = datetime.strptime(valid_args["date_1"], date_format) <= datetime.strptime(valid_args["date_2"], date_format)

        ##DATE ORDER CHECK
        if range_check:
            valid_args["date_1"] = datetime.strptime(valid_args["date_1"], date_format)
            valid_args["date_2"] = datetime.strptime(valid_args["date_2"], date_format)
        else:
            sys.exit(f'Start date is bigger than end date: {valid_args["date_1"]} > {valid_args["date_2"]}')

    else:  #FORMAT VALIDATION FAILED
        if not any(checks): #BOTH DATES FAILED
            sys.exit(f'Both dates are not of valid format (YYYY-MM-DD required): {valid_args["date_1"]} {valid_args["date_2"]}')
        elif not checks[0]: #date_1 FAILED
            sys.exit(f'End date is not of valied format (YYYY-MM-DD required): {valid_args["date_1"]}')
        else: #date_2 FAILED
            sys.exit(f'End date is not of valied format (YYYY-MM-DD required): {valid_args["date_2"]}')


def find_report_files(arguments:dict, walk_path:str, context_path_map:dict, pid_format:str="COV", pipeline_outdir_format:str = "NMRL", by_date:bool=True):
    '''
    Given dictionary of validated arguments, a path to a directory to look in (parent folder) and 
    context_path_map dictionary where every key represents a file extension to be included in the report,
    fills the context_path_map with paths to files if folder contains pipeline_outdir_format as substring and
    file matches the format and contains pid format as substring.  
    If by_date set to False, search by id list is performed. Returns None.
    '''
    if by_date:
        for (root,dirs,_) in os.walk(walk_path, topdown=True): #FINDING PATH TO EACH FOLDER THAT MATCHES TIME CRITERIA
            for dir in dirs: 
                if pipeline_outdir_format not in dir: #SCAN ONLY PROPERLY NAMED FOLDERS
                    continue          
                else: 
                    date = datetime.strptime(dir.split("-")[1],"%Y_%m_%d") #EXTRACT DATE FROM FOLDER NAME
                    if arguments['date_1'] <= date <= arguments['date_2']: #DATE IN SEARCH RANGE
                        for file in os.listdir(os.path.join(root,dir)): #CHECK FILES IN THAT FOLDER
                            for context in context_path_map: #AGAINST EVERY TYPE REQUIRED IN CONTEXT MAP
                                if context in file and pid_format in file: #TYPE AND ID IS MATCHED
                                    context_path_map[context].append(os.path.join(root,dir,file))
            break
    else:
        id_list = list(pd.read_csv(arguments["id_list_path"], header=None).iloc[:,0]) #CONVERT COLUMN OF IDS TO LIST
        for (root,dirs,_) in os.walk(covid_output_path, topdown=True): #FINDING PATH TO EACH FILE THAT MATCHES TIME CRITERIA
            for dir in sorted(dirs)[::-1]: 
                if pipeline_outdir_format not in dir: #SCAN ONLY PROPERLY NAMED FOLDERS
                    continue          
                else: 
                    for file in os.listdir(os.path.join(root,dir)):
                        for id in id_list:
                         #CHECK FILES IN THAT FOLDER
                            if f'_{id}_' in file or f'_{id}.' in file:
                                for context in context_path_map:
                                    if context in file:
                                        context_path_map[context].append(os.path.join(root, file))
                                break
                    if all(len(context_path_map[context])==len(id_list) for context in context_path_map): break
            break


def copy_files(file_path:str, source_files_path:str):
    '''Copy wrapper to use in multiprocessing: copy if not already there or different content'''
    #HTTPS://STACKOVERFLOW.COM/QUESTIONS/36821178/HOW-TO-SHUTIL-COPYFILE-ONLY-IF-FILE-DIFFER/36821211
    if not os.path.exists(f'{source_files_path}/{file_path.split("/")[-1]}') or not filecmp.cmp(file_path, f'{source_files_path}/{file_path.split("/")[-1]}'):
        shutil.copy(file_path, f'{source_files_path}/')


def copy_files_parallel(path_list:list, source_files_path:str, procs:int=6, progress_bar:bool=True):
    '''Function to copy files using multiprocessing'''
    with concurrent.futures.ProcessPoolExecutor(max_workers=procs) as executor:
        results = [executor.submit(copy_files, file_path, source_files_path) for file_path in path_list]
        processed_count = 0 #TO VIEW PROGRESS
        for _ in concurrent.futures.as_completed(results):
            processed_count += 1 #counting processed files
            if progress_bar: printProgressBar(processed_count, len(path_list), prefix = 'Progress:', suffix = 'Complete', length = 50)


def validate_context_file_paths(context_map:dict):
    '''
    Given dictionary mapping file format to full file path, returns None if all file formats are mapped to the same number of paths.
    Else returns a dict of context_map keys for which number of paths is less than the maximum number of paths found in context_map.
    '''
    max_path_count = len(max(list(context_map.values()), key=len))
    missing_files = {context:len(context_map[context]) for context in context_map if len(context_map[context]) != max_path_count}
    if not missing_files: return
    return missing_files


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    from https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


def readGenome(combined_fasta_path:str):
    '''
    Given path to multifasta, parses a file and returns python dictionary mapping each sequence header to the corresponding sequence.
    '''
    genome_dict = {}
    with open(combined_fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            if line[0] == '>': # IF LINE IS A HEADER (NEW HEADER IS REACHED)
                fasta_header = line # STORE HEADER UNTILE NEXT HEADER IS REACHED
                genome_dict[fasta_header] = '' # INITIALIZE AN EMPTY SEQUENCE MAPPED TO HEADER
            if line[0] != '>': # IF LINE IS A PART OF THE GENOME SEQUENCE
                genome_dict[fasta_header] += line.rstrip() # ADD THE LINE TO THE STRING THAT IS MAPPED TO THE LAST REACHED HEADER
    return genome_dict


def get_nt_counting_stats(combined_fasta_path:str, skip:bool=False, headers_changed:bool = True, column_list:list = ["receiving_lab_sample_id","genome_length","genome_N_percentage","genome_GC_content"]):
    '''
    Given path to multifasta, converts it to python dictionary.
    For each sequence in the dictionary, extracts sample id, sequence length, 
    N base content (%), GC-content (%). Returns the dataframe containing the 
    results for each sequence as rows and column_list contents as columns.
    '''

    if skip: #IF SKIPPING OPTION IS CHOSEN
        print("Statistics calculation option skipped.")
        return None #RETURN NONE FOR TESTING PURPOSES
    elif not skip: #IF SKIPPING OPTION IS NOT CHOSEN
        sequence_stats_df = pd.DataFrame(columns=column_list) #INITIALIZING DATAFRAME TO STORE THE RESULTING STATISTICS
        genome_dict_items = readGenome(combined_fasta_path).items() #GET HEADER:SEQUENCE PAIRS FROM READGENOME() OUTPUT
        for header, sequence in genome_dict_items: #LOOPING THROUGH EACH HEADER-SEQUENCE PAIR IN GENOME DICT
            sequence_length = len(sequence) #CALCULATING GENOME LENGTH
            if sequence_length > 0: 
                gc_content = round(100 * (sequence.count('G') + sequence.count('C')) / len(sequence), 2) #CALCULATING GC-CONTENT FOR EACH SEQUENCE
                n_content = round(100 * (sequence.count('N') / len(sequence)), 2)
            else: 
                gc_content = 0
                n_content = 0
            if re.search(r'[A-Z]{2}[0-9]{4}\.B([0-9]{2}|[0-9])', header): new_header = re.search(r'[A-Z]{2}[0-9]{4}\.B([0-9]{2}|[0-9])', header).group(0)
            elif headers_changed: new_header = header.split("/")[2] #EXTRACTING SAMPLE ID FROM FASTA HEADER TO BE USED LATER IN THE COLUMN MAPPING PROCESS
            else: new_header = header.split("_")[1].strip() #EXTRACTING SAMPLE ID FROM FASTA HEADER TO BE USED LATER IN THE COLUMN MAPPING PROCESS IF FASTA HEADERS WERE NOT CHANGED
            sequence_stats_df = sequence_stats_df.append({column_list[0]:new_header, column_list[1]:sequence_length, column_list[2]:n_content, column_list[3]:gc_content}, ignore_index=True, sort=False)
            #ADDING RESULT ROW TO THE DATAFRAME
        return sequence_stats_df # RETURN RESULTING DATAFRAME


def csv_info_extractor(file_path:str, filter_list:list, f_value:float):
    """
    Parses a csv file, exctract mutation information
    based on set of provided filters.
    """
    sample_id = file_path.split('/')[-1].split("_", 1)[1]
    sample_id = sample_id[:len(sample_id) - 8]
    processing_id = file_path.split('/')[-1][:9]
    result_row = {"SAMPLE_ID": sample_id, "processing_id": processing_id}
    try:
        df = pd.read_csv(file_path)
        for key in filter_list.keys(): 
            if use_f_filter:
                mut_df = df.loc[
                    (df["AMINO_ACID_CHANGE"].isin(filter_list[key])) & (df["FREQUENCY"] >= f_value)
                ]
            elif use_p_filter:
                mut_df = df.loc[
                    (df["AMINO_ACID_CHANGE"].isin(filter_list[key])) & (df["P_ERR_MUT_CALL"] < p_value)
                ]
            if len(mut_df["P_ERR_MUT_CALL"]) == len(filter_list[key]): #IF THE NUMBER OF EXTRACTED (UNIQUE) ROWS MATCHES THE NUMBER OF MUTATION_NAMES FOR A GIVEN FILTER_NAME
                result_row[key] = str(1) #ASSUME THAT THE MUTATION SPECIFIED BY THE FILTER WAS FOUND - ADD FILTER_NAME:1 PAIR TO THE RESULT_ROW DICT
            else: result_row[key] = str(round(len(mut_df["P_ERR_MUT_CALL"])/len(filter_list[key]), 2)) #IF THE NUMBER OF EXTRACTED (UNIQUE) ROWS DOES NOT MATCH THE NUMBER OF MUTATION_NAMES FOR A GIVEN FILTER_NAME - ADD FILTER_NAME:%MATCH PAIR TO THE RESULT_ROW DICT
    except:
        print(f'WARNING: failed to parse {file_path}')
        for key in filter_list.keys():
            result_row[key] = str(0)
    return result_row


def depth_info_extractor(file_path):
    '''
    Helper function to extract average and median coverage, 
    sequencing_date and sequencing lab for each sample.
    Returns result dictionary where sample is identified by SAMPLE_ID column.
    '''
    #ID EXTRACTION
    path_split = file_path.split('/')
    sample_id = path_split[-1].split("_", 1)[1]
    sample_id = sample_id[:len(sample_id) - 14]
    seq_date = path_split[-2].split("-")[1].replace("_","-")
    seq_lab = path_split[-2].split("-")[0].replace("_","(")+")"
    processing_id = file_path.split('/')[-1][:9]
    
    #BY-DEFAULT ASSUMING THAT COVERAGE FILE IS EMPTY (NO READS MAPPED - COVERAGE 0)
    result_row = {
        'SAMPLE_ID':sample_id, 
        "seq_date":seq_date, 
        'seq_institution':seq_lab, 
        "processing_id":processing_id, 
        'AVERAGE_COVERAGE':0, 
        'MEDIAN_COVERAGE':0
    }

    if not os.stat(file_path).st_size == 0:
        data = pd.read_csv(file_path, delimiter='\t').iloc[:,2]
        result_row['AVERAGE_COVERAGE'], result_row['MEDIAN_COVERAGE'] = sum(data)/len(data), data.median()
    return result_row
    

def mapped_info_extractor(file_path):
    '''
    Helper function to extract samtools flagstat report information for each sample.
    Returns result dictionary where sample is identified by SAMPLE_ID column.
    '''

    with open(file_path, "r+") as file: data = file.readlines()
    sample_id = file_path.split('/')[-1].split("_", 1)[1]
    sample_id = sample_id[:len(sample_id) - 18]
    processing_id = file_path.split('/')[-1][:9]
    try:
        result_row = {'SAMPLE_ID':sample_id, "processing_id":processing_id, 'TOTAL_READS':data[0].split(" ")[0], 'READS_MAPPED':data[4].split(" ")[0], 'MAPPED_FRACTION':round(int(data[4].split(" ")[0])/int(data[0].split(" ")[0]),2)}
    except IndexError:
        sys.exit(f'ERROR: empty file {file_path}')
    return result_row


def run_extractors_parallel(context_map:dict, filter_list:dict, f_value:float):
    '''
    Runs extractor functions using multiprocessing and assembles the results in single dataframe.
    Returns pandas dataframe. 
    Contect_map - dictionary that maps file types to file paths, 
    filter_list - dictionary that maps arbitrary name to mutation name, 
    f_value - frequency threshold to reduce false-negative mutation detections.
    '''

    #DETECTING SPECIFIC MUTATIONS FROM FILTERS
    result_df = pd.DataFrame() #INIT EMPTY DATAFRAME FOR THE REPORT
    with concurrent.futures.ProcessPoolExecutor() as executor: # APPLYING THE EXTRACTOR FUNCTION IN-PARALLEL ON DIFFERENT CORES
        results = [executor.submit(csv_info_extractor, file_path, filter_list, f_value) for file_path in context_map['ann.csv']] #SUBMITTING FUNCTION CALLS TO DIFFERENT PROCESSES
        processed_count = 0
        for f in concurrent.futures.as_completed(results): #COLLECTING PROCESSING RESULTS TO THE DATAFRAME
            result_df = result_df.append(f.result(), ignore_index=True) #ADD THE RESULT_ROW TO THE REPORT DF
            processed_count += 1
            printProgressBar(processed_count, len(context_map['ann.csv']), prefix = 'Extracting mutation stats:', suffix = 'Complete', length = 50)
    
#EXTRACTING COVERAGE DATA
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(depth_info_extractor, file_path) for file_path in context_map['seq_depth.txt']] #SUBMITTING FUNCTION CALLS TO DIFFERENT PROCESSES
        coverage_df = pd.DataFrame()
        processed_count = 0
        for f in concurrent.futures.as_completed(results): #COLLECTING PROCESSING RESULTS TO THE DATAFRAME
            coverage_df = coverage_df.append(f.result(), ignore_index=True) #ADD THE RESULT_ROW TO THE REPORT DF
            processed_count += 1
            printProgressBar(processed_count, len(context_map['seq_depth.txt']), prefix = 'Extracting coverage stats:', suffix = 'Complete', length = 50)
    result_df = pd.merge(result_df.applymap(str), coverage_df.applymap(str), how="left", on="SAMPLE_ID")

#EXTRACTING MAPPING STATISTICS
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(mapped_info_extractor, file_path) for file_path in context_map['mapped_report.txt']] #SUBMITTING FUNCTION CALLS TO DIFFERENT PROCESSES
        mapped_df = pd.DataFrame()
        processed_count = 0
        for f in concurrent.futures.as_completed(results): #COLLECTING PROCESSING RESULTS TO THE DATAFRAME
            mapped_df = mapped_df.append(f.result(), ignore_index=True) #ADD THE RESULT_ROW TO THE REPORT DF
            processed_count += 1
            printProgressBar(processed_count, len(context_map['mapped_report.txt']), prefix = 'Extracting mapping stats:', suffix = 'Complete', length = 50)

#ASSEMBLING AND FORMATTING RESULT DATAFRAME
    result_df = pd.merge(result_df.applymap(str), mapped_df.applymap(str), how="left", on="SAMPLE_ID")    
    column_reorder_list = ['SAMPLE_ID', 'processing_id'] + list(filter_list.keys()) + ['AVERAGE_COVERAGE', 'MAPPED_FRACTION', 'READS_MAPPED', 'TOTAL_READS', 'MEDIAN_COVERAGE', 'seq_date', 'seq_institution']
    result_df = result_df[column_reorder_list] #REORDERING COLUMNS BASED IN ORDER OF FILTERS IN REPORT_FILTERS.TXT FILE
    result_df.rename(columns={'SAMPLE_ID':"receiving_lab_sample_id"}, inplace=True) #FOR MERGING PURPOSES
    return result_df


def add_meta_pango(df:pd.DataFrame):
    '''Adds metadata to the raw result_df produced by extractor runner function. Returns resulting pandas dataframe.'''

    #ADDING METADATA AND PANGOLIN LINEAGES
    if all_excel_dump_df is not None: #IF METADATA MINING OPTION SELECTED - ADD METADATA TO THE MUTATION REPORT
        df = pd.merge(df.applymap(str), all_excel_dump_df.applymap(str), how="left", on="receiving_lab_sample_id")
    if sequence_stats_df is not None: #IF STATISTICS CALCULATION OPTION SELECTED - ADD SEQUENCE STATISTICS DATA TO THE MUTATION REPORT
        df = pd.merge(df.applymap(str), sequence_stats_df.applymap(str), how="left", on="receiving_lab_sample_id")
    if pango_report_path is not None: #IF PANGOLIN TYPING OPTION SELECTED - ADD PANGOLIN LINEAGE DATA TO THE MUTATION REPORT
        pango_report_df= pd.read_csv(pango_report_path) # READ PANGOLIN REPORT TO A DATAFRAME
        if re.match(r'[A-Z]{2}[0-9]{4}\.B', pango_report_df['taxon'][0]): pango_id_dict = {taxon:taxon for taxon in pango_report_df['taxon']}# MAP TAXOT TO SAMPLE ID IN A DICT
        elif rename_fasta_header: pango_id_dict = {taxon:taxon.split("/")[-2] for taxon in pango_report_df["taxon"]} # OPTION FOR GISAID-ACCEPTED HEADERS
        else: pango_id_dict = {taxon:taxon.split("_")[1] for taxon in pango_report_df["taxon"]} # MAP TAXON TO SAMPLE ID IN A DICT IF FASTA HEADERS WERE NOT CHANGED
        pango_report_df["receiving_lab_sample_id"] = pango_report_df["taxon"].map(pango_id_dict)
        pango_report_df = pango_report_df[["receiving_lab_sample_id", 'lineage']] #KEEP ONLY LINEAGE MATCHED AGAINST SAMPLE IDS
        df = pd.merge(df.applymap(str), pango_report_df.applymap(str), how="left", on="receiving_lab_sample_id")
    return df


def add_static_values(df:pd.DataFrame, 
    zero_fill_list:list = [
        'sequencing_notes', 
        'result_notes', 
        'sub_lineage'
    ], 
    static_fill_dict:dict = {
        'analysis_institution':analysis_institution, 
        'sequencing_platform':'illumina',
        'analysis_pipeline':pipeline_string
    }):
    '''
    Given list of columns that are to be filled with 0-s (zero_fill_list)
    and dictionary of columns that are to be filled with constants (static_fill_dict),
    fills the columns accordingly. In addition adds special values to rows defined in special case section of the function.
    Returns pandas dataframe.
    '''

    #ADDING STATIC VALUES
    static_fill_dict = {
        'analysis_institution':analysis_institution, 
        'sequencing_platform':'illumina',
        'analysis_pipeline':pipeline_string
    }

    for col_name in zero_fill_list: df[col_name] = np.zeros(len(df['receiving_lab_sample_id']))
    for col_name in static_fill_dict:
        df.loc[df.processing_id != 'Z_BMC', col_name] = static_fill_dict[col_name]    
        # result_df[col_name] = np.chararray(result_df['receiving_lab_sample_id'].shape, itemsize=len(static_fill_dict[col_name])+1).tostring()
    
    #SPECIAL CASE
    df.loc[df.processing_id == 'Z_BMC', 'sequencing_platform'] = 'MGI'
    df.loc[df.processing_id == 'Z_BMC', 'analysis_institution'] = 'BMC'
    return df


def parse_directory_tree(df:pd.DataFrame, report_path:str, raw_folder_path:str):
    '''
    Given path to the report folder and raw folder, parses the directory structure to extract relevant metadata.
    Adds metadata to the df based on receiving_lab_sample_id. Returns tuple of 3 pandas dataframes in order:
    df - input dataframe annotated with analysis_date & analysis batch; nmrl_samples - dataframe where each analysed
    sample that was sequenced in NMRL is mapped to its fastq file path (one for each sample); eurofins_samples - 
    same as nmrl_samples but for samples sequenced by eurofins. 
    '''

    # ANALYSIS_DATE & ANALYSIS_BATCH_ID
    analysis_batch_id = report_path.split('/')[-1]
    analysis_date = analysis_batch_id.split('_')[1]
    df['analysis_date'] = np.chararray(df['receiving_lab_sample_id'].shape, itemsize=len(analysis_date)+1).tostring()
    df['analysis_batch_id'] = np.chararray(df['receiving_lab_sample_id'].shape, itemsize=len(analysis_batch_id)+1).tostring()
    df.loc[df.processing_id != 'Z_BMC', ['analysis_date','analysis_batch_id']] = [analysis_date,analysis_batch_id]

    #PREPARING FILES TO USE BASH UTILITIES
    sample_ids = df['receiving_lab_sample_id']
    sample_ids.to_csv(f'{report_path}/sample_ids.csv',header=False,index=False)

    #FINDNG FULL PATHS TO FASTQ FILES USING BASH
    if not os.path.isfile(f'{report_path}/path_list.csv') or os.stat(f'{report_path}/path_list.csv').st_size == 0:
        os.system(f'ls -R {raw_folder_path} | grep 1.fastq.gz > {report_path}/fastq_files.csv')
        os.system(f'cat {report_path}/sample_ids.csv | while IFS="," read a ; do cat {report_path}/fastq_files.csv | grep "$a" ; done > {report_path}/sample_file_list.csv')
        os.system(f'cat {report_path}/sample_file_list.csv | xargs -n1 -P12 -I% find {raw_folder_path}/ -type f -name % > {report_path}/path_list.csv')
        os.system(f'rm {report_path}/sample_file_list.csv {report_path}/fastq_files.csv {report_path}/sample_ids.csv')
    sample_paths = pd.read_csv(f'{report_path}/path_list.csv', header=None)
    sample_frame = pd.DataFrame({'id':[],'path':[]})

    #MAPPING RUN FOLDER TO SAMPLE ID
    for id in sample_ids:
        paths = sample_paths.loc[sample_paths[sample_paths.columns[0]].str.contains(id)].reset_index() #GET LIST OF PATHS THAT CONTAIN SAMPLE ID
        paths.rename(columns={sample_paths.columns[0]:'path'}, inplace=True) #PROVIDE CORRECT COLUMN NAME
        paths["id"] = [id for _ in range(len(paths.path))] #ADD ID
        paths.drop('index', inplace=True, axis=1) #REMOVE INDEX COLUMN
        sample_frame = sample_frame.append(paths) #COMBINE IN ONE DF
    # os.system(f'rm {report_path}/path_list.csv')

    #GETTING RUN FOLDER NAMES FROM FULL PATHS
    sample_frame[['raw','run_folder','fastq']] = sample_frame['path'].str.rsplit('/', n=2, expand=True) #GENERATE NEEDED COLUMNS
    sample_frame.drop('raw', inplace=True, axis=1) #REMOVE USELESS COLUMN
    sample_frame.drop('fastq', inplace=True, axis=1) #REMOVE USELESS COLUMN
    eurofins_samples = sample_frame[~sample_frame['run_folder'].str.contains('nmrl')]
    nmrl_samples = sample_frame[sample_frame['run_folder'].str.contains('nmrl')]
    # sample_frame.to_csv('/mnt/home/jevgen01/nmrl/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/sample_frame.csv', header=True, index=False)
    #SPLIT RUN FOLDER NAME TO COLUMNS AND DROP SAMPLE DUPLICATES
    if len(eurofins_samples) > 0:
        eurofins_samples[['seq_date','1','2','3','4']] = eurofins_samples['run_folder'].str.rsplit('-', n=4, expand=True)
    if len(nmrl_samples) > 0:
        nmrl_samples[['seq_date','lab','kit','instrument', 'seq_mode', 'primer']] = nmrl_samples['run_folder'].str.rsplit('-', n=5, expand=True)
    eurofins_samples.drop_duplicates(subset='id', keep='first', inplace=True)
    nmrl_samples.drop_duplicates(subset='id', keep='first', inplace=True)
    # nmrl_samples.to_csv('/mnt/home/jevgen01/nmrl/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/nmrl_samples.csv', header=True, index=False)
    
    return df, eurofins_samples, nmrl_samples


def add_tree_info(result_df:pd.DataFrame, nmrl_samples:pd.DataFrame, eurofins_samples:pd.DataFrame):
    '''Given result dataframe, dataframe obtained from paths to samples sequenced in NMRL and
    dataframe obtained from paths to samples sequences in eurofins, performs data formatting and metadata extraction.
    Returns result dataframe annotated with run-related metadata.
    '''

    #INIT EMPTY COLUMNS TO STORE INFORMATION
    blanks = ["used_sequencing_run_ids", "used_batch_ids","used_batch_ids", "library_prep_method", "analysis_pipeline_notes"]
    for col_name in blanks: result_df[col_name] = [None for _ in range(len(result_df))]

    #ADDING METADATA TO EACH SAMPLE ID IN REPORT
    for id in result_df['receiving_lab_sample_id']:
        reported_samples = result_df.loc[result_df['receiving_lab_sample_id'] == id] #GET ALL ROWS THAT CONTAIN THIS ID
        if len(reported_samples) == 1:  #IF ID IS UNIQUE
            reported_samples = reported_samples.reset_index(drop=True) #TO GET 0-BASED INDEXES

    #GET METADATA FOR THE SAMPLE
            if reported_samples['seq_institution'][0] == 'NMRL(LIC)':  #IF SEQUENCING INSTITUTION IS DETERMINED AS NMRL
                id_match = nmrl_samples.loc[nmrl_samples.id == id].reset_index(drop=True)
                run_folder = id_match['run_folder'][0]
                batch_id = id_match['seq_date'][0]
                used_batch_ids = batch_id
                library_prep_method = id_match['kit'][0]
                analysis_pipeline_notes = id_match['primer'][0]
            else:  #IF NOT NMRL
                id_match = eurofins_samples.loc[eurofins_samples.id == id].reset_index(drop=True)
                run_folder = id_match['run_folder'][0]
                batch_id = id_match['seq_date'][0]
                library_prep_method = 'eurofins in-house'
                used_batch_ids = f'Eurofins_{batch_id}'
                analysis_pipeline_notes = 'arctic v4'

    #ADD RESULTS TO THE DATAFRAME (ID CAN BE USED AS UNIQUE ROW IDENTIFIER)
            result_map = {
                'used_batch_ids':used_batch_ids,
                'used_sequencing_run_ids':run_folder,
                'library_prep_method':library_prep_method,
                'analysis_pipeline_notes':analysis_pipeline_notes
                }
            for col_name in result_map: result_df.loc[result_df.receiving_lab_sample_id == id, col_name] = result_map[col_name]

#IF ID IS NOT UNIQUE (E.G. MORE THAN ONE SAMPLE WITH THE SAME ID SEQUENCED BY DIFFERENT/SAME INSTITUTIONS) (TO BE TESTED)
        elif len(reported_samples) > 1:
            # reported_samples.to_csv("/mnt/home/jevgen01/nmrl/cov_analysis/SARS-CoV2_assembly/subscripts/downstream/reported_samples.csv", header=True, index=False)
            for i in reported_samples.index: #ROW INDEX IS USED AS UNIQUE ROW IDENTIFIER (INDECES OF RESULT_DF ARE KEPT)
                if reported_samples.iloc[i]['seq_institution'] == 'NMRL(LIC)':   #IF SEQUENCING INSTITUTION IS DETERMINED AS NMRL
                    run_folder = nmrl_samples.loc[nmrl_samples.id == id]['run_folder'][0]
                    batch_id = nmrl_samples.loc[nmrl_samples.id == id]['seq_date'][0]
                    used_batch_ids = f'{batch_id}'
                    library_prep_method = nmrl_samples.loc[nmrl_samples.id == id]['kit'][0]
                    analysis_pipeline_notes = 'arctic v3'
                else: #IF NOT NMRL
                    run_folder = eurofins_samples.loc[eurofins_samples.id == id]['run_folder'][0]
                    batch_id = eurofins_samples.loc[eurofins_samples.id == id]['seq_date'][0]
                    library_prep_method = 'eurofins in-house'
                    used_batch_ids = f'Eurofins_{batch_id}'
                    analysis_pipeline_notes = 'arctic v4'
#ADD RESULTS TO THE DATAFRAME (INDEX CAN BE USED AS UNIQUE ROW IDENTIFIER)
                result_map = {
                'used_batch_ids':used_batch_ids,
                'used_sequencing_run_ids':run_folder,
                'library_prep_method':library_prep_method,
                'analysis_pipeline_notes':analysis_pipeline_notes
                }
                for col_name in result_map: result_df.at[result_df.index == i, col_name] = result_map[col_name]
    return result_df


def mutation_report_generator(output_path:str, f_value:float):
    """
    Generate mutation report using helper functions to extract and process metadata and pipeline results.
    """

#EXTRACTING DATA FROM PIPELINE REPORTS, ADDING METADATA AND PANGO LINEAGES
    result_df = run_extractors_parallel(context_path_map, filter_list, f_value)
    result_df = add_meta_pango(result_df)
    result_df.drop_duplicates(subset=['receiving_lab_sample_id'],inplace=True)

#ADDING STATIC VALUES
    result_df = add_static_values(result_df)

#EXTRACTING METADATA FROM DIRECTORY STRUCTURE
    tree_data = parse_directory_tree(result_df, report_path, raw_folder_path)
    result_df, eurofins_samples, nmrl_samples = tree_data[0], tree_data[1], tree_data[2]

#ADDING EXTRACTED INFO TO THE RESULTS AND GENERATING REPORT
    result_df = add_tree_info(result_df, nmrl_samples, eurofins_samples)
    result_df.to_csv(f'{output_path}/pipeline_report.csv', header=True, index=False, encoding='utf-8-sig') #GENERATING THE REPORT




#####################
#RUNNING THE SCRIPT
#####################

if __name__ == "__main__":
    
#ARGUMENT DECLARATION
    parser = argparse.ArgumentParser(description='A script to generate a mutation report given range of dates or list of sample ids.') #ARGPARSER OBJECT TO PROVIDE COMMAND-LINE FUNCTIONALITY
    parser.add_argument('-d1', '--start_date', metavar='\b', help = 'A starting sequencing date of the reference interval (YYYY-MM-DD).', default=None, required=False)
    parser.add_argument('-d2', '--end_date', metavar='\b', help = 'An ending sequencing date of the reference interval (YYYY-MM-DD).', default=None, required=False)
    parser.add_argument('-l', '--name_list', metavar='\b', help = 'Path to list of file names to lookup in folders', default=None, required=False)
    parser.add_argument('-o', '--output_dir', metavar='\b', help = 'Path to the folder where report folder should be created', default=report_folder_path, required=False)
    parser.add_argument('-m', '--metadata', metavar='\b', help = 'Path to the metadata table', default=default_metadata_path, required=False)
    parser.add_argument('-s', '--skip_pango', help = "Flag to skip pangolin typing.", action='store_true')

#ARGUMENT PRE-ROCESSING
    args = parse_arguments(parser)
    valid_args = validate_arguments(args)
    if valid_args["date_1"]: validate_dates(valid_args)


#SEARCHING FOR RELEVANT FILES & FOLDERS
    if valid_args["date_1"]: find_report_files(valid_args, covid_output_path, context_path_map=context_path_map) #SEARCHING BY DATE RANGE
    else: find_report_files(valid_args, covid_output_path, context_path_map=context_path_map, by_date=False) #SEARCHING BY SAMPLE ID


#VERIFYING PROCESSING COMPLETION
    file_count_check = validate_context_file_paths(context_path_map)
    if file_count_check: #PRINT THE FILE FORMATS THAT FAILED WITH CORRESPONDING FILE COUNT
        [print(context, file_count_check[context]) for context in file_count_check]
        sys.exit('\nERROR: Some samples lack input files!\nPlease reprocess the samples or modify the query so that all required files exist.')


#CREATING REPORT FOLDER AND COPYING FILES
    if valid_args['id_list_path'] is not None: #CREATE NEW REPORT_FOLDER
        source_files_path = f'{report_folder_path}/report_{datetime.now().date().strftime("%Y-%m-%d")}_{valid_args["id_list_path"]}/source_files' #PATH TO FOLDER WHERE SOURCE FILES SHOULD BE SAVED
        report_path = str(pathlib.Path(source_files_path).parents[0]) #PATH TO FOLDER WHERE REPORTS SHOULD BE SAVED
        pathlib.Path(source_files_path).mkdir(parents=True,exist_ok=True)
    else:
        source_files_path = f'{report_folder_path}/report_{datetime.now().date().strftime("%Y-%m-%d")}_{valid_args["date_1"].strftime("%Y-%m-%d")}_{valid_args["date_2"].strftime("%Y-%m-%d")}/source_files'
        report_path = str(pathlib.Path(source_files_path).parents[0]) #PATH TO FOLDER WHERE REPORTS SHOULD BE SAVED
        pathlib.Path(source_files_path).mkdir(parents=True, exist_ok=True)


    #COPY FILES REQUIRED FOR REPORT TO THE SOURCE_FILES FOLDER UNDER CREATED REPORTS FOLDER
    for context in context_path_map:
        print(f'Copying {context} files:')
        copy_files_parallel(context_path_map[context],source_files_path)


##################
#DATA PROCESSING
##################

#READING FILTER FILE
    with open(filter_path, "r+") as filter_file:
        #SAVES FILTER_NAME:FILTER_SET PAIR TO A DICTIONARY
        filter_list = {filter.split(" ")[0]:set(filter.split(" ")[1].split(",")) for filter in filter_file.read().strip().split("\n")} 

#READING METADATA IF ALLOWED BY CONTROLLERS
    if not skip_excel_mining:
        all_excel_dump_df = pd.read_csv(valid_args["metadata_file_path"])
    else: all_excel_dump_df = None

# RUNNING PANGOLIN IF ALLOWED BY CONTROLLERS
    if not valid_args['skip_pango']: # PROVIDING CORRECT PATH FOR RUN_PANGOLIN.SH & USING SUBPROCESS TO RUN IT WITH OPTIONS FROM CONTROL VARIABLE
        #RUNNING THE COMMAND
        print(f'Submitting pangolin typing job to HPC.')
        log_path = f'{log_folder_path}{datetime.now().strftime("%Y-%m-%d-%H-%M")}_pangolin.log'
        command = ['qsub', "-o", log_path, "-e", log_path, "-F", f"{1} {source_files_path}", f"{subprocess_path}run_pangolin.sh"]
        subprocess.check_call(command)
        #WAITING FOR JOB TO COMPLETE
        while not os.path.isfile(f'{source_files_path}/{timestr_fasta}_lineage_report.csv'): time.sleep(30) #whould be a good place for asyncio
        #PATHS TO REPORTS
        pango_report_path = f'{source_files_path}/{timestr_fasta}_lineage_report.csv'
        combined_path = [path for path in os.listdir(source_files_path) if 'combined.fasta' in path][0]
        combined_fasta_path = f'{source_files_path}/{combined_path}'
    elif valid_args['skip_pango']:
        print(f'Pangolin typing skipped.')
        pango_report_path = f'{report_path}/{timestr_fasta}_lineage_report.csv'
        combined_path = [path for path in os.listdir(source_files_path) if 'combined.fasta' in path][0]
        combined_fasta_path = f'{source_files_path}/{combined_path}'



#GENERATING MUTATION REPORT
    if not skip_sequence_stats: sequence_stats_df = get_nt_counting_stats(combined_fasta_path, skip_sequence_stats, rename_fasta_header) 
    mutation_report_generator(source_files_path, f_value)

#PERFORMING CLEANUP
    os.system(f'mv {source_files_path}/pipeline_report.csv {report_path}')
    if not valid_args['skip_pango']: os.system(f'mv {source_files_path}/{timestr_fasta}_lineage_report.csv {report_path}')
