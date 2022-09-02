#IMPORTS
import sys, os, pandas as pd, re, time, datetime, concurrent.futures, subprocess, argparse, shutil, pathlib, filecmp, numpy as np
from datetime import datetime
from pathlib import Path

parser = argparse.ArgumentParser(description='A script to generate a mutation report given range of dates or list of sample ids.') #ARGPARSER OBJECT TO PROVIDE COMMAND-LINE FUNCTIONALITY
parser.add_argument('-d1', '--start_date', metavar='\b', help = 'A starting sequencing date of the reference interval (YYYY-MM-DD).', default=None, required=False)
parser.add_argument('-d2', '--end_date', metavar='\b', help = 'An ending sequencing date of the reference interval (YYYY-MM-DD).', default=None, required=False)
parser.add_argument('-l', '--name_list', metavar='\b', help = 'Path to list of file names to lookup in folders', default=None, required=False)
parser.add_argument('-o', '--output_dir', metavar='\b', help = 'Path to the folder where report folder should be created', default='/home/groups/nmrl/cov_analysis/reports', required=False)
parser.add_argument('-m', '--metadata', metavar='\b', help = 'Path to the metadata table', default='/home/groups/nmrl/cov_analysis/reports', required=True)
parser.add_argument('-q', '--multiqc', help = 'Export contents of summary file in the SISdb format', action="store_true")

#PATHS & SETUPS
subprocess_path = f'/home/groups/nmrl/cov_analysis/fastq_processing/subprocess'
report_folder_path = f'/mnt/home/groups/nmrl/cov_analysis/reports'
raw_folder_path = f'/mnt/home/groups/nmrl/cov_analysis/raw'
timestr_fasta = time.strftime("%m_%d_%Y")
filter_path = f'/home/groups/nmrl/cov_analysis/fastq_processing/' # DOWNSTREAM DIRECTORY
f_value = 65 #SET MAX ACCEPTABLE FREQUENCY VALUE TO CONFIRM THAT MUTATION IS DETECTED
p_value = 0.05 #SET MAX ACCEPTABLE P VALUE

#CONTROLLERS - DEV MODE
skip_excel_mining = False
run_pangolin = True
skip_sequence_stats = False
rename_fasta_header = True
use_p_filter = False
use_f_filter = True
request_new = True

#PARSING ARGUMENTS
if len(sys.argv)==1: #IF NO COMMAND-LINE ARGUMENTS PROVIDED - DISPLAY HELP AND STOP SCRIPT EXCECUTION
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args() #ARGS LIST FROM COMMAND-LINE INPUT

#VALIDATING ARGUMENTS
if (args.start_date == None or args.end_date == None) and args.name_list == None: #IF NO DATE AND NO SAMPLE ID LIST PROVIDED
    print("If using dates, both d1 and d2 should be provided. If one or both date is missing, a csv file with sample ids should be provided.")
    parser.print_help(sys.stderr)
    sys.exit(1)

date_1 = args.start_date
date_2 = args.end_date
id_list_path = args.name_list
metadata_file_path= args.metadata

if (args.start_date != None and args.end_date != None) and args.name_list != None: #IF SAMPLE ID LIST AND BOTH SEQUENCING DATES ARE PROVIDED
    print("WARNING: Both sample id list path and sequencing date range were provided. Using sample id list to calculate mutation statistics.")
    date_1 = None
    date_2 = None

def check_date(date_string):
    '''The function is used to verify that provided dates are in correct format.'''
    check = True if re.match("[0-9]{4}-[0-9]{2}-[0-9]{2}", date_string) else False
    return check  

#STRINGS TO LOOK FOR
processing_id = "COV"
mutation_file = "ann.csv"
coverage_file = 'seq_depth.txt'
mapping_file = 'mapped_report.txt'
consensus_file = 'consensus.fasta'
fastp_file = 'fastp.json'
cutadapt_file = 'cutadapt_log.txt'
ivar_file = 'ivar_log.txt'
screen= "R1_001_screen.html"
inhouse_header = 'COV'

#EXTRACTING PATHS TO FILES
csv_list = []
coverage_list = []
mapped_list = []
fasta_list = []
ivar_list = []
cutadapt_list = []
fastp_list = []
if request_new:
    screen_list = []

if date_1 != None: #IF SELECTING BY DATE RANGE
    if check_date(date_1) and check_date(date_2):
        start = datetime.strptime(date_1, "%Y-%m-%d")
        end = datetime.strptime(date_2, "%Y-%m-%d")
        if start <= end:
            for (root,dirs,files) in os.walk("/home/groups/nmrl/cov_analysis/covid_output/", topdown=True): #FINDING PATH TO EACH FILE THAT MATCHES TIME CRITERIA
                for dir in dirs: 
                    if "NMRL" not in dir: #SCAN ONLY PROPERLY NAMED FOLDERS
                        continue                    
                    else: 
                        date = datetime.strptime(dir.split("-")[1],"%Y_%m_%d")
                        if start <= date <= end:
                            for file in os.listdir(os.path.join(root,dir)):
                                if mutation_file in file and processing_id in file:
                                    csv_list.append(os.path.join(root,dir,file))
                                elif coverage_file in file and processing_id in file:
                                    coverage_list.append(os.path.join(root,dir,file))
                                elif mapping_file in file and processing_id in file:
                                    mapped_list.append(os.path.join(root,dir,file))
                                elif consensus_file in file and processing_id in file:
                                    fasta_list.append(os.path.join(root,dir,file))
                                elif cutadapt_file in file and processing_id in file:
                                    cutadapt_list.append(os.path.join(root,dir,file))
                                elif ivar_file in file and processing_id in file:
                                    ivar_list.append(os.path.join(root,dir,file))
                                elif fastp_file in file and processing_id in file:
                                    fastp_list.append(os.path.join(root,dir,file))
                                if request_new:
                                    if screen in file:
                                        screen_list.append(os.path.join(root,dir,file))
        else:
            sys.exit(f'Start date is bigger than end date: {date_1} > {date_2}')

    else: 
        if not check_date(date_1):
            sys.exit(f'Start date is not of valid format (YYYY-MM-DD required): {date_1}')
        else:
            sys.exit(f'End date is not of valied format (YYYY-MM-DD required): {date_2}')

else: #IF SELECTING BY SAMPLE LIST
    id_list = list(pd.read_csv(id_list_path, header=None).iloc[:,0]) #CONVERT COLUMN OF IDS TO LIST
    for (root,dirs,files) in os.walk("/home/groups/nmrl/cov_analysis/covid_output/", topdown=True): #FINDING PATH TO EACH FILE THAT MATCHES TIME CRITERIA
            for file in files:
                if mutation_file in file and processing_id in file:
                    if inhouse_header in file: #IF INHOUSE ID IN FILE NAME
                       if file.split(".",1)[0][10:] in id_list:
                            csv_list.append(os.path.join(root, file))
                    else: #IF INHOUSE ID NOT IN FILE NAME
                        if file.split(".",1)[0] in id_list:
                            csv_list.append(os.path.join(root, file))
                elif coverage_file in file and processing_id in file:
                    if inhouse_header in file: #IF INHOUSE ID IN FILE NAME
                       if file[10:len(file)-14] in id_list:
                            coverage_list.append(os.path.join(root, file))
                    else: #IF INHOUSE ID NOT IN FILE NAME
                        if file[:len(file)-14] in id_list:
                            coverage_list.append(os.path.join(root, file))
                elif mapping_file in file and processing_id in file:
                    if inhouse_header in file: #IF INHOUSE ID IN FILE NAME
                       if file[10:len(file)-18] in id_list:
                            mapped_list.append(os.path.join(root, file))
                    else: #IF INHOUSE ID NOT IN FILE NAME
                        if file[:len(file)-18] in id_list:
                            mapped_list.append(os.path.join(root, file))
                elif consensus_file in file and processing_id in file:
                    if inhouse_header in file: #IF INHOUSE ID IN FILE NAME
                        if file[10:len(file)-16] in id_list:
                            fasta_list.append(os.path.join(root, file))
                    else: #IF INHOUSE ID NOT IN FILE NAME
                        if file[:len(file)-16] in id_list:
                            fasta_list.append(os.path.join(root, file))
                elif cutadapt_file in file and processing_id in file:
                    if inhouse_header in file: #IF INHOUSE ID IN FILE NAME
                       if file[10:len(file)-17] in id_list:
                            cutadapt_list.append(os.path.join(root, file))
                    else: #IF INHOUSE ID NOT IN FILE NAME
                        if file[:len(file)-17] in id_list:
                            cutadapt_list.append(os.path.join(root, file))
                elif ivar_file in file and processing_id in file:
                    if inhouse_header in file: #IF INHOUSE ID IN FILE NAME
                       if file[10:len(file)-13] in id_list:
                            ivar_list.append(os.path.join(root, file))
                    else: #IF INHOUSE ID NOT IN FILE NAME
                        if file[:len(file)-13] in id_list:
                            ivar_list.append(os.path.join(root, file))
                elif fastp_file in file and processing_id in file:
                    if inhouse_header in file: #IF INHOUSE ID IN FILE NAME
                        if file[10:len(file)-11] in id_list:
                            fastp_list.append(os.path.join(root, file))
                    else: #IF INHOUSE ID NOT IN FILE NAME
                        if file[:len(file)-11] in id_list:
                            fastp_list.append(os.path.join(root, file))
                if request_new:
                    if screen in file: #FOR FASTQSCREEN FILES
                        if file[:2] == "C0":
                            if file.split("_",3)[-2] in id_list:
                                screen_list.append(os.path.join(root,dir,file))
                        else:
                            if file.split("_",1)[0] in id_list:
                                screen_list.append(os.path.join(root,dir,file))
                        screen_length = len(screen_list)

#VERIFYING DATA INTEGRITY
if not request_new:
    screen_length = len(csv_list)
else:
    screen_length = len(screen_list)

if not len(csv_list) == len(coverage_list) == len(mapped_list) == len(fasta_list) == screen_length: #CHECKING IF SOME SAMPLES LACK ANY OF THE REPORT FILES
    print(f'WARNING:!!! some samples lack one or more report file(s) !!! Details below:')
    id_list_csv = [path.split("/")[-1].split(".",1)[0][10:] if inhouse_header in path else path.split("/")[-1].split(".",1)[0] for path in csv_list]
    id_list_fasta = [path.split("/")[-1][10:len(path.split("/")[-1])-16] if inhouse_header in path else path.split("/")[-1][:len(path.split("/")[-1])-16] for path in fasta_list]
    id_list_coverage = [path.split("/")[-1][10:len(path.split("/")[-1])-14] if inhouse_header in path else path.split("/")[-1][:len(path.split("/")[-1])-14] for path in coverage_list]
    id_list_mapped = [path.split("/")[-1][10:len(path.split("/")[-1])-18] if inhouse_header in path else path.split("/")[-1][:len(path.split("/")[-1])-18] for path in mapped_list]
    id_list_fastp = [path.split("/")[-1][10:len(path.split("/")[-1])-11] if inhouse_header in path else path.split("/")[-1][:len(path.split("/")[-1])-11] for path in fastp_list]
    id_list_cutadapt = [path.split("/")[-1][10:len(path.split("/")[-1])-17] if inhouse_header in path else path.split("/")[-1][:len(path.split("/")[-1])-17] for path in cutadapt_list]
    id_list_ivar = [path.split("/")[-1][10:len(path.split("/")[-1])-13] if inhouse_header in path else path.split("/")[-1][:len(path.split("/")[-1])-13] for path in ivar_list]
    if request_new:
        id_list_screen = list(set(path.split("/")[-1].split("_",3)[-2] if path.split("/")[-1][:2] == "CO" else path.split("/")[-1].split("_",1)[0] for path in screen_list))
        lists = [id_list_csv, id_list_fasta, id_list_coverage, id_list_mapped, id_list_fastp, id_list_cutadapt, id_list_ivar, id_list_screen]
    else:
        lists = [id_list_csv, id_list_fasta, id_list_coverage, id_list_mapped, id_list_fastp, id_list_cutadapt, id_list_ivar]
    max_len = max([len(i) for i in lists])
    if request_new:
        lists_to_check = {mutation_file:id_list_csv,consensus_file:id_list_fasta,coverage_file:id_list_coverage,mapping_file:id_list_mapped, ivar_file:id_list_ivar, cutadapt_file:id_list_cutadapt, fastp_file:id_list_fastp, screen:id_list_screen}
    else:
        lists_to_check = {mutation_file:id_list_csv,consensus_file:id_list_fasta,coverage_file:id_list_coverage,mapping_file:id_list_mapped, ivar_file:id_list_ivar, cutadapt_file:id_list_cutadapt, fastp_file:id_list_fastp}
    template_list = None
    for lst in lists: #FINDING LIST WITH MOST SAMPLE IDS INCLUDED
        if len(lst) == max_len:
            template_list = lst
            break
    for lst_1 in lists_to_check.keys(): #CHECKING IF LISTS MATCH - SORTING ENFORCES EQUALITY OF INDEXES
        if sorted(lists_to_check[lst_1]) == sorted(template_list): #IF CORRESPONDING ELEMENTS OF THE LISTS ARE EQUAL & LENGTHS ARE EQUAL
            continue
        else:
            for id in sorted(template_list): #IF A QUERY LIST LACKS SOME ID - CORRESPONDING FILE IS MISSING
                if not id in lists_to_check[lst_1]:
                    print(f'\tSample {id}: {lst_1} file not found!')

    sys.exit('\nERROR: Some samples lack input files!\n\n\tPlease reprocess the samples or modify the query so that all required files exist.\n')

#COPY FASTA FILES TO THE REPORT FOLDER
if args.name_list is not None: #CREATE NEW REPORT_FOLDER
    source_files_path = f'{report_folder_path}/report_{"{0:%Y-%m-%d}".format(datetime.now())}_{args.name_list}/source_files' #PATH TO FOLDER WHERE SOURCE FILES SHOULD BE SAVED
    report_path = str(pathlib.Path(source_files_path).parents[0]) #PATH TO FOLDER WHERE REPORTS SHOULD BE SAVED
    pathlib.Path(source_files_path).mkdir(parents=True,exist_ok=True)
else:
    source_files_path = f'{report_folder_path}/report_{"{0:%Y-%m-%d}".format(datetime.now())}_{date_1}_{date_2}/source_files'
    report_path = str(pathlib.Path(source_files_path).parents[0]) #PATH TO FOLDER WHERE REPORTS SHOULD BE SAVED
    pathlib.Path(source_files_path).mkdir(parents=True, exist_ok=True)


def copy_files(file_path):
    '''Copy wrapper to use in multiprocessing: copy if not already there or different content'''
    #HTTPS://STACKOVERFLOW.COM/QUESTIONS/36821178/HOW-TO-SHUTIL-COPYFILE-ONLY-IF-FILE-DIFFER/36821211
    if not os.path.exists(f'{source_files_path}/{file_path.split("/")[-1]}') or not filecmp.cmp(file_path, f'{source_files_path}/{file_path.split("/")[-1]}'):
        shutil.copy(file_path, f'{source_files_path}/')
        return file_path.split("/")[-1],True
    else:
        return file_path.split("/")[-1],False

def copy_files_parallel(path_list, file_type):
    '''Function to copy files of given type using multiprocessing'''
    with concurrent.futures.ProcessPoolExecutor() as executor:
        start_time = time.time() #TO ESTIMATE THE TIME IT TOOK TO PROCESS THE FILES
        results = [executor.submit(copy_files, file_path) for file_path in path_list] 
        for f in concurrent.futures.as_completed(results):
            result = f.result()
            if result[1]:
                print(f'{result[0]} - COPIED TO {source_files_path} - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.now())}') #PRINTING THE TIME WHEN THE PROCESSING WAS FINISHED FOR A GIVEN SAMPLE
            elif not result[1]:
                print(f'{result[0]} - FOUND IN {source_files_path} - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.now())}') #PRINTING THE TIME WHEN THE PROCESSING WAS FINISHED FOR A GIVEN SAMPLE
        print(f"\nCopy finished for {len(path_list)} {file_type} in --- %s seconds ---\n" % round(time.time() - start_time, 2))

#USING MULTIPROCESSING TO COPY FASTA FILES TO THE REPORT DIR
copy_files_parallel(fasta_list, consensus_file)
copy_files_parallel(coverage_list, coverage_file)
copy_files_parallel(mapped_list, mapping_file)
copy_files_parallel(csv_list, mutation_file)
copy_files_parallel(ivar_list, ivar_file)
copy_files_parallel(cutadapt_list, cutadapt_file)
copy_files_parallel(fastp_list, fastp_file)
if request_new:
    copy_files_parallel(screen_list, screen)
    
# READING METADATA IF ALLOWED BY CONTROLLERS
if not skip_excel_mining:
    all_excel_dump_df = pd.read_csv(metadata_file_path)
else: all_excel_dump_df = None

# RUNNING PANGOLIN IF ALLOWED BY CONTROLLERS
if run_pangolin: # PROVIDING CORRECT PATH FOR RUN_PANGOLIN.SH & USING SUBPROCESS TO RUN IT WITH OPTIONS FROM CONTROL VARIABLE
    os.chdir(subprocess_path), subprocess.check_call(["./run_pangolin.sh", str(1), source_files_path]) 
elif not run_pangolin: # IF PANGOLIN TYPING OPTION SKIPPED - CHECK IF PANGOLIN UDPATE IS REQUIRED
    os.chdir(subprocess_path), subprocess.check_call(["./run_pangolin.sh", str(0), source_files_path]) 
# PATHS TO PANGOLIN OUTPUT
combined_path = [path for path in os.listdir(source_files_path) if 'combined.fasta' in path][0]
combined_fasta_path = f'{source_files_path}/{combined_path}'
pango_report_path = f'{source_files_path}/{timestr_fasta}_lineage_report.csv'

#READING FILTER FILE
with open(f"{filter_path}/resources/report_filters.txt", "r+") as filter_file:
    #SAVES FILTER_NAME:FILTER_SET PAIR TO A DICTIONARY
    filter_list = {filter.split(" ")[0]:set(filter.split(" ")[1].split(",")) for filter in filter_file.read().strip().split("\n")} 

#DATA PROCESSING FUNCTIONS
def readGenome(combined_fasta_path): # READING COMBINED FASTA FILE INTO GENOME_DICT DICTIONARY, MAPPING HEADER (CONTAINING UNIQUE SAMPLE ID) TO THE GENOME SEQUENCE STRING
    genome_dict = {} # INITIALIZE A DICT TO STORE HEADER:SEQUENCE PAIRS
    with open(combined_fasta_path, 'r') as fasta_file: # OPEN FASTA FILE FOR READING
        for line in fasta_file: # GO THROUGH THE LINES OF THE FASTA FILE
            if line[0] == '>': # IF LINE IS A HEADER (NEW HEADER IS REACHED)
                fasta_header = line # STORE HEADER UNTILE NEXT HEADER IS REACHED
                genome_dict[fasta_header] = '' # INITIALIZE AN EMPTY SEQUENCE MAPPED TO HEADER
            if not line[0] == '>': # IF LINE IS A PART OF THE GENOME SEQUENCE
                genome_dict[fasta_header] += line.rstrip() # ADD THE LINE TO THE STRING THAT IS MAPPED TO THE LAST REACHED HEADER
    return genome_dict # DICTIONARY CONTAINING HEADER:GENOME PAIRS FROM COMBINED FASTA FILE


def get_nt_counting_stats(combined_fasta_path, skip=False, headers_changed = True):
    if skip: # IF SKIPPING OPTION IS CHOSEN
        print("Statistics calculation option skipped.")
        return None # RETURN NONE FOR TESTING PURPOSES
    elif not skip: # IF SKIPPING OPTION IS NOT CHOSEN
        sequence_stats_df = pd.DataFrame(columns=["receiving_lab_sample_id","genome_length","genome_N_percentage","genome_GC_content"]) # INITIALIZING DATAFRAME TO STORE THE RESULTING STATISTICS
        genome_dict_items = readGenome(combined_fasta_path).items() # GET HEADER:SEQUENCE PAIRS FROM READGENOME() OUTPUT
        for header, sequence in genome_dict_items: # LOOPING THROUGH EACH HEADER-SEQUENCE PAIR IN GENOME DICT
            sequence_length = len(sequence) # CALCULATING GENOME LENGTH
            gc_content = round(100 * (sequence.count('G') + sequence.count('C')) / len(sequence), 2) # CALCULATING GC-CONTENT FOR EACH SEQUENCE
            n_content = round(100 * (sequence.count('N') / len(sequence)), 2) # CALCULATING INVALID BASE CONTENT FOR EACH SEQUENCE
            if re.search(r'[A-Z]{2}[0-9]{4}\.B([0-9]{2}|[0-9])', header): new_header = re.search(r'[A-Z]{2}[0-9]{4}\.B([0-9]{2}|[0-9])', header).group(0)
            elif headers_changed: new_header = header.split("/")[2] # EXTRACTING SAMPLE ID FROM FASTA HEADER TO BE USED LATER IN THE COLUMN MAPPING PROCESS
            else: new_header = header.split("_")[1].strip() # EXTRACTING SAMPLE ID FROM FASTA HEADER TO BE USED LATER IN THE COLUMN MAPPING PROCESS IF FASTA HEADERS WERE NOT CHANGED
            sequence_stats_df = sequence_stats_df.append({"receiving_lab_sample_id":new_header, "genome_length":sequence_length, "genome_N_percentage":n_content, "genome_GC_content":gc_content}, ignore_index=True)
            #ADDING RESULT ROW TO THE DATAFRAME
        return sequence_stats_df # RETURN RESULTING DATAFRAME


def csv_info_extractor(file_path, filter_list, f_value):
    """The function is used to exctract mutation data from vcf-based csv reports based on set of provided filters"""
    try:
        sample_id = file_path.split('/')[-1].split("_", 1)[1]
        sample_id = sample_id[:len(sample_id) - 8]
        processing_id = file_path.split('/')[-1][:9]
        df = pd.read_csv(file_path)
        result_row = {"SAMPLE_ID":sample_id, "processing_id":processing_id}
    except:
        print(file_path)
        sys.exit(1)
    for key in filter_list.keys(): 
        if use_f_filter:
            mut_df = df.loc[(df["AMINO_ACID_CHANGE"].isin(filter_list[key])) & (df["FREQUENCY"] >= f_value)]
        elif use_p_filter:
            mut_df = df.loc[(df["AMINO_ACID_CHANGE"].isin(filter_list[key])) & (df["P_ERR_MUT_CALL"] < p_value)]
        if len(mut_df["P_ERR_MUT_CALL"]) == len(filter_list[key]): #IF THE NUMBER OF EXTRACTED (UNIQUE) ROWS MATCHES THE NUMBER OF MUTATION_NAMES FOR A GIVEN FILTER_NAME
            result_row[key] = str(1) #ASSUME THAT THE MUTATION SPECIFIED BY THE FILTER WAS FOUND - ADD FILTER_NAME:1 PAIR TO THE RESULT_ROW DICT
        else: result_row[key] = str(round(len(mut_df["P_ERR_MUT_CALL"])/len(filter_list[key]), 2)) #IF THE NUMBER OF EXTRACTED (UNIQUE) ROWS DOES NOT MATCH THE NUMBER OF MUTATION_NAMES FOR A GIVEN FILTER_NAME - ADD FILTER_NAME:%MATCH PAIR TO THE RESULT_ROW DICT
    return result_row


def depth_info_extractor(file_path):
    data = pd.read_csv(file_path, delimiter='\t').iloc[:,2]
    path_split = file_path.split('/')
    sample_id = path_split[-1].split("_", 1)[1]
    sample_id = sample_id[:len(sample_id) - 14]
    seq_date = path_split[-2].split("-")[1].replace("_","-")
    seq_lab = path_split[-2].split("-")[0].replace("_","(")+")"
    processing_id = file_path.split('/')[-1][:9]
    result_row = {'SAMPLE_ID':sample_id, "seq_date":seq_date, 'seq_institution':seq_lab, "processing_id":processing_id, 'AVERAGE_COVERAGE':sum(data)/len(data), 'MEDIAN_COVERAGE':data.median()}
    return result_row
    

def mapped_info_extractor(file_path):
    with open(file_path, "r+") as file: data = file.readlines()
    sample_id = file_path.split('/')[-1].split("_", 1)[1]
    sample_id = sample_id[:len(sample_id) - 18]
    processing_id = file_path.split('/')[-1][:9]
    try:
        result_row = {'SAMPLE_ID':sample_id, "processing_id":processing_id, 'TOTAL_READS':data[0].split(" ")[0], 'READS_MAPPED':data[4].split(" ")[0], 'MAPPED_FRACTION':round(int(data[4].split(" ")[0])/int(data[0].split(" ")[0]),2)}
    except IndexError:
        print(sample_id)
        print('Mapped_report_empty_error')
        sys.exit(1)
    return result_row


def mutation_report_generator(csv_list, coverage_list, output_path, f_value, all_excel_dump_df=None, sequence_stats_df=None, pango_report_path=None):
    """The function is used to generate mutation report based on data extracted from mutations reports for individual samples."""
    result_df = pd.DataFrame() #INIT EMPTY DATAFRAME FOR THE REPORT
    with concurrent.futures.ProcessPoolExecutor() as executor: # APPLYING THE EXTRACTOR FUNCTION IN-PARALLEL ON DIFFERENT CORES
        start_time = time.time() #TO ESTIMATE THE TIME IT TOOK TO PROCESS THE FILES
        results = [executor.submit(csv_info_extractor, file_path, filter_list, f_value) for file_path in csv_list] #SUBMITTING FUNCTION CALLS TO DIFFERENT PROCESSES
        for f in concurrent.futures.as_completed(results): #COLLECTING PROCESSING RESULTS TO THE DATAFRAME
            print(f'{f.result()["SAMPLE_ID"]} - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.now())}') #PRINTING THE TIME WHEN THE PROCESSING WAS FINISHED FOR A GIVEN SAMPLE
            result_df = result_df.append(f.result(), ignore_index=True) #ADD THE RESULT_ROW TO THE REPORT DF
        print(f"Processed {len(csv_list)} mutation reports in --- %s seconds ---" % round(time.time() - start_time, 2)) #THE TIME IT TOOK TO PROCESS ALL THE CSV FILES FROM FILE_PATH_SET    
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        start_time = time.time() #TO ESTIMATE THE TIME IT TOOK TO PROCESS THE FILES
        results = [executor.submit(depth_info_extractor, file_path) for file_path in coverage_list] #SUBMITTING FUNCTION CALLS TO DIFFERENT PROCESSES
        coverage_df = pd.DataFrame()
        for f in concurrent.futures.as_completed(results): #COLLECTING PROCESSING RESULTS TO THE DATAFRAME
            print(f'{f.result()["SAMPLE_ID"]} - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.now())}') #PRINTING THE TIME WHEN THE PROCESSING WAS FINISHED FOR A GIVEN SAMPLE
            coverage_df = coverage_df.append(f.result(), ignore_index=True) #ADD THE RESULT_ROW TO THE REPORT DF
        print(f"Processed {len(coverage_list)} coverage reports in --- %s seconds ---" % round(time.time() - start_time, 2))  
    result_df = pd.merge(result_df.applymap(str), coverage_df.applymap(str), how="left", on="SAMPLE_ID")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        start_time = time.time() #TO ESTIMATE THE TIME IT TOOK TO PROCESS THE FILES
        results = [executor.submit(mapped_info_extractor, file_path) for file_path in mapped_list] #SUBMITTING FUNCTION CALLS TO DIFFERENT PROCESSES
        mapped_df = pd.DataFrame()
        for f in concurrent.futures.as_completed(results): #COLLECTING PROCESSING RESULTS TO THE DATAFRAME
            print(f'{f.result()["SAMPLE_ID"]} - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.now())}') #PRINTING THE TIME WHEN THE PROCESSING WAS FINISHED FOR A GIVEN SAMPLE
            mapped_df = mapped_df.append(f.result(), ignore_index=True) #ADD THE RESULT_ROW TO THE REPORT DF
        print(f"Processed {len(mapped_list)} mapping reports in --- %s seconds ---" % round(time.time() - start_time, 2))
    result_df = pd.merge(result_df.applymap(str), mapped_df.applymap(str), how="left", on="SAMPLE_ID")
    
    column_reorder_list = ['SAMPLE_ID', 'processing_id'] + list(filter_list.keys()) + ['AVERAGE_COVERAGE', 'MAPPED_FRACTION', 'READS_MAPPED', 'TOTAL_READS', 'MEDIAN_COVERAGE', 'seq_date', 'seq_institution']
    result_df = result_df[column_reorder_list] #REORDERING COLUMNS BASED IN ORDER OF FILTERS IN REPORT_FILTERS.TXT FILE
    result_df.rename(columns={'SAMPLE_ID':"receiving_lab_sample_id"}, inplace=True) #FOR MERGING PURPOSES

    if all_excel_dump_df is not None: #IF METADATA MINING OPTION SELECTED - ADD METADATA TO THE MUTATION REPORT
        result_df = pd.merge(result_df.applymap(str), all_excel_dump_df.applymap(str), how="left", on="receiving_lab_sample_id")
    if sequence_stats_df is not None: #IF STATISTICS CALCULATION OPTION SELECTED - ADD SEQUENCE STATISTICS DATA TO THE MUTATION REPORT
        result_df = pd.merge(result_df.applymap(str), sequence_stats_df.applymap(str), how="left", on="receiving_lab_sample_id")
    if pango_report_path is not None: #IF PANGOLIN TYPING OPTION SELECTED - ADD PANGOLIN LINEAGE DATA TO THE MUTATION REPORT
        pango_report_df= pd.read_csv(pango_report_path) # READ PANGOLIN REPORT TO A DATAFRAME
        if re.match(r'[A-Z]{2}[0-9]{4}\.B', pango_report_df['taxon'][0]): pango_id_dict = {taxon:taxon for taxon in pango_report_df['taxon']}# MAP TAXOT TO SAMPLE ID IN A DICT
        elif rename_fasta_header: pango_id_dict = {taxon:taxon.split("/")[-2] for taxon in pango_report_df["taxon"]} # OPTION FOR GISAID-ACCEPTED HEADERS
        else: pango_id_dict = {taxon:taxon.split("_")[1] for taxon in pango_report_df["taxon"]} # MAP TAXON TO SAMPLE ID IN A DICT IF FASTA HEADERS WERE NOT CHANGED
        pango_report_df["receiving_lab_sample_id"] = pango_report_df["taxon"].map(pango_id_dict)
        pango_report_df = pango_report_df[["receiving_lab_sample_id", 'lineage']] #KEEP ONLY LINEAGE MATCHED AGAINST SAMPLE IDS
        result_df = pd.merge(result_df.applymap(str), pango_report_df.applymap(str), how="left", on="receiving_lab_sample_id")
    result_df.drop_duplicates(subset=['receiving_lab_sample_id'],inplace=True)

    ###EXPORT_CODE
    ###ADD COLUMN - FILL WITH CONSTANT
    #SEQUENCING_NOTES - NONE - TO BE FILLED MANUALLY
    result_df['sequencing_notes'] = np.zeros(len(result_df['receiving_lab_sample_id']))
    #RESULT_NOTES - NONE - TO BE FILLED MANUALLY
    result_df['result_notes'] = np.zeros(len(result_df['receiving_lab_sample_id']))
    #SUB_LINEAGE - NONE FOR COVID - MLST FOR BACT
    result_df['sub_lineage'] = np.zeros(len(result_df['receiving_lab_sample_id']))
    #ANALYSIS_PIPELINE - BWA 0.7.17-R1198-DIRTY MEM - PIPELINE-SPECIFIC
    pipeline_string = 'bwa 0.7.17-r1198-dirty mem'
    result_df['analysis_pipeline'] = np.chararray(result_df['receiving_lab_sample_id'].shape, itemsize=len(pipeline_string)+1).tostring()
    result_df.loc[result_df.processing_id != 'Z_BMC', 'analysis_pipeline'] = pipeline_string

    # ANALYSIS_INSTITUTION - DEFAULT NMRL
    an_inst = 'NMRL'
    result_df['analysis_institution'] = np.chararray(result_df['receiving_lab_sample_id'].shape, itemsize=len(an_inst)+1).tostring()
    result_df.loc[result_df.processing_id != 'Z_BMC', 'analysis_institution'] = an_inst
    result_df.loc[result_df.processing_id == 'Z_BMC', 'analysis_institution'] = 'BMC'

    #SEQUENCING_PLATFORM - ILLUMINA/MGI - DEDUCE FROM LIBRARY PREPARATION
    result_df['sequencing_platform'] = np.chararray(result_df['receiving_lab_sample_id'].shape, itemsize=len('illumina')+1).tostring()
    result_df.loc[result_df.processing_id != 'Z_BMC', 'sequencing_platform'] = 'illumina'
    result_df.loc[result_df.processing_id == 'Z_BMC', 'sequencing_platform'] = 'MGI'

    ###EXTRACT FROM REPORT FOLDER NAMES
    # ANALYSIS_DATE & ANALYSIS_BATCH_ID
    analysis_batch_id = report_path.split('/')[-1]
    analysis_date = analysis_batch_id.split('_')[1]
    result_df['analysis_date'] = np.chararray(result_df['receiving_lab_sample_id'].shape, itemsize=len(analysis_date)+1).tostring()
    result_df['analysis_batch_id'] = np.chararray(result_df['receiving_lab_sample_id'].shape, itemsize=len(analysis_batch_id)+1).tostring()
    result_df.loc[result_df.processing_id != 'Z_BMC', ['analysis_date','analysis_batch_id']] = [analysis_date,analysis_batch_id]

    ###EXTRACT FROM OUTPUT FOLDER NAMES
    #GET LIST OF SAMPLE
    sample_ids = result_df['receiving_lab_sample_id']
    sample_ids.to_csv(f'{report_path}/sample_ids.csv',header=False,index=False)

    #FIND PATHS TO FASTQ FILES USING BASH UTILITIES
    os.system(f'ls -R {raw_folder_path} | grep 1.fastq.gz > {report_path}/fastq_files.csv')
    os.system(f'cat {report_path}/sample_ids.csv | while IFS="," read a ; do cat {report_path}/fastq_files.csv | grep "$a" ; done > {report_path}/sample_file_list.csv')
    os.system(f'cat {report_path}/sample_file_list.csv | xargs -n1 -P12 -I% find {raw_folder_path}/ -type f -name % > {report_path}/path_list.csv')
    os.system(f'rm {report_path}/sample_file_list.csv {report_path}/fastq_files.csv {report_path}/sample_ids.csv')
    sample_paths = pd.read_csv(f'{report_path}/path_list.csv', header=None)
    sample_frame = pd.DataFrame({'id':[],'path':[]})

    #MAP RUN FOLDER TO SAMPLE ID
    for id in sample_ids:
        paths = sample_paths.loc[sample_paths[sample_paths.columns[0]].str.contains(id)].reset_index() #GET LIST OF PATHS THAT CONTAIN SAMPLE ID
        paths.rename(columns={sample_paths.columns[0]:'path'}, inplace=True) #PROVIDE CORRECT COLUMN NAME
        paths["id"] = [id for _ in range(len(paths.path))] #ADD ID
        paths.drop('index', inplace=True, axis=1) #REMOVE INDEX COLUMN
        sample_frame = sample_frame.append(paths) #COMBINE IN ONE DF
    os.system(f'rm {report_path}/path_list.csv')

    #KEEP RUN FOLDER FROM PATHS
    sample_frame[['raw','run_folder','fastq']] = sample_frame['path'].str.rsplit('/', n=2, expand=True) #GENERATE NEEDED COLUMNS
    sample_frame.drop('raw', inplace=True, axis=1) #REMOVE USELESS COLUMN
    sample_frame.drop('fastq', inplace=True, axis=1) #REMOVE USELESS COLUMN
    eurofins_samples = sample_frame[~sample_frame['run_folder'].str.contains('nmrl')]
    nmrl_samples = sample_frame[sample_frame['run_folder'].str.contains('nmrl')]

    #SPLIT TEXT TO COLUMNS AND DROP SAMPLE DUPLICATES
    if len(eurofins_samples) > 0:
        eurofins_samples[['seq_date','1','2','3','4']] = eurofins_samples['run_folder'].str.rsplit('-', n=4, expand=True)
    if len(nmrl_samples) > 0:
        nmrl_samples[['seq_date','lab','kit','instrument', 'seq_mode', 'primer']] = nmrl_samples['run_folder'].str.rsplit('-', n=5, expand=True)
    eurofins_samples.drop_duplicates(subset='id', keep='first', inplace=True)
    nmrl_samples.drop_duplicates(subset='id', keep='first', inplace=True)

    #ADD EXTRACTED INFO TO RESULT DF
    #INIT EMPTY COLUMNS TO STORE INFORMATION
    result_df["used_sequencing_run_ids"] = [None for _ in range(len(result_df))]
    result_df["used_batch_ids"] = [None for _ in range(len(result_df))]
    result_df["library_prep_method"] = [None for _ in range(len(result_df))]
    result_df["analysis_pipeline_notes"] = [None for _ in range(len(result_df))]

    #FOR EACH SAMPLE ID IN REPORT
    for id in result_df['receiving_lab_sample_id']:
        reported_samples = result_df.loc[result_df['receiving_lab_sample_id'] == id] #GET ALL ROWS THAT CONTAIN THIS ID
        if len(reported_samples) == 1:  #IF ID IS UNIQUE
            reported_samples = reported_samples.reset_index(drop=True) #TO GET 0 INDEX
            if reported_samples['seq_institution'][0] == 'NMRL(LIC)':  #IF SEQUENCING INSTITUTION IS DETERMINED AS NMRL
                    run_folder = nmrl_samples.loc[nmrl_samples.id == id]['run_folder'][0]
                    batch_id = nmrl_samples.loc[nmrl_samples.id == id]['seq_date'][0]
                    used_batch_ids = f'{batch_id}'
                    library_prep_method = nmrl_samples.loc[nmrl_samples.id == id]['kit'][0]
                    analysis_pipeline_notes = nmrl_samples.loc[nmrl_samples.id == id]['primer'][0]
            else:  #IF NOT NMRL
                run_folder = eurofins_samples.loc[eurofins_samples.id == id]['run_folder'][0]
                batch_id = eurofins_samples.loc[eurofins_samples.id == id]['seq_date'][0]
                library_prep_method = 'eurofins in-house'
                used_batch_ids = f'Eurofins_{batch_id}'
                analysis_pipeline_notes = 'arctic v4'
            #ADD RESULTS TO THE DATAFRAME (ID CAN BE USED AS UNIQUE ROW IDENTIFIER)
            result_df.loc[result_df.receiving_lab_sample_id == id, 'used_batch_ids'] = used_batch_ids
            result_df.loc[result_df.receiving_lab_sample_id == id, 'used_sequencing_run_ids'] = run_folder
            result_df.loc[result_df.receiving_lab_sample_id == id, 'library_prep_method'] = library_prep_method
            result_df.loc[result_df.receiving_lab_sample_id == id, 'analysis_pipeline_notes'] = analysis_pipeline_notes
        
        elif len(reported_samples) > 1: #IF ID IS NOT UNIQUE (E.G. MORE THAN ONE SAMPLE WITH THE SAME ID SEQUENCED BY DIFFERENT INSTITUTIONS) (TO BE TESTED)
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
                #ADD RESULTS TO THE DATAFRAME (BY ROW INDEX)
                result_df.at[result_df.index == i, 'used_batch_ids'] = used_batch_ids
                result_df.at[result_df.index == i, 'used_sequencing_run_ids'] = run_folder
                result_df.at[result_df.index == i, 'library_prep_method'] = library_prep_method
                result_df.at[result_df.index == i, 'analysis_pipeline_notes'] = analysis_pipeline_notes
    
    result_df.to_csv(f'{output_path}/{"{0:%Y-%m-%d-%H-%M}".format(datetime.now())}_ft_{f_value}_mutation_report.csv', header=True, index=False, encoding='utf-8-sig') #GENERATING THE REPORT


#GENERATING SEQUENCE STATISTICS WITH OPTION TO SKIP AND PROCESS >LVA000_2104047246 AND >HCOV-19/LATVIA/2103098474/2021 TYPE FASTA HEADERS
if not skip_sequence_stats: sequence_stats_df = get_nt_counting_stats(combined_fasta_path, skip_sequence_stats, rename_fasta_header) 
mutation_report_generator(csv_list, coverage_list, source_files_path, f_value, all_excel_dump_df, sequence_stats_df, pango_report_path)

#MOVE REPORTS TO THE PARENT DIRECTORY
os.system(f'mv {source_files_path}/{"{0:%Y-%m-%d-%H-%M}".format(datetime.now())}_ft_{f_value}_mutation_report.csv {report_path}')
os.system(f'mv {source_files_path}/{timestr_fasta}_lineage_report.csv {report_path}')
if args.multiqc:
    multiqc_source_path = source_files_path.replace('/mnt/home/jevgen01/nmrl/cov_analysis/', '')
    multiqc_report_path = report_path.replace('/mnt/home/jevgen01/nmrl/cov_analysis/', '')
    os.system(f'cd ~/nmrl/cov_analysis; module load singularity; singularity run ~/nmrl/image_files/fastq_screen.sif multiqc --interactive {multiqc_source_path} -o {multiqc_report_path}')


#BASIC LOGIC
#-READ FILTERS FROM REPORT_FILTERS.TXT
#-GET LIST OF CSV FILES
#-INIT AN EMPTY DF TO STORE TOTAL RESULTS - RESULT DF
#-ITERATE OVER THE CSV FILES
#-FOR EACH FILE
# - GET SAMPLE ID FROM FILE NAME
# - EXTRACT CONTENT TO DF
# - LOOK FOR MUTATIONS BASED ON FILTERS
# - COMBINE MUTATION INFORMATION IN TESSY REPORT TABLE ROW FORMAT
# - ADD ROW TO THE RESULT DF
#-WHEN ALL FILES ARE PROCESSED - GENERATE A CSV FROM RESULT DF



