#!/mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/tools/rbase_env/bin/python
import pandas as pd, sys, shutil, re
import argparse
import os,  time

db_file = [f'/mnt/home/groups/nmrl/cov_analysis/analysis_history/{file}' for file in os.listdir('/mnt/home/groups/nmrl/cov_analysis/analysis_history/') if 'summary_file' in file][0] #LOOKUP FOR THE DATABASE FILE IN THE FOLDER WHERE SCRIPT IS LOCATED
log_path = '/mnt/home/groups/nmrl/cov_analysis/analysis_history/database_log_files'
backup_path = '/mnt/home/groups/nmrl/cov_analysis/analysis_history/backup'
print(f'Current database file: {db_file}')
cur_time = time.strftime("%d_%m_%Y") #TO TIMESTAMP UPDATES OF SUMMARY FILES IN FILE NAME

#CMD ARGUMENTS & SCRIPT USAGE MESSAGES
parser = argparse.ArgumentParser(description='A script to update database file with new batch data.') 
parser.add_argument('-d', '--dupl', metavar='\b', help = 'Full path to the csv files containing samples with ids that are already in db file', default=None, required=False)
parser.add_argument('-p', '--pipe', metavar='\b', help = 'Full path to the mutation_report.csv', default=None, required=False)
parser.add_argument('-e', '--export', help = 'Export contents of summary file in the SISdb format', action="store_true")



print('RUNNING WITH -p FLAG')
print('INFO: The -b flag is used to integrate new batch report data into db file.')
print('INFO: The script will automatically check if db file and batch report file contain matching sample ids.')
print('INFO: The user will be prompted to confirm db file changing and backup file generation steps in the command line.')
print('INFO: Duplicates.csv file is generated by this script if matching sample ids were found in new file and db file.')
print('INFO: added_unique_records.csv file is generated by this script following corresponding user selection.\n')

print('RUNNING WITH -d FLAG:')
print('INFO: The -d flag is used to integrate duplicated values into db file.')
print('INFO: Before running the script on the duplicates.csv file, the file should be manually edited by filling "processed" column as specified below.')
print('INFO: values 0, 1, 2 should be entered into "processed" column of duplicates.csv to control the processing of duplicated data.')
print('INFO: 0 - the record from db file should be replaced with record from batch report file.')
print('INFO: 1 - the record from batch report file should be added to db file as new records (make sure to change the id, otherwise duplicated records will be created).')
print('INFO: 2 - the record batch report file should be dropped.')
print('INFO: removed_duplicates.csv will be generated if value 2 was entered into "processed" column for any sample in duplicates.csv')


#IF SCRIPT IS RUN WITHOUT ARGUMENTS - PRINT HELP MESSAGE & EXIT
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()


###EXPORT OPTION
if args.export:
    #READ SUMMARY FILE
    summary_file = pd.read_csv(db_file)
    #CONVERT DATAFRAME INTO SISdb format
    rename_dict = {
        "receiving_lab_sample_id":"alt_sample_id",
        "processing_id":"analysis_id",
        "seq_institution":"sequencing_institution",
        "lineage":"lineage",
        "genome_length":"assembly_length",
        "genome_N_percentage":"coverage_percentage",
        "genome_GC_content":"GC_content",
        "AVERAGE_COVERAGE":"average_coverage",
        "MAPPED_FRACTION":"mapped_fraction",
        "READS_MAPPED":"reads_mapped",
        "TOTAL_READS":"total_reads",
        "seq_date":"sequencing_date",
        "testing_lab_sample_id":"sample_id",
        "MEDIAN_COVERAGE":"median_coverage",
        "used_batch_ids":"used_batch_ids",
        "used_sequencing_run_ids":"used_sequencing_run_ids",
        "sub_lineage":"sub_lineage",
        "result_notes":"result_notes",
        "analysis_pipeline_notes":"analysis_pipeline_notes",
        "analysis_institution":"analysis_institution",
        "analysis_date":"analysis_date",
        "analysis_batch_id":"analysis_batch_id",
        "analysis_pipeline":"analysis_pipeline",
        "sequencing_platform":"sequencing_platform",
        "library_prep_method":"library_prep_method",
        "sequencing_notes":"sequencing_notes",
        "testing_lab":"sample_origin_lab",
        "age":"age",
        "gender":"sex",
        "sample_type":"sample_type",
        "sampling_date":"sample_collection_date"
    }

    column_reorder = [
    "alt_sample_id",
    "analysis_id",
    "sequencing_institution",
    "lineage","assembly_length",
    "coverage_percentage",
    "GC_content",
    "sample_collection_date",
    "sample_origin_lab",
    "age","sex",
    "sample_type",
    "average_coverage",
    "mapped_fraction",
    "reads_mapped",
    "total_reads",
    "sequencing_date",
    "sample_id",
    "median_coverage",
    "used_batch_ids",
    "used_sequencing_run_ids",
    "sub_lineage",
    "result_notes",
    "analysis_pipeline_notes",
    "analysis_institution",
    "analysis_date",
    "analysis_batch_id",
    "analysis_pipeline",
    "sequencing_platform",
    "library_prep_method",
    "sequencing_notes"
    ]

    summary_file.rename(columns=rename_dict, inplace=True)
    #COMPUTE COVERAGE FROM N%
    summary_file['coverage_percentage'] = round((1-summary_file['coverage_percentage']/100),2)
    #REPLACE COVERAGE VALUES WITH 0 WHERE ASSEMBLY LENGTH IS 0
    summary_file.loc[summary_file.assembly_length == 0, 'coverage_percentage'] = 0.0
    #EXPORT DATAFRAME WITH NAME SISdb_yyyy_mm_dd_hh_mm_ss.csv (backup export to db_log_files)
    sisdb_export = summary_file[column_reorder]
    file_name = f'SISdb_{time.strftime("%Y_%m_%d_%H:%M:%S")}.csv'
    sisdb_export.to_csv(file_name, header=True, index=False)
    os.system(f'cp {file_name} database_log_files/')
    



if args.pipe: #IF SCRIPT IS USED TO ADD DATA FOR A NEW BATCH TO DATABASE FILE
    pipe_data = pd.read_csv(args.pipe) #READ PIPELINE REPORT
    pipe_data['testing_lab_sample_id'] = pipe_data['receiving_lab_sample_id'] #TO MATCH COLUMN USED TO STORE DOUBLE-LABELLED IDS (ARTIFACT)
    pipe_data['sample_type'] = [0 for _ in range(len(pipe_data))] #TO MATCH COLUMN STORING SAMPLE TYPE (ARTIFACT)  

    db_frame = pd.read_csv(db_file).applymap(str) #NORMALIZE THE DB FILE DATAFRAME TO STRING TYPE
    batch_data = pipe_data[db_frame.columns] #EXTRACTING & REORDERING RELEVANT COLUMNS
    batch_data.fillna(0, inplace=True) #STANDARDIZE EMPTY RECORDS
    batch_data.drop_duplicates(subset="receiving_lab_sample_id", keep='first', inplace=True) #REMOVE RECORDS WHERE SAMPLE ID IS DUPLICATED
    bmc_records = db_frame[db_frame.processing_id == "Z_BMC"] #SAVING TO REINSERT AT WRITING STAGE
    db_frame = db_frame[db_frame.processing_id != "Z_BMC"] #EXCLUDING BMC DATA FROM SEARCH TO AVOID DUPLICATES
    db_frame['receiving_lab_sample_id_1'] = db_frame.receiving_lab_sample_id #REPEAT SAMPLE ID COLUMN TO USE IN INDEXING IN DB FILE DATAFRAME
    
    reind_df_db = db_frame.set_index('receiving_lab_sample_id_1') #CREATE A COPY OF DB FILE DATAFRAME USING SAMPLE IDS AS INDEXES
    batch_data['receiving_lab_sample_id_1'] = batch_data.receiving_lab_sample_id #REPEAT SAMPLE ID COLUMN TO USE IN INDEXING IN BATCH DATAFRAME
    reind_df_batch = batch_data.set_index('receiving_lab_sample_id_1') #CREATE A COPY OF BATCH DATAFRAME USING SAMPLE IDS AS INDEXES IN COMPARISON BETWEEN DB FILE AND BATCH FILE

    batch_data.age = batch_data.age.astype(str).astype(float).astype(int).astype(str) #NORMALIZE AGE COLUMN IN BATCH DATAFRAME (NEEDED TO AVOID NUMERIC TYPE MISMATCH FOR IDENTICAL VALUES DURING COMPARISON (E.G. 10(INT) == 10.0(FLOAT) => False))
    db_frame.age = db_frame.age.astype(float).astype(int).astype(str) #NORMALIZE AGE COLUMN IN DB FILE DATAFRAME
    
    
    match = reind_df_batch.isin(reind_df_db) #IF ALL ROW-COLUMN COMBINATIONS FROM BATCH DATAFRAME WITH ROW-COLUMN COMBINATIONS FROM DB FILE DATAFRAME
    
    match.index = match.index.astype(str) #NORMALIZE INDEXES AS TYPE STRING
    batch_data.index = batch_data.index.astype(int) #NORMALIZE INDEXES OF BATCH FRAME AS TYPE INT TO LOOK FOR DUPLICATES
    db_frame.index = db_frame.index.astype(int) #NORMALIZE INDEXES OF DB FILE FRAME AS TYPE INT TO LOOK FOR DUPLICATES
    duplicates = False #INITIALIZING TO DEFAULT STATE - NO DUPLICATES EXPECTED BY-DEFAULT BETWEEN DB FILE AND BATCH FILE
    
    if any(match.receiving_lab_sample_id): #IF THERE IS A MATCH IN SAMPLE ID BETWEEN DB FILE AND BATCH FILE
        duplicates = True #SET TO TRUE - PROCESS DUPLICATES
        samples = match.index #GET ALL SAMPLE IDS IN THE MATCH MATRIX
        condition_dup_smpl = match['receiving_lab_sample_id'] == True #CONDITION TO FIND SAMPLE IDS THAT ARE FOUND BOTH IN DB FILE AND BATCH FILE
        
        condition_unq_smpl = match['receiving_lab_sample_id'] == False #CONDITION TO FIND SAMPLE IDS THAT ARE FOUND ONLY IN BATCH FILE
        unique_samples = list(map(str, list(samples[condition_unq_smpl]))) #EXTRACTING SAMPLE IDS OF SAMPLES UNIQUE TO BATCH FILE
        duplicate_samples = list(map(str, list(samples[condition_dup_smpl]))) #EXTRACTING IDS OF SAMPLES THAT ARE FOUND BOTH IN DB FILE AND BATCH FILE
        
        indexes = batch_data.index #GET INDEXES OF BATCH FILE
        condition_unq_idx = batch_data['receiving_lab_sample_id'].isin(unique_samples) #CONDITION TO EXTRACT BATCH FILE INDEXES OF SAMPLES UNIQUE TO BATCH FILE
        condition_dup_idx = batch_data['receiving_lab_sample_id'].isin(duplicate_samples) #CONDITION TO EXTRACT BATCH FILE INDEXES OF SAMPLES THAT ARE FOUND BOTH IN DB FILE AND BATCH FILE
        
        unique_indexes = list(indexes[condition_unq_idx]) #EXTRACTING BATCH FILE INDEXES OF SAMPLES UNIQUE TO BATCH FILE
        duplicate_indexes = list(indexes[condition_dup_idx]) #EXTRACTING BATCH FILE INDEXES OF SAMPLES THAT ARE FOUND BOTH IN DB FILE AND BATCH FILE
        duplicate_frame = batch_data.iloc[duplicate_indexes] #GENERATING FRAME WITH DATA FOR SAMPLES THAT ARE FOUND BOTH IN DB FILE AND BATCH FILE
        
        unique_records = batch_data.drop(duplicate_indexes) #KEEPING ONLY DATA ON SAMPLES UNIQUE TO BATCH FILE IN THE BATCH FILE FRAME
        unique_records.drop(['receiving_lab_sample_id_1'], axis=1, inplace=True) #REMOVING COLUMN USED FOR INDEXING
        match = match.drop(unique_samples) #REMOVING DUPLICATE ENTRIES FROM MATCH MATRIX AND SAVING IT TO BE USED IN DUPLICATE PROCESSING
    
    else: #IF THERE IS NO MATCH IN SAMPLE ID BETWEEN DB FILE AND BATCH FILE 
        unique_records = batch_data #CONSIDER ALL ENTRIES IN BATCH FILE UNIQUE

    #INPUT OPTIONS TO APPEND UNIQUE ENTRIES FROM BATCH FILE TO DB FILE
    update_unique_to_df = 'y' #input('Would you like to add unique entries into DB file? (Y\y)/(N/n): ')
    if update_unique_to_df in ['y','Y']:
        print('Unique records were appended to DB file.')
        record_unique = 'y'#input('Would you like to save added entries in separate csv file? (Y\y)/(N/n): ')
        
        if record_unique in ['y','Y']:
            print('Saving added unique records into separate csv file.')
            shutil.move(f'{db_file}', f'{backup_path}/{os.path.basename(db_file)}') #BACKING-UP DATABASE FILE TO REVERT UPDATE IF NEEDED
            unique_records.to_csv(f'{log_path}/added_unique_records_{cur_time}.csv', header=True, index=False)
        
        elif record_unique in ['n','N']:
            print('Added records:\n')
            print(unique_records)
        
        else: 
            print('Invalid user input on unique records save.')
            sys.exit(1)
        
        db_frame = db_frame.append(unique_records) #APPENDING UNIQUE RECORDS TO DB FILE FRAME
        db_frame.drop(['receiving_lab_sample_id_1'], axis=1, inplace=True) #REMOVING COLUMN USED IN INDEX SEARCH
        db_file_new = db_file.replace(re.search(r'[0-9]{2}_[0-9]{2}_[0-9]{4}',db_file).group(0),cur_time)
        db_frame = db_frame.sort_values(by=['processing_id'], ascending=True)
        db_frame = bmc_records.append(db_frame) #ADDING BMC RECORDS BACK TO THE DF
        db_frame.to_csv(db_file_new,header=True, index=False) #UPDATING DB FILE

    elif update_unique_to_df in ['n','N']:
        print(unique_records)
        print('Unique records were NOT added to DB file.')
    
    else: 
        print('Invalid user input on db file update.')
        sys.exit(1)

    if duplicates: #IF THERE IS A MATCH IN SAMPLE ID BETWEEN DB FILE AND BATCH FILE
        duplicate_frame.drop(['receiving_lab_sample_id_1'], axis=1, inplace=True) #REMOVE REDUNDANT INDEXING COLUMN FROM DUPLICATE FRAME
        duplicate_frame.reset_index(drop=True, inplace=True) #RESET INDEXES FOR DUPLICATE FRAME
        match.reset_index(drop=True, inplace=True) #RESET INDEXES FOR MATCH MATRIX (REQUIRED TO MATCH INDEXES OF DUPLICATE FRAME AND MATCH MATRIX AFTER UNIQUE ENTRIES WERE REMOVED FROM MATCH MATRIX)
        
        dup_indexes = [idx for idx in match.index if all(match.iloc[idx])] #GENERATE LIST OF INDEXES FOR ROWS THAT FULLY MATCH ROWS FROM DB FILE (IDENTICAL VALUES IN ALL COLUMNS)
        duplicate_frame = duplicate_frame.drop(dup_indexes) #REMOVING FULL MATCHES BY INDEX
        sample_ids = list(duplicate_frame.receiving_lab_sample_id) #GENERATING LIST OF SAMPLE IDS THAT REMAIN AFTER FULL MATHCES WERE EXCLUDED
        db_duplicates = db_frame[db_frame['receiving_lab_sample_id'].isin(sample_ids)] #GENERATING NEW FRAME FROM DB FILE CONTAINING ONLY SAMPLES WITH NOT-FULLY MATCHED ROWS (VALUES IS SOME COLUMNS DIFFER, BUT IDS ARE THE SAME)
        db_duplicates = db_duplicates.add_suffix('_DB') #ADDING SUFFIX TO COLUMNS EXTRACTED FROM DB FILE
        db_duplicates = db_duplicates.rename(columns={'receiving_lab_sample_id_DB':'receiving_lab_sample_id'}) #NORMALISING COLUMN NAME TO MERGE UPON
        merge_product = duplicate_frame.merge(db_duplicates, on='receiving_lab_sample_id', how='outer') #COMBINE VALUES EXTRACTED FROM DB FILE WITH VALUES EXTRACTED FROM BATCH FILE (KEEP BOTH VALUES)
        
        print('\nThe following duplicate rows are found to contain mismatches for the same id in the db file.')
        print(merge_product)

        if not merge_product.empty: #IF THERE ARE ANY PARTIAL MATCHES BETWEEN DB FILE AND BATCH FILE FOR THE SAME ID - PROCESS DUPLICATE FRAME
            merge_product = merge_product[sorted(merge_product.columns, reverse=True)] #SORT COLUMNS SO THAT _DB-SUFFIXED COLUMNS ARE NEAR COLUMNS FROM BATCH FILE
            merge_product['process'] = [None for _ in merge_product.index] #ADD PROCESS COLUMN TO BE USED IN FURTHER PROCESSING OF DUPLICATE VALUES
            duplicates_to_csv = 'y' #input('Would you like to save list of duplicates as csv? (Y\y)/(N/n): ')
            
            if duplicates_to_csv in ['y','Y']:
                print('Saving duplicates to csv.')
                merge_product.to_csv(f'{log_path}/duplicates_{cur_time}.csv', header=True, index=False) #SAVE DUPLICATE FRAME TO CSV IF PROMPTED BY THE USER
                sys.exit(0)
            elif duplicates_to_csv in ['n','N']:
                print('Duplicate record export skipped.')
                sys.exit(0)
            else: 
                print('Invalid user input on duplicate records save.')
                sys.exit(1)
            

if args.dupl: #IF DUPLICATE PROCESSING FLAG WAS SET
    db_frame = pd.read_csv(db_file).applymap(str) #READ DB FILE FILE AND NORMALIZE AS TYPE STRING
    input_dupl = pd.read_csv(args.dupl) #READ DUPLICATES FILE AND NORMALIZE AS TYPE STRING

    if 'process' not in input_dupl.columns: #IF THERE IS NO PROCESS COLUMN IN DUPLICATES FILE, ABORT - COLUMN IS REQUIRED TO PROCESS DUPLICATE VALUES
        print('No process column in provided file.')
        sys.exit(1)
    
    columns = [col for col in input_dupl.columns if 'DB' not in col] #GET LIST OF COLUMNS FROM BATCH FILE
    input_dupl=input_dupl[columns] #KEEP ONLY COLUMNS FROM BATCH FILE IN DUPLICATES FRAME

    if 2 in input_dupl['process']: #VALUE IN 'PROCESS' COLUMN == 2 - DATA IS TO BE DISCARDED
        print('Exporting list of dropped duplicates.')
        removed_duplicates = input_dupl[input_dupl['process'] == 2]
        if len(removed_duplicates) > 0:
            removed_duplicates.to_csv(f'{log_path}/removed_duplicates_{cur_time}.csv',header=True, index=False) #GENERATE A CSV CONTAINING DISCARDED DATA FOR BACKUP PURPOSES
        input_dupl=input_dupl[input_dupl['process'] != 2] #DROPPING ROWS BASED ON 'PROCESS' COLUMN VALUE
    
    replace = input_dupl[input_dupl['process'] == 0] #VALUE IN 'PROCESS' COLUMN == 0 - DATA IN DB FILE IS TO BE REPLACED WITH DATA FROM DUPLICATES FILE
    append = input_dupl[input_dupl['process'] == 1] #VALUE IN 'PROCESS' COLUMN == 1 - DATA FROM DUPLICATES FILE IS TO BE ADDED AS NEW DATA TO THE DB FILE (SAMPLE ID CHANGING OPTION)
    replace = replace[[col for col in replace.columns if col != 'process']] #REMOVING 'PROCESS' COLUMN FROM REPLACE FRAME
    
    append = append[[col for col in append.columns if col != 'process']] #REMOVING 'PROCESS' COLUMN FROM APPEND FRAME
    print('The following data will replace data from DB.')
    print(replace)
    
    replace.receiving_lab_sample_id = replace.receiving_lab_sample_id.astype(str) #NORMALISING SAMPLE IDS AS TYPE STRING IN REPLACE FRAME
    sample_ids = list(replace['receiving_lab_sample_id']) #GETTING SAMPLE IDS FROM REPLACE FRAME
    index = db_frame.index #GETTING INDEXES FROM DB FILE FRAME
    
    condition = db_frame["receiving_lab_sample_id"].isin(sample_ids) #CONDITION TO EXTRACT INDEXES FROM DB FILE FRAME
    rep_indexes = index[condition].to_list() #GETTING LIST OF INDEXES WHERE REPLACEMENT WAS NEEDED
    db_frame.drop(rep_indexes, inplace=True) #REMOVING ROWS THAT REQUIRED REPLACEMENT BY INDEXES
    replace.sort_values(by=['processing_id'], inplace=True) #ORDER RECORDS BEFORE ADDING TO DATABASE
    db_frame = db_frame.append(replace) #APPENDING ROWS FROM REPLACE FRAME TO DB FILE FRAME
    print('The following data will be appended to the database.')
    print(append)
    
    db_file_new = db_file.replace(re.search(r'[0-9]{2}_[0-9]{2}_[0-9]{4}',db_file).group(0),cur_time)
    append.sort_values(by=['processing_id'], inplace=True) #ORDER RECORDS BEFORE ADDING TO DATABASE
    db_frame = db_frame.append(append) #APPENDING VALUES FROM APPEND FRAME TO DB FILE FRAME
    db_frame.to_csv(db_file_new,header=True, index=False) #UPDATING DB FILE
    os.system(f'mv {args.dupl} {log_path}/{args.dupl}')

###CLEANUP
os.system('chmod -R 775 ./*')
