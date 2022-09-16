import sys, pandas as pd, os, time, concurrent.futures, numpy as np

#READ SUMMARY DATA
summary_data_path = sys.argv[2] #PATH TO SUMMARY STATISTICS FILE
summary_data = pd.read_csv(summary_data_path) #READ METADATA FROM SUMMARY FILE
summary_data = summary_data[['receiving_lab_sample_id','sampling_date']] #SELECT RELEVANT COLUMNS
summary_data = summary_data[summary_data['sampling_date'].str.match(r'[0-9]{4}-[0-9]{2}-[0-9]{2}')==True] #FILTER INVALID VALUES BASED ON REGEX FOR DATE FORMAT
output_folder = '/mnt/home/groups/nmrl/cov_analysis/mutation_heatmap/mut_heatmap_data/'

def process_mutations(file_path):
    df = pd.read_csv(file_path) #READ EACH REPORT
    try:
        filtered_df = df[(df['FREQUENCY']>65) & (df['GENE'] == 'S') & (df['ANNOTATION'] != 'synonymous_variant')] #FILTER MUTATIONS BASED ON FREQUENCY
        allowed_mutations = list(set(filtered_df.AMINO_ACID_CHANGE)) #GET LIST OF UNIQUE MUTATIONS
    except KeyError:
        print(file_path)
        return 1
    # ADD UNIQUE MUTATIONS TO THE LIST OF ALL MUTATIONS
    return allowed_mutations

def count_mutation(file_path):
    file_name = file_path.split("/")[-1] #GET FILE NAME FROM FULL PATH
    if "COV" in file_name:
        id = file_name[:len(file_name)-8].split("_")[1] #EXTRACT REAL ID FROM FILE NAME
    else:
        id = file_name[:len(file_name)-8] #EXTRACT REAL ID FROM FILE NAME
    df = pd.read_csv(file_path) #READ REPORT INTO DF
    found_mut = df[df['AMINO_ACID_CHANGE'].isin(mutation_set)] #EXTRACT MUTATIONS FROM GIVEN REPORT
    found_mut = list(found_mut['AMINO_ACID_CHANGE']) #CONVERT MUTATIONS TO LIST 
    unfound_mut = mutation_data[~mutation_data['Mutations'].isin(found_mut)] #EXTRACT MUTATIONS THAT WERE NOT OBSERVED IN GIVEN SAMPLE
    total_mut = found_mut + list(unfound_mut['Mutations']) #GET ORDERED LIST OF FOUND AND UNFOUND MUTATIONS
    found_mut = list(map(lambda x:[1],found_mut)) #REPLACE MUTATION NAME WITH BINARY COUNT - 1 FOR FOUND MUTATIONS
    unfound_mut = list(map(lambda x:[0],list(unfound_mut['Mutations']))) #REPLACE MUTATION NAME WITH BINARY COUNT - 0 FOR UNFOUND MUTATIONS
    total_data = found_mut + unfound_mut #GET ORDERED LIST OF BINARY DATA
    data = pd.DataFrame.from_dict(dict(zip(total_mut,total_data)), orient='columns') #CONVERT TO DATAFRAME
    data.index = [id] #SET INDEX TO SAMPLE ID
    return data


#PARSING MUTATION REPORTS
mutation_set = [] #EMPTY LIST TO STORE MUTATIONS
folder_path = sys.argv[1] #PATH TO FOLDER CONTAINING MUTATION REPORTS
file_list = os.listdir(folder_path) #LIST OF ALL MUTATIONS REPORTS
print(f'Total files to process: {len(file_list)}')
total_time = time.time()
with concurrent.futures.ProcessPoolExecutor() as executor: #APPLY FUNCTION USING MULTIPROCESSING
        results = [executor.submit(process_mutations, f'{folder_path}/{file}') for file in file_list] #SUBMITTING FUNCTION CALL TO DIFFERENT PROCESSES
        for f in concurrent.futures.as_completed(results): #COLLECTING RESULTS TO A LIST
            result = f.result()
            if result != 1: 
                mutation_set+=f.result()
            else: #IF FILE FORMAT IS INCORRECT (SOME COLUMNS MISSING)
                continue
mutation_set=[mut for mut in set(mutation_set) if str(mut) != 'nan'] #FILTER NAN VALUES AND KEEP ONLY UNIQUE MUTATIONS

batch_data = pd.DataFrame() #TO STORE RAW MUTATION DATA IN BINARY FORMAT FOR EACH SAMPLE
mutation_data = pd.DataFrame({'Mutations':mutation_set}) #CONVERT LIST TO DF TO USE PANDAS MATCHING TOOLS
with concurrent.futures.ProcessPoolExecutor() as executor: #APPLY FUNCTION USING MULTIPROCESSING
        results = [executor.submit(count_mutation, f'{folder_path}/{file}') for file in file_list] #SUBMITTING FUNCTION CALL TO DIFFERENT PROCESSES
        for f in concurrent.futures.as_completed(results): #COLLECTING RESULTS TO A DATAFRAME
            batch_data = pd.concat([batch_data, f.result()]) 
batch_data['receiving_lab_sample_id'] = batch_data.index #SET SAMPLE ID AS INDEX FOR MERGING
batch_data['receiving_lab_sample_id'] = batch_data['receiving_lab_sample_id'].apply(str) #CONVERT SAMPLE ID TO STR FOR MERGING
summary_data['receiving_lab_sample_id'] = summary_data['receiving_lab_sample_id'].apply(str) #CONVERT SAMPLE ID TO STR FOR MERGING
batch_data = pd.merge(batch_data, summary_data, how="left", on="receiving_lab_sample_id") #ADDING SAMPLING DATE COLUMN
batch_data.drop_duplicates(subset=['receiving_lab_sample_id'],keep='first', inplace=True) #REMOVE DUPLICATES
batch_data.to_csv(f'{output_folder}mutation_binary_table.csv', header=True, index=False) #SAVE DATA AS CSV
del batch_data['receiving_lab_sample_id'] #REMOVE UNUSED COLUMN
n_by_date = batch_data.groupby("sampling_date").sum()/batch_data.groupby("sampling_date").count() #COMPUTE FREQUENCY OF ALL MUTATIONS
c_by_date = batch_data.groupby("sampling_date").count() #COUNT DATA

n_by_date = n_by_date.T #TRANSPOSE DF
c_by_date = c_by_date.T #COUNTS TRANSPOSED
n_by_date = n_by_date.stack().reset_index() #STACK COLUMNS TO USE IN HEATMAP
c_by_date = c_by_date.stack().reset_index() #STACK COLUMNS TO USE IN HEATMAP TOOLTIP
n_by_date.columns = ['Mutation','sampling_date','Fill'] #RENAME COLUMNS
c_by_date.columns = ['Mutation','sampling_date','Count'] #RENAME COLUMNS
n_by_date=n_by_date[n_by_date['Fill']>=0] #FILTER MUTATIONS WITH FREQUENCY LOWER THAN SET
n_by_date['Count'] = c_by_date["Count"]
n_by_date.to_csv(f'{output_folder}mutation_frequency_table.csv', header=True, index=False) #SAVE DATA AS CSV
print(f'Finished mutation frequency computations in {round(time.time() - total_time, 2)} seconds.')

print('Starting filtering by 2 consecutive dates.')
consec_filter = 1 #FILTER SIZE - 1
unique_dates = sorted(n_by_date['sampling_date'].unique())

#GENERATING FILTER RANGE LISTS
filter_list = []
for i in range(0,len(unique_dates),consec_filter): #"SLIDING WINDOW" APPROACH - ITERATE OVER SORTED (ASC) LIST OF UNIQUE DATES
    if i+consec_filter < len(unique_dates): #IF LIST CONTAINS DATES TO INCLUDE INTO FILTER (TO AVOID OUT-OF-BOUNDS)
        filter_list.append([unique_dates[i+l] for l in range(consec_filter+1)]) #ADD CURRENT DATE AND NEXT n DATES (n = number of consecutive dates to observe mutations)
filter_array = np.array(filter_list, np.datetime64) #TO NUMPY ARRAY TO SPEED-UP FILTERING
ordered_blocks = [n_by_date[(n_by_date['sampling_date'].astype(np.datetime64).isin(filter_array[i]))] for i in range(len(filter_array))] #SPLIT RAW DATA INTO BLOCKS ORDERED BY CONSECUTIVE DATES

output_df = pd.DataFrame() #TO STORE RESULTS
for block in ordered_blocks: 
    mutation_filter = block[block['Fill'] == 0]['Mutation'].unique() #BLOCKS ARE DEFINED BY CONSECUTIVE DATES - IF ANY MUTATION HAS FREQUENCY 0, IT HAS TO BE EXCLUDED FROM THE GIVEN BLOCK
    filtered_blocks = block[~block['Mutation'].isin(mutation_filter)] #EXCLUDING MUTATIONS THAT WERE NOT OBSERVED IN A GIVEN BLOCK FOR n CONSECUTIVE DATES*
    #IF MUTATION OCCURS FOR n CONSECUTIVE DAYS IN ANOTHER BLOCK, IT WILL BE INCLUDED INTO THE RESULTING DATAFRAME
    output_df = output_df.append(filtered_blocks) #STORING FILTERING RESULTS

#DUPLICATES ARISE FROM OVERLAP, WHICH IS NEEDED TO CHECK ALL CONSECUTIVE DATE
output_df.drop_duplicates(keep='last').to_csv(f'{output_folder}filtered_mutation_table.csv', header=True, index=False)

#CLEAUP
print('Adding permissions & performing cleanup')
os.system(f'mv {output_folder}mutāciju_apkopojums* {output_folder}backup/')

#GENERATING MUTATION PLOT
os.system(f'Rscript generate_mutation_heatmap.R {output_folder}filtered_mutation_table.csv')
os.system(f'chmod -R 775 {output_folder}*')