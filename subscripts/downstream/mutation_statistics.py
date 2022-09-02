import pandas as pd, sys, os, datetime, argparse, datetime, re


parser = argparse.ArgumentParser(description='A script to calculate mutation statistics given range of dates or list of sample ids.') #argparser object to provide command-line functionality
parser.add_argument('-d1', '--start_date', metavar='\b', help = 'A starting sequencing date of the reference interval (YYYY-MM-DD).', default=None, required=False)
parser.add_argument('-d2', '--end_date', metavar='\b', help = 'An ending sequencing date of the reference interval (YYYY-MM-DD).', default=None, required=False)
parser.add_argument('-l', '--name_list', metavar='\b', help = 'Path to list of file names to lookup in folders', default=None, required=False)
parser.add_argument('-o', '--out_dir', metavar='\b', help = 'Path to output folder', default="./", required=False)

if len(sys.argv)==1: #if no command-line arguments provided - display help and stop script excecution
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args() #args list from command-line input


if (args.start_date == None or args.end_date == None) and args.name_list == None: #IF NO DATE AND NO SAMPLE ID LIST PROVIDED
    print("If using dates, both d1 and d2 should be provided. If one or both date is missing, a csv file with sample ids should be provided.")
    parser.print_help(sys.stderr)
    sys.exit(1)

date_1 = args.start_date
date_2 = args.end_date
id_list_path = args.name_list
out_dir = args.out_dir if args.out_dir[-1] == '/' else f"{args.out_dir}/"

if (args.start_date != None and args.end_date != None) and args.name_list != None: #IF SAMPLE ID LIST AND BOTH SEQUENCING DATES ARE PROVIDED
    print("WARNING: Both sample id list path and sequencing date range were provided. Using sample id list to calculate mutation statistics.")
    date_1 = None
    date_2 = None

def check_date(date_string):
    '''The function is used to verify that provided dates are in correct format.'''
    check = True if re.match("[0-9]{4}-[0-9]{2}-[0-9]{2}", date_string) else False
    return check    

file_list = []

if date_1 != None: #IF SELECTING BY DATE RANGE
    if check_date(date_1) and check_date(date_2):
        start = datetime.datetime.strptime(date_1, "%Y-%m-%d")
        end = datetime.datetime.strptime(date_2, "%Y-%m-%d") 
        for (root,dirs,files) in os.walk("/home/groups/nmrl/cov_analysis/covid_output/", topdown=True): #FINDING PATH TO EACH FILE THAT MATCHES TIME CRITERIA
            for dir in dirs: 
                if "NMRL" not in dir:
                     continue                    
                else: 
                    date = datetime.datetime.strptime(dir.split("-")[1],"%Y_%m_%d")
                    if start <= date <= end:
                        for file in os.listdir(os.path.join(root,dir)):
                            if ".ann.csv" in file:
                                file_list.append(os.path.join(root,dir,file))
else: #IF SELECTING BY SAMPLE LIST
    id_list = list(pd.read_csv(id_list_path, header=0).iloc[:,0]) #CONVERT COLUMN OF IDS TO LIST
    for (root,dirs,files) in os.walk("/home/groups/nmrl/cov_analysis/covid_output/", topdown=True): #FINDING PATH TO EACH FILE THAT MATCHES TIME CRITERIA
            for file in files:
                if '.ann.csv' in file:
                    if 'COV' in file: #IF INHOUSE ID IN FILE NAME
                       if file.split(".",1)[0][10:] in id_list:
                            file_list.append(os.path.join(root, file))  
                    else: #IF INHOUSE ID NOT IN FILE NAME
                        if file.split(".",1)[0] in id_list:
                            file_list.append(os.path.join(root, file))  
                

def calculate_statistics(file_list):
    valid_mut_dict = dict()
    valid_ann_dict = dict()
    sample_count = 0
    for file in file_list:
        if 'ann.csv' in file or '_variants.tsv' in file:
            sample_count +=1
            if len(valid_mut_dict) == 0:
                if 'ann.csv' in file:
                    df = pd.read_csv(file).applymap(str)
                    cur_mut_set = list(df.AMINO_ACID_CHANGE)
                    cur_nt_list = list(df.MUTATION)
                    valid_nt_dict = dict(zip(cur_mut_set, cur_nt_list))
                    cur_gene_list = list(df.GENE)
                    valid_gene_dict = dict(zip(cur_mut_set, cur_gene_list))
                    cur_ann_list = list(df.ANNOTATION)
                    valid_ann_dict = dict(zip(cur_mut_set, cur_ann_list))
                    while 'nan' in cur_mut_set: cur_mut_set.remove('nan')
                    valid_mut_dict = dict(zip(cur_mut_set, [1 for _ in range(len(cur_mut_set))]))
                elif '_variants.tsv' in file:
                    df = pd.read_csv(file, sep='\t').applymap(str)
                    cur_mut_set = list(df.AMINOACID_CHANGE)
                    cur_nt_list = list(df.MUTATION)
                    valid_nt_dict = dict(zip(cur_mut_set, cur_nt_list))
                    cur_gene_list = list(df.GENE)
                    valid_gene_dict = dict(zip(cur_mut_set, cur_gene_list))
                    cur_ann_list = list(df.VARIANT_TYPE)
                    valid_ann_dict = dict(zip(cur_mut_set, cur_ann_list))
                    while 'nan' in cur_mut_set: cur_mut_set.remove('nan')
                    valid_mut_dict = dict(zip(cur_mut_set, [1 for _ in range(len(cur_mut_set))]))
            else:
                if 'ann.csv' in file:
                    df = pd.read_csv(file).applymap(str)
                    cur_mut_set = df.AMINO_ACID_CHANGE
                    for mut in cur_mut_set:
                        if mut in valid_mut_dict.keys():
                            valid_mut_dict[mut] += 1
                        elif mut not in valid_mut_dict.keys() and mut != 'nan':
                            valid_mut_dict[mut] = 1
                        if mut not in valid_ann_dict.keys():
                            condition = df['AMINO_ACID_CHANGE'] == mut
                            valid_ann_dict[mut] = df.iloc[df.index[condition], 3].values[0]
                        if mut not in valid_gene_dict.keys():
                            condition = df['AMINO_ACID_CHANGE'] == mut
                            valid_gene_dict[mut] = df.iloc[df.index[condition], 1].values[0]
                        if mut not in valid_nt_dict.keys():
                            condition = df['AMINO_ACID_CHANGE'] == mut
                            valid_nt_dict[mut] = df.iloc[df.index[condition], 0].values[0]
                elif '_variants.tsv' in file:
                    df = pd.read_csv(file, sep='\t').applymap(str)
                    cur_mut_set = df.AMINOACID_CHANGE
                    for mut in cur_mut_set:
                        if mut in valid_mut_dict.keys():
                            valid_mut_dict[mut] += 1
                        elif mut not in valid_mut_dict.keys() and mut != 'nan':
                            valid_mut_dict[mut] = 1
                        if mut not in valid_ann_dict.keys():
                            condition = df['AMINOACID_CHANGE'] == mut
                            valid_ann_dict[mut] = df.iloc[df.index[condition], 5].values[0]
                        if mut not in valid_gene_dict.keys():
                            condition = df['AMINOACID_CHANGE'] == mut
                            valid_gene_dict[mut] = df.iloc[df.index[condition], 4].values[0]
                        if mut not in valid_nt_dict.keys():
                            condition = df['AMINOACID_CHANGE'] == mut
                            valid_nt_dict[mut] = df.iloc[df.index[condition], 1].values[0]

    df = pd.DataFrame(data={'Mutation':list(valid_mut_dict.keys()), "Number_of_samples":list(valid_mut_dict.values())})
    df['Annotation'] = df['Mutation'].map(valid_ann_dict)
    df['Nt_change'] = df['Mutation'].map(valid_nt_dict)
    df['Gene'] = df['Mutation'].map(valid_gene_dict)
    df['%_of_samples'] = round(100*(df['Number_of_samples'] / sample_count),2)
    df = df[['Nt_change', 'Mutation', 'Number_of_samples', '%_of_samples', 'Gene', 'Annotation']]
    df = df[df['Mutation'] != '.']
    return df

total_df = calculate_statistics(file_list)
if args.name_list is not None:
    report_name = args.name_list.split('/')[-1]
    report_name = report_name[:len(report_name) - 4]
    total_df.to_csv(f"{out_dir}{datetime.date.today().strftime('%Y-%m-%d')}_{report_name}_mutation_statistics.csv", index=False)
else:
    total_df.to_csv(f"{out_dir}{datetime.date.today().strftime('%Y-%m-%d')}_{date_1}_{date_2} mutation_statistics.csv", index=False)
