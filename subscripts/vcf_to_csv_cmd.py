
#!/usr/bin/python3
import sys, allel, pandas as pd, datetime, re

one_letter = {
    'Val':'V',
    'Ile':'I',
    'Leu':'L',
    'Glu':'E',
    'Gln':'Q',
    'Asp':'D',
    'Asn':'N',
    'His':'H',
    'Trp':'W',
    'Phe':'F',
    'Tyr':'Y',
    'Arg':'R',
    'Lys':'K',
    'Ser':'S',
    'Thr':'T',
    'Met':'M',
    'Ala':'A',
    'Gly':'G',
    'Pro':'P',
    'Cys':'C'
    }

def aa_rename(aa_change): 
    '''The function is used to replace all 3-letter-encoded aa to 1-letter encoding in the string provided.'''
    for key in one_letter.keys(): aa_change = re.sub(key,one_letter[key],aa_change) #replace 3-letter codes for 1-letter code in aa_change string
    return aa_change #edited string

def vcf_to_csv(path_to_vcf,output_path=None):
    '''The function is used to generate csv report from annotated vcf file.'''

    callset = allel.vcf_to_dataframe(path_to_vcf, fields="*")  #reading vcf file using scikit-allel package
    ann_names = 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'.split("|") #annotation data to be extracted from vcf
    ann_data = [result.split('|') for result in callset['ANN']] #annotation data fields found in vcf file
    t_ann_data = {col:dat for col, dat in zip(ann_names,list(map(list, zip(*ann_data))))} #map data from ann_data to ann_name in a dict
    keep_list = ' Gene_Name | Annotation | HGVS.p | ERRORS / WARNINGS / INFO'.split("|") # keep data from provided fields
    old_list = 'MUTATION| Gene_Name | Annotation | HGVS.p |P_ERR_MUT_CALL| ERRORS / WARNINGS / INFO'.split("|") #renaming vcf fields for improved readability of mutation report
    new_list = 'MUTATION|GENE|ANNOTATION|AMINO_ACID_CHANGE|P_ERR_MUT_CALL|ERRORS/WARNINGS/INFO'.split("|")
    file = path_to_vcf.split('/')[-1][:-3] #extracting file name from file path

    output_data = pd.DataFrame.from_dict(t_ann_data) #generate df from vcf data
    output_data = output_data.drop(columns=[col for col in output_data if col not in keep_list]) #drop irrelevant columns
    output_data['MUTATION']=callset['REF']+callset['POS'].astype(str)+callset['ALT_1'] #edit mutation data from vcf to be in human-readable format
    output_data['P_ERR_MUT_CALL'] = callset['QUAL'].apply(lambda x:str(10**(x/(-10)))) #compute error probability from phred-score for mutation call
    output_data = output_data['MUTATION| Gene_Name | Annotation | HGVS.p |P_ERR_MUT_CALL| ERRORS / WARNINGS / INFO'.split("|")] #drop columns used only in calculation
    output_data = output_data.rename(columns={old:new for old,new in zip(old_list, new_list)}) #rename columns using new_list
    output_data['AMINO_ACID_CHANGE'] = output_data['AMINO_ACID_CHANGE'].apply(lambda x:x[2:]) #remove scikit-allel data from aa_change column
    output_data['AMINO_ACID_CHANGE'] = output_data['AMINO_ACID_CHANGE'].apply(lambda x:aa_rename(x)) #change 3-letter aa format to 1-letter aa format
    output_data['COVERAGE'] = callset['DP'] #add coverage depth column
    output_data['FREQUENCY'] = round(100*callset['AO_1']/callset['DP'],2) #compute frequency column
    output_data = output_data['MUTATION|GENE|AMINO_ACID_CHANGE|ANNOTATION|COVERAGE|FREQUENCY|P_ERR_MUT_CALL|ERRORS/WARNINGS/INFO'.split("|")] #reorder dataframe columns
    output_data.to_csv(f'{output_path}', header = True, index = False) #generating mutation report for the sample
    
    return f'{file}vcf - processed - {"{0:%Y-%m-%d %H:%M:%S}".format(datetime.datetime.now())}' #stdoutput upon completion

path_to_vcf, output_path = sys.argv[1], sys.argv[2]
print(vcf_to_csv(path_to_vcf, output_path)) #stdoutput
