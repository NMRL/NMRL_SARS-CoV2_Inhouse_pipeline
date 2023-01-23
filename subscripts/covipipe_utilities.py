
import pandas as pd, gzip, re, os, json, sys, subprocess as sp
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.dirname(Path(__file__).absolute())))
from subscripts.src.utilities import Housekeeper as hk

class covipipe_housekeeper(hk):
    '''Class extends the standard housekeeper class to implement functions required by specific pipeline'''

    @staticmethod
    def map_replace_column(df:pd.DataFrame, new_dict:dict, column_to_map:str ,column_to_replace:str):
        '''
        Given pandas dataframe,a new_dict dictionary, column to map the new_dict and a column to 
        replace with the new dict, performs the replacement and returns resulting dataframe.
        '''
        df['temp_col'] = df[column_to_map].map(new_dict) #adding numbered sample ids to the sample sheet to be used by downstream module
        del df[column_to_replace] #removing old sample ids
        df.rename(columns={'temp_col':column_to_replace}, inplace=True) #adding new sample ids
        return df

    
    @staticmethod
    def read_plaintext_file(path_to_file:str):
        '''
        Given path to file, reads the contents of the file and returns a list of strings. 
        Plaintext file expected.
        '''
        with open(path_to_file, "r+") as file: 
            contents = file.readlines()
            return list(map(str.strip, contents))


    @staticmethod
    def overwrite_plaintext_file(path_to_file:str, text:str):
        '''
        Given path to a file and a text string, overwrites the contents of the file with text.
        '''
        with open(path_to_file, "w+") as file: file.write(text)

    
    @staticmethod
    def check_fastq_integrity(path_to_file:str):
        '''
        Given path to a fastq file returns True if file is intact, otherwise returns False.
        '''
        try: 
            gzip.open(path_to_file, "rt").read()
        except EOFError:
            return False, path_to_file
        return True, path_to_file



    @staticmethod
    def name_formatter(path_to_fastq:str):
        '''Helper function to remove prefix and suffix patters from fastq file name before processing.'''
        prefix_patterns = [r'CO-[0-9]{5}_LVA[0-9]{3}_'] #Eurofins prefixes
        suffix_patterns = [r'_lib[0-9]{6}'] #Eurofins suffixes
        new_path = path_to_fastq
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
        return path_to_fastq, new_path #save to forward map


    @staticmethod
    def renamer(old_path:str, new_path:str):
        '''Helper function to rename fastq files to illumina format.'''
        try: #attempt renaming
            os.rename(old_path, new_path)
            return f'{old_path} {new_path} OK\n'#return info about successful renaming
        except OSError: #if failed
            return f'{old_path} {new_path} FAILED\n'#return info about successful renaming


class Wrapper():
    '''Toolkit class to store wrapper methods for different tools'''

    _config_dict = hk.read_yaml("./config_files/yaml/config_modular.yaml")
    _fastp_version = sp.run(f'module load singularity && singularity run {_config_dict["fastq_sif"]} fastp --version', stderr=sp.PIPE, shell=True).stderr.decode('utf-8').strip()
    _fastqc_version = sp.run(f'module load singularity && singularity run {_config_dict["multiqc_sif"]} fastq_screen --version 2> /dev/null', stdout=sp.PIPE, shell=True).stdout.decode('utf=8').strip()
    _ivar_version = sp.run(f'module load singularity && singularity run {_config_dict["fastq_sif"]} ivar -v 2> /dev/null', stdout=sp.PIPE, shell=True).stdout.decode('utf=8').split('\n')[0]
    _qualimap_version = sp.run(f'module load singularity && singularity run {_config_dict["qualimap_sif_path"]} qualimap bamqc --version 2> /dev/null', stdout=sp.PIPE, shell=True).stdout.decode('utf-8').split('\n')[3]
    _samtools_version = sp.run(f'module load singularity && singularity run {_config_dict["fastq_sif"]} samtools --version 2> /dev/null', stdout=sp.PIPE, shell=True).stdout.decode('utf=8').split('\n')[0]


    @staticmethod
    def parse_fastp(sample_id:str, path_to_report:str) -> dict:
        '''Serializes report as python dictionary'''
        mri_map = {
            'before_filtering':         ['summary','before_filtering'],
            'after_filtering':          ['summary', 'after_filtering'],
            'filtering_results':        ['filtering_result'],
            'duplication':              ['duplication', 'rate'],
            'insert_size':              ['insert_size'],
            'adapter_cutting':          ['adapter_cutting'],
            'read1_before_filtering':   ['read1_before_filtering'],
            'read2_before_filtering':   ['read2_before_filtering'],
            'read1_after_filtering':    ['read1_after_filtering'],
            'read2_after_filtering':    ['read2_after_filtering'],
            'command':                  ["command"]
        }
        report = hk.read_json_dict(path_to_report)
        data = {mri:hk.find_in_nested_dict(report, mri_map[mri]) for mri in mri_map}

        #filter read_before/after_filtering data
        tags = [
            'total_reads',
            'total_bases',
            'q20_bases',
            'q30_bases',
            'total_cycles',
            'quality_curves',
            'content_curves'
            ]
        filters = ['read1_before_filtering', 'read2_before_filtering', 'read1_after_filtering', 'read2_after_filtering']
        for mri in filters:
            record = dict(zip(tags,[data[mri][tag] for tag in tags]))
            data[mri] = record

        #get sequencing mode string
        seq_mode = f'paired end ({data["read1_before_filtering"]["total_cycles"]} cycles + {data["read2_before_filtering"]["total_cycles"]})'
        data['sequencing_mode'] = seq_mode
        data['fastp_version'] = Wrapper._fastp_version

        with open(f'{os.path.abspath(os.path.dirname(path_to_report))}/{sample_id}_fastp_et.json', 'w+') as f: json.dump(data,f, indent=4)


    @staticmethod
    def parse_fqscreen(sample_id:str, path_to_report:str, target_organism:str='SarsCoV2') -> dict:
        '''Serializes report as python dictionary'''
        df = pd.read_csv(path_to_report, sep='\t')
        new_header = df.iloc[0]
        df.columns = new_header
        df = df[:-1]
        data = {df.index[i][0]:[df.index[i][4],df.index[i][5]] for i in range(1,len(df)) if df.index[i][0] == target_organism}
        data['fastq_screen_version'] = Wrapper._fastqc_version
        with open(f'{os.path.abspath(os.path.dirname(path_to_report))}/{sample_id}_fqscreen_et.json', 'w+') as f: json.dump(data,f, indent=4)


    @staticmethod
    def parse_ivar(sample_id:str, path_to_report:str) -> dict:
        '''Serializes report as python dictionary'''
        # try:
        df = pd.read_csv(path_to_report, sep='\r\n', engine='python')[19:-5].reset_index(drop=True)
        if 'Results:' not in df.iloc[0].values:
            new_header = df.iloc[0].values[0].split('\t')
            df.columns=df.iloc[0].values
            df[new_header] = df[df.columns[0]].str.split('\t', expand = True)
        else:
            new_header = df.iloc[1].values[0].split('\t')
            df.columns=df.iloc[1].values
            df[new_header] = df[df.columns[0]].str.split('\t', expand = True)
        data = df[new_header][1:].set_index('Primer Name').to_dict()
        data['ivar_version'] = Wrapper._ivar_version
        with open(f'{os.path.abspath(os.path.dirname(path_to_report))}/{sample_id}_ivar_et.json', 'w+') as f: json.dump(data,f, indent=4)

                
    @staticmethod
    def parse_qualimap(sample_id:str, path_to_report_cov:str, path_to_report_qual:str) -> dict:
        '''Serializes report as python dictionary'''
        gen_cov_hist = pd.read_csv(path_to_report_cov, sep='\t')
        gen_cov_hist = {gen_cov_hist.columns[0]: list(gen_cov_hist[gen_cov_hist.columns[0]]), gen_cov_hist.columns[1]: list(gen_cov_hist[gen_cov_hist.columns[1]])}
        gen_mapq_hist = pd.read_csv(path_to_report_qual, sep='\t')
        gen_mapq_hist = {gen_mapq_hist.columns[0]: list(gen_mapq_hist[gen_mapq_hist.columns[0]]), gen_mapq_hist.columns[1]: list(gen_mapq_hist[gen_mapq_hist.columns[1]])}
        data = {**gen_cov_hist, **gen_mapq_hist}
        data['qualimap_version'] = Wrapper._qualimap_version
        with open(f'{os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(path_to_report_cov))))}/{sample_id}_qmap_et.json', 'w+') as f: json.dump(data,f, indent=4)


    @staticmethod
    def parse_samtools(sample_id:str, path_to_report:str) -> dict:
        '''Serializes report as python dictionary'''
        with open(path_to_report, 'r+') as f:
            lines = f.readlines()
        data = {'reads_mapped':lines[4].strip().split(' ')[0]}
        data['samtools_version'] = Wrapper._samtools_version
        with open(f'{os.path.abspath(os.path.dirname(path_to_report))}/{sample_id}_samtools_et.json', 'w+') as f: json.dump(data,f, indent=4)


    @staticmethod
    def combine_jsons(sample_id:str, output_path:str, report_paths:list) -> None:
        data = {}
        for path in report_paths:
            with open(path, 'r+') as f:
                data[os.path.basename(path)] = json.load(f)
        
        with open(f'{output_path}{sample_id}_combined.json', 'w+') as f:
            json.dump(data,f, indent=4)
        
        