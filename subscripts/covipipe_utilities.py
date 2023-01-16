from subscripts.src.utilities import Housekeeper as hk
import pandas as pd, gzip, re, os

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

    @staticmethod
    def parse_fastp(path_to_report:str) -> dict:
        '''Serializes report as python dictionary'''
        mri_map = {
            'failed_avq':   ["filtering_result", "low_quality_reads"],
            'failed_msb':   ["filtering_result", "too_many_N_reads"],
            'failed_cpx':   ["filtering_result", "low_complexity_reads"],
            'failed-len':   ["filtering_result", "too_short_reads"],
            'bcaf_histr1':  ["read1_after_filtering", "content_curves"],
            'bqaf_histr1':  ["read1_after_filtering", "quality_curves"],
            'bcaf_histr2':  ["read2_after_filtering", "content_curves"],
            'bqaf_histr2':  ["read2_after_filtering", "quality_curves"],
            'rpf_r1':       ["read1_after_filtering","total_reads"],
            'rpf_r2':       ["read2_after_filtering","total_reads"],
            'avlnaf_r1':    ["summary", "after_filtering", "read1_mean_length"],
            'avlnaf_r2':    ["summary", "after_filtering", "read2_mean_length"],
            'gc_af':        ["summary", "after_filtering", "gc_content"],
            'q30_af':       ["summary", "after_filtering", "q30_rate"],
            'cmd':          ["command"]
        }
        report = hk.read_json_dict(path_to_report)
        return {mri:hk.find_in_nested_dict(report, mri_map[mri]) for mri in mri_map}
            
    @staticmethod
    def parse_fqscreen(path_to_report:str, target_organism:str='SarsCoV2') -> dict:
        '''Serializes report as python dictionary'''
        df = pd.read_csv(path_to_report, sep='\t')
        new_header = df.iloc[0]
        df.columns = new_header
        df = df[:-1]
        data = {df.index[i][0]:[df.index[i][4],df.index[i][5]] for i in range(1,len(df)) if df.index[i][0] == target_organism}
        return data

    @staticmethod
    def parse_ivar(path_to_report:str) -> dict:
        '''Serializes report as python dictionary'''
        df = pd.read_csv(path_to_report, sep='\n')[19:-5]
        new_header = df.iloc[0].get_values()[0].split('\t')
        df[new_header] = df[df.columns[0]].str.split('\t', expand = True)
        data = df[new_header][1:].set_index('Primer Name').to_dict()
        return data

    @staticmethod
    def parse_qualimap(path_to_report_cov:str, path_to_report_qual:str) -> dict:
        '''Serializes report as python dictionary'''
        gen_cov_hist = pd.read_csv(path_to_report_cov, sep='\t')
        gen_cov_hist = {gen_cov_hist.columns[0]: list(gen_cov_hist[gen_cov_hist.columns[0]]), gen_cov_hist.columns[1]: list(gen_cov_hist[gen_cov_hist.columns[1]])}
        gen_mapq_hist = pd.read_csv(path_to_report_qual, sep='\t')
        gen_mapq_hist = {gen_mapq_hist.columns[0]: list(gen_mapq_hist[gen_mapq_hist.columns[0]]), gen_mapq_hist.columns[1]: list(gen_mapq_hist[gen_mapq_hist.columns[1]])}
        data = {**gen_cov_hist, **gen_mapq_hist}
        return data

    @staticmethod
    def parse_samtools(path_to_report:str) -> dict:
        '''Serializes report as python dictionary'''
        with open(path_to_report, 'r+') as f:
            lines = f.readlines()
        return {'reads_mapped':lines[4].strip().split(' ')[0]}

if __name__ == "__main__":
    print('''
Module to define tool-specific wrappers to be used from snakefile.
    ''')