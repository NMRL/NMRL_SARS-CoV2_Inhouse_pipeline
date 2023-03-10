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
            if re.search(prefix,new_path) is not None: #if pattern is detected in file name
                new_path = re.sub(prefix,"",new_path) #replace
                break #done with prefixes
        #replace suffix
        for suffix in suffix_patterns:
            if re.search(suffix,new_path) is not None: #if pattern is detected in file name
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

