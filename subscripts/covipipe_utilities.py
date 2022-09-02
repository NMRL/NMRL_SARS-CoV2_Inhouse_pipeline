from subscripts.src.utilities import Housekeeper as hk
import pandas as pd

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

