from src.utilities import Housekeeper as hk
import pandas as pd

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