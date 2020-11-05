import os, subprocess, shutil
import pandas as pd
import numpy as np
import src.utilities as utils

class QualityChecker(object):
    '''Examines directory of MAGs for completeness, contamination, and number of bins to find candidate CMAGs'''
    def __init__(self, mag_dir, output_dir, completeness_val=95.0, contam_value=5.0, num_contigs=10):
        '''Inits.
        mag_dir: Directory containing MAGs of interest
        output_dir: Where to place high quality MAGs
        completeness_val (float): Minimum completeness to qualify as high quality. For CMAGs, this is 95.0
        contam_val (float): Maximum value for contamination to qualify as high quality. For CMAGs, this is 10.0 - currently set to 5.0
        num_contigs (int): Maximum number of contigs to qualify as high quality. For CMAGs, default is 10.
        '''
        self.mag_dir = mag_dir
        self.output_dir = output_dir
        self.completeness_val = completeness_val
        self.contam_val = contam_value
        self.num_contigs = num_contigs
            
    def run_checkm(self, outloc):
        '''Runs CheckM on a given set of MAGs'''
        checkm_cmd = ['checkm', 'lineage_wf', '-x', 'fa', '-f', outloc, '--tab_table', self.mag_dir, self.output_dir]
        subprocess.run(checkm_cmd)     
    
    def count_contigs(self):
        contig_count_dict = {}
        for mag in os.listdir(self.mag_dir):
            mag_path = os.path.join(self.mag_dir, mag)
            mag_name = os.path.basename(mag)
            num_contigs = 0
            with open(mag_path, 'r') as magfile:
                for line in magfile:
                    if line.startswith('>'):
                        num_contigs += 1
            contig_count_dict[mag_name] = num_contigs
        contig_count_df = pd.DataFrame(list(contig_count_dict.items()), columns = ['Bin Id', 'Num Contigs'])
        return contig_count_df

    def summarize_completeness_contam(self, contig_count_df):
        '''Adds column to CheckM result file denoting:
        1. Bins as high (HQ), medium (MQ), or low (LQ) quality
        2. Number of contigs in a bin
        '''
        checkm_df = pd.read_csv(self.checkm_outpath, sep='\t')
        wanted_cols = checkm_df[['Bin Id', 'Completeness', 'Contamination']]
        df_with_contig_nums = pd.merge(wanted_cols, contig_count_df, on = 'Bin Id', how = 'outer')
        df_with_contig_nums['Quality'] = np.where((df_with_contig_nums['Completeness'] >= self.completeness_val) & 
                                                  (df_with_contig_nums['Contamination'] <= self.contam_val) &
                                                  (df_with_contig_nums['Num Contigs'] <= self.num_contigs), 'HQ',
                                                  (np.where((df_with_contig_nums['Completeness'] < self.completeness_val) &
                                                            (df_with_contig_nums['Contamination'] > self.contam_val) &
                                                            (df_with_contig_nums['Num Contigs'] > self.num_contigs), 'LQ', 'MQ')))
        return df_with_contig_nums
    
    def df_to_csv(self, df, filename):
        output_path = os.path.join(self.output_dir, filename)
        df.to_csv(filename)
        
    def copy_hqmags_out(self, summary_df):
        hq_only = summary_df[summary_df['Quality'] == 'HQ']
        hq_names = hq_only['Bin Id'].tolist()
        for mag in os.listdir(self.mag_dir):
            mag_name = os.path.splitext(mag)[0]
            for hqname in hq_names:
                if mag_name == hqname:
                    shutil.copyfile(mag, os.path.join(self.output_dir, mag))
                    
    def run(self):
        checkm_outloc = os.path.join(self.output_dir, "_checkm.tab")
        complete_summary_outloc = os.path.join(self.output_dir, "_summary.csv")
        counts = self.count_contigs()
        summary = self.summarize_completeness_contam(counts)
        self.copy_hqmags_out(summary)
        self.df_to_csv(summary, complete_summary_outloc)
        
        