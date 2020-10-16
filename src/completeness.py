import os, subprocess, shutil, logging
import pandas as pd
import numpy as np

class CompletenessChecker(object):
    
    def __init__(self, mag_dir, output_dir, completeness_val=95.0, contam_value=5.0):
        '''Checks all MAGs for completeness/contamination and quantifies number of contigs
        mag_dir: Directory containing MAGs of interest
        output_dir: Where to place output files from CheckM
        completeness_val (float): Minimum completeness to qualify as high quality. For CMAGs, this is 95.0
        contam_val (float): Maximum value for contamination to qualify as high quality. For CMAGs, this is 10.0 - currently set to 5.0
        '''
        self.mag_dir = os.path.abspath(mag_dir)
        self.output_dir = os.path.abspath(output_dir)
        self.completeness_val = completeness_val
        self.contam_val = contam_value
        
    def create_logger(self):
        self.logger = logging.getLogger(__file__)
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('completeness.log')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

    def create_output_dir(self):
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
            
    def run_checkm(self):
        '''Runs CheckM on a given set of MAGs'''
        self.logger.info('Checking MAG completeness and contamination using CheckM')
        self.checkm_outpath = os.path.join(self.output_dir, 'checkm_out.tab')
        checkm_cmd = ['checkm', 'lineage_wf', '-x', 'fa', '-f', self.checkm_outpath, '--tab_table', self.mag_dir, self.output_dir]
        subprocess.run(checkm_cmd)
        
    # Commented out because this probably should be handled separately, since this class is just to check completeness
    #def place_highqual_bins(self, checkm_df):
    #    self.high_quality_bins = []
    #    for index, row in checkm_df.iterrows():
    #        if row['Completeness'] >= self.completeness_val and row['Contamination'] <= self.contam_val:
    #            self.high_quality_bins.append(row['Bin Id'])
    #    for bin in self.high_quality_bins:
    #        for mag in os.listdir(self.mag_dir):
    #            mag_path = os.path.join(self.mag_dir, mag)
    #            shutil.copyfile(mag_path, os.path.join(self.output_dir, mag))        
    
    def count_contigs(self, mag):
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
        self.logger.info('Checking MAG quality using ' + self.completeness_val + ' as Completeness cutoff and ' +
                         self.contam_val + ' as Contamination cutoff')
        checkm_df = pd.read_csv(self.checkm_outpath, sep='\t')
        conditions = [
            checkm_df['Completeness'] >= self.completeness_val and checkm_df['Contamination'] <= self.contam_val,
            checkm_df['Completeness'] >= self.completeness_val and checkm_df['Contamination'] >= self.contam_val,
            checkm_df['Completeness'] <= self.completeness_val and checkm_df['Contamination'] <= self.contam_val,
            checkm_df['Completeness'] <= self.completeness_val and checkm_df['Contamination'] >= self.contam_val
        ]
        values = ['HQ', 'MQ', 'MQ', 'LQ']
        checkm_df['Quality'] = np.select(conditions, values)
        wanted_cols = checkm_df[['Bin Id', 'Completeness', 'Contamination']]
        df_with_contig_nums = pd.merge(wanted_cols, contig_count_df, on = 'Bin Id', how = 'outer')
        return df_with_contig_nums
    
    def df_to_csv(self, df, filename):
        output_path = os.path.join(self.output_dir, filename)
        df.to_csv(filename)