import os, subprocess, shutil
import pandas as pd
from collections import Counter
import src.utilities as utils

class MagMobilome(object):
    '''Checks MAG for eukaryotic mobile elements.
    This problem was recently highlighted in https://msphere.asm.org/content/5/6/e00854-20#ref-8,
    which suggests that eukaryotic sequences are an urgent problem for most MAGs. This uses
    TREP via MMSeqs2 to search for contigs with these elements.'''
    def __init__(self, mag, outdir, tmp_dir, score=100):
        self.mag = mag
        self.outdir = outdir
        self.mag_name = os.path.splitext(os.path.basename(mag))[0]
        self.tmp_dir = tmp_dir
        self.score = score
        
    def create_mag_db(self):
        createdb_cmd = ['mmseqs', 'createdb', self.mag, os.path.join(self.tmp_dir, mag_name + "_db")]
        subprocess.run(createdb_cmd)
        
    def search_trep(self, trep_db):
        search_cmd = ['mmseqs', 'search', os.path.join(self.tmp_dir, mag_name + "_db"), 
                      trep_db, os.path.join(self.tmp_dir, mag_name + "_search_db"), self.tmp_dir, '--search-type', '3']
        subprocess.run(search_cmd)
        convert_cmd = ['mmseqs', 'convertalis', os.path.join(self.tmp_dir, mag_name + "_db"),
                       trep_db, os.path.join(self.tmp_dir, mag_name + "search_db"), 
                       os.path.join(self.tmp_dir, mag_name + "_hits.tab")]
        
    def parse_trep(self):
        mmseqs_cols = ['query_id', 'target_id', 'seq_identity', 'ali_len', 'num_mismatch',
                       'num_gaps', 'dom_start_query', 'dom_end_query',
                       'dom_start_target', 'dom_end_target',
                       'e_value', 'bit_score']
        search_df = pd.read_csv(os.path.join(self.tmp_dir, mag_name + "_hits.tab"),
                                sep='\t',
                                header=None,
                                names = mmseqs_cols)
        real_hits = search_df[search_df['bit_score'] > 30]
        output_df = real_hits[['query_id', 
                               'target_id', 
                               'bit_score']].rename(columns = {'query_id': 'Contig Name',
                                                               'target_id': 'TREP ID',
                                                               'bit_score': 'Bit Score'})
        return output_df
    
    def identify_erroneous_contigs(self, trep_df):
        err_df = trep_df[trep_df['Contig Name'] >= self.score]
        return err_df
    
    def write_err_df(self, err_df, output_file):
        '''Writes DataFrames to file.'''
        err_df.to_csv(output_file, index = False, header = True)
        
    def run(self):
        outfile_loc = os.path.join(self.outdir, self.mag_name + "_err_euk.csv")
        
        
    
        
        
    