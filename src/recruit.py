import os, subprocess, logging
from Bio import SeqIO, SearchIO, Seq
import pandas as pd

class MagRecruiter(object):
    '''Search metagenome assembly for sections to recruit
    NOTE: This should work after one or two rounds of extension, since newly extended scaffolds 
        can match to previously unmatched contigs
    '''
    def __init__(self, extended_scaffold_fasta, metagenome_assembly, output_dir):
        self.input_seqfile = extended_scaffold_fasta
        self.search_seqfile = metagenome_assembly
        self.output_dir = output_dir
        self.tmp_dir = os.path.join(self.output_dir, "tmp")
        
    def create_mmseqs_db(self, seqfile):
        seqfile_name = os.path.splitext(os.path.basename(seqfile))[0]
        mmseqs_cmd = ['mmseqs', 'createdb', seqfile, seqfile_name + "_db"]
        if not os.path.exists(seqfile + "_db"):
            subprocess.run(mmseqs_cmd)
            
    def search_assembly(self, input_db, assembly_db, output_db):
        self.searched_db = os.path.join(self.output_dir, output_db)
        search_cmd = ['mmseqs', 'search', input_db, assembly_db, self.searched_db, self.tmp_dir, '--search-type', '3']
        if not os.path.exists(self.searched_db):
            subprocess.run(search_cmd)
    
    def convert_searched(self, input_db, assembly_db):
        self.outfile = os.path.join(self.output_dir, "search_output.tab")
        convert_cmd = ['mmseqs', 'convertalis', input_db, assembly_db, self.searched_db, self.outfile]
        if not os.path.exists(self.outfile):
            subprocess.run(convert_cmd) 
        
    def get_mag_scaffold_names(self):
        scaffold_names = []
        for seqrecord in SeqIO.parse(self.input_seqfile):
            scaffold_names.append(seqrecord.id)
        return scaffold_names
    
    def parse_search_output(self, scaffold_names):
        search_df = pd.read_csv(self.outfile, sep='\t', 
                                header=None, 
                                names = ['query_seq','target_seq','seqid',
                                         'ali_len', 'num_mismatch', 'num_gap_openings',
                                         'query_dom_start', 'query_dom_end', 
                                         'target_dom_start', 'target_dom_end',
                                         'e_value', 'bitscore'])
        noself_df = search_df[~search_df['target_seq'].isin(scaffold_names)]
        return noself_df
    
    def identify_links(self, noself_df):
        pass # TODO: Find scaffols that align with the query sequences