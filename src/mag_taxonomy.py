import os, subprocess
import pandas as pd
from Bio import SeqIO, SearchIO
import src.camag_utilities as utils

class MagTaxonomy(object):
    '''Assigns individual contigs taxonomy using MMSeqs2'''
    def __init__(self, mag, output_path, database_path):
        self.mag = os.path.abspath(mag)
        self.output_dir = output_path
        self.mag_name = os.path.splitext(os.path.basename(mag))[0]
        self.database_path = database_path
        
    def create_outputs(self):
        utils.create_dir(self.output_dir)
        utils.create_dir(os.path.join(self.output_dir, "tmp"))
        self.tmp_dir = os.path.join(self.output_dir, "tmp")
        
    def create_mag_mmseqsdb(self):
        createdb_cmd = ['mmseqs', 'createdb', self.mag, os.path.join(self.output_dir, self.mag_name + '_db'), self.tmp_dir]
        subprocess.run(createdb_cmd)
    
    def get_contig_taxonomy(self):
        taxonomy_cmd = ['mmseqs', 'taxonomy', os.path.join(self.output_dir, self.mag_name + '_db'), self.database_path, os.path.join(self.output_dir, self.mag_name + '_tax'), self.tmp_dir, '--merge-query', 1, '--remove-tmp-files', '--tax-lineage', 1]
        subprocess.run(taxonomy_cmd)
        maketsv_cmd = ['mmseqs', 'createtsv', os.path.join(self.output_dir, self.mag_name + '_db'), os.path.join(self.output_dir, self.mag_name + '_tax'), os.path.join(self.output_dir, 'tax.tsv')]
        subprocess.run(maketsv_cmd)
        
    def parse_taxonomy(self):
        tax = pd.read_csv(os.path.join(self.output_dir, 'tax.tsv'), sep = '\t', header=None, names=['Contig','Acc','Cat','LCA','Full Tax'])
        tax = tax.join(tax['Full Tax'].str.split(';', expand = True))
        return tax
        
        