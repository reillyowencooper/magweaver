import os, subprocess
import pandas as pd
from Bio import SeqIO, SearchIO
from collections import Counter
import src.camag_utilities as utils

class MagTaxonomy(object):
    '''Assigns individual contigs taxonomy using MMSeqs2'''
    def __init__(self, mag, output_path, database_path):
        self.mag = os.path.abspath(mag)
        self.mag_name = os.path.splitext(os.path.basename(mag))[0]
        self.output_dir = output_path
        self.mag_name = os.path.splitext(os.path.basename(mag))[0]
        self.database_path = database_path # NOTE: Probably best to use SwissProt for now
        
    def create_outputs(self):
        utils.create_dir(self.output_dir)
        utils.create_dir(os.path.join(self.output_dir, "tmp"))
        self.tmp_dir = os.path.join(self.output_dir, "tmp")
        
    def predict_cds(self):
        self.aa_path = os.path.join(self.output_dir, self.mag_name + ".faa")
        if not os.path.exists(self.aa_path):
            utils.predict_cds(self.mag, self.aa_path)
        
    def create_mag_mmseqsdb(self):
        self.mag_db = os.path.join(self.output_dir, self.mag_name + "_db")
        createdb_cmd = ['mmseqs', 'createdb', self.aa_path, self.mag_db]
        subprocess.run(createdb_cmd)
    
    def get_contig_taxonomy(self):
        taxonomy_cmd = ['mmseqs', 'taxonomy', self.mag_db, self.database_path, os.path.join(self.output_dir, self.mag_name + '_tax'), self.tmp_dir, '--merge-query', 1, '--remove-tmp-files', '--tax-lineage', 1]
        subprocess.run(taxonomy_cmd)
        maketsv_cmd = ['mmseqs', 'createtsv', os.path.join(self.output_dir, self.mag_name + '_db'), os.path.join(self.output_dir, self.mag_name + '_tax'), os.path.join(self.output_dir, self.mag_name + 'tax.tsv')]
        subprocess.run(maketsv_cmd)
        
    def parse_taxonomy(self):
        tax = pd.read_csv(os.path.join(self.output_dir, 'tax.tsv'), sep = '\t', header=None, names=['Contig','Acc','Cat','LCA','Full Tax'])
        tax = tax.join(tax['Full Tax'].str.split(';', expand = True)).drop(['Acc', 'Cat', 'Full Tax'], axis = 1)
        tax[['Contig Name','ORF']] = tax['Contig'].str.rsplit('_', 1, expand=True)
        tax = tax.drop(['Contig'], axis = 1)
        tax = tax.drop(tax.iloc[:, 10:36], axis = 1)
        tax = tax[tax['LCA'] != 'unclassified']
        return tax
    
    def get_consensus_taxonomy(self, tax_df):
        contig_tax = {}
        contig_names = tax_df['Contig Name'].unique().tolist()
        for contig in contig_names:
            contig_taxlist = []
            for index, row in tax_df.iterrows():
                if row['Contig Name'] == contig:
                    rowlist = [item for item in row.tolist() if item is not None]
                    for item in rowlist:
                        if item.startswith('p_'):
                            contig_taxlist.append(item)
            contig_phylum = Counter(contig_taxlist).most_common(1)[0][0]
            contig_tax[contig] = contig_phylum
        tax_df = pd.DataFrame(list(contig_tax.items()), columns = ['Contig', 'Phylum'])
        return tax_df
        
    def remove_mmseqs_files(self):
        for filename in os.listdir(self.output_dir):
            if not filename == self.mag_name + "tax.tsv":
                os.remove(filename)
                
    
        
        