import os, subprocess
import pandas as pd
import numpy as np

class MagClassifier(object):
    
    def __init__(self, mag_dir, database_dir, output_dir):
        self.mag_dir = os.path.abspath(mag_dir)
        self.database_dir = database_dir
        self.output_dir = output_dir
        self.cat_db = os.path.join(self.database_dir, "cat_db")
        
    def create_output_folder(self):
        output_path = os.path.join(self.mag_dir, self.output_dir)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        self.output_dir = output_path
        
    def run_cat(self):
        cat_cmd = ['CAT', 'bins', '-b', self.mag_dir, '-d', self.cat_db, '-t', self.output_dir]
        subprocess.run(cat_cmd)
        cat_names_cmd = ['CAT', 'add_names', '-i', os.path.join(self.output_dir, "out.BAT.bin2classification.txt"), '-o', os.path.join(self.output_dir, 'cat_classification.tab'), '-t', self.output_dir]
        subprocess.run(cat_names_cmd)
        
    def parse_cat_taxonomy(self):
        cat_tax_filepath = os.path.join(self.output_dir, 'cat_classification.tab')
        # Coded assuming the CAT file looks like it does in the docs - will check once database is downloaded
        cat_tax = pd.from_csv(cat_tax_filepath, sep = "\t")
        bin_identity_dict = {}
        for index, row in cat_tax.iterrows():
            if row['lineage'] is NaN:
                bin_identity_dict[row['bin']] = 'Unidentified'
            else:
                kingdom = row['classification'].split(';')[0]
                bin_identity_dict[row['bin']] = kingdom
        return bin_identity_dict
    
    def place_bins_in_tax_folders(self):
        # Put CAT-id'd bins into respective bacteria/archaea/eukaryote filters
        
    def classify_remaining_bins(self):
        # Use other DBs to identify the rest of the bins not able to be identified by CAT/BAT
        
        
        