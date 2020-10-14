import os, subprocess, shutil
import pandas as pd
import numpy as np

class MagClassifier(object):
    
    def __init__(self, mag_dir, database_dir, output_dir):
        self.mag_dir = os.path.abspath(mag_dir)
        self.database_dir = database_dir
        self.output_dir = output_dir
        self.cat_db_dir = os.path.join(self.cat_db, "cat_db/2020-06-18_CAT_database")
        self.cat_db_tax = os.path.join(self.cat_db, "cat_db/2020-06-18_taxonomy")
        
    def create_output_folder(self):
        '''Generates folder to output classified bins'''
        output_path = os.path.join(self.mag_dir, self.output_dir)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        self.output_dir = output_path
        
    def run_cat(self):
        '''Runs CAT/BAT on a set of input bins, then adds taxonomic identity to those classified bins'''
        cat_cmd = ['CAT', 'bins', '-b', self.mag_dir, '-d', self.cat_db_dir, '-t', self.cat_db_tax, '-s', '.fa', '--I_know_what_Im_doing', '--top', 6]
        subprocess.run(cat_cmd)
        cat_names_cmd = ['CAT', 'add_names', '-i', "out.BAT.bin2classification.txt", '-o', 'cat_classification.tab', '-t', self.cat_db_tax]
        subprocess.run(cat_names_cmd)
        shutil.move('cat_classification.tab', os.path.join(self.output_dir, 'cat_classification.tab'))
        
    def parse_cat_taxonomy(self):
        '''Parses the output CAT/BAT file to find taxonomic identity for each bin
        NOTE: Coded assuming the CAT file looks like it does in the docs - will check once database is downloaded
        '''
        cat_tax_filepath = os.path.join(self.output_dir, 'cat_classification.tab')
        bin_identity_dict = {}
        bin_list = []
        with open(cat_tax_filepath, 'r') as classification:
            for line in classification:
                if not line.startswith('#'):
                    line = line.split()
                    bin_list.append(line)
        for bin in bin_list:
            mag_name = str(bin[0])
            kingdom = str(bin[17])
            bin_identity_dict[mag_name] = kingdom
        return bin_identity_dict
    
    def place_bins_in_tax_folders(self, bin_identity_dict):
        '''Places CAT-id-ed bins into folders for bacteria, archaea, or eukaryotes'''
        bac_bins_dir = os.path.join(self.output_dir, "bacteria_bins")
        arc_bins_dir = os.path.join(self.output_dir, "archaea_bins")
        euk_bins_dir = os.path.join(self.output_dir, "eukaryota_bins")
        for bindir in [bac_bins_dir, arc_bins_dir, euk_bins_dir]:
            if not os.path.exists(bindir):
                os.mkdir(bindir)
        for mag_name, identity in bin_identity_dict.items():
            for mag in os.listdir(self.mag_dir):
                mag_name_indir = os.path.basename(mag)
                if mag_name == mag_name_indir:
                    if identity == "Bacteria":
                        shutil.move(os.path.join(self.mag_dir, mag), os.path.join(bac_bins_dir, mag))
                    elif identity == "Archaea":
                        shutil.move(os.path.join(self.mag_dir, mag), os.path.join(arc_bins_dir, mag))
                    elif identity == "Eukaryota":
                        shutil.move(os.path.join(self.mag_dir, mag), os.path.join(euk_bins_dir, mag))
                    else:
                        continue
    
    #def remove_cat_outputs(self):
        # Remove the extra files that CAT/BAT generates, no need to have them around
        
    #def classify_remaining_bins(self):
        # Attempt to clasify remaining bins not classified by CAT/BAT, not sure how to do this yet or if needed
        
    def workhorse(self):
        '''Workhorse function that performs folder and file handling, CAT/BAT running and parsing
        NOTE: Will eventually include processing of unclassified bins and viruses
        '''
        self.create_output_folder()
        self.run_cat()
        bin_ids = self.parse_cat_taxonomy()
        self.place_bins_in_tax_folders(bin_ids)