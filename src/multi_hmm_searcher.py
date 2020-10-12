import subprocess, os, shutil, csv
from Bio import SearchIO, SeqIO
import pandas as pd

class SearchSuiteHMM(object):
    
    def __init__(self, mag_dir, hmm_list, hmm_dir, e_value_path, output_dir):
        '''Searches a set of metagenome-assembled genomes (MAGs) against a given set of profile HMMs
        mag_dir: Path to directory containing MAGs of interest
        hmm_list: List of HMMs to search against
        hmm_dir: Directory containing HMMs to search. Structured this way if only a subset of an HMM library is wanted
        e_value_path: Path to file containing a two-column csv file of HMM name, e-value cutoff to search
        output_dir: Where to place output files
        '''
        self.mag_dir = os.path.abspath(mag_dir)
        self.hmms = hmm_list
        self.hmm_dir = os.path.abspath(hmm_dir)
        self.output_dir = output_dir
        self.e_value = e_value_path
        self.mag_paths = []
        self.mag_hmm_dirs = []
        self.cds_paths = []
        
    def create_output_folder(self):
        '''Prepares output folder for file deposit, nested within the MAG dir
        '''
        output_path = os.path.join(self.mag_dir, self.output_dir)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        self.output_dir = output_path
        
    def get_mags(self):
        for mag in os.listdir(self.mag_dir):
            if mag.endswith('.fa') or mag.endswith('.fna'):
                mag_path = os.path.join(self.mag_dir, mag)
                self.mag_paths.append(mag_path)
    
    def get_hmms(self):
        for hmm in os.listdir(self.hmm_dir):
            hmm_name = os.path.splitext(os.path.basename(hmm))[0]
            if hmm_name in self.hmms:
                hmm_path = os.path.join(self.hmm_dir, hmm_name)
                self.hmm_paths.append(hmm_path)
    
    def get_evals(self):
        eval_df = pd.read_csv(self.e_value, names = ["hmm", "e_value"])
        self.hmm_evals = dict(zip(eval_df.hmm, eval_df.e_value))
        
    def predict_cds(self):
        '''Predicts coding sequences in each MAG using Prodigal'''
        for mag in self.mag_paths:
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            output_path = os.path.join(self.output_dir, mag_name + ".faa")
            prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta']
            if not os.path.exists(mag_name + ".faa"):
                subprocess.run(prodigal_cmd)
            self.cds_list.append(output_path)
            
    def search_hmms(self):
        for mag_cds in self.cds_paths:
            mag_name = os.path.splitext(os.path.basename(mag_cds))[0]
            mag_hmmsearch_outpath = os.path.join(self.output_dir, mag_name)
            if not os.path.exists(mag_hmmsearch_outpath):
                os.mkdir(mag_hmmsearch_outpath)
            for hmm in self.hmm_paths:
                hmm_name = os.path.splitext(os.path.basename(hmm))[0]
                for key, value in self.hmm_evals.items():
                    if hmm_name == key:
                        e_value = value
                output_path = os.path.join(mag_hmmsearch_outpath, hmm_name + ".tab")
                hmmsearch_cmd = ['hmmsearch', '-E', e_value, '--tblout', output_path, hmm, mag_cds]
                if not os.path.exists(output_path):
                    subprocess.run(hmmsearch_cmd)
                self.mag_hmm_dirs.append(mag_hmmsearch_outpath)
    
    def parse_hmms(self):
        mag_dict = {}
        for mag_hmm_dir in self.mag_hmm_dirs:
            hit_list = []
            mag_name = os.path.basename(mag_hmm_dir)
            for hmmtbl in os.listdir(mag_hmm_dir):
                hmm_dict = {}  
                hmm_name = os.path.splitext(os.path.basename(hmmtbl))[0]
                with open(hmmtbl, "r") as tblout:
                    for hit_res in SearchIO.parse(tblout, 'hmmer3-tab'):
                        hits = hit_res.hits
                        if len(hits) > 0:
                            for i in range(0, len(hits)):
                                hit_id = hits[i].id
                                hit_list.append(hit_id)
                    hmm_dict[hmm_name] = hit_list
                mag_dict[mag_name] = hmm_dict
        return mag_dict
    
    def create_mag_hmm_abundance_df(self, mag_hmm_dict):
        # Do something with the dict created above