import os, subprocess, logging

def create_dir(filepath):
    '''Checks if a directory exists then creates it if not'''
    if not os.path.exists(filepath):
        os.mkdir(filepath)
        
def run_hmmsearch(aa_file, e_value, output_path, hmm):
    '''Runs HMMsearch, given an input file, an HMM, whatever E-value you want to use'''
    hmmsearch_cmd = ['hmmsearch', '-E', e_value, '--tblout', output_path, hmm, aa_file]
    if not os.path.exists(output_path):
        subprocess.run(hmmsearch_cmd)