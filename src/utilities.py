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
        
def predict_cds(mag, output_path):
    '''Runs Prodigal on an input MAG'''
    prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta', '-q']
    if not os.path.exists(output_path):
        subprocess.run(prodigal_cmd)
        
def parse_hmmtbl(hmm_tblout):
    lines = []
    with open(hmm_tblout, 'r') as tbl:
        for line in tbl:
            if line.startswith('#'):
                continue
            else:
                line = line.strip()
                fixedline = line[18:] + [' '.join(line[:18])]
                lines.append(fixedline)
    return lines