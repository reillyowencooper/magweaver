import os, subprocess, shutil
from collections import defaultdict
from Bio import SearchIO, SeqIO
import src.utilities as utils
import pandas as pd


class MagSCG(object):
    '''Searches a MAG for single copy genes, flags redundant HMMs and their contigs.
    This uses the Bacteria HMMs from Albertsen et al., 2013, found at https://github.com/MadsAlbertsen/multi-metagenome/blob/master/R.data.generation/essential.hmm
    NOTE: This is a bacteria-specific HMM set, should probably add Archaea at some point
    '''
    def __init__(self, mag, scg_db, outdir):
        self.mag = mag
        self.mag_name = os.path.splitext(os.path.basename(mag))[0]
        self.scg_db = scg_db
        self.outdir = outdir
        self.tmp = os.path.join(self.outdir, "tmp")
        self.evalue = "1e-15"
        
    def create_tmp(self):
        utils.create_dir(self.tmp)
        
    def get_aas(self):
        self.aas = os.path.join(self.tmp, self.mag_name + ".faa")
        utils.predict_cds(self.mag, self.aas)
        
    def search_scgs(self):
        self.essential = os.path.join(self.tmp, self.mag_name + "_essential.tab")
        utils.run_hmmsearch(self.aas, self.evalue, self.bac, self.scg_db)
        
    def find_potential_redundancy(self):
        hmm_hits = defaultdict(int)
        essential = utils.parse_hmmtbl(self.esssential)
        for line in essential:
            hmm_acc = line[3]
            hmm_hits[hmm_acc] += 1
        flagged_contigs = defaultdict(list)
        flagged_hmms = []
        for hmm_acc, num_hits in hmm_hits.items():
            if num_hits > 1:
                flagged_hmms.append(hmm_acc)
        for hmm in flagged_hmms:
            for line in essential:
                hmm_acc = line[3]
                contig = line[0].rsplit('_',1)[0]
                if hmm_acc == hmm:
                    flagged_contigs[hmm_acc].append(contig)
        err_df = pd.DataFrame(list(flagged_contigs.items()), columns = ['Single Copy Gene', 'Contig Names'])
        return err_df
    
    def write_df(self, err_df, outfile):
        err_df.to_csv(outfile, index = False, header = True)
        
    def remove_tmp(self):
        shutil.rmtree(self.tmp)
        
    def run(self):
        outfile_loc = os.path.join(self.outdir, self.mag_name + "_err_scg.csv")
        self.create_tmp()
        self.get_aas()
        self.search_scgs()
        err_df = self.find_potential_redundancy()
        self.write_df(err_df, outfile_loc)
        self.remove_tmp()
                
                        
                        
                        