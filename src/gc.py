import os, subprocess, shutil, logging
import pandas as pd
import numpy as np
from Bio import SeqUtils, SearchIO, SeqIO

class MagGC(object):
    '''Evalutes a single MAG for GC content'''
    def __init__(self, mag, outdir):
        '''Inits.
        mag: str
            Path to MAG of interest
        outdir: str
            Directory to output erroneous contig file
        '''
        self.mag = mag
        self.outdir = outdir
        
    def create_logger(self):
        self.logger = logging.getLogger(__file__)
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('gc.log')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        
    def retrieve_contig_len(self):
        '''Creates dict of contig name: contig length'''
        length_dict = {}
        self.logger.info('Retrieving contig lengths')
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            seqid = seqrecord.id
            seqlen = len(seqrecord.seq)
            length_dict[seqid] = seqlen
        return length_dict
    
    def retrieve_contig_gc(self):
        '''Creates a dict of contig name: contig GC content'''
        gc_dict = {}
        self.logger.info('Retrieving contig GC content')
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            seqid = seqrecord.id
            seqgc = SeqUtils.GC(seqrecord.seq)
            gc_dict[seqid] = seqgc
            
    def find_mean_gc(self, length_dict, gc_dict):
        '''Finds weighted mean of GC content across MAG
        length_dict: dict
            Dictionary of contig name: contig length
        gc_dict: dict
            Dictionary of contig name: contig GC content
            
        NOTE: Chose to do it this way because longer contigs are more likely to be assembled correctly OR
            more important for binning, so their GC content has priority
        '''
        self.logger.info('Calculating weighted GC content mean')
        lengths = [length for length in length_dict.values()]
        gc = [gc for gc in gc_dict.values()]
        mean = np.average(gc, weights = lengths)
        return mean
    
    def find_std_gc(self, gc_dict):
        '''Finds standard deviation of GC content across contigs'''
        self.logger.info('Finding GC standard deviation')
        gc = [gc for gc in gc_dict.values()]
        std = np.std(gc)
        return std
    
    def contig_diff(self, mean_gc, gc_dict):
        '''Creates a dict of contig name: contig GC difference from weighted mean
        mean_gc: float
            Value for weighted GC mean across MAG
        gc_dict: dict
            Dictionary of contig name: contig GC content
        '''
        self.logger.info('Calculating deviation from weighted mean per contig')
        diff_dict = {}
        for contig, gc in gc_dict.items():
            diff_dict[contig] = abs(mean_gc - gc)
            
    def identify_erroneous_contigs(self, diff_dict, mean_gc, std_gc):
        '''Examines each contig for GC signature outside of standard deviation
        diff_dict: dict
            Dictionary of contig name: contig GC difference from weighted mean
        mean_gc: float
            Value for weighted GC mean across MAG
        std_gc: float
            Value for standard deviation of GC across MAG
        Returns:
            DataFrame of Contig, GC difference from weighted mean
        '''
        self.logger.info('Identifying contigs with GC content outside of standard deviation')
        erroneous_contigs = {}
        min_gc = mean_gc - std_gc
        max_gc = mean_gc + std_gc
        for contig, gc_diff in diff_dict.items():
            if (gc_diff > max_gc) or (gc_diff < min_gc):
                erroneous_contigs[contig] = gc_diff
        err_df = pd.DataFrame(list(erroneous_contigs.items()), columns = ['Contig', 'GC Diff'])
        return err_df
    
    def write_df(self, err_df, outfile):
        '''Writes DataFrame to file'''
        self.logger.info('Writing erroneous contig dataframe to ' + outfile)
        err_df.to_csv(outfile, index = False, header = True)
        
    def run(self):
        outfile_loc = os.path.join(self.outdir, os.path.splitext(mag)[0] + "_err_gc.csv")
        length_dict = self.retrieve_contig_len()
        gc_dict = self.retrieve_contig_gc()
        mean_gc = self.find_mean_gc(length_dict, gc_dict)
        std_gc = self.find_std_gc(gc_dict)
        diff_dict = self.contig_diff(mean_gc, gc_dict)
        err_df = self.identify_erroneous_contigs(diff_dict, mean_gc, std_gc)
        self.write_df(err_df, outfile_loc)
        
        