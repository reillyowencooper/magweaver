import os, subprocess, shutil, sys
import pandas as pd
import numpy as np
from Bio import SeqUtils, SeqIO
from src.file_handling import ReadNucFasta

class MagGC(object):
    '''Evalutes a single MAG for GC content'''
    def __init__(self, mag, outdir):
        '''Inits.
        mag: str
            Path to MAG of interest
        outdir: str
            Directory to output erroneous contig file
        '''
        self.mag = ReadNucFasta(mag).fasta
        self.len_dict = ReadNucFasta(mag).retrieve_contig_len()
        self.outdir = outdir
    
    def retrieve_contig_gc(self):
        '''Creates a dict of contig name: contig GC content'''
        gc_dict = {}
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            gc_dict[seqrecord.id] = SeqUtils.GC(seqrecord.seq)
        return gc_dict
            
    def find_mean_gc(self, length_dict, gc_dict):
        '''Finds weighted mean of GC content across MAG
        length_dict: dict
            Dictionary of contig name: contig length
        gc_dict: dict
            Dictionary of contig name: contig GC content
            
        NOTE: Chose to do it this way because longer contigs are more likely to be assembled correctly OR
            more important for binning, so their GC content has priority
        '''
        lengths = [length for length in length_dict.values()]
        gc = [gc for gc in gc_dict.values()]
        mean = np.average(gc, weights = lengths)
        return mean
    
    def find_std_gc(self, gc_dict):
        '''Finds standard deviation of GC content across contigs'''
        gc = [gc for gc in gc_dict.values()]
        std = np.std(gc)
        return std
            
    def identify_erroneous_contigs(self, gc_dict, mean_gc, std_gc):
        '''Examines each contig for GC signature outside of 2x standard deviation
        diff_dict: dict
            Dictionary of contig name: contig GC difference from weighted mean
        mean_gc: float
            Value for weighted GC mean across MAG
        std_gc: float
            Value for standard deviation of GC across MAG
        Returns:
            DataFrame of Contig, GC content for erroneous contigs
        '''
        erroneous_contigs = {}
        min_gc = mean_gc - 2*std_gc
        max_gc = mean_gc + 2*std_gc
        for contig, gc in gc_dict.items():
            if (gc > max_gc) or (gc < min_gc):
                erroneous_contigs[contig] = gc
        err_df = pd.DataFrame(list(erroneous_contigs.items()), columns = ['Contig', 'GC Content'])
        return err_df
    
    def write_df(self, err_df, outfile):
        '''Writes DataFrame to file'''
        err_df.to_csv(outfile, index = False, header = True)
        
    def run(self):
        outfile_loc = os.path.join(self.outdir, os.path.splitext(self.mag)[0] + "_err_gc.csv")
        gc_dict = self.retrieve_contig_gc()
        mean_gc = self.find_mean_gc(self.len_dict, gc_dict)
        std_gc = self.find_std_gc(gc_dict)
        err_df = self.identify_erroneous_contigs(gc_dict, mean_gc, std_gc)
        self.write_df(err_df, outfile_loc)
        
        