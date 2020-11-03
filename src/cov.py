import os, subprocess, shutil
import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO
import src.utilities as utils

# NOTE: Might move BAM file handling to a separate file, because I have a feeling I'll need to implement mapping in other places

class MagCov(object):
    '''Evalutes read coverage in a MAG.'''
    def __init__(self, mag, forward_reads, reverse_reads, outdir, tmp_dir):
        '''Inits.
        mag: str
            Path to MAG of interest
        outdir: str
            Directory to output erroneous contig file
        '''
        self.mag = mag
        self.forward_reads = forward_reads
        self.reverse_reads = reverse_reads
        self.outdir = outdir
        self.tmp_dir = tmp_dir
        self.mag_name = os.path.splitext(os.path.basename(self.mag))[0]
        
    def index_mag(self):
        self.index_loc = os.path.join(self.tmp_dir, os.path.basename(self.mag))
        if not os.path.exists(self.index_loc):
            shutil.copyfile(self.mag, self.index_loc)
        index_cmd = ['bwa', 'index', self.index_loc]
        subprocess.run(index_cmd)
    
    def map_reads(self):
        self.bam_loc = os.path.join(self.tmp_dir, self.mag_name + ".bam")
        full_cmd = "bwa mem " + self.tmp_mag_loc + " " + self.forward_reads + " " + self.reverse_reads + " | samtools sort -o " + self.bam_loc + " -"
        subprocess.call(full_cmd, shell = True) # This is not how I want to implement, but best way to get to BAM immediately for right now
        
    def sort_bam(self):
        '''Sorts BAM file for further steps'''
        self.sorted_bam = os.path.join(self.tmp_dir, self.mag_name + "_sorted.bam")
        if not os.path.exists(self.sorted_bam):
            pysam.sort("-o", self.sorted_bam, self.bam_file)
            
    def get_contig_coverage(self):
        '''Gets per-contig mean coverage
        Returns:
        cov_dict: dict
            Dictionary of contig name:mean coverage
        '''
        depth = pysam.depth(self.sorted_bam)
        lines_for_df = [line.split('\t') for line in depth.split('\n')]
        cov_df = pd.DataFrame(lines_for_df, columns = ['Contig', 'Position', 'Depth'], dtype = float).dropna()
        contig_cov_df = cov_df.groupby(['Contig'])['Depth'].mean().reset_index(name = 'Mean Coverage')
        cov_dict = dict(zip(contig_cov_df['Contig'], contig_cov_df['Mean Coverage']))
        return cov_dict
    
    def retrieve_contig_len(self):
        '''Creates a dict of contig name:contig length'''
        length_dict = {}
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            length_dict[seqrecord.id] = len(seqrecord.seq)
        return length_dict
    
    def find_mean_cov(self, cov_dict, len_dict):
        '''Finds weighted mean of coverage across MAG.
        NOTE: This is a little dicey because it's an average of averages, not sure how to deal with that
        cov_dict: dict
            Dictionary of contig name: contig coverage
        len_dict: dict
            Dictionary of contig name: contig length
        Returns:
        mean: float
            Weighted mean of coverage across MAG
        '''
        lengths = [length for length in len_dict.values()]
        cov = [cov for cov in cov_dict.values()]
        mean = np.average(cov, weights = lengths)
        return mean
    
    def find_std_cov(self, cov_dict):
        '''Finds standard deviation of coverage across MAG.'''
        cov = [cov for cov in cov_dict.values()]
        std = np.std(cov)
        return std
    
    def identify_erroneous_contigs(self, mean_cov, std_cov, cov_dict):
        '''Identifies contigs with erroneous mean coverage in MAG.
        mean_cov: float
            Weighted mean of coverage across MAG
        std_cov: float
            Standard deviation of coverage across MAG
        cov_dict: dict
            Dictionary of contig name:contig mean coverage
        Returns:
            DataFrame of Contig, Coverage for erroneous contigs
        '''
        erroneous_contigs = {}
        min_cov = mean_cov - std_cov
        max_cov = mean_cov + std_cov
        for contig, cov in cov_dict.items():
            if (cov > max_cov) or (cov < min_cov):
                erroneous_contigs[contig] = cov
        err_df = pd.DataFrame(list(erroneous_contigs.items()), columns = ['Contig', 'Coverage'])
        return err_df
    
    def write_df(self, err_df, outfile):
        '''Writes DataFrame to file'''
        err_df.to_csv(outfile, index = False, header = True)
        
    def run(self):
        outfile_loc = os.path.join(self.outdir, os.path.splitext(self.mag)[0] + "_err_cov.csv")
        self.index_mag()
        self.map_reads()
        self.sort_bam()
        cov_dict = self.get_contig_coverage()
        length_dict = self.retrieve_contig_len()
        mean_cov = self.find_mean_cov(cov_dict, length_dict)
        std_cov = self.find_std_cov(cov_dict)
        err_df = self.identify_erroneous_contigs(mean_cov, std_cov, cov_dict)
        self.write_df(err_df, outfile_loc)
    