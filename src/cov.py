import os, subprocess, shutil, logging
import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO, SearchIO

class MagCov(object):
    '''Evalutes read coverage in a MAG.'''
    def __init__(self, mag, bam, outdir):
        '''Inits.
        mag: str
            Path to MAG of interest
        bam: str
            Path to BAM file of reads mapped to MAG
        outdir: str
            Directory to output erroneous contig file
        '''
        self.mag = mag
        self.bam = bam
        self.outdir = outdir
        
    def create_logger(self):
        self.logger = logging.getLogger(__file__)
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('cov.log')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        
    def sort_bam(self):
        '''Sorts BAM file for further steps'''
        self.logger.info('Sorting input BAM')
        self.sorted_bam = os.path.join(self.outdir, os.path.splitext(bam)[0] + "_sorted.bam")
        if not os.path.exists(self.sorted_bam):
            pysam.sort("-o", self.sorted_bam, self.bam_file)
            
    def get_contig_coverage(self):
        '''Gets per-contig mean coverage
        Returns:
        cov_dict: dict
            Dictionary of contig name:mean coverage
        '''
        self.logger.info('Calculating coverage depth for ' + self.mag + ' using ' + self.bam)
        depth = pysam.depth(self.sorted_bam)
        lines_for_df = [line.split('\t') for line in depth.split('\n')]
        cov_df = pd.DataFrame(lines_for_df, columns = ['Contig', 'Position', 'Depth'], dtype = float).dropna()
        self.logger.info('Calculating mean coverage per contig')
        contig_cov_df = cov_df.groupby(['Contig'])['Depth'].mean().reset_index(name = 'Mean Coverage')
        cov_dict = dict(zip(contig_cov_df['Contig'], contig_cov_df['Mean Coverage']))
        return cov_dict
    
    def retrieve_contig_len(self):
        '''Creates a dict of contig name:contig length'''
        self.logger.info('Retrieving contig lengths')
        length_dict = {}
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            seqid = seqrecord.id
            seqlen = len(seqrecord.seq)
            length_dict[seqid] = seqlen
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
        self.logger.info('Calculating weighted mean coverage across MAG')
        lengths = [length for length in len_dict.values()]
        cov = [cov for cov in cov_dict.values()]
        mean = np.average(cov, weights = lengths)
        return mean
    
    def find_std_cov(self, cov_dict):
        '''Finds standard deviation of coverage across MAG.'''
        self.logger.info('Calculating coverage standard deviation across MAG')
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
        self.logger.info('Flagging contigs with erroneous coverage')
        erroneous_contigs = {}
        min_cov = mean_cov - std_cov
        max_cov = mean_cov + std_cov
        for contig, cov in cov_dict.items():
            if (cov > max_cov) or (cov < min_cov):
                erroneous_contigs[contig] = cov
        err_df = pd.DataFrame(list(erroneous_contigs.items()), column = ['Contig', 'Coverage'])
        return err_df
    
    def write_df(self, err_df, outfile):
        '''Writes DataFrame to file'''
        self.logger.info('Writing erroneous contig data to ' + outfile)
        err_df.to_csv(outfile, index = False, header = True)
        
    def run(self):
        outfile_loc = os.path.join(self.outdir, os.path.splitext(mag)[0] + "_err_cov.csv")
        self.create_logger()
        self.sort_bam()
        cov_dict = self.get_contig_coverage()
        length_dict = self.retrieve_contig_len()
        mean_cov = self.find_mean_cov(cov_dict, length_dict)
        std_cov = self.find_std_cov(cov_dict)
        err_df = self.identify_erroneous_contigs(mean_cov, std_cov, cov_dict)
        self.write(err_df, outfile_loc)
    