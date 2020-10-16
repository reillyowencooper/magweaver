import os, subprocess, shutil, logging
import pandas as pd
import numpy as np
import src.camag_utilities as utils
from Bio import SearchIO, SeqIO, SeqUtils
import pysam

# TODO: Add logging

class MagAligner(object):
    '''Given a MAG and a set of forward/reverse reads, aligns reads and generates stats'''
    def __init__(self, mag, forward_reads, reverse_reads, keep_bam=False):
        self.mag = os.path.abspath(mag)
        self.mag_name = os.path.splitext(os.path.basename(mag))[0]
        self.forward_reads = os.path.abspath(forward_reads)
        self.reverse_reads = os.path.abspath(reverse_reads)
        self.tmp_dir = os.path.join(os.path.dirname(self.mag), "tmp")
        self.keep_bam = keep_bam
        
    def make_tmp_dir(self):
        utils.create_dir(self.tmp_dir)
        
    def index_mag(self):
        index_cmd = ['bwa', 'index', self.mag]
        subprocess.run(index_cmd)
        
    def align_reads(self):
        self.sam_file = os.path.join(self.tmp_dir, self.mag_name + ".sam")
        bwamem_cmd = ['bwa', 'mem', self.mag, self.forward_reads, self.reverse_reads]
        if not os.path.exists(self.sam_file):
            with open(self.sam_file, "w") as outfile:
                subprocess.run(bwamem_cmd, stdout=outfile)
            
    def sam_to_bam(self):
        self.bam_file = os.path.join(self.tmp_dir, self.mag_name + ".bam")
        if not os.path.exists(self.bam_file):
            pysam.view("-bS", "-o", self.bam_file, self.sam_file, catch_stdout = False)
    
    def sort_bam(self):
        self.sorted_bam = os.path.join(self.tmp_dir, self.mag_name + "_sorted.bam")
        if not os.path.exists(self.sorted_bam):
            pysam.sort("-o", self.sorted_bam, self.bam_file)
            
    def get_num_contigs(self):
        contig_headers = []
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            contig_headers.append(seqrecord.id)
        num_contigs = len(contig_headers)
        return num_contigs
        
    def get_contig_gc(self):
        contig_gc = {}
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            contig_name = seqrecord.id
            contig_gc_content = SeqUtils.GC(seqrecord.seq)
            contig_gc[contig_name] = contig_gc_content
        gc_df = pd.DataFrame(list(contig_gc.items()), columns = ['Contig', 'GC Content'])
        return gc_df
        
    def get_contig_coverage_depth(self):
        depth = pysam.depth(self.sorted_bam)
        lines_for_df = [line.split('\t') for line in depth.split('\n')]
        coverage_df = pd.DataFrame(lines_for_df, columns=['Contig','Position','Depth'], dtype=float).dropna()
        contig_cov_mean = coverage_df.groupby(['Contig'])['Depth'].mean().reset_index(name = "Mean Coverage")
        return contig_cov_mean

    def get_contig_tetranucleotide_freq(self):
        contig_tetramer_df_list = []
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            contig_name = seqrecord.id
            contig_seq = str(seqrecord.seq)
            contig_tetramers = self.count_tetramers(contig_seq, contig_name)
            contig_tetramer_df_list.append(contig_tetramers)
        contig_tetramers = pd.concat(contig_tetramer_df_list)
        return contig_tetramers
        
    def count_tetramers(self, seq, seqname):
        tetramers = {}
        for i in range(len(seq) - 3):
            tetramer = seq[i:i+4]
            if tetramer in tetramers:
                tetramers[tetramer] += 1
            else:
                tetramers[tetramer] = 1
        total_tetramers = float(sum(tetramers.values()))
        for tetra in tetramers.keys():
            tetramers[tetra] = float(tetramers[tetra])/total_tetramers
        tetra_df = pd.DataFrame(tetramers, index = [0])
        tetra_df['Contig'] = seqname
        return tetra_df
    
    def build_complete_stats_df(self, tetra_df, depth_df, gc_df):
        tetra_depth_df = tetra_df.merge(depth_df, how = "outer", on = "Contig")
        full_df = tetra_depth_df.merge(gc_df, how = "outer", on = "Contig")
        return full_df
        
    def write_files(self):
        if self.keep_bam == True:
            output_path = os.path.splitext(self.mag)[0] + '_sorted.bam'
            shutil.move(self.sorted_bam, output_path)
    
    def remove_tmp(self):
        index_suffixes = ('.amb', '.ann', '.bwt', '.pac', '.sa')
        for filename in os.path.dirname(self.mag):
            if filename.endswith(index_suffixes):
                os.remove(filename)
        shutil.rmtree(self.tmp_dir)
