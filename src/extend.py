import os, subprocess, shutil, logging
import pandas as pd
import numpy as np
import src.utilities as utils
from Bio import SearchIO, SeqIO, SeqUtils
import pysam

# TODO: Add logging

class MagExtender(object):
    '''Given a MAG and a set of forward/reverse reads, aligns reads and attempt to extend at ends of each contig'''
    def __init__(self, mag, forward_reads, reverse_reads):
        self.mag = mag
        self.mag_name = os.path.splitext(os.path.basename(mag))[0]
        self.forward_reads = forward_reads
        self.reverse_reads = reverse_reads
        self.tmp_dir = os.path.join(os.path.dirname(self.mag), "tmp")
        
        
    def make_tmp_dir(self):
        utils.create_dir(self.tmp_dir)
        
    def split_mag(self):
        '''Splits MAG bin into individual contigs'''
        self.contig_fas = []
        for seqrecord in SeqIO.parse(self.mag):
            seqid = seqrecord.id
            sequence = seqrecord.seq
            with open(os.path.join(self.tmp_dir, seqid), 'w') as cfile:
                cfile.write('>' + str(id) + '\n' + str(sequence))
                self.contig_fas.append(os.path.join(self.tmp_dir, seqid))
            
    def index_contig(self):
        '''Indexes each contig separately using BWA'''
        for cfile in self.contig_fas:
            index_cmd = ['bwa', 'index', cfile]
            subprocess.run(index_cmd)    
        
    def align_reads_to_contig(self, contig_file):
        samfile = os.path.join(self.tmp_dir, os.path.basename(contig_file) + ".sam")
        bwamem_cmd = ['bwa', 'mem', contig_file, self.forward_reads, self.reverse_reads]
        if not os.path.exists(samfile):
            with open(samfile, "w'") as outfile:
                    subprocess.run(bwamem_cmd, stdout=outfile)
                       
    def extract_all_mapped(self, samfile):
        sam_name = os.path.splitext(os.path.basename(samfile))[0]
        # Reads with both mates mapping
        both_outfile = os.path.join(self.tmp_dir, sam_name + "_both.bam")
        both_sorted_outfile = os.path.join(self.tmp_dir, sam_name + "_both_sorted.bam")
        pysam.view('-bS', '-f', '12', '-o', both_outfile, samfile, catch_stdout = False)
        pysam.sort('-o', both_sorted_outfile, both_outfile)
        # Reads with only forward mapping
        read_only_outfile = os.path.join(self.tmp_dir, sam_name + "_read.bam")
        read_sorted_outfile = os.path.join(self.tmp_dir, sam_name + "_read_sorted.bam")
        pysam.view('-bS', '-f', '8', '-F', '4', '-o', read_only_outfile, samfile, catch_stdout = False)
        pysam.sort('-o', read_sorted_outfile, read_only_outfile)
        # Reads with only mate mapping
        mate_only_outfile = os.path.join(self.tmp_dir, sam_name + "_mate.bam")
        mate_sorted_outfile = os.path.join(self.tmp_dir, sam_name + "_mate_sorted.bam")
        pysam.view('-bS','-f', '4', '-F', '8', '-o', read_only_outfile, samfile, catch_stdout = False)
        pysam.sort('-o', mate_sorted_outfile, mate_only_outfile)  

    def merge_bams(self, both, mate1, mate2):
        merged_outfile = os.path.join(self.tmp_dir, self.mag_name + "_merged.bam")
        pysam.merge(merged_outfile, both, mate1, mate2)
            
    def write_bam_to_fq(self, bam, outfile):
        bedtools_cmd = ['bamToFastq', '-i', bam, '-fq', outfile + '_forward.fq', '-fq2', outfile + '_reverse.fq']
        subprocess.run(bedtools_cmd)
    
    def assemble_reads(self, forward, reverse):
        pass # TODO: Create contig assembly

        
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
