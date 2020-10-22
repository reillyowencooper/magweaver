import os, subprocess, shutil, logging
import pandas as pd
import numpy as np
import src.utilities as utils
from Bio import SearchIO, SeqIO, SeqUtils, Seq
import pysam

# TODO: Add logging

class MagExtender(object):
    '''Given a MAG and a set of forward/reverse reads, aligns reads and attempt to extend at ends of each scaffold'''
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
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            seqid = seqrecord.id
            sequence = seqrecord.seq
            with open(os.path.join(self.tmp_dir, seqid), 'w') as cfile:
                cfile.write('>' + str(seqid) + '\n' + str(sequence))
                self.contig_fas.append(os.path.join(self.tmp_dir, seqid))
            
    def index_contig(self):
        '''Indexes each contig separately using Bowtie2'''
        for cfile in self.contig_fas:
            index_cmd = ['bowtie2-build', cfile, cfile]
            subprocess.run(index_cmd)    
        
    def align_reads_to_contig(self, contig_file):
        samfile = os.path.join(self.tmp_dir, os.path.basename(contig_file) + ".sam")
        bwamem_cmd = ['bwa', 'mem', contig_file, self.forward_reads, self.reverse_reads]
        if not os.path.exists(samfile):
            with open(samfile, "w") as outfile:
                    subprocess.run(bwamem_cmd, stdout=outfile)
        bamfile = os.path.join(self.tmp_dir, os.path.basename(contig_file) + ".bam")
        if not os.path.exists(bamfile):
            pysam.view('-bS', '-o', bamfile, samfile, catch_stdout=False)
        bamfile_sorted = os.path.join(self.tmp_dir, os.path.basename(contig_file) + "_sorted.bam")
        pysam.sort('-o', bamfile_sorted, bamfile)
        pysam.index(bamfile_sorted)
        
    # Working with BAM files and trying to turn into Fastq has to be the worst thing in the history of the universe maybe ever

        
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
