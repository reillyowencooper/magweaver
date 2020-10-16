import os, subprocess, shutil, logging
import pandas as pd
import numpy as np
from Bio import SeqIO, SearchIO

def MagContamChecker(object):
    '''Examines a MAG for contamination using coverage, GC content, tetranucleotide'''
    def __init__(self, 
                 mag, 
                 forward_reads, 
                 reverse_reads, 
                 sorted_bam, 
                 mag_stats,
                 uniprot_mmdb):
        self.mag = os.path.abspath(mag)
        self.forward_reads = os.path.abspath(forward_reads)
        self.reverse_reads = os.path.abspath(reverse_reads)
        self.sorted_bam = os.path.abspath(sorted_bam)
        self.mag_stats = mag_stats
        self.uniprot_mmdb = uniprot_mmdb
    
    def create_mag_mmseqs(self):
        # Create MMseqs2 db from input MAG
        pass
    
    def find_contig_taxonomy(self):
        # Search input MAG db against Uniprot50, assign top hit at currently undecided taxonomic rank
        pass
    
    def check_gc(self):
        # Check contig GC content vs mean
        pass
    
    def check_cov(self):
        # Check contig coverage vs mean
        pass
    
    def check_taxonomy(self):
        # Check contig taxonomy vs mean
        pass
    
    def check_tetranucleotide(self):
        # Check contig tetra vs mean (???)
        pass
    
    def summarize_similarity(self):
        # Highligh contigs that don't look similar
        pass
    