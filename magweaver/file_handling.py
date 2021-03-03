import sys, os, shutil
from distutils.spawn import find_executable
from Bio import SeqIO, Seq

def does_file_exist(file_path):
    actual_path = os.path.abspath(file_path)
    if not os.path.exists(actual_path):
        raise IOError("'%s' does not seem to exist" % file_path)
    return True

def does_program_exist(program_name):
    if find_executable(program_name) is not None:
        return True
    else:
        raise IOError("'%s' is not installed in your path" % program_name)

class ReadNucFasta(object):
    '''Class for reading and dealing with nucleotide FASTA files'''
    def __init__(self, fasta):
        self.fasta = fasta
        
    def retrieve_contigs(self):
        '''Reads a FASTA and returns contig name: sequence dict'''
        contig_dict = {}
        for seqrecord in SeqIO.parse(self.fasta, "fasta"):
            contig_dict[seqrecord.id] = str(seqrecord.seq)
        return contig_dict
    
    def retrieve_contig_len(self):
        '''Reads a FASTA and returns contig name: contig lenght dict'''
        length_dict = {}
        for seqrecord in SeqIO.parse(self.fasta, "fasta"):
            length_dict[seqrecord.id] = len(seqrecord.seq)
        return length_dict
    
    def retrieve_contig_revcomp(self):
        '''For a nucleotide FASTA, returns a dict of contig name: reverse complement seq dict'''
        revcomp_dict = {}
        for seqrecord in SeqIO.parse(self.fasta, "fasta"):
            revcomp_dict[seqrecord.id] = Seq.reverse_complement(seqrecord.seq)
        return revcomp_dict
    
class BamOps(object):
    pass # Will add functionality later probably