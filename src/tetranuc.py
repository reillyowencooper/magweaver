import os, subprocess, shutil, logging
import pandas as pd
import numpy as np
from Bio import SeqIO
from sklearn import decomposition


def MagTetranuc(object):
    
    def __init__(self, mag, outdir):
        self.mag = mag
        self.outdir = outdir
        
    def retrieve_contig_len(self):
        length_dict = {}
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            seqid = seqrecord.id
            seqlen = len(seqrecord.seq)
            length_dict[seqid] = seqlen
        return length_dict
        
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
    
    def get_tetranucleotide_freq(self):
        contig_tetramer_df_list = []
        for seqrecord in SeqIO.parse(self.mag, "fasta"):
            contig_name = seqrecord.id
            contig_seq = str(seqrecord.seq)
            contig_tetramers = self.count_tetramers(contig_seq, contig_name)
            contig_tetramer_df_list.append(contig_tetramers)
        contig_tetramers = pd.concat(contig_tetramer_df_list)
        return contig_tetramers
    
    def run_pca(self, contig_tetramer_df):
        pca = decomposition.PCA(n_components = 1)
        pca.fit(contig_tetramer_df)
        first_axis = pca.components_[0]
        contig_names = contig_tetramer_df['Contig'].tolist()
        contig_component_dict = dict(zip(contig_names, first_axis))
        return contig_component_dict
    
    def find_mean_component(self, contig_component_dict, len_dict):
        lengths = [length for length in len_dict.values()]
        components = [component for component in contig_component_dict.values()]
        mean = np.average(components, weights = lengths)
        return mean
    
    def find_std_component(self, contig_component_dict):
        components = [component for component in contig_component_dict.values()]
        std = np.std(components)
        return std
    
    def identify_erroneous_contigs(self, contig_component_dict, mean, std):
        erroneous_contigs = {}
        min_comp = mean - std
        max_comp = mean + std
        for contig, comp in contig_component_dict.items():
            if (comp > max_comp) or (comp < min_comp):
                erroneous_contigs[contig] = comp
        err_df = pd.DataFrame(list(erroneous_contigs.items()), columns = ['Contig', 'PCA Value'])
        return err_df
    
    def write_df(self, err_df, outfile):
        err_df.to_csv(outfile, index = False, header = True)
        
    def run(self):
        outfile = os.path.join(self.outdir, os.path.splitext(self.mag)[0] + "_err_tetra.csv")
        len_dict = self.retrieve_contig_len()
        contig_tetramers = self.get_tetranucleotide_freq()
        pca_component_dict = self.run_pca(contig_tetramers)
        mean = self.find_mean_component(pca_component_dict, len_dict)
        std = self.find_std_component(pca_component_dict)
        err_df = self.identify_erroneous_contigs(pca_component_dict, mean, std)
        self.write_df(err_df, outfile)