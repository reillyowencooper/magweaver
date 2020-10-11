#!/usr/bin/env python3

import subprocess, os, shutil, argparse
from Bio import SearchIO, SeqIO
import pandas as pd

# TODO: Something is funky with the file locations - prodigal and HMMsearch both just output directly to the mag_dir instead of the output_dir like they're supposed to
# TODO: Everything is working up until the extract_seqs step of opening the multiple sequence FASTA - not sure what's happening

class SearchSingleHMM(object):
    '''Searches a suite of metagenome-assembled genomes to find an HMM of interest'''
    def __init__(self, mag_dir, output_dir, hmm, msa_option, phylo_option, evalue=1e-10):
        self.mag_dir = os.path.abspath(mag_dir)
        self.output_temp = output_dir
        self.hmm = os.path.abspath(hmm)
        self.msa_option = msa_option.lower()
        self.phylo_option = phylo_option.lower()
        self.evalue = str(evalue)
        self.mag_list = []
        self.cds_list = []
        self.hmm_list = []
        
    def get_mags(self):
        for mag in os.listdir(self.mag_dir):
            if mag.endswith('.fa') or mag.endswith('.fna'):
                mag_path = os.path.join(self.mag_dir, mag)
                self.mag_list.append(mag_path)
            
    def create_output_folder(self):
        output_path = os.path.join(self.mag_dir, self.output_temp)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        self.output_dir = output_path
            
    def predict_cds(self):
        for mag in self.mag_list:
            mag_name = os.path.splitext(mag)[0]
            output_path = os.path.join(self.output_dir, mag_name + ".faa")
            prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta']
            if not os.path.exists(output_path):
                subprocess.run(prodigal_cmd)
            self.cds_list.append(output_path)
            
    def search_hmm(self):
        for mag in self.cds_list:
            mag_name = os.path.splitext(mag)[0]
            output_path = os.path.join(self.output_dir, mag_name + ".tab")
            hmmsearch_cmd = ['hmmsearch', '-E', self.evalue, '--tblout', output_path, self.hmm, mag]
            if not os.path.exists(output_path):
                subprocess.run(hmmsearch_cmd)
            self.hmm_list.append(output_path)
            
    def parse_hmmoutput(self):
        hit_ids = {}
        for mag in self.hmm_list:
            mag_name = os.path.splitext(mag)[0]
            hit_list = []
            with open(mag, "r") as tblout:
                for hit_res in SearchIO.parse(tblout, 'hmmer3-tab'):
                    hits = hit_res.hits
                    if len(hits) > 0:
                        for i in range(0, len(hits)):
                            hit_id = hits[i].id
                            hit_list.append(hit_id)
            hit_ids[mag_name] = hit_list
        return hit_ids
    
    def extract_seqs(self, hit_ids):
        self.mag_alignment_path = os.path.join(self.output_dir, 'hmm_mag_hits.fa')
        for mag in self.cds_list:
            mag_name = os.path.splitext(mag)[0]
            contig_seqs = {}
            for seqrecord in SeqIO.parse(mag, 'fasta'):
                contig_seqs[seqrecord.id] = seqrecord.seq
            hit_subset = {key: value for key, value in hit_ids.items() if key == mag_name}
            contig_id_list = [i for i in hit_subset.values()]
            contig_names = []
            for list in contig_id_list:
                for item in list:
                    contig_names.append(item)
            for contig_name in contig_names:
                for seqrecord_id, seqrecord_seq in contig_seqs.items():
                    if contig_name == seqrecord_id:
                        with open(self.mag_alignment_path, 'r') as output:
                            output.write(">" + mag_name + " | " + seqrecord_id + "\n" + str(seqrecord_seq) + "\n")
                                

    def align_seqs(self):
        self.msa_path = os.path.join(self.output_dir, os.path.splitext(self.hmm)[0] + '_msa.aln')
        if self.msa_option == "muscle":
            muscle_cmd = ['muscle', '-in', self.mag_alignment_path, '-clwout', self.msa_path]
            subprocess.run(muscle_cmd)
        elif self.msa_option == "clustal":
            clustal_cmd = ['clustalo', '-i', self.mag_alignment_path, '-o', self.msa_path]
            subprocess.run(clustal_cmd)
            
    def build_phylogeny(self):
        hmm_name = os.path.splitext(self.hmm)[0]
        self.phylogenetic_tree = os.path.join(self.output_dir, hmm_name + "_tree.tre")
        # Still need to figure out how to get PhyML output in correct location
        #if self.phylo_option == "phyml":
            #phyml_cmd = ['phyml', '-i', self.msa_path, '-d', 'aa', '-b', '-1']
            #subprocess.run(phyml_cmd)
        if self.phylo_option == "raxml":
            raxml_cmd = ['raxmlHPC', '-f', 'a', '-p', 99999, '-m', 'PROTGAMMAAUTO', '-x', 99999, '-N', 100, '-s', self.msa_path, '-n', self.phylogenetic_tree]
            subprocess.run(raxml_cmd)
        elif self.phylo_option == "fasttree":
            fasttree_cmd = ['FastTree', '-out', self.phylogenetic_tree, self.msa_path]
            subprocess.run(fasttree_cmd)
        shutil.copy(self.phylogenetic_tree, self.output_dir)
            
    def create_hmm_abundance_df(self):
        abundance_df = pd.DataFrame()
        for hmmout in self.hmm_list:
            mag_name = os.path.splitext(hmmout)[0]
            with open(hmmout, 'r') as tblout:
                lines_for_df = list()
                for line in tblout:
                    if line.startswith('#'):
                        continue
                else:
                    line = line.split()
                    line_fixed = line[:18] + [' '.join(line[18:])]
                    lines_for_df.append(line_fixed)
            hit_dict = {}
            for line in lines_for_df:
                contig_name = str(line[0])
                hit_evalue = float(line[4])
                if hit_evalue >= self.evalue:
                    hit_dict[mag_name] = contig_name, hit_evalue
            abundance_df = abundance_df.append(hit_dict, ignore_index = True)
        abundance_df = abundance_df.rename(index = {0: "mag", 1: "contig"})
        return abundance_df
    
    def create_abundance_summary(self, abundance_df):
        abundance_by_mag = abundance_df.groupby(["mag"]).count()
        output_filename = os.path.join(self.output_dir, "hmm_abundance_by_mag.csv")
        abundance_by_mag.to_csv(output_filename)        
        
    def workhorse(self):
        self.get_mags()
        self.create_output_folder()
        #self.create_intermediate_dirs()
        self.predict_cds()
        self.search_hmm()
        hitids = self.parse_hmmoutput()
        self.extract_seqs(hitids)
        self.align_seqs()
        self.build_phylogeny()
        ab_df = self.create_hmm_abundance_df()
        self.create_abundance_summary(ab_df)
        

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir')
    parser.add_argument('-o', '--output_dir')
    parser.add_argument('-p', '--profile_hmm')
    parser.add_argument('-m', '--msa_option')
    parser.add_argument('-t', '--tree_option')
    parser.add_argument('-e', '--e_value')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    magsearch = SearchSingleHMM(args.input_dir, args.output_dir, args.profile_hmm, args.msa_option, args.tree_option, args.e_value)
    magsearch.workhorse()

if __name__ == "__main__":
    main()