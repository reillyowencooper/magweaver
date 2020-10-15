import subprocess, os, shutil, csv, logging
from Bio import SearchIO, SeqIO
import pandas as pd
import src.camag_utilities as utils


class SearchSingleHMM(object):
    '''Searches a suite of metagenome-assembled genomes (MAGs) to find an HMM of interest
    mag_dir: Directory containing fasta files for each MAG to examine
    output_dir: Directory to output all generated files (nested within mag_dir)
    hmm: The profile HMM to search MAGs against
    msa_option: Muscle or Clustal Omega, personal preference
    phylo_option: RaxML or FastTree, personal preference
    e-value: E-value cutoff for HMMsearch, will not keep hits that have a higher e-value than this
    '''
    def __init__(self, mag_dir, output_dir, hmm, msa_option, phylo_option, evalue=1e-10):
        self.mag_dir = os.path.abspath(mag_dir)
        self.output_dir = output_dir
        self.hmm = os.path.abspath(hmm)
        self.msa_option = msa_option.lower()
        self.phylo_option = phylo_option.lower()
        self.evalue = str(evalue)
        self.mag_list = []
        self.cds_list = []
        self.hmm_list = []
    
    def create_logger(self):
        self.logger = logging.getLogger(__file__)
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('single_hmm_search.log')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        
    def get_mags(self):
        '''Gets a list of MAGs from mag_dir
        '''
        for mag in os.listdir(self.mag_dir):
            if mag.endswith('.fa') or mag.endswith('.fna'):
                mag_path = os.path.join(self.mag_dir, mag)
                self.mag_list.append(mag_path)
            
    def create_output_folder(self):
        '''Prepares output folder for file deposit
        '''
        utils.create_dir(self.output_dir)
            
    def predict_cds(self):
        '''Predicts coding sequences in each MAG using Prodigal'''
        for mag in self.mag_list:
            mag_name = os.path.splitext(mag)[0]
            self.logger.info('Predicting coding regions for ' + mag_name)
            output_path = os.path.join(self.output_dir, mag_name + ".faa")
            prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta', '-q']
            if not os.path.exists(mag_name + ".faa"):
                subprocess.run(prodigal_cmd)
            self.cds_list.append(output_path)
            
    def search_hmm(self):
        '''Searches coding regions for given HMM using Hmmsearch'''
        for mag in self.cds_list:
            mag_name = os.path.splitext(mag)[0]
            self.logger.info('Searching ' + mag_name + ' for HMM')
            output_path = os.path.join(self.output_dir, mag_name + ".tab")
            utils.run_hmmsearch(mag, self.evalue, output_path, self.hmm)
            self.hmm_list.append(output_path)
            
    def parse_hmmoutput(self):
        '''Parses HMM table output and creates dict of MAG:contig names'''
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
        '''Takes in dict of MAG:contig names, then writes a FASTA file with all genome FASTA hits and MAG-specific headers'''
        self.mag_alignment_path = os.path.join(self.mag_dir, 'hmm_mag_hits.fa')
        if not os.path.exists(self.mag_alignment_path):
            seqfile = open(self.mag_alignment_path, "x")
            self.logger.info('Building multiple sequence FASTA with all HMM hits from MAGs')
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
                            with open(self.mag_alignment_path, 'a') as output:
                                output.write(">" + os.path.basename(mag) + "_" + seqrecord_id + "\n" + str(seqrecord_seq) + "\n")
                                
    def align_seqs(self):
        '''Aligns sequence file using either MUSCLE or Clustal Omega'''
        hmm_name = os.path.splitext(os.path.basename(self.hmm))[0]
        self.msa_path = os.path.join(self.mag_dir, hmm_name + '_msa.aln')
        if not os.path.exists(self.msa_path):
            if self.msa_option == "muscle":
                self.logger.info('Generating multiple sequence alignment using MUSCLE')
                muscle_cmd = ['muscle', '-in', self.mag_alignment_path, '-clwout', self.msa_path]
                subprocess.run(muscle_cmd)
            elif self.msa_option == "clustal":
                self.logger.info('Generating multiple sequence alignment using Clustal Omega')
                clustal_cmd = ['clustalo', '-i', self.mag_alignment_path, '-o', self.msa_path]
                subprocess.run(clustal_cmd)
            
    def build_phylogeny(self):
        '''Builds phylogenetic tree using RAxML or FastTree (and PhyML once I figure it out)'''
        hmm_name = os.path.splitext(os.path.basename(self.hmm))[0]
        self.phylogenetic_tree = os.path.join(self.mag_dir, hmm_name + "_tree.tree")
        if not os.path.exists(self.phylogenetic_tree):
            # Still need to figure out how to get PhyML output in correct location
            #if self.phylo_option == "phyml":
                #phyml_cmd = ['phyml', '-i', self.msa_path, '-d', 'aa', '-b', '-1']
                #subprocess.run(phyml_cmd)
            if self.phylo_option == "raxml":
                self.logger.info('Generating phylogeny with RAxML')
                raxml_cmd = ['raxmlHPC', '-f', 'a', '-p', 99999, '-m', 'PROTGAMMAAUTO', '-x', 99999, '-N', 100, '-s', self.msa_path, '-n', self.phylogenetic_tree]
                subprocess.run(raxml_cmd)
            elif self.phylo_option == "fasttree":
                self.logger.info('Generating phylogeny with FastTree')
                fasttree_cmd = ['FastTree', '-out', self.phylogenetic_tree, self.msa_path]
                subprocess.run(fasttree_cmd)
            
    def create_hmm_abundance_csv(self):
        '''Generates a table of MAG:HMM occurrence from Hmmsearch table'''
        self.logger.info('Summarizing number of HMM hits across MAGs')
        mag_abundance_dict = {}
        abundance_df = pd.DataFrame()
        for hmmout in self.hmm_list:
            mag_name = os.path.basename(hmmout)
            with open(hmmout, 'r') as tblout:
                contigs = list()
                for line in tblout:
                    if line.startswith('#'):
                        continue
                    else:
                        line = line.split()
                        line_fixed = line[:18] + [' '.join(line[18:])]
                        contig_name = str(line_fixed[0])
                        contigs.append(contig_name)
            mag_abundance_dict[mag_name] = len(contigs)
        output_filename = os.path.join(self.mag_dir, "hmm_abundance_by_mag.csv")
        with open(output_filename, 'w') as magout:
            w = csv.writer(magout)
            w.writerows(mag_abundance_dict.items())
                 
    def move_files_to_output(self):
        '''Cleans up mag_dir and puts output files from each step in the correct place'''
        self.logger.info('Arranging output files')
        hmm_name = os.path.splitext(os.path.basename(self.hmm))[0]
        prodigal_dir = os.path.join(self.output_dir, "predicted_cds")
        hmm_dir = os.path.join(self.output_dir, "hmm_hits")
        msa_dir = os.path.join(self.output_dir, "msa")
        if not os.path.exists(prodigal_dir):
            os.mkdir(prodigal_dir)
        if not os.path.exists(hmm_dir):
            os.mkdir(hmm_dir)
        if not os.path.exists(msa_dir):
            os.mkdir(msa_dir)
        for filename in os.listdir(self.mag_dir):
            if filename.endswith(".faa"):
                shutil.move(os.path.join(self.mag_dir, filename), os.path.join(prodigal_dir, filename))
            elif filename.endswith(".tab"):
                shutil.move(os.path.join(self.mag_dir, filename), os.path.join(hmm_dir, filename))
        shutil.move(os.path.join(self.mag_dir, hmm_name + '_msa.aln'), os.path.join(self.mag_dir, msa_dir, hmm_name + '_msa.aln'))
        shutil.copyfile(os.path.join(self.mag_dir, hmm_name + '_tree.tree'), os.path.join(self.output_dir, hmm_name + '_tree.tre'))
        shutil.move(os.path.join(self.mag_dir, hmm_name + '_tree.tree'), os.path.join(self.mag_dir, msa_dir, hmm_name + '_tree.tre'))
        shutil.move(os.path.join(self.mag_dir, "hmm_abundance_by_mag.csv"), os.path.join(self.output_dir, "hmm_abundance_by_mag.csv"))
        
    def workhorse(self):
        self.create_logger()
        self.get_mags()
        self.create_output_folder()
        self.predict_cds()
        self.search_hmm()
        hitids = self.parse_hmmoutput()
        self.extract_seqs(hitids)
        self.align_seqs()
        self.build_phylogeny()
        self.create_hmm_abundance_csv()
        self.move_files_to_output()
