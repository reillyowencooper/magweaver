# Creates classes for individual contigs and for entire MAGs
# Runs basic contig metrics on individual contigs (Contig), including GC/tetramers/length/reverse complement
# Runs contig metrics MAG-wide for efficiency (MAG) , like taxonomy/coverage/core conserved genes/mobilome content
# Runs MAG-wide statistics (MAG), like concensus LCA, tetranucleotide deviation (PCA), standard deviation and weighted mean of GC content, coverage based on contig length

import os, subprocess, shutil
import pysam
import pandas as pd
import numpy as np
from Bio import SeqUtils, SeqIO, Seq
from collections import defaultdict, Counter
from sklearn import decomposition

BASEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCG_DB = os.path.join(BASEPATH, "data/essential.hmm.gz")
TAX_DB = os.path.join(BASEPATH, "data/swissprot")
TREP_DB = os.path.join(BASEPATH, "data/trep")

class Contig(object):
    '''Defines a contig, evaluates contig-wise metrics assuming input is a BioPython SeqRecord'''
    def __init__(self, contig):
        self.contig = contig
        self.contig_dict = self.create_contig_dict()
        
    def create_contig_dict(self):
        contig_dict = defaultdict()
        seq_id = self.contig.id
        seq = self.contig.seq
        length = len(self.contig.seq)
        revcomp = Seq.reverse_complement(self.contig.seq)
        gc = SeqUtils.GC(self.contig.seq)
        tetra = self.count_tetramers(self.contig.seq)
        tetra_proportion = self.tetramer_proportion(tetra)
        contig_dict[seq_id] = {"seq": seq,
                               "len": length,
                               "revcomp": revcomp,
                               "gc": gc,
                               "tetra": tetra,
                               "tetra_proportion": tetra_proportion}
        return contig_dict
        
    # --------- Tetranucleotide frequencies in a sequence --------       
    
    def count_tetramers(self, seq):
        sequence_str = str(seq)
        tetramers = {}
        for i in range(len(sequence_str) - 3):
                tetramer = sequence_str[i:i+4]
                if tetramer in tetramers:
                    tetramers[tetramer] += 1
                else:
                    tetramers[tetramer] = 1
        return tetramers
    
    def tetramer_proportion(self, tetramer_dict):
        total_tetra = float(sum(tetramer_dict.values()))
        tetramer_prop = tetramer_dict
        for tetra in tetramer_prop.keys():
            tetramer_prop[tetra] = float(tetramer_prop[tetra])/total_tetra
        return total_tetra
      
        
class Mag(object):
    '''Defines a metagenome-assembled genome, evaluates for useful metrics'''
    def __init__(self, mag_fasta_loc, mg_forward_reads_loc, mg_reverse_reads_loc, tmp_dir):
        self.mag_name = os.path.basename(mag_fasta_loc)
        self.mag = self.create_seq_dict(mag_fasta_loc)
        self.freads = mg_forward_reads_loc
        self.rreads = mg_reverse_reads_loc
        self.tmp_dir = tmp_dir
        
    # ~~~~~~~~~~~ Create initial sequence dictionary ~~~~~~~~~~~~~~ 
    
    def create_seq_dict(self, fasta):
        '''Sets up sequence dictionary for all contigs in '''
        complete_contig_dict = defaultdict()
        for seqrecord in SeqIO.parse(fasta, "fasta"):
            contig = Contig(seqrecord)
            contig_dict = contig.contig_dict
            complete_contig_dict.update(contig_dict)
        return complete_contig_dict
    
    # ~~~~~~~~~~~~ Contig statistics performed MAG-wide ~~~~~~~~~~~~~~
    
    def craft_mag(self):
        mag = self.create_seq_dict(self.mag)
        self.index_mag()
        self.map_reads()
        mag = self.count_coverage(mag)
        self.predict_cds()
        self.run_hmmsearch()
        mag = self.create_scg_dict(mag)
        self.create_mag_db()
        self.run_taxonomy()
        mag = self.create_tax_dict(mag)
        self.search_trep()
        mag = self.create_mobilome_dict(mag)
        mag = self.create_pca_dict(mag)
        self.mag_contigs = mag

    # --------- Read coverage --------      
    
    def index_mag(self):
        self.index_loc = os.path.join(self.tmp_dir, self.mag_name)
        if not os.path.exists(self.index_loc):
            shutil.copyfile(self.mag, self.index_loc)
        index_cmd = ["bwa", "index", self.index_loc]
        subprocess.run(index_cmd)
        
    def map_reads(self):
        self.bam = os.path.join(self.tmp_dir, self.mag_name + ".bam")
        map_cmd = "bwa mem" + self.index_loc + " " + self.freads + " " + self.rreads + " | samtools sort -o " + self.bam + " -"
        if not os.path.exists(self.bam):
            subprocess.call(map_cmd, shell = True) # Only way I could figure out how to implement with the pipe to samtools
    
    def count_coverage(self, contig_dict):
        depth = pysam.depth(self.bam)
        lines_for_df = [line.split('\t') for line in depth.split('\n')]
        cov_df = pd.DataFrame(lines_for_df, columns = ['Contig', 'Position', 'Depth'], dtype = float).dropna()
        contig_cov_df = cov_df.groupby(['Contig'])['Depth'].mean().reset_index(name = 'Mean Coverage')
        cov_dict = dict(zip(contig_cov_df['Contig'], contig_cov_df['Mean Coverage']))
        
        for seqid, seqcontents in contig_dict.items():
            for contig, cov in cov_dict.items():
                if seqid == contig:
                    seqcontents["cov"] = cov
        
        return contig_dict

    # --------- Single copy genes --------      
    
    def predict_cds(self):
        self.cds = os.path.join(self.tmp_dir, self.mag_name + ".faa")
        prodigal_cmd = ["prodigal", "-i", self.mag, "-a", self.cds, "-p", "meta", "-q"]
        if not os.path.exists(self.cds):
            subprocess.run(prodigal_cmd)
            
    def run_hmmsearch(self):
        self.hmm = os.path.join(self.tmp_dir, self.mag_name + ".tab")
        e_value = "1e-15"
        scg_db = SCG_DB
        hmmsearch_cmd = ["hmmsearch", "-E", e_value, "--tblout", self.hmm, scg_db, self.cds]
        if not os.path.exists(self.hmm):
            subprocess.run(hmmsearch_cmd)
    
    def parse_hmmtbl(self, hmmtbl):
        lines = []
        with open(hmmtbl, 'r') as tbl:
            for line in tbl:
                if line.startswith("#"):
                    continue
                else:
                    line = line.split()
                    lines.append(line)
        return lines
                
    def create_scg_dict(self, contig_dict):
        hmmtbl = self.parse_hmmtbl(self.hmm)
        scg_dict = defaultdict(dict)
        for hit in hmmtbl:
            contig = hit[0].split("_")[0]
            hmm_acc = hit[2]
            e_val = hit[4]
            bitscore = float(hit[5])
            start_pos = int(hit[19])
            end_pos = int(hit[21])
            hit_len = abs(start_pos - end_pos)
            if hit[23] == "-1":
                strand = "R"
            else:
                strand = "F"
            scg_dict[contig][hmm_acc] = {"e_val": e_val,
                                         "bitscore": bitscore,
                                         "start_pos": start_pos,
                                         "end_pos": end_pos,
                                         "hit_len": hit_len,
                                         "strand": strand}
        
        for seqid, seqcontent in contig_dict.items():
            for contig, hmm in scg_dict.items():
                if seqid == contig:
                    seqcontent["scg"] = hmm
        
        return contig_dict
            
    # --------- Taxonomy --------      
            
    def create_mag_db(self):
        self.db = os.path.join(self.tmp_dir, self.mag_name + "_db")
        create_db_cmd = ["mmseqs", "createdb", self.cds, self.db]
        subprocess.run(create_db_cmd)
    
    def run_taxonomy(self):
        self.tax = os.path.join(self.tmp_dir, self.mag_name + "_taxdb")
        self.taxtsv = os.path.join(self.tmp_dir, self.mag_name + "_tax.tsv")
        gettax_cmd = ["mmseqs", "taxonomy", self.db, TAX_DB, self.tax, self.tmp_dir, "--merge-query", "1", "--remove-tmp-files", "--tax-lineage", "1"]
        tsv_cmd = ["mmseqs", "createtsv", self.db, self.tax, self.taxtsv]
        subprocess.run(gettax_cmd)
        subprocess.run(tsv_cmd)
        
    def create_tax_dict(self, contig_dict):
        tax = pd.read_csv(self.taxtsv, sep = '\t', header=None, names=['Contig','Acc','Cat','LCA','Full Tax'])
        tax = tax.join(tax['Full Tax'].str.split(';', expand = True)).drop(['Acc', 'Cat', 'Full Tax'], axis = 1)
        tax[['Contig Name','ORF']] = tax['Contig'].str.rsplit('_', 1, expand=True)
        tax = tax.drop(tax.iloc[:, 10:36], axis = 1)
        tax = tax[tax['LCA'] != 'unclassified']
        contig_names = tax["Contig Name"].unique().tolist()
        contig_taxonomy = defaultdict()
        for contig in contig_names:
            contig_phylist = []
            contig_classlist = []
            contig_orderlist = []
            contig_famlist = []
            contig_genuslist = []
            contig_specieslist = []
            for index, row in tax.iterrows():
                if row["Contig Name"] == contig:
                    rowlist = [item for item in row.tolist() if item is not None]
                    for item in rowlist:
                        if item.startswith("p_"):
                            contig_phylist.append(item)
                        elif item.startswith("c_"):
                            contig_classlist.append(item)
                        elif item.startswith("o_"):
                            contig_orderlist.append(item)
                        elif item.startswith("f_"):
                            contig_famlist.append(item)
                        elif item.startswith("g_"):
                            contig_genuslist.append(item)
                        elif item.startswith("s_"):
                            contig_specieslist.append(item)
            contig_taxonomy[contig] = {"phylum": self.count_tax_hits(contig_phylist),
                                  "class": self.count_tax_hits(contig_classlist),
                                  "order": self.count_tax_hits(contig_orderlist),
                                  "family": self.count_tax_hits(contig_famlist),
                                  "genus": self.count_tax_hits(contig_genuslist),
                                  "species": self.count_tax_hits(contig_specieslist)}
        
        for seqid, seqcontents in contig_dict.items():
            for contig, taxonomy in contig_taxonomy.items():
                if seqid == contig:
                    seqcontents["tax"] = taxonomy
        
        return contig_dict
                               
    def count_tax_hits(self, tax_list):
        contig_tax = Counter(tax_list).most_common(1)[0][0]
        if contig_tax is None:
            return "Unknown"
        else:
            return contig_tax
            
    # --------- Mobilome --------      
                
    def search_trep(self):
        self.trep = os.path.join(self.tmp_dir, self.mag_name + "_trepdb")
        self.trep_tsv = os.path.join(self.tmp_dir, self.mag_name + "_trep.tsv")
        search_cmd = ["mmseqs", "search", self.db, TREP_DB, self.trep, self.tmp_dir, "--search-type", "3"]
        convert_cmd = ["mmseqs", "convertalis", self.db, TREP_DB, self.trep, self.trep_tsv]
        subprocess.run(search_cmd)
        subprocess.run(convert_cmd)
        
    def parse_mmseqs_search(self, mmseqs_tsv):
        lines = []
        with open(mmseqs_tsv, "r") as stsv:
            for line in stsv:
                line = line.split()
                lines.append(line)
        return lines
    
    def create_mobilome_dict(self, contig_dict):
        trep_res = self.parse_mmseqs_search(self.trep_tsv)
        mobilome_dict = defaultdict(dict)
        for hit in trep_res:
            contig = hit[0].split("_")[0]
            trep = hit[1]
            e_value = hit[11]
            bitscore = float(hit[12])
            start_pos = int(hit[7])
            end_pos = int(hit[8])
            hit_len = abs(start_pos - end_pos)
            mobilome_dict[contig][trep] = {"e_val": e_value,
                                         "bitscore": bitscore,
                                         "start_pos": start_pos,
                                         "end_pos": end_pos,
                                         "hit_len": hit_len}
        
        for seqid, seqcontent in contig_dict.items():
            for contig, trep in mobilome_dict.items():
                if seqid == contig:
                    seqcontent['mob'] = trep
                    
        return contig_dict

    # --------- PCA --------      
    
    def create_pca_dict(self, contig_dict):
        pca = decomposition.PCA(n_components = 1)
        tetramers = []
        for seqcontents in contig_dict.values():
            tetramer_list = seqcontents["tetra_proportion"]
            tetramers.append(tetramer_list)
        tetramer_df = pd.concat(tetramers)
        pca.fit(tetramer_df)
        first_axis = pca.components_[0]
        contig_names = list(contig_dict.keys())
        pca_dict = dict(zip(contig_names, first_axis))
        
        for seqid, seqcontents in contig_dict.items():
            for contig, pc in pca_dict.items():
                if seqid == contig:
                    seqcontents["pca"] = pc
                    
        return contig_dict

    # ~~~~~~~~~~~ Run MAG-wide MAG statistics ~~~~~~~~~~~~
    
    def craft_summary_dict(self, contig_dict):
        summary_dict = defaultdict()
        summary_dict["total_length"] = sum(x["len"] for x in contig_dict.values())
        gc = self.find_mean_std(contig_dict, "gc")
        summary_dict["gc_mean"] = gc[0]
        summary_dict["gc_std"] = gc[1]
        cov = self.find_mean_std(contig_dict, "cov")
        summary_dict["cov_mean"] = cov[0]
        summary_dict["cov_std"] = cov[1]
        pc = self.find_mean_std(contig_dict, "pca")
        summary_dict["pca_mean"] = pc[0]
        summary_dict["pca_std"] = pc[1]
        comp_red = self.find_completeness_redundancy(contig_dict)
        summary_dict["completeness"] = comp_red[0]
        summary_dict["redundancy"] = comp_red[1]
        tax = self.find_consensus_taxonomy(contig_dict)
        summary_dict["phylum"] = tax[0]
        summary_dict["class"] = tax[1]
        
        self.summary_dict = summary_dict
    
    # Generic function to query MAG contig dict
    def find_mean_std(self, contig_dict, query_for_dict):
        query_list = []
        length_list = []
        for seqcontents in contig_dict.values():
            query_val = seqcontents[query_for_dict]
            query_list.append(query_val)
            length = seqcontents["len"]
            length_list.append(length)
        mean_val = np.average(query_list, weights = length_list)
        std_val = np.std(query_list)
        return mean_val, std_val
    
    # SCG completeness and redundancy
    def find_completeness_redundancy(self, contig_dict):
        total_scgs = 111
        scg_dict = {}
        for seqcontents in contig_dict.values():
            contig_scg_dict = seqcontents["scg"]
            for hmm in contig_scg_dict.keys():
                if hmm not in scg_dict.keys():
                    scg_dict[hmm] = 1
                else:
                    scg_dict[hmm] += 1
                            
        duplicate_hits = {hmm: hits for hmm, hits in scg_dict.items() if hits > 1}
        
        completeness = len(scg_dict) / total_scgs
        redundancy = len(duplicate_hits) / total_scgs
        
        return completeness, redundancy
    
    # Consensus taxonomy at the phylum and class ranks
    def find_consensus_taxonomy_phylum(self, contig_dict):
        phylum_dict = defaultdict(int)
        class_dict = defaultdict(int)
        for seqcontents in contig_dict.values():
            phylum_dict[seqcontents["tax"]["phylum"]] += 1
            class_dict[seqcontents["tax"]["class"]] += 1
        most_common_phylum = max(phylum_dict, key = phylum_dict.get)
        most_common_class = max(class_dict, key = class_dict.get)
        return most_common_phylum, most_common_class
