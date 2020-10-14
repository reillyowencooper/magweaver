import os, subprocess, shutil, logging
import pandas as pd
from Bio import SearchIO, SeqIO


class MagAnnotator(object):
    
    def __init__(self, mag_dir, database_dir, e_value, output_dir):
        self.mag_dir = os.path.abspath(mag_dir)
        self.mag_type_dirs = {"bacteria": os.path.join(self.mag_dir, "bacteria_bins"),
                              "archaea": os.path.join(self.mag_dir, "archaea_bins"),
                              "eukaryota": os.path.join(self.mag_dir, "eukaryota_bins"),
                              "virus": os.path.join(self.mag_dir, "virus_bins")}
        self.database_dir = os.path.abspath(database_dir)
        self.hmm_dbs = {"kofam": os.path.join(self.database_dir, "kofam.hmm"),
                    "amrfinder": os.path.join(self.database_dir, "amrfinder.hmm"),
                    "vog": os.path.join(self.database_dir, "vog.hmm"),
                    "pfam": os.path.join(self.database_dir, "pfam.hmm")}
        self.mmseqs_dbs = {"uniref": os.path.join(self.database_dir, "uniref_mmseqs.db")}
        self.e_value = e_value
        self.output_dir = output_dir
        self.mag_list = []
        
    def create_logger(self):
        self.logger = logging.getLogger(__file__)
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('annotate_mags.log')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
          
    def prepare_output_dirs(self):
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        self.tmp_dir = os.path.join(self.output_dir, "tmp")
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
            
    def get_mags(self):
        self.bac_list = []
        self.arc_list = []
        self.euk_list = []
        self.vir_list = []
        for mag_type, loc in self.mag_type_dirs.items():
            for mag in os.listdir(loc):
                if mag_type == "bacteria":
                    self.bac_list.append(os.path.join(loc, mag))
                elif mag_type == "archaea":
                    self.arc_list.append(os.path.join(loc, mag))
                elif mag_type == "eukaryota":
                    self.euk_list.append(os.path.join(loc, mag))
                elif mag_type == "virus":
                    self.vir_list.append(os.path.join(loc, mag))
                    
    def predict_cds(self):
        # Doing it this way in case viruses/archaea have a better ORF predictor than Prodigal
        self.bac_cds = []
        self.arc_cds = []
        self.euk_cds = []
        self.vir_cds = []
        for mag in self.bac_list:
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            self.logger.info('Predicting coding regions for ' + mag_name + ' using Prodigal')
            output_path = os.path.join(self.tmp_dir, mag_name + ".faa")
            prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta', '-q']
            subprocess.run(prodigal_cmd)
            self.bac_cds.append(output_path)
        for mag in self.arc_list:
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            self.logger.info('Predicting coding regions for ' + mag_name + ' using Prodigal')
            output_path = os.path.join(self.tmp_dir, mag_name + ".faa")
            prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta', '-q']
            subprocess.run(prodigal_cmd)
            self.arc_cds.append(output_path)
        for mag in self.vir_list:
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            self.logger.info('Predicting coding regions for ' + mag_name + ' using Prodigal')
            output_path = os.path.join(self.tmp_dir, mag_name + ".faa")
            prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta', '-q']
            subprocess.run(prodigal_cmd)
            self.vir_cds.append(output_path)
        for mag in self.euk_list:
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            self.logger.info('Predicting coding regions for ' + mag_name + ' using MetaEuk, this may take a while')
            output_path = os.path.join(self.tmp_dir, mag_name + ".faa")
            metaeuk_cmd = ['metaeuk', 'easy-predict', mag, self.dbs['uniref'], mag_name, self.output_dir]            subprocess.run(prodigal_cmd)
            self.euk_cds.append(output_path)
            
    def run_hmmsearch(self, hmm_db, input_cd, output_path):
        hmmsearch_cmd = ['hmmsearch', '-E', str(self.e_value),
                         '--tblout', output_path, hmm_db, input_cd]
        subprocess.run(hmmsearch_cmd)
    
    def parse_hmmtbl(self, hmmtbl):
        contig_ids = []
        with open(hmmtbl, 'r') as tblout:
            lines_for_df = list()
            for line in tblout:
                if line.startswith('#'):
                    continue
                else:
                    line = line.split()
                    line_fixed = line[:18] + [' '.join(line[18:])]
                    lines_for_df.append(line_fixed)
        for line in lines_for_df:
            contig_name = str(line[0])
            contig_ids.append(contig_name)
        return contig_ids
    
    def reduce_input_fasta(self, contig_ids, input_cd, output_fasta):
        with open(output_fasta, "w") as reduced:
            with SeqIO.parse(open(input_cd), "fasta") as aas:
                for record in aas:
                    if record.id not in contig_ids:
                        SeqIO.write([record], reduced, "fasta")
    
    def search_hmmdbs(self):
        for mag in self.bac_cds:
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            output_path = os.path.join(self.tmp_dir, mag_name + '_kofam.tab')
            self.run_hmmsearch(self.hmm_dbs["kofam"], mag, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            tmp_seqfile_kofam = os.path.join(self.tmp_dir, "tmp_seqfile_kofam.faa")
            self.reduce_input_fasta(hit_contigs, mag, tmp_seqfile_kofam)
            output_path = os.path.join(self.tmp_dir, mag_name + '_pfam.tab')
            self.run_hmmsearch(self.hmm_dbs["pfam"], tmp_seqfile_kofam, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            tmp_seqfile_pfam = os.path.join(self.tmp_dir, "tmp_seqfile_pfam.faa")
            self.reduce_input_fasta(hit_contigs, tmp_seqfile_kofam, tmp_seqfile_pfam)
            output_path = os.path.join(self.tmp_dir, mag_name + '_amr.tab')
            self.run_hmmsearch(self.hmm_dbs["amrfinder"], tmp_seqfile_pfam, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            unannotated_seqfile = os.path.join(self.tmp_dir, mag_name + "_unannotated.faa")
            self.reduce_input_fasta(hit_contigs, tmp_seqfile_pfam, unannotated_seqfile)
            os.remove(tmp_seqfile_kofam)
            os.remove(tmp_seqfile_pfam)
            # TODO: Make dict of all HMM hits together to assign to sequence headers
            # TODO: Now do this for archaea, viruses (need to use VOG)
            
    
    
    def workhorse(self):
        self.create_logger()
        self.prepare_output_dirs()
        self.get_mags()
        
        
    
        