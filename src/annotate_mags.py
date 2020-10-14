import os, subprocess, shutil, logging, csv
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
    
    def combine_hmmtbls(self, list_of_hmmtbls, output_name):
        lines_for_output = []
        for hmmtbl in list_of_hmmtbls:
            with open(hmmtbl, 'r') as tblout:
                for line in tblout:
                    if line.startswith('#'):
                        continue
                    else:
                        line = line.split()
                        line_fixed = line[:18] + [' '.join(line[18:])]
                        lines_for_output.append(line_fixed)
        with open(output_name, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(lines_for_output)
                        
    def search_hmmdbs(self):
        self.bac_unannot = []
        self.arc_unannot = []
        self.vir_unannot = []
        self.bac_hmmannot = []
        self.arc_hmmannot = []
        self.vir_hmmannot = []
        kofam_tmpfile = os.path.join(self.tmp_dir, "tmp_seqfile_kofam.faa")
        pfam_tmpfile = os.path.join(self.tmp_dir, "tmp_seqfile_pfam.faa")
        for mag in self.bac_cds:
            hmm_files = []
            files_to_remove = []
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            output_path = os.path.join(self.tmp_dir, mag_name + '_kofam.tab')
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["kofam"], mag, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            files_to_remove.append(kofam_tmpfile)
            self.reduce_input_fasta(hit_contigs, mag, kofam_tmpfile)
            output_path = os.path.join(self.tmp_dir, mag_name + '_pfam.tab')
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["pfam"], kofam_tmpfile, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            files_to_remove.append(pfam_tmpfile)
            self.reduce_input_fasta(hit_contigs, kofam_tmpfile, pfam_tmpfile)
            output_path = os.path.join(self.tmp_dir, mag_name + '_amr.tab')
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["amrfinder"], pfam_tmpfile, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            unannotated_seqfile = os.path.join(self.tmp_dir, mag_name + "_unannotated.faa")
            self.bac_unannot.append(unannotated_seqfile)
            self.reduce_input_fasta(hit_contigs, pfam_tmpfile, unannotated_seqfile)
            self.combine_hmmtbls(hmm_files, os.path.join(tmp_dir, mag_name + '_hmms.tab'))
            self.bac_hmmannot.append(os.path.join(tmp_dir, mag_name + '_hmms.tab'))
            for filename in files_to_remove:
                os.remove(filename)
        for mag in self.arc_cds:
            hmm_files = []
            files_to_remove = []
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            output_path = os.path.join(self.tmp_dir)
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["kofam"], mag, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            files_to_remove.append(kofam_tmpfile)
            self.reduce_input_fasta(hit_contigs, mag, kofam_tmpfile)
            output_path = os.path.join(self.tmp_dir, mag_name + '_pfam.tab')
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["pfam"], kofam_tmpfile, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            files_to_remove.append(pfam_tmpfile)
            self.reduce_input_fasta(hit_contigs, kofam_tmpfile, pfam_tmpfile)
            output_path = os.path.join(self.tmp_dir, mag_name + '_amr.tab')
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["amrfinder"], pfam_tmpfile, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            unannotated_seqfile = os.path.join(self.tmp_dir, mag_name + "_unannotated.faa")
            self.arc_unannot.append(unannotated_seqfile)
            self.reduce_input_fasta(hit_contigs, pfam_tmpfile, unannotated_seqfile)
            self.combine_hmmtbls(hmm_files, os.path.join(tmp_dir, mag_name + '_hmms.tab'))
            self.arc_hmmannot.append(os.path.join(tmp_dir, mag_name + '_hmms.tab'))
            for filename in files_to_remove:
                os.remove(filename)
        for mag in self.vir_cds:
            hmm_files = []
            files_to_remove = []
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            output_path = os.path.join(self.tmp_dir)
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["kofam"], mag, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            files_to_remove.append(kofam_tmpfile)
            self.reduce_input_fasta(hit_contigs, mag, kofam_tmpfile)
            output_path = os.path.join(self.tmp_dir, mag_name + '_pfam.tab')
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["pfam"], kofam_tmpfile, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            files_to_remove.append(pfam_tmpfile)
            self.reduce_input_fasta(hit_contigs, kofam_tmpfile, pfam_tmpfile)
            output_path = os.path.join(self.tmp_dir, mag_name + '_vog.tab')
            hmm_files.append(output_path)
            files_to_remove.append(output_path)
            self.run_hmmsearch(self.hmm_dbs["vog"], pfam_tmpfile, output_path)
            hit_contigs = self.parse_hmmtbl(output_path)
            unannotated_seqfile = os.path.join(self.tmp_dir, mag_name + "_unannotated.faa")
            self.vir_unannot.append(unannotated_seqfile)
            self.reduce_input_fasta(hit_contigs, pfam_tmpfile, unannotated_seqfile)
            self.combine_hmmtbls(hmm_files, os.path.join(tmp_dir, mag_name + '_hmms.tab'))
            self.vir_hmmannot.append(os.path.join(tmp_dir, mag_name + '_hmms.tab'))
            for filename in files_to_remove:
                os.remove(filename)          
            # TODO: Figure out how eukaryotes work (life tip, just in general)
    
    def run_mmseqs(self, input_aa, list_to_append_output):
        mag_name = os.path.splitext(os.path.basename(input_aa))[0]
        output_path = os.path.join(self.tmp_dir, mag_name + '.db')
        create_mmseqs_db_cmd = ['mmseqs', 'createdb', input_aa, output_path]
        list_to_append_output.append(output_path)
        subprocess.run(create_mmseqs_db_cmd)
        
    def convert_aa_to_mmseqs(self):
        self.bac_mmseqs = []
        self.arc_mmseqs = []
        self.vir_mmseqs = []
        for mag in self.bac_unannot:
            self.run_mmseqs(mag, self.bac_mmseqs)
        for mag in self.arc_unannot:
            self.run_mmseqs(mag, self.arc_mmseqs)
        for mag in self.vir_mmseqs:
            self.run_mmseqs(mag, self.vir_mmseqs)
            
    def run_mmseq_search(self, input_db):
        # Run MMseqs on each unannotated MAG
        pass
    
    def workhorse(self):
        self.create_logger()
        self.prepare_output_dirs()
        self.get_mags()
        
        
    
        