import os, subprocess, logging, shutil
import pandas as pd
from Bio import SeqIO, SearchIO
import src.camag_utilities as utils

# Sequentially search genome against databases in the following order:
# 1: Kofam
# 2: AMRfinder
# 3: VOGdb
# 4: Pfam
# 5: Uniref50

# NOTE: I made this early on, and it doesn't work. Remaking in annotate_mags.py so it actually works

class Annotator(object):
    
    def __init__(self, mag_dir, output_dir, database_dir, e_value, mag_type):
        self.mag_dir = os.path.abspath(mag_dir)
        self.output_dir = os.path.abspath(output_dir)
        self.database_dir = os.path.abspath(database_dir)
        self.e_value = e_value
        self.mag_type = mag_type
        self.hmm_dbs = {"kofam": os.path.join(self.database_dir, "kofam.hmm"),
                    "amrfinder": os.path.join(self.database_dir, "amrfinder.hmm"),
                    "vog": os.path.join(self.database_dir, "vog.hmm"),
                    "pfam": os.path.join(self.database_dir, "pfam.hmm")
                    }
        self.mmseqs_dbs = {"uniref": os.path.join(self.database_dir, "uniref_mmseqs.db")}
        self.mag_list = []
        
    def create_logger(self):
        self.logger = logging.getLogger(__file__)
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('single_hmm_search.log')
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        
    def get_mags(self):
        for mag in os.listdir(self.mag_dir):
            mag_name = os.path.basename(mag)
            if mag.endswith('.fa') or mag.endswith('.fna'):
                mag_path = os.path.join(self.mag_dir, mag)
                self.mag_list.append(mag_path)
                self.logger.info('Preparing ' + mag_name + ' for ' + self.mag_type + ' annotation')
            else:
                self.logger.info('Omitting ' + mag_name + ' due to incorrect filetype')
                
    def create_output_folders(self):
        self.bac_output = os.path.join(self.output_dir, "bacteria_bin_annotation")
        self.arc_output = os.path.join(self.output_dir, "archaea_bin_annotation")
        self.euk_output = os.path.join(self.output_dir, "eukaryota_bin_annotation")
        self.vir_output = os.path.join(self.output_dir, "virus_bin_annotation")
        for annotype in [self.bac_output, self.arc_output, self.euk_output, self.vir_output]:
            utils.create_output_dir(annotype)
            
    def predict_cds(self):
        for mag in self.mag_list:
            mag_name = os.path.splitext(os.path.basename(mag))[0]
            if self.mag_type == "Eukaryota":
                self.logger.info('Predicting coding regions for ' + mag_name + ' using MetaEuk')
                output_path = os.path.join(self.euk_output, mag_name + '.faa')
                metaeuk_cmd = ['metaeuk', 'easy-predict', mag, self.dbs['uniref'], mag_name, self.output_dir]
                if not os.path.exists(output_path):
                    subprocess.run(metaeuk_cmd)
                else:
                    self.logger.info('Prediction file already exists, moving on')
                self.euk_cds[mag_name] = output_path
            else:
                self.logger.info('Predicting coding regions for ' + mag_name + ' using Prodigal')
                if self.mag_type == "Bacteria":
                    output_path = os.path.join(self.bac_output, mag_name + '.faa')
                    self.bac_cds[mag_name] = output_path
                elif self.mag_type == "Archaea":
                    output_path = os.path.join(self.arc_output, mag_name + '.faa')
                    self.arc_cds[mag_name] = output_path
                elif self.mag_type == "Virus":
                    output_path = os.path.join(self.vir_output, mag_name + '.faa')
                    self.vir_cds[mag_name] = output_path
                prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta', '-q']
                if not os.path.exists(output_path):
                    subprocess.run(prodigal_cmd)
                else:
                    self.logger.info('Prediction file already exists, moving on')
                
    def search_hmms(self):
        # TODO: Write function that combines the following three in a nice way
        # Basically, create a way to search each MAG against each HMM database sequentially,
        # reducing search space each time by taking out CDS that have already hit to something in an earlier HMM database
        # Bacteria first
        for mag_name, cd_loc in self.bac_cds.items():
            
    def parse_hmmoutput(self, hmm_tbl):
        hit_ids = {}
        name = os.path.splitext(os.path.basename(hmm_tbl))[0]
        with open(hmm_tbl, 'r') as tblout:
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
            hit_name = str(line[2])
            hit_evalue = float(line[4])
            if hit_evalue <= float(self.e_value):
                hit_ids[name] = contig_name
        return hit_ids
    
    def reduce_search(self, mag_path, hit_dict, output_path):
        mag_name = os.path.splitext(os.path.basename(mag_path))
        for cd_file in self.cds_list:
            cd_name = os.path.splitext(os.path.basename(cd_file))[0]
            hit_dict_magonly = {k:v for k,v in hit_dict.items() if k == cd_name}
            for mag in self.mag_list:
                mag_name = os.path.splitext(os.path.basename(mag))[0]
                if mag_name == cd_name:
                    with open(output_path, "w") as filtered:
                        with SeqIO.parse(open(cd_file), "fasta") as aa_file:
                            for record in aa_file:
                                if record.id not in hit_dict_magonly.values():
                                    SeqIO.write([record], filtered, "fasta")
                                    
    def search_hmms(self, hmm_db_name, hmm_db_loc):
        # TODO: Rethink this and see if there's a less clunky way to write
        for mag_name, cd_path in self.bac_cds:
            # Run against KoFam first
            db_loc = self.hmm_dbs["kofam"]
            output_path = os.path.join(self.bac_output, mag_name + "_kofam.tab")
            hmmsearch_cmd = ['hmmsearch', '-E', str(self.e_value), '--tblout', output_path, db_loc, cd_path]
            subprocess.run(hmmsearch_cmd)
            hits = self.parse_hmmoutput(output_path)
            reduced_path = os.path.join(self.bac_output, mag_name + "_reduced.fa")
            self.reduce_search(hits, reduced_path)
            # Run reduced FASTA against Pfam
            db_loc = self.hmm_dbs["pfam"]
            output_path = os.path.join(self.bac_output, mag_name + "_pfam.tab")
            hmmsearch2_cmd = ['hmmsearch', '-E', str(self.e_value), '--tblout', output_path, db_loc, reduced_path]
            subprocess.run(hmmsearch2_cmd)
            hits = self.parse_hmmoutput(output_path)
            os.remove(reduced_path)
            self.reduce_search(hits, reduced_path)
            # Run against AMRfinder
            db_loc = self.hmm_dbs["amrfinder"]
            output_path = os.path.join(self.bac_output, mag_name + "_amr.tab")
            subprocess.run(hmmsearch2_cmd)
            hits = self.parse_hmmoutput(output_path)
            self.reduce_search()
            
            for db_name, db_loc in self.hmm_dbs:
                self.logger.info('Searching ' + mag_name + ' against ' + db_name)
                
                
                if not os.path.exists(output_path):
                    subprocess.run(hmmsearch_cmd)

    

class MagAnnotator(object):

    def __init__(self, input_mag, output_dir, database_dir, evalue, bitscore):
        self.mag = os.path.abspath(input_mag)
        self.mag_name = os.path.splitext(input_mag)[0]
        self.output_dir = os.path.abspath(output_dir)
        self.database_dir = os.path.abspath(database_dir)
        self.evalue = str(evalue)
        self.bitscore = bitscore
        self.dbs = {"kofam": os.path.join(self.database_dir, "kofam.hmm"),
                    "amrfinder": os.path.join(self.database_dir, "amrfinder.hmm"),
                    "vog": os.path.join(self.database_dir, "vog.hmm"),
                    "pfam": os.path.join(self.database_dir, "pfam.hmm"),
                    "uniref": os.path.join(self.database_dir, "uniref_mmseqs.db")}

    def prodigal_search(self):
        print('Finding ORFs')
        aa_name = os.path.join(self.output_dir, self.mag_name + '.faa')
        prodigal_cmd = ['prodigal', '-i', self.mag, '-a', aa_name, '-p', 'meta']
        subprocess.run(prodigal_cmd)
        self.aa_filepath = aa_name

    def hmm_search(self, db_name, db_loc, aa_file):
        print('Searching HMMs of ' + db_name)
        domtblout_filepath = os.path.join(self.output_dir, self.mag_name + '_' + db_name + '.tab')
        hmmsearch_cmd = ['hmmsearch', '-E', self.evalue, '--domtblout', domtblout_filepath, '--noali', db_loc, aa_file]
        subprocess.run(hmmsearch_cmd)
        return domtblout_filepath

    def domtblout_parser(self, domtblout_loc):
        print("Parsing HMMsearch output")
        with open(os.path.join(self.output_dir, domtblout_loc), 'r') as domtblout:
            lines_for_df = list()
            for line in domtblout:
                if line.startswith('#'):
                    continue
                else:
                    line = line.split()
                    line_fixed = line[:22] + [' '.join(line[22:])]
                    lines_for_df.append(line_fixed)
        hit_dict = {}
        for line in lines_for_df:
            contig = str(line[0])
            hit_ko = str(line[3])
            hit_bitscore = float(line[7])
            if hit_bitscore >= self.bitscore:
                hit_dict[contig] = hit_ko, hit_bitscore
        return hit_dict

    def reduce_search(self, hit_dict):
        print("Reducing input fasta based on prior HMMsearch results")
        already_hit_contigs = hit_dict.keys()
        filtered_mag_loc = os.path.join(self.output_dir, os.path.splitext(self.aa_filepath)[0] + '_filtered.faa')
        aa_file = SeqIO.parse(open(self.aa_filepath), "fasta")
        with open(filtered_mag_loc, "w") as filtered:
            for record in aa_file:
                if record.id not in already_hit_contigs:
                    SeqIO.write([record], filtered, "fasta")

    def convert_aa_to_mmseqs(self, input_fasta):
        print("Converting remaining ORFs to MMSeqs2 db")
        output_db = os.path.join(os.output_dir, os.path.splitext(self.mag)[0] + '_mmseqs.db')
        createdb_cmd = ['mmseqs', 'createdb', input_fasta, output_db]
        subprocess.run(createdb_cmd)
        return output_db

    def search_uniref(self, input_db):
        print("Searching UniRef50 db")
        input_db_name = os.path.splitext(input_db)[0]
        output_db_name = os.path.join(self.output_dir, input_db_name + '_uniref.db')
        output_besthit_db_name = os.path.join(self.output_dir, input_db_name + '_uniref_besthit.db')
        mmseqs_search_cmd = ['mmseqs', 'search', input_db, self.dbs["uniref"], output_db_name, self.output_dir]
        subprocess.run(mmseqs_search_cmd)
        get_besthit_cmd = ['mmseqs', 'filterdb', output_db_name, output_besthit_db_name, '--extract-lines', '1']
        subprocess.run(get_besthit_cmd)
        tab_output_name = os.path.join(self.output_dir, input_db_name + '_uniref.tab')
        create_output_cmd = ['mmseqs', 'convertalis', input_db, self.dbs["uniref"], output_besthit_db_name, tab_output_name]
        subprocess.run(create_output_cmd)
        hit_dict = {}
        with open(tab_output_name, 'r') as hits:
            for line in hits:
                contig = str(line[0])
                identifier = str(line[1])
                bitscore = float(line[11])
            if bitscore >= self.bitscore:
                hit_dict[contig] = identifier, bitscore
        return hit_dict

    def complete_search(self):
        self.prodigal_search()
        kofam_search = self.hmm_search("kofam", self.dbs["kofam"], self.aa_filepath)
        total_hits = {}
        kofam_hitdict = self.domtblout_parser(kofam_search)
        total_hits.update(kofam_hitdict)
        self.reduce_search(total_hits)
        filtered_fasta = os.path.splitext(self.aa_filepath)[0] + "_filtered.faa"
        amrfinder_search = self.hmm_search("amrfinder", self.dbs["amrfinder"], filtered_fasta)
        amrfinder_hitdict = self.domtblout_parser(amrfinder_search)
        total_hits.update(amrfinder_hitdict)
        os.remove(filtered_fasta)
        self.reduce_search(total_hits)
        vog_search = self.hmm_search("vog", self.dbs["kofam"], filtered_fasta)
        vog_hitdict = self.domtblout_parser(vog_search)
        total_hits.update(vog_hitdict)
        os.remove(filtered_fasta)
        self.reduce_search(total_hits)
        pfam_search = self.hmm_search("pfam", self.dbs["pfam"], filtered_fasta)
        pfam_hitdict = self.domtblout_parser(pfam_search)
        total_hits.update(pfam_hitdict)
        os.remove(filtered_fasta)
        self.reduce_search(total_hits)
        remaining_aas = self.convert_aa_to_mmseqs(filtered_fasta)
        uniref_hitdict = self.search_uniref(remaining_aas)
        total_hits.update(uniref_hitdict)
        os.remove(filtered_fasta)
        self.reduce_search(total_hits)
        return total_hits