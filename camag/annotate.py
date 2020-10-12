import os
import subprocess
from Bio import SeqIO

# Sequentially search genome against databases in the following order:
# 1: Kofam
# 2: AMRfinder
# 3: VOGdb
# 4: Pfam
# 5: Uniref50

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