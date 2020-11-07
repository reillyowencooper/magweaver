import os, subprocess, shutil, logging, sys
from collections import defaultdict, Counter
import pandas as pd
from Bio import SeqIO
from src.cov import MagCov
from src.gc import MagGC
from src.tax import MagTaxonomy
from src.tetranuc import MagTetranuc
from src.scg import MagSCG
import src.utilities as utils
import src.file_handling as filehandling

class Decontaminator(object):
    def __init__(self, mag, 
                 forward_reads, 
                 reverse_reads, 
                 swissprot_db, 
                 scg_db,
                 outdir):
        self.mag = mag
        self.forward_reads = forward_reads
        self.reverse_reads = reverse_reads
        self.swissprot_db = swissprot_db
        self.scg_db = scg_db
        self.mag_name = os.path.splitext(os.path.basename(self.mag))[0]
        self.outdir = outdir
        self.tmp_dir = os.path.join(self.outdir, "tmp")

    def get_mag_gc(self):
        logging.info('Flagging contigs with erroneous GC content in ' + self.mag_name)
        gc = MagGC(self.mag, self.outdir)
        gc.run()
        gc_df = pd.read_csv(os.path.join(self.outdir, self.mag_name + "_err_gc.csv"))
        gc_contigs = gc_df['Contig'].tolist()
        self.gc = gc
        return gc_contigs
    
    def get_mag_cov(self):
        logging.info('Flagging contigs with erroneous read coverage content in ' + self.mag_name)
        cov = MagCov(self.mag, self.forward_reads, self.reverse_reads, self.outdir, self.tmp_dir)
        cov.run()
        cov_df = pd.read_csv(os.path.join(self.outdir, self.mag_name + "_err_cov.csv"))
        cov_contigs = cov_df['Contig'].tolist()
        self.cov = cov
        return cov_contigs
    
    def get_mag_tax(self):
        logging.info('Flagging contigs with erroneous taxonomy in ' + self.mag_name)
        tax = MagTaxonomy(self.mag, self.outdir, self.tmp_dir)
        tax.run(self.swissprot_db, self.tmp_dir)
        tax_df = pd.read_csv(os.path.join(self.outdir, self.mag_name + "_err_tax.csv"))
        tax_contigs = tax_df['Contig'].tolist()
        self.tax = tax
        return tax_contigs
    
    def get_mag_tetra(self):
        logging.info('Flagging contigs with erroneous tetranucleotide frequencies in ' + self.mag_name)
        tetra = MagTetranuc(self.mag, self.outdir)
        tetra.run()
        tetra_df = pd.read_csv(os.path.join(self.outdir, self.mag_name + "_err_tetra.csv"))
        tetra_contigs = tetra_df['Contig'].tolist()
        self.tetra = tetra
        return tetra_contigs
    
    def rank_suspicion(self, gc_contigs, cov_contigs, tax_contigs, tetra_contigs):
        logging.info('Ranking contigs by suspicious content for ' + self.mag_name)
        all_lists = gc_contigs + cov_contigs + tax_contigs + tetra_contigs
        sus_dict = dict(Counter(all_lists))
        sus_df = pd.DataFrame(list(sus_dict.items()), columns = ['Contig', 'Suspicion Level'])
        return sus_dict, sus_df
    
    def find_high_sus_contigs(self, sus_dict, minimum_suspicion):
        logging.info('Flagging contigs with suspicion greater than or equal to ' + str(minimum_suspicion) + ' in ' + self.mag_name)
        high_sus_contigs = []
        for contig_name, sus_rank in sus_dict.items():
            if sus_rank >= minimum_suspicion:
                high_sus_contigs.append(contig_name)
        return high_sus_contigs
    
    def filter_mag(self, high_sus_contigs):
        logging.info('Filtering ' + self.mag_name + ' of contigs over minimum suspicion')
        mag_contigs = filehandling.ReadNucFasta(self.mag).retrieve_contigs()
        filtered_contigs = {contig_name: contig_seq for contig_name, contig_seq in mag_contigs.items() if not contig_name in high_sus_contigs}
        return filtered_contigs
    
    def write_fasta(self, contig_dict):
        outfile = os.path.join(self.outdir, self.mag_name + "_decontaminated.fa")
        logging.info('Writing filtered version of ' + self.mag_name + ' to ' + outfile)
        for contig_name, contig_seq in contig_dict.items():
            with open(outfile, 'a+') as outmag:
                outmag.write('>' + str(contig_name) + '\n' + str(contig_seq) + '\n')
    
    def check_scgs(self, contig_dict):
        pass # Do something here with the single-copy genes
        
    def compare_after_decontam(self, filtered_contig_dict):
        start_mag = filehandling.ReadNucFasta(self.mag).retrieve_contigs()
        filtered_mag = filtered_contig_dict
        num_removed_contigs = len(start_mag) - len(filtered_mag)
        seqlength_start_mag = 0
        seqlength_filtered_mag = 0
        for contig_seq in start_mag.values():
            seqlength_start_mag += len(contig_seq)
        for contig_seq in filtered_mag.values():
            seqlength_filtered_mag += len(contig_seq)
        seqlength_removed = seqlength_start_mag - seqlength_filtered_mag
        return seqlength_start_mag, num_removed_contigs, seqlength_removed
        
    def run(self, minimum_suspicion_limit):
        utils.create_dir(self.outdir)
        utils.create_dir(self.tmp_dir)
        gc = self.get_mag_gc()
        cov = self.get_mag_cov()
        tax = self.get_mag_tax()
        tetra = self.get_mag_tetra()
        suspicious_contigs, suspicious_df = self.rank_suspicion(gc, cov, tax, tetra)
        very_suspicious = self.find_high_sus_contigs(suspicious_contigs, minimum_suspicion_limit)
        filtered_contigs = self.filter_mag(very_suspicious)
        self.write_fasta(filtered_contigs)
        suspicious_df.to_csv(os.path.join(self.outdir, self.mag_name + "_contig_suspicion.csv"), index = False, header = True)
        stats = self.compare_after_decontam(filtered_contigs)
        print('After filtering ' + self.mag_name + ', ' + str(stats[1]) + ' contigs were removed, removing ' + str(stats[2]) + ' from an initial length of ' + str(stats[0]))
        utils.remove_tmp_dir(self.tmp_dir)        

def decontaminate_wrapper(swissprot_db, scg_db, **kwargs):
    decontam = Decontaminator(kwargs['input_mag'],
                              kwargs['forward_reads'],
                              kwargs['reverse_reads'],
                              swissprot_db,
                              scg_db,
                              kwargs['output_dir'])
    decontam.run(kwargs['suspicion'])