import os, subprocess, shutil, logging
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from src.cov import MagCov
from src.gc import MagGC
from src.tax import MagTaxonomy
from src.tetranuc import MagTetranuc
import src.utilities as utils

# TODO: Weave together coverage, GC, taxonomy, and tetranucleotide

def create_tmp_dir(tmpdir):
    utils.create_dir("tmp")

def check_tax_db():
    pass # TODO: Write func to check Swissprot MMSeqs is present

def rank_suspicion(gc_list, cov_list, tax_list, tetra_list):
    '''Weights suspicion of contigs based on how many lists they appear in'''
    # This is probably the dumbest way to rank contigs, 
    # but I'm tired and can't think of anything better
    sus_dict = defaultdict(int)
    checked_items = []
    for item in gc_list:
        sus_dict[item] += 1
        if item in cov_list:
            sus_dict[item] += 1
        elif item in tax_list:
            sus_dict[item] += 1
        elif item in tetra_list:
            sus_dict[item] += 1
        checked_items.append(item)
    updated_cov_list = [x for x in cov_list if x not in checked_items]
    for item in updated_cov_list:
        sus_dict[item] += 1
        if item in tax_list:
            sus_dict[item] += 1
        elif item in tetra_list:
            sus_dict[item] += 1
        checked_items.append(item)
    updated_tax_list = [x for x in tax_list if x not in checked_items]
    for item in updated_tax_list:
        sus_dict[item] += 1
        if item in tetra_list:
            sus_dict[item] += 1
        checked_items.append(item)
    updated_tetra_list = [x for x in tetra_list if x not in checked_items]
    for item in updated_tetra_list:
        sus_dict[item] += 1
    return sus_dict

def get_consensus_contigs(mag, bam, taxdb,
                          tmpdir, outdir):
    mag_name = os.path.splitext(os.path.basename(mag))[0]
    aa_file = mag_name + ".faa"
    utils.predict_cds(mag, aa_file)
    gc = MagGC(mag, outdir)
    cov = MagCov(mag, bam, outdir)
    tax = MagTaxonomy(aa_file, outdir)
    tetra = MagTetranuc(mag, outdir)
    gc.run()
    cov.run()
    tax.run(taxdb, tmpdir)
    tetra.run()
    gc_df = pd.read_csv(mag_name + "_err_gc.csv")
    gc_contigs = gc_df['Contig'].tolist()
    cov_df = pd.read_csv(mag_name + "_err_cov.csv")
    cov_contigs = cov_df['Contig'].tolist()
    tax_df = pd.read_csv(mag_name + "_err_tax.csv")
    tax_contigs = tax_df['Contig'].tolist()
    tetra_df = pd.read_csv(mag_name + "err_tetra.csv")
    tetra_contigs = tetra_df['Contig'].tolist()
    ranked_contigs = rank_suspicion(gc_contigs, cov_contigs,
                                    tax_contigs, tetra_contigs)
    return ranked_contigs
    
def remove_sus_contigs(mag, ranked_contig_dict, outdir, minimum_suspicion=3):
    outfile = os.path.join(outdir, os.path.basename(mag))
    for contig_name, sus_rank in ranked_contig_dict.items():
        if sus_rank >= minimum_suspicion:
            for seqrecord in SeqIO.parse(mag, "fasta"):
                seq_name = seqrecord.id
                if not seq_name == contig_name:
                    with open(outfile, 'w') as newmag:
                        newmag.write('>' + str(seq_name) + '\n' + str(seqrecord.seq) + '\n')
 
def main():
    pass