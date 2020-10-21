import os, subprocess, shutil, logging
import pandas as pd
from src.cov import MagCov
from src.gc import MagGC
from src.tax import MagTaxonomy
from src.tetranuc import MagTetranuc
import src.utilities as utils

# TODO: Weave together coverage, GC, taxonomy, and tetranucleotide

def create_tmp_dir(tmpdir):
    utils.create_dir("tmp")

def check_tax_db():
    pass

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
    cov_df = pd.read_csv(mag_name + "_err_cov.csv")
    tax_df = pd.read_csv(mag_name + "_err_tax.csv")
    tetra_df = pd.read_csv(mag_name + "err_tetra.csv")
    # Get dataframe of contgs that are erroneous in all 4 categories, 3, 2, 1
    # return that df
    
def remove_sus_contigs(mag, outdir):
    pass # Write MAG fasta to file, excluding sus contigs
 
def main():
    pass