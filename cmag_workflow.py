import os, shutil, subprocess, logging, argparse
import pandas as pd
import numpy as np
from src.completeness import CompletenessChecker
from src.mag_aligner import MagAligner
from src.mag_taxonomy import MagTaxonomy
from src.database_handling import DatabaseDownloader

def parse_args():
    parser = argparse.ArgumentParser(description = 'Complete workflow for processing MAG bins into CMAGs')
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('--mag_dir', help = 'Directory with MAGs of interest')
    required.add_argument('--output_dir', help = 'Directory to output files')
    required.add_argument('--assembly', help = 'Path to metagenome assembly for MAGs of interest')
    required.add_argument('--forward', help = 'Path to all forward reads used in assembly')
    required.add_argument('--reverse', help = 'Path to all reverse reads used in assembly')
    optional.add_argument('--completeness', help = 'Minimum completeness value for MAG to qualify as high quality. Default = 95')
    optional.add_argument('--contamination', help = 'Maximum contamination value for MAG to qualify as high quality. Default = 5')
    optional.add_argument('--keep_bam', help = 'Keep sorted BAM for each paired read-MAG alignment, default FALSE')
    args = parser.parse_args()
    return args

def place_hiqual_bins(mag_dir, checkm_df, completeness, contamination, output_dir):
    hq_bins = []
    for index, row in checkm_df.iterrows():
        if row['Completeness'] >= completeness and row['Contamination'] <= contamination:
            bin_id = row['Bin Id']
            for mag in os.listdir(mag_dir):
                if os.path.basename(mag) == bin_id:
                    shutil.copyfile(os.path.join(mag_dir, mag), os.path.join(output_dir, mag))
    
def prep_nt(database_dir):
    if not os.path.exists(database_dir):
        os.mkdir(database_dir)
    databaseHandler = DatabaseDownloader(database_dir)
    databaseHandler.retrieve_swissprot_mmseqs()
    
def main():
    args = parse_args()
    prep_nt("databases")
    completor = CompletenessChecker(args.mag_dir,
                                    args.output_dir,
                                    args.completeness,
                                    args.contamination)
    completor.create_logger()
    completor.create_output_dir()
    completor.run_checkm()
    contig_count = completor.count_contigs()
    final_df = completor.summarize_completeness_contam(contig_count)
    place_hiqual_bins(args.mag_dir,
                      final_df,
                      args.completeness,
                      args.contamination,
                      args.output_dir)
    for filename in os.listdir(args.output_dir):
        aligned = MagAligner(os.path.join(args.output_dir, filename),
                             args.forward,
                             args.reverse,
                             args.keep_bam)
        aligned.make_tmp_dir()
        aligned.index_mag()
        aligned.align_reads()
        aligned.sam_to_bam()
        aligned.sort_bam()
        aligned.write_files()
        num = aligned.get_num_contigs()
        gc = aligned.get_contig_gc()
        depth = aligned.get_contig_coverage_depth()
        tetras = aligned.get_contig_tetranucleotide_freq()
        full = aligned.build_complete_stats_df(tetras, depth, gc)
        mean_depth = depth['Mean Coverage'].mean()
        sd_depth = depth['Mean Coverage'].std()
        mean_gc = gc['GC Content'].mean()
        sd_gc = gc['GC Content'].std()
        swissprot_db = os.path.join()
        mag_tax = MagTaxonomy(os.path.join(args.output_dir, filename),
                               args.output_dir,
                               "databases/nt/")
        mag_tax.create_outputs()
        mag_tax.create_mag_mmseqsdb()
        mag_tax.get_contig_taxonomy()
        taxonomy = mag_tax.parse_taxonomy()
        # TODO: Write function that defines what a weird contig looks like
        weird_contigs = []
        
        
    
