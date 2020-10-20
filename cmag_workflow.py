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

# TODO: Combine tax, gc, cov, tetranuc, quality here
        
        
    
