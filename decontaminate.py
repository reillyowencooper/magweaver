import argparse
from time import sleep
from src.mag_aligner import MagAligner

# TODO: Once remove_contam is finished, add to script

def parse_args():
    parser = argparse.ArgumentParser(description = 'Workflow for examining contigs in a bin for potential contamination')
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('--mag', help = 'Path to input MAG')
    required.add_argument('--forward', help = 'Path to forward reads')
    required.add_argument('--reverse', help = 'Path to reverse reads')
    optional.add_argument('--keep_bam', help = 'Specify whether to keep sorted output BAM')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    mag_aligned = MagAligner(args.mag, args.forward, args.reverse, args.keep_bam)
    mag_aligned.make_tmp_dir()
    mag_aligned.index_mag()
    mag_aligned.align_reads()
    mag_aligned.sam_to_bam()
    mag_aligned.sort_bam()
    mag_aligned.write_files()
    num = mag_aligned.get_num_contigs()
    gc = mag_aligned.get_contig_gc()
    contig_depth = mag_aligned.get_contig_coverage_depth()
    tetras = mag_aligned.get_contig_tetranucleotide_freq()
    full = mag_aligned.build_complete_stats_df(tetras, contig_depth, gc)
    
if __name__ == "__main__":
    main()