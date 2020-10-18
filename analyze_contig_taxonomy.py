#!/usr/bin/env python3

import argparse
from src.tax import MagTaxonomy

def parse_args():
    parser = argparse.ArgumentParser(description='Get contig-level taxonomy for a MAG')
    required = parser.add_argument_group("Required arguments")
    required.add_argument('--input_mag', help = "Path to predicted CDS file for MAG")
    required.add_argument('--search_db', default = "databases/swissprot/swissprot", help = "Path to MMSeqs SwissProt DB")
    required.add_argument('--output_dir', help = "Directory to place output files")
    required.add_argument('--tmp_dir', default= "tmp", help = "Directory to store temporary files")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    if not os.path.exists(args["tmp_dir"]):
        os.mkdir(args["tmp_dir"])
    magtax = MagTaxonomy(args["output_dir"])
    magtax.run(args["input_mag"],
               args["search_db"],
               args["tmp_dir"])
    
if __name__ == "__main__":
    main()