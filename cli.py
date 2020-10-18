#!/usr/bin/env python3

import argparse
import sys
from src import tax

def cli():
    parser = argparse.ArgumentParser(description = 'Create CMAGs from binned MAGs',
                                     help = 'Create CMAGs from binned MAGs',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument()
    subparsers = parser.add_subparsers()
    tax_parser = subparsers.add_parser("contig-tax",
                                        help="Identify contig taxonomy within a MAG",
                                        description="Identify contig taxonomy within a MAG",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        )
    tax.fetch_args(tax_parser)
    args = vars(parser.parse_args())
    args["func"](args)