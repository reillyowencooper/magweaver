import argparse
from src import tax

def cli():
    parser = argparse.ArgumentParser(description = 'Create CMAGs from binned MAGs',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    tax_parser = subparsers.add_parser("contig-tax",
                                       help = "Identify contig-level taxonomy within a MAG",
                                       )
    tax.fetch_args(tax_parser)
    args = vars(parser.parse_args())
    args["func"](args)