#!/usr/bin/env python3

import argparse
from src import tax

def cli():
    parser = argparse.ArgumentParser(description = 'Create CMAGs from binned MAGs',
                                     help = 'Create CMAGs from binned MAGs',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument()
