import click
import os
from decontaminate import decontaminator
from prepare import prepare_databases

# TODO: Put relevant filepaths into a JSON file
BASEPATH = os.path.dirname(os.path.abspath(__file__))
TMP_DIR = os.path.join(BASEPATH, "tmp")
OUT_DIR = os.path.join(BASEPATH, "results")
DB_DIR = os.path.join(BASEPATH, "data")

DB_PATHS = {"data": DB_DIR,
            "taxonomy": os.path.join(DB_DIR, "swissprot/swissprot"),
            "scg": os.path.join(DB_DIR, "essential.hmm"),
            "trep": os.path.join(DB_DIR, "trep/trep")}

@click.command()
def preparedb():
    prepare_databases(DB_PATHS, TMP_DIR)

@click.command()
@click.option('--mag', help = 'Path to the FASTA file of a MAG')
@click.option('--fwd', help = 'Path to forward reads used to create MAG')
@click.option('--rev', help = 'Path to reverse reads used to create MAG')
@click.option('--threads', help = 'Number of threads to use')
@click.option('--sus', default = 2, help = 'Minimum number of erroneous features to flag a contig')
@click.option('--out', help = 'Path to output filtered MAG FASTA')
@click.option('--tmp', default = TMP_DIR, help = 'Path to use as temporary directory')          
def decontam(mag, fwd, rev, threads, sus, out, tmp):
    decontaminator(mag, fwd, rev, threads, sus, out, tmp)




if __name__ == '__main__':
    decontam()