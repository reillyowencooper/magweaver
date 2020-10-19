import os, subprocess, shutil, logging
import pandas as pd

# TODO: Move count_tetramers, get_contig_tetranucleotide_freq here
# TODO: Use PCA to get clusters, then flag contigs that don't cluster