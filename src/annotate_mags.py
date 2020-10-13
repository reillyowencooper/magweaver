import os, subprocess, shutil, logging
import pandas as pd
from Bio import SearchIO, SeqIO


class MagAnnotator(object):
    
    def __init__(self):
        