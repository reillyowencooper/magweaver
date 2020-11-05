import os, subprocess, gzip, shutil, logging
import src.file_handling as filehandling
import src.utilities as utils

class DataBaseHandler(object):
    
    def __init__(self):
        self.database_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
        self.dbs = {}
        self.tmp_dir = os.path.join(self.database_dir, "tmp")
        
    def retrieve_swissprot(self):
        logging.info('Downloading MMSeqs2 SwissProt database')
        uniref_mmseqs_loc = os.path.join(self.database_dir, 'swissprot')
        utils.create_dir(uniref_mmseqs_loc)
        swissprot_loc = os.path.join(uniref_mmseqs_loc, 'swissprot')
        getdb_cmd = ['mmseqs', 'databases', 'UniProtKB/Swiss-Prot', swissprot_loc, self.tmp_dir]
        subprocess.run(getdb_cmd)
        gettaxdb_cmd = ['mmseqs', 'createtaxdb', swissprot_loc, self.tmp_dir]
        subprocess.run(gettaxdb_cmd)
        index_cmd = ['mmseqs', 'createindex', swissprot_loc, self.tmp_dir]
        subprocess.run(index_cmd)
        self.dbs['swissprot'] = swissprot_loc
        
    def unzip_scg(self):
        logging.info('Unzipping single-copy gene HMMs')
        scg_loc = os.path.join(self.database_dir, 'essential.hmm.gz')
        gunzip_cmd = ['gunzip', scg_loc]
        subprocess.run(gunzip_cmd)
        self.dbs['scg'] = scg_loc
        
    def run(self):
        self.retrieve_swissprot()
        self.unzip_scg()