import os, subprocess, gzip, shutil, logging
from Bio import SeqIO
import src.file_handling as filehandling
import src.utilities as utils

## TODO: COMPLETELY REVAMP DATABASE HANDLING

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
        self.dbs['swissprot'] = swissprot_loc
        if not os.path.exists(swissprot_loc):
            getdb_cmd = ['mmseqs', 'databases', 'UniProtKB/Swiss-Prot', swissprot_loc, self.tmp_dir]
            subprocess.run(getdb_cmd)
            gettaxdb_cmd = ['mmseqs', 'createtaxdb', swissprot_loc, self.tmp_dir]
            subprocess.run(gettaxdb_cmd)
            index_cmd = ['mmseqs', 'createindex', swissprot_loc, self.tmp_dir]
            subprocess.run(index_cmd)
          
    def retrieve_trep(self):
        logging.info('Downloading eukaryotic repetitive sequences from TREP')
        self.dbs['trep_fasta'] = os.path.join(self.database_dir, "trep.fasta")
        # Keep bacterial transposable elements!
        bac_treps = ('>DTX_Bumb', '>DTX_Bcre', '>DTX_Ecol', '>DTX_Sphi') 
        if not os.path.exists(self.dbs['trep_fasta']):
            trep_gz_loc = os.path.join(self.database_dir, 'trep_all.fasta.gz')
            retrieve_cmd = ['wget', '-O', trep_gz_loc, 'https://botserv2.uzh.ch/kelldata/trep-db/downloads/trep-db_complete_Rel-19.fasta.gz']
            subprocess.run(retrieve_cmd)
            gunzip_cmd = ['gunzip', trep_gz_loc]
            subprocess.run(gunzip_cmd)
        for seqrecord in SeqIO.parse(os.path.join(self.database_dir, 'trep_all.fasta'), 'fasta'):
            with open(self.dbs['trep_fasta'], 'a+') as trep:
                if not seqrecord.id.startswith(bac_treps):
                    trep.write('>' + str(seqrecord.id) + '\n' + str(seqrecord.seq) + '\n')
       
    def unzip_scg(self):
        logging.info('Unzipping single-copy gene HMMs')
        self.dbs['scg'] = os.path.join(self.database_dir, 'essential.hmm')
        if not os.path.exists(self.dbs['scg']):
            scg_loc = os.path.join(self.database_dir, 'essential.hmm.gz')
            gunzip_cmd = ['gunzip', scg_loc]
            subprocess.run(gunzip_cmd)
           
    def create_trep_db(self):
        logging.info('Converting TREP database to MMSeqs database')
        
        trep_mmseqs_loc = os.path.join(self.database_dir, 'trep')
        self.dbs['trep_db'] = os.path.join(trep_mmseqs_loc, 'trep')
        if not os.path.exists(self.dbs['trep_db']):
            utils.create_dir(trep_mmseqs_loc)
            createdb_cmd = ['mmseqs', 'createdb', self.dbs['trep_fasta'], os.path.join(trep_mmseqs_loc, 'trep')]
            subprocess.run(createdb_cmd)
        
    def run(self):
        self.retrieve_swissprot()
        self.retrieve_trep()
        self.unzip_scg()
        self.create_trep_db()
        return self.dbs