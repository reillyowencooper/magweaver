import os
import subprocess

def get_taxdb(data_dir, taxdb_dir, tmp_dir):
    "Get MMSeqs prebuilt SwissProt database for taxonomy and index"
    taxdb_loc = os.path.join(data_dir, taxdb_dir)
    get_db_cmd = ["mmseqs", "databases", "UniProtKB/Swiss-Prot", taxdb_loc, tmp_dir]
    subprocess.run(get_db_cmd)
    create_taxdb_cmd = ["mmseqs", "createtaxdb", taxdb_loc, tmp_dir]
    subprocess.run(create_taxdb_cmd)
    index_db_cmd = ["mmseqs", "createindex", taxdb_loc, tmp_dir]
    subprocess.run(index_db_cmd)
    
def unzip_scgs(data_dir):
    "Unzip the SCG HMMs that come built into package"
    scg_zip_loc = os.path.join(data_dir, "essential.hmm.gz")
    gunzip_cmd = ["gunzip", scg_zip_loc]
    subprocess.run(gunzip_cmd)
    
def get_trepdb(data_dir, trep_dir, tmp_dir):
    "Download the TREP database, convert to MMSeqs database, and place"
    trepdb_loc = os.path.join(data_dir, trep_dir)
    tmp_loc = os.path.join(tmp_dir, "trep.gz")
    get_db_cmd = ["wget", "-O", tmp_loc, "https://botserv2.uzh.ch/kelldata/trep-db/downloads/trep-db_complete_Rel-19.fasta.gz"]
    subprocess.run(get_db_cmd)
    create_trepdb_cmd = ["mmseqs", "createdb", tmp_loc, trepdb_loc]
    subprocess.run(create_trepdb_cmd)
    
# TODO: Incorporate custom HMM database input
