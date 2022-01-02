import os
import subprocess
from magweaver.database_handling import get_taxdb, unzip_scgs, get_trepdb

def prepare_databases(db_dict, tmp_dir):
    "Cycles through database dictionary and runs downloader or formatter if path does not exist"
    if not os.path.exists(db_dict["taxonomy"]):
        get_taxdb(db_dict["data"], db_dict["taxonomy"], tmp_dir)
    if not os.path.exists(db_dict["scg"]):
        unzip_scgs(db_dict["data"])
    if not os.path.exists(db_dict["trep"]):
        get_trepdb(db_dict["data"], db_dict["trep"], tmp_dir)
    else:
        print("All databases ready to use")