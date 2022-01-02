import os, subprocess, shutil

def download_file(file_url, output_loc = None, verbose = True):
    '''Downloads a file, prints to stdout if no output_loc specified'''
    url_str = str(file_url)
    if verbose:
        print(f'Downloading file from {url_str}')
    if output_loc is not None:
        run_external()
        
def run_external(external_cmd, shell = False, verbose = True, keep_stdout = False):
    '''Basic external shell command function, can run via shell or not'''
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if not keep_stdout:
        subprocess.run(external_cmd, shell = shell, 
                       stdout = stdout, stderr = stderr)
    else:
        return subprocess.run(external_cmd, shell = shell, 
                              stdout = subprocess.PIPE, 
                              stderr = stderr).stdout.decode(errors = "ignore")
    
def create_dir(filepath):
    '''Checks if a directory exists then creates it if not'''
    if not os.path.exists(filepath):
        os.mkdir(filepath)
        
def run_hmmsearch(aa_file, e_value, output_path, hmm):
    '''Runs HMMsearch, given an input file, an HMM, whatever E-value you want to use'''
    hmmsearch_cmd = ['hmmsearch', '-E', e_value, '--tblout', output_path, hmm, aa_file]
    if not os.path.exists(output_path):
        subprocess.run(hmmsearch_cmd)
        
def predict_cds(mag, output_path):
    '''Runs Prodigal on an input MAG'''
    prodigal_cmd = ['prodigal', '-i', mag, '-a', output_path, '-p', 'meta', '-q']
    if not os.path.exists(output_path):
        subprocess.run(prodigal_cmd)
        
def parse_hmmtbl(hmm_tblout):
    lines = []
    with open(hmm_tblout, 'r') as tbl:
        for line in tbl:
            if line.startswith('#'):
                continue
            else:
                line = line.split()
                lines.append(line)
    return lines

def remove_tmp_dir(tmp_dir):
    shutil.rmtree(tmp_dir)

logger = logging.getLogger(__name__)

def gunzip_file(filepath, outpath = None):
    """Gunzips a .gz file

    Args:
        filepath (str): Path to gzip'd file
        outpath (str, optional): Path to output file. Defaults to None, which writes to filename of gzip'd file without the .gz.

    Returns:
        [str]: Path to output file
    """
    if filepath.endswith('.gz'):
        if outpath is None:
            fname = os.path.splitext(os.path.basename(filepath))[0]
        else:
            fname = fname
        if os.path.exists(fname):
            raise FileExistsError(fname)
        try:
            with gzip.open(filepath, 'rb') as inf:
                with open(fname, 'wb') as outf:
                    shutil.copyfileobj(inf, outf)
                    logger.debug('Gunzipped ' + filepath + ' to ' + fname)
                    os.remove(filepath)
            return fname
        except IOError:
            logger.debug('Gunzipping ' + filepath + ' failed, likely due to corrupt input file')