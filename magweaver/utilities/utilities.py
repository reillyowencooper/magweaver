#!/usr/bin/env python3

import os, gzip, tarfile, shutil, logging

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