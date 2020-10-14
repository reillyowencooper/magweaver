import os, subprocess, shutil
import pandas as pd

class CompletenessChecker(object):
    def __init__(self, mag_dir, output_dir, completeness_val):
        self.mag_dir = os.path.abspath(mag_dir)
        self.output_dir = os.path.abspath(output_dir)
        self.completeness_val = float(completeness_val)
        self.contam_val = float(5.0)
        
    def create_output_dir(self):
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
            
    def run_checkm(self):
        checkm_cmd = ['checkm', 'lineage_wf', '-x', 'fa', '--tab_table', self.mag_dir, self.output_dir]
        subprocess.run(checkm_cmd)
        
    def parse_checkm(self):
        # Parse checkm for completeness and contamination, store MAG name, completeness/contamination in dict
        pass
    
    def place_highqual_bins(self, mag_dict):
        # Put bins > completeness, < contamination into output_dir for further work
        pass
    
    def summarize_completeness(self):
        # Summarize CheckM results in file
        pass
        