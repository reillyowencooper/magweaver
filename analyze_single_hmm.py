import argparse
from src.single_hmm_searcher import SearchSingleHMM

def parse_args():
    parser = argparse.ArgumentParser(description = "Search a directory of MAGs against a single profile HMM, construct an abundance table and a phylogenetic tree")
    parser.add_argument('--mag_dir', help = 'Directory containing your input MAGs')
    parser.add_argument('--output_dir', help = 'Directory to place output files and folders')
    parser.add_argument('--hmm', help = 'Path to your profile HMM to search against')
    parser.add_argument('--msa_opt', help = 'Choose MUSCLE or Clustal Omega for multiple sequence alignment')
    parser.add_argument('--phylo_opt', help = 'Choose RAxML or FastTree for phylogenetic tree construction')
    parser.add_argument('--e_value', help = 'The E-value cutoff for your HMM', default=1e-10)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    mag_examiner = SearchSingleHMM()
    mag_examiner.workhorse()


if __name__ == "__main__":
    main()