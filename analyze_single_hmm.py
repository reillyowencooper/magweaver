import argparse
from src.single_hmm_searcher import SearchSingleHMM

def parse_args():
    parser = argparse.ArgumentParser(description = "Search a directory of MAGs against a single profile HMM, construct an abundance table and a phylogenetic tree")
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--mag_dir', help = 'Directory containing your input MAGs')
    required.add_argument('--output_dir', help = 'Directory to place output files and folders')
    required.add_argument('--hmm', help = 'Path to your profile HMM to search against')
    required.add_argument('--msa_opt', help = 'Choose MUSCLE or Clustal Omega for multiple sequence alignment')
    required.add_argument('--phylo_opt', help = 'Choose RAxML or FastTree for phylogenetic tree construction')
    required.add_argument('--e_value', help = 'The E-value cutoff for your HMM', default=1e-10)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    mag_examiner = SearchSingleHMM(args.mag_dir,
                                   args.output_dir,
                                   args.hmm,
                                   args.msa_opt,
                                   args.phylo_opt,
                                   args.e_value)
    mag_examiner.workhorse()


if __name__ == "__main__":
    main()