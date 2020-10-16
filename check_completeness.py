import argparse
from src.completeness import CompletenessChecker

def parse_args():
    parser = argparse.ArgumentParser(description = "Quantify completion and contamination of MAGs using CheckM")
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('--mag_dir', help = 'Directory containing your input MAGs')
    required.add_argument('--output_dir', help = 'Directory to output summary')
    optional.add_argument('--completeness', help = 'Minimum completeness to qualify as High Quality, default 95')
    optional.add_argument('--contamination', help = 'Maximum contamination to qualify as High Quality, default 5')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    completor = CompletenessChecker(args.mag_dir,
                                       args.output_dir,
                                       args.completeness,
                                       args.contamination)
    completor.create_logger()
    completor.create_output_dir()
    completor.run_checkm()
    contig_count = completor.count_contigs()
    final_df = completor.summarize_completeness_contam(contig_count)
    completor.df_to_csv(final_df)

if __name__ == "__main__":
    main()