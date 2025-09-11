# Import libraries
import argparse
import pandas as pd
import hedgehog

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process and filter panda network.")

    # Parse args for the main function
    parser.add_argument('prior_file', type=str, help='Path to edge list file.')
    parser.add_argument('input_file', type=str, help='Path to edge list file.')
    parser.add_argument('output_file', type=str, help='Path to output file.')
    parser.add_argument('--delimiter', type=str, help='Delimiter used in the edge list files')
    parser.add_argument('--filtering_method', type=str, help='Filtering method to use. Options: prior, hedgehog, both')

    return parser.parse_args()

def main():
    args = parse_arguments()
    print("Filtering PANDA network using method:", args.filtering_method)
    if args.filtering_method not in ['prior', 'hedgehog', 'both', 'none']:
        raise ValueError("Invalid filtering method. Choose from: prior, hedgehog, both")
    
    if args.filtering_method == 'prior':
        # process panda result into edgelist
        hedgehog.filter_panda.filter_panda(args.prior_file, args.input_file, args.output_file, delimiter = args.delimiter, prior_only = True)
    elif args.filtering_method == 'hedgehog':
        hedgehog.filter_panda.filter_panda(args.prior_file, args.input_file, args.output_file, delimiter = args.delimiter, prior_only = False)
    elif args.filtering_method == 'both':
        # process panda result into edgelist
        hedgehog.filter_panda.filter_panda(args.prior_file, args.input_file, args.output_file.replace(".txt", "_prior.txt"), delimiter = args.delimiter, prior_only = True)
        hedgehog.filter_panda.filter_panda(args.prior_file, args.input_file, args.output_file.replace(".txt", "_hedgehog.txt"), delimiter = args.delimiter, prior_only = False)

if __name__ == "__main__":
    main()
