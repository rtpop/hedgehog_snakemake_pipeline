# Import libraries
import argparse
import pandas as pd
import eland

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process and filter panda network.")

    # Parse args for the main function
    parser.add_argument('prior_file', type=str, help='Path to edge list file.')
    parser.add_argument('input_file', type=str, help='Path to edge list file.')
    parser.add_argument('output_file', type=str, help='Path to output file.')
    parser.add_argument('--delimiter', type=str, help='Delimiter used in the edge list files')
    parser.add_argument('--prior_only', type=bool, help='Whether or not to filter only to the prior or perform eland filtering.')
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    # process panda result into edgelist
    eland.filter_panda.filter_panda(args.prior_file, args.input_file, args.output_file, delimiter = args.delimiter, prior_only = args.prior_only)

if __name__ == "__main__":
    main()
