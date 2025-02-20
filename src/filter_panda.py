# Import libraries
import argparse
import pandas as pd
from eland import filter_panda, process_panda

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process and filter panda network.")
    subparsers = parser.add_subparsers(dest='func', help='sub-command help')

    # Create a parser for 'process_edge_list' function
    parser_process_edge_list = subparsers.add_parser('process_edge_list', help='Process edge list')
    parser_process_edge_list.add_argument('input_file', type=str, help='Path to edge list file.')
    parser_process_edge_list.add_argument('output_file', type=str, help='Path to output file.')

    # Create a parser for 'filter_panda' function
    parser_filter_panda = subparsers.add_parser('filter_panda', help='Filter panda network')
    parser_filter_panda.add_argument('prior_file', type=str, help='Path to motif prior edge list file.')
    parser_filter_panda.add_argument('delimiter', type=str, help='Delimiter used in the edge list files')
    parser_filter_panda.add_argument('edgelist_file', type=str, help='Path to output file.')
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    # process panda result into edgelist
    process_panda.process_edge_list(args.input_file, args.output_file)
    
    # filter panda network
    fil_edges = filter_panda.filter_panda(args.prior_file, args.output_file, args.delimiter)

    # save edgelist to file
    fil_edges.to_csv(args.edgelist_file, sep=args.delimiter, index=False)

if __name__ == "__main__":
    main()
