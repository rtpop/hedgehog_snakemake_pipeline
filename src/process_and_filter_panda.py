# Import libraries
import argparse
import pandas as pd
import eland

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process and filter panda network.")
    subparsers = parser.add_subparsers(dest='func', help='sub-command help')

    # Create a parser for 'process_edge_list' function
    parser_process_edge_list = subparsers.add_parser('process_edge_list', help='Process edge list')
    parser_process_edge_list.add_argument('input_file', type=str, help='Path to edge list file.')
    parser_process_edge_list.add_argument('output_file', type=str, help='Path to output file.')
    parser_process_edge_list.add_argument('delimiter', type=str, help='Delimiter used in the edge list files')

    # Create a parser for 'filter_panda' function
    parser_filter_panda = subparsers.add_parser('filter_panda', help='Filter panda network')
    parser_filter_panda.add_argument('prior_file', type=str, help='Path to motif prior edge list file.')
    parser_process_edge_list.add_argument('edgelist', type=str, help='Path to output file.')
    parser_filter_panda.add_argument('delimiter', type=str, help='Delimiter used in the edge list files')
    
    # create subparser for to_csv
    parser_to_csv = subparsers.add_parser('to_csv', help='print to csv')
    parser_to_csv.add_argument('filtered_net', type=str, help='Path to output file.')
    parser_to_csv.add_argument('delimiter', type=str, help='Delimiter used in the edge list files')
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    # process panda result into edgelist
    eland.process_edge_list(args.input_file, args.output_file, sep = args.delimiter)
    
    # filter panda network
    fil_edges = eland.filter_panda.filter_panda(args.prior_file, args.edgelist, args.delimiter)

    # save edgelist to file
    fil_edges.to_csv(args.filtered_net, sep=",", index=False)

if __name__ == "__main__":
    main()
