from netZooPy import sambar
import argparse

def arg_parser():
    parser = argparse.ArgumentParser(description='Run sambar')
    parser.add_argument('--mut_file', type=str, help='File with mutation data')
    parser.add_argument('--esize_file', type=str, help='File with gene exon size')
    parser.add_argument('--can_genes', type=int, help='File with cancer associated genes')
    parser.add_argument('--gmt-file', type=str, help='GMT file to be used for desparsification.')
    parser.add_argument('--output_dir', type=str, help='Output directory')
    return parser.parse_args()

def main():
    args = arg_parser()
    sambar.sambar(args.mut_file, args.esize_file, args.can_genes, args.gmt_file)