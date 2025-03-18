from netZooPy import sambar
import argparse
import os

def arg_parser():
    parser = argparse.ArgumentParser(description='Run sambar')
    parser.add_argument('--mut-file', type=str, help='File with mutation data')
    parser.add_argument('--esize-file', type=str, help='File with gene exon size')
    parser.add_argument('--can-genes', type=str, help='File with cancer associated genes')
    parser.add_argument('--gmt-file', type=str, help='GMT file to be used for desparsification.')
    parser.add_argument('--output-dir', type=str, help='Output directory')
    return parser.parse_args()

def main():
    args = arg_parser()
    print(args)
    os.chdir(args.output_dir)
    
    # change input filenamesto account for changed directory
    args.mut_file = os.path.join("../../../../../" + args.mut_file)
    args.esize_file = os.path.join("../../../../../" + args.esize_file)
    args.can_genes = os.path.join("../../../../../" + args.can_genes)
    args.gmt_file = os.path.join("../../../../../" + args.gmt_file)
    
    sambar.sambar(args.mut_file, args.esize_file, args.can_genes, args.gmt_file)
    
if __name__ == '__main__':
    main()
