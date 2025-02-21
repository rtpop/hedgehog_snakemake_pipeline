import argparse
import os
from ELAND import bihidef

def parse_args():
    parser = argparse.ArgumentParser("Run bihidef")
    parser.add_argument("edgelist_file", required=True, type = str, help="Input file")
    parser.add_argument("max_res", type = float, help="The maximum resolution to use.")
    parser.add_argument("max_communities", type = int, help="The maximum number of communities to find.")
    parser.add_argument("output_dir", type = str, help="Output directory")
    parser.add_argument("output_prefix_reg", type = str, help="Output prefix for regulator communities.")
    parser.add_argument("output_prefix_tar", type = str, help="Output prefix for target communities.")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # run bihidef
    bihidef.bihidef(args.edgelist_file, args.max_res, args.max_communities, args.output_dir, args.output_prefix_reg, args.output_prefix_tar)