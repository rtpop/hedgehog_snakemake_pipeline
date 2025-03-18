import argparse
import os
import bihidef
import sys

def parse_args():
    parser = argparse.ArgumentParser("Run bihidef")
    parser.add_argument("edgelist_file", type = str, help="Input file")
    parser.add_argument("--max_res", type = float, help="The maximum resolution to use.")
    parser.add_argument("--comm_mult", type = float, help="The maximum number of communities to find.")
    parser.add_argument("--output_dir", type = str, help="Output directory.")
    parser.add_argument("--output_prefix_reg", type = str, help="Output prefix for regulator communities.")
    parser.add_argument("--output_prefix_tar", type = str, help="Output prefix for target communities.")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # change to output directory
    os.chdir(args.output_dir)
    
    # edit input file path to be relative to output directory
    args.edgelist_file = os.path.join("../../../../../", args.edgelist_file)

    # run bihidef
    bihidef.bihidef(filename = args.edgelist_file, maxres = args.max_res, comm_mult = args.comm_mult, oR= args.output_prefix_reg, oT = args.output_prefix_tar)

if __name__ == "__main__":
    main()