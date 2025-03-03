import argparse
from eland import process_bihidef

def parese_args():
    parser = argparse.ArgumentParser("Select communities")
    parser.add_argument("communities", type = str, help="Input file")
    parser.add_argument("output_gmt", type = str, help="File name for the gmt file to be saved.")
    parser.add_argument("--logs", type = str, help="File name for log file detailing the selected communities.")
    parser.add_argument("--min_size", type = int, help="Minimum size of the community to be selected.")
    parser.add_argument("--max_size", type = int, help="Maximum size of the community to be selected.")
    return parser.parse_args()

def main():
    args = parese_args()
    
    # Select communities
    communities = process_bihidef.select_communities(args.communities, args.min_size, args.max_size, args.logs)
    
    # Save selected communities to a gmt file    
    process_bihidef.gmt_from_bihidef(communities, args.output_gmt)
    
if __name__ == "__main__":
    main()
