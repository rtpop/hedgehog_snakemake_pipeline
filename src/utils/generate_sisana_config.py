import yaml
import argparse

def generate_sisana_params(exp_file, motif_file, ppi_file, number, outdir, processed_paths, method, pandafilepath, output_filename):
    data = {
        'preprocess': {
            'exp_file': exp_file,
            'motif_file': motif_file,
            'ppi_file': ppi_file,
            'number': number,
            'outdir': outdir
        },
        'generate': {
            'processed_paths': processed_paths,
            'method': method,
            'pandafilepath': pandafilepath
        }
    }

    with open(output_filename, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SiSaNA params YAML file")
    parser.add_argument('--exp', required=True, help='Path to expression file')
    parser.add_argument('--motif', required=True, help='Path to motif file')
    parser.add_argument('--ppi', required=True, help='Path to PPI file')
    parser.add_argument('--number', type=int, required=True, help='Number of samples a gene must be expressed in')
    parser.add_argument('--outdir', required=True, help='Output directory for preprocess stage')
    parser.add_argument('--processed_paths', required=True, help='Path to processed data paths YAML file')
    parser.add_argument('--method', required=True, help='Method to use (panda or lioness)')
    parser.add_argument('--pandafilepath', required=True, help='Path to PANDA output file')
    parser.add_argument('--output', required=True, help='Output YAML file name')

    args = parser.parse_args()

    generate_sisana_params(
        exp_file=args.exp,
        motif_file=args.motif,
        ppi_file=args.ppi,
        number=args.number,
        outdir=args.outdir,
        processed_paths=args.processed_paths,
        method=args.method,
        pandafilepath=args.pandafilepath,
        output_filename=args.output
    )