import yaml

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