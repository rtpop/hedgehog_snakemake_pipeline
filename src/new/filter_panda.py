# Import libraries
print("Importing libraries", flush=True)
import pandas as pd
import hedgehog


prior_file = snakemake.input.prior
input_file = snakemake.input.panda
output_file = snakemake.output.filtered_panda
delimiter = snakemake.params.delimiter
fil = snakemake.params.filtering_method

def main(prior_file, input_file, output_file, delimiter, filtering_method):
    print("Filtering PANDA network using method:", filtering_method, flush=True)
    if filtering_method not in ['prior', 'hedgehog', 'both', 'none']:
        raise ValueError("Invalid filtering method. Choose from: prior, hedgehog, both")

    if filtering_method == 'prior':
        # process panda result into edgelist
        hedgehog.filter_panda.filter_panda(prior_file, input_file, output_file[0], delimiter=delimiter, prior_only=True)
    elif filtering_method == 'hedgehog':
        hedgehog.filter_panda.filter_panda(prior_file, input_file, output_file[0], delimiter=delimiter, prior_only=False)
    elif filtering_method == 'both':
        if len(output_file) != 2:
            raise ValueError("For filtering_method 'both', provide two output files: prior and hedgehog.")
        # process panda result into edgelist
        hedgehog.filter_panda.filter_panda(prior_file, input_file, output_file[0], delimiter=delimiter, prior_only=True)
        hedgehog.filter_panda.filter_panda(prior_file, input_file, output_file[1], delimiter=delimiter, prior_only=False)
    elif filtering_method == 'none':
        # do nothing
        pass

if __name__ == "__main__":
    main(prior_file, input_file, output_file, delimiter, fil)