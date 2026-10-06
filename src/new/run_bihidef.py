import os
import sys
from pathlib import Path
import bihidef

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "old" / "utils"))
from measure_resources import measure_resources_context

# diagnostics
from condor import condor_object

_original_qscores = condor_object.qscores


def diagnostic_qscores(self, *args, **kwargs):
    try:
        return _original_qscores(self, *args, **kwargs)
    except KeyError as exc:
        community = exc.args[0]

        # Only diagnose missing community contributions.
        if community in self.Qcol_lookup:
            raise

        n_reg = int((self.reg_memb["community"] == community).sum())
        n_tar = int((self.tar_memb["community"] == community).sum())

        raise RuntimeError(
            f"Qscore failed: community {community} is missing from "
            f"Qcol_lookup; regulators={n_reg}, targets={n_tar}"
        ) from exc


condor_object.qscores = diagnostic_qscores
# end diagnostics

net = snakemake.input.net
gene_communities = snakemake.output.gene_communities
max_communities = snakemake.params.max_communities
max_res = snakemake.params.max_resolution
output_prefix_reg = snakemake.params.output_prefix_reg
output_prefix_tar = snakemake.params.output_prefix_tar
filtering_method = snakemake.params.filtering_method
out_dir = snakemake.params.outdir
ncores = snakemake.params.ncores
resource_log = os.path.abspath(snakemake.params.resource_log)

def main(net, gene_communities, max_communities, max_res, output_prefix_reg, output_prefix_tar, filtering_method, out_dir, ncores):
    
    selected_file = next((f for f in net if filtering_method.lower() in f.lower()), None)
    if selected_file is None:
        raise ValueError("No valid input file found for filtering method: {}".format(filtering_method))

    selected_file = os.path.abspath(selected_file)
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)

    with measure_resources_context(resource_log):
        bihidef.bihidef(
            filename=selected_file,
            maxres=max_res,
            comm_mult=max_communities,
            oR=output_prefix_reg,
            oT=output_prefix_tar,
            processes=ncores,
        )

if __name__ == "__main__":
    main(net, gene_communities, max_communities, max_res, output_prefix_reg, output_prefix_tar, filtering_method, out_dir, ncores)