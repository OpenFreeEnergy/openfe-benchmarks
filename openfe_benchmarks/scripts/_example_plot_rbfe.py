from openfe_benchmarks.data._results_utils import build_femap_from_relative_results
import json
from gufe.tokenization import JSON_HANDLER
from cinnabar import plotting
import pathlib

RESULTS_FILE = "../results/2026-03-18-openmm-840-qa-testing/computational_results.json"
OUTPUT_DIR = "outputs"

def main():
    """
    An example script which can load the calculated DDG values from RBFE calculations and plot vs experimental data for each system.

    Notes:
        - Does not plot the DG values
    """

    # load the results file
    results = json.load(open(RESULTS_FILE), cls=JSON_HANDLER.decoder)
    # check we have DDG values
    if "DDG" not in results:
        raise ValueError(f"Results file {RESULTS_FILE} does not contain 'DDG' values, cannot plot")

    # build FEMaps and load with experimental data
    femaps_by_system = build_femap_from_relative_results(results=results["DDG"])

    output_dir = pathlib.Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    # for each system plot the RBFE results
    for (system_group, system_name), femap in femaps_by_system.items():
        leg_graph = femap.to_legacy_graph()
        plotting.plot_DDGs(
            graph=leg_graph,
            title=f"{system_group}-{system_name}",
            figsize=5,
            scatter_kwargs={"s": 20, "marker": "o"},
            filename=(output_dir / f"{system_group}_{system_name}_DDG.png").as_posix()
        )

if __name__ == "__main__":
    main()