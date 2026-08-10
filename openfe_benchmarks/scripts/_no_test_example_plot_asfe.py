"""Plot ASFEs, only for cinnabar >= 0.6.1, which is not compatible with current env"""

import pathlib
import json
import bz2

from gufe.tokenization import JSON_HANDLER
from cinnabar import plotting

from openfe_benchmarks.scripts._results_utils import build_femap_from_absolute_results

RESULTS_FILE = "../results/2026-08-06-openff-2.3.0-solvation_set_freesolv/computational_results.json.bz2"
OUTPUT_DIR = "outputs"


def _load_results(results_file: str) -> dict:
    results_path = pathlib.Path(results_file)

    if not results_path.exists():
        raise FileNotFoundError(f"Could not find results file: {results_path}")

    open_func = bz2.open if "bz2" in results_file else open

    with open_func(results_path, "rt") as handle:
        return json.load(handle, cls=JSON_HANDLER.decoder)


def main():
    """
    An example script which can load the calculated DG values from ASFE calculations and plot vs experimental solvation data.

    This script creates plots comparing the computed absolute solvation free energies (DG)
    to experimental values for each benchmark system.
    """

    # load the results file, whether compressed or not
    results = _load_results(RESULTS_FILE)
    if "dg" not in results:
        raise ValueError(
            f"Results file {RESULTS_FILE} does not contain 'dg' values, cannot plot"
        )

    # build FEMaps and load with experimental data
    femaps_by_system = build_femap_from_absolute_results(results=results["dg"])

    output_dir = pathlib.Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    # for each system plot the ASFE results compared to experimental data
    for (system_group, system_name), femap in femaps_by_system.items():
        plotting.plot_DGs(
            femap,
            source="Computational",
            title=f"{system_group}-{system_name}",
            figsize=5,
            scatter_kwargs={"s": 20, "marker": "o"},
            filename=(output_dir / f"{system_group}_{system_name}_DG.png").as_posix(),
        )


if __name__ == "__main__":
    main()
