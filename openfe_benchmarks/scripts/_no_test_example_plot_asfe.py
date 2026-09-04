"""Example script demonstrating BenchmarkResults API for ASFE plotting.

This script showcases the new BenchmarkResults features:
- Loading results using get_benchmark_results()
- Filtering results with filter_results()
- Using lazy dg_femaps() method for automatic FEMap generation
- Generating plots using cinnabar

Simply edit the SUBMISSION_ID variable to plot a different submission.

Note: Requires cinnabar >= 0.6.1 for ASFE plotting support.
"""

import pathlib

from cinnabar import plotting

from openfe_benchmarks.results import get_benchmark_results

# Edit this to plot a different submission
SUBMISSION_ID = "2026-08-06-openff-2.3.0-solvation_set_freesolv"
OUTPUT_DIR = "outputs"


def main():
    """Load ASFE results and plot vs experimental solvation data for each system."""

    print(f"Loading submission: {SUBMISSION_ID}")
    results = get_benchmark_results(SUBMISSION_ID)
    print(f"Loaded: {results.title}")
    print(f"Calculation type: {results.calculation_type}")
    print(f"Tags: {', '.join(results.tags)}")

    # Verify this is an ASFE submission
    if "dg" not in results.raw_results:
        raise ValueError(
            f"Submission {SUBMISSION_ID} does not contain 'dg' values. "
            f"This script is for ASFE calculations only."
        )

    print("\nGenerating FEMaps...")
    femaps_by_system = results.dg_femaps()
    print(f"Generated {len(femaps_by_system)} FEMaps:")
    for (system_group, system_name), femap in femaps_by_system.items():
        n_values = len(femap.graph.nodes)
        print(f"  - {system_group}/{system_name}: {n_values} calculations")

    output_dir = pathlib.Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating plots in {output_dir}/...")
    for (system_group, system_name), femap in femaps_by_system.items():
        output_file = output_dir / f"{system_group}_{system_name}_DG.png"
        plotting.plot_DGs(
            femap,
            source=results.submission_id,
            title=f"{system_group}-{system_name}",
            figsize=5,
            scatter_kwargs={"s": 20, "marker": "o"},
            filename=output_file.as_posix(),
        )
        print(f"  Created: {output_file}")

    print(f"\nComplete! Generated {len(femaps_by_system)} plots.")


if __name__ == "__main__":
    main()
