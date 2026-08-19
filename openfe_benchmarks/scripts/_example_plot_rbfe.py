"""Example script demonstrating BenchmarkResults API for RBFE plotting.

This script showcases the new BenchmarkResults features:
- Loading results using get_benchmark_results()
- Filtering results with filter_results()
- Using lazy .ddg_femaps property for automatic FEMap generation
- Generating plots using cinnabar

Simply edit the SUBMISSION_ID variable to plot a different submission.
"""

import pathlib

from cinnabar import plotting

from openfe_benchmarks.results import get_benchmark_results

# Edit this to plot a different submission
SUBMISSION_ID = "2026-03-18-openmm-840-qa-testing"
OUTPUT_DIR = "outputs"


def main():
    """Load RBFE results and plot vs experimental data for each system."""

    # 1. Load results using the new get_benchmark_results() factory function
    print(f"Loading submission: {SUBMISSION_ID}")
    results = get_benchmark_results(SUBMISSION_ID)
    print(f"Loaded: {results.title}")
    print(f"Calculation type: {results.calculation_type}")
    print(f"Tags: {', '.join(results.tags)}")

    # Verify this is an RBFE submission
    if "ddg" not in results.raw_results:
        raise ValueError(
            f"Submission {SUBMISSION_ID} does not contain 'ddg' values. "
            f"This script is for RBFE calculations only."
        )

    print("\nGenerating FEMaps...")
    femaps_by_system = results.ddg_femaps
    print(f"Generated {len(femaps_by_system)} FEMaps:")
    for (system_group, system_name), femap in femaps_by_system.items():
        n_edges = femap.n_edges
        print(f"  - {system_group}/{system_name}: {n_edges} edges")

    output_dir = pathlib.Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating plots in {output_dir}/...")
    for (system_group, system_name), femap in femaps_by_system.items():
        output_file = output_dir / f"{system_group}_{system_name}_DDG.png"
        legacy_graph = femap.to_legacy_graph()
        plotting.plot_DDGs(
            graph=legacy_graph,
            title=f"{system_group}-{system_name}",
            figsize=5,
            scatter_kwargs={"s": 20, "marker": "o"},
            filename=output_file.as_posix(),
        )
        print(f"  Created: {output_file}")

    print(f"\nComplete! Generated {len(femaps_by_system)} plots.")


if __name__ == "__main__":
    main()
