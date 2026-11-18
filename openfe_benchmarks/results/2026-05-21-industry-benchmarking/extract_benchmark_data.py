import pandas as pd
from openff.units import unit
from gufe.tokenization import JSON_HANDLER
import json
from cinnabar import FEMap
import bz2
from pathlib import Path
from collections import defaultdict

from openfe_benchmarks.data import get_data_by_system_name
from openfe_benchmarks.results import BenchmarkResults

def main():
    """
    Create a submission from the industry benchmarking results csv file.
    """
    results_df = pd.read_csv("https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/refs/heads/main/industry_benchmarks/analysis/processed_results/combined_pymbar3_edge_data.csv")

    gathered_results = {"dg": [], "ddg": []}

    # track the system groups and the system names under them
    benchmark_systems = defaultdict(set)

    # some ligands were renamed in openfe-benchmarks after the industry benchmarking
    name_conversions = {
        "41-flip": "41 flip",
        "38-flip": "38 flip",
        "40-flip": "40 flip",
        "30-flip": "30 flip",
        "43-flip": "43 flip",
        "47-flip": "47 flip",
        "48-flip": "48 flip",
        "46-flip": "46 flip",
        "36o": "36 out",
        "37o": "37 out",
        "38o": "38 out",
        "28o": "28 out",
        "39o": "39 out",
        "CHEMBL3402757_6.5_redocked": "CHEMBL3402757_6.5 redocked",
        "CHEMBL3402758_10_redocked": "CHEMBL3402758_10 redocked",
        "CHEMBL3402760_1_redocked": "CHEMBL3402760_1 redocked",
        "CHEMBL3402762_1_redocked": "CHEMBL3402762_1 redocked",
        "CHEMBL3402756_2.7_redocked": "CHEMBL3402756_2.7 redocked",
        "CHEMBL3402759_5.7_redocked": "CHEMBL3402759_5.7 redocked",
        "CHEMBL3402761_1_redocked": "CHEMBL3402761_1 redocked",
        "Example-22": "Example 22",
        "Example-23": "Example 23",
        "Example-14": "Example 14",
        "Example-9": "Example 9",
        "SHP099-1-Example-7": "SHP099-1 Example 7",
        "Example-28": "Example 28",
        "Example-24": "Example 24",
        "Example-26": "Example 26",
        "Example-6": "Example 6",
        "Example-1": "Example 1",
        "Example-30": "Example 30",
        "Example-8": "Example 8",
        "Example-29": "Example 29",
        "Example-2": "Example 2",
        "Example-25": "Example 25",
        "Example-4": "Example 4",
        "Example-3": "Example 3",
        "Example-27": "Example 27",
        "Example-5": "Example 5",
        "9-flip": "9 flip",
    }


    # loop over all the systems in the results dataframe and extract the ddg values
    for system_group in results_df["system group"].unique():
        system_group_df = results_df[results_df["system group"] == system_group]

        system_names = system_group_df["system name"].unique()
        for system_name in system_names:
            print(f"Processing {system_group}:{system_name}")
            system_df = system_group_df[system_group_df["system name"] == system_name]

            # Make an femap for this dataset
            system_femap = FEMap()

            # this group was renamed
            if system_group == "scaffold_hopping_set":
                system_group = "scaffold_hopping"
            if system_group == "scaffold_hopping" and system_name == "factor_xa":
                system_name = "fxa"

            benchmark_systems[system_group].add(system_name)

            # load the experimental data for this system to check that the ligands are present
            ref_data = get_data_by_system_name(system_group, system_name)
            experimental_file = ref_data.reference_data["experimental_binding_data"]
            with open(experimental_file) as f:
                experimental_data = json.load(f, cls=JSON_HANDLER.decoder)

            for _, row in system_df.iterrows():

                ligand_a = name_conversions.get(row["ligand_A"], row["ligand_A"])
                ligand_b = name_conversions.get(row["ligand_B"], row["ligand_B"])
                # check that both the ligands have experimental data not one system has intermediate ligands
                if ligand_a not in experimental_data and "intermediate" not in ligand_a:
                    raise ValueError(f"Ligand {ligand_a} not found in experimental data for {system_group}/{system_name}")
                if ligand_b not in experimental_data and "intermediate" not in ligand_b:
                    raise ValueError(f"Ligand {ligand_b} not found in experimental data for {system_group}/{system_name}")

                # add the data to the femap and to the gathered results
                result_data = {
                    "system_group": system_group,
                    "system_name": system_name,
                    "ligand_a": ligand_a,
                    "ligand_b": ligand_b,
                    "ddg": row["DDG_bind (kcal/mol)"] * unit.kilocalories_per_mole,
                    "ddg_uncertainty": row["stdev DDG_bind (kcal/mol)"] * unit.kilocalories_per_mole,
                    "dgs_complex": [row[f"complex_repeat_{i}_DG (kcal/mol)"] * unit.kilocalories_per_mole for i in range(3)],
                    "dgs_solvent": [row[f"solvent_repeat_{i}_DG (kcal/mol)"] * unit.kilocalories_per_mole for i in range(3)],
                    "complex_mbar_errors": [row[f"complex_repeat_{i}_dDG (kcal/mol)"] * unit.kilocalories_per_mole for i in range(3)],
                    "solvent_mbar_errors": [row[f"solvent_repeat_{i}_dDG (kcal/mol)"] * unit.kilocalories_per_mole for i in range(3)],
                    "complex_smallest_mbar_overlaps": [row[f"complex_repeat_{i}_smallest_overlap"] for i in range(3)],
                    "solvent_smallest_mbar_overlaps": [row[f"solvent_repeat_{i}_smallest_overlap"] for i in range(3)],
                }
                gathered_results["ddg"].append(result_data)

                system_femap.add_relative_calculation(
                    labelA=ligand_a,
                    labelB=ligand_b,
                    value=result_data["ddg"],
                    uncertainty=result_data["ddg_uncertainty"],
                )

            # generate the absolute values for this system
            system_femap.generate_absolute_values()
            abs_df = system_femap.get_absolute_dataframe()

            # extract the absolute values and add them to the gathered results
            for _, row in abs_df.iterrows():
                entry_data = {
                    "ligand": row["label"],
                    "dg": row["DG (kcal/mol)"] * unit.kilocalories_per_mole,
                    "dg_uncertainty": row["uncertainty (kcal/mol)"]
                    * unit.kilocalories_per_mole,
                    "system_group": system_group,
                    "system_name": system_name,
                    "source": row["source"],
                }
                gathered_results["dg"].append(entry_data)

    # write out the data to a bz2-compressed JSON file
    output_file = "computational_results.json.bz2"
    with bz2.open(output_file, "wt") as w:
        json.dump(gathered_results, w, cls=JSON_HANDLER.encoder, indent=4)

    tags = set(["rbfe", "openff-2.2.0", "am1bcc_at", "openfe"] + sorted(benchmark_systems.keys()) + sorted([name for names in benchmark_systems.values() for name in names]))

    # write out the data for the yaml file based on the other submission data
    submission_data = BenchmarkResults(
        submission_id="2026-05-21-industry-benchmarking",
        title="OpenFE RBFE Industry Benchmarking Results - 2026-05-21-industry-benchmarking",
        summary="RBFE results from the OpenFE industry benchmarking dataset.",
        tags=list(tags),
        calculation_type="rbfe",
        authors=[{"name": "Josh Horton"}],
        date="2026-05-21",
        openfe_version="1.0.1",
        openmm_version="8.1.1",
        openff_toolkit_version="0.15.2",
        mapper="KartografAtomMapper 1.0.1 (LSA)",
        forcefield=["ff14SB", "phosaa10", "tip3p_HFE_multivalent", "tip3p_standard"],
        small_molecule_forcefield="openff-2.2.0",
        partial_charges="am1bcc_at",
        benchmark_data={
            "source_repository": "https://github.com/OpenFreeEnergy/IndustryBenchmarks2024",
        },
        results=output_file,
        archive={"doi": "https://zenodo.org/communities/openfreeenergy/records?q=&l=list&p=1&s=10&sort=newest", "archive_provider": "zenodo"},
        license="CC-BY-4.0",
        protocol_settings=[
            {
                "protocol": "RelativeHybridTopologyProtocol",
                "timestep": "4.0 fs",
                "temperature": "298.15 K",
                "pressure": "1 bar",
                "forcefields": ["ff14SB", "phosaa10", "tip3p_HFE_multivalent", "tip3p_standard"],
                "small_molecule_forcefield": "openff-2.2.0",
                "partial_charges": "am1bcc_at",
                "equilibration_time": "1 ns",
                "production_time": "5 ns",
                "vacuum_equilibration_time": "",
                "vacuum_production_time": "",
                "solvent_equilibration_time": "",
                "solvent_production_time": "",
                "lambda_functions": "default",
                "lambda_windows": "11",
                "lambda_schedule": "",
            }
        ],
    )
    submission_data.write_submission_yaml(Path("submission.yaml"))

if __name__ == "__main__":
    main()