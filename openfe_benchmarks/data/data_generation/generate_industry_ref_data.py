"""
For the given industry system group ie JACS/MERCK and system name ie TYK2 extract the reference data for the ligands in
the system and save it to the output directory.

The reference data will be saved to JSON using gufe which allows us to include units. The structure is a dictionary
where each key is the ligand name and the value is a dict of the experimental DG and the associated uncertainty.
"""
import click
import pathlib
from gufe.tokenization import JSON_HANDLER
import pandas as pd
from openff.units import unit
import json


@click.command()
@click.option("--system-group", type=str, required=True, help="The system group of system names used to determine the edges, i.e., JACS/MERCK.")
@click.option("--system-name", type=str, required=True, help="The system name used to determine the edges, i.e., TYK2.")
@click.option("--out-dir", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path), required=True, help="The output dir name")
def main(system_group: str, system_name: str, out_dir: pathlib.Path):
    """
    Extract the reference data for the given industry benchmark system.

    Parameters
    ----------
    system_group : str
        The industry system group ie JACS/MERCK used to find the reference data.
    system_name : str
        The industry system name ie TYK2 used to find the reference data.
    out_dir : pathlib.Path
        The output dir name where the extracted reference data will be saved.

    Raises
    ------
    ValueError
        If the system group or system name is not found in the reference data.
    """
    # load the ref data from the industry benchmark results stored on github
    ref_dg_data = pd.read_csv("https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/refs/tags/v1.0.0/industry_benchmarks/analysis/schrodinger_21_4_results/combined_schrodinger_dg.csv")
    available_groups = ref_dg_data["system group"].unique().tolist()
    if system_group not in available_groups:
        raise ValueError(f"System group {system_group} not found. Available groups: {available_groups}")
    # check the system name
    group_data = ref_dg_data[ref_dg_data["system group"] == system_group].reset_index(drop=True)
    available_names = group_data["system name"].unique().tolist()
    if system_name not in available_names:
        raise ValueError(f"System name {system_name} not found in group {system_group}. Available names: {available_names}")
    # filter to the given system name
    system_data = group_data[group_data["system name"] == system_name].reset_index(drop=True)
    print(f"Extracted reference data for system group {system_group} and system name {system_name} with {len(system_data)} ligands.")
    ref_dataset = {}
    for _, row in system_data.iterrows():
        ligand_name = row["Ligand name"]
        dg = row["Exp. dG (kcal/mol)"] * unit.kilocalories_per_mole
        dg_uncertainty = row["Exp. dG error (kcal/mol)"] * unit.kilocalories_per_mole
        ref_dataset[ligand_name] = {"dg": dg, "uncertainty": dg_uncertainty}

    # save the extracted reference data to the output directory as json using gufe serialization
    out_file = out_dir / "experimental_binding_data.json"
    with open(out_file, "w") as f:
        json.dump(ref_dataset, f, indent=4, cls=JSON_HANDLER.encoder)
    print(f"Saved extracted reference data to {out_file}")

if __name__ == "__main__":
    main()