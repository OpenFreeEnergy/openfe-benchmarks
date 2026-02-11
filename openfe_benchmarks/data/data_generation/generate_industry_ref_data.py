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
from openff.toolkit import Molecule

# name conversions extracted from the industry benchmark files, these will be reversed later
name_conversions = {
    "41 flip": "41-flip",
    "40 flip": "40-flip",
    "38 flip": "38-flip",
    "30 flip": "30-flip",
    "43 flip": "43-flip",
    "47 flip": "47-flip",
    "48 flip": "48-flip",
    "46 flip": "46-flip",
    "36 out": "36o",
    "37 out": "37o",
    "38 out": "38o",
    "39 out": "39o",
    "28 out": "28o",
    "CHEMBL3402756_2.7 redocked": "CHEMBL3402756_2.7_redocked",
    "CHEMBL3402757_6.5 redocked": "CHEMBL3402757_6.5_redocked",
    "CHEMBL3402758_10 redocked": "CHEMBL3402758_10_redocked",
    "CHEMBL3402760_1 redocked": "CHEMBL3402760_1_redocked",
    "CHEMBL3402762_1 redocked": "CHEMBL3402762_1_redocked",
    "CHEMBL3402759_5.7 redocked": "CHEMBL3402759_5.7_redocked",
    "CHEMBL3402761_1 redocked": "CHEMBL3402761_1_redocked",
    "Example 22": "Example-22",
    "Example 23": "Example-23",
    "Example 14": "Example-14",
    "Example 9": "Example-9",
    "SHP099-1 Example 7": "SHP099-1-Example-7",
    "Example 28": "Example-28",
    "Example 24": "Example-24",
    "Example 26": "Example-26",
    "Example 6": "Example-6",
    "Example 1": "Example-1",
    "Example 30": "Example-30",
    "Example 8": "Example-8",
    "Example 29": "Example-29",
    "Example 2": "Example-2",
    "Example 25": "Example-25",
    "Example 4": "Example-4",
    "Example 3": "Example-3",
    "Example 27": "Example-27",
    "Example 5": "Example-5",
    "9 flip": "9-flip",
}


@click.command()
@click.option("--system-group", type=str, required=True, help="The system group of system names used to determine the edges, i.e., JACS/MERCK.")
@click.option("--system-name", type=str, required=True, help="The system name used to determine the edges, i.e., TYK2.")
@click.option("--out-dir", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path), required=True, help="The output dir name")
@click.option("--input-sdf", type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path), required=False, help="The input sdf file containing the ligands from this network, used to extract identifiers.")
def main(system_group: str, system_name: str, out_dir: pathlib.Path, input_sdf: pathlib.Path):
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
    input_sdf
        The input sdf file containing the ligands from this network, used to extract identifiers which will be stored in the reference data.

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
    # load the ligands with openff toolkit to extract the identifiers
    molecules = Molecule.from_file(input_sdf.as_posix())
    molecule_by_name = {molecule.name: molecule for molecule in molecules}
    # reverse the name conversions to match the names in the input sdf
    reversed_name_conversions = dict((value, key) for key, value in name_conversions.items())
    ref_dataset = {}
    for _, row in system_data.iterrows():
        ligand_name = reversed_name_conversions.get(row["Ligand name"], row["Ligand name"])
        if ligand_name not in molecule_by_name:
            # this ligand was not included in the benchmark dataset so skip it
            continue
        ligand_data = {
            "dg": row["Exp. dG (kcal/mol)"] * unit.kilocalories_per_mole,
            "canonical_smiles": molecule_by_name[ligand_name].to_smiles(isomeric=True),
            # do we need to be tautomer specific?
            "inchikey": molecule_by_name[ligand_name].to_inchikey(fixed_hydrogens=True)
        }
        dg_uncertainty = row["Exp. dG error (kcal/mol)"] * unit.kilocalories_per_mole
        if dg_uncertainty.m != 0.0:
            ligand_data["uncertainty"] = dg_uncertainty
        ref_dataset[ligand_name] = ligand_data

    # save the extracted reference data to the output directory as json using gufe serialization
    out_file = out_dir / "experimental_binding_data.json"
    with open(out_file, "w") as f:
        json.dump(ref_dataset, f, indent=4, cls=JSON_HANDLER.encoder)
    print(f"Saved extracted reference data to {out_file}")

if __name__ == "__main__":
    main()