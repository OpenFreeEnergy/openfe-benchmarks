"""
For the given industry system group ie JACS/MERCK and system name ie TYK2 generate the LOMAP networks and atom mappings.
"""

import click
import pathlib
from kartograf import KartografAtomMapper
from gufe import LigandNetwork, SmallMoleculeComponent
import pandas as pd
from rdkit import Chem
from importlib.metadata import version

# some ligands reported in the industry benchmarks have names that differ from those in the input sdf which were edited by the submitters to be more file-system friendly
# we provide a mapping here to convert the names accordingly
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
@click.option(
    "--system-group",
    type=str,
    required=True,
    help="The industry system group of system names used to determine the edges, i.e., JACS/MERCK.",
)
@click.option(
    "--system-name",
    type=str,
    required=True,
    help="The industry system name used to determine the edges, i.e., TYK2.",
)
@click.option(
    "--input-sdf",
    type=click.Path(
        exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path
    ),
    required=False,
    help="The input sdf file containing the ligands to be used in this network.",
)
@click.option(
    "--out-dir",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path
    ),
    required=True,
    help="The output dir name",
)
def main(
    system_group: str, system_name: str, input_sdf: pathlib.Path, out_dir: pathlib.Path
):
    """
    Generate LOMAP networks and atom mappings for the given industry system group and name.

    Parameters
    ----------
    system_group : str
        The industry system group ie JACS/MERCK used to determine the edges.
    system_name : str
        The industry system name ie TYK2 used to determine the edges.
    input_sdf : pathlib.Path
        The input sdf file containing the ligands to be used in this network.
    out_dir : pathlib.Path
        The output dir name where the generated LOMAP network will be saved.

    Notes
    -----
    - The edges will be extracted from the industry benchmark results
    - New mappings will be generated using Kartograf with the same default settings as used in the industry benchmarks
    - The generated LOMAP network will be saved to the given output directory as 'industry_benchmarks_network.json'

    Raises
    ------
    ValueError
        If the system group or system name is not found in the reference data.
    RuntimeError
        If the ligands in the input SDF do not match the ligands in the reference edges, this is checked by name.
    """
    # load the ref data stored on github
    all_edges = pd.read_csv(
        "https://raw.githubusercontent.com/OpenFreeEnergy/IndustryBenchmarks2024/refs/tags/v1.0.0/industry_benchmarks/analysis/processed_results/combined_pymbar3_edge_data.csv",
        dtype={
            "ligand_A": str,
            "ligand_B": str,
            "system group": str,
            "system name": str,
        },
    )
    available_groups = all_edges["system group"].unique().tolist()
    if system_group not in available_groups:
        raise ValueError(
            f"System group {system_group} not found. Available groups: {available_groups}"
        )
    # check the system name
    group_edges = all_edges[all_edges["system group"] == system_group].reset_index(
        drop=True
    )
    available_names = group_edges["system name"].unique().tolist()
    if system_name not in available_names:
        raise ValueError(
            f"System name {system_name} not found in group {system_group}. Available names: {available_names}"
        )
    # filter to the given system name
    system_edges = group_edges[group_edges["system name"] == system_name].reset_index(
        drop=True
    )

    # load the ligands to generate the mappings and strip charges if present
    with Chem.SDMolSupplier(str(input_sdf), removeHs=False) as suppl:
        ligands = []
        for mol in suppl:
            if mol is not None:
                smc = SmallMoleculeComponent(mol)

                # check for charges and strip if present
                off_mol = smc.to_openff()
                if off_mol.partial_charges is not None:
                    print("Stripping partial charges from ligand:", smc.name, smc.key)
                    # set the charge to None
                    off_mol.partial_charges = None
                    # remove the provenance info
                    _ = off_mol.properties.pop("charge_provenance")
                    # remove the dprop info used for rdkit conversions
                    _ = off_mol.properties.pop("atom.dprop.PartialCharge")
                    smc = SmallMoleculeComponent.from_openff(off_mol)

                ligands.append(smc)

    # create a mapping of ligand name to ligand component, use the converted names to match the edges in the reference data
    ligands_by_name = {
        name_conversions.get(ligand.name, ligand.name): ligand for ligand in ligands
    }

    # check that the ligands in the edges are present in the input sdf
    edge_ligand_names = set(system_edges["ligand_A"]).union(
        set(system_edges["ligand_B"])
    )
    input_ligand_names = set(ligands_by_name.keys())
    if edge_ligand_names != input_ligand_names:
        raise RuntimeError(
            f"Ligands in edges do not match ligands in input sdf. Edge ligands: {edge_ligand_names}, Input ligands: {input_ligand_names}"
        )

    # generate the mappings using kartograf
    mapper = KartografAtomMapper(map_hydrogens_on_hydrogens_only=True)
    # workaround missing __version__ attribute in kartograf
    mapper_provenance = {
        "mapper_settings": mapper.to_dict(),
        "mapper_version": version("kartograf"),
    }
    edges = []
    print(
        f"generating mappings for system group {system_group}, system name {system_name}..."
    )
    for _, row in system_edges.iterrows():
        ligand_a = ligands_by_name[row["ligand_A"]]
        ligand_b = ligands_by_name[row["ligand_B"]]
        mapping = next(mapper.suggest_mappings(ligand_a, ligand_b))
        edges.append(mapping.with_annotations(annotations=mapper_provenance))

    # create the LOMAP network
    network = LigandNetwork(edges=edges)
    # make sure the network has the same number of edges as what we expected
    if len(network.edges) != len(system_edges):
        raise RuntimeError(
            f"Generated network has {len(network.edges)} edges, but expected {len(system_edges)} edges."
        )

    # save the network
    out_path = out_dir / "industry_benchmarks_network.json"
    network.to_json(out_path)
    print(
        f"Industry benchmarks LOMAP style network with kartograf mapping saved to {out_path}"
    )


if __name__ == "__main__":
    main()
