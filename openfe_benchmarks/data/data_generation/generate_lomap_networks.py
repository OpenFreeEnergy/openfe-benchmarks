"""
For the given system generate the LOMAP networks and atom mappings.
"""

import click
import pathlib
from kartograf import KartografAtomMapper
from gufe import SmallMoleculeComponent
import pandas as pd
import openfe
import lomap
from openfe.setup.ligand_network_planning import generate_lomap_network
from rdkit import Chem
from importlib.metadata import version


@click.command()
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
def main(input_sdf: pathlib.Path, out_dir: pathlib.Path):
    """
    Generate LOMAP networks and atom mappings for the given ligands.

    Parameters
    ----------
    input_sdf : pathlib.Path
        The input sdf file containing the ligands to be used in this network.
    out_dir : pathlib.Path
        The output dir name where the generated LOMAP network will be saved.

    Notes
    -----
    - New mappings will be generated using Kartograf with the same default settings as used in the industry benchmarks
    - The generated LOMAP network will be saved to the given output directory as 'lomap_network.json'
    """
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

    # generate the mappings using kartograf
    mapper = KartografAtomMapper(map_hydrogens_on_hydrogens_only=True)

    print(
        f"generating mappings for system group {system_group}, system name {system_name}..."
    )
    # Generate the LOMAP network
    network = generate_lomap_network(
        ligands=ligands,
        scorer=openfe.lomap_scorers.default_lomap_score,
        mappers=[mapper],
    )
    network_provenance = {
        "network_method": "LOMAP network",
        "mapper_settings": mapper.to_dict(),
        "mapper_version": version("kartograf"),
        "scorer": "default_lomap_score",
        "lomap_version": version("lomap2"),
    }

    for edge in network.edges:
        edge.annotations.update(network_provenance)

    # save the network
    out_path = out_dir / "lomap_network.json"
    network.to_json(out_path)
    print(
        f"LOMAP network with kartograf mapping saved to {out_path}"
    )


if __name__ == "__main__":
    main()
