import json
import requests

from pathlib import Path
from openff.toolkit import Molecule
from openff.units import unit
from openff.toolkit.utils.toolkit_registry import ToolkitRegistry
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
from gufe.tokenization import JSON_HANDLER

from openfe_benchmarks.scripts.utils import process_sdf
import openfe_benchmarks

toolkit_registry = ToolkitRegistry([RDKitToolkitWrapper()])


def main():
    """
    Extract the reference experimental solvation free energy data from the FreeSolv dataset and save it as a JSON file in a slightly different format with units.

    InchiKeys are derived from structures in data/benchmark_systems/solvation_set/freesolv/ligands.sdf

    Notes
    -----
    - to be consistent with the bfe files we use dg and uncertainty for the experimental data and its associated uncertainty.
    - we store the original reference and notes for each entry, as well as the iupac name and the original smiles string. We also store the inchikey and inchi for each molecule to allow matching if we lose the name.

    """
    # tagged to version 0.52
    url = "https://raw.githubusercontent.com/MobleyLab/FreeSolv/refs/tags/v0.52/database.json"
    response = requests.get(url)
    data = response.json()

    # get sdf structures
    package_root = Path(openfe_benchmarks.__file__).parent.parent
    sdf_path = (
        package_root
        / "openfe_benchmarks/data/benchmark_systems/solvation_set/freesolv/ligands.sdf"
    )
    ligands = process_sdf(str(sdf_path), return_dict=True)

    # water is the only solvent in the dataset
    water = Molecule.from_smiles("O")
    water_inchi = water.to_inchi(fixed_hydrogens=True)
    water_inchikey = water.to_inchikey(fixed_hydrogens=True)

    ref_data = {}
    for name, entry in data.items():
        # make the molecule from the provided smiles
        off_mol = ligands[name].to_openff()
        # use the solute solvent names as the keys
        ref_data[f"{name},water"] = {
            "dg": entry["expt"] * unit.kilocalories_per_mole,
            "uncertainty": entry["d_expt"] * unit.kilocalories_per_mole,
            "reference": entry["expt_reference"],
            "solute_name": name,
            "solvent_name": "water",
            "solute_iupac": entry["iupac"],
            "notes": entry["notes"],
            "solute_smiles": entry["smiles"],
            # store extra identifiers which can be used to match the molecule if we lose the name
            "solute_inchikey": off_mol.to_inchikey(
                fixed_hydrogens=True, toolkit_registry=toolkit_registry
            ),
            "solute_inchi": off_mol.to_inchi(
                fixed_hydrogens=True, toolkit_registry=toolkit_registry
            ),
            "solvent_smiles": "O",
            "solvent_inchikey": water_inchikey,
            "solvent_inchi": water_inchi,
        }

    out_path = (
        package_root
        / "openfe_benchmarks/data/benchmark_systems/solvation_set/freesolv/experimental_solvation_free_energy_data.json"
    )
    with open(out_path, "w") as f:
        json.dump(ref_data, f, cls=JSON_HANDLER.encoder, indent=4)


if __name__ == "__main__":
    main()
