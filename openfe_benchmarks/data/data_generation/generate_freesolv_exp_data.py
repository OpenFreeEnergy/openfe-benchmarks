import json
import requests
from openff.toolkit import Molecule
from openff.units import unit
from gufe.tokenization import JSON_HANDLER


def main():
    """
    Extract the reference experimental solvation free energy data from the FreeSolv dataset and save it as a JSON file in a slightly different format with units.

    Notes
    -----
    - to be consistent with the bfe files we use dg and uncertainty for the experimental data and its associated uncertainty.
    - we store the original reference and notes for each entry, as well as the iupac name and the original smiles string. We also store the inchikey and inchi for each molecule to allow matching if we lose the name.

    """
    # tagged to version 0.52
    url = "https://raw.githubusercontent.com/MobleyLab/FreeSolv/refs/tags/v0.52/database.json"
    response = requests.get(url)
    data = response.json()
    # water is the only solvent in the dataset
    water = Molecule.from_smiles("O")
    water_inchi = water.to_inchi(fixed_hydrogens=True)
    water_inchikey = water.to_inchikey(fixed_hydrogens=True)

    ref_data = {}
    for name, entry in data.items():
        # make the molecule from the provided smiles
        off_mol = Molecule.from_smiles(entry["smiles"], allow_undefined_stereo=True)
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
            "solute_inchikey": off_mol.to_inchikey(fixed_hydrogens=True),
            "solute_inchi": off_mol.to_inchi(fixed_hydrogens=True),
            "solvent_smiles": "O",
            "solvent_inchikey": water_inchikey,
            "solvent_inchi": water_inchi,
        }

    with open("experimental_solvation_free_energy_data.json", "w") as f:
        json.dump(ref_data, f, cls=JSON_HANDLER.encoder, indent=4)


if __name__ == "__main__":
    main()
