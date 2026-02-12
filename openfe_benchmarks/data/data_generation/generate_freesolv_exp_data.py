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

    ref_data = {}
    for name, entry in data.items():
        # make the molecule from the provided smiles
        off_mol = Molecule.from_smiles(entry["smiles"], allow_undefined_stereo=True)
        ref_data[name] = {
            "dg": entry["expt"] * unit.kilocalories_per_mole,
            "uncertainty": entry["d_expt"] * unit.kilocalories_per_mole,
            "reference": entry["expt_reference"],
            "iupac": entry["iupac"],
            "notes": entry["notes"],
            "smiles": entry["smiles"],
            # store extra identifiers which can be used to match the molecule if we lose the name
            "inchikey": off_mol.to_inchikey(fixed_hydrogens=True),
            "inchi": off_mol.to_inchi(fixed_hydrogens=True),
        }

    with open("experimental_solvation_free_energy.json", "w") as f:
        json.dump(ref_data, f, cls=JSON_HANDLER.encoder, indent=4)



if __name__ == "__main__":
    main()