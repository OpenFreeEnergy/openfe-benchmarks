import json
import csv
import os

import click
import pathlib

from openff.toolkit import Molecule
from openff.units import unit
from gufe.tokenization import JSON_HANDLER

NAMES_TO_SMILES = json.load(open("mnsol-name-to-smiles.json", "r"))


@click.command()
@click.option(
    "--mnsol-alldata",
    type=click.Path(exists=True, dir_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to MNSol_alldata.txt",
)
def main(mnsol_alldata: pathlib.Path):
    """
    Parse MNSol_alldata.txt and write two JSON files:

    - ``experimental_solvation_free_energy_data.json`` — experimental ΔG_solv
      (kcal/mol) with uncertainties (0.2 for neutrals).
    - ``systems_data.json`` — molecule metadata (SMILES, InChI, InChIKey) without
      energetic data; only written if the file does not already exist.

    Neutral, non-radical, non-dimer entries whose solute and solvent names appear
    in ``mnsol-name-to-smiles.json`` are included; all others are skipped with a
    printed message.

    Example
    -------
    .. code-block:: bash

        python generate_mnsol_data.py --mnsol-alldata /path/to/MNSol_alldata.txt
    """

    exp_filename = "../benchmark_systems/solvation_set/mnsol_neutral/experimental_solvation_free_energy_data.json"
    sys_filename = "../benchmark_systems/solvation_set/mnsol_neutral/systems_data.json"

    sys_data = {}
    exp_data = {}
    with mnsol_alldata.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if not row or not row.get("FileHandle"):
                continue

            key = f"mnsol-{int(row.get('No.')):04d}"
            solute_name = row.get("SoluteName").strip().strip('"')
            solvent_name = row.get("Solvent").strip().strip('"')
            group = row.get("Level1").strip().strip('"')
            charge = row.get("Charge").strip().strip('"')
            dg = float(row.get("DeltaGsolv").strip().strip('"'))
            uncertainty = 0.2 if charge == 0 else 3

            if solute_name not in NAMES_TO_SMILES:
                print(f"Skipping {key}: SMILES for solute name not found")
                continue
            if solvent_name not in NAMES_TO_SMILES:
                print(f"Skipping {key}: SMILES for solvent name not found")
                continue
            if "radical" in solute_name:
                print(f"Skipping {key}: Solute is a radical")
                continue
            if group == "14":
                print(f"Skipping {key}: Solute is a dimer")
                continue
            if charge != "0":
                print(f"Skipping {key}: Charge is not zero")
                continue

            offmol_solute = Molecule.from_smiles(
                NAMES_TO_SMILES[solute_name], allow_undefined_stereo=True
            )
            offmol_solvent = Molecule.from_smiles(
                NAMES_TO_SMILES[solvent_name], allow_undefined_stereo=True
            )

            exp_data[key] = {
                "mnsol No.": key.split("-")[1],
                "dg": dg * unit.kilocalories_per_mole,
                "uncertainty": uncertainty * unit.kilocalories_per_mole,
                "reference": "https://doi.org/10.13020/3eks-j059",
                "notes": ";".join(f"{k}={v}" for k, v in list(row.items())[:12]),
                "solute_charge": charge,
                "solute_smiles": offmol_solute.to_smiles(explicit_hydrogens=True),
                "solute_inchikey": offmol_solute.to_inchikey(fixed_hydrogens=True),
                "solute_inchi": offmol_solute.to_inchi(fixed_hydrogens=True),
                "solvent_smiles": offmol_solvent.to_smiles(explicit_hydrogens=True),
                "solvent_inchikey": offmol_solute.to_inchikey(fixed_hydrogens=True),
                "solvent_inchi": offmol_solute.to_inchi(fixed_hydrogens=True),
            }
            sys_data[key] = {
                "mnsol No.": key.split("-")[1],
                "reference": "https://doi.org/10.13020/3eks-j059",
                "notes": ";".join(f"{k}={v}" for k, v in list(row.items())[:12]),
                "solute_charge": charge,
                "solute_smiles": offmol_solute.to_smiles(explicit_hydrogens=True),
                "solute_inchikey": offmol_solute.to_inchikey(fixed_hydrogens=True),
                "solute_inchi": offmol_solute.to_inchi(fixed_hydrogens=True),
                "solvent_smiles": offmol_solvent.to_smiles(explicit_hydrogens=True),
                "solvent_inchikey": offmol_solute.to_inchikey(fixed_hydrogens=True),
                "solvent_inchi": offmol_solute.to_inchi(fixed_hydrogens=True),
            }

    print(len(sys_data))
    with open(exp_filename, "w") as f:
        json.dump(exp_data, f, cls=JSON_HANDLER.encoder, indent=4)
    if not os.path.isfile(sys_filename):
        with open(sys_filename, "w") as f:
            json.dump(sys_data, f, cls=JSON_HANDLER.encoder, indent=4)


if __name__ == "__main__":
    main()
