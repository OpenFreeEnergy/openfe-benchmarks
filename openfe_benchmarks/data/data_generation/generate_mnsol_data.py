import json
import csv
import os

import click
import pathlib
from rdkit.Chem import SDWriter
from tqdm import tqdm

from openff.toolkit import Molecule
from openff.units import unit
from gufe.tokenization import JSON_HANDLER
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import AllChem

NAMES_TO_SMILES = json.load(open("mnsol-name-to-smiles.json", "r"))

RDLogger.DisableLog("rdApp.*")


def _best_conformer_rdmol(offmol: "Molecule"):
    """
    Generate conformers and return the OpenFF ``Molecule`` with only the
    lowest-energy conformer retained.

    Up to 10 conformers are generated via RDKit's ETKDGv3 (through
    ``RDKitToolkitWrapper``) with a 0.25 Å RMS pruning cutoff. Each conformer
    is assigned OpenFF 2.3.0 SMIRNOFF parameters and energy-minimised with
    OpenMM's ``LocalEnergyMinimizer`` on the Reference platform. The conformer
    with the lowest post-minimisation potential energy is kept; all others are
    discarded.

    Parameters
    ----------
    offmol : openff.toolkit.Molecule
        Input molecule without conformers. Modified in-place.

    Returns
    -------
    openff.toolkit.Molecule
        The same molecule object with ``_conformers`` replaced by a single
        lowest-energy conformer.
    """
    best_energy = float("inf")
    best_conformer = None
    rdmol = offmol.to_rdkit()
    params = AllChem.ETKDGv3()
    params.pruneRmsThresh = 0.25
    AllChem.EmbedMultipleConfs(rdmol, numConfs=10, params=params)
    for conf_id in range(rdmol.GetNumConformers()):
        try:
            AllChem.UFFOptimizeMolecule(rdmol, confId=conf_id)
            energy = AllChem.UFFGetMoleculeForceField(
                rdmol, confId=conf_id
            ).CalcEnergy()
        except Exception:
            continue
        if energy < best_energy:
            best_energy = energy
            best_conformer = conf_id

    new_rdmol = Chem.Mol(rdmol)
    new_rdmol.RemoveAllConformers()
    conf = (
        rdmol.GetConformer(best_conformer)
        if best_conformer is not None
        else rdmol.GetConformer(0)
    )
    new_conf = Chem.rdchem.Conformer(conf)
    new_rdmol.AddConformer(new_conf, assignId=True)
    return new_rdmol


@click.command()
@click.option(
    "--mnsol-alldata",
    type=click.Path(exists=True, dir_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to MNSol_alldata.txt",
)
def main(mnsol_alldata: pathlib.Path):
    """
    Parse MNSol_alldata.txt and write output files for the ``mnsol_neutral``
    benchmark set.

    Output files written
    --------------------
    - ``experimental_solvation_free_energy_data.json`` — experimental ΔG_solv
      (kcal/mol) with uncertainties (0.2 kcal/mol for neutrals), one entry per
      MNSol row that passes all filters.
    - ``systems_data.json`` — molecule metadata (SMILES, InChI, InChIKey)
      without energetic data; only written if the file does not already exist.
    - ``ligands.sdf`` — 3-D coordinates for every unique solute and solvent
      molecule, one lowest-energy conformer each (OpenFF 2.3.0 + OpenMM
      minimised); only written if the file does not already exist.

    Filters applied
    ---------------
    Only neutral (charge == 0), non-radical, non-dimer entries whose solute
    and solvent names both appear in ``mnsol-name-to-smiles.json`` are
    included. All skipped entries are printed to stdout.

    Parameters
    ----------
    mnsol_alldata : pathlib.Path
        Path to ``MNSol_alldata.txt`` from the MNSol v2012 database
        (https://doi.org/10.13020/3eks-j059).

    Example
    -------
    .. code-block:: bash

        python generate_mnsol_data.py --mnsol-alldata /path/to/MNSol_alldata.txt
    """

    exp_filename = "../benchmark_systems/solvation_set/mnsol_neutral/experimental_solvation_free_energy_data.json"
    sys_filename = "../benchmark_systems/solvation_set/mnsol_neutral/systems_data.json"
    sdf_filename = "../benchmark_systems/solvation_set/mnsol_neutral/ligands.sdf"
    flag_sdf = not os.path.isfile(sdf_filename)

    sys_data = {}
    exp_data = {}
    molecules = {}
    skip_molecules = []
    with mnsol_alldata.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        total_lines = sum(1 for _ in fh) - 1  # subtract header
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        for row in tqdm(reader, total=total_lines):
            if not row or not row.get("FileHandle"):
                continue

            key = f"mnsol-{int(row.get('No.')):04d}"
            solute_name = row.get("SoluteName").strip().strip('"')
            solvent_name = row.get("Solvent").strip().strip('"')
            group = row.get("Level1").strip().strip('"')
            charge = row.get("Charge").strip().strip('"')
            dg = float(row.get("DeltaGsolv").strip().strip('"'))
            uncertainty = 0.2 if charge == 0 else 3

            if solute_name in skip_molecules or solvent_name in skip_molecules:
                continue
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
            if solute_name == "water":
                print(f"Skipping {key}: Water")
                continue

            offmol_solute = Molecule.from_smiles(
                NAMES_TO_SMILES[solute_name],
                allow_undefined_stereo=True,
                name=solute_name,
            )
            if flag_sdf and solute_name not in molecules:
                try:
                    rdmol = _best_conformer_rdmol(offmol_solute)
                except Exception as e:
                    print(
                        f"Failed to find a conformer for solute: {solute_name}, {str(e)}"
                    )
                    skip_molecules.append(solute_name)
                    continue
                molecules[solute_name] = rdmol
            offmol_solvent = Molecule.from_smiles(
                NAMES_TO_SMILES[solvent_name],
                allow_undefined_stereo=True,
                name=solvent_name,
            )
            if flag_sdf and solvent_name not in molecules:
                try:
                    rdmol = _best_conformer_rdmol(offmol_solvent)
                except Exception as e:
                    print(
                        f"Failed to find a conformer for solvent: {solvent_name}, {str(e)}"
                    )
                    skip_molecules.append(solvent_name)
                    continue
                molecules[solvent_name] = rdmol

            exp_data[key] = {
                "mnsol No.": key.split("-")[1],
                "dg": dg * unit.kilocalories_per_mole,
                "uncertainty": uncertainty * unit.kilocalories_per_mole,
                "reference": "https://doi.org/10.13020/3eks-j059",
                "notes": ";".join(f"{k}={v}" for k, v in list(row.items())[:12]),
                "solute_name": solute_name,
                "solvent_name": solvent_name,
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
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "solute_smiles": offmol_solute.to_smiles(explicit_hydrogens=True),
                "solute_inchikey": offmol_solute.to_inchikey(fixed_hydrogens=True),
                "solute_inchi": offmol_solute.to_inchi(fixed_hydrogens=True),
                "solvent_smiles": offmol_solvent.to_smiles(explicit_hydrogens=True),
                "solvent_inchikey": offmol_solute.to_inchikey(fixed_hydrogens=True),
                "solvent_inchi": offmol_solute.to_inchi(fixed_hydrogens=True),
            }

    print(f"There are {len(sys_data)} systems, and {len(molecules)} molecules")
    with open(exp_filename, "w") as f:
        json.dump(exp_data, f, cls=JSON_HANDLER.encoder, indent=4)
    if not os.path.isfile(sys_filename):
        with open(sys_filename, "w") as f:
            json.dump(sys_data, f, cls=JSON_HANDLER.encoder, indent=4)

    if flag_sdf:
        with SDWriter(sdf_filename) as writer:
            for molecule in molecules.values():
                writer.write(molecule)


if __name__ == "__main__":
    main()
