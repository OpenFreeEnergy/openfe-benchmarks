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
    Generate conformers for an OpenFF ``Molecule`` and return an
    RDKit molecule containing only the lowest-energy conformer.

    Up to 10 conformers are generated via RDKit's ETKDGv3 algorithm with a 0.25 Å
    RMS pruning cutoff. Each conformer is optimized using RDKit's UFF force field.
    The conformer with the lowest post-optimization potential energy is kept; all
    others are discarded.

    Parameters
    ----------
    offmol : openff.toolkit.Molecule
        Input molecule without conformers.

    Returns
    -------
    rdkit.Chem.Mol
        A new RDKit ``Mol`` object containing a single lowest‑energy
        conformer.
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
    - ``ligands.sdf`` — 3‑D coordinates for every unique solute and solvent
      molecule, one lowest‑energy conformer each (OpenFF 2.3.0 + OpenMM
      minimised); only written if the file does not already exist. When this
      file is generated each record is stamped with a
      ``conformer_provenance`` SDF property containing a JSON object with
      keys ``rdkit_version``, ``openff_toolkit_version``,
      ``conformer_generation_method``, ``optimization_method``,
      ``rms_pruning_threshold`` and ``num_conformers_generated``.

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
    sdf_filename = "../benchmark_systems/solvation_set/mnsol_neutral/ligands.sdf"
    flag_sdf = not os.path.isfile(sdf_filename)

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
            # Experimental uncertainties for neutral molecules are set to 0.2 kcal/mol,
            # following the recommendation in the MNSol documentation.
            uncertainty = 0.2 if charge == "0" else 3

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
                    rdmol_solute = _best_conformer_rdmol(offmol_solute)
                except Exception as e:
                    print(
                        f"Failed to find a conformer for solute: {solute_name}, {str(e)}"
                    )
                    skip_molecules.append(solute_name)
                    continue

            offmol_solvent = Molecule.from_smiles(
                NAMES_TO_SMILES[solvent_name],
                allow_undefined_stereo=True,
                name=solvent_name,
            )
            if flag_sdf and solvent_name not in molecules:
                try:
                    rdmol_solvent = _best_conformer_rdmol(offmol_solvent)
                except Exception as e:
                    print(
                        f"Failed to find a conformer for solvent: {solvent_name}, {str(e)}"
                    )
                    skip_molecules.append(solvent_name)
                    continue

            if flag_sdf and solute_name not in molecules:
                molecules[solute_name] = rdmol_solute
            if flag_sdf and solvent_name not in molecules:
                molecules[solvent_name] = rdmol_solvent

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
                "solvent_inchikey": offmol_solvent.to_inchikey(fixed_hydrogens=True),
                "solvent_inchi": offmol_solvent.to_inchi(fixed_hydrogens=True),
            }

    print(f"There are {len(exp_data)} systems, and {len(molecules)} molecules")
    with open(exp_filename, "w") as f:
        json.dump(exp_data, f, cls=JSON_HANDLER.encoder, indent=4)

    if flag_sdf:
        provenance = {
            "rdkit_version": Chem.rdBase.rdkitVersion,
            "openff_toolkit_version": Molecule.__module__,
            "conformer_generation_method": "ETKDGv3",
            "optimization_method": "UFF",
            "rms_pruning_threshold": 0.25,
            "num_conformers_generated": 10,
        }
        with SDWriter(sdf_filename) as writer:
            for molecule in molecules.values():
                molecule.SetProp("conformer_provenance", json.dumps(provenance))
                writer.write(molecule)


if __name__ == "__main__":
    main()
