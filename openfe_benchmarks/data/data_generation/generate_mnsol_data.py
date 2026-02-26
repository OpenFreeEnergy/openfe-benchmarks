import json
import csv
import os
import copy
from collections import defaultdict

import click
import pathlib
from rdkit.Chem import SDWriter
import openmm
import openmm.unit
from tqdm import tqdm

from openff.toolkit import Molecule
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
from openff.units import unit
from gufe.tokenization import JSON_HANDLER
from openff.toolkit.typing.engines.smirnoff import ForceField

NAMES_TO_SMILES = json.load(open("mnsol-name-to-smiles.json", "r"))


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
    offmol.generate_conformers(
        n_conformers=10,
        rms_cutoff=0.25 * unit.angstrom,
        toolkit_registry=RDKitToolkitWrapper(),
    )

    ff = ForceField("openff-2.3.0.offxml")
    best_energy = float("inf") * openmm.unit.kilojoules_per_mole
    for conf_idx in range(len(offmol.conformers)):
        #        try:
        offmol_copy = copy.deepcopy(offmol)
        offmol_copy._conformers = [offmol.conformers[conf_idx]]
        positions = offmol_copy.conformers[0]

        interchange = ff.create_interchange(offmol_copy.to_topology())
        openmm_system = interchange.to_openmm_system(combine_nonbonded_forces=False)
        context = openmm.Context(
            openmm_system,
            openmm.VerletIntegrator(0.1 * openmm.unit.femtoseconds),
            openmm.Platform.getPlatformByName("Reference"),
        )
        context.setPositions(
            (positions * openmm.unit.angstrom).in_units_of(openmm.unit.nanometer),
        )
        openmm.LocalEnergyMinimizer.minimize(
            context=context,
            tolerance=10,
            maxIterations=0,
        )
        energy = context.getState(getEnergy=True).getPotentialEnergy()

        if energy < best_energy:
            best_energy = energy
            best_conformer = offmol_copy.conformers[0]

    offmol._conformers = [best_conformer]
    return offmol


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

    sys_data = {}
    exp_data = {}
    molecules = {}
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
            formula = row.get("Formula").strip().strip('"')
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
            if "si" in formula.lower():
                print(f"Skipping {key}: Silicon in molecule")
                continue

            offmol_solute = Molecule.from_smiles(
                NAMES_TO_SMILES[solute_name],
                allow_undefined_stereo=True,
                name=solute_name,
            )
            if solute_name not in molecules:
                molecules[solute_name] = offmol_solute
            offmol_solvent = Molecule.from_smiles(
                NAMES_TO_SMILES[solvent_name],
                allow_undefined_stereo=True,
                name=solvent_name,
            )
            if solvent_name not in molecules:
                molecules[solvent_name] = offmol_solvent

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

    print(f"There are {len(sys_data)} systems, and {len(molecules)} molecules")
    with open(exp_filename, "w") as f:
        json.dump(exp_data, f, cls=JSON_HANDLER.encoder, indent=4)
    if not os.path.isfile(sys_filename):
        with open(sys_filename, "w") as f:
            json.dump(sys_data, f, cls=JSON_HANDLER.encoder, indent=4)

    errors = defaultdict(list)
    if not os.path.isfile(sdf_filename):
        with SDWriter(sdf_filename) as writer:
            for molecule in tqdm(molecules.values()):
                try:
                    writer.write(_best_conformer_rdmol(molecule).to_rdkit())
                except Exception as e:
                    errors[str(e)[:80]].append(molecule.name)
        for err, tmp in errors.items():
            print(err, tmp)


if __name__ == "__main__":
    main()
