"""Generate the ligands.sdf and experimental_solvation_free_energy_data.json for the MNSol dataset.

One can generate experimental_solvation_free_energy_data.json if you have access to the MNSol dataset
in compliance with the dataset license.

Note that the molecules themselves are generated from SMILES strings, with conformers produced from RDKit.
The molecular metadata in experimental_solvation_free_energy_data.json must then be generated from a RDKit
molecule after the conformer stereochemistry is detected, like what will be done when the SDF file is imported
into the openfe toolkit.
"""

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

from openfe_benchmarks.scripts import utils as ofebu

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parent
NAMES_TO_SMILES = json.loads(
    (PACKAGE_ROOT / "mnsol-name-to-smiles.json").read_text(encoding="utf-8")
)

RDLogger.DisableLog("rdApp.*")


def _specificity_key(name: str) -> tuple[int, int, int, str]:
    return (int(name[0].isdigit()), name.count("-"), len(name), name)


def _canonical_smiles(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    return Chem.MolToSmiles(mol, canonical=True)


def _canonicalize_names(
    names_to_smiles: dict[str, str],
) -> tuple[dict[str, str], dict[str, str], set[str], set[str]]:
    name_to_molecule: dict[str, Molecule] = {}
    for name, smiles in names_to_smiles.items():
        try:
            name_to_molecule[name] = Molecule.from_smiles(
                smiles, allow_undefined_stereo=True, name=name
            )
        except Exception:
            continue

    clusters: list[list[str]] = []
    for name, mol in name_to_molecule.items():
        for cluster in clusters:
            representative = name_to_molecule[cluster[0]]
            if Molecule.are_isomorphic(mol, representative, return_atom_map=False)[0]:
                cluster.append(name)
                break
        else:
            clusters.append([name])

    name_to_canonical: dict[str, str] = {}
    for cluster in clusters:
        canonical_name = max(cluster, key=_specificity_key)
        for name in cluster:
            name_to_canonical[name] = canonical_name

    canonical_name_to_smiles = {
        canonical_name: names_to_smiles[canonical_name]
        for canonical_name in set(name_to_canonical.values())
    }

    octane_signature = _canonical_smiles("CCCCCCCC")
    water_signature = _canonical_smiles("O")
    octane_names = {
        name_to_canonical[name]
        for name, smiles in names_to_smiles.items()
        if _canonical_smiles(smiles) == octane_signature
    }
    water_names = {
        name_to_canonical[name]
        for name, smiles in names_to_smiles.items()
        if _canonical_smiles(smiles) == water_signature
    }

    return name_to_canonical, canonical_name_to_smiles, octane_names, water_names


(
    NAME_CANONICAL_NAME,
    CANONICAL_NAME_TO_SMILES,
    OCTANE_CANONICAL_NAMES,
    WATER_CANONICAL_NAMES,
) = _canonicalize_names(NAMES_TO_SMILES)


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
    flag_no_sdf = not os.path.isfile(sdf_filename)

    if flag_no_sdf:
        molecules: dict[str, Chem.rdchem.Mol] = {}
    else:
        molecules_by_smiles = {
            k: m._rdkit
            for k, m in ofebu.process_sdf(
                sdf_filename, return_dict=True, key_type="smiles"
            ).items()
        }
        molecules = {}
        for canonical_name, smiles in CANONICAL_NAME_TO_SMILES.items():
            try:
                offmol = Molecule.from_smiles(
                    smiles,
                    allow_undefined_stereo=True,
                    name=canonical_name,
                )
            except Exception:
                continue
            explicit_smiles = offmol.to_smiles(explicit_hydrogens=True)
            if explicit_smiles in molecules_by_smiles:
                molecules[canonical_name] = molecules_by_smiles[explicit_smiles]

    exp_data = {}
    skip_molecules: list[str] = []
    skipped_systems = []
    with mnsol_alldata.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        total_lines = sum(1 for _ in fh) - 1  # subtract header
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        for row in tqdm(reader, total=total_lines):
            if not row or not row.get("FileHandle"):
                continue

            index = f"{int(row.get('No.')):04d}"
            key = f"mnsol-{index}"
            raw_solute_name = row.get("SoluteName").strip().strip('"')
            raw_solvent_name = row.get("Solvent").strip().strip('"')
            group = row.get("Level1").strip().strip('"')
            charge = row.get("Charge").strip().strip('"')
            dg = float(row.get("DeltaGsolv").strip().strip('"'))
            # Experimental uncertainties for neutral molecules are set to 0.2 kcal/mol,
            # following the recommendation in the MNSol documentation.
            uncertainty = 0.2 if charge == "0" else 3

            if raw_solute_name not in NAMES_TO_SMILES:
                skipped_systems.append(index)
                continue
            if raw_solvent_name not in NAMES_TO_SMILES:
                skipped_systems.append(index)
                continue
            if "radical" in raw_solute_name:
                skipped_systems.append(index)
                continue
            if group == "14":
                skipped_systems.append(index)
                continue
            if charge != "0":
                skipped_systems.append(index)
                continue

            solute_name = NAME_CANONICAL_NAME.get(raw_solute_name)
            solvent_name = NAME_CANONICAL_NAME.get(raw_solvent_name)
            if solute_name is None or solvent_name is None:
                skipped_systems.append(index)
                continue
            if solute_name in skip_molecules or solvent_name in skip_molecules:
                skipped_systems.append(index)
                continue

            # Filter on source SMILES (before 3-D assignment bakes in arbitrary stereo).
            # A racemic solvent measured experimentally as a racemate must not be
            # modelled as a specific enantiopure solvent in simulation.
            solvent_source_smiles = NAMES_TO_SMILES[solvent_name]
            solvent_rdmol = Chem.MolFromSmiles(solvent_source_smiles, sanitize=True)
            if solvent_rdmol is not None:
                Chem.AssignStereochemistry(solvent_rdmol, cleanIt=True, force=True)
                _chiral = Chem.FindMolChiralCenters(
                    solvent_rdmol, includeUnassigned=True, useLegacyImplementation=False
                )
                _stereo = Chem.FindPotentialStereo(solvent_rdmol)
                _undefined = any(tag == "?" for _, tag in _chiral) or any(
                    s.specified == Chem.rdchem.StereoSpecified.Unspecified
                    for s in _stereo
                )
                if _undefined:
                    skipped_systems.append(index)
                    continue

            if (
                solute_name in OCTANE_CANONICAL_NAMES
                and solvent_name in WATER_CANONICAL_NAMES
            ) or (
                solute_name in WATER_CANONICAL_NAMES
                and solvent_name in OCTANE_CANONICAL_NAMES
            ):
                skipped_systems.append(index)
                continue

            offmol_solute = Molecule.from_smiles(
                NAMES_TO_SMILES[solute_name],
                allow_undefined_stereo=True,
                name=solute_name,
            )
            offmol_solvent = Molecule.from_smiles(
                NAMES_TO_SMILES[solvent_name],
                allow_undefined_stereo=True,
                name=solvent_name,
            )

            if flag_no_sdf:
                if solute_name not in molecules:
                    try:
                        rdmol_solute = _best_conformer_rdmol(offmol_solute)
                        Chem.AssignStereochemistryFrom3D(
                            rdmol_solute
                        )  # for comparable inchi and smiles to SDF
                    except Exception as e:
                        print(
                            f"Failed to find a conformer for solute: {solute_name}, {str(e)}"
                        )
                        skip_molecules.append(solute_name)
                        continue

                if solvent_name not in molecules:
                    try:
                        rdmol_solvent = _best_conformer_rdmol(offmol_solvent)
                        Chem.AssignStereochemistryFrom3D(
                            rdmol_solvent
                        )  # for comparable inchi and smiles to SDF
                    except Exception as e:
                        print(
                            f"Failed to find a conformer for solvent: {solvent_name}, {str(e)}"
                        )
                        skip_molecules.append(solvent_name)
                        continue

                if solute_name not in molecules:
                    molecules[solute_name] = rdmol_solute
                if solvent_name not in molecules:
                    molecules[solvent_name] = rdmol_solvent

            if solute_name not in molecules or solvent_name not in molecules:
                continue

            # Recreate the offmols like importing the SDF will to ensure comparable metadata
            offmol_solute = Molecule.from_rdkit(molecules[solute_name])
            offmol_solvent = Molecule.from_rdkit(molecules[solvent_name])

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

    print(f"The following systems were skipped MNSol No.: {skipped_systems}")
    print(f"There are {len(exp_data)} systems, and {len(molecules)} molecules")
    with open(exp_filename, "w") as f:
        json.dump(exp_data, f, cls=JSON_HANDLER.encoder, indent=4)

    if flag_no_sdf:
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
