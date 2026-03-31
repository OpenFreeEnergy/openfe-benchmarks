"""Build aligned MNSol and FreeSolv benchmark subsets with overlap-first seeding.

Compared with the openff-sage subsets
(https://github.com/openforcefield/openff-sage), this script changes the MNSol
selection order and objective. This script first detects MNSol-FreeSolv solute
overlap using OpenFF isomorphism and uses those overlapping solutes to seed
selection. In the openff-sage subsets, MNSol subset construction is not driven
by prior overlap with FreeSolv.

For MNSol, this script limits each solute to at most two selected
solute/solvent systems. Selection uses RDKit Morgan fingerprints and Tanimoto
distance in two stages: solvent-fingerprint diversity is used when choosing the
seeded solvent rows per overlapping solute, and solute-fingerprint diversity is
used during subsequent environment-coverage expansion. This differs from
openff-sage subset selection, which relies on OpenEye fingerprints via
SelectSubstances.

After overlap-seeded MNSol selection, this script performs a second pass to
fill missing checkmol chemical-environment coverage (up to a per-environment
target) subject to the per-solute cap and distance-based ranking.

For FreeSolv, this script applies the same chemical filtering rules used here
for MNSol-style curation, then performs overlap-seeded selection and an
independent chemical-environment fill step. In contrast, the openff-sage
FreeSolv subset is primarily an overlap-membership filter against MNSol
solutes and does not perform an additional independent chemical-space fill
stage.
"""

import csv
import json
import pathlib
from collections import defaultdict
from pathlib import Path
from typing import Any

import click
import requests
from gufe.tokenization import JSON_HANDLER
from openff.toolkit import Molecule
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem
from tqdm import tqdm
from yammbs import checkmol
from yammbs.checkmol import ChemicalEnvironment

import openfe_benchmarks

RDLogger.DisableLog("rdApp.*")
PACKAGE_ROOT = Path(openfe_benchmarks.__file__).parent.parent

_MN_SOL_JSON_PATH = pathlib.Path(__file__).parent / "mnsol-name-to-smiles.json"
with _MN_SOL_JSON_PATH.open("r", encoding="utf-8") as _f:
    NAMES_TO_SMILES = json.load(_f)

SMIRKS_FILTERS = [
    # Long chain alkane / ether
    "-".join(["[#6X4,#8X2]"] * 10),
    # 1,3 carbonyls with at least one ketone carbonyl
    "[#6](=[#8])-[#6](-[#1])(-[#1])-[#6](=[#8])-[#6]",
    # Things which dissociate
    "[H]-[Cl]",
    "[H]-[Br]",
]

SMIRKS_FILTER_MOLS = [Chem.MolFromSmarts(s) for s in SMIRKS_FILTERS]
ALLOWED_ELEMENTS = {"C", "O", "N", "Cl", "Br", "H", "S", "P"}

CHEMICAL_ENVIRONMENTS = [
    ChemicalEnvironment.Alkane,
    ChemicalEnvironment.Alkene,
    ChemicalEnvironment.Alcohol,
    ChemicalEnvironment.CarbonylHydrate,
    ChemicalEnvironment.Hemiacetal,
    ChemicalEnvironment.Acetal,
    ChemicalEnvironment.Hemiaminal,
    ChemicalEnvironment.Aminal,
    ChemicalEnvironment.Thioacetal,
    ChemicalEnvironment.CarboxylicAcidEster,
    ChemicalEnvironment.Ether,
    ChemicalEnvironment.Aldehyde,
    ChemicalEnvironment.Ketone,
    ChemicalEnvironment.Aromatic,
    ChemicalEnvironment.CarboxylicAcidPrimaryAmide,
    ChemicalEnvironment.CarboxylicAcidSecondaryAmide,
    ChemicalEnvironment.CarboxylicAcidTertiaryAmide,
    ChemicalEnvironment.PrimaryAmine,
    ChemicalEnvironment.SecondaryAmine,
    ChemicalEnvironment.TertiaryAmine,
    ChemicalEnvironment.Cyanate,
    ChemicalEnvironment.Isocyanate,
    ChemicalEnvironment.Heterocycle,
    ChemicalEnvironment.AlkylFluoride,
    ChemicalEnvironment.ArylFluoride,
    ChemicalEnvironment.AlkylChloride,
    ChemicalEnvironment.ArylChloride,
    ChemicalEnvironment.AlkylBromide,
    ChemicalEnvironment.ArylBromide,
    ChemicalEnvironment.Thiol,
    ChemicalEnvironment.Thioaldehyde,
    ChemicalEnvironment.Thioketone,
    ChemicalEnvironment.Thioether,
    ChemicalEnvironment.Disulfide,
    ChemicalEnvironment.Thiourea,
    ChemicalEnvironment.Thiocyanate,
    ChemicalEnvironment.Isothiocyanate,
    ChemicalEnvironment.ThiocarboxylicAcid,
    ChemicalEnvironment.ThiocarboxylicAcidEster,
    ChemicalEnvironment.SulfonicAcidEster,
    ChemicalEnvironment.Sulfone,
    ChemicalEnvironment.PhosphonicAcid,
    ChemicalEnvironment.PhosphoricAcid,
    ChemicalEnvironment.PhosphoricAcidEster,
]


def contains_any_filter_smirks(smiles: str) -> bool:
    """Return True if the SMILES matches any SMIRKS/SMARTS filter pattern."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        return False

    if mol is None:
        return False

    # Make hydrogens explicit so patterns containing [#1]/[H] can match.
    mol = Chem.AddHs(mol)

    return any(
        pattern is not None and mol.HasSubstructMatch(pattern)
        for pattern in SMIRKS_FILTER_MOLS
    )


def has_undefined_stereochemistry(smiles: str) -> bool:
    """Return True if a sanitized RDKit molecule has unspecified stereochemistry."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        return False

    if mol is None:
        return False

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    chiral_centers = Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    )
    if any(tag == "?" for _, tag in chiral_centers):
        return True

    return any(
        stereo_info.specified == Chem.rdchem.StereoSpecified.Unspecified
        for stereo_info in Chem.FindPotentialStereo(mol)
    )


def has_only_allowed_elements(
    smiles: str, allowed_elements: set[str] = ALLOWED_ELEMENTS
) -> bool:
    """Return True if the molecule contains only the allowed element symbols."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        return False

    if mol is None:
        return False

    mol = Chem.AddHs(mol)
    return all(atom.GetSymbol() in allowed_elements for atom in mol.GetAtoms())


def get_checkmol_environments(smiles: str) -> set[ChemicalEnvironment]:
    """Return the set of checkmol chemical environments for a SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        return set()

    if mol is None:
        return set()

    candidate_functions = (
        "find_chemical_environments",
        "get_chemical_environments",
        "chemical_environments",
        "analyze_functional_groups",
    )
    for fn_name in candidate_functions:
        fn = getattr(checkmol, fn_name, None)
        if fn is None:
            continue

        for arg in (smiles, mol):
            try:
                result = fn(arg)
            except Exception:
                continue

            if result is None:
                continue

            environments: set[ChemicalEnvironment] = set()
            if isinstance(result, dict):
                items = [k for k, v in result.items() if v]
            else:
                items = list(result)

            for item in items:
                if isinstance(item, ChemicalEnvironment):
                    environments.add(item)
                elif isinstance(item, str):
                    try:
                        environments.add(ChemicalEnvironment[item])
                    except Exception:
                        try:
                            environments.add(ChemicalEnvironment(item))
                        except Exception:
                            pass

            return environments

    return set()


def canonicalize_smiles(smiles: str) -> str | None:
    """Return an RDKit canonical isomeric SMILES for matching."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        return None

    if mol is None:
        return None

    return Chem.MolToSmiles(mol, isomericSmiles=True)


def make_openff_molecule(smiles: str) -> Molecule | None:
    """Build an OpenFF Molecule for isomorphism checks."""
    try:
        return Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    except Exception:
        return None


def molecules_are_isomorphic(mol_a: Molecule, mol_b: Molecule) -> bool:
    """Check OpenFF isomorphism while being tolerant to toolkit API differences."""
    return Molecule.are_isomorphic(mol_a, mol_b, return_atom_map=False)[0]


def solute_fingerprint(smiles: str):
    """Return a Morgan fingerprint used for diversity tie-breaking."""
    if smiles is None:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        return None

    if mol is None:
        return None

    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def select_diverse_keys(
    candidate_keys: list[str],
    already_selected_keys: set[str],
    fp_by_key: dict[str, Any],
    order_map: dict[str, int],
    n_to_select: int,
) -> list[str]:
    """Select keys maximizing minimum distance to already selected entries."""
    available = [k for k in candidate_keys if k not in already_selected_keys]
    chosen: list[str] = []

    while available and len(chosen) < n_to_select:
        reference_keys = [*already_selected_keys, *chosen]

        if not reference_keys:
            best_key = min(available, key=lambda k: order_map[k])
        else:

            def min_distance_to_references(key: str) -> float:
                fp = fp_by_key.get(key)
                if fp is None:
                    return 0.0

                distances = []
                for ref_key in reference_keys:
                    ref_fp = fp_by_key.get(ref_key)
                    if ref_fp is None:
                        distances.append(0.0)
                    else:
                        similarity = DataStructs.TanimotoSimilarity(fp, ref_fp)
                        distances.append(1.0 - similarity)

                return min(distances) if distances else 0.0

            best_key = max(
                available,
                key=lambda k: (min_distance_to_references(k), -order_map[k]),
            )

        chosen.append(best_key)
        available.remove(best_key)

    return chosen


def fill_environment_coverage(
    data: dict[str, dict],
    seed_keys: set[str],
    n_per_environment: int,
    max_rows_per_solute: int | None = None,
    environment_mode: str = "solvent_solute_pairs",
) -> tuple[dict[str, dict], list[str], list[str]]:
    """Fill environment coverage using a SelectSubstances-like loop.

    The selection iterates over either:

    - solvent/solute environment tuples (``environment_mode='solvent_solute_pairs'``)
    - solute-only environments (``environment_mode='solute_only'``)

    It then selects up to ``n_per_environment`` systems per target using
    fingerprint diversity and removes newly selected systems from the candidate
    pool before moving to the next target.
    """
    target_environment_names = [env.name for env in CHEMICAL_ENVIRONMENTS]
    ordered_keys = list(data)

    def build_environment_tuples() -> list[tuple[str, str]]:
        return [
            *[(env_name, env_name) for env_name in target_environment_names],
            *[
                pair
                for i, env_i in enumerate(target_environment_names)
                for pair in [
                    (env_i, env_j) for env_j in target_environment_names[i + 1 :]
                ]
            ],
        ]

    def matches_environment_tuple(
        record: dict, environment_tuple: tuple[str, str]
    ) -> bool:
        solvent_envs = set(record.get("solvent_env", []))
        solute_envs = set(record.get("solute_env", []))
        env_a, env_b = environment_tuple

        if env_a == env_b:
            return env_a in solvent_envs and env_a in solute_envs

        return (env_a in solvent_envs and env_b in solute_envs) or (
            env_b in solvent_envs and env_a in solute_envs
        )

    environment_tuples = build_environment_tuples()
    if environment_mode not in {"solvent_solute_pairs", "solute_only"}:
        raise ValueError(
            "environment_mode must be 'solvent_solute_pairs' or 'solute_only'"
        )

    fp_by_key = {key: solute_fingerprint(data[key]["solute_smiles"]) for key in data}
    order_map = {key: i for i, key in enumerate(ordered_keys)}
    key_to_envs = {key: set(data[key].get("solute_env", [])) for key in data}
    key_to_solute = {key: data[key]["solute_name"] for key in data}

    selected: set[str] = {key for key in seed_keys if key in data}
    component_pool_keys = [key for key in ordered_keys if key not in selected]
    coverage_additions: list[str] = []
    solute_counts = defaultdict(int)
    for key in selected:
        solute_counts[key_to_solute[key]] += 1

    def count_for_environment(env_name: str) -> int:
        return sum(1 for key in selected if env_name in key_to_envs[key])

    def count_for_environment_tuple(environment_tuple: tuple[str, str]) -> int:
        return sum(
            1
            for key in selected
            if matches_environment_tuple(data[key], environment_tuple)
        )

    if environment_mode == "solute_only":
        for env_name in target_environment_names:
            needed = n_per_environment - count_for_environment(env_name)
            if needed <= 0:
                continue

            candidate_keys = [
                key
                for key in component_pool_keys
                if env_name in key_to_envs[key]
                and (
                    max_rows_per_solute is None
                    or solute_counts[key_to_solute[key]] < max_rows_per_solute
                )
            ]

            picked = select_diverse_keys(
                candidate_keys,
                selected,
                fp_by_key,
                order_map,
                needed,
            )

            selected.update(picked)
            coverage_additions.extend(picked)
            for key in picked:
                solute_counts[key_to_solute[key]] += 1
            if picked:
                picked_set = set(picked)
                component_pool_keys = [
                    key for key in component_pool_keys if key not in picked_set
                ]

    else:
        for environment_tuple in environment_tuples:
            needed = n_per_environment - count_for_environment_tuple(environment_tuple)
            if needed <= 0:
                continue

            candidate_keys = [
                key
                for key in component_pool_keys
                if matches_environment_tuple(data[key], environment_tuple)
                and (
                    max_rows_per_solute is None
                    or solute_counts[key_to_solute[key]] < max_rows_per_solute
                )
            ]

            picked = select_diverse_keys(
                candidate_keys,
                selected,
                fp_by_key,
                order_map,
                needed,
            )

            selected.update(picked)
            coverage_additions.extend(picked)
            for key in picked:
                solute_counts[key_to_solute[key]] += 1
            if picked:
                picked_set = set(picked)
                component_pool_keys = [
                    key for key in component_pool_keys if key not in picked_set
                ]

    missing_environments = [
        env_name
        for env_name in target_environment_names
        if count_for_environment(env_name) < n_per_environment
    ]

    selected_data = {key: data[key] for key in sorted(selected)}
    return selected_data, sorted(coverage_additions), missing_environments


def build_filtered_mnsol_records(
    mnsol_alldata: pathlib.Path,
) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Load and filter MNSol rows with the same gates as the existing script."""
    data: dict[str, dict] = {}
    skipped_systems = defaultdict(list)
    skip_molecules = ["water"]

    with mnsol_alldata.open("r", encoding="utf-8", errors="replace") as fh:
        total_lines = sum(1 for _ in fh)
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        for row in tqdm(reader, total=total_lines):
            if not row or not row.get("FileHandle"):
                continue

            index_raw = row.get("No.")
            if index_raw is None:
                continue
            index = f"{int(index_raw):04d}"
            key = f"mnsol-{index}"
            solute_name = row.get("SoluteName", "").strip().strip('"')
            solvent_name = row.get("Solvent", "").strip().strip('"')
            group = row.get("Level1", "").strip().strip('"')
            charge = row.get("Charge", "").strip().strip('"')

            checks = [
                (solute_name in skip_molecules, "solute in skip_molecules"),
                (solvent_name in skip_molecules, "solvent in skip_molecules"),
                (solute_name not in NAMES_TO_SMILES, "no solute smiles"),
                (solvent_name not in NAMES_TO_SMILES, "no solvent smiles"),
                ("radical" in solute_name, "solute is a radical"),
                (group == "14", "explicit solvent"),
                (charge != "0", "charged system"),
                (solute_name == solvent_name, "solute and solvent are the same"),
            ]
            skip_reason = next(
                (reason for condition, reason in checks if condition), None
            )
            if skip_reason:
                skipped_systems[skip_reason].append(index)
                continue

            solute_smiles = NAMES_TO_SMILES[solute_name]
            solvent_smiles = NAMES_TO_SMILES[solvent_name]
            checks = [
                (
                    contains_any_filter_smirks(solute_smiles),
                    "solute has disqualifying smirks",
                ),
                (
                    contains_any_filter_smirks(solvent_smiles),
                    "solvent has disqualifying smirks",
                ),
                (
                    has_undefined_stereochemistry(solvent_smiles),
                    "solvent has undefined stereochemistry",
                ),
                (
                    has_undefined_stereochemistry(solute_smiles),
                    "solute has undefined stereochemistry",
                ),
                (
                    not has_only_allowed_elements(solvent_smiles),
                    "solvent elements outside chemical space",
                ),
                (
                    not has_only_allowed_elements(solute_smiles),
                    "solute elements outside chemical space",
                ),
            ]
            skip_reason = next(
                (reason for condition, reason in checks if condition), None
            )
            if skip_reason:
                skipped_systems[skip_reason].append(index)
                continue

            data[key] = {
                "mnsol No.": key.split("-")[1],
                "solute_name": solute_name,
                "solvent_name": solvent_name,
                "solute_smiles": solute_smiles,
                "solvent_smiles": solvent_smiles,
                "solute_charge": charge,
                "solute_env": sorted(
                    env.name for env in get_checkmol_environments(solute_smiles)
                ),
                "solvent_env": sorted(
                    env.name for env in get_checkmol_environments(solvent_smiles)
                ),
            }

    return data, skipped_systems


def build_filtered_freesolv_records() -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Load and filter FreeSolv rows with the same gates as the existing script."""
    url = "https://raw.githubusercontent.com/MobleyLab/FreeSolv/refs/tags/v0.52/database.json"
    response = requests.get(url)
    response.raise_for_status()
    freesolv = response.json()

    data: dict[str, dict] = {}
    skipped_systems = defaultdict(list)
    for name, entry in freesolv.items():
        offmol_solute = make_openff_molecule(entry["smiles"])
        if offmol_solute is None:
            skipped_systems["invalid openff solute"].append(name)
            continue

        charge = sum(a.GetFormalCharge() for a in offmol_solute.to_rdkit().GetAtoms())
        checks = [
            (charge != 0, "charged system"),
            (
                contains_any_filter_smirks(entry["smiles"]),
                "solute has disqualifying smirks",
            ),
            (
                has_undefined_stereochemistry(entry["smiles"]),
                "solute has undefined stereochemistry",
            ),
            (
                not has_only_allowed_elements(entry["smiles"]),
                "solute elements outside chemical space",
            ),
        ]
        skip_reason = next((reason for condition, reason in checks if condition), None)
        if skip_reason:
            skipped_systems[skip_reason].append(name)
            continue

        key = f"freesolv-{name}"
        data[key] = {
            "freesolv_id": name,
            "solute_name": name,
            "solvent_name": "water",
            "solute_smiles": entry["smiles"],
            "solvent_smiles": "O",
            "solute_charge": charge,
            "solute_env": sorted(
                env.name for env in get_checkmol_environments(entry["smiles"])
            ),
            "solvent_env": sorted(env.name for env in get_checkmol_environments("O")),
        }

    return data, skipped_systems


def find_overlapping_solutes(
    mnsol_data: dict[str, dict], freesolv_data: dict[str, dict]
) -> tuple[set[str], set[str]]:
    """Find overlapping solutes across MNSol and FreeSolv using OpenFF isomorphism."""
    mnsol_solute_to_smiles = {
        record["solute_name"]: record["solute_smiles"] for record in mnsol_data.values()
    }
    freesolv_solute_to_smiles = {
        record["solute_name"]: record["solute_smiles"]
        for record in freesolv_data.values()
    }

    mnsol_mols = {
        name: make_openff_molecule(smiles)
        for name, smiles in mnsol_solute_to_smiles.items()
    }
    freesolv_mols = {
        name: make_openff_molecule(smiles)
        for name, smiles in freesolv_solute_to_smiles.items()
    }

    freesolv_by_canonical = defaultdict(list)
    for name, smiles in freesolv_solute_to_smiles.items():
        canonical = canonicalize_smiles(smiles)
        if canonical is not None:
            freesolv_by_canonical[canonical].append(name)

    overlapping_mnsol_solutes: set[str] = set()
    overlapping_freesolv_solutes: set[str] = set()

    for mnsol_name, mnsol_smiles in mnsol_solute_to_smiles.items():
        mnsol_mol = mnsol_mols.get(mnsol_name)
        if mnsol_mol is None:
            continue

        canonical = canonicalize_smiles(mnsol_smiles)
        candidate_freesolv_names = []
        if canonical is not None:
            candidate_freesolv_names = freesolv_by_canonical.get(canonical, [])

        if not candidate_freesolv_names:
            candidate_freesolv_names = list(freesolv_solute_to_smiles)

        for freesolv_name in candidate_freesolv_names:
            freesolv_mol = freesolv_mols.get(freesolv_name)
            if freesolv_mol is None:
                continue

            if molecules_are_isomorphic(mnsol_mol, freesolv_mol):
                overlapping_mnsol_solutes.add(mnsol_name)
                overlapping_freesolv_solutes.add(freesolv_name)

    return overlapping_mnsol_solutes, overlapping_freesolv_solutes


def find_multi_solvent_mnsol_solutes(mnsol_data: dict[str, dict]) -> set[str]:
    """Return MNSol solutes that appear with more than one solvent."""
    solvents_by_solute = defaultdict(set)
    for record in mnsol_data.values():
        solvents_by_solute[record["solute_name"]].add(record["solvent_name"])

    return {
        solute_name
        for solute_name, solvents in solvents_by_solute.items()
        if len(solvents) > 1
    }


def select_overlap_mnsol_seed_keys(
    mnsol_data: dict[str, dict],
    overlapping_mnsol_solutes: set[str],
    max_rows_per_solute: int,
) -> set[str]:
    """Select up to ``max_rows_per_solute`` overlap rows per MNSol solute.

    For each overlapping solute, rows are selected to maximize solvent
    diversity using solvent fingerprint distance.
    """
    ordered_keys = sorted(mnsol_data, key=lambda k: int(k.split("-")[1]))
    order_map = {key: i for i, key in enumerate(ordered_keys)}
    solvent_fp_by_key = {
        key: solute_fingerprint(mnsol_data[key]["solvent_smiles"])
        for key in ordered_keys
    }

    keys_by_solute = defaultdict(list)
    for key in ordered_keys:
        solute_name = mnsol_data[key]["solute_name"]
        if solute_name in overlapping_mnsol_solutes:
            keys_by_solute[solute_name].append(key)

    selected_seed_keys: set[str] = set()

    for solute_name in sorted(keys_by_solute):
        candidate_keys = keys_by_solute[solute_name]
        if len(candidate_keys) <= max_rows_per_solute:
            selected_seed_keys.update(candidate_keys)
            continue

        selected_for_solute = select_diverse_keys(
            candidate_keys,
            set(),
            solvent_fp_by_key,
            order_map,
            max_rows_per_solute,
        )
        selected_seed_keys.update(selected_for_solute)

    return selected_seed_keys


@click.command()
@click.option(
    "--mnsol-alldata",
    type=click.Path(exists=True, dir_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to MNSol_alldata.txt",
)
@click.option(
    "--n-per-environment",
    default=2,
    show_default=True,
    type=int,
    help="Target count per chemical environment when filling coverage.",
)
def main(mnsol_alldata: pathlib.Path, n_per_environment: int):
    """Create aligned FreeSolv/MNSol subsets using overlap-only seeds,
    then fill chemical-environment coverage.

    The seed selection follows the spirit of curate-fsolv-test-set.py by using
    only MNSol/FreeSolv overlapping solutes. Additional structures are then
    added to improve CHEMICAL_ENVIRONMENTS coverage.
    """

    mnsol_filtered_path = (
        PACKAGE_ROOT
        / "openfe_benchmarks/data/benchmark_systems/solvation_set/mnsol_neutral/subset_openff_filtered.json"
    )
    mnsol_subset_path = (
        PACKAGE_ROOT
        / "openfe_benchmarks/data/benchmark_systems/solvation_set/mnsol_neutral/subset_openff_small.json"
    )
    freesolv_filtered_path = (
        PACKAGE_ROOT
        / "openfe_benchmarks/data/benchmark_systems/solvation_set/freesolv/subset_openff_filtered.json"
    )
    freesolv_subset_path = (
        PACKAGE_ROOT
        / "openfe_benchmarks/data/benchmark_systems/solvation_set/freesolv/subset_openff_small.json"
    )

    mnsol_data, mnsol_skips = build_filtered_mnsol_records(mnsol_alldata)
    freesolv_data, freesolv_skips = build_filtered_freesolv_records()

    overlapping_mnsol_solutes, overlapping_freesolv_solutes = find_overlapping_solutes(
        mnsol_data, freesolv_data
    )

    mnsol_seed_keys = select_overlap_mnsol_seed_keys(
        mnsol_data,
        overlapping_mnsol_solutes,
        max_rows_per_solute=2,
    )
    freesolv_seed_keys = {
        key
        for key, record in freesolv_data.items()
        if record["solute_name"] in overlapping_freesolv_solutes
    }

    mnsol_subset, mnsol_added_for_coverage, mnsol_missing = fill_environment_coverage(
        mnsol_data,
        mnsol_seed_keys,
        n_per_environment=n_per_environment,
        max_rows_per_solute=2,
        environment_mode="solvent_solute_pairs",
    )
    freesolv_subset, freesolv_added_for_coverage, freesolv_missing = (
        fill_environment_coverage(
            freesolv_data,
            freesolv_seed_keys,
            n_per_environment=n_per_environment,
            environment_mode="solute_only",
        )
    )

    print(f"Filtered MNSol systems: {len(mnsol_data)}")
    print(f"Filtered FreeSolv systems: {len(freesolv_data)}")
    print(f"Overlapping MNSol solutes: {len(overlapping_mnsol_solutes)}")
    print(f"Overlapping FreeSolv solutes: {len(overlapping_freesolv_solutes)}")

    print("MNSol skip summary:")
    for reason, systems in sorted(mnsol_skips.items()):
        print(f"    {reason}: {len(systems)}")

    print("FreeSolv skip summary:")
    for reason, systems in sorted(freesolv_skips.items()):
        print(f"    {reason}: {len(systems)}")

    print(
        "MNSol subset: "
        f"seed={len(mnsol_seed_keys)} total={len(mnsol_subset)} "
        f"coverage_additions={len(mnsol_added_for_coverage)} "
        f"missing_env={len(mnsol_missing)}"
    )
    print(
        "FreeSolv subset: "
        f"seed={len(freesolv_seed_keys)} total={len(freesolv_subset)} "
        f"coverage_additions={len(freesolv_added_for_coverage)} "
        f"missing_env={len(freesolv_missing)}"
    )

    if mnsol_missing:
        print(f"MNSol missing environments: {', '.join(mnsol_missing)}")
    if freesolv_missing:
        print(f"FreeSolv missing environments: {', '.join(freesolv_missing)}")

    with open(str(mnsol_filtered_path), "w", encoding="utf-8") as f:
        json.dump(mnsol_data, f, cls=JSON_HANDLER.encoder, indent=4)
    with open(str(mnsol_subset_path), "w", encoding="utf-8") as f:
        json.dump(mnsol_subset, f, cls=JSON_HANDLER.encoder, indent=4)

    with open(str(freesolv_filtered_path), "w", encoding="utf-8") as f:
        json.dump(freesolv_data, f, cls=JSON_HANDLER.encoder, indent=4)
    with open(str(freesolv_subset_path), "w", encoding="utf-8") as f:
        json.dump(freesolv_subset, f, cls=JSON_HANDLER.encoder, indent=4)


if __name__ == "__main__":
    main()
