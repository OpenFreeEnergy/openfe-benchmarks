"""Build aligned MNSol and FreeSolv benchmark subsets with Tanimoto-guided selection.

Output format
-------------
This script produces two output files per dataset, each containing structural identifiers:
  - ``subset_openff_filtered.json``  — the full filtered pool (all valid records)
  - ``subset_openff_small.json``     — the curated subset chosen by the algorithm below

Each record is keyed by system identifier (e.g. "mobley_1017962,water") and contains only
solute_* and solvent_* structural identifiers:
  - ``solute_smiles``, ``solute_inchi``, ``solute_inchikey``, ``solute_name``, ``solute_iupac``
  - ``solvent_smiles``, ``solvent_inchi``, ``solvent_inchikey``, ``solvent_name``

Excluded (looked up from original experimental JSON instead):
  - Experimental data: dG, reference, uncertainty, notes
  - Metadata: solute_charge
  - Runtime metadata: solute_env, solvent_env

No experimental data is replicated in subset files. Original experimental records can be
looked up from the source JSON files if needed.

MNSol subset selection
----------------------
Phase 1 — FreeSolv-overlap seeding
    Solutes confirmed to overlap with FreeSolv (via OpenFF SMILES isomorphism)
    are guaranteed inclusion. Up to ``max_rows_per_solute`` rows per solute are
    chosen by greedy Tanimoto MaxMin on concatenated 4096-bit pair fingerprints
    to maximise solvent diversity.

Phase 2 — checkmol environment coverage
    For each checkmol ``ChemicalEnvironment`` still below ``n_per_environment``
    representatives, additional rows are added from the full filtered pool via
    greedy Tanimoto MaxMin.

Phase 3 — uniform Tanimoto fill
    Remaining budget (``max_subset_size`` minus current size) is filled by
    greedy MaxMin over concatenated 4096-bit (solute ‖ solvent) Morgan
    fingerprints, maximizing pair-space coverage.

Phase 4 — gap rebalancing
    Up to ``max_swaps`` swap iterations move points from dense Tanimoto clusters
    into sparse regions. Each swap removes the selected point closest to any
    other selected neighbour and replaces it with the unselected point that
    covers the largest remaining gap, subject to:
      - per-environment hard floor of ``n_per_environment - 1``
      - per-solute cap of ``max_rows_per_solute``
    Swaps are scored by (cap_penalty, env_scarcity, gap_distance). A tabu set
    prevents oscillation.

FreeSolv subset selection
-------------------------
FreeSolv has only water as solvent, so environment coverage uses solute
environments only. MNSol-overlapping solutes are seeded first; remaining budget
is filled by ``fill_environment_coverage`` targeting ``n_per_environment`` per
environment.

Differences from openff-sage subsets
-------------------------------------
The openff-sage <= 2.30 MNSol subset does not use Tanimoto-guided selection or gap
rebalancing; its FreeSolv subset is primarily an MNSol overlap-membership filter
without an independent chemical-space fill stage.
"""

import json
import pathlib
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any

import click
from gufe.tokenization import JSON_HANDLER
from openff.toolkit import Molecule
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem
import numpy as np
from yammbs import checkmol
from yammbs.checkmol import ChemicalEnvironment

import openfe_benchmarks

RDLogger.DisableLog("rdApp.*")
# Collect noisy third-party warnings and print them once at the end.
_warning_log: list[warnings.WarningMessage] = []

PACKAGE_ROOT = Path(openfe_benchmarks.__file__).parent.parent
MNSOL_SOLVENTS_PER_SOLUTE = 3

DEFAULT_MNSOL_JSON_PATH = (
    PACKAGE_ROOT
    / "openfe_benchmarks/data/benchmark_systems/solvation_set/mnsol_neutral/experimental_solvation_free_energy_data.json"
)
DEFAULT_FREESOLV_JSON_PATH = (
    PACKAGE_ROOT
    / "openfe_benchmarks/data/benchmark_systems/solvation_set/freesolv/experimental_solvation_free_energy_data.json"
)

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

# ChemicalEnvironment values outside OpenFF chemical space
# (require elements absent from ALLOWED_ELEMENTS) and therefore excluded:
#   Fluorine (F):    AlkylFluoride, ArylFluoride, AcylFluoride
#   Iodine (I):      AlkylIodide, ArylIodide, AcylIodide
#   Boron (B):       BoronicAcidDeriv, BoronicAcid, BoronicAcidEster
#   Organometallics: Organometallic, Organolithium, Organomagnesium
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
    ChemicalEnvironment.CarboxylicAcid,
    ChemicalEnvironment.CarboxylicAcidEster,
    ChemicalEnvironment.Ether,
    ChemicalEnvironment.Aldehyde,
    ChemicalEnvironment.Ketone,
    ChemicalEnvironment.Imine,
    ChemicalEnvironment.Aromatic,
    ChemicalEnvironment.CarboxylicAcidPrimaryAmide,
    ChemicalEnvironment.CarboxylicAcidSecondaryAmide,
    ChemicalEnvironment.CarboxylicAcidTertiaryAmide,
    ChemicalEnvironment.PrimaryAmine,
    ChemicalEnvironment.SecondaryAmine,
    ChemicalEnvironment.TertiaryAmine,
    ChemicalEnvironment.Nitrile,
    ChemicalEnvironment.Cyanate,
    ChemicalEnvironment.Isocyanate,
    ChemicalEnvironment.NitroCompound,
    ChemicalEnvironment.Heterocycle,
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
    ChemicalEnvironment.Sulfoxide,
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
    """Return an RDKit canonical SMILES for matching."""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        return None

    if mol is None:
        return None

    return Chem.MolToSmiles(mol)


def make_openff_molecule(smiles: str) -> Molecule | None:
    """Build an OpenFF Molecule for isomorphism checks."""
    try:
        return Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    except Exception:
        return None


def to_cmiles(smiles: str) -> str | None:
    """Return CMILES for a SMILES string."""

    try:
        rdmol = Chem.MolFromSmiles(smiles)
        if rdmol is None:
            return None
        rdmol = Chem.AddHs(rdmol)
        return Chem.MolToSmiles(rdmol, canonical=True)
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


def pair_fingerprint(solute_smiles: str, solvent_smiles: str):
    """Return a concatenated solute/solvent Morgan fingerprint for a system pair.

    Solute bits occupy positions 0–2047 and solvent bits occupy positions
    2048–4095 of a 4096-bit vector. Concatenation preserves solute/solvent
    role information (unlike OR-combination, where (solute A + solvent B) and
    (solute B + solvent A) are indistinguishable). This is the standard
    approach for multi-component systems (reaction fingerprints, mixture QSAR).
    """
    fp_solute = solute_fingerprint(solute_smiles)
    fp_solvent = solute_fingerprint(solvent_smiles)
    if fp_solute is None or fp_solvent is None:
        return None
    combined = DataStructs.ExplicitBitVect(4096)
    for bit in fp_solute.GetOnBits():
        combined.SetBit(bit)  # positions 0–2047
    for bit in fp_solvent.GetOnBits():
        combined.SetBit(2048 + bit)  # positions 2048–4095
    return combined


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

    below_target_environments = [
        env_name
        for env_name in target_environment_names
        if count_for_environment(env_name) < n_per_environment
    ]

    selected_data = {key: data[key] for key in sorted(selected)}
    return selected_data, sorted(coverage_additions), below_target_environments


def build_filtered_records_from_experimental_json(
    reference_json: pathlib.Path,
    skip_inchikey: set[str] | None = None,
) -> tuple[dict[str, dict], dict[str, list[str]]]:
    """Load validated experimental solvation reference data from a JSON file."""
    data: dict[str, dict] = {}
    skipped_systems = defaultdict(list)
    skip_inchikey = skip_inchikey or set()

    with reference_json.open("r", encoding="utf-8") as fh:
        raw = json.load(fh)

    if not isinstance(raw, dict):
        raise ValueError(
            f"Experimental reference JSON {reference_json} must contain a top-level object"
        )

    for key, record in raw.items():
        if not isinstance(record, dict):
            skipped_systems["invalid record"].append(str(key))
            continue

        if skip_inchikey:
            solute_inchikey = record.get("solute_inchikey")
            solvent_inchikey = record.get("solvent_inchikey")
            if solute_inchikey in skip_inchikey or solvent_inchikey in skip_inchikey:
                skipped_systems["skip inchikey"].append(str(key))
                continue

        missing_fields = [
            field
            for field in (
                "solute_name",
                "solvent_name",
                "solute_smiles",
                "solvent_smiles",
            )
            if not isinstance(record.get(field), str)
        ]
        if missing_fields:
            skipped_systems["missing fields"].append(str(key))
            continue

        solute_smiles = record["solute_smiles"]
        solvent_smiles = record["solvent_smiles"]

        if record["solute_inchikey"] == record["solvent_inchikey"]:
            skipped_systems["solute and solvent are the same"].append(str(key))
            continue

        smiles_checks = [
            (
                contains_any_filter_smirks(solute_smiles),
                "solute has disqualifying smirks",
            ),
            (
                contains_any_filter_smirks(solvent_smiles),
                "solvent has disqualifying smirks",
            ),
            (
                has_undefined_stereochemistry(solute_smiles),
                "solute has undefined stereochemistry",
            ),
            (
                has_undefined_stereochemistry(solvent_smiles),
                "solvent has undefined stereochemistry",
            ),
            (
                not has_only_allowed_elements(solute_smiles),
                "solute elements outside chemical space",
            ),
            (
                not has_only_allowed_elements(solvent_smiles),
                "solvent elements outside chemical space",
            ),
        ]
        skip_reason = next(
            (reason for condition, reason in smiles_checks if condition), None
        )
        if skip_reason:
            skipped_systems[skip_reason].append(str(key))
            continue

        data[str(key)] = {
            **record,
            "solute_name": record["solute_name"],
            "solvent_name": record["solvent_name"],
            "solute_smiles": solute_smiles,
            "solvent_smiles": solvent_smiles,
            "solute_charge": record.get("solute_charge", 0),
        }

    return data, skipped_systems


def tag_records_with_checkmol_environments(data: dict[str, dict]) -> None:
    """Annotate each record in *data* with checkmol environment lists (in place).

    Adds ``solute_env`` and ``solvent_env`` keys to each record: lists of
    ``ChemicalEnvironment.name`` strings drawn from ``CHEMICAL_ENVIRONMENTS``.
    Deduplicates checkmol calls across records sharing the same SMILES.
    """
    target_env_names = {env.name for env in CHEMICAL_ENVIRONMENTS}
    smiles_to_envs: dict[str, list[str]] = {}

    def _compute(smiles: str) -> list[str]:
        if smiles not in smiles_to_envs:
            envs = get_checkmol_environments(smiles)
            smiles_to_envs[smiles] = [
                e.name for e in envs if e.name in target_env_names
            ]
        return smiles_to_envs[smiles]

    for record in data.values():
        record["solute_env"] = _compute(record["solute_smiles"])
        record["solvent_env"] = _compute(record["solvent_smiles"])


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


def rebalance_tanimoto_coverage(
    selected: set[str],
    all_valid_keys: list[str],
    fp_by_key: dict[str, Any],
    key_to_envs: dict[str, set[str]],
    key_to_solute: dict[str, str],
    n_per_environment: int,
    max_rows_per_solute: int,
    max_swaps: int = 500,
) -> set[str]:
    """Swap over-represented points for under-represented ones to fill Tanimoto gaps.

    After the initial selection phases there can be residual dense clusters
    (e.g. many similar alkane/solvent pairs grouped together) leaving sparse
    regions of pair chemical space uncovered.  This function addresses that
    by iterating up to *max_swaps* times:

    1. Precompute the full pairwise Tanimoto distance matrix for all valid keys
       (feasible at ~2000 records).  Gap distances and crowdedness scores are
       derived from this matrix rather than from a 2D UMAP projection, giving
       exact Tanimoto distances consistent with the selection metric in Phases
       1–3.

    2. Build the list of the 50 unselected candidates with the largest minimum
       Tanimoto distance to any selected point (the biggest gap candidates).

    3. For each gap candidate p_in (largest gap first), score every selected
       point as a potential swap-out using a weighted penalty:

         penalty(env) = 10^(n_per_env - n_actual + 1)

       where n_actual is the current count of that environment in the selection.
       This gives near-zero penalty when an environment is well above the target
       and grows steeply as n_actual approaches n_per_env.  An absolute hard
       floor prevents any environment from dropping below n_per_env - 1.
       Ties are broken by nearest-neighbour Tanimoto distance (most crowded =
       smallest nn_dist preferred for removal).

    4. Cap handling: the per-solute cap is treated as a penalty rather than a
       hard gate on the swap-out search.  For each (p_in, p_out) pair the
       combined score is (cap_penalty, env_scarcity, nn_dist):

         cap_penalty = 0  if p_out is same solute as p_in (count-neutral), or
                            p_in's solute is still below cap (no violation)
         cap_penalty = 1  if adding p_in would exceed cap by exactly one
                            (cross-solute fallback, used only when
                            same-solute and cap-neutral options are worse)
         hard-blocked     if adding p_in would exceed cap by more than one

       This means same-solute swaps are always preferred, but cross-solute
       swaps that bring a new (solute, solvent) pair to a sparse region are
       permitted at the cost of temporarily putting one solute at cap+1.

    5. If no valid swap-out exists for the current gap candidate, continue to
       the next gap candidate rather than stopping the whole phase.

    6. Stall detection: if _MAX_STALL consecutive outer iterations all fail to
       find any swap, stop early (the remaining gaps cannot be filled under the
       current constraints).

    The set size is unchanged; only its distribution in pair chemical space
    improves.
    """
    selected = set(selected)  # work on a copy
    target_env_names = [env.name for env in CHEMICAL_ENVIRONMENTS]

    env_counts: dict[str, int] = {
        env: sum(1 for k in selected if env in key_to_envs.get(k, set()))
        for env in target_env_names
    }
    solute_counts: defaultdict[str, int] = defaultdict(int)
    for k in selected:
        solute_counts[key_to_solute[k]] += 1

    unselected = [k for k in all_valid_keys if k not in selected]

    # Precompute pairwise Tanimoto distance matrix for all valid keys.
    # Gap and crowdedness scores are read from this matrix each swap iteration.
    # At ~2000 records this is ~16 MB (float32) — feasible in memory.
    all_fps = [fp_by_key[k] for k in all_valid_keys]
    key_to_idx = {k: i for i, k in enumerate(all_valid_keys)}
    n_all = len(all_valid_keys)
    sim_matrix = np.zeros((n_all, n_all), dtype=np.float32)
    for i, fp in enumerate(all_fps):
        sim_matrix[i] = DataStructs.BulkTanimotoSimilarity(fp, all_fps)
    dist_matrix = 1.0 - sim_matrix
    np.fill_diagonal(dist_matrix, np.inf)  # self-distance excluded from nn/gap lookups

    def _swap_out_penalty(cand_key: str, p_in_envs: set[str]) -> float:
        """Scarcity-weighted penalty for removing cand_key from the selection.

        Returns the maximum scarcity score over all tracked environments covered
        by cand_key:

            scarcity(env) = n_per_environment / n_actual   (in [0, 1])

        Approaches 1 when an environment is barely at target; near 0 when
        abundant.  This correctly dominates the sort: rare-env-covering points
        are last to be removed.

        Returns inf if the swap would drop any environment below the hard floor
        of n_per_environment - 1, absolutely blocking that candidate.
        """
        max_pen = 0.0
        for env in key_to_envs.get(cand_key, set()):
            if env not in env_counts:
                continue
            n_actual = env_counts[env]
            if n_actual == 0:
                continue
            n_after = n_actual - 1 + (1 if env in p_in_envs else 0)
            if n_after < n_per_environment - 1:
                return float("inf")  # hard floor: never drop below n-1
            pen = n_per_environment / n_actual  # in (0, 1]; higher = scarcer
            if pen > max_pen:
                max_pen = pen
        return max_pen

    # Report gap statistics before rebalancing.
    sel_init_idx = [key_to_idx[k] for k in selected]
    unsel_init_idx = [key_to_idx[k] for k in unselected] if unselected else None
    if unsel_init_idx is not None:
        d_init = dist_matrix[np.ix_(unsel_init_idx, sel_init_idx)].min(axis=1)
        print(
            f"  Phase 4 start: max_gap={d_init.max():.3f}  "
            f"mean_gap={d_init.mean():.3f}  "
            f"p90_gap={np.percentile(d_init, 90):.3f}"
        )

    swaps_made = 0
    same_solute_swaps = 0
    cross_solute_swaps = 0
    blocked_all_inf = 0  # gap candidates where all swap-out options are env-blocked
    swap_gaps: list[float] = []  # gap distance at the moment of each successful swap
    swap_env_pens: list[float] = []  # env_penalty of the chosen swap-out
    swap_cap_pens: list[float] = []  # cap_penalty (0 or 1) of the chosen swap-out
    stall_count = 0
    _MAX_STALL = 15  # stop when this many consecutive outer iterations all fail
    # Tabu: recently-swapped-out keys cannot be swapped back in for _TABU_TTL
    # iterations, preventing pointless oscillation between equivalent structures.
    _TABU_TTL = 5
    tabu: dict[str, int] = {}  # key -> swap_iter at which it expires

    for swap_iter in range(max_swaps):
        if not unselected:
            break

        sel_list = list(selected)
        sel_idx = [key_to_idx[k] for k in sel_list]
        unsel_idx = [key_to_idx[k] for k in unselected]

        # Expire tabu entries whose TTL has passed.
        tabu = {k: exp for k, exp in tabu.items() if exp > swap_iter}
        tabu_set = set(tabu)

        # Gap: min Tanimoto distance from each unselected point to any selected point.
        gap_dists = dist_matrix[np.ix_(unsel_idx, sel_idx)].min(axis=1)
        for ti, tk in enumerate(unselected):
            if tk in tabu_set:
                gap_dists[ti] = 0.0  # suppress tabu keys from gap ranking

        # NN distance for selected points: smallest = most crowded = best to remove.
        sel_pair_dists = dist_matrix[np.ix_(sel_idx, sel_idx)]
        nn_dist = sel_pair_dists.min(axis=1)

        swap_done = False

        # Try gap candidates from largest to smallest (capped at 50 for speed).
        for gap_rank, gap_idx in enumerate(np.argsort(-gap_dists)[:50]):
            if gap_dists[gap_idx] < 1e-6:
                break  # no real gap left

            p_in = unselected[gap_idx]
            p_in_envs = key_to_envs.get(p_in, set())
            p_in_solute = key_to_solute[p_in]
            p_in_count = solute_counts[p_in_solute]

            # Score each candidate: (cap_penalty, env_penalty, nn_dist).
            # cap_penalty=0 → same-solute (neutral) or p_in not at cap (no violation).
            # cap_penalty=1 → cross-solute swap that would put p_in's solute at cap+1.
            # cap_penalty=inf → would exceed cap+1 (hard block, never allowed).
            # env_penalty = scarcity of the env being removed (lower = easier to remove).
            # nn_dist = crowdedness tie-breaker (smaller = more redundant in Tanimoto pair space).
            scored: list[tuple[float, float, float, int]] = []
            n_inf = 0
            worst_env_blocking: str = ""
            for i, k in enumerate(sel_list):
                k_solute = key_to_solute[k]
                if k_solute == p_in_solute:
                    cap_pen = 0.0  # same-solute: count unchanged
                elif p_in_count < max_rows_per_solute:
                    cap_pen = 0.0  # p_in not yet at cap: cross-solute fine
                elif p_in_count == max_rows_per_solute:
                    cap_pen = 1.0  # would go to cap+1: allowed as fallback
                else:
                    continue  # already over cap+1: hard block

                env_pen = _swap_out_penalty(k, p_in_envs)
                if env_pen < float("inf"):
                    scored.append((cap_pen, env_pen, nn_dist[i], i))
                else:
                    n_inf += 1
                    if not worst_env_blocking:
                        for env in key_to_envs.get(k, set()):
                            if env not in env_counts:
                                continue
                            n_actual = env_counts[env]
                            n_after = n_actual - 1 + (1 if env in p_in_envs else 0)
                            if n_after < n_per_environment - 1:
                                worst_env_blocking = f"{env}(count={n_actual},floor={n_per_environment - 1})"
                                break

            if not scored:
                blocked_all_inf += 1
                continue

            scored.sort()
            best_cap_pen, best_env_pen, best_nn, best_out_idx = scored[0]
            p_out = sel_list[best_out_idx]
            p_out_envs = key_to_envs.get(p_out, set())
            p_out_solute = key_to_solute[p_out]
            is_same_solute = p_out_solute == p_in_solute

            # Accumulate per-swap statistics.
            swap_gaps.append(float(gap_dists[gap_idx]))
            swap_env_pens.append(best_env_pen)
            swap_cap_pens.append(best_cap_pen)

            # Execute swap.
            selected.discard(p_out)
            selected.add(p_in)
            unselected.pop(gap_idx)
            unselected.append(p_out)

            for env in p_out_envs:
                if env in env_counts:
                    env_counts[env] -= 1
            for env in p_in_envs:
                if env in env_counts:
                    env_counts[env] += 1
            if is_same_solute:
                same_solute_swaps += 1
            else:
                cross_solute_swaps += 1
                solute_counts[p_out_solute] -= 1
                solute_counts[p_in_solute] += 1

            swaps_made += 1
            stall_count = 0
            swap_done = True
            tabu[p_out] = (
                swap_iter + _TABU_TTL
            )  # block p_out from re-entry for TTL iters
            break  # recompute gaps with updated selection

        if not swap_done:
            stall_count += 1
            if stall_count >= _MAX_STALL:
                print(
                    f"  Phase 4: stalled after {_MAX_STALL} consecutive no-swap iterations."
                )
                break

    # Report gap statistics after rebalancing.
    if unselected:
        sel_final_idx = [key_to_idx[k] for k in selected]
        unsel_final_idx = [key_to_idx[k] for k in unselected]
        d_final = dist_matrix[np.ix_(unsel_final_idx, sel_final_idx)].min(axis=1)
        print(
            f"  Phase 4 end:   max_gap={d_final.max():.3f}  "
            f"mean_gap={d_final.mean():.3f}  "
            f"p90_gap={np.percentile(d_final, 90):.3f}"
        )

    _a = np.array(swap_gaps) if swap_gaps else np.array([0.0])
    _e = np.array(swap_env_pens) if swap_env_pens else np.array([0.0])
    cap1_swaps = int(sum(1 for c in swap_cap_pens if c > 0))
    print(
        f"  Phase 4: {swaps_made}/{max_swaps} swaps "
        f"(same-solute={same_solute_swaps}, cross-solute={cross_solute_swaps}, "
        f"cap+1={cap1_swaps}) | "
        f"skipped: all_env_blocked={blocked_all_inf}\n"
        f"    gap at swap:  min={_a.min():.3f}  mean={_a.mean():.3f}  "
        f"max={_a.max():.3f}  p90={np.percentile(_a, 90):.3f}\n"
        f"    env_pen:      min={_e.min():.3f}  mean={_e.mean():.3f}  "
        f"max={_e.max():.3f}  p90={np.percentile(_e, 90):.3f}"
    )
    if swaps_made > 0 and cross_solute_swaps == 0:
        # Every swap replaced one solvent partner for a solute already in the
        # subset.  This typically means the unselected pool is saturated at the
        # per-solute cap, leaving no room to introduce genuinely new solutes.
        # Adjustable parameters that change this outcome:
        #   --max-rows-per-solute N   raise the per-solute cap (currently
        #                             {max_rows_per_solute}); a higher value
        #                             allows cross-solute swaps but reduces
        #                             solute diversity.
        #   --max-subset-size N       increase the subset budget so more novel
        #                             solutes can be added in Phases 2/3 before
        #                             Phase 4 is reached.
        #   --max-swaps N             more swap iterations (currently
        #                             {max_swaps}) give Phase 4 more chances to
        #                             find rare cross-solute opportunities.
        print(
            f"  Phase 4 note: all {swaps_made} swaps were same-solute "
            f"(solvent-partner changes only). This usually means the unselected pool "
            f"is saturated at the per-solute cap={max_rows_per_solute}. "
            f"To allow cross-solute swaps, consider:\n"
            f"    --max-rows-per-solute  raise from {max_rows_per_solute} "
            f"(e.g. 4 or 5)\n"
            f"    --max-subset-size      raise from current value to add more\n"
            f"                           unique solutes in earlier phases"
        )
    return selected


def build_mnsol_subset(
    data: dict[str, dict],
    overlapping_mnsol_solutes: set[str],
    n_per_environment: int,
    max_size: int,
    max_rows_per_solute: int,
    max_swaps: int = 500,
) -> tuple[dict[str, dict], list[str]]:
    """Build MNSol subset with uniform Tanimoto pair-space coverage.

    Selection proceeds in four phases, all using greedy Tanimoto MaxMin on
    concatenated 4096-bit (solute ‖ solvent) Morgan fingerprints.
    SAGE data is never consulted.

    Phase 1 — Seed (FreeSolv alignment)
        For each FreeSolv-overlapping solute, add up to *max_rows_per_solute*
        rows choosing solvents by Tanimoto MaxMin so diverse solvents are
        included first (needed for transfer free energy estimates).

    Phase 2 — Environment coverage
        For each checkmol target environment with fewer than *n_per_environment*
        representatives, add the most Tanimoto-distant eligible candidate that
        respects the per-solute cap.

    Phase 3 — Uniform fill
        Add records up to *max_size* using Tanimoto MaxMin, still respecting
        the per-solute cap so multi-solvent solutes are not over-represented.

    Phase 4 — Rebalance
        Swap over-represented points from dense clusters for unselected points
        in sparse regions, without reducing any environment count below
        *n_per_environment* or violating the per-solute cap.

    Returns
    -------
    selected_data : dict
        Subset records keyed by MNSol key.
    below_target_environments : list[str]
        Target environment names that have fewer than *n_per_environment*
        representatives after selection. This includes both environments with
        zero representatives (absent) and environments where all available
        representatives were selected but the pool had fewer than
        *n_per_environment* (pool-limited).
    """
    # --- Build pair fingerprints ---
    ordered_keys = sorted(data, key=lambda k: int(k.split("-")[1]))
    order_map = {key: i for i, key in enumerate(ordered_keys)}

    valid_keys: list[str] = []
    pair_fps: dict[str, Any] = {}
    for key in ordered_keys:
        fp = pair_fingerprint(data[key]["solute_smiles"], data[key]["solvent_smiles"])
        if fp is None:
            continue
        valid_keys.append(key)
        pair_fps[key] = fp

    if not valid_keys:
        return {}, [env.name for env in CHEMICAL_ENVIRONMENTS]

    key_to_solute = {key: data[key]["solute_name"] for key in valid_keys}
    # Phases 2–4 track both solute and solvent environments so that the
    # coverage fill and swap penalty apply to the full pair chemistry.
    # FreeSolv uses solute-only tracking (fill_environment_coverage) because
    # FreeSolv has only water as a solvent and solvent coverage is not meaningful.
    key_to_envs = {
        key: set(data[key].get("solute_env", []))
        | set(data[key].get("solvent_env", []))
        for key in valid_keys
    }

    selected: set[str] = set()
    solute_counts: defaultdict[str, int] = defaultdict(int)

    def can_add(key: str) -> bool:
        return (
            key not in selected
            and solute_counts[key_to_solute[key]] < max_rows_per_solute
        )

    # --- Phase 1: Seed with FreeSolv-overlapping solutes ---
    keys_by_overlap_solute: dict[str, list[str]] = defaultdict(list)
    for key in valid_keys:
        if key_to_solute[key] in overlapping_mnsol_solutes:
            keys_by_overlap_solute[key_to_solute[key]].append(key)

    print(
        f"  Phase 1: seeding {len(keys_by_overlap_solute)} FreeSolv-overlap solutes..."
    )
    for solute_name in sorted(keys_by_overlap_solute):
        if len(selected) >= max_size:
            break
        candidates = [k for k in keys_by_overlap_solute[solute_name] if can_add(k)]
        n = min(max_rows_per_solute, max_size - len(selected))
        picked = select_diverse_keys(candidates, selected, pair_fps, order_map, n)
        for k in picked:
            selected.add(k)
            solute_counts[key_to_solute[k]] += 1

    # --- Phase 2: Fill environment coverage ---
    target_env_names = [env.name for env in CHEMICAL_ENVIRONMENTS]

    def env_count(env_name: str) -> int:
        return sum(1 for k in selected if env_name in key_to_envs.get(k, set()))

    print(f"  Phase 2: filling environment coverage ({n_per_environment} per env)...")
    coverage_count_before = len(selected)
    for env_name in target_env_names:
        if len(selected) >= max_size:
            break
        needed = n_per_environment - env_count(env_name)
        if needed <= 0:
            continue
        candidates = [
            k
            for k in valid_keys
            if env_name in key_to_envs.get(k, set()) and can_add(k)
        ]
        n = min(needed, max_size - len(selected))
        picked = select_diverse_keys(candidates, selected, pair_fps, order_map, n)
        for k in picked:
            selected.add(k)
            solute_counts[key_to_solute[k]] += 1

    coverage_additions = len(selected) - coverage_count_before

    # --- Phase 3: Tanimoto MaxMin fill up to max_size ---
    if len(selected) < max_size:
        candidates = [k for k in valid_keys if can_add(k)]
        remaining = max_size - len(selected)
        print(f"  Phase 3: Tanimoto MaxMin fill ({remaining} slots remaining)...")
        picked = select_diverse_keys(
            candidates, selected, pair_fps, order_map, remaining
        )
        for k in picked:
            selected.add(k)
            solute_counts[key_to_solute[k]] += 1

    # --- Phase 4: Rebalance gaps in Tanimoto pair space ---
    print("  Phase 4: rebalancing Tanimoto pair-space coverage...")
    selected = rebalance_tanimoto_coverage(
        selected,
        valid_keys,
        pair_fps,
        key_to_envs,
        key_to_solute,
        n_per_environment=n_per_environment,
        max_rows_per_solute=max_rows_per_solute,
        max_swaps=max_swaps,
    )

    below_target_environments = [
        env_name
        for env_name in target_env_names
        if sum(1 for k in selected if env_name in key_to_envs.get(k, set()))
        < n_per_environment
    ]

    print(
        f"\nMNSol subset: total={len(selected)} "
        f"coverage_additions={coverage_additions} below_target_env={len(below_target_environments)}"
    )

    _pool_env_counts = {
        env_name: sum(1 for k in valid_keys if env_name in key_to_envs.get(k, set()))
        for env_name in below_target_environments
    }
    _mnsol_absent = sorted(
        e for e in below_target_environments if _pool_env_counts[e] == 0
    )
    _mnsol_pool_limited = sorted(
        e
        for e in below_target_environments
        if 0 < _pool_env_counts[e] < n_per_environment
    )
    _mnsol_budget_limited = sorted(
        e for e in below_target_environments if _pool_env_counts[e] >= n_per_environment
    )
    print(
        f"  below-target env breakdown — absent from pool: {len(_mnsol_absent)}  "
        f"pool-limited (all available selected, but <{n_per_environment} in pool): {len(_mnsol_pool_limited)}  "
        f"budget-limited (>={n_per_environment} in pool, not all selected): {len(_mnsol_budget_limited)}"
    )
    if _mnsol_absent:
        print(f"    absent: {', '.join(_mnsol_absent)}")
    if _mnsol_pool_limited:
        print(f"    pool-limited: {', '.join(_mnsol_pool_limited)}")
    if _mnsol_budget_limited:
        print(f"    budget-limited: {', '.join(_mnsol_budget_limited)}")

    selected_data = {key: data[key] for key in sorted(selected)}
    return selected_data, below_target_environments


@click.command()
@click.option(
    "--mnsol-json",
    type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path),
    default=DEFAULT_MNSOL_JSON_PATH,
    show_default=True,
    help="Path to MNSol experimental_solvation_free_energy_data.json",
)
@click.option(
    "--freesolv-json",
    type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path),
    default=DEFAULT_FREESOLV_JSON_PATH,
    show_default=True,
    help="Path to FreeSolv experimental_solvation_free_energy_data.json",
)
@click.option(
    "--n-per-environment",
    default=3,
    show_default=True,
    type=int,
    help="Target count per checkmol environment (phases 2 of MNSol selection and FreeSolv fill).",
)
@click.option(
    "--max-subset-size",
    default=400,
    show_default=True,
    type=int,
    help="Hard ceiling on MNSol subset size.",
)
@click.option(
    "--max-rows-per-solute",
    default=MNSOL_SOLVENTS_PER_SOLUTE,
    show_default=True,
    type=int,
    help=(
        "Maximum rows per unique solute in the MNSol subset. "
        "Keeping this >= 2 preserves multi-solvent pairs for transfer free energy."
    ),
)
@click.option(
    "--max-swaps",
    default=500,
    show_default=True,
    type=int,
    help=(
        "Maximum number of swap iterations in Phase 4. Each swap moves one point "
        "from a dense cluster to the largest remaining gap in Tanimoto pair space, "
        "subject to environment-coverage and per-solute-cap constraints."
    ),
)
def main(
    mnsol_json: pathlib.Path,
    freesolv_json: pathlib.Path,
    n_per_environment: int,
    max_subset_size: int,
    max_rows_per_solute: int,
    max_swaps: int,
):
    """Create aligned FreeSolv/MNSol subsets.

    Writes four JSON files:
      subset_openff_filtered.json  (MNSol full filtered pool)
      subset_openff_small.json     (MNSol curated subset)
      subset_openff_filtered.json  (FreeSolv full filtered pool)
      subset_openff_small.json     (FreeSolv curated subset)

    MNSol selection — four phases, SAGE data never consulted:
      Phase 1 — seed FreeSolv-overlapping solutes (Tanimoto MaxMin over solvents).
      Phase 2 — fill checkmol environment coverage to n_per_environment (Tanimoto MaxMin).
      Phase 3 — fill remaining budget with uniform Tanimoto MaxMin over pair space.
      Phase 4 — swap dense-cluster points into Tanimoto gaps (up to max_swaps iterations).

    FreeSolv selection seeds from MNSol-overlapping solutes then fills
    environment coverage independently (solute environments only; solvent is always water).
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

    if not pathlib.Path(mnsol_json).exists():
        raise FileNotFoundError(
            f"MNSol data file not found: {mnsol_json}\n"
            "This file is not included in the repository because MNSol is a "
            "licensed database. Licensed MNSol users can regenerate it with:\n"
            "  openfe_benchmarks/data/data_generation/generate_mnsol_data.py"
        )

    mnsol_data, mnsol_skips = build_filtered_records_from_experimental_json(
        mnsol_json, skip_inchikey=["XLYOFNOQVPJJNP-UHFFFAOYNA-N"]
    )
    freesolv_data, freesolv_skips = build_filtered_records_from_experimental_json(
        freesolv_json,
    )

    overlapping_mnsol_solutes, overlapping_freesolv_solutes = find_overlapping_solutes(
        mnsol_data, freesolv_data
    )

    mnsol_unique_solutes = {rec["solute_name"] for rec in mnsol_data.values()}
    freesolv_unique_solutes = {rec["solute_name"] for rec in freesolv_data.values()}

    print("=== Pool summary ===")
    print(
        f"MNSol:    {len(mnsol_data)} systems, {len(mnsol_unique_solutes)} unique solutes"
    )
    print(
        f"FreeSolv: {len(freesolv_data)} systems, {len(freesolv_unique_solutes)} unique solutes"
    )
    print(
        f"Solute overlap (SMILES isomorphism): "
        f"{len(overlapping_mnsol_solutes)} MNSol / "
        f"{len(overlapping_freesolv_solutes)} FreeSolv"
    )

    print("\nMNSol skip summary:")
    for reason, systems in sorted(mnsol_skips.items()):
        print(f"    {reason}: {len(systems)}")
    print("FreeSolv skip summary:")
    for reason, systems in sorted(freesolv_skips.items()):
        print(f"    {reason}: {len(systems)}")

    print("\nTagging records with checkmol environments...")
    tag_records_with_checkmol_environments(mnsol_data)
    tag_records_with_checkmol_environments(freesolv_data)

    print("\nBuilding MNSol subset (Tanimoto-guided, no SAGE data)...")
    mnsol_subset, mnsol_below_target = build_mnsol_subset(
        mnsol_data,
        overlapping_mnsol_solutes,
        n_per_environment=n_per_environment,
        max_size=max_subset_size,
        max_rows_per_solute=max_rows_per_solute,
        max_swaps=max_swaps,
    )

    freesolv_seed_keys = {
        key
        for key, record in freesolv_data.items()
        if record["solute_name"] in overlapping_freesolv_solutes
    }
    freesolv_subset, freesolv_added_for_coverage, freesolv_below_target = (
        fill_environment_coverage(
            freesolv_data,
            freesolv_seed_keys,
            n_per_environment=n_per_environment,
            environment_mode="solute_only",
        )
    )

    print(
        f"FreeSolv subset: total={len(freesolv_subset)} "
        f"coverage_additions={len(freesolv_added_for_coverage)} "
        f"below_target_env={len(freesolv_below_target)}\n"
    )
    if mnsol_below_target:
        print(
            f"MNSol environments below target (<{n_per_environment} representatives): {', '.join(mnsol_below_target)}"
        )
    if freesolv_below_target:
        print(
            f"FreeSolv environments below target (<{n_per_environment} representatives): {', '.join(freesolv_below_target)}"
        )

    _freesolv_key_to_envs = {
        key: set(freesolv_data[key].get("solute_env", [])) for key in freesolv_data
    }
    _freesolv_pool_env_counts = {
        env_name: sum(1 for k in freesolv_data if env_name in _freesolv_key_to_envs[k])
        for env_name in freesolv_below_target
    }
    _freesolv_absent = sorted(
        e for e in freesolv_below_target if _freesolv_pool_env_counts[e] == 0
    )
    _freesolv_pool_limited = sorted(
        e
        for e in freesolv_below_target
        if 0 < _freesolv_pool_env_counts[e] < n_per_environment
    )
    _freesolv_budget_limited = sorted(
        e
        for e in freesolv_below_target
        if _freesolv_pool_env_counts[e] >= n_per_environment
    )
    print(
        f"FreeSolv env breakdown — absent from pool: {len(_freesolv_absent)}  "
        f"pool-limited (all available selected, but <{n_per_environment} in pool): {len(_freesolv_pool_limited)}  "
        f"budget-limited (>={n_per_environment} in pool, not all selected): {len(_freesolv_budget_limited)}"
    )
    if _freesolv_absent:
        print(f"  absent: {', '.join(_freesolv_absent)}")
    if _freesolv_pool_limited:
        print(f"  pool-limited: {', '.join(_freesolv_pool_limited)}")
    if _freesolv_budget_limited:
        print(f"  budget-limited: {', '.join(_freesolv_budget_limited)}")

    # --- Final statistics ---
    # Pool row counts per solute (from the full filtered pool)
    pool_solute_counts: dict[str, int] = {}
    for rec in mnsol_data.values():
        sname = rec["solute_name"]
        pool_solute_counts[sname] = pool_solute_counts.get(sname, 0) + 1

    # MNSol solute repetition in subset
    mnsol_subset_solute_counts: dict[str, int] = {}
    for rec in mnsol_subset.values():
        sname = rec["solute_name"]
        mnsol_subset_solute_counts[sname] = mnsol_subset_solute_counts.get(sname, 0) + 1

    n_unique_mnsol = len(mnsol_subset_solute_counts)
    n_repeated = sum(1 for c in mnsol_subset_solute_counts.values() if c > 1)
    repeat_counts = sorted(mnsol_subset_solute_counts.values(), reverse=True)
    max_repeats = repeat_counts[0] if repeat_counts else 0
    dist: dict[int, int] = {}
    for c in repeat_counts:
        dist[c] = dist.get(c, 0) + 1

    # For each selection-count bucket, classify solutes by why they weren't selected more:
    #   data-limited  → selected == pool rows (no more rows exist)
    #   at cap        → selected == max_rows_per_solute (cap prevented more; pool has additional rows)
    #   selection-limited → selected < pool rows AND selected < cap (more available but not chosen)
    bucket_breakdown: dict[int, dict[str, int]] = {}
    for sname, sel_count in mnsol_subset_solute_counts.items():
        pool_count = pool_solute_counts.get(sname, sel_count)
        if sel_count not in bucket_breakdown:
            bucket_breakdown[sel_count] = {
                "data_limited": 0,
                "at_cap": 0,
                "selection_limited": 0,
            }
        if sel_count >= pool_count:
            bucket_breakdown[sel_count]["data_limited"] += 1
        elif sel_count >= max_rows_per_solute:
            bucket_breakdown[sel_count]["at_cap"] += 1
        else:
            bucket_breakdown[sel_count]["selection_limited"] += 1

    # MNSol/FreeSolv solute overlap in the final subsets
    # Use the SMILES-isomorphism-verified overlap sets (name spaces differ between datasets)
    mnsol_subset_solutes = {rec["solute_name"] for rec in mnsol_subset.values()}
    freesolv_subset_solutes = {rec["solute_name"] for rec in freesolv_subset.values()}
    # overlapping_mnsol_solutes: MNSol-name keys confirmed to match a FreeSolv molecule
    # overlapping_freesolv_solutes: FreeSolv-ID keys confirmed to match an MNSol molecule
    overlap_in_mnsol_subset = mnsol_subset_solutes & overlapping_mnsol_solutes
    overlap_in_freesolv_subset = freesolv_subset_solutes & overlapping_freesolv_solutes

    print("\n=== Final statistics ===")
    print(
        f"MNSol subset:   {len(mnsol_subset)} systems, "
        f"{n_unique_mnsol} unique solutes, "
        f"{n_repeated} solutes appear more than once (max {max_repeats}x)"
    )
    print(f"  Solute repetition (cap={max_rows_per_solute}):")
    for sel_count, n_solutes in sorted(dist.items()):
        bd = bucket_breakdown.get(sel_count, {})
        data_lim = bd.get("data_limited", 0)
        at_cap = bd.get("at_cap", 0)
        sel_lim = bd.get("selection_limited", 0)
        parts = []
        if data_lim:
            parts.append(f"{data_lim} data-limited (pool exhausted)")
        if at_cap:
            parts.append(f"{at_cap} at cap (pool has more)")
        if sel_lim:
            parts.append(f"{sel_lim} selection-limited (pool has more, under cap)")
        breakdown_str = f"  [{', '.join(parts)}]" if parts else ""
        print(f"    {sel_count}x: {n_solutes} solutes{breakdown_str}")

    print(
        f"FreeSolv subset: {len(freesolv_subset)} systems, "
        f"{len(freesolv_subset_solutes)} unique solutes"
    )
    print(
        f"Solute overlap (MNSol subset ∩ FreeSolv pool): "
        f"{len(overlap_in_mnsol_subset)} / {n_unique_mnsol} MNSol solutes "
        f"({100 * len(overlap_in_mnsol_subset) / n_unique_mnsol:.1f}%)"
    )
    print(
        f"Solute overlap (FreeSolv subset ∩ MNSol pool): "
        f"{len(overlap_in_freesolv_subset)} / {len(freesolv_subset_solutes)} FreeSolv solutes "
        f"({100 * len(overlap_in_freesolv_subset) / max(len(freesolv_subset_solutes), 1):.1f}%)"
    )
    if overlap_in_mnsol_subset:
        print(
            f"  MNSol overlapping solutes: {', '.join(sorted(overlap_in_mnsol_subset))}"
        )

    def _extract_system_keys(data: dict[str, dict]) -> dict[str, dict]:
        """Extract only system-identifying fields (solute_* and solvent_* keys, excluding runtime metadata).

        This keeps the subset files focused on structural identifiers (SMILES, InChI, InChIKey, names, IUPAC)
        while avoiding duplication of experimental data (dG, reference, uncertainty, etc.)
        and runtime metadata (solute_env, solvent_env, solute_charge) that can be looked up
        from the original experimental JSON files.
        """
        excluded = {"solute_env", "solvent_env", "solute_charge"}
        return {
            k: {
                f: v
                for f, v in rec.items()
                if (f.startswith("solute_") or f.startswith("solvent_"))
                and f not in excluded
            }
            for k, rec in data.items()
        }

    with open(str(mnsol_filtered_path), "w", encoding="utf-8") as f:
        json.dump(
            _extract_system_keys(mnsol_data), f, cls=JSON_HANDLER.encoder, indent=4
        )
    with open(str(mnsol_subset_path), "w", encoding="utf-8") as f:
        json.dump(
            _extract_system_keys(mnsol_subset), f, cls=JSON_HANDLER.encoder, indent=4
        )

    with open(str(freesolv_filtered_path), "w", encoding="utf-8") as f:
        json.dump(
            _extract_system_keys(freesolv_data), f, cls=JSON_HANDLER.encoder, indent=4
        )
    with open(str(freesolv_subset_path), "w", encoding="utf-8") as f:
        json.dump(
            _extract_system_keys(freesolv_subset), f, cls=JSON_HANDLER.encoder, indent=4
        )


if __name__ == "__main__":
    with warnings.catch_warnings(record=True) as _warning_log:
        warnings.simplefilter("always")
        main(standalone_mode=False)
    if _warning_log:
        print(f"\n--- {len(_warning_log)} suppressed warning(s) ---")
        seen: set[str] = set()
        for w in _warning_log:
            key = f"{w.category.__name__}: {w.message}"
            if key not in seen:
                seen.add(key)
                print(
                    f"  [{w.category.__name__}] {w.message} ({w.filename}:{w.lineno})"
                )
