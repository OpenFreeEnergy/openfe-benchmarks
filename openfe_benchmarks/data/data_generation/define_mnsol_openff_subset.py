import json
import csv
from collections import defaultdict

from pathlib import Path
import click
import pathlib
from tqdm import tqdm
from rdkit import RDLogger
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

from gufe.tokenization import JSON_HANDLER
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


def _down_select_representative_subset(
    data: dict[str, dict], n_per_environment: int = 2
) -> tuple[dict[str, dict], list[str]]:
    """Select up to ``n_per_environment`` maximally diverse solutes per
    (solvent, chemical-environment pair), using Morgan fingerprint Tanimoto
    distance for diversity.

    Returns the down-selected data dict and a list of removed keys.
    """
    target_environment_names = [env.name for env in CHEMICAL_ENVIRONMENTS]

    def build_target_environment_pairs() -> list[tuple[str, str]]:
        pairs = [(env, env) for env in target_environment_names]
        for i, env_i in enumerate(target_environment_names):
            for env_j in target_environment_names[i + 1 :]:
                pairs.append((env_i, env_j))
        return pairs

    target_environment_pairs = build_target_environment_pairs()

    def matches_environment_pair(
        solvent_envs: set[str], solute_envs: set[str], environment_pair: tuple[str, str]
    ) -> bool:
        env_a, env_b = environment_pair
        if env_a == env_b:
            return env_a in solvent_envs and env_a in solute_envs
        return (env_a in solvent_envs and env_b in solute_envs) or (
            env_b in solvent_envs and env_a in solute_envs
        )

    def solute_fingerprint(solute_name: str):
        smiles = NAMES_TO_SMILES.get(solute_name)
        if smiles is None:
            return None
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=True)
        except Exception:
            return None
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

    def parse_key_index(key: str) -> int:
        return int(key.split("-")[1])

    def select_diverse_keys(
        candidate_keys: list[str],
        already_selected_keys: set[str],
        fp_by_key: dict[str, object],
        n_to_select: int,
    ) -> list[str]:
        available = [k for k in candidate_keys if k not in already_selected_keys]
        chosen: list[str] = []

        while available and len(chosen) < n_to_select:
            reference_keys = [*already_selected_keys, *chosen]

            if not reference_keys:
                best_key = min(available, key=parse_key_index)
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
                    key=lambda k: (min_distance_to_references(k), -parse_key_index(k)),
                )

            chosen.append(best_key)
            available.remove(best_key)

        return chosen

    keys_by_solvent = defaultdict(list)
    for key, record in data.items():
        keys_by_solvent[record["solvent_name"]].append(key)

    selected_keys: set[str] = set()

    for _, solvent_keys in keys_by_solvent.items():
        ordered_keys = sorted(solvent_keys, key=lambda k: int(k.split("-")[1]))
        key_to_envs: dict[str, set[str]] = {}
        fp_by_key: dict[str, object] = {}
        for key in ordered_keys:
            key_to_envs[key] = set(data[key].get("solute_env", []))
            fp_by_key[key] = solute_fingerprint(data[key]["solute_name"])

        solvent_envs = (
            set(data[ordered_keys[0]].get("solvent_env", [])) if ordered_keys else set()
        )
        selected_for_solvent: set[str] = set()

        for environment_pair in target_environment_pairs:
            candidate_keys = [
                key
                for key in ordered_keys
                if key not in selected_for_solvent
                and matches_environment_pair(
                    solvent_envs, key_to_envs[key], environment_pair
                )
            ]

            selected_for_pair = select_diverse_keys(
                candidate_keys,
                selected_for_solvent,
                fp_by_key,
                n_per_environment,
            )
            selected_for_solvent.update(selected_for_pair)

        selected_keys.update(selected_for_solvent)

    selected_data = {key: data[key] for key in sorted(selected_keys)}
    removed_keys = [key for key in sorted(data) if key not in selected_data]
    return selected_data, removed_keys


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
    # Unassigned tetrahedral stereocenters
    chiral_centers = Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    )
    if any(tag == "?" for _, tag in chiral_centers):
        return True
    # Unassigned stereochemistry reported by RDKit (e.g., double bonds)
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
    """Return the set of ``ChemicalEnvironment`` values detected by
    ``yammbs.checkmol`` for the given SMILES. Returns an empty set if the
    SMILES is invalid or no checkmol detection function is found.
    """
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


@click.command()
@click.option(
    "--mnsol-alldata",
    type=click.Path(exists=True, dir_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to MNSol_alldata.txt",
)
def main(mnsol_alldata: pathlib.Path):
    """Filter MNSol_alldata.txt and write ``subset_openff.json`` for the
    ``mnsol_neutral`` benchmark set, keyed by ``mnsol-XXXX``.

    Rows are excluded if the solute or solvent is water, absent from
    ``mnsol-name-to-smiles.json``, contains 'radical', is an explicit-solvent
    entry (Level1==14), is charged, is a self-solvation pair, matches a
    disqualifying SMIRKS (long chains, 1,3-dicarbonyls, or dissociating molecules, i.e., HBr or HCl), has
    undefined stereochemistry, or contains elements outside
    ``ALLOWED_ELEMENTS``. The passing set is then down-selected to 2 diverse
    entries per (solvent, chemical-environment pair). Skip counts are printed
    to stdout.
    """

    subset_filename_filtered = (
        PACKAGE_ROOT
        / "openfe_benchmarks/data/benchmark_systems/solvation_set/mnsol_neutral/subset_openff_filtered.json"
    )
    subset_filename_small = (
        PACKAGE_ROOT
        / "openfe_benchmarks/data/benchmark_systems/solvation_set/mnsol_neutral/subset_openff_small.json"
    )
    data = {}
    skip_molecules = ["water"]

    skipped_systems = defaultdict(list)
    with mnsol_alldata.open("r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        total_lines = sum(1 for _ in fh)
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter="\t")
        for row in tqdm(reader, total=total_lines):
            if not row or not row.get("FileHandle"):
                continue

            index = f"{int(row.get('No.')):04d}"
            key = f"mnsol-{index}"
            solute_name = row.get("SoluteName").strip().strip('"')
            solvent_name = row.get("Solvent").strip().strip('"')
            group = row.get("Level1").strip().strip('"')
            charge = row.get("Charge").strip().strip('"')

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

            checks = [
                (
                    contains_any_filter_smirks(NAMES_TO_SMILES[solute_name]),
                    "solute has disqualifying smirks",
                ),
                (
                    contains_any_filter_smirks(NAMES_TO_SMILES[solvent_name]),
                    "solvent has disqualifying smirks",
                ),
                (
                    has_undefined_stereochemistry(NAMES_TO_SMILES[solvent_name]),
                    "solvent has undefined stereochemistry",
                ),
                (
                    has_undefined_stereochemistry(NAMES_TO_SMILES[solute_name]),
                    "solute has undefined stereochemistry",
                ),
                (
                    not has_only_allowed_elements(NAMES_TO_SMILES[solvent_name]),
                    "solvent elements outside chemical space",
                ),
                (
                    not has_only_allowed_elements(NAMES_TO_SMILES[solute_name]),
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
                "solute_charge": charge,
                "solute_env": sorted(
                    env.name
                    for env in get_checkmol_environments(NAMES_TO_SMILES[solute_name])
                ),
                "solvent_env": sorted(
                    env.name
                    for env in get_checkmol_environments(NAMES_TO_SMILES[solvent_name])
                ),
            }

    # Down Selection: representative subset per solvent (n_per_environment=2).
    data_subset, removed_keys = _down_select_representative_subset(
        data, n_per_environment=2
    )
    for key in removed_keys:
        skipped_systems["down selection: representative subset"].append(
            key.split("-")[1]
        )

    print(f"There were skipped {sum(len(x) for x in skipped_systems.values())} systems")
    for reason, systems in skipped_systems.items():
        print(f"    {reason}: {len(systems)}")
    print(f"Down-selected filtered systems from {len(data)} to {len(data_subset)}")
    with open(str(subset_filename_filtered), "w") as f:
        json.dump(data, f, cls=JSON_HANDLER.encoder, indent=4)
    with open(str(subset_filename_small), "w") as f:
        json.dump(data_subset, f, cls=JSON_HANDLER.encoder, indent=4)


if __name__ == "__main__":
    main()
