#!/usr/bin/env python3
"""Prepare submission metadata from OpenFE archives using BenchmarkResults."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date
import glob as glob_module
import json
import logging
from pathlib import Path
import re
import sys
import textwrap
from typing import cast
import yaml

from pint import Quantity

from gufe import (
    AlchemicalNetwork,
    ProteinComponent,
    SmallMoleculeComponent,
    SolventComponent,
)
from gufe.archival import AlchemicalArchive
from gufe.transformations.transformation import Transformation

from openfe_benchmarks.data import BenchmarkIndex
from openfe_benchmarks.results import BenchmarkResults
from openfe_benchmarks.results._benchmark_results import Archive, LiteralStr
from openfe_benchmarks.scripts.utils import load_archive, load_alchemical_network
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openfe.protocols.openmm_septop import SepTopProtocol
from openfe.protocols.openmm_afe import AbsoluteSolvationProtocol
from pontibus.protocols.relative import HybridTopProtocol
from pontibus.protocols.solvation import ASFEProtocol

logger = logging.getLogger(__name__)

# map the protocol to a calculation type for the purposes of metadata aggregation and summary
_PROTOCOL_MAPPING = {
    RelativeHybridTopologyProtocol: "rbfe",
    HybridTopProtocol: "rbfe",
    SepTopProtocol: "rbfe",
    ASFEProtocol: "rbfe",
    AbsoluteSolvationProtocol: "rbfe",
}


@dataclass(frozen=True)
class _ModeSpec:
    key: str
    detect_prefixes: tuple[str, ...]
    detect_default: bool
    use_mapping_ligands: bool
    rbfe_like: bool
    ligand_count_min: int
    ligand_count_max: int
    summary_mode_label: str
    summary_sentence_builder: Callable[[int, int, int], str]
    simulation_setting_keys: tuple[tuple[str, str, str], ...]


def _rbfe_summary_sentence(
    n_transformations: int, n_ligands: int, _n_solvents: int
) -> str:
    return f"The submission contains {n_transformations} edges across {n_ligands} unique ligands."


def _asfe_summary_sentence(
    n_transformations: int, n_ligands: int, n_solvents: int
) -> str:
    return (
        f"The submission contains {n_transformations} edges across {n_ligands} unique solutes and "
        f"{n_solvents} unique solvents."
    )


_MODE_SPECS: dict[str, _ModeSpec] = {
    "rbfe": _ModeSpec(
        key="rbfe",
        detect_prefixes=("complex_", "solvent_"),
        detect_default=False,
        use_mapping_ligands=True,
        rbfe_like=True,
        ligand_count_min=1,
        ligand_count_max=2,
        summary_mode_label="RBFE",
        summary_sentence_builder=_rbfe_summary_sentence,
        simulation_setting_keys=(
            ("simulation_settings", "equilibration_time", "production_time"),
        ),
    ),
    "asfe": _ModeSpec(
        key="asfe",
        detect_prefixes=(),
        detect_default=True,
        use_mapping_ligands=False,
        rbfe_like=False,
        ligand_count_min=1,
        ligand_count_max=1,
        summary_mode_label="ASFE",
        summary_sentence_builder=_asfe_summary_sentence,
        simulation_setting_keys=(
            (
                "vacuum_simulation_settings",
                "vacuum_equilibration_time",
                "vacuum_production_time",
            ),
            (
                "solvent_simulation_settings",
                "solvent_equilibration_time",
                "solvent_production_time",
            ),
        ),
    ),
}


def _mode_spec(mode: str) -> _ModeSpec:
    if mode not in _MODE_SPECS:
        raise ValueError(f"Unsupported calculation mode: {mode}")
    return _MODE_SPECS[mode]


def _as_obj_dict(value: object) -> dict[str, object]:
    if isinstance(value, dict):
        return cast(dict[str, object], value)
    return {}


def _add_str_value_with_keys(
    values: list[tuple[str, list[str]]],
    value: str,
    keys: list[str],
) -> None:
    for current_value, current_keys in values:
        if current_value == value:
            for key in keys:
                if key not in current_keys:
                    current_keys.append(key)
            return
    values.append((value, list(keys)))


def _add_tuple_value_with_keys(
    values: list[tuple[tuple[str, ...], list[str]]],
    value: tuple[str, ...],
    keys: list[str],
) -> None:
    for current_value, current_keys in values:
        if current_value == value:
            for key in keys:
                if key not in current_keys:
                    current_keys.append(key)
            return
    values.append((value, list(keys)))


def _add_protocol_value_with_keys(
    values: list[tuple[dict[str, str | list[str]], list[str]]],
    value: dict[str, str | list[str]],
    keys: list[str],
) -> None:
    for current_value, current_keys in values:
        if current_value == value:
            for key in keys:
                if key not in current_keys:
                    current_keys.append(key)
            return
    values.append((value, list(keys)))


def _quantity_to_text(value: object) -> str:
    if isinstance(value, dict) and "unit" in value:
        quantity = Quantity(value["val"], value["unit"])
        return f"{quantity:#~}"
    return str(value)


def _extract_charge_provenance(
    component_dict: dict[str, object],
) -> dict[str, object] | None:
    """Extract and parse charge_provenance JSON from component molprops when available."""
    molprops_raw = component_dict.get("molprops")
    if not isinstance(molprops_raw, dict):
        return None

    charge_provenance = molprops_raw.get("charge_provenance")
    if not isinstance(charge_provenance, str) or not charge_provenance:
        return None

    try:
        parsed = json.loads(charge_provenance)
    except (json.JSONDecodeError, TypeError):
        return None

    if isinstance(parsed, dict):
        return cast(dict[str, object], parsed)
    return None


def _normalize_charge_method_from_provenance(provenance: dict[str, object]) -> str:
    """Build normalized charge-method tag from parsed charge_provenance metadata."""
    method = str(provenance.get("charge_method", "")).lower().strip()
    if not method:
        return ""

    if "nagl" in method:
        nagl_model = str(provenance.get("nagl_model", "")).strip()
        if nagl_model:
            nagl_model = nagl_model.split("/")[-1].split("\\")[-1]
            return f"nagl_{nagl_model}"
        return "nagl_off"

    if "am1bccelf10" in method or "elf10" in method:
        return "am1bccelf10_oe"

    if "am1bcc" in method:
        backend = str(provenance.get("off_toolkit_backend", "")).lower().strip()
        if backend == "ambertools":
            return "am1bcc_at"
        if backend == "openeye":
            return "am1bcc_oe"
        return "TODO"

    return re.sub(r"[^a-z0-9._-]+", "_", method).strip("_")


def _partial_charge_from_transformation(
    trans: Transformation,
    mode_spec: _ModeSpec,
) -> str:
    """Extract normalized ligand/cofactor partial-charge method(s) from component provenance."""
    methods: set[str] = set()

    for state_key in ("stateA", "stateB"):
        chemical_system = getattr(trans, state_key)
        if not chemical_system:
            continue

        for label, component in chemical_system.components.items():
            component_type_name = type(component).__name__
            is_small_molecule = isinstance(component, SmallMoleculeComponent) or (
                component_type_name == "SmallMoleculeComponent"
            )
            if not is_small_molecule:
                continue

            label_text = str(label).lower()
            if mode_spec.rbfe_like:
                if "solvent" in label_text:
                    continue
            else:
                if "solute" not in label_text and "solvent" in label_text:
                    continue

            component_to_dict = getattr(component, "to_dict", None)
            if not callable(component_to_dict):
                continue

            provenance = _extract_charge_provenance(component_to_dict())
            if provenance is None:
                continue

            method = _normalize_charge_method_from_provenance(provenance)
            if method:
                methods.add(method)

    if not methods:
        return ""

    return "/".join(sorted(methods))


def _looks_like_serialized_forcefield(value: str) -> bool:
    trimmed = value.lstrip()
    return trimmed.startswith("<?xml") or "<SMIRNOFF" in value or "<ForceField" in value


def _normalize_forcefield_label(value: str) -> str:
    if _looks_like_serialized_forcefield(value):
        return ""
    return value.strip()


def _load_network(
    input_path: Path,
) -> tuple[AlchemicalArchive | AlchemicalNetwork, str]:
    try:
        return load_archive(input_path), "alchemicalarchive"
    except Exception:
        try:
            return load_alchemical_network(input_path), "alchemicalnetwork"
        except Exception as exc:
            raise ImportError(
                "Could not import as AlchemicalArchive or AlchemicalNetwork: "
                f"{input_path}"
            ) from exc


def _transformation_refs(
    network_obj: AlchemicalArchive | AlchemicalNetwork,
    network_mode: str,
) -> list[Transformation]:
    if network_mode == "alchemicalarchive":
        archive_obj = cast(AlchemicalArchive, network_obj)
        return [
            cast(Transformation, trans)
            for trans, _ in archive_obj.transformation_results
        ]
    network = cast(AlchemicalNetwork, network_obj)
    return [cast(Transformation, trans) for trans in network.edges]


def _network_key(
    network_obj: AlchemicalArchive | AlchemicalNetwork, network_mode: str
) -> str:
    if network_mode == "alchemicalarchive":
        return str(cast(AlchemicalArchive, network_obj).network.key)
    return str(cast(AlchemicalNetwork, network_obj).key)


def _detect_mode(transformations: list[Transformation]) -> str:
    """
    Use the protocol type to detect the calculation mode (RBFE vs ASFE) for the submission.
    """
    # grab the first one here we are assuming that all transformations in the archive are of the same type
    protocol_cls = transformations[0].protocol.__class__
    mode_key = _PROTOCOL_MAPPING.get(protocol_cls, None)

    if mode_key is None:
        raise ValueError(
            "Unable to detect calculation from transformation protocol. "
            f"Observed protocol class: {protocol_cls}. "
            f"Known classes: {''.join(cls.__name__ for cls in _PROTOCOL_MAPPING.keys())}."
        )
    return mode_key


def _get_mapping_annotations(trans: Transformation) -> dict[str, object]:
    """Return mapping annotations when present; otherwise an empty dict."""
    mapping = trans.mapping
    if mapping is None:
        return {}
    if hasattr(mapping, "annotations"):
        annotations = getattr(mapping, "annotations")
        if isinstance(annotations, dict):
            return annotations
    return {}


def _get_alchemical_ligands(trans: Transformation) -> set[object]:
    """Return alchemical ligand components for RBFE-style mappings when available
    TODO update to look for SMC differences for SepTop style calculations
    """
    mapping = trans.mapping
    if mapping is None:
        return set()
    component_a = getattr(mapping, "componentA", None)
    component_b = getattr(mapping, "componentB", None)
    if component_a is None or component_b is None:
        return set()
    return {component_a, component_b}


def _default_submission_id(network_key: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", network_key.lower()).strip("-")
    return f"{date.today().isoformat()}-{slug}"


def _normalize_submission_date(
    value: date | str | None, submission_id: str | None = None
) -> str:
    if value is not None:
        if isinstance(value, date):
            return value.isoformat()
        try:
            return date.fromisoformat(value).isoformat()
        except ValueError as exc:
            raise ValueError("submission_date must be ISO 8601 YYYY-MM-DD") from exc

    # No explicit date given: fall back to the ISO date embedded in the
    # submission_id, so the two stay consistent by construction.
    if submission_id is not None:
        try:
            return date.fromisoformat(submission_id[:10]).isoformat()
        except ValueError:
            pass

    return date.today().isoformat()


def _dedupe_preserve_order(items: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for item in items:
        value = item.strip()
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def _generate_title(
    mode: str, systems: list[tuple[str, str]], submission_id: str
) -> str:
    n_systems = len(systems)
    groups = sorted({group for group, _ in systems})
    mode_upper = mode.upper()

    if n_systems == 0:
        return f"OpenFE {mode_upper} Benchmark - {submission_id}"

    if len(groups) == 1:
        group = groups[0]
        names = [name for _, name in systems]
        if n_systems <= 3:
            return (
                f"OpenFE {mode_upper} - {group} - {', '.join(names)} - {submission_id}"
            )
        return f"OpenFE {mode_upper} - {group} ({n_systems} systems) - {submission_id}"

    if n_systems <= 3:
        desc = ", ".join(f"{group}/{name}" for group, name in systems)
        return f"OpenFE {mode_upper} - {desc} - {submission_id}"

    unique_names = len({name for _, name in systems})
    return (
        f"OpenFE {mode_upper} - Multi-group Benchmark "
        f"({len(groups)} groups, {unique_names} systems) - {submission_id}"
    )


def _component_name(component_dict: dict[str, object]) -> str:
    molprops = component_dict.get("molprops") or {}
    if isinstance(molprops, dict) and molprops.get("ofe-name"):
        return str(molprops["ofe-name"])
    for key in ("name", "smiles", "solvent_molecule"):
        if component_dict.get(key):
            return str(component_dict[key])
    return "unknown"


def _infer_system_group_name(
    trans: Transformation,
    override_group: str | None,
    override_name: str | None,
) -> tuple[str, str]:
    annotations = _get_mapping_annotations(trans)
    original_group = annotations.get("system_group")
    original_name = annotations.get("system_name")

    if override_group and original_group and override_group != original_group:
        raise ValueError(
            f"Transformation '{trans.name}' annotation system_group='{original_group}' "
            f"conflicts with override '{override_group}'"
        )
    if override_name and original_name and override_name != original_name:
        raise ValueError(
            f"Transformation '{trans.name}' annotation system_name='{original_name}' "
            f"conflicts with override '{override_name}'"
        )

    system_group = override_group or original_group
    system_name = override_name or original_name

    if system_name and not system_group:
        try:
            index = BenchmarkIndex()
            all_systems = index.list_system_names_by_tag()
            matching_groups = [
                group for group, name in all_systems if name == system_name
            ]
            if matching_groups:
                system_group = matching_groups[0]
        except Exception:
            system_group = None

    if not system_group:
        system_group = "TODO"
    if not system_name:
        system_name = "TODO"

    return str(system_group), str(system_name)


def _extract_system_components(
    trans: Transformation,
    mode_spec: _ModeSpec,
) -> dict[str, set[str] | list[str]]:
    solvents: set[str] = set()
    proteins: set[str] = set()
    ligands: list[str] = []
    cofactors: set[str] = set()

    alchemical_ligands: set[object] = set()
    if mode_spec.use_mapping_ligands:
        alchemical_ligands = _get_alchemical_ligands(trans)

    for state_key in ("stateA", "stateB"):
        chemical_system = getattr(trans, state_key)
        if not chemical_system:
            continue
        for label, component in chemical_system.components.items():
            name = _component_name(component.to_dict())
            if isinstance(component, SolventComponent):
                solvents.add(name)
            elif isinstance(component, ProteinComponent):
                proteins.add(name)
            elif isinstance(component, SmallMoleculeComponent):
                if not mode_spec.use_mapping_ligands:
                    if name not in ligands:
                        ligands.append(name)
                elif component in alchemical_ligands:
                    if name not in ligands:
                        ligands.append(name)
                elif "solvent" not in label:
                    cofactors.add(name)

    if not (mode_spec.ligand_count_min <= len(ligands) <= mode_spec.ligand_count_max):
        raise ValueError(
            f"{mode_spec.summary_mode_label} transformation has invalid ligand count: "
            f"{trans.name} -> {ligands}"
        )

    return {
        "solvents": solvents,
        "proteins": proteins,
        "ligands": ligands,
        "cofactors": cofactors,
    }


def _make_edge_key(
    network_key: str,
    system_group: str,
    system_name: str,
    mode_spec: _ModeSpec,
    components: dict[str, set[str] | list[str]],
) -> str:
    ligands = cast(list[str], components["ligands"])
    ligand_start = ligands[0]
    ligand_final = "none"
    if mode_spec.rbfe_like and len(ligands) > 1:
        ligand_final = ligands[1]

    return (
        f"{network_key} {system_group}-{system_name}: "
        f"ligand_start={ligand_start}, ligand_final={ligand_final}, "
        f"solvent={components['solvents'] or 'none'}, "
        f"cofactors={components['cofactors'] or 'none'}, "
        f"protein={components['proteins'] or 'none'}"
    )


def _extract_protocol_settings(
    protocol_obj: object | None,
    mode_spec: _ModeSpec,
) -> dict[str, str | list[str]]:
    if protocol_obj is None:
        return {
            "protocol": "TODO",
            "protocol_library": "TODO",
            "notes": "Protocol settings unavailable in archive.",
        }

    protocol_name = str(type(protocol_obj)).rstrip("'>").split(".")[-1]
    module_name = type(protocol_obj).__module__
    library = module_name.split(".")[0] if module_name else "TODO"

    settings_obj = getattr(protocol_obj, "settings", None)
    model_dump = getattr(settings_obj, "model_dump", None)
    settings = _as_obj_dict(model_dump() if callable(model_dump) else {})
    payload: dict[str, str | list[str]] = {
        "protocol": protocol_name,
        "protocol_library": library,
        "notes": "",
    }

    if not settings:
        payload["notes"] = (
            "Protocol class found, but detailed settings were unavailable."
        )
        return payload

    # Keep a full-settings fingerprint so aggregation does not collapse subtly
    # different protocols into one entry.
    payload["_settings_fingerprint"] = json.dumps(
        settings,
        sort_keys=True,
        default=str,
    )

    integrator = _as_obj_dict(settings.get("integrator_settings") or {})
    thermo = _as_obj_dict(settings.get("thermo_settings") or {})
    lambda_settings = _as_obj_dict(settings.get("lambda_settings") or {})

    if integrator.get("timestep") is not None:
        payload["timestep"] = _quantity_to_text(integrator["timestep"])
    if thermo.get("temperature") is not None:
        payload["temperature"] = _quantity_to_text(thermo["temperature"])
    if thermo.get("pressure") is not None:
        payload["pressure"] = _quantity_to_text(thermo["pressure"])

    payload["lambda_functions"] = str(lambda_settings.get("lambda_functions", ""))
    if lambda_settings.get("lambda_windows") is not None:
        payload["lambda_windows"] = str(lambda_settings["lambda_windows"])

    ff_settings = _as_obj_dict(
        settings.get("forcefield_settings")
        or settings.get("solvent_forcefield_settings")
        or settings.get("vacuum_forcefield_settings")
        or {}
    )
    if ff_settings:
        payload["small_molecule_forcefield"] = _normalize_forcefield_label(
            str(ff_settings.get("small_molecule_forcefield") or "")
        )
        forcefields_raw = ff_settings.get("forcefields")
        forcefields = forcefields_raw if isinstance(forcefields_raw, list) else []
        normalized = []
        for ff in forcefields:
            label = _normalize_forcefield_label(str(ff).split("/")[-1].split(".")[0])
            if label:
                normalized.append(label)
        payload["forcefields"] = sorted(set(normalized))

    payload["partial_charges"] = "TODO"  # partial charges come from BenchmarkData

    for source_key, eq_key, prod_key in mode_spec.simulation_setting_keys:
        sim = _as_obj_dict(settings.get(source_key) or {})
        if sim.get("equilibration_length") is not None:
            payload[eq_key] = _quantity_to_text(sim["equilibration_length"])
        if sim.get("production_length") is not None:
            payload[prod_key] = _quantity_to_text(sim["production_length"])

    return payload


@dataclass
class _SystemRecord:
    system_group: str
    system_name: str
    network_key: str
    ligands: set[str] = field(default_factory=set)
    proteins: set[str] = field(default_factory=set)
    cofactors: set[str] = field(default_factory=set)
    solvents: set[str] = field(default_factory=set)


@dataclass
class _Metadata:
    mode: str
    network_mode: str
    n_transformations: int = 0
    network_keys: list[str] = field(default_factory=list)
    systems: dict[tuple[str, str], _SystemRecord] = field(default_factory=dict)
    system_order: list[tuple[str, str]] = field(default_factory=list)
    openfe_version: list[tuple[str, list[str]]] = field(default_factory=list)
    openmm_version: list[tuple[str, list[str]]] = field(default_factory=list)
    openff_toolkit_version: list[tuple[str, list[str]]] = field(default_factory=list)
    pontibus_version: list[tuple[str, list[str]]] = field(default_factory=list)
    mapper: list[tuple[str, list[str]]] = field(default_factory=list)
    forcefield: list[tuple[tuple[str, ...], list[str]]] = field(default_factory=list)
    small_molecule_forcefield: list[tuple[str, list[str]]] = field(default_factory=list)
    partial_charges: list[tuple[str, list[str]]] = field(default_factory=list)
    protocol_libraries: list[tuple[str, list[str]]] = field(default_factory=list)
    protocol_settings: list[tuple[dict[str, str | list[str]], list[str]]] = field(
        default_factory=list
    )


def _resolve_input_paths(
    input_files: Path | list[Path] | str | None,
    systems: list[tuple[str, str, str | Path]]
    | tuple[tuple[str, str, str | Path], ...]
    | None,
    system_group: str | None,
    system_name: str | None,
) -> tuple[list[Path], dict[Path, tuple[str | None, str | None]]]:
    if systems is not None:
        if input_files is not None:
            raise ValueError("Cannot specify both input_files and systems.")
        if system_group is not None or system_name is not None:
            raise ValueError(
                "When systems is provided, do not pass system_group/system_name; "
                "use per-item system tuples."
            )
        entries: list[tuple[str, str, Path]] = []
        for item in systems:
            if not isinstance(item, (list, tuple)) or len(item) != 3:
                raise ValueError(
                    "Each systems item must be (system_group, system_name, archive_path)."
                )
            group, name, archive_path = item
            entries.append((str(group).strip(), str(name).strip(), Path(archive_path)))
        if not entries:
            raise ValueError("At least one systems entry is required.")
        paths = [entry[2] for entry in entries]
        overrides: dict[Path, tuple[str | None, str | None]] = {
            entry[2].resolve(): (entry[0], entry[1]) for entry in entries
        }
        return paths, overrides

    if input_files is None:
        raise ValueError("At least one input file must be provided")

    if isinstance(input_files, str):
        matched = glob_module.glob(input_files, recursive=True)
        if not matched:
            raise ValueError(f"No files matched glob pattern: {input_files}")
        paths = [Path(path) for path in sorted(matched)]
    elif isinstance(input_files, Path):
        paths = [input_files]
    else:
        paths = list(input_files)

    if not paths:
        raise ValueError("At least one input file must be provided")

    return paths, {}


def _collapse_value_keys(
    values: list[tuple[str, list[str]]] | list[tuple[tuple[str, ...], list[str]]],
    label: str,
    keys_label: str = "edges",
) -> str | list[str] | list[dict[str, str | list[str]]]:
    if not values:
        return "TODO"

    def _normalize(value: str | tuple[str, ...]) -> str | list[str]:
        if isinstance(value, tuple):
            return [str(item) for item in value]
        return str(value)

    if len(values) == 1:
        return _normalize(values[0][0])

    collapsed: list[dict[str, str | list[str]]] = []
    for value, keys in sorted(values, key=lambda item: (len(item[1]), str(item[0]))):
        entry = {label: _normalize(value)}
        if keys:
            sorted_keys = sorted(str(key) for key in keys)
            entry[keys_label] = sorted_keys[:5] + (
                ["etc."] if len(sorted_keys) > 5 else []
            )
        collapsed.append(entry)
    return collapsed


def _flatten_value_keys(
    values: list[tuple[str, list[str]]] | list[tuple[tuple[str, ...], list[str]]],
) -> list[str]:
    """Return sorted unique string values from value-key pairs for summaries."""
    flattened: set[str] = set()
    for value, _ in values:
        if isinstance(value, (list, tuple, set)):
            for item in value:
                text = str(item).strip()
                if text:
                    flattened.add(text)
        else:
            text = str(value).strip()
            if text:
                flattened.add(text)
    return sorted(flattened)


def _build_protocol_payload(
    protocol_settings: list[tuple[dict[str, str | list[str]], list[str]]],
) -> list[dict[str, str | list[str]]]:
    if not protocol_settings:
        return []

    def _parsed_settings_fingerprint(settings: dict[str, str | list[str]]) -> object:
        raw = settings.get("_settings_fingerprint")
        if not isinstance(raw, str) or not raw:
            return {}
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            return {}

    def _format_path(path: list[str]) -> str:
        dotted = ".".join(path)
        if dotted.endswith(".val"):
            return dotted[:-4]
        return dotted

    def _diff_settings(base: object, other: object) -> list[tuple[str, object, object]]:
        diffs: list[tuple[str, object, object]] = []

        def _recurse(path: list[str], left: object, right: object) -> None:
            if type(left) is not type(right):
                diffs.append((_format_path(path), left, right))
                return

            if isinstance(left, dict):
                left_keys = set(left.keys())
                right_keys = set(right.keys()) if isinstance(right, dict) else set()
                for key in sorted(left_keys | right_keys):
                    left_has = key in left
                    right_has = isinstance(right, dict) and key in right
                    if not left_has and right_has:
                        diffs.append(
                            (_format_path(path + [str(key)]), None, right[key])
                        )
                    elif left_has and not right_has:
                        diffs.append((_format_path(path + [str(key)]), left[key], None))
                    elif right_has:
                        _recurse(path + [str(key)], left[key], right[key])
                return

            if isinstance(left, list):
                if not isinstance(right, list) or left != right:
                    diffs.append((_format_path(path), left, right))
                return

            if left != right:
                diffs.append((_format_path(path), left, right))

        _recurse([], base, other)
        return diffs

    output: list[dict[str, str | list[str]]] = []
    ordered = sorted(
        enumerate(protocol_settings),
        key=lambda item: (len(item[1][1]), str(item[1][0]), item[0]),
    )
    primary_index = max(
        range(len(protocol_settings)),
        key=lambda i: (len(protocol_settings[i][1]), -i),
    )
    primary_settings = protocol_settings[primary_index][0]
    primary_full_settings = _parsed_settings_fingerprint(primary_settings)

    traversal_order = [primary_index] + [
        idx for idx, _ in ordered if idx != primary_index
    ]

    for index in traversal_order:
        settings, keys = protocol_settings[index]
        entry = {k: v for k, v in settings.items() if not k.startswith("_")}
        sorted_keys = sorted(str(key) for key in keys)

        notes_lines: list[str] = []
        existing_notes = entry.get("notes")
        if isinstance(existing_notes, str) and existing_notes.strip():
            notes_lines.append(existing_notes.strip())

        if index != primary_index:
            other_full_settings = _parsed_settings_fingerprint(settings)
            diffs = _diff_settings(primary_full_settings, other_full_settings)
            if diffs:
                notes_lines.append("Detailed protocol settings differ:")
                for path, base_value, other_value in diffs:
                    notes_lines.append(f"- {path}: {base_value!r} -> {other_value!r}")

        if len(sorted_keys) == 0:
            notes_lines.append("Applies to all edges")
        else:
            notes_lines.append(f"Applies to {len(sorted_keys)} edges:")
            for edge_key in sorted_keys[:5]:
                notes_lines.append(f"- {edge_key}")
            if len(sorted_keys) > 5:
                notes_lines.append("- etc.")

        entry["notes"] = LiteralStr("\n".join(notes_lines))
        output.append(entry)
    return output


def _build_content_summary(
    metadata: _Metadata,
    mode_spec: _ModeSpec,
    used_alchemiscale: bool,
) -> str:
    forcefields = _flatten_value_keys(metadata.forcefield)
    if not forcefields or forcefields == ["TODO"]:
        ff_text = "an unspecified force field (TODO)"
    elif len(forcefields) == 1:
        ff_text = forcefields[0]
    else:
        ff_text = "/".join(forcefields)

    charge_methods = _flatten_value_keys(metadata.partial_charges)
    if not charge_methods or charge_methods == ["TODO"]:
        charge_text = "an unspecified partial charge method (TODO)"
    elif len(charge_methods) == 1:
        charge_text = charge_methods[0]
    else:
        charge_text = "/".join(charge_methods)

    small_molecule_ffs = _flatten_value_keys(metadata.small_molecule_forcefield)
    if not small_molecule_ffs or small_molecule_ffs == ["TODO"]:
        small_molecule_ff_text = "an unspecified small molecule force field (TODO)"
    elif len(small_molecule_ffs) == 1:
        small_molecule_ff_text = small_molecule_ffs[0]
    else:
        small_molecule_ff_text = "/".join(small_molecule_ffs)

    systems_by_group: dict[str, list[str]] = defaultdict(list)
    for group, name in metadata.system_order:
        systems_by_group[group].append(name)

    group_parts = []
    for group in sorted(systems_by_group):
        system_names = sorted(systems_by_group[group])
        group_parts.append(f"{group}: {', '.join(system_names)}")
    systems_desc = ", ".join(group_parts)

    all_ligands = set()
    all_solvents = set()
    for system in metadata.systems.values():
        all_ligands.update(system.ligands)
        all_solvents.update(system.solvents)

    if mode_spec.rbfe_like:
        summary = (
            f"This submission describes the {mode_spec.summary_mode_label} benchmark covering {systems_desc} "
            f"prepared with {ff_text} for proteins and solvents, and {small_molecule_ff_text} with {charge_text} "
            "for ligands, solutes, and cofactors. "
            f"{mode_spec.summary_sentence_builder(metadata.n_transformations, len(all_ligands), len(all_solvents))}"
        )
    else:
        summary = (
            f"This submission describes the {mode_spec.summary_mode_label} benchmark covering {systems_desc} "
            f"prepared with {ff_text} for solvents and {charge_text} for solutes and cofactors. "
            f"{mode_spec.summary_sentence_builder(metadata.n_transformations, len(all_ligands), len(all_solvents))}"
        )

    if used_alchemiscale:
        summary += " Results are derived from archived Alchemiscale workflow data."

    return summary


def _make_tags(metadata: _Metadata, user_tags: str) -> list[str]:
    tags: list[str] = [metadata.mode, metadata.network_mode]

    for value, _ in metadata.forcefield:
        values = value if isinstance(value, (list, tuple, set)) else [value]
        for ff in values:
            label = str(ff).strip()
            if label:
                tags.append(label)

    for value, _ in metadata.small_molecule_forcefield:
        label = str(value).strip()
        if label:
            tags.append(label)

    for group, name in metadata.system_order:
        tags.extend([group, name])

    for value, _ in metadata.partial_charges:
        label = str(value).strip()
        if label:
            tags.append(label)

    for value, _ in metadata.protocol_libraries:
        label = str(value).strip()
        if label and label != "TODO":
            tags.append(label)

    tags.extend(tag.strip() for tag in user_tags.split(",") if tag.strip())
    return _dedupe_preserve_order(tags)


def _build_benchmark_data(
    systems: dict[tuple[str, str], _SystemRecord],
) -> dict[str, object]:
    out: dict[str, object] = {
        "source_repository": "https://github.com/OpenFreeEnergy/openfe-benchmarks"
    }
    grouped: dict[str, dict[str, str]] = defaultdict(dict)
    for system in systems.values():
        grouped[system.system_group][system.system_name] = system.network_key
    for group in sorted(grouped):
        out[group] = {name: grouped[group][name] for name in sorted(grouped[group])}
    return out


def _make_zenodo_description(
    benchmark_results: BenchmarkResults,
    network_mode: str,
    used_alchemiscale: bool,
) -> str:
    payload = benchmark_results.to_submission_dict()

    if network_mode == "alchemicalarchive":
        source_description = "AlchemicalArchive"
    elif network_mode == "alchemicalnetwork":
        source_description = "AlchemicalNetwork"
    else:
        source_description = "OpenFE archive"

    workflow_text = "OpenFE"
    if used_alchemiscale:
        workflow_text += " and Alchemiscale"

    authors = payload.get("authors", [])
    author_names = [
        entry.get("name", "")
        for entry in authors
        if isinstance(entry, dict) and isinstance(entry.get("name"), str)
    ]

    tags = payload.get("tags", [])
    tag_text = (
        ", ".join(str(tag) for tag in tags) if isinstance(tags, list) else str(tags)
    )

    benchmark_data = payload.get("benchmark_data", {})
    system_lines: list[str] = []
    provenance_lines: list[str] = []
    network_lines: list[str] = []
    if isinstance(benchmark_data, dict):
        source_repo = benchmark_data.get("source_repository")
        if isinstance(source_repo, str) and source_repo.strip():
            provenance_lines.append(f"- source_repository: {source_repo}")
        for group, group_data in sorted(benchmark_data.items()):
            if group == "source_repository":
                continue
            if isinstance(group_data, dict):
                names = sorted(str(name) for name in group_data)
                system_lines.append(f"- {group}: {', '.join(names)}")
                for name in names:
                    key_value = group_data.get(name)
                    network_lines.append(f"- {key_value}: {group}/{name}")

    systems_block = "\n".join(system_lines) if system_lines else "- TODO"
    network_block = "\n".join(network_lines) if network_lines else "- TODO"
    provenance_block = (
        "\n".join(provenance_lines) if provenance_lines else "- source_repository: TODO"
    )

    protocol_settings = payload.get("protocol_settings", [])
    protocol_count = (
        len(protocol_settings) if isinstance(protocol_settings, list) else 0
    )
    if isinstance(protocol_settings, list):
        protocol_yaml_block = yaml.safe_dump(
            {"protocol_settings": protocol_settings},
            sort_keys=False,
            default_flow_style=False,
            allow_unicode=False,
        ).rstrip()
    else:
        protocol_yaml_block = "protocol_settings: []"

    benchmark_data_yaml_block = yaml.safe_dump(
        {"benchmark_data": benchmark_data if isinstance(benchmark_data, dict) else {}},
        sort_keys=False,
        default_flow_style=False,
        allow_unicode=False,
    ).rstrip()

    repository_reference = (
        "https://github.com/OpenFreeEnergy/openfe-benchmarks/tree/main/openfe_benchmarks/results/"
        f"{benchmark_results.submission_id}"
    )

    software_versions_block = "\n".join(
        [
            f"- openfe_version: {payload.get('openfe_version', 'TODO')}",
            f"- openmm_version: {payload.get('openmm_version', 'TODO')}",
            f"- openff_toolkit_version: {payload.get('openff_toolkit_version', 'TODO')}",
            f"- pontibus_version: {payload.get('pontibus_version', 'TODO')}",
        ]
    )

    recommended_descriptors_block = "\n".join(
        [
            f"- partial_charges: {payload.get('partial_charges', 'TODO')}",
            f"- mapper: {payload.get('mapper', 'TODO')}",
            f"- forcefield: {payload.get('forcefield', 'TODO')}",
            f"- small_molecule_forcefield: {payload.get('small_molecule_forcefield', 'TODO')}",
        ]
    )

    lines = [
        f"# {benchmark_results.title}",
        "",
        "## Overview",
        (
            f"{benchmark_results.calculation_type.upper()} benchmark results prepared from "
            f"{source_description} JSON file(s) generated with {workflow_text}."
        ),
        "",
        benchmark_results.summary,
        "",
        "## Submission Snapshot",
        f"- Submission ID: {benchmark_results.submission_id}",
        f"- Date: {benchmark_results.date}",
        f"- Results file: {benchmark_results.results}",
        f"- License: {benchmark_results.license}",
        f"- Authors: {', '.join(author_names) if author_names else 'TODO'}",
        f"- Tags: {tag_text if tag_text else 'TODO'}",
        "",
        "## Systems Covered",
        systems_block,
        "",
        "## Repository Reference",
        repository_reference,
        "",
        "## Software Versions",
        software_versions_block,
        "",
        "## Alchemical Network Keys",
        network_block,
        "",
        "## Recommended Descriptors",
        recommended_descriptors_block,
        "",
        "## BenchmarkData Provenance",
        provenance_block,
        "",
        "```yaml",
        benchmark_data_yaml_block,
        "```",
        "",
        "## Protocol Settings",
        f"- Protocol settings entries: {protocol_count}",
        "",
        "```yaml",
        protocol_yaml_block,
        "```",
    ]

    return "\n".join(lines) + "\n"


def _update_metadata_from_transformation(
    metadata: _Metadata,
    trans: Transformation,
    network_key: str,
    override_group: str | None,
    override_name: str | None,
    mode_spec: _ModeSpec,
) -> None:
    system_group, system_name = _infer_system_group_name(
        trans, override_group, override_name
    )
    system_key = (system_group, system_name)

    if system_key not in metadata.systems:
        metadata.systems[system_key] = _SystemRecord(
            system_group, system_name, network_key
        )
        metadata.system_order.append(system_key)

    system_record = metadata.systems[system_key]
    components = _extract_system_components(trans, mode_spec)
    system_record.solvents.update(components["solvents"])
    system_record.proteins.update(components["proteins"])
    system_record.cofactors.update(components["cofactors"])
    system_record.ligands.update(components["ligands"])

    edge_key = _make_edge_key(
        network_key, system_group, system_name, mode_spec, components
    )

    annotations = _get_mapping_annotations(trans)

    for key, value in annotations.items():
        value_str = str(value)
        if "openfe" in key:
            _add_str_value_with_keys(metadata.openfe_version, value_str, [key])
        if "openmm" in key:
            _add_str_value_with_keys(metadata.openmm_version, value_str, [key])
        if "openff" in key and "toolkit" in key:
            _add_str_value_with_keys(metadata.openff_toolkit_version, value_str, [key])
        if "pontibus" in key:
            _add_str_value_with_keys(metadata.pontibus_version, value_str, [key])

    mapper_settings = annotations.get("mapper_settings")
    mapper_version = annotations.get("mapper_version")
    if mode_spec.rbfe_like and isinstance(mapper_settings, dict) and mapper_version:
        mapper_name = str(mapper_settings.get("__qualname__", "TODO")).split(".")[-1]
        mapping_algorithm = str(mapper_settings.get("_mapping_algorithm", "TODO"))
        mapper_value = f"{mapper_name} {mapper_version} ({mapping_algorithm})"
        _add_str_value_with_keys(metadata.mapper, mapper_value, [edge_key])

    protocol_settings = _extract_protocol_settings(trans.protocol, mode_spec)
    _add_protocol_value_with_keys(
        metadata.protocol_settings, protocol_settings, [edge_key]
    )

    if protocol_settings.get("forcefields"):
        forcefields = protocol_settings["forcefields"]
        if isinstance(forcefields, list):
            _add_tuple_value_with_keys(
                metadata.forcefield, tuple(forcefields), [edge_key]
            )
    if protocol_settings.get("small_molecule_forcefield"):
        small_molecule_ff = protocol_settings["small_molecule_forcefield"]
        if isinstance(small_molecule_ff, str):
            _add_str_value_with_keys(
                metadata.small_molecule_forcefield, small_molecule_ff, [edge_key]
            )
    partial_charges = _partial_charge_from_transformation(trans, mode_spec)
    if not partial_charges:
        protocol_partial_charges = protocol_settings.get("partial_charges")
        if isinstance(protocol_partial_charges, str):
            partial_charges = protocol_partial_charges
    if partial_charges:
        _add_str_value_with_keys(metadata.partial_charges, partial_charges, [edge_key])
    if protocol_settings.get("protocol_library"):
        protocol_library = protocol_settings["protocol_library"]
        if isinstance(protocol_library, str):
            _add_str_value_with_keys(
                metadata.protocol_libraries, protocol_library, [edge_key]
            )


def _collect_metadata(
    input_paths: list[Path],
    system_overrides: dict[Path, tuple[str | None, str | None]],
    system_group: str | None,
    system_name: str | None,
) -> tuple[_Metadata, _ModeSpec]:
    metadata: _Metadata | None = None
    selected_mode_spec: _ModeSpec | None = None

    for path in input_paths:
        resolved = path.resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"Input file not found: {resolved}")

        network_obj, network_mode = _load_network(resolved)
        transformations = _transformation_refs(network_obj, network_mode)
        mode = _detect_mode(transformations)
        mode_spec = _mode_spec(mode)
        network_key = _network_key(network_obj, network_mode)

        if metadata is None:
            metadata = _Metadata(mode=mode, network_mode=network_mode)
            selected_mode_spec = mode_spec
        else:
            if metadata.mode != mode:
                raise ValueError("Mixed ASFE/RBFE input files are not supported")
            if metadata.network_mode != network_mode:
                raise ValueError(
                    "Mixed AlchemicalArchive/AlchemicalNetwork inputs are not supported"
                )

        metadata.network_keys.append(network_key)
        metadata.n_transformations += len(transformations)

        override_group, override_name = system_overrides.get(
            resolved,
            (system_group, system_name),
        )
        for trans in transformations:
            _update_metadata_from_transformation(
                metadata,
                trans,
                network_key,
                override_group,
                override_name,
                mode_spec,
            )

    if metadata is None or selected_mode_spec is None:
        raise ValueError("No metadata could be extracted from input files")
    return metadata, selected_mode_spec


def _apply_overrides(
    metadata: _Metadata,
    forcefields: list[str] | str | None,
    small_molecule_forcefield: str | None,
    openfe_version: str | None,
    openmm_version: str | None,
    openff_toolkit_version: str | None,
    pontibus_version: str | None,
) -> None:
    if forcefields is not None:
        labels = [forcefields] if isinstance(forcefields, str) else list(forcefields)
        cleaned = [_normalize_forcefield_label(str(label)) for label in labels]
        cleaned = [label for label in cleaned if label]
        if cleaned:
            metadata.forcefield = [(tuple(cleaned), ["override"])]

    if small_molecule_forcefield:
        metadata.small_molecule_forcefield = [(small_molecule_forcefield, ["override"])]

    if openfe_version is not None:
        metadata.openfe_version = [(openfe_version, ["override"])]
    if openmm_version is not None:
        metadata.openmm_version = [(openmm_version, ["override"])]
    if openff_toolkit_version is not None:
        metadata.openff_toolkit_version = [(openff_toolkit_version, ["override"])]
    if pontibus_version is not None:
        metadata.pontibus_version = [(pontibus_version, ["override"])]


def process_network(
    input_files: Path | list[Path] | str | None = None,
    systems: list[tuple[str, str, str | Path]]
    | tuple[tuple[str, str, str | Path], ...]
    | None = None,
    output_dir: Path = Path("."),
    submission_id: str | None = None,
    tags: str = "",
    author: list[str] | None = None,
    license: str = "CC-BY-4.0",
    used_alchemiscale: bool = True,
    summary_suffix: str | None = None,
    results_file: str = "computational_results.json.bz2",
    submission_date: date | str | None = None,
    system_group: str | None = None,
    system_name: str | None = None,
    forcefields: list[str] | str | None = None,
    small_molecule_forcefield: str | None = None,
    openfe_version: str | None = None,
    openmm_version: str | None = None,
    openff_toolkit_version: str | None = None,
    pontibus_version: str | None = None,
) -> tuple[Path, Path]:
    """Generate submission metadata artifacts from OpenFE archive inputs."""
    out_dir = output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    results_path = out_dir / results_file
    if not results_path.exists():
        raise FileNotFoundError(
            f"Required file '{results_file}' not found in output directory: {out_dir}"
        )

    input_paths, system_overrides = _resolve_input_paths(
        input_files=input_files,
        systems=systems,
        system_group=system_group,
        system_name=system_name,
    )

    metadata, mode_spec = _collect_metadata(
        input_paths=input_paths,
        system_overrides=system_overrides,
        system_group=system_group,
        system_name=system_name,
    )

    _apply_overrides(
        metadata,
        forcefields=forcefields,
        small_molecule_forcefield=small_molecule_forcefield,
        openfe_version=openfe_version,
        openmm_version=openmm_version,
        openff_toolkit_version=openff_toolkit_version,
        pontibus_version=pontibus_version,
    )

    summary = _build_content_summary(
        metadata,
        mode_spec=mode_spec,
        used_alchemiscale=used_alchemiscale,
    )
    if summary_suffix:
        summary = summary.rstrip() + " " + summary_suffix.strip()

    systems_for_title = list(metadata.system_order)
    final_submission_id = submission_id or _default_submission_id(
        "_".join(metadata.network_keys)
    )
    title = _generate_title(metadata.mode, systems_for_title, final_submission_id)
    final_tags = _make_tags(metadata, user_tags=tags)

    mapper_value: str | None = None
    small_molecule_ff_value: str | None = None
    if mode_spec.rbfe_like:
        mapper_raw = _collapse_value_keys(metadata.mapper, "mapper")
        if isinstance(mapper_raw, str):
            mapper_value = mapper_raw

        small_molecule_ff_raw = _collapse_value_keys(
            metadata.small_molecule_forcefield,
            "small_molecule_forcefield",
        )
        if isinstance(small_molecule_ff_raw, str):
            small_molecule_ff_value = small_molecule_ff_raw

    forcefield_raw = _collapse_value_keys(metadata.forcefield, "forcefield")
    forcefield_value: list[str] | str | None = None
    if isinstance(forcefield_raw, str):
        forcefield_value = forcefield_raw
    elif isinstance(forcefield_raw, list):
        if all(isinstance(item, str) for item in forcefield_raw):
            forcefield_value = cast(list[str], forcefield_raw)

    benchmark_results = BenchmarkResults(
        submission_id=final_submission_id,
        title=title,
        summary=summary,
        tags=final_tags,
        calculation_type=metadata.mode,
        authors=[{"name": name} for name in (author or ["TODO add author name"])],
        date=_normalize_submission_date(submission_date, final_submission_id),
        results=results_file,
        archive=Archive(
            doi="TODO add DOI",
            archive_provider="TODO add archive provider",
        ),
        license=license,
        openfe_version=str(_collapse_value_keys(metadata.openfe_version, "version")),
        openmm_version=str(_collapse_value_keys(metadata.openmm_version, "version")),
        openff_toolkit_version=str(
            _collapse_value_keys(metadata.openff_toolkit_version, "version")
        ),
        partial_charges=str(
            _collapse_value_keys(metadata.partial_charges, "partial_charges")
        ),
        benchmark_data=cast(dict[str, object], _build_benchmark_data(metadata.systems)),
        protocol_settings=cast(
            list[dict[str, object]], _build_protocol_payload(metadata.protocol_settings)
        ),
        mapper=mapper_value,
        forcefield=forcefield_value,
        small_molecule_forcefield=small_molecule_ff_value,
        pontibus_version=str(
            _collapse_value_keys(metadata.pontibus_version, "version")
        ),
    )

    submission_yaml_path = out_dir / "submission.yaml"
    zenodo_description_path = out_dir / "zenodo_description.md"

    benchmark_results.write_submission_yaml(submission_yaml_path)
    zenodo_description_path.write_text(
        _make_zenodo_description(
            benchmark_results=benchmark_results,
            network_mode=metadata.network_mode,
            used_alchemiscale=used_alchemiscale,
        )
    )

    logger.info(f"Processed {len(input_paths)} input file(s)")
    logger.info(f"Detected mode: {metadata.mode}")
    logger.info(f"Submission YAML: {submission_yaml_path}")

    return submission_yaml_path, zenodo_description_path


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate submission.yaml and zenodo_description.md from OpenFE JSON archives",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """
            Examples:
              # Single archive file
              %(prog)s archive.json.bz2

              # Multiple archive files
              %(prog)s archive1.json.bz2 archive2.json.bz2 --output-dir ./results

              # Glob pattern
              %(prog)s "networks/*/*.json" --output-dir ./results

              # Multiple glob patterns
              %(prog)s "charge_changes/*/*.json" "jacs_set/*/*.json"

              # Full example with all options
              %(prog)s "networks/*/*.json" \\
                  --output-dir ./output \\
                  --submission-id "2026-06-03-tyk2-rbfe" \\
                  --tags "openfe,rbfe,tyk2" \\
                  --author "Jane Doe" \\
                  --author "John Smith" \\
                  --license "CC-BY-4.0"
            """
        ),
    )

    parser.add_argument(
        "input_patterns",
        type=str,
        nargs="+",
        metavar="INPUT",
        help="One or more file paths or glob patterns (e.g., 'networks/*/*.json'). "
        "Glob patterns support * and ** wildcards.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory for submission.yaml and zenodo_description.md (default: current directory)",
    )
    parser.add_argument(
        "-s",
        "--submission-id",
        type=str,
        default=None,
        help="Submission ID (default: auto-generated from date and network key)",
    )
    parser.add_argument(
        "-t",
        "--tags",
        type=str,
        default="",
        help="Comma-separated extra tags to append (default: none)",
    )
    parser.add_argument(
        "-a",
        "--author",
        type=str,
        action="append",
        dest="authors",
        help="Author name (can be specified multiple times)",
    )
    parser.add_argument(
        "-l",
        "--license",
        type=str,
        default="CC-BY-4.0",
        help="License identifier (default: CC-BY-4.0)",
    )
    parser.add_argument(
        "--submission-date",
        type=str,
        required=True,
        help="Submission date in ISO 8601 format (YYYY-MM-DD). This date is required for submission.yaml.",
    )
    parser.add_argument(
        "--no-alchemiscale",
        action="store_true",
        help="Indicate that Alchemiscale was NOT used to generate the results",
    )
    parser.add_argument(
        "--summary-suffix",
        type=str,
        default=None,
        help="Additional text to append to the auto-generated summary",
    )
    parser.add_argument(
        "-r",
        "--results-file",
        type=str,
        default="computational_results.json",
        help="Name of the results file in output directory (default: computational_results.json)",
    )
    parser.add_argument(
        "--system-group",
        type=str,
        default=None,
        help="Benchmark set name (e.g., 'jacs_set', 'solvation_set'); overrides values from transformation annotations",
    )
    parser.add_argument(
        "--system-name",
        type=str,
        default=None,
        help="System name (e.g., 'tyk2', 'hsp90'); overrides values from transformation annotations",
    )
    parser.add_argument(
        "--forcefields",
        type=str,
        action="append",
        default=None,
        help="Custom protein/solvent force field labels when the archive uses serialized force field contents. Repeat to add multiple values.",
    )
    parser.add_argument(
        "--openfe-version",
        type=str,
        default=None,
        help="Override the OpenFE version written into submission metadata.",
    )
    parser.add_argument(
        "--openmm-version",
        type=str,
        default=None,
        help="Override the OpenMM version written into submission metadata.",
    )
    parser.add_argument(
        "--openff-toolkit-version",
        type=str,
        default=None,
        help="Override the OpenFF Toolkit version written into submission metadata.",
    )
    parser.add_argument(
        "--pontibus-version",
        type=str,
        default=None,
        help="Override the Pontibus version written into submission metadata.",
    )
    parser.add_argument(
        "--small-molecule-forcefield",
        type=str,
        default=None,
        help="Custom small molecule force field label when the archive uses serialized force field contents.",
    )

    return parser


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()

    all_files: list[Path] = []
    for pattern in args.input_patterns:
        matched = glob_module.glob(pattern, recursive=True)
        if matched:
            all_files.extend(Path(path) for path in sorted(matched))
        else:
            all_files.append(Path(pattern))

    if not all_files:
        logger.error("No input files found")
        return 1

    unique_files: list[Path] = []
    seen: set[Path] = set()
    for path in all_files:
        if path not in seen:
            seen.add(path)
            unique_files.append(path)

    process_network(
        input_files=unique_files,
        output_dir=args.output_dir,
        submission_id=args.submission_id,
        tags=args.tags,
        author=args.authors,
        license=args.license,
        used_alchemiscale=not args.no_alchemiscale,
        summary_suffix=args.summary_suffix,
        results_file=args.results_file,
        submission_date=args.submission_date,
        system_group=args.system_group,
        system_name=args.system_name,
        forcefields=args.forcefields,
        small_molecule_forcefield=args.small_molecule_forcefield,
        openfe_version=args.openfe_version,
        openmm_version=args.openmm_version,
        openff_toolkit_version=args.openff_toolkit_version,
        pontibus_version=args.pontibus_version,
    )

    logger.info("Successfully generated submission metadata")
    return 0


if __name__ == "__main__":
    sys.exit(main())
