#!/usr/bin/env python3
"""Prepare benchmark submission artifacts from AlchemicalArchive or AlchemicalNetwork JSON files.

This module generates `submission.yaml` and `zenodo_description.md` from one or more JSON archives
exported by OpenFE/Alchemiscale. Supports single files, lists of files, or glob patterns with 
potentially different protocol settings.

Example (Python API):

    from pathlib import Path
    from openfe_benchmarks.scripts.prepare_metadata_submission import process_network

    # Using explicit file list
    process_network(
        input_files=[Path("archive1.json.bz2"), Path("archive2.json.bz2")],
        output_dir=Path("."),
        submission_id="2026-04-15-example",
        tags="openfe,alchemicalarchive",
        author=["Jane Doe"],
        license="CC-BY-4.0",
    )
    
    # Using glob pattern
    process_network(
        input_files="networks/*/*.json",
        output_dir=Path("."),
        submission_id="2026-04-15-example",
        tags="openfe,alchemicalarchive",
        author=["Jane Doe"],
        license="CC-BY-4.0",
    )

Example (CLI):

    # Explicit files
    python prepare_metadata_submission.py archive1.json.bz2 archive2.json.bz2 \\
        --output-dir ./output \\
        --submission-id "2026-04-15-example" \\
        --tags "openfe,alchemicalarchive" \\
        --author "Jane Doe" \\
        --license "CC-BY-4.0"
    
    # Glob pattern
    python prepare_metadata_submission.py "networks/*/*.json" \\
        --output-dir ./output \\
        --submission-id "2026-04-15-example"
"""

from __future__ import annotations

import argparse
import bz2
import glob as glob_module
import json
import re
import sys
import textwrap
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

from openfe_benchmarks.data import BenchmarkIndex


@dataclass
class ProtocolSettingsInfo:
    """Container for protocol settings with source metadata."""

    settings: dict[str, str]
    source_file: str
    benchmark_set: str
    benchmark_system: str
    network_key: str


@dataclass
class SystemInfo:
    """Per-system information extracted from transformations."""

    system_group: str
    system_name: str
    n_transformations: int
    ligands: set[str] = field(default_factory=set)
    proteins: set[str] = field(default_factory=set)
    solutes: set[str] = field(default_factory=set)
    solvents: set[str] = field(default_factory=set)
    max_atoms: int = 0
    files: set[str] = field(default_factory=set)


@dataclass
class AutoMetadata:
    openfe_version: str = ""
    openmm_version: str = ""
    openff_toolkit_version: str = ""
    forcefield: str = ""
    partial_charges: str = ""
    benchmark_data_set: str = ""
    benchmark_system: str = ""
    protocol_settings_list: list[ProtocolSettingsInfo] = field(default_factory=list)
    system_info_list: list[SystemInfo] = field(default_factory=list)


def _require_ref_key(ref: dict[str, Any]) -> str:
    key = ref.get(":gufe-key:")
    if not key:
        raise KeyError(f"Expected :gufe-key: reference, got: {ref}")
    return key


def _get_network_key(
    archive_obj: dict[str, Any] | None, network_obj: dict[str, Any] | None
) -> str:
    if archive_obj is not None:
        network_ref = archive_obj.get("network")
        if isinstance(network_ref, dict):
            return _require_ref_key(network_ref)

    if network_obj is not None:
        name = network_obj.get("name")
        if isinstance(name, str) and name:
            return name
        return "NA"

    raise ValueError("Could not determine network key from input file")


def _open_json_file(path: Path):
    if path.suffix == ".bz2":
        return bz2.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def _load_token_table(
    input_path: Path,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any] | None, dict[str, Any] | None]:
    with _open_json_file(input_path) as f:
        token_table = json.load(f)

    if isinstance(token_table, dict):
        token_table = list(token_table.items())

    if not isinstance(token_table, list):
        raise ValueError("Unsupported JSON token table format")

    by_key: dict[str, dict[str, Any]] = {}
    archive_obj: dict[str, Any] | None = None
    network_obj: dict[str, Any] | None = None

    for item in token_table:
        if not (isinstance(item, list) and len(item) == 2):
            continue
        table_key, payload = item
        if not isinstance(payload, dict):
            continue

        gufe_key = payload.get(":gufe-key:") or table_key
        if isinstance(gufe_key, str):
            by_key[gufe_key] = payload

        qualname = payload.get("__qualname__")
        if qualname == "AlchemicalArchive":
            archive_obj = payload
        elif qualname == "AlchemicalNetwork":
            network_obj = payload

    if archive_obj is None and network_obj is None:
        raise ValueError(
            "Could not find AlchemicalArchive or AlchemicalNetwork object in token table"
        )

    return by_key, archive_obj, network_obj


def _transformation_refs(
    archive_obj: dict[str, Any] | None, network_obj: dict[str, Any] | None
) -> list[Any]:
    if archive_obj is not None:
        return [
            item[0]
            for item in archive_obj.get("transformation_results", [])
            if isinstance(item, list) and len(item) == 2
        ]
    if network_obj is not None:
        edges = network_obj.get("edges") or []
        if isinstance(edges, list):
            return edges
    return []


def _detect_mode(
    by_key: dict[str, dict[str, Any]],
    archive_obj: dict[str, Any] | None,
    network_obj: dict[str, Any] | None,
) -> str:
    names: list[str] = []
    for transformation_ref in _transformation_refs(archive_obj, network_obj)[:20]:
        transformation = by_key[_require_ref_key(transformation_ref)]
        name = str(transformation.get("name") or "")
        if name:
            names.append(name)

    if any(n.startswith("complex_") or n.startswith("solvent_") for n in names):
        return "rbfe"
    return "asfe"


def _slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")


def _default_submission_id(network_key: str) -> str:
    return f"{date.today().isoformat()}-{_slugify(network_key)}"


def _generate_title(
    mode: str,
    benchmark_sets: list[str],
    systems: list[str],
    submission_id: str,
) -> str:
    """
    Generate a descriptive title for the submission.

    Format rules:
    - Single set, 1-3 systems: "OpenFE RBFE - set - sys1, sys2, sys3 - submission_id"
    - Single set, 4+ systems: "OpenFE RBFE - set (N systems) - submission_id"
    - 2-3 sets, any systems: "OpenFE RBFE - set1, set2 (N systems) - submission_id"
    - 4+ sets: "OpenFE RBFE - Multi-set Benchmark (N sets, M systems) - submission_id"
    """
    mode = mode.upper()
    n_sets = len(benchmark_sets)
    n_systems = len(systems)

    if n_sets == 0:
        # Fallback if no benchmark set detected
        return f"OpenFE {mode} Benchmark - {submission_id}"

    if n_sets == 1:
        set_name = benchmark_sets[0]
        if n_systems <= 3:
            # List system names
            systems_str = ", ".join(systems)
            return f"OpenFE {mode} - {set_name} - {systems_str} - {submission_id}"
        else:
            # Use count
            return f"OpenFE {mode} - {set_name} ({n_systems} systems) - {submission_id}"

    if n_sets <= 3:
        # List set names with system count
        sets_str = ", ".join(benchmark_sets)
        return f"OpenFE {mode} - {sets_str} ({n_systems} systems) - {submission_id}"

    # Many sets - use multi-set notation
    return f"OpenFE {mode} - Multi-set Benchmark ({n_sets} sets, {n_systems} systems) - {submission_id}"


def _iter_nested_items(obj: Any) -> list[tuple[str, Any]]:
    items: list[tuple[str, Any]] = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            items.append((str(k), v))
            items.extend(_iter_nested_items(v))
    elif isinstance(obj, list):
        for v in obj:
            items.extend(_iter_nested_items(v))
    return items


def _quantity_to_text(value: Any) -> str:
    if isinstance(value, dict) and "magnitude" in value and "unit" in value:
        return f"{value['magnitude']} {value['unit']}"
    return str(value)


def _extract_system_info_from_mapping(
    by_key: dict[str, dict[str, Any]], transformation_ref: Any
) -> tuple[str, str]:
    """
    Extract system_group and system_name from LigandAtomMapping annotations.

    Returns:
        (system_group, system_name) tuple, or ("", "") if not found
    """
    transformation = by_key.get(_require_ref_key(transformation_ref), {})
    mapping_ref = transformation.get("mapping")
    if not mapping_ref:
        return ("", "")

    mapping = by_key.get(_require_ref_key(mapping_ref), {})
    if mapping.get("__qualname__") != "LigandAtomMapping":
        return ("", "")

    annotations = mapping.get("annotations")
    if isinstance(annotations, str):
        try:
            annotations = json.loads(annotations)
        except json.JSONDecodeError:
            return ("", "")

    if not isinstance(annotations, dict):
        return ("", "")

    system_group = annotations.get("system_group", "")
    system_name = annotations.get("system_name", "")
    return (str(system_group), str(system_name))


def _infer_benchmark_data_set_system(
    *, by_key: dict[str, dict[str, Any]], mode: str, archive_stem: str, network_key: str
) -> tuple[str, str]:
    """Infer benchmark set and system from file contents using BenchmarkIndex.

    This searches for any known benchmark set or system name in the file's
    metadata (filename, network key, and JSON contents).

    Returns:
        (benchmark_set, system_name) tuple, or ("", "") if not found
    """
    blob = json.dumps(list(by_key.values())).lower()
    search_space = " ".join([blob, archive_stem.lower(), network_key.lower()])

    # Get all known benchmark sets and systems from the index
    index = BenchmarkIndex()
    benchmark_sets = index.list_benchmark_sets()

    # Check for benchmark set matches
    benchmark_set = ""
    for set_name in benchmark_sets:
        if set_name.lower() in search_space:
            benchmark_set = set_name
            break

    # Check for system name matches within the found set (or all sets if no set found)
    system = ""
    sets_to_check = [benchmark_set] if benchmark_set else benchmark_sets

    for set_name in sets_to_check:
        try:
            systems = index.list_systems_by_benchmark_set(set_name)
            for system_name in systems:
                if system_name.lower() in search_space:
                    system = system_name
                    if not benchmark_set:
                        benchmark_set = set_name
                    break
            if system:
                break
        except ValueError:
            # Skip if benchmark set doesn't exist
            continue

    return benchmark_set, system


def _extract_sim_times(settings_block: dict[str, Any]) -> tuple[str, str]:
    equilibration = settings_block.get("equilibration_length")
    production = settings_block.get("production_length")
    return _quantity_to_text(
        equilibration
    ) if equilibration is not None else "", _quantity_to_text(
        production
    ) if production is not None else ""


def _build_protocol_settings(
    protocol_obj: dict[str, Any] | None, mode: str
) -> dict[str, str]:
    if not protocol_obj:
        return {
            "protocol": "unknown",
            "notes": "Protocol settings unavailable in archive payload.",
        }

    # Detect protocol name from the object
    protocol_name = (
        protocol_obj.get("__qualname__") or protocol_obj.get("qualname") or "unknown"
    )
    out: dict[str, str] = {"protocol": protocol_name}

    settings = protocol_obj.get("settings") or {}

    integrator_settings = settings.get("integrator_settings") or {}
    if isinstance(integrator_settings, dict):
        timestep = integrator_settings.get("timestep")
        if timestep is not None:
            out["timestep"] = _quantity_to_text(timestep)

    thermo_settings = settings.get("thermo_settings") or {}
    if isinstance(thermo_settings, dict):
        temperature = thermo_settings.get("temperature")
        pressure = thermo_settings.get("pressure")
        if temperature is not None:
            out["temperature"] = _quantity_to_text(temperature)
        if pressure is not None:
            out["pressure"] = _quantity_to_text(pressure)

    lambda_settings = settings.get("lambda_settings") or {}
    if isinstance(lambda_settings, dict):
        lambda_windows = lambda_settings.get("lambda_windows")
        if lambda_windows is not None:
            out["lambda_windows"] = str(lambda_windows)
        else:
            lambda_counts: list[str] = []
            for lambda_key in ("lambda_elec", "lambda_vdw", "lambda_restraints"):
                values = lambda_settings.get(lambda_key)
                if isinstance(values, list):
                    lambda_counts.append(f"{lambda_key}:{len(values)}")
            if lambda_counts:
                out["lambda_schedule"] = ", ".join(lambda_counts)

    # Protocol-specific handling: RBFE typically has a single simulation block;
    # ASFE commonly has separate vacuum and solvent simulation settings.
    if mode == "rbfe":
        sim = settings.get("simulation_settings") or {}
        if isinstance(sim, dict):
            eq, prod = _extract_sim_times(sim)
            if eq:
                out["equilibration_time"] = eq
            if prod:
                out["production_time"] = prod
    elif mode == "asfe":
        for prefix, key in (
            ("vacuum", "vacuum_simulation_settings"),
            ("solvent", "solvent_simulation_settings"),
        ):
            sim = settings.get(key) or {}
            if not isinstance(sim, dict):
                continue
            eq, prod = _extract_sim_times(sim)
            if eq:
                out[f"{prefix}_equilibration_time"] = eq
            if prod:
                out[f"{prefix}_production_time"] = prod
    else:
        ValueError(
            f"Calculation type {mode} is not yet supported. Add capability to `_build_protocol_settings`"
        )

    if len(out) == 1:
        out["notes"] = "Protocol class found, but detailed settings were unavailable."

    return out


def _render_protocol_settings_yaml(
    protocol_settings_list: list[ProtocolSettingsInfo],
) -> str:
    """Render protocol settings as YAML, grouping by unique setting combinations."""
    if not protocol_settings_list:
        return "protocol_settings:\n  - protocol: unknown\n    notes: Protocol settings unavailable."

    lines = ["protocol_settings:"]

    # Group identical protocol settings and collect source metadata
    settings_groups: dict[str, list[ProtocolSettingsInfo]] = defaultdict(list)

    for info in protocol_settings_list:
        # Create a hashable key from the settings only
        settings_key = json.dumps(info.settings, sort_keys=True)
        settings_groups[settings_key].append(info)

    preferred_order = [
        "production_time",
        "equilibration_time",
        "vacuum_production_time",
        "vacuum_equilibration_time",
        "solvent_production_time",
        "solvent_equilibration_time",
        "timestep",
        "temperature",
        "pressure",
        "lambda_windows",
        "lambda_schedule",
        "notes",
    ]

    for _, (settings_key, info_list) in enumerate(settings_groups.items()):
        settings = info_list[0].settings
        protocol_name = settings.get("protocol", "")
        lines.append(f"  - protocol: {protocol_name}")

        # Add count and source information
        if len(protocol_settings_list) > 1:
            lines.append(f"    count: {len(info_list)}")

            # Collect unique benchmark systems
            systems = sorted(
                set(
                    f"{info.benchmark_set}/{info.benchmark_system}"
                    for info in info_list
                    if info.benchmark_set or info.benchmark_system
                )
            )
            if systems:
                systems_str = ", ".join(systems)
                lines.append(f"    systems: {json.dumps(systems_str)}")

            # Collect source files (just filenames, not full paths)
            files = sorted(set(Path(info.source_file).name for info in info_list))
            if files and len(files) <= 5:  # Only show if not too many
                files_str = ", ".join(files)
                lines.append(f"    files: {json.dumps(files_str)}")
            elif files:
                lines.append(f"    files: {json.dumps(f'{len(files)} files')}")

        # Add protocol settings
        for key in preferred_order:
            if key in settings:
                lines.append(f"    {key}: {json.dumps(str(settings[key]))}")

        for key in sorted(
            k for k in settings if k not in set(preferred_order) | {"protocol"}
        ):
            lines.append(f"    {key}: {json.dumps(str(settings[key]))}")

    return "\n".join(lines)


def _resolve_payload(
    by_key: dict[str, dict[str, Any]], ref: Any
) -> tuple[str | None, dict[str, Any] | None]:
    if isinstance(ref, dict):
        ref_key = ref.get(":gufe-key:")
        if isinstance(ref_key, str):
            return ref_key, by_key.get(ref_key)
    return None, None


def _component_name(component: dict[str, Any], component_key: str) -> str:
    molprops = component.get("molprops") or {}
    if isinstance(molprops, dict):
        ofe_name = molprops.get("ofe-name")
        if ofe_name:
            return str(ofe_name)

    if component.get("name"):
        return str(component.get("name"))

    if component.get("smiles"):
        return str(component.get("smiles"))

    if component.get("solvent_molecule"):
        return str(component.get("solvent_molecule"))

    return component_key


def _component_atoms(component: dict[str, Any]) -> int:
    atoms = component.get("atoms")
    if isinstance(atoms, list):
        return len(atoms)
    return 0


def _build_content_summary(
    by_key: dict[str, dict[str, Any]],
    archive_objs: list[dict[str, Any]],
    network_objs: list[dict[str, Any]],
    mode: str,
    benchmark_data_set: str,
    forcefield: str,
    partial_charges: str,
    used_alchemiscale: bool = True,
    source_files: list[str] = None,
) -> tuple[str, list[SystemInfo]]:
    """
    Build content summary and extract per-system information.

    Parameters
    ----------
    archive_objs: List of AlchemicalArchive objects
    network_objs: List of AlchemicalNetwork objects

    Returns:
        (summary_text, list of SystemInfo objects)
    """
    if source_files is None:
        source_files = []

    repeat_counts: list[int] = []
    transformation_refs: list[Any] = []

    # Collect transformation refs from all archives and networks
    for archive_obj in archive_objs:
        transformation_results = archive_obj.get("transformation_results", [])
        trans_refs = [
            item[0]
            for item in transformation_results
            if isinstance(item, list) and len(item) == 2
        ]
        transformation_refs.extend(trans_refs)
        repeat_counts.extend(
            [
                len(item[1])
                for item in transformation_results
                if isinstance(item, list) and len(item) == 2
            ]
        )

    for network_obj in network_objs:
        edges = network_obj.get("edges") or []
        if isinstance(edges, list):
            transformation_refs.extend(edges)

    # Per-system tracking
    system_data: dict[
        tuple[str, str], SystemInfo
    ] = {}  # key: (system_group, system_name)

    visited_systems_for_cofactors: set[str] = set()
    all_cofactors: set[str] = set()
    systems_with_cofactors: set[str] = set()

    transformation_count = len(transformation_refs)

    for transformation_ref in transformation_refs:
        _, transformation = _resolve_payload(by_key, transformation_ref)
        if not transformation:
            continue

        # Extract system info from mapping annotations
        system_group, system_name = _extract_system_info_from_mapping(
            by_key, transformation_ref
        )

        # Initialize SystemInfo if first time seeing this system
        system_key = (system_group, system_name)
        if system_key not in system_data:
            system_data[system_key] = SystemInfo(
                system_group=system_group,
                system_name=system_name,
                n_transformations=0,
            )

        system_info = system_data[system_key]
        system_info.n_transformations += 1

        for state_key in ("stateA", "stateB"):
            cs_key, chemical_system = _resolve_payload(
                by_key, transformation.get(state_key)
            )
            if not chemical_system:
                continue

            components = chemical_system.get("components") or {}
            if not isinstance(components, dict):
                continue

            system_atoms = 0
            local_cofactors: set[str] = set()

            for label, component_ref in components.items():
                comp_key, component = _resolve_payload(by_key, component_ref)
                if not component:
                    continue

                label_l = str(label).lower()
                qualname = str(component.get("__qualname__") or "")
                comp_name = _component_name(component, comp_key or "unknown")
                comp_atoms = _component_atoms(component)
                system_atoms += comp_atoms

                if mode == "asfe":
                    if "solvent" in label_l or "solventcomponent" in qualname.lower():
                        system_info.solvents.add(comp_name)
                    elif "solute" in label_l or qualname == "SmallMoleculeComponent":
                        system_info.solutes.add(comp_name)
                elif mode == "rbfe":
                    if "protein" in label_l or qualname == "ProteinComponent":
                        system_info.proteins.add(comp_name)
                    elif "ligand" in label_l:
                        system_info.ligands.add(comp_name)
                    elif "cofactor" in label_l:
                        all_cofactors.add(comp_name)
                        local_cofactors.add(comp_name)
                    elif (
                        qualname == "SmallMoleculeComponent"
                        and "solvent" not in label_l
                    ):
                        # Non-solvent small molecules that are not explicit ligands are treated as cofactors.
                        all_cofactors.add(comp_name)
                        local_cofactors.add(comp_name)
                else:
                    ValueError(
                        f"Calculation type {mode} is not yet supported. Add capability to `_build_content_summary`"
                    )

            system_info.max_atoms = max(system_info.max_atoms, system_atoms)

            if (
                mode == "rbfe"
                and cs_key
                and cs_key not in visited_systems_for_cofactors
            ):
                visited_systems_for_cofactors.add(cs_key)
                if local_cofactors:
                    systems_with_cofactors.add(cs_key)

    # Add source files to system_info
    for system_info in system_data.values():
        system_info.files = set(source_files)

    field_info = forcefield or "an unspecified force field"
    charge_info = partial_charges or "unspecified partial charges"

    # Group systems by benchmark set for explicit listing
    sets_to_systems: dict[str, list[str]] = defaultdict(list)
    for si in system_data.values():
        if si.system_group and si.system_name:
            sets_to_systems[si.system_group].append(si.system_name)

    # Sort systems within each set
    for systems_list in sets_to_systems.values():
        systems_list.sort()

    unique_sets = sorted(sets_to_systems.keys())

    # Build descriptive subject line
    if len(unique_sets) == 0:
        subject = benchmark_data_set or "benchmark"
    elif len(unique_sets) == 1:
        subject = unique_sets[0]
    else:
        subject = ", ".join(unique_sets)

    # Build explicit systems description: "set1: sys1, sys2; set2: sys3, sys4"
    if len(unique_sets) > 1:
        set_descriptions = [
            f"{set_name}: {', '.join(sets_to_systems[set_name])}"
            for set_name in unique_sets
        ]
        systems_desc = "; ".join(set_descriptions)
    elif len(unique_sets) == 1:
        systems_desc = ", ".join(sets_to_systems[unique_sets[0]])
    else:
        systems_desc = f"{len(system_data)} systems"

    # Count totals across all systems
    total_ligands = sum(len(si.ligands) for si in system_data.values())
    total_proteins = sum(len(si.proteins) for si in system_data.values())
    total_solutes = sum(len(si.solutes) for si in system_data.values())
    total_solvents = sum(len(si.solvents) for si in system_data.values())
    max_atoms_overall = max((si.max_atoms for si in system_data.values()), default=0)

    # Build summary
    if mode == "rbfe":
        cofactor_list = ", ".join(sorted(all_cofactors)) if all_cofactors else "none"
        if len(system_data) > 1:
            summary_parts = [
                f"This submission describes the {subject} RBFE benchmark ({systems_desc}) prepared with {field_info} and {charge_info}.",
                f"The submission contains {transformation_count} transformations, {total_ligands} unique ligands, and {total_proteins} unique proteins.",
            ]
        else:
            summary_parts = [
                f"This submission describes the {subject} RBFE benchmark prepared with {field_info} and {charge_info}.",
                f"The network contains {transformation_count} transformations across {total_ligands} unique ligands and {total_proteins} unique proteins.",
            ]
        if systems_with_cofactors:
            summary_parts.append(
                f"{len(systems_with_cofactors)} systems include cofactors ({cofactor_list})."
            )
    else:
        if len(system_data) > 1:
            summary_parts = [
                f"This submission describes the {subject} ASFE benchmark ({systems_desc}) prepared with {field_info} and {charge_info}.",
                f"The submission contains {transformation_count} transformations, {total_solutes} unique solutes, and {total_solvents} unique solvents.",
            ]
        else:
            summary_parts = [
                f"This submission describes the {subject} ASFE benchmark prepared with {field_info} and {charge_info}.",
                f"The archive contains {transformation_count} transformations across {total_solutes} unique solutes and {total_solvents} unique solvents.",
            ]

    summary_parts.append(
        f"The largest simulated chemical system contains {max_atoms_overall} atoms."
    )
    if used_alchemiscale:
        summary_parts.append(
            "Results are derived from archived Alchemiscale workflow data."
        )
    summary_text = " ".join(summary_parts)
    return textwrap.fill(summary_text, width=100), list(system_data.values())


def _extract_auto_metadata(
    *,
    by_key: dict[str, dict[str, Any]],
    mode: str,
    archive_path: Path,
    network_key: str,
    archive_stem: str,
) -> AutoMetadata:
    metadata = AutoMetadata()
    metadata.benchmark_data_set, metadata.benchmark_system = (
        _infer_benchmark_data_set_system(
            by_key=by_key,
            mode=mode,
            archive_stem=archive_stem,
            network_key=network_key,
        )
    )

    protocol_obj: dict[str, Any] | None = None

    for payload in by_key.values():
        qualname = payload.get("__qualname__")
        if qualname in {
            "RelativeHybridTopologyProtocol",
            "AbsoluteSolvationProtocol",
            "ASFEProtocol",
        }:
            protocol_obj = payload

            settings = payload.get("settings") or {}
            forcefield_settings = (
                settings.get("forcefield_settings")
                or settings.get("solvent_forcefield_settings")
                or settings.get("vacuum_forcefield_settings")
                or {}
            )
            if not metadata.forcefield and isinstance(forcefield_settings, dict):
                metadata.forcefield = str(
                    forcefield_settings.get("small_molecule_forcefield") or ""
                )
                if not metadata.forcefield:
                    ffs = forcefield_settings.get("forcefields")
                    if isinstance(ffs, list) and ffs:
                        preferred = next(
                            (
                                ff
                                for ff in ffs
                                if isinstance(ff, str) and "openff" in ff
                            ),
                            None,
                        )
                        if preferred:
                            metadata.forcefield = str(preferred).replace(".offxml", "")

            partial_charge_settings = settings.get("partial_charge_settings") or {}
            if not metadata.partial_charges and isinstance(
                partial_charge_settings, dict
            ):
                method = partial_charge_settings.get("partial_charge_method")
                nagl_model = partial_charge_settings.get("nagl_model")
                if method and nagl_model:
                    metadata.partial_charges = f"{method} ({nagl_model})"
                elif method:
                    metadata.partial_charges = str(method)

        for key, value in _iter_nested_items(payload):
            if (
                not metadata.openmm_version
                and key == "openmm_version"
                and isinstance(value, str)
            ):
                metadata.openmm_version = value
            if (
                not metadata.openfe_version
                and key == "openfe_version"
                and isinstance(value, str)
            ):
                metadata.openfe_version = value
            if (
                not metadata.openff_toolkit_version
                and key == "openff_toolkit_version"
                and isinstance(value, str)
            ):
                metadata.openff_toolkit_version = value

        molprops = payload.get("molprops")
        if isinstance(molprops, dict):
            charge_provenance = molprops.get("charge_provenance")
            if isinstance(charge_provenance, str):
                try:
                    charge_provenance = json.loads(charge_provenance)
                except json.JSONDecodeError:
                    charge_provenance = None
            if isinstance(charge_provenance, dict):
                if not metadata.openfe_version and isinstance(
                    charge_provenance.get("openfe_version"), str
                ):
                    metadata.openfe_version = charge_provenance["openfe_version"]
                if not metadata.openff_toolkit_version and isinstance(
                    charge_provenance.get("openff_toolkit_version"), str
                ):
                    metadata.openff_toolkit_version = charge_provenance[
                        "openff_toolkit_version"
                    ]
                if not metadata.partial_charges:
                    charge_method = charge_provenance.get("charge_method")
                    nagl_model = charge_provenance.get("nagl_model")
                    if charge_method and nagl_model:
                        metadata.partial_charges = f"{charge_method} ({nagl_model})"
                    elif charge_method:
                        metadata.partial_charges = str(charge_method)

    # Build protocol settings for this archive with source metadata
    protocol_settings = _build_protocol_settings(protocol_obj, mode)
    protocol_info = ProtocolSettingsInfo(
        settings=protocol_settings,
        source_file=str(archive_path),
        benchmark_set=metadata.benchmark_data_set,
        benchmark_system=metadata.benchmark_system,
        network_key=network_key,
    )
    metadata.protocol_settings_list = [protocol_info]

    return metadata


def _normalize_partial_charge_info(partial_charges: str) -> str:
    value = partial_charges.strip()
    if not value:
        return ""

    lower = value.lower()

    # Canonical openfe-benchmarks style: nagl_<model>.pt
    model_match = re.search(
        r"(openff-gnn-am1bcc-[0-9.]+\.0\.pt|openff-gnn-am1bcc-[0-9.]+\.pt)", lower
    )
    if "nagl" in lower and model_match:
        model = model_match.group(1)
        return f"nagl_{model}"

    if lower in {"am1bcc", "am1-bcc"}:
        return "am1bcc"

    normalized = re.sub(r"[^a-z0-9._-]+", "_", lower).strip("_")
    return normalized


def _make_tags(
    *, mode: str, forcefield: str, partial_charge_tag: str, user_keywords: list[str]
) -> list[str]:
    tags: list[str] = []
    tags.append(mode)
    if forcefield:
        tags.append(forcefield)
    if partial_charge_tag:
        tags.append(partial_charge_tag)
    tags.extend(user_keywords)

    # Deduplicate while preserving order.
    out: list[str] = []
    seen: set[str] = set()
    for tag in tags:
        t = tag.strip()
        if not t or t in seen:
            continue
        seen.add(t)
        out.append(t)
    return out


def _yaml_block(text: str, indent_spaces: int = 2) -> str:
    indent = " " * indent_spaces
    lines = text.splitlines() or [""]
    return "\n".join(f"{indent}{line}" for line in lines)


def _make_submission_yaml(
    *,
    submission_id: str,
    title: str,
    summary: str,
    tags: list[str],
    authors: list[str],
    openfe_version: str,
    openmm_version: str,
    openff_toolkit_version: str,
    forcefield: str,
    partial_charges: str,
    benchmark_data_set: str,
    benchmark_system: str,
    archive_doi: str,
    archive_provider: str,
    license_name: str,
    protocol_settings_list: list[ProtocolSettingsInfo],
    network_key_to_systems: dict[str, list[str]],
    results_file: str,
) -> str:
    if not authors:
        authors = ["TODO add author name"]

    tags_yaml = ", ".join(tags)
    authors_yaml = "\n".join(f"  - name: {name}" for name in authors)
    protocol_settings_yaml = _render_protocol_settings_yaml(protocol_settings_list)

    # Render network_key to systems mapping
    network_keys_yaml = ""
    if network_key_to_systems:
        network_keys_yaml = "\n# Network keys to systems mapping\nnetwork_keys:\n"
        for network_key in sorted(network_key_to_systems.keys()):
            systems = ", ".join(network_key_to_systems[network_key])
            network_keys_yaml += f"  {json.dumps(network_key)}: {json.dumps(systems)}\n"

    return f"""# REQUIRED: unique, kebab-case identifier for this submission
submission_id: {submission_id}

# REQUIRED: short descriptive title
title: {title}

# REQUIRED: short descriptive summary (1-2 sentences)
summary: |
{_yaml_block(summary, 2)}

# REQUIRED: list of submission tags
tags: [{tags_yaml}]

# REQUIRED: list of contributing authors (name, affiliation; ORCID optional)
authors:
{authors_yaml}

# REQUIRED: publication/submission date (ISO 8601)
date: {date.today().isoformat()}

# REQUIRED: OpenFE/OpenMM/toolkit versions used to produce gathered reports
openfe_version: {openfe_version}
openmm_version: {openmm_version}
openff_toolkit_version: {openff_toolkit_version}

# Recommended descriptors
forcefield: {forcefield}
partial_charges: {partial_charges}

# BenchmarkData provenance (from openfe-benchmarks planning script)
benchmark_data:
  source_repository: https://github.com/OpenFreeEnergy/openfe-benchmarks
  set: {benchmark_data_set}
  system: {benchmark_system}
{network_keys_yaml}
# REQUIRED: results file
results: {results_file}

# REQUIRED: long-term archive pointer (at least doi or url)
archive:
  doi: {archive_doi}
  archive_provider: {archive_provider}

# REQUIRED: license for the submission
license: {license_name}

# RECOMMENDED / OPTIONAL metadata for protocol settings
{protocol_settings_yaml}
"""


def _make_zenodo_description(
    *,
    title: str,
    archive_filename: str,
    mode: str,
    content_summary: str,
    openfe_version: str,
    openmm_version: str,
    openff_toolkit_version: str,
    forcefield: str,
    partial_charges: str,
    benchmark_data_set: str,
    benchmark_system: str,
    license_name: str,
    protocol_settings_list: list[ProtocolSettingsInfo],
    has_archive_objects: bool,
    used_alchemiscale: bool,
    system_info_list: list[SystemInfo] = None,
    network_key_to_systems: dict[str, list[str]] = None,
) -> str:
    content_kind = "ASFE" if mode == "asfe" else "RBFE"
    openmm_display = openmm_version or "<not found in alchemical network>"

    # Determine the source type for the overview text
    if has_archive_objects:
        source_description = "AlchemicalArchive"
    else:
        source_description = "AlchemicalNetwork"

    # Build workflow description
    workflow_text = "OpenFE"
    if used_alchemiscale:
        workflow_text += " and Alchemiscale"

    if system_info_list is None:
        system_info_list = []

    # Group protocol settings by unique combinations
    settings_groups: dict[str, list[ProtocolSettingsInfo]] = defaultdict(list)
    for info in protocol_settings_list:
        settings_key = json.dumps(info.settings, sort_keys=True)
        settings_groups[settings_key].append(info)

    # Build protocol settings section
    protocol_blocks: list[str] = []
    for idx, (settings_key, info_list) in enumerate(settings_groups.items()):
        protocol_settings = info_list[0].settings

        if len(settings_groups) > 1:
            protocol_blocks.append(f"\n### Protocol Settings Group {idx + 1}")

        # Add system information
        systems = sorted(
            set(
                f"{info.benchmark_set}/{info.benchmark_system}"
                for info in info_list
                if info.benchmark_set or info.benchmark_system
            )
        )
        if systems:
            protocol_blocks.append(
                f"**Systems ({len(info_list)} networks):** {', '.join(systems)}"
            )

        protocol_lines = "\n".join(
            [f"- {k}: {v}" for k, v in protocol_settings.items()]
        )
        if not protocol_lines:
            protocol_lines = "- notes: Protocol settings unavailable"
        protocol_blocks.append(protocol_lines)

    protocol_section = "\n".join(protocol_blocks)

    # Per-system details section removed per user request
    system_details_section = ""

    # Build network keys to systems mapping section
    # Only show for AlchemicalArchive files (network_key_to_systems will be populated)
    # For AlchemicalNetwork files, the "name" field is not a meaningful network key
    network_keys_section = ""
    if network_key_to_systems:
        network_keys_lines = []
        for network_key_item in sorted(network_key_to_systems.keys()):
            systems = ", ".join(network_key_to_systems[network_key_item])
            network_keys_lines.append(f"  - {network_key_item}: {systems}")
        if network_keys_lines:
            network_keys_section = "\n- network keys:\n" + "\n".join(network_keys_lines)

    return f"""# {title}

## Description

## Overview

{content_kind} benchmark results prepared from {source_description} JSON file(s) generated with {workflow_text}.

{content_summary}

## Software versions

- openfe_version: {openfe_version}
- openmm_version: {openmm_display}
- openff_toolkit_version: {openff_toolkit_version}

## Recommended descriptors

- forcefield: {forcefield}
- partial_charges: {partial_charges}

## BenchmarkData provenance

- source_repository: https://github.com/OpenFreeEnergy/openfe-benchmarks
- benchmark_set: {benchmark_data_set}
- systems: {benchmark_system}

## Protocol settings

{protocol_section}

- archive file: {archive_filename}{network_keys_section}
{system_details_section}

## Rights

- License: {license_name}
"""


def process_network(
    input_files: Path | list[Path] | str,
    output_dir: Path = Path("."),
    submission_id: str | None = None,
    tags: str = "openfe,alchemicalarchive",
    author: list[str] | None = None,
    license: str = "CC-BY-4.0",
    used_alchemiscale: bool = True,
    summary_suffix: str | None = None,
    results_file: str = "computational_results.json",
) -> tuple[Path, Path]:
    """Generate submission metadata from one or more archived OpenFE JSON networks.

    Parameters
    ----------
    input_files:
        Path to a single AlchemicalArchive/AlchemicalNetwork JSON file, a list
        of such files, or a glob pattern string (e.g., "networks/*/*.json").
        Supported extensions are `.json`, `.bz2`, and `.json.bz2`.
        When multiple files are provided, protocol settings from each are collected
        and grouped in the output.
    output_dir:
        Directory where `submission.yaml` and `zenodo_description.md` will be
        written. Defaults to the current working directory.
    submission_id:
        Optional identifier to use in `submission.yaml`. If omitted, a default
        value is generated from the current date and network key.
    tags:
        Comma-separated list of additional tags to include in the submission
        metadata. The generated tag list also always includes the detected
        `mode` (either ``asfe`` or ``rbfe``), the resolved forcefield string,
        and normalized partial charge information.
    author:
        Optional list of author entries for the submission YAML. Each entry is
        treated as a raw string and written to the `authors` section.
    license:
        License string to write into the submission metadata.
    used_alchemiscale:
        Whether Alchemiscale was used to generate the results. If True, the
        description will mention Alchemiscale. Defaults to True.
    summary_suffix:
        Optional text to append to the auto-generated summary.
    results_file:
        Name of the results file to reference in submission.yaml and validate
        exists in output_dir. Defaults to 'computational_results.json'.

    Notes
    -----
    The generated submission YAML will include placeholder values for archive
    DOI and provider, which are intentionally not propagated into the Zenodo
    description output.

    Returns
    -------
    tuple[Path, Path]
        Paths to the generated `submission.yaml` and `zenodo_description.md`.
    """
    # Normalize input to list and expand glob patterns
    if isinstance(input_files, str):
        # Glob pattern
        matched_files = glob_module.glob(input_files, recursive=True)
        if not matched_files:
            raise ValueError(f"No files matched glob pattern: {input_files}")
        input_paths = [Path(f) for f in sorted(matched_files)]
    elif isinstance(input_files, Path):
        input_paths = [input_files]
    else:
        input_paths = input_files

    if not input_paths:
        raise ValueError("At least one input file must be provided")

    # Validate all input files exist
    for input_path in input_paths:
        resolved_path = input_path.resolve()
        if not resolved_path.exists():
            raise FileNotFoundError(f"Input file not found: {resolved_path}")

    out_dir = output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Check for required results file
    results_path = out_dir / results_file
    if not results_path.exists():
        raise FileNotFoundError(
            f"Required file '{results_file}' not found in output directory: {out_dir}"
        )

    # Process all input files and collect metadata
    all_metadata: list[AutoMetadata] = []
    all_by_key: dict[str, dict[str, Any]] = {}
    all_archive_objs: list[dict[str, Any]] = []
    all_network_objs: list[dict[str, Any]] = []
    all_network_keys: list[str] = []
    modes: set[str] = set()

    for input_path in input_paths:
        resolved_path = input_path.resolve()

        by_key, archive_obj, network_obj = _load_token_table(resolved_path)
        mode = _detect_mode(by_key, archive_obj, network_obj)
        modes.add(mode)

        network_key = _get_network_key(archive_obj, network_obj)
        all_network_keys.append(network_key)

        archive_stem = resolved_path.name
        for suffix in (".json.bz2", ".bz2", ".json"):
            if archive_stem.endswith(suffix):
                archive_stem = archive_stem[: -len(suffix)]
                break

        metadata = _extract_auto_metadata(
            by_key=by_key,
            mode=mode,
            archive_path=resolved_path,
            network_key=network_key,
            archive_stem=archive_stem,
        )
        all_metadata.append(metadata)

        # Accumulate objects for summary generation
        all_by_key.update(by_key)
        if archive_obj:
            all_archive_objs.append(archive_obj)
        if network_obj:
            all_network_objs.append(network_obj)

    # Check consistency
    if len(modes) > 1:
        raise ValueError(
            f"Mixed modes detected across input files: {modes}. All files must be either ASFE or RBFE."
        )

    mode = modes.pop()

    # Merge metadata from all files
    merged_metadata = AutoMetadata()

    # Collect protocol settings from all files
    for metadata in all_metadata:
        merged_metadata.protocol_settings_list.extend(metadata.protocol_settings_list)

        # Use first non-empty value for scalar fields
        if not merged_metadata.openfe_version and metadata.openfe_version:
            merged_metadata.openfe_version = metadata.openfe_version
        if not merged_metadata.openmm_version and metadata.openmm_version:
            merged_metadata.openmm_version = metadata.openmm_version
        if (
            not merged_metadata.openff_toolkit_version
            and metadata.openff_toolkit_version
        ):
            merged_metadata.openff_toolkit_version = metadata.openff_toolkit_version
        if not merged_metadata.forcefield and metadata.forcefield:
            merged_metadata.forcefield = metadata.forcefield
        if not merged_metadata.partial_charges and metadata.partial_charges:
            merged_metadata.partial_charges = metadata.partial_charges
        if not merged_metadata.benchmark_data_set and metadata.benchmark_data_set:
            merged_metadata.benchmark_data_set = metadata.benchmark_data_set
        if not merged_metadata.benchmark_system and metadata.benchmark_system:
            merged_metadata.benchmark_system = metadata.benchmark_system

    # Use merged data for outputs
    openfe_version = merged_metadata.openfe_version
    openmm_version = merged_metadata.openmm_version
    openff_toolkit_version = merged_metadata.openff_toolkit_version
    forcefield = merged_metadata.forcefield
    partial_charges_raw = merged_metadata.partial_charges
    partial_charge_tag = _normalize_partial_charge_info(partial_charges_raw)
    partial_charges = partial_charge_tag or partial_charges_raw

    # Build content summary from combined data
    # Get list of source file names
    source_file_names = [p.name for p in input_paths]

    content_summary, system_info_list = _build_content_summary(
        all_by_key,
        all_archive_objs,
        all_network_objs,
        mode,
        merged_metadata.benchmark_data_set,
        forcefield,
        partial_charges,
        used_alchemiscale,
        source_file_names,
    )

    # Append additional summary text if provided
    if summary_suffix:
        content_summary = content_summary.rstrip() + " " + summary_suffix.strip()

    # Store system_info_list in merged_metadata
    merged_metadata.system_info_list = system_info_list

    # Build network_key to systems mapping (only for AlchemicalArchive files)
    # For AlchemicalNetwork files, the "name" field is not a meaningful network key
    network_key_to_systems: dict[str, list[str]] = defaultdict(list)
    has_archive_objects = len(all_archive_objs) > 0

    # Override benchmark_data_set and benchmark_system from system_info_list if available
    if system_info_list:
        # Use system_group from system info (more reliable than string matching)
        system_groups = set(
            si.system_group for si in system_info_list if si.system_group
        )
        if len(system_groups) == 1:
            merged_metadata.benchmark_data_set = system_groups.pop()
        elif len(system_groups) > 1:
            # Multiple groups - list them
            merged_metadata.benchmark_data_set = ", ".join(sorted(system_groups))

        # List all system names
        system_names = sorted(
            set(si.system_name for si in system_info_list if si.system_name)
        )
        if system_names:
            merged_metadata.benchmark_system = ", ".join(system_names)

        # Update protocol_settings_list with correct system info
        # For each input file, extract system info from first transformation
        file_to_system: dict[str, tuple[str, str]] = {}

        for idx, input_path in enumerate(input_paths):
            by_key, archive_obj, network_obj = _load_token_table(input_path)
            transformation_refs = _transformation_refs(archive_obj, network_obj)
            if transformation_refs:
                system_group, system_name = _extract_system_info_from_mapping(
                    by_key, transformation_refs[0]
                )
                file_to_system[str(input_path)] = (system_group, system_name)

                # Build network_key to systems mapping only for AlchemicalArchive files
                if has_archive_objects and archive_obj:
                    network_key = all_network_keys[idx]
                    system_path = (
                        f"{system_group}/{system_name}"
                        if system_group and system_name
                        else system_name
                    )
                    if (
                        system_path
                        and system_path not in network_key_to_systems[network_key]
                    ):
                        network_key_to_systems[network_key].append(system_path)

        # Update each protocol settings info with correct system data
        for protocol_info in merged_metadata.protocol_settings_list:
            # protocol_info.source_file is the full path
            if protocol_info.source_file in file_to_system:
                group, name = file_to_system[protocol_info.source_file]
                protocol_info.benchmark_set = group
                protocol_info.benchmark_system = name

    # Extract these AFTER overriding from system_info_list
    benchmark_data_set = merged_metadata.benchmark_data_set
    benchmark_system = merged_metadata.benchmark_system

    # Calculate total transformations across all files
    total_transformations = 0
    for archive_obj in all_archive_objs:
        total_transformations += len(_transformation_refs(archive_obj, None))
    for network_obj in all_network_objs:
        if not any(
            net_obj == network_obj
            for net_obj in [a.get("network") for a in all_archive_objs if a]
        ):
            total_transformations += len(_transformation_refs(None, network_obj))

    # Use primary network key for submission ID
    primary_network_key = all_network_keys[0] if all_network_keys else "unknown"

    # Generate a descriptive title
    submission_id_str = submission_id or _default_submission_id(primary_network_key)

    # Extract unique benchmark sets and systems for title
    unique_sets = []
    if system_info_list:
        unique_sets = sorted(
            set(si.system_group for si in system_info_list if si.system_group)
        )
    if not unique_sets and benchmark_data_set:
        unique_sets = [s.strip() for s in benchmark_data_set.split(",")]

    unique_systems = []
    if system_info_list:
        unique_systems = sorted(
            set(si.system_name for si in system_info_list if si.system_name)
        )
    elif benchmark_system:
        unique_systems = [s.strip() for s in benchmark_system.split(",")]

    title = _generate_title(mode, unique_sets, unique_systems, submission_id_str)

    submission_yaml_filename = "submission.yaml"
    zenodo_description_filename = "zenodo_description.md"

    submission_yaml_path = out_dir / submission_yaml_filename
    zenodo_description_path = out_dir / zenodo_description_filename

    submission_id = submission_id or _default_submission_id(primary_network_key)
    tags_list = [k.strip() for k in tags.split(",") if k.strip()]
    tags_final = _make_tags(
        mode=mode,
        forcefield=forcefield,
        partial_charge_tag=partial_charges,
        user_keywords=tags_list,
    )

    submission_yaml_text = _make_submission_yaml(
        submission_id=submission_id,
        title=title,
        summary=content_summary,
        tags=tags_final,
        authors=author or [],
        openfe_version=openfe_version,
        openmm_version=openmm_version,
        openff_toolkit_version=openff_toolkit_version,
        forcefield=forcefield,
        partial_charges=partial_charges,
        benchmark_data_set=benchmark_data_set,
        benchmark_system=benchmark_system,
        archive_doi="TODO add DOI",
        archive_provider="TODO add archive provider",
        license_name=license,
        protocol_settings_list=merged_metadata.protocol_settings_list,
        network_key_to_systems=network_key_to_systems,
        results_file=results_file,
    )
    submission_yaml_path.write_text(submission_yaml_text)

    # For Zenodo description, list all input files
    archive_filenames = ", ".join(p.name for p in input_paths)

    zenodo_description_text = _make_zenodo_description(
        title=title,
        archive_filename=archive_filenames,
        mode=mode,
        content_summary=content_summary,
        openfe_version=openfe_version,
        openmm_version=openmm_version,
        openff_toolkit_version=openff_toolkit_version,
        forcefield=forcefield,
        partial_charges=partial_charges,
        benchmark_data_set=benchmark_data_set,
        benchmark_system=benchmark_system,
        license_name=license,
        protocol_settings_list=merged_metadata.protocol_settings_list,
        has_archive_objects=has_archive_objects,
        used_alchemiscale=used_alchemiscale,
        system_info_list=system_info_list,
        network_key_to_systems=network_key_to_systems,
    )
    zenodo_description_path.write_text(zenodo_description_text)

    print(f"Processed {len(input_paths)} input file(s)")
    print(f"Detected mode: {mode}")
    print(f"Submission YAML: {submission_yaml_path}")
    print(f"Zenodo description: {zenodo_description_path}")

    return submission_yaml_path, zenodo_description_path


def main():
    """CLI entry point for prepare_metadata_submission."""
    parser = argparse.ArgumentParser(
        description="Generate submission.yaml and zenodo_description.md from OpenFE JSON archives",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
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
        """),
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
        default="openfe,alchemicalarchive",
        help="Comma-separated tags (default: 'openfe,alchemicalarchive')",
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

    args = parser.parse_args()

    # Expand glob patterns and collect all matching files
    all_files: list[Path] = []
    for pattern in args.input_patterns:
        matched = glob_module.glob(pattern, recursive=True)
        if matched:
            all_files.extend(Path(f) for f in sorted(matched))
        else:
            # If no glob match, treat as literal file path
            all_files.append(Path(pattern))

    if not all_files:
        print("Error: No input files found", file=sys.stderr)
        return 1

    # Remove duplicates while preserving order
    seen = set()
    unique_files = []
    for f in all_files:
        if f not in seen:
            seen.add(f)
            unique_files.append(f)

    try:
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
        )
        print("\n✓ Successfully generated submission metadata")
        return 0
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
