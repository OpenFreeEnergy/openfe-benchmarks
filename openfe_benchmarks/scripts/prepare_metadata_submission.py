#!/usr/bin/env python3
"""Prepare benchmark submission artifacts from an AlchemicalArchive or AlchemicalNetwork JSON file.

This module generates `submission.yaml` and `zenodo_description.md` from a JSON archive
exported by OpenFE/Alchemiscale. It no longer exposes a CLI; instead, call
`process_network(...)` directly from Python.

Example:

    from pathlib import Path
    from openfe_benchmarks.scripts.prepare_metadata_submission import process_network

    process_network(
        input_file=Path("archive.json.bz2"),
        output_dir=Path("."),
        submission_id="2026-04-15-example",
        keywords="openfe,alchemicalarchive",
        author=["Jane Doe"],
        license="CC-BY-4.0",
    )
"""

from __future__ import annotations

import bz2
import json
import re
import statistics
import textwrap
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any


@dataclass
class AutoMetadata:
    openfe_version: str = ""
    openmm_version: str = ""
    openff_toolkit_version: str = ""
    forcefield: str = ""
    partial_charges: str = ""
    network_descriptor: str = ""
    benchmark_data_set: str = ""
    benchmark_system: str = ""
    protocol_settings: dict[str, str] = field(default_factory=dict)


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


def _guess_network_descriptor(archive_stem: str, network_key: str) -> str:
    if archive_stem.startswith(f"{network_key}-"):
        return archive_stem[len(network_key) + 1 :]
    if archive_stem.startswith("AlchemicalNetwork-"):
        parts = archive_stem.split("-", maxsplit=2)
        if len(parts) == 3:
            return parts[2]
    return archive_stem


def _infer_benchmark_data_set_system(
    *, by_key: dict[str, dict[str, Any]], mode: str, archive_stem: str, network_key: str
) -> tuple[str, str]:
    blob = json.dumps(list(by_key.values())).lower()
    descriptor = _guess_network_descriptor(archive_stem, network_key).lower()
    search_space = " ".join(
        [blob, descriptor, archive_stem.lower(), network_key.lower()]
    )

    system = ""
    for candidate in ("freesolv", "tyk2", "mnsol"):
        if candidate in search_space:
            system = candidate
            break

    benchmark_set = ""
    if "jacs_set" in search_space or "jacs" in search_space:
        benchmark_set = "jacs_set"
    elif "solvation_set" in search_space or mode == "asfe":
        benchmark_set = "solvation_set"

    if not system and mode == "asfe":
        system = "freesolv"

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
    protocol_name = (
        "AbsoluteSolvationProtocol"
        if mode == "asfe"
        else "RelativeHybridTopologyProtocol"
    )
    out: dict[str, str] = {"protocol": protocol_name}

    if not protocol_obj:
        out["notes"] = "Protocol settings unavailable in archive payload."
        return out

    settings = protocol_obj.get("settings") or {}

    repeats = settings.get("protocol_repeats")
    if repeats is not None:
        out["repeats"] = str(repeats)

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
    else:
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

    if len(out) == 1:
        out["notes"] = "Protocol class found, but detailed settings were unavailable."

    return out


def _render_protocol_settings_yaml(protocol_settings: dict[str, str]) -> str:
    lines = ["protocol_settings:"]
    protocol_name = protocol_settings.get("protocol", "")
    lines.append(f"  - protocol: {protocol_name}")

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
        "repeats",
        "lambda_windows",
        "lambda_schedule",
        "notes",
    ]

    for key in preferred_order:
        if key in protocol_settings:
            lines.append(f"    {key}: {json.dumps(str(protocol_settings[key]))}")

    for key in sorted(
        k for k in protocol_settings if k not in set(preferred_order) | {"protocol"}
    ):
        lines.append(f"    {key}: {json.dumps(str(protocol_settings[key]))}")

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


def _repeat_stats_summary(repeat_counts: list[int]) -> str:
    if not repeat_counts:
        return "repeats stats unavailable"

    mean_repeats = statistics.fmean(repeat_counts)
    median_repeats = statistics.median(repeat_counts)
    dist: dict[int, int] = {}
    for count in repeat_counts:
        dist[count] = dist.get(count, 0) + 1
    dist_text = ", ".join(f"{k}:{v}" for k, v in sorted(dist.items()))
    return (
        f"min={min(repeat_counts)}, median={median_repeats:g}, mean={mean_repeats:.2f}, "
        f"max={max(repeat_counts)}, distribution={{ {dist_text} }}"
    )


def _build_content_summary(
    by_key: dict[str, dict[str, Any]],
    archive_obj: dict[str, Any] | None,
    network_obj: dict[str, Any] | None,
    mode: str,
    benchmark_data_set: str,
    forcefield: str,
    partial_charges: str,
) -> str:
    transformation_results = (
        archive_obj.get("transformation_results", []) if archive_obj else []
    )
    repeat_counts: list[int] = []
    max_atoms_per_system = 0

    solutes: set[str] = set()
    solvents: set[str] = set()

    ligands: set[str] = set()
    proteins: set[str] = set()
    cofactors: set[str] = set()
    systems_with_cofactors: set[str] = set()

    visited_systems_for_cofactors: set[str] = set()

    if archive_obj is not None:
        transformation_refs = [
            item[0]
            for item in transformation_results
            if isinstance(item, list) and len(item) == 2
        ]
        repeat_counts = [
            len(item[1])
            for item in transformation_results
            if isinstance(item, list) and len(item) == 2
        ]
    else:
        transformation_refs = _transformation_refs(archive_obj, network_obj)

    transformation_count = len(transformation_refs)

    for transformation_ref in transformation_refs:
        _, transformation = _resolve_payload(by_key, transformation_ref)
        if not transformation:
            continue

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
                        solvents.add(comp_name)
                    elif "solute" in label_l or qualname == "SmallMoleculeComponent":
                        solutes.add(comp_name)
                else:
                    if "protein" in label_l or qualname == "ProteinComponent":
                        proteins.add(comp_name)
                    elif "ligand" in label_l:
                        ligands.add(comp_name)
                    elif "cofactor" in label_l:
                        cofactors.add(comp_name)
                        local_cofactors.add(comp_name)
                    elif (
                        qualname == "SmallMoleculeComponent"
                        and "solvent" not in label_l
                    ):
                        # Non-solvent small molecules that are not explicit ligands are treated as cofactors.
                        cofactors.add(comp_name)
                        local_cofactors.add(comp_name)

            max_atoms_per_system = max(max_atoms_per_system, system_atoms)

            if (
                mode == "rbfe"
                and cs_key
                and cs_key not in visited_systems_for_cofactors
            ):
                visited_systems_for_cofactors.add(cs_key)
                if local_cofactors:
                    systems_with_cofactors.add(cs_key)

    repeats_text = _repeat_stats_summary(repeat_counts)

    subject = benchmark_data_set or "benchmark"
    field_info = forcefield or "an unspecified force field"
    charge_info = partial_charges or "unspecified partial charges"
    if mode == "rbfe":
        cofactor_list = ", ".join(sorted(cofactors)) if cofactors else "none"
        summary_parts = [
            f"This submission describes the {subject} RBFE benchmark prepared with {field_info} and {charge_info}.",
            f"The network contains {transformation_count} transformations across {len(ligands)} unique ligands and {len(proteins)} unique proteins.",
        ]
        if systems_with_cofactors:
            summary_parts.append(
                f"{len(systems_with_cofactors)} systems include cofactors ({cofactor_list})."
            )
    else:
        summary_parts = [
            f"This submission describes the {subject} ASFE benchmark prepared with {field_info} and {charge_info}.",
            f"The archive contains {transformation_count} transformations across {len(solutes)} unique solutes and {len(solvents)} unique solvents.",
        ]

    summary_parts.append(
        f"The largest simulated chemical system contains {max_atoms_per_system} atoms, and repeat counts per transformation are {repeats_text}."
    )
    summary_parts.append(
        "Results are derived from archived Alchemiscale workflow data."
    )
    summary_text = " ".join(summary_parts)
    return textwrap.fill(summary_text, width=100)


def _extract_auto_metadata(
    *,
    by_key: dict[str, dict[str, Any]],
    mode: str,
    archive_path: Path,
    network_key: str,
    archive_stem: str,
) -> AutoMetadata:
    metadata = AutoMetadata()
    metadata.network_descriptor = _guess_network_descriptor(archive_stem, network_key)
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

    metadata.protocol_settings = _build_protocol_settings(protocol_obj, mode)
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
    network_descriptor: str,
    benchmark_data_set: str,
    benchmark_system: str,
    archive_doi: str,
    archive_provider: str,
    license_name: str,
    protocol_settings: dict[str, str],
) -> str:
    if not authors:
        authors = ["TODO: add author name"]

    tags_yaml = ", ".join(tags)
    authors_yaml = "\n".join(f"  - name: {name}" for name in authors)
    protocol_settings_yaml = _render_protocol_settings_yaml(protocol_settings)

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
network: {network_descriptor}

# BenchmarkData provenance (from openfe-benchmarks planning script)
benchmark_data:
  source_repository: https://github.com/OpenFreeEnergy/openfe-benchmarks
  set: {benchmark_data_set}
  system: {benchmark_system}

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
    summary: str,
    archive_filename: str,
    network_key: str,
    mode: str,
    content_summary: str,
    openfe_version: str,
    openmm_version: str,
    openff_toolkit_version: str,
    forcefield: str,
    partial_charges: str,
    network_descriptor: str,
    tags: list[str],
    benchmark_data_set: str,
    benchmark_system: str,
    n_transformations: int,
    submission_yaml_file: str,
    license_name: str,
    protocol_settings: dict[str, str],
) -> str:
    content_kind = "ASFE" if mode == "asfe" else "RBFE"
    openmm_display = openmm_version or "<not found in archive>"
    protocol_lines = (
        "\n".join(
            [f"- {k}: {v}" for k, v in protocol_settings.items() if k != "protocol"]
        )
        or "- notes: Protocol settings unavailable"
    )

    return f"""# {title}

## Description

## Overview

{content_kind} benchmark results prepared from an AlchemicalArchive generated with OpenFE and Alchemiscale.

{content_summary}

## Software versions

- openfe_version: {openfe_version}
- openmm_version: {openmm_display}
- openff_toolkit_version: {openff_toolkit_version}

## Recommended descriptors

- forcefield: {forcefield}
- partial_charges: {partial_charges}
- network: {network_descriptor}

## BenchmarkData provenance

- source_repository: https://github.com/OpenFreeEnergy/openfe-benchmarks
- set: {benchmark_data_set}
- system: {benchmark_system}

## Protocol settings

- protocol: {protocol_settings.get("protocol", "unknown")}
{protocol_lines}

- archive file: {archive_filename}
- network key: {network_key}

## Contents

### Data files

- {submission_yaml_file}: submission metadata
- zenodo_description.md: Zenodo metadata description

### Network summary

- total transformations in archive: {n_transformations}

## Simulation details

- generation workflow: network prepared with plan_*.py scripts from openfe-benchmarks and then archived from Alchemiscale

## Changelog

- Generated by prepare_archive_submission.py on {date.today().isoformat()}

## Rights

- License: {license_name}
"""


def process_network(
    input_file: Path,
    output_dir: Path = Path("."),
    submission_id: str | None = None,
    keywords: str = "openfe,alchemicalarchive",
    author: list[str] | None = None,
    license: str = "CC-BY-4.0",
) -> tuple[Path, Path]:
    """Generate submission metadata from an archived OpenFE JSON network.

    Parameters
    ----------
    input_file:
        Path to the AlchemicalArchive or AlchemicalNetwork JSON file. Supported
        extensions are `.json`, `.bz2`, and `.json.bz2`.
    output_dir:
        Directory where `submission.yaml` and `zenodo_description.md` will be
        written. Defaults to the current working directory.
    submission_id:
        Optional identifier to use in `submission.yaml`. If omitted, a default
        value is generated from the current date and network key.
    keywords:
        Comma-separated list of additional tags to include in the submission
        metadata. The generated tag list also always includes the detected
        `mode` (either ``asfe`` or ``rbfe``), the resolved forcefield string,
        and normalized partial charge information.
    author:
        Optional list of author entries for the submission YAML. Each entry is
        treated as a raw string and written to the `authors` section.
    license:
        License string to write into the submission metadata.

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
    input_path = input_file.resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    out_dir = output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    by_key, archive_obj, network_obj = _load_token_table(input_path)
    mode = _detect_mode(by_key, archive_obj, network_obj)

    network_key = _get_network_key(archive_obj, network_obj)

    archive_stem = input_path.name
    for suffix in (".json.bz2", ".bz2", ".json"):
        if archive_stem.endswith(suffix):
            archive_stem = archive_stem[: -len(suffix)]
            break

    auto_metadata = _extract_auto_metadata(
        by_key=by_key,
        mode=mode,
        archive_path=input_path,
        network_key=network_key,
        archive_stem=archive_stem,
    )

    openfe_version = auto_metadata.openfe_version
    openmm_version = auto_metadata.openmm_version
    openff_toolkit_version = auto_metadata.openff_toolkit_version
    forcefield = auto_metadata.forcefield
    partial_charges_raw = auto_metadata.partial_charges
    partial_charge_tag = _normalize_partial_charge_info(partial_charges_raw)
    partial_charges = partial_charge_tag or partial_charges_raw
    network_descriptor = auto_metadata.network_descriptor

    content_summary = _build_content_summary(
        by_key,
        archive_obj,
        network_obj,
        mode,
        auto_metadata.benchmark_data_set,
        forcefield,
        partial_charges,
    )

    benchmark_data_set = auto_metadata.benchmark_data_set
    benchmark_system = auto_metadata.benchmark_system

    # Requested behavior: always leave title for manual curation.
    title = "TODO: add title"

    submission_yaml_filename = "submission.yaml"
    zenodo_description_filename = "zenodo_description.md"

    submission_yaml_path = out_dir / submission_yaml_filename
    zenodo_description_path = out_dir / zenodo_description_filename

    submission_id = submission_id or _default_submission_id(network_key)
    keywords_list = [k.strip() for k in keywords.split(",") if k.strip()]
    tags = _make_tags(
        mode=mode,
        forcefield=forcefield,
        partial_charge_tag=partial_charges,
        user_keywords=keywords_list,
    )

    submission_yaml_text = _make_submission_yaml(
        submission_id=submission_id,
        title=title,
        summary=content_summary,
        tags=tags,
        authors=author or [],
        openfe_version=openfe_version,
        openmm_version=openmm_version,
        openff_toolkit_version=openff_toolkit_version,
        forcefield=forcefield,
        partial_charges=partial_charges,
        network_descriptor=network_descriptor,
        benchmark_data_set=benchmark_data_set,
        benchmark_system=benchmark_system,
        archive_doi="TODO: add DOI",
        archive_provider="TODO: add archive provider",
        license_name=license,
        protocol_settings=auto_metadata.protocol_settings,
    )
    submission_yaml_path.write_text(submission_yaml_text)

    zenodo_description_text = _make_zenodo_description(
        title=title,
        summary=content_summary,
        archive_filename=input_path.name,
        network_key=network_key,
        mode=mode,
        content_summary=content_summary,
        openfe_version=openfe_version,
        openmm_version=openmm_version,
        openff_toolkit_version=openff_toolkit_version,
        forcefield=forcefield,
        partial_charges=partial_charges,
        network_descriptor=network_descriptor,
        tags=tags,
        benchmark_data_set=benchmark_data_set,
        benchmark_system=benchmark_system,
        n_transformations=len(_transformation_refs(archive_obj, network_obj)),
        submission_yaml_file=submission_yaml_filename,
        license_name=license,
        protocol_settings=auto_metadata.protocol_settings,
    )
    zenodo_description_path.write_text(zenodo_description_text)

    print(f"Input file: {input_path}")
    print(f"Detected mode: {mode}")
    print(f"Submission YAML: {submission_yaml_path}")
    print(f"Zenodo description: {zenodo_description_path}")

    return submission_yaml_path, zenodo_description_path
