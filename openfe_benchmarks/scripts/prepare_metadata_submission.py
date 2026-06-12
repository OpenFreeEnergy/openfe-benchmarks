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

import os
import argparse
import ast
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
import warnings
import pprint

from pint import UnitRegistry

from gufe.archival import AlchemicalArchive
from gufe import AlchemicalNetwork
from gufe.transformations.transformation import Transformation

from openfe_benchmarks.data import BenchmarkIndex

ureg = UnitRegistry()


def _add_value_with_keys(
    list_obj: list[tuple[Any, list[str]]],
    value: Any,
    keys: list[str],
) -> None:
    for existing_value, existing_keys in list_obj:
        if existing_value == value:
            for key in keys:
                if key not in existing_keys:
                    existing_keys.append(key)
            return
    list_obj.append((value, list(keys)))


@dataclass
class ProtocolSettingsInfo:
    """Container for protocol settings with source metadata."""

    calculation_mode: str
    protocol: str
    full_protocol_settings: str
    timestep: str
    temperature: str
    pressure: str
    lambda_functions: str
    small_molecule_forcefield: str
    forcefields: set[str]
    partial_charges: str
    lambda_windows: str = ""
    lambda_schedule: str = ""
    notes: str = ""
    # for rbfe
    equilibration_time: str | None = None
    production_time: str | None = None
    # for asfe
    vacuum_equilibration_time: str | None = None
    vacuum_production_time: str | None = None
    solvent_equilibration_time: str | None = None
    solvent_production_time: str | None = None

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ProtocolSettingsInfo):
            return NotImplemented
        return (
            self.calculation_mode == other.calculation_mode
            and self.protocol == other.protocol
            and self.notes == other.notes
            and self.full_protocol_settings == other.full_protocol_settings
            and self.timestep == other.timestep
            and self.temperature == other.temperature
            and self.pressure == other.pressure
            and self.lambda_functions == other.lambda_functions
            and self.lambda_windows == other.lambda_windows
            and self.lambda_schedule == other.lambda_schedule
            and self.small_molecule_forcefield == other.small_molecule_forcefield
            and self.forcefields == other.forcefields
            and self.partial_charges == other.partial_charges
            and self.equilibration_time == other.equilibration_time
            and self.production_time == other.production_time
            and self.vacuum_equilibration_time == other.vacuum_equilibration_time
            and self.vacuum_production_time == other.vacuum_production_time
            and self.solvent_equilibration_time == other.solvent_equilibration_time
            and self.solvent_production_time == other.solvent_production_time
        )


@dataclass
class SystemInfo:
    """Per-system information extracted from edges."""

    benchmark_set: str
    benchmark_system: str
    calculation_mode: str
    source_file: str
    network_key: str
    ligands: set[str] = field(default_factory=set)
    proteins: set[str] = field(default_factory=set)
    cofactors: set[str] = field(default_factory=set)
    solvents: set[str] = field(default_factory=set)
    files: set[str] = field(default_factory=set)
    openfe_version: list[tuple[str, list[str]]] = field(default_factory=list)
    openmm_version: list[tuple[str, list[str]]] = field(default_factory=list)
    openff_toolkit_version: list[tuple[str, list[str]]] = field(default_factory=list)
    mapper: list[tuple[str, list[str]]] = field(default_factory=list)
    protocol_settings_list: list[tuple[ProtocolSettingsInfo, list[str]]] = field(
        default_factory=list
    )

    def make_key(
        self,
        network_key,
        ligand_start,
        cofactors,
        solvent,
        ligand_final=None,
        protein=None,
    ):
        if self.calculation_mode == "rbfe":
            return f"{network_key} {self.benchmark_set}-{self.benchmark_system}: ligand_start={ligand_start}, ligand_final={ligand_final}, solvent={solvent or 'none'}, cofactors={cofactors or 'none'}, protein={protein or 'none'}"
        elif self.calculation_mode == "asfe":
            if protein is not None or ligand_final is not None:
                warnings.warn("ASFEs do not use final ligand or protein information.")
            return f"{network_key} {self.benchmark_set}-{self.benchmark_system}: ligand_start={ligand_start}, solvent={solvent or 'none'}, cofactors={cofactors or 'none'}"
        else:
            raise ValueError(
                "Set the calculation mode to a supported value: 'rbfe', 'asfe'"
            )

    def add_version_setting(self, attribute, value, key):
        """Add a version attribute to the appropriate list

        Parameters
        ----------
        attribute : str
            Attribute of SystemInfo, one of openmm_version, openfe_version, or openff_toolkit_version
        value : str
            Version string
        key : str
            String representing the calculation run with this version
        """
        _add_value_with_keys(getattr(self, attribute), value, [key])

    def add_protocol_settings(self, protocol_settings: ProtocolSettingsInfo, key):
        """Add or update protocol settings with associated transformation key.

        Stores unique ProtocolSettingsInfo objects with a list of edge keys
        that use that protocol configuration.
        """
        _add_value_with_keys(self.protocol_settings_list, protocol_settings, [key])


@dataclass
class AutoMetadata:
    calculation_mode: str = ""
    network_key: str = ""
    n_transformations: int = 0
    benchmark_sets_systems: list[tuple] = field(default_factory=list)
    system_info_dict: dict[tuple, SystemInfo] = field(default_factory=dict)
    openfe_version: list[tuple[str, list[str]]] = field(default_factory=list)
    openmm_version: list[tuple[str, list[str]]] = field(default_factory=list)
    openff_toolkit_version: list[tuple[str, list[str]]] = field(default_factory=list)
    mapper: list[tuple[str, list[str]]] = field(default_factory=list)
    protocols: list[tuple[str, list[str]]] = field(default_factory=list)
    forcefield: list[tuple[str, list[str]]] = field(default_factory=list)
    small_molecule_force_field: list[tuple[str, list[str]]] = field(
        default_factory=list
    )
    partial_charges: list[tuple[str, list[str]]] = field(default_factory=list)
    protocol_settings_list: list[tuple[ProtocolSettingsInfo, list[str]]] = field(
        default_factory=list
    )

    def update_from_system_info(self) -> None:
        """Update aggregated fields from the contained SystemInfo entries."""
        self.benchmark_sets_systems = list(self.system_info_dict.keys())

        self.openfe_version = []
        self.openmm_version = []
        self.openff_toolkit_version = []
        self.mapper = []
        self.protocols = []
        self.forcefield = []
        self.small_molecule_force_field = []
        self.partial_charges = []
        self.protocol_settings_list = []

        for system_info in self.system_info_dict.values():
            for version, keys in system_info.openfe_version:
                _add_value_with_keys(self.openfe_version, version, keys)
            for version, keys in system_info.openmm_version:
                _add_value_with_keys(self.openmm_version, version, keys)
            for version, keys in system_info.openff_toolkit_version:
                _add_value_with_keys(self.openff_toolkit_version, version, keys)
            for mapper_info, keys in system_info.mapper:
                _add_value_with_keys(self.mapper, mapper_info, keys)

            for protocol_settings, keys in system_info.protocol_settings_list:
                if protocol_settings.protocol:
                    _add_value_with_keys(
                        self.protocols, protocol_settings.protocol, keys
                    )
                if protocol_settings.forcefields:
                    _add_value_with_keys(
                        self.forcefield, protocol_settings.forcefields, keys
                    )
                if protocol_settings.small_molecule_forcefield:
                    _add_value_with_keys(
                        self.small_molecule_force_field,
                        protocol_settings.small_molecule_forcefield,
                        keys,
                    )
                if protocol_settings.partial_charges:
                    _add_value_with_keys(
                        self.partial_charges,
                        protocol_settings.partial_charges,
                        keys,
                    )

                for existing_settings, existing_keys in self.protocol_settings_list:
                    if existing_settings == protocol_settings:
                        for key in keys:
                            if key not in existing_keys:
                                existing_keys.append(key)
                        break
                else:
                    self.protocol_settings_list.append((protocol_settings, list(keys)))


def _load_network(
    input_path: Path,
) -> AlchemicalNetwork | AlchemicalArchive:
    try:
        archive = AlchemicalArchive.from_json(file=input_path)
        alchemical_network = archive.network
        mode = "alchemicalarchive"
    except Exception:
        try:
            alchemical_network = AlchemicalNetwork.from_json(file=input_path)
            mode = "alchemicalnetwork"
        except Exception:
            raise ImportError(
                f"Could not import file neither an AlchemicalArchive nor AlchemicalNetwork: {input_path}"
            )

    return alchemical_network, mode


def _get_network_key(
    network_obj: AlchemicalArchive | AlchemicalNetwork,
    mode: str,
) -> str:
    if mode == "alchemicalarchive":
        return network_obj.network.key
    elif mode == "alchemicalnetwork":
        return network_obj.key
    else:
        raise ValueError(
            f"Network mode must be either 'alchemical network' or 'alchemicalarchive', not {mode},"
        )


def _transformation_refs(
    network_obj: AlchemicalArchive | AlchemicalNetwork,
    mode: str,
) -> list[Any]:
    if mode == "alchemicalarchive":
        return network_obj.transformation_results
    elif mode == "alchemicalnetwork":
        return network_obj.edges
    else:
        raise ValueError(
            f"Network mode must be either 'alchemical network' or 'alchemicalarchive', not {mode},"
        )


def _detect_calc_mode(
    network_obj: AlchemicalArchive | AlchemicalNetwork,
    mode: str,
) -> str:
    names: list[str] = []

    names = [trans.name for trans in _transformation_refs(network_obj, mode)]
    if names:
        # !!!! NoteHere !!! what is hard coded to detect calculation type?
        if any(n.startswith("complex_") or n.startswith("solvent_") for n in names):
            return "rbfe"
        return "asfe"


def _slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")


def _default_submission_id(network_key: str) -> str:
    return f"{date.today().isoformat()}-{_slugify(network_key)}"


def _generate_title(
    mode: str,
    benchmark_set_systems: list[tuple(str, str)],
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
    n_set_systems = len(benchmark_set_systems)

    if n_set_systems == 0:
        # Fallback if no benchmark set detected
        return f"OpenFE {mode} Benchmark - {submission_id}"

    if len(set([x[0] for x in benchmark_set_systems])) == 1:
        set_name = benchmark_set_systems[0][0]
        if n_set_systems <= 3:
            # List system names
            systems_str = ", ".join([x[1] for x in benchmark_set_systems])
            return f"OpenFE {mode} - {set_name} - {systems_str} - {submission_id}"
        else:
            # Use count
            return f"OpenFE {mode} - {set_name} ({n_set_systems} systems) - {submission_id}"

    if n_set_systems <= 3:
        return f"OpenFE {mode} - {', '.join([f'{x}/{y}' for x, y in benchmark_set_systems])} - {submission_id}"

    # Many sets - use multi-set notation
    return f"OpenFE {mode} - Multi-set Benchmark ({len(set([x[0] for x in benchmark_set_systems]))} sets, {len(set([x[1] for x in benchmark_set_systems]))} systems) - {submission_id}"


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
    # For pint quantities that have been processed by pydantic model_dump()
    if isinstance(value, dict) and "unit" in value:
        q = value["val"] * ureg.parse_expression(value["unit"])
        return f"{q:#~}"
    else:
        return str(value)


def _infer_benchmark_data_set_system(
    trans: Transformation,
) -> tuple[str, str]:
    """Infer benchmark set and system from Transformation contents using BenchmarkIndex.

    This searches for any known benchmark set or system name in the transformation mapping metadata.

    Returns:
        (benchmark_set, system_name) tuple, or ("", "") if not found
    """
    benchmark_set = trans.mapping.annotations.get("system_group", None)
    system = trans.mapping.annotations.get("system_name", None)

    if benchmark_set is None and system is not None:
        # Get all known benchmark sets and systems from the index
        index = BenchmarkIndex()
        benchmark_sets_systems = index.list_systems_by_tag()
        benchmark_set = [x[0] for x in benchmark_sets_systems if x[1] == system]
        if benchmark_set:
            benchmark_set = benchmark_set[0]  # just take the first one

    if system is None or benchmark_set is None:
        raise ValueError(
            f"Benchmark set / system could not be found for the transformation, {trans.name}. The following was found: {benchmark_set} / {system}. See valid combinations with `index = ofebm.BenchmarkIndex(); index.list_systems_by_tag()`"
        )

    return benchmark_set, system


def _extract_sim_times(settings_block: dict[str, Any]) -> tuple[str, str]:
    equilibration = settings_block.get("equilibration_length")
    production = settings_block.get("production_length")
    return _quantity_to_text(
        equilibration
    ) if equilibration is not None else "", _quantity_to_text(
        production
    ) if production is not None else ""


def _build_protocol_settings(protocol_obj, calc_mode) -> dict[str, str | set(str)]:
    if not protocol_obj:
        return {
            "protocol": "unknown",
            "notes": "Protocol settings unavailable in archive.",
        }

    # Detect protocol name from the object
    protocol_name = str(type(protocol_obj)).rstrip("'>").split(".")[-1]
    out: dict[str, str] = {"protocol": protocol_name, "calculation_mode": calc_mode}
    settings = protocol_obj.settings.model_dump()

    if not settings:
        out["notes"] = "Protocol class found, but detailed settings were unavailable."
    else:
        out["full_protocol_settings"] = pprint.pformat(settings)

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
        out["lambda_functions"] = lambda_settings.get("lambda_functions", "")
        lambda_windows = lambda_settings.get("lambda_windows")
        if lambda_windows is not None:
            out["lambda_windows"] = str(lambda_windows)
        else:
            lambda_counts: list[str] = []
            for lambda_key, values in lambda_settings.items():
                if lambda_key not in ["lambda_functions", "lambda_windows"]:
                    continue
                lambda_counts.append(f"{lambda_key}:{len(values)}")
            if lambda_counts:
                out["lambda_schedule"] = ", ".join(lambda_counts)

    forcefield_settings = (
        settings.get("forcefield_settings")
        or settings.get("solvent_forcefield_settings")
        or settings.get("vacuum_forcefield_settings")
        or {}
    )
    if forcefield_settings:
        out["small_molecule_forcefield"] = str(
            forcefield_settings.get("small_molecule_forcefield") or ""
        )
        ffs = forcefield_settings.get("forcefields")
        if isinstance(ffs, list) and ffs:
            out["forcefields"] = set(
                sorted([os.path.splitext(ff.split("/")[1])[0] for ff in ffs])
            )

    partial_charge_settings = settings.get("partial_charge_settings") or {}
    if partial_charge_settings:
        out["partial_charges"] = _normalize_partial_charge_info(partial_charge_settings)

    # Protocol-specific handling: RBFE typically has a single simulation block;
    # ASFE commonly has separate vacuum and solvent simulation settings.
    if calc_mode == "rbfe":
        sim = settings.get("simulation_settings") or {}
        if isinstance(sim, dict):
            eq, prod = _extract_sim_times(sim)
            if eq:
                out["equilibration_time"] = eq
            if prod:
                out["production_time"] = prod
    elif calc_mode == "asfe":
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
            f"Calculation type {calc_mode} is not yet supported. Add capability to `_build_protocol_settings`"
        )

    return out


def _component_name(component) -> str:
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

    return "unknown"


def _get_system_info(trans, calc_mode) -> dict[str, set | list]:
    # Per-system tracking
    for state_key in ("stateA", "stateB"):
        chemical_system = getattr(trans, state_key)
        if not chemical_system:
            continue

        solvents = set()
        proteins = set()
        ligands = list()
        cofactors = set()
        for label, component in chemical_system.components.items():
            qualname = str(type(component)).rstrip("'>").split(".")[-1]

            component = component.to_dict()
            comp_name = _component_name(component)

            if calc_mode == "asfe":
                if "solvent" in label or "solventcomponent" in qualname.lower():
                    solvents.add(comp_name)
                elif "solute" in label or qualname == "SmallMoleculeComponent":
                    ligands.append(comp_name)
            elif calc_mode == "rbfe":
                if "protein" in label or qualname == "ProteinComponent":
                    proteins.add(comp_name)
                elif "ligand" in label:
                    if comp_name not in ligands:
                        ligands.append(comp_name)
                elif "cofactor" in label:
                    cofactors.add(comp_name)
                elif qualname == "SmallMoleculeComponent" and "solvent" not in label:
                    # Non-solvent small molecules that are not explicit ligands are treated as cofactors.
                    cofactors.add(comp_name)
            else:
                ValueError(
                    f"Calculation type {calc_mode} is not yet supported. Add capability to `_build_content_summary`"
                )

        return {
            "solvents": solvents,
            "proteins": proteins,
            "ligands": ligands,
            "cofactors": cofactors,
        }


def _extract_auto_metadata(
    network_obj: AlchemicalArchive | AlchemicalNetwork,
    network_mode: str,
    source_file: str,
) -> AutoMetadata:
    metadata = AutoMetadata()
    metadata.network_key = _get_network_key(network_obj, network_mode)
    metadata.calculation_mode = _detect_calc_mode(network_obj, network_mode)

    transformations = _transformation_refs(network_obj, network_mode)
    metadata.n_transformations = len(transformations)
    for trans in transformations:
        benchmark_set_system = _infer_benchmark_data_set_system(trans)
        if benchmark_set_system not in metadata.system_info_dict:
            metadata.system_info_dict[benchmark_set_system] = SystemInfo(
                *benchmark_set_system,
                metadata.calculation_mode,
                source_file,
                metadata.network_key,
            )

        system_info = _get_system_info(trans, metadata.calculation_mode)
        for key, value in system_info.items():
            current_attr = getattr(metadata.system_info_dict[benchmark_set_system], key)
            if isinstance(current_attr, set):
                current_attr.update(value)
            elif isinstance(current_attr, list):
                current_attr.extend(value)
            else:
                raise TypeError(
                    f"Unsupported system_info attribute type for {key}: {type(current_attr)}"
                )

        if metadata.calculation_mode == "rbfe":
            if len(system_info["ligands"]) == 0 or len(system_info["ligands"]) > 2:
                raise ValueError(
                    f"Transformation detects a count other than one or two ligands: network_key={metadata.network_key}, set/system: {benchmark_set_system}, transformation: {trans.name}, ligands: {system_info['ligands']}"
                )
            ligand_start = system_info["ligands"][0]
            ligand_final = (
                system_info["ligands"][1] if len(system_info["ligands"]) > 1 else "none"
            )
        elif metadata.calculation_mode == "asfe":
            if len(system_info["ligands"]) < 1 or len(system_info["ligands"]) > 1:
                raise ValueError(
                    f"Transformation detects a count other than one ligand: network_key={metadata.network_key}, set/system: {benchmark_set_system}, transformation: {trans.name}, ligands: {system_info['ligands']}"
                )
            ligand_start = system_info["ligands"][0]
            ligand_final = "none"

        key = metadata.system_info_dict[benchmark_set_system].make_key(
            metadata.network_key,
            ligand_start,
            system_info["cofactors"],
            system_info["solvents"],
            ligand_final=ligand_final,
            protein=system_info["proteins"],
        )

        protocol_info = ProtocolSettingsInfo(
            **_build_protocol_settings(trans.protocol, metadata.calculation_mode)
        )
        metadata.system_info_dict[benchmark_set_system].add_protocol_settings(
            protocol_info, key
        )

        annotations = trans.mapping.annotations

        # Extract mapper info if available (Option 1: concatenated string)
        if "mapper_settings" in annotations and "mapper_version" in annotations:
            mapper_settings = annotations.get("mapper_settings")
            mapper_version = annotations.get("mapper_version", "unknown")
            if isinstance(mapper_settings, dict):
                mapper_name = mapper_settings.get("__qualname__", "unknown").split(".")[
                    -1
                ]
                mapping_algorithm = mapper_settings.get("_mapping_algorithm", "unknown")
                mapper_str = f"{mapper_name} {mapper_version} ({mapping_algorithm})"
                metadata.system_info_dict[benchmark_set_system].add_version_setting(
                    "mapper", mapper_str, key
                )

        for annotation_key, value in annotations.items():
            if "openmm" in annotation_key:
                metadata.system_info_dict[benchmark_set_system].add_version_setting(
                    "openmm_version", value, annotation_key
                )
            if "openfe" in annotation_key:
                metadata.system_info_dict[benchmark_set_system].add_version_setting(
                    "openfe_version", value, annotation_key
                )
            if "openff" in annotation_key and "toolkit" in annotation_key:
                metadata.system_info_dict[benchmark_set_system].add_version_setting(
                    "openff_toolkit_version", value, annotation_key
                )

    metadata.update_from_system_info()

    return metadata


def _normalize_partial_charge_info(partial_charge_settings: dict) -> str:
    """Normalize partial charge settings to standardized method tags.

    Maps protocol charge method names to the standard method names used in
    openfe_benchmarks.data.data_generation.charge_molecules:
    - am1bcc_at (AM1BCC with AmberTools)
    - am1bcc_oe (AM1BCC with OpenEye)
    - am1bccelf10_oe (AM1BCC ELF10 with OpenEye)
    - nagl_off (NAGL with OpenFF Toolkit)

    For nagl_off, appends the model name if available.

    Parameters
    ----------
    partial_charge_settings : dict
        Protocol partial charge settings dict containing 'partial_charge_method',
        optionally 'off_toolkit_backend', and optionally 'nagl_model'.

    Returns
    -------
    str
        Normalized method tag, e.g., "nagl_off_openff-gnn-am1bcc-1.0.0.pt" or "am1bccelf10_oe".
    """
    if not partial_charge_settings or not isinstance(partial_charge_settings, dict):
        return ""

    method = partial_charge_settings.get("partial_charge_method", "").lower().strip()
    if not method:
        return ""

    # Map method names to standardized tags matching charge_molecules.py
    if "nagl" in method:
        # nagl_off with optional model
        nagl_model = partial_charge_settings.get("nagl_model", "").strip()
        if nagl_model:
            # Extract just the filename if it's a path
            nagl_model = nagl_model.split("/")[-1].split("\\")[-1]
            return f"nagl_off_{nagl_model}"
        return "nagl_off"
    elif "am1bccelf10" in method or "elf10" in method:
        return "am1bccelf10_oe"
    elif "am1bcc" in method:
        # Check toolkit backend to determine if AmberTools or OpenEye
        backend = partial_charge_settings.get("off_toolkit_backend", "").lower().strip()
        if backend == "ambertools":
            return "am1bcc_at"
        elif backend == "openeye":
            return "am1bcc_oe"
        else:
            raise ValueError("Unknown charge backend")
    else:
        # Fallback: normalize any other method name
        normalized = re.sub(r"[^a-z0-9._-]+", "_", method).strip("_")
        return normalized


def _make_tags(
    *,
    mode: str,
    forcefield: list[tuple],
    partial_charge_tag: list[tuple],
    benchmark_data: list[tuple],
    user_keywords: list[str],
) -> list[str]:
    tags: list[str] = []
    tags.append(mode)
    if forcefield:
        tags.extend(list(set(ff for ff_set, _ in forcefield for ff in ff_set)))
    if benchmark_data:
        tags.extend(list(set(y for x in benchmark_data for y in x)))
    if partial_charge_tag:
        tags.extend(list(set(x[0] for x in partial_charge_tag)))
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


def _build_content_summary(
    metadata: AutoMetadata,
    used_alchemiscale: bool = True,
) -> str:
    """
    Build content summary and extract per-system information.

    Parameters
    ----------
    metadata : AutoMetadata
        Consilidated

    Returns:
        (summary_text, list of SystemInfo objects)
    """

    field_info = "/".join(set(ff for ff_set, _ in metadata.forcefield for ff in ff_set))
    if not field_info:
        field_info = "an unspecified force field"

    charge_info = "/".join(set(x[0] for x in metadata.partial_charges))
    if not charge_info:
        charge_info = "an unspecified partial charges"

    # Group systems by benchmark set for explicit listing
    sets_to_systems: dict[str, list[str]] = defaultdict(list)
    for system_group, system_name in metadata.benchmark_sets_systems:
        sets_to_systems[system_group].append(system_name)

    # Sort systems within each set
    for systems_list in sets_to_systems.values():
        systems_list.sort()

    unique_sets = sorted(sets_to_systems.keys())

    # Build descriptive subject line
    if len(unique_sets) == 0:
        subject = "benchmark"
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
        systems_desc = f"{len(metadata.benchmark_sets_systems)} edges"

    # Count totals across all edges
    all_structures = {
        "ligands": set(),
        "proteins": set(),
        "cofactors": set(),
        "solvents": set(),
    }
    systems_with_cofactors = []
    for si in metadata.system_info_dict.values():
        for key in all_structures.keys():
            all_structures[key].update(getattr(si, key))
        if si.cofactors:
            systems_with_cofactors.append(f"{si.benchmark_set}/{si.benchmark_system}")

    # Build summary
    if metadata.calculation_mode == "rbfe":
        cofactor_list = (
            ", ".join(sorted(all_structures["cofactors"]))
            if all_structures["cofactors"]
            else "none"
        )
        if len(unique_sets) > 1:
            summary_parts = [
                f"This submission describes the {subject} RBFE benchmark ({systems_desc}) prepared with {field_info} and {charge_info}.",
                f"The submission contains {metadata.n_transformations} edges, {len(all_structures['ligands'])} unique ligands, and {len(all_structures['proteins'])} unique proteins.",
            ]
        else:
            summary_parts = [
                f"This submission describes the {subject} RBFE benchmark prepared with {field_info} and {charge_info}.",
                f"The network contains {metadata.n_transformations} edges across {len(all_structures['ligands'])} unique ligands and {len(all_structures['proteins'])} unique proteins.",
            ]
        if systems_with_cofactors:
            summary_parts.append(
                f"{len(systems_with_cofactors)} systems include cofactors ({cofactor_list})."
            )
    else:
        if len(unique_sets) > 1:
            summary_parts = [
                f"This submission describes the {subject} ASFE benchmark ({systems_desc}) prepared with {field_info} and {charge_info}.",
                f"The submission contains {metadata.n_transformations} edges, {len(all_structures['ligands'])} unique solutes, and {len(all_structures['solvents'])} unique solvents.",
            ]
        else:
            summary_parts = [
                f"This submission describes the {subject} ASFE benchmark prepared with {field_info} and {charge_info}.",
                f"The archive contains {metadata.n_transformations} edges across {len(all_structures['ligands'])} unique solutes and {len(all_structures['solvents'])} unique solvents.",
            ]

    if used_alchemiscale:
        summary_parts.append(
            "Results are derived from archived Alchemiscale workflow data."
        )
    summary_text = " ".join(summary_parts)
    return textwrap.fill(summary_text, width=100)


def _render_protocol_settings_yaml(
    protocol_settings_list: list[tuple[ProtocolSettingsInfo, list[str]]],
) -> str:
    """Take a list of alchemical protocols pairs with strings identifying systems that use it

    All keys in the ProtocolSettingsInfo class are listed except for ``full_protocol_settings``.

    If only one protocol is used, it is labeled as the submission protocol and the notes specify "Applies to all
    edges". If more than one protocol is present, the protocol that represents the largest number of systems is
    listed last and notes specify as "All remaining edges". The other protocols are listed with notes containing
    the list of identifying strings.

    Parameters
    ----------
    protocol_settings_list : list[tuple[ProtocolSettingsInfo, list[str]]]
        List of unique protocol settings paired with system identifiers that use that protocol.

    Returns
    -------
    str
        Output yaml section
    """
    if not protocol_settings_list:
        return "protocol_settings: []\n"

    def _format_value(value: Any) -> str:
        """Format a value for YAML. Quantity fields (with units) are rendered unquoted."""
        if value is None:
            return ""
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, (int, float)):
            return str(value)
        return json.dumps(str(value))

    def _format_identifier(identifier: Any) -> str:
        if identifier is None:
            return "None"
        if isinstance(identifier, str):
            return identifier
        if isinstance(identifier, (list, tuple, set)):
            return ", ".join(str(item) for item in identifier)
        return str(identifier)

    ordered_settings = sorted(
        enumerate(protocol_settings_list),
        key=lambda item: (len(item[1]), item[0]),
    )

    primary_index = max(
        range(len(protocol_settings_list)),
        key=lambda i: (len(protocol_settings_list[i][1]), -i),
    )
    primary_settings, _ = protocol_settings_list[primary_index]

    def _parse_full_protocol_settings(value: str) -> Any:
        try:
            return ast.literal_eval(value)
        except Exception:
            return None

    def _format_path(path: list[str]) -> str:
        return ".".join(path)

    def _compare_full_protocol_settings(
        base: ProtocolSettingsInfo, other: ProtocolSettingsInfo
    ) -> list[tuple[str, Any, Any]]:
        base_obj = _parse_full_protocol_settings(base.full_protocol_settings)
        other_obj = _parse_full_protocol_settings(other.full_protocol_settings)
        diffs: list[tuple[str, Any, Any]] = []

        def recurse(path: list[str], a: Any, b: Any) -> None:
            if type(a) is not type(b):
                diffs.append((_format_path(path), a, b))
                return
            if isinstance(a, dict):
                for key in sorted(set(a) | set(b)):
                    if key not in a:
                        diffs.append((_format_path(path + [key]), None, b[key]))
                    elif key not in b:
                        diffs.append((_format_path(path + [key]), a[key], None))
                    else:
                        recurse(path + [key], a[key], b[key])
            elif isinstance(a, list):
                if a != b:
                    diffs.append((_format_path(path), a, b))
            else:
                if a != b:
                    diffs.append((_format_path(path), a, b))

        if isinstance(base_obj, dict) and isinstance(other_obj, dict):
            recurse([], base_obj, other_obj)
        return diffs

    def _full_protocol_setting_notes(
        base: ProtocolSettingsInfo, other: ProtocolSettingsInfo
    ) -> list[str]:
        diffs = _compare_full_protocol_settings(base, other)
        if not diffs:
            return []
        notes: list[str] = [
            "Detailed protocol settings differ:",
        ]
        for path, base_value, other_value in diffs:
            notes.append(f"- {path}: {base_value!r} -> {other_value!r}")
        return notes

    output_lines = ["protocol_settings:"]
    multiple_protocols = len(protocol_settings_list) > 1
    field_names = [
        "protocol",
        "timestep",
        "temperature",
        "pressure",
        "forcefields",
        "small_molecule_forcefield",
        "partial_charges",
        "equilibration_time",
        "production_time",
        "vacuum_equilibration_time",
        "vacuum_production_time",
        "solvent_equilibration_time",
        "solvent_production_time",
        "lambda_functions",
        "lambda_windows",
        "lambda_schedule",
        "notes",
    ]

    primary_order = [primary_index] + [
        idx for idx, _ in ordered_settings if idx != primary_index
    ]
    for index in primary_order:
        protocol_settings, identifiers = protocol_settings_list[index]
        is_primary = index == primary_index
        sorted_ids = sorted(_format_identifier(item) for item in identifiers)

        if is_primary:
            if multiple_protocols:
                notes_lines = [f"Applies to {len(sorted_ids)} edges:"] + [
                    f"- {item}" for item in sorted_ids[:5]
                ]
                if len(notes_lines) - 1 != len(sorted_ids):
                    notes_lines.append("- etc.")
                notes = "\n".join(notes_lines)
                notes_is_multiline = True
            else:
                notes = "Applies to all edges"
                notes_is_multiline = False
        else:
            notes_lines = _full_protocol_setting_notes(
                primary_settings, protocol_settings
            )
            trans_lines = [f"Applies to {len(sorted_ids)} edges:"] + [
                f"- {item}" for item in sorted_ids[:5]
            ]
            if len(trans_lines) - 1 != len(sorted_ids):
                trans_lines.append("- etc.")
            notes = "\n".join(notes_lines + trans_lines)
            notes_is_multiline = True

        output_lines.append(
            "  - protocol: " + _format_value(protocol_settings.protocol)
        )
        for field_name in field_names:
            if field_name == "protocol":
                continue
            if field_name == "notes":
                if notes_is_multiline:
                    output_lines.append("    notes: |")
                    for line in notes.splitlines():
                        output_lines.append("      " + line)
                else:
                    output_lines.append("    notes: " + _format_value(notes))
                continue
            # Special handling for forcefields: render as JSON array
            if field_name == "forcefields":
                ff_value = getattr(protocol_settings, field_name)
                if isinstance(ff_value, (list, tuple, set)) and ff_value:
                    items = [json.dumps(str(x)) for x in sorted(ff_value)]
                    if is_primary:
                        output_lines.append(f"    {field_name}: [{', '.join(items)}]")
                continue
            if is_primary:
                output_lines.append(
                    f"    {field_name}: {_format_value(getattr(protocol_settings, field_name))}"
                )
            elif field_name == "protocol":
                output_lines.append(
                    f"    {field_name}: {_format_value(getattr(protocol_settings, field_name))}"
                )

        # For non-primary protocols, only protocol and notes are listed.

    return "\n".join(output_lines) + "\n"


def _render_keyed_values_yaml(
    section_name: str,
    value_keys: list[tuple[Any, list[str]]],
    value_label: str = "value",
    keys_label: str = "edges",
) -> str:
    """Render simple value-with-systems metadata into YAML.

    Parameters
    ----------
    section_name:
        YAML section name.
    value_keys:
        List of (value, system identifiers) pairs.
    value_label:
        Label to use for the scalar value.
    keys_label:
        Label to use for the identifying system keys.

    Returns
    -------
    str
        Output yaml section.
    """
    if not value_keys:
        return f"{section_name}: TODO"
    if len(value_keys) == 1:
        if isinstance(value_keys[0][0], str):
            return f"{section_name}: {json.dumps(str(value_keys[0][0]))}"
        elif isinstance(value_keys[0][0], (list, tuple, set)):
            items = [json.dumps(str(x)) for x in sorted(value_keys[0][0])]
            return f"{section_name}: [{', '.join(items)}]"
        else:
            raise ValueError(f"Unknown value type to print: {value_keys[0][0]}")

    ordered_settings = sorted(
        enumerate(value_keys),
        key=lambda item: (len(item[1]), item[0]),
    )

    lines = [f"{section_name}:"]
    for _, (value, keys) in ordered_settings:
        if isinstance(value, str):
            lines.append(f"  - {value_label}: {json.dumps(str(value))}")
        elif isinstance(value, (list, tuple, set)):
            items = [json.dumps(str(x)) for x in sorted(value)]
            lines.append(f"  - {value_label}: [{', '.join(items)}]")
        else:
            raise ValueError(f"Unknown value type to print: {value_label}")
        if keys:
            lines.append(f"    {keys_label}:")
            for i, key in enumerate(sorted(keys)):
                if i < 5:
                    lines.append(f"      - {json.dumps(str(key))}")
                else:
                    lines.append("      - etc.")
                    break

    return "\n".join(lines)


def _render_benchmark_system_yaml(system_info_dict: dict[tuple, SystemInfo]) -> str:
    """Render a list of SystemInfo objects into a list of network keys and the benchmark systems they contain

    Parameters
    ----------
    system_info_dict : dict[tuple, SystemInfo]
        List of unique protocol settings paired with system identifiers that use that protocol.

    Returns
    -------
    str
        Output yaml section
    """

    network_breakdown = defaultdict(lambda: defaultdict(str))
    for si in system_info_dict.values():
        network_breakdown[si.benchmark_set][si.benchmark_system] = si.network_key

    benchmark_yaml = """
# BenchmarkData provenance (from openfe-benchmarks planning script) with associated network key
benchmark_data:
  source_repository: https://github.com/OpenFreeEnergy/openfe-benchmarks
"""

    for benchmark_set in sorted(network_breakdown):
        benchmark_yaml += f"  {json.dumps(benchmark_set)}:\n"
        for benchmark_system, network_key in sorted(
            network_breakdown[benchmark_set].items()
        ):
            benchmark_yaml += f"    {json.dumps(benchmark_system)}: {network_key}\n"

    return benchmark_yaml


def _make_submission_yaml(
    metadata: AutoMetadata,
    submission_id: str,
    title: str,
    summary: str,
    tags: list[str],
    authors: list[str],
    archive_doi: str,
    archive_provider: str,
    license_name: str,
    results_file: str,
) -> str:
    if not authors:
        authors = ["TODO add author name"]

    tags_yaml = ", ".join(tags)
    authors_yaml = "\n".join(f"  - name: {name}" for name in authors)
    protocol_settings_yaml = _render_protocol_settings_yaml(
        metadata.protocol_settings_list
    )
    benchmark_system_yaml = _render_benchmark_system_yaml(metadata.system_info_dict)
    openfe_version_yaml = _render_keyed_values_yaml(
        "openfe_version", metadata.openfe_version, "version", "edges"
    )
    openmm_version_yaml = _render_keyed_values_yaml(
        "openmm_version", metadata.openmm_version, "version", "edges"
    )
    openff_toolkit_version_yaml = _render_keyed_values_yaml(
        "openff_toolkit_version", metadata.openff_toolkit_version, "version", "edges"
    )
    mapper_yaml = _render_keyed_values_yaml(
        "mapper", metadata.mapper, "mapper", "edges"
    )
    forcefield_yaml = _render_keyed_values_yaml(
        "forcefield", metadata.forcefield, "forcefield", "edges"
    )
    partial_charges_yaml = _render_keyed_values_yaml(
        "partial_charges", metadata.partial_charges, "partial_charges", "edges"
    )

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
{openfe_version_yaml}
{openmm_version_yaml}
{openff_toolkit_version_yaml}
{mapper_yaml}
{forcefield_yaml}
{partial_charges_yaml}
{benchmark_system_yaml}

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
    metadata: AutoMetadata,
    network_mode: str,
    title: str,
    archive_filename: str,
    mode: str,
    content_summary: str,
    license_name: str,
    used_alchemiscale: bool,
) -> str:
    content_kind = "ASFE" if mode == "asfe" else "RBFE"

    # Determine the source type for the overview text
    if network_mode == "alchemicalarchive":
        source_description = "AlchemicalArchive"
    elif network_mode == "alchemicalnetwork":
        source_description = "AlchemicalNetwork"
    else:
        source_description = "OpenFE archive"

    # Build workflow description
    workflow_text = "OpenFE"
    if used_alchemiscale:
        workflow_text += " and Alchemiscale"

    protocol_settings_yaml = _render_protocol_settings_yaml(
        metadata.protocol_settings_list
    )
    benchmark_system_yaml = _render_benchmark_system_yaml(metadata.system_info_dict)
    openfe_version_yaml = _render_keyed_values_yaml(
        "openfe_version", metadata.openfe_version, "version", "edges"
    )
    openmm_version_yaml = _render_keyed_values_yaml(
        "openmm_version", metadata.openmm_version, "version", "edges"
    )
    openff_toolkit_version_yaml = _render_keyed_values_yaml(
        "openff_toolkit_version",
        metadata.openff_toolkit_version,
        "version",
        "edges",
    )
    forcefield_yaml = _render_keyed_values_yaml(
        "forcefield", metadata.forcefield, "forcefield", "edges"
    )
    partial_charges_yaml = _render_keyed_values_yaml(
        "partial_charges",
        metadata.partial_charges,
        "partial_charges",
        "edges",
    )

    # Build network keys to systems mapping section
    network_keys_section = ""
    network_breakdown = defaultdict(list)
    for si in metadata.system_info_dict.values():
        network_breakdown[si.network_key].append(
            f"{si.benchmark_set}/{si.benchmark_system}"
        )

    if network_breakdown:
        network_keys_lines = []
        for network_key_item, systems in sorted(network_breakdown.items()):
            network_keys_lines.append(
                f"  - {network_key_item}: {', '.join(sorted(set(systems)))}"
            )
        network_keys_section = "## Alchemical Network Keys:\n" + "\n".join(
            network_keys_lines
        )

    return f"""# {title}

## Description

## Overview

{content_kind} benchmark results prepared from {source_description} JSON file(s) generated with {workflow_text}.

{content_summary}

## Software versions

{openfe_version_yaml}
{openmm_version_yaml}
{openff_toolkit_version_yaml}

{network_keys_section}

## Recommended descriptors

{forcefield_yaml}
{partial_charges_yaml}

{benchmark_system_yaml}

## Protocol settings

{protocol_settings_yaml}

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

    out_dir = output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Check for required results file
    results_path = out_dir / results_file
    if not results_path.exists():
        raise FileNotFoundError(
            f"Required file '{results_file}' not found in output directory: {out_dir}"
        )

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

    # Process all input files and collect metadata
    all_network_objs: list[dict[str, Any]] = []
    modes: set[str] = set()
    all_network_keys: list[str] = []
    all_metadata: list[AutoMetadata] = []
    for input_path in input_paths:
        resolved_path = input_path.resolve()
        network_obj, network_mode = _load_network(resolved_path)
        all_network_objs.append(network_obj)

        metadata = _extract_auto_metadata(network_obj, network_mode, str(resolved_path))
        modes.add(metadata.calculation_mode)
        all_network_keys.append(metadata.network_key)
        all_metadata.append(metadata)

    # Check consistency
    if len(modes) > 1:
        raise ValueError(
            f"Mixed modes detected across input files: {modes}. All files must be either ASFE or RBFE."
        )
    mode = modes.pop()

    # Merge metadata from all files
    merged_metadata = AutoMetadata()
    merged_metadata.calculation_mode = mode
    for metadata in all_metadata:
        merged_metadata.n_transformations += metadata.n_transformations
        if not merged_metadata.network_key:
            merged_metadata.network_key = metadata.network_key
        else:
            merged_metadata.network_key += f", {metadata.network_key}"
        merged_metadata.benchmark_sets_systems.extend(metadata.benchmark_sets_systems)

        # Use first non-empty value for scalar fields
        for key in [
            "openmm_version",
            "openfe_version",
            "openff_toolkit_version",
            "forcefield",
            "partial_charges",
            "small_molecule_force_field",
            "protocols",
            "protocol_settings_list",
            "mapper",
        ]:
            for value, keys in getattr(metadata, key):
                _add_value_with_keys(getattr(merged_metadata, key), value, keys)

        if any(
            [
                x in merged_metadata.system_info_dict
                for x in metadata.system_info_dict.keys()
            ]
        ):
            raise ValueError(
                f"System is already documented: {[x for x in metadata.system_info_dict.keys() if x in merged_metadata.system_info_dict]}"
            )

        merged_metadata.system_info_dict.update(metadata.system_info_dict)

    # Build content summary from combined data
    # Get list of source file names
    content_summary = _build_content_summary(
        merged_metadata,
        used_alchemiscale,
    )

    # Append additional summary text if provided
    if summary_suffix:
        content_summary = textwrap.fill(
            content_summary.rstrip() + " " + summary_suffix.strip(), width=100
        )

    sets_to_systems: dict[str, list[str]] = defaultdict(list)
    for system_group, system_name in merged_metadata.benchmark_sets_systems:
        sets_to_systems[system_group].append(system_name)

    # Generate a descriptive title
    submission_id = submission_id or _default_submission_id("_".join(all_network_keys))

    title = _generate_title(mode, merged_metadata.benchmark_sets_systems, submission_id)

    submission_yaml_filename = "submission.yaml"
    zenodo_description_filename = "zenodo_description.md"

    submission_yaml_path = out_dir / submission_yaml_filename
    zenodo_description_path = out_dir / zenodo_description_filename

    tags_list = [k.strip() for k in tags.split(",") if k.strip()]
    tags_final = _make_tags(
        mode=mode,
        forcefield=merged_metadata.forcefield,
        partial_charge_tag=merged_metadata.partial_charges,
        benchmark_data=merged_metadata.benchmark_sets_systems,
        user_keywords=tags_list,
    )

    submission_yaml_text = _make_submission_yaml(
        merged_metadata,
        submission_id=submission_id,
        title=title,
        summary=content_summary,
        tags=tags_final,
        authors=author or [],
        archive_doi="TODO add DOI",
        archive_provider="TODO add archive provider",
        license_name=license,
        results_file=results_file,
    )
    submission_yaml_path.write_text(submission_yaml_text)

    # For Zenodo description, list all input files
    archive_filenames = ", ".join(p.name for p in input_paths)

    zenodo_description_text = _make_zenodo_description(
        merged_metadata,
        network_mode=network_mode,
        title=title,
        archive_filename=archive_filenames,
        mode=mode,
        content_summary=content_summary,
        license_name=license,
        used_alchemiscale=used_alchemiscale,
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


if __name__ == "__main__":
    sys.exit(main())
