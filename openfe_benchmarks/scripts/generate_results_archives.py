import json
import click
import pathlib
from collections import defaultdict
import logging

import numpy as np
import bz2

from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openfe.protocols.openmm_septop import SepTopProtocol
from gufe.tokenization import JSON_HANDLER
from gufe import ProteinComponent, SolventComponent
from gufe.archival import AlchemicalArchive
from openff.units import unit
from cinnabar import FEMap
from pontibus.protocols.relative import HybridTopProtocol
from pontibus.protocols.solvation import ASFEProtocol

logger = logging.getLogger(__name__)

MIN_ALLOWED_REPEATS = 3  # Less than this results in an error


def _load_archive(archive_path: pathlib.Path):
    """
    Load an AlchemicalArchive from a bz2-compressed JSON archive.

    Parameters
    ----------
    archive_path : pathlib.Path
        Path to the .json.bz2 archive file

    Returns
    -------
    AlchemicalArchive
        The deserialized alchemical archive
    """

    archive_path = pathlib.Path(archive_path)

    if str(archive_path).endswith(".bz2"):
        with bz2.open(archive_path, "rt") as f:
            json_content = f.read()
        archive = AlchemicalArchive.from_json(content=json_content)
    else:
        archive = AlchemicalArchive.from_json(str(archive_path))

    return archive


def _extract_hybrid_topology_rfe_data(transformation, results):
    """We assume that the inputs were prepared with the _example_plan_rbfe.py script and have the system name and group stored in the mapping."""
    if transformation.stateA.contains(ProteinComponent):
        phase = "complex"
    elif transformation.stateA.contains(SolventComponent):
        phase = "solvent"
    else:
        phase = "TODO"

    ligand_a_name = transformation.mapping.componentA.name
    ligand_b_name = transformation.mapping.componentB.name
    protocol_result = transformation.protocol.gather(results)
    individual_estimates = protocol_result.get_individual_estimates()
    overlap_matrices = protocol_result.get_overlap_matrices()
    replica_transition_statistics = protocol_result.get_replica_transition_statistics()
    equilibration_iterations = protocol_result.equilibration_iterations()
    production_iterations = protocol_result.production_iterations()

    return {
        "ligand_a": ligand_a_name,
        "ligand_b": ligand_b_name,
        "system_group": transformation.mapping.annotations.get("system_group", "TODO"),
        "system_name": transformation.mapping.annotations.get("system_name", "TODO"),
        "estimate": protocol_result.get_estimate(),
        "estimate_error": protocol_result.get_uncertainty(),
        "phase": phase,
        "individual_estimates": [e[0] for e in individual_estimates],
        "individual_mbar_errors": [e[1] for e in individual_estimates],
        "smallest_mbar_overlaps": [
            np.diagonal(om["matrix"], offset=1).min() for om in overlap_matrices
        ],
        "smallest_replica_mixing": [
            np.diagonal(rts["matrix"], offset=1).min()
            for rts in replica_transition_statistics
        ],
        "equilibration_iterations": equilibration_iterations,
        "production_iterations": production_iterations,
    }


def _extract_asfe_data(transformation, results):
    """We assume that the inputs were prepared with the _example_plan_asfe.py script and the system components are in the transformation name"""
    solute, solvent = transformation.name.split(",")
    protocol_result = transformation.protocol.gather(results)
    individual_estimates = protocol_result.get_individual_estimates()
    dgs_solvent = [e[0] for e in individual_estimates["solvent"]]
    solvent_mbar_errors = [e[1] for e in individual_estimates["solvent"]]
    dgs_vacuum = [e[0] for e in individual_estimates["vacuum"]]
    vacuum_mbar_errors = [e[1] for e in individual_estimates["vacuum"]]
    overlap_matrices = protocol_result.get_overlap_matrices()
    replica_mixing_matrices = protocol_result.get_replica_transition_statistics()
    equilibration_iterations = protocol_result.equilibration_iterations()
    production_iterations = protocol_result.production_iterations()

    return {
        "solute": solute,
        "solvent": solvent,
        "estimate": protocol_result.get_estimate(),
        "estimate_error": protocol_result.get_uncertainty(),
        "dgs_solvent": dgs_solvent,
        "solvent_mbar_errors": solvent_mbar_errors,
        "dgs_vacuum": dgs_vacuum,
        "vacuum_mbar_errors": vacuum_mbar_errors,
        "solvent_smallest_mbar_overlaps": [
            np.diagonal(om["matrix"], offset=1).min()
            for om in overlap_matrices["solvent"]
        ],
        "vacuum_smallest_mbar_overlaps": [
            np.diagonal(om["matrix"], offset=1).min()
            for om in overlap_matrices["vacuum"]
        ],
        "solvent_smallest_replica_mixing": [
            np.diagonal(rts["matrix"], offset=1).min()
            for rts in replica_mixing_matrices["solvent"]
        ],
        "vacuum_smallest_replica_mixing": [
            np.diagonal(rts["matrix"], offset=1).min()
            for rts in replica_mixing_matrices["vacuum"]
        ],
        "solvent_equilibration_iterations": equilibration_iterations["solvent"],
        "vacuum_equilibration_iterations": equilibration_iterations["vacuum"],
        "solvent_production_iterations": production_iterations["solvent"],
        "vacuum_production_iterations": production_iterations["vacuum"],
    }


def _extract_septop_rbfe_data(transformation, results):
    """We assume that the inputs were prepared with the _example_plan_septop.py script and the system components are in the transformation name."""
    protocol_result = transformation.protocol.gather(results)
    individual_estimates = protocol_result.get_individual_estimates()
    ligand_a_name = transformation.mapping.componentA.name
    ligand_b_name = transformation.mapping.componentB.name
    dgs_solvent = [e[0] for e in individual_estimates["solvent"]]
    solvent_mbar_errors = [e[1] for e in individual_estimates["solvent"]]
    dgs_complex = [e[0] for e in individual_estimates["complex"]]
    complex_mbar_errors = [e[1] for e in individual_estimates["complex"]]
    overlap_matrices = protocol_result.get_overlap_matrices()
    replica_mixing_matrices = protocol_result.get_replica_transition_statistics()
    equilibration_iterations = protocol_result.equilibration_iterations()
    production_iterations = protocol_result.production_iterations()

    return {
        "ligand_a": ligand_a_name,
        "ligand_b": ligand_b_name,
        "system_group": transformation.mapping.annotations.get("system_group", "TODO"),
        "system_name": transformation.mapping.annotations.get("system_name", "TODO"),
        "ddg": protocol_result.get_estimate(),
        "ddg_uncertainty": protocol_result.get_uncertainty(),
        "dgs_solvent": dgs_solvent,
        "solvent_mbar_errors": solvent_mbar_errors,
        "dgs_complex": dgs_complex,
        "complex_mbar_errors": complex_mbar_errors,
        "solvent_smallest_mbar_overlaps": [
            np.diagonal(om["matrix"], offset=1).min()
            for om in overlap_matrices["solvent"]
        ],
        "complex_smallest_mbar_overlaps": [
            np.diagonal(om["matrix"], offset=1).min()
            for om in overlap_matrices["complex"]
        ],
        "solvent_smallest_replica_mixing": [
            np.diagonal(rts["matrix"], offset=1).min()
            for rts in replica_mixing_matrices["solvent"]
        ],
        "complex_smallest_replica_mixing": [
            np.diagonal(rts["matrix"], offset=1).min()
            for rts in replica_mixing_matrices["complex"]
        ],
        # for analytical corrections there is no uncertainty
        "standard_state_correction_complex_A": [
            e[0] for e in individual_estimates["standard_state_correction_complex_A"]
        ],
        "standard_state_correction_complex_B": [
            e[0] for e in individual_estimates["standard_state_correction_complex_B"]
        ],
        "standard_state_correction_solvent": [
            e[0] for e in individual_estimates["standard_state_correction_solvent"]
        ],
        "complex_production_iterations": production_iterations["complex"],
        "solvent_production_iterations": production_iterations["solvent"],
        "complex_equlibration_iterations": equilibration_iterations["complex"],
        "solvent_equlibration_iterations": equilibration_iterations["solvent"],
    }


def _combine_hybrid_topology_results(results: list[dict]) -> list[dict]:
    """For the given list of hybrid topology results combine compatible legs of the calculations to form a single estimate in line with other protocols.

    Notes
    -----
    - We currently only support combining solvent and complex legs into relative binding free energies
    """
    # group the results by a key made from the ligand names
    results_by_key = defaultdict(list)

    for result in results:
        key = (result["ligand_a"], result["ligand_b"])
        results_by_key[key].append(result)

    processed_results = []
    for paired_results in results_by_key.values():
        # construct a single output for this transformation
        results_by_phase = dict((r["phase"], r) for r in paired_results)
        if (phases := set(results_by_phase.keys())) != {"solvent", "complex"}:
            raise RuntimeError(
                f"Only relative binding free energy calculations are currently supported with the hybrid topology found phases: {phases}"
            )

        combined_data = {
            "ligand_a": results_by_phase["solvent"]["ligand_a"],
            "ligand_b": results_by_phase["solvent"]["ligand_b"],
            "system_group": results_by_phase["solvent"]["system_group"],
            "system_name": results_by_phase["solvent"]["system_name"],
            "ddg": results_by_phase["complex"]["estimate"]
            - results_by_phase["solvent"]["estimate"],
            "ddg_uncertainty": np.sqrt(
                results_by_phase["complex"]["estimate_error"] ** 2
                + results_by_phase["solvent"]["estimate_error"] ** 2
            ),
            "dgs_complex": results_by_phase["complex"]["individual_estimates"],
            "dgs_solvent": results_by_phase["solvent"]["individual_estimates"],
            "complex_mbar_errors": results_by_phase["complex"][
                "individual_mbar_errors"
            ],
            "solvent_mbar_errors": results_by_phase["solvent"][
                "individual_mbar_errors"
            ],
            "complex_smallest_mbar_overlaps": results_by_phase["complex"][
                "smallest_mbar_overlaps"
            ],
            "solvent_smallest_mbar_overlaps": results_by_phase["solvent"][
                "smallest_mbar_overlaps"
            ],
            "complex_smallest_replica_mixing": results_by_phase["complex"][
                "smallest_replica_mixing"
            ],
            "solvent_smallest_replica_mixing": results_by_phase["solvent"][
                "smallest_replica_mixing"
            ],
            "complex_equilibration_iterations": results_by_phase["complex"][
                "equilibration_iterations"
            ],
            "solvent_equilibration_iterations": results_by_phase["solvent"][
                "equilibration_iterations"
            ],
            "complex_production_iterations": results_by_phase["complex"][
                "production_iterations"
            ],
            "solvent_production_iterations": results_by_phase["solvent"][
                "production_iterations"
            ],
        }
        processed_results.append(combined_data)
    return processed_results


def _extract_results_from_archive(alchemical_archive):
    """Extract results from an AlchemicalArchive.

    Returns a dictionary keyed by the protocol class with lists of result dictionaries
    for each transformation.
    """
    extraction_functions = {
        RelativeHybridTopologyProtocol: _extract_hybrid_topology_rfe_data,
        HybridTopProtocol: _extract_hybrid_topology_rfe_data,
        ASFEProtocol: _extract_asfe_data,
        SepTopProtocol: _extract_septop_rbfe_data,
        # TODO add support for openfe ASFE
        # TODO add support for openfe ABFE
    }

    raw_results = defaultdict(list)

    for transformation, dag_results_list in alchemical_archive.transformation_results:
        if len(dag_results_list) < MIN_ALLOWED_REPEATS:
            raise ValueError(
                f"Transformation {transformation.key} is does not meet minimum number of repeats requirement. Must at least {MIN_ALLOWED_REPEATS}."
            )
        protocol_cls = transformation.protocol.__class__
        extract_func = extraction_functions.get(protocol_cls)
        result_data = extract_func(transformation, dag_results_list)
        raw_results[protocol_cls].append(result_data)

    if RelativeHybridTopologyProtocol in raw_results:
        raw_results[RelativeHybridTopologyProtocol] = _combine_hybrid_topology_results(
            raw_results[RelativeHybridTopologyProtocol]
        )
    if HybridTopProtocol in raw_results:
        raw_results[HybridTopProtocol] = _combine_hybrid_topology_results(
            raw_results[HybridTopProtocol]
        )

    return raw_results


def _parse_systems_input(systems):
    """Normalize systems iterable into a list of tuples.

    Expected shape is an iterable of (system_group, system_name, archive_path).
    """
    parsed = []
    for item in systems:
        if not isinstance(item, (list, tuple)) or len(item) != 3:
            raise ValueError(
                "Each item in systems must be a tuple/list of "
                "(system_group, system_name, archive_path)"
            )
        system_group, system_name, archive_path = item
        if system_group is None or system_name is None:
            raise ValueError(
                "Each systems entry must include a non-empty system_group and system_name"
            )
        parsed.append((system_group, system_name, pathlib.Path(archive_path)))
    if not parsed:
        raise ValueError("At least one system must be provided in systems")
    return parsed


def run_generate_results(
    archive=None,
    output_dir=None,
    system_group=None,
    system_name=None,
    systems=None,
):
    """Generate computational results from one or more alchemical archives.

    Parameters
    ----------
    archive : Path, optional
        Path to the alchemical archive file (.json.bz2)
    output_dir : Path, optional
        Directory to write the results JSON to
    system_group : str, optional
        Benchmark set name (e.g., 'jacs_set', 'solvation_set')
    system_name : str, optional
        System name (e.g., 'tyk2', 'hsp90')
    systems : iterable of tuple, optional
        Iterable of (system_group, system_name, archive_path) entries to process
        multiple archives together. Replaces archive, system_group, system_name
    """

    # Convert string paths to Path objects if needed
    if output_dir is not None and isinstance(output_dir, str):
        output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    archive_results_by_protocol = defaultdict(list)
    if systems is not None:
        if archive is not None:
            raise ValueError(
                "Cannot specify both archive and systems. Use one input style only."
            )
        if system_group is not None or system_name is not None:
            raise ValueError(
                "When providing systems, do not also pass system_group or system_name. "
                "Use per-system overrides in systems instead."
            )

        system_entries = _parse_systems_input(systems)
        for system_group_entry, system_name_entry, archive_path in system_entries:
            logger.info(f"Loading archive:{archive_path}")
            alchemical_archive = _load_archive(archive_path)
            logger.info("Extracting results from archive")
            archive_results = _extract_results_from_archive(alchemical_archive)
            for protocol_cls, results in archive_results.items():
                for result in results:
                    result["system_group"] = system_group_entry
                    result["system_name"] = system_name_entry
                archive_results_by_protocol[protocol_cls].extend(results)
    else:
        if archive is None:
            raise ValueError("archive must be provided when systems is not used")

        # Load from archive
        logger.info(f"Loading archive:{archive}")
        alchemical_archive = _load_archive(archive)

        logger.info("Extracting results from archive")
        archive_results_by_protocol = _extract_results_from_archive(alchemical_archive)

    # raise an error if we find specific mixes of protocols
    if ASFEProtocol in archive_results_by_protocol and (
        HybridTopProtocol in archive_results_by_protocol
        or RelativeHybridTopologyProtocol in archive_results_by_protocol
    ):
        raise RuntimeError(
            f"Found a mix of ASFE and Hybrid Topology protocols in the archive. This is not supported. Found protocols: {list(archive_results_by_protocol.keys())}"
        )

    # Initialize output structure
    gathered_results = {"dg": [], "ddg": []}

    final_system_group = None
    final_system_name = None
    if systems is None:
        # Get system group and system name from the first result we can find
        for p_rs in archive_results_by_protocol.values():
            data = p_rs[0]
            extracted_system_group = data.get("system_group", "TODO")
            extracted_system_name = data.get("system_name", "TODO")
            break

        # Use CLI options if provided, otherwise use extracted values
        final_system_group = system_group if system_group else extracted_system_group
        final_system_name = system_name if system_name else extracted_system_name
        logger.info(f"System: {final_system_group} / {final_system_name}")
    else:
        logger.info("System groups and names are set per archive from systems input")

    # workout what we are going to do based on the found protocols
    if ASFEProtocol in archive_results_by_protocol:
        # if its just absolute values add them to the DGs and be done
        for result in archive_results_by_protocol[ASFEProtocol]:
            if final_system_group is not None:
                result["system_group"] = final_system_group
            if final_system_name is not None:
                result["system_name"] = final_system_name
            gathered_results["dg"].append(result)

    else:
        # we have some relative calculations add them to the ddg results and then build the DGs from the ddgs
        for results in archive_results_by_protocol.values():
            for result in results:
                if final_system_group is not None:
                    result["system_group"] = final_system_group
                if final_system_name is not None:
                    result["system_name"] = final_system_name
                gathered_results["ddg"].append(result)

        # now build the DGs from the ddgs grouped by system
        ddgs_by_system = defaultdict(list)
        for result in gathered_results["ddg"]:
            system_key = (
                result.get("system_group", "TODO"),
                result.get("system_name", "TODO"),
            )
            ddgs_by_system[system_key].append(result)

        for (system_group_key, system_name_key), system_ddgs in ddgs_by_system.items():
            fe_map = FEMap()
            for result in system_ddgs:
                fe_map.add_relative_calculation(
                    labelA=result["ligand_a"],
                    labelB=result["ligand_b"],
                    value=result["ddg"],
                    uncertainty=result["ddg_uncertainty"],
                )

            if not fe_map.check_weakly_connected():
                logger.warning(
                    "Could not generate DGs for system %s / %s because the relative map is not weakly connected.",
                    system_group_key,
                    system_name_key,
                )
                continue

            try:
                fe_map.generate_absolute_values()
                abs_df = fe_map.get_absolute_dataframe()
                for _, row in abs_df.iterrows():
                    entry_data = {
                        "ligand": row["label"],
                        "dg": row["DG (kcal/mol)"] * unit.kilocalories_per_mole,
                        "dg_uncertainty": row["uncertainty (kcal/mol)"]
                        * unit.kilocalories_per_mole,
                        "system_group": system_group_key,
                        "system_name": system_name_key,
                        "source": row["source"],
                    }
                    gathered_results["dg"].append(entry_data)
            except Exception as e:
                logger.warning(
                    "Could not generate absolute values (DG) for the alchemical map of system %s / %s: %s. "
                    "This may occur when uncertainties are NaN (single replicate data).",
                    system_group_key,
                    system_name_key,
                    e,
                )

    # write out the data to a json file
    output_file = output_dir / "computational_results.json"
    with open(output_file, "w") as w:
        json.dump(gathered_results, w, cls=JSON_HANDLER.encoder, indent=4)

    logger.info(f"Writing results to: {output_file}")
    logger.info(
        f"Done! Found {len(gathered_results['ddg'])} DDG entries and {len(gathered_results['dg'])} DG entries"
    )


@click.command()
@click.option(
    "--archive",
    help="Path to alchemical archive (.json.bz2 file; can be specified multiple times)",
    type=click.Path(
        exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path
    ),
    required=True,
)
@click.option(
    "--output-dir",
    help="Directory to write the results JSON to",
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    required=True,
)
@click.option(
    "--system-group",
    help="Benchmark set name (e.g., 'jacs_set', 'solvation_set'); overrides value from network annotations",
    type=str,
    default=None,
)
@click.option(
    "--system-name",
    help="System name (e.g., 'tyk2', 'hsp90'); overrides value from network annotations",
    type=str,
    default=None,
)
def main(archive, output_dir, system_group, system_name):
    """CLI wrapper for run_generate_results.

    Gather transformation results and compute DDG/DG values.
    """
    run_generate_results(
        archive=archive,
        output_dir=output_dir,
        system_group=system_group,
        system_name=system_name,
    )


if __name__ == "__main__":
    main()
