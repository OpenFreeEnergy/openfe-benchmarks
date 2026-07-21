import json
import click
from gufe.tokenization import JSON_HANDLER
from gufe import AlchemicalNetwork, ProteinComponent
from openff.units import unit
import pathlib
from cinnabar import FEMap
from collections import defaultdict
import numpy as np
import bz2
import logging
import tempfile

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
    from gufe.archival import AlchemicalArchive

    with bz2.open(archive_path, "rt") as f:
        json_content = f.read()

    # Write to temporary file for AlchemicalArchive.from_json()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tmp:
        tmp.write(json_content)
        tmp_path = tmp.name

    try:
        archive = AlchemicalArchive.from_json(tmp_path)
        return archive
    finally:
        pathlib.Path(tmp_path).unlink()


def _extract_results_from_archive(alchemical_archive):
    """
    Extract results from an AlchemicalArchive.

    Returns a dictionary keyed by (ligand_a, ligand_b) with lists of (phase, result_dict)
    """
    raw_results = defaultdict(list)

    for transformation, dag_results_list in alchemical_archive.transformation_results:
        ligand_a_name = transformation.mapping.componentA.name
        ligand_b_name = transformation.mapping.componentB.name

        # Determine phase from transformation state
        if transformation.stateA.contains(ProteinComponent):
            phase = "complex"
        else:
            phase = "solvent"

        key = (ligand_a_name, ligand_b_name)

        # Extract estimate from the DAG result
        if dag_results_list:
            dag_result = dag_results_list[0]
            if dag_result.terminal_protocol_unit_results:
                term_result = dag_result.terminal_protocol_unit_results[0]
                outputs = term_result.outputs

                # Extract unit_estimate and error
                estimate = None
                estimate_error = None
                mbar_overlap_matrix = None
                replica_mixing_matrix = None

                if "unit_estimate" in outputs:
                    estimate_qty = outputs["unit_estimate"]
                    estimate = estimate_qty.magnitude

                if "unit_estimate_error" in outputs:
                    error_qty = outputs["unit_estimate_error"]
                    estimate_error = error_qty.magnitude

                if "unit_mbar_overlap" in outputs:
                    mbar_dict = outputs["unit_mbar_overlap"]
                    if isinstance(mbar_dict, dict) and "matrix" in mbar_dict:
                        mbar_overlap_matrix = mbar_dict["matrix"]

                if "replica_exchange_statistics" in outputs:
                    replica_dict = outputs["replica_exchange_statistics"]
                    if isinstance(replica_dict, dict) and "matrix" in replica_dict:
                        replica_mixing_matrix = replica_dict["matrix"]

                result_dict = {
                    "phase": phase,
                    "estimate": estimate,
                    "estimate_error": estimate_error,
                    "mbar_overlap": mbar_overlap_matrix,
                    "replica_mixing": replica_mixing_matrix,
                    "transformation": transformation,
                }
                raw_results[key].append((phase, result_dict))

    return raw_results


def _extract_results_from_files(results_dirs):
    """
    Extract results from result files on disk.

    Returns a dictionary keyed by (ligand_a, ligand_b) with lists of (phase, result_dict)
    """
    import zstandard as zstd
    from openfecli.commands.gather import _get_names, _get_type

    def _get_simulation_key(result: dict) -> tuple[tuple[str, str], str]:
        lig_a_name, lig_b_name = _get_names(result)
        phase = _get_type(result)
        return (lig_a_name, lig_b_name), phase

    raw_results = defaultdict(list)

    for result_dir in results_dirs:
        # search for the results json files
        for result_file in result_dir.glob("*.json.zst"):
            with open(result_file, "rb") as f:
                dctx = zstd.ZstdDecompressor()
                with dctx.stream_reader(f) as reader:
                    result = json.load(reader, cls=JSON_HANDLER.decoder)
                    # make a key for this result
                    key, phase = _get_simulation_key(result)
                    raw_results[key].append((phase, result))

    return raw_results


def run_generate_results(
    archive=None,
    network=None,
    network_key=None,
    results_dir=None,
    output_dir=None,
    system_group=None,
    system_name=None,
):
    """Generate computational results from alchemical archives or networks.

    Parameters
    ----------
    archive : tuple of Path, optional
        Paths to alchemical archive files (.json.bz2)
    network : tuple of Path, optional
        Paths to alchemical network JSON files
    network_key : tuple of str, optional
        Scope keys of networks from alchemiscale
    results_dir : tuple of Path, optional
        Directories containing transformation results
    output_dir : Path, optional
        Directory to write the results JSON to
    system_group : str, optional
        Benchmark set name (e.g., 'jacs_set', 'solvation_set')
    system_name : str, optional
        System name (e.g., 'tyk2', 'hsp90')

    Notes
    -----
    Can accept input from either:
    - One or more alchemical archives (archive) with embedded results
    - One or more network files (network) with separate results directories (results_dir)
    - One or more network keys (network_key) from alchemiscale with separate results directories (results_dir)
    """
    # Convert tuples to lists if needed, handle None values
    archive = tuple(archive) if archive else ()
    network = tuple(network) if network else ()
    network_key = tuple(network_key) if network_key else ()
    results_dir = tuple(results_dir) if results_dir else ()

    # Convert string paths to Path objects if needed
    if output_dir is not None and isinstance(output_dir, str):
        output_dir = pathlib.Path(output_dir)

    # Validate input combinations
    input_sources = sum([bool(archive), bool(network), bool(network_key)])
    if input_sources > 1:
        raise ValueError(
            "Specify only one of: archive, network, or network_key (not combinations)"
        )
    if input_sources == 0:
        raise ValueError("Must specify one of: archive, network, or network_key")

    if archive:
        # Load from archives
        nl = "\n"
        logger.info(
            f"Loading {len(archive)} archive(s):\n{nl.join([x.name for x in archive])}"
        )
        alchemical_archives = [_load_archive(arch) for arch in archive]
        network_obj = alchemical_archives[0].network

        # Verify all archives have the same network
        for i, arch in enumerate(alchemical_archives[1:], start=1):
            if arch.network != network_obj:
                raise ValueError(
                    f"Archive {i} has a different network than the first archive. "
                    "All archives must have the same network."
                )

        logger.info("Extracting results from archives")
        raw_results = defaultdict(list)
        for arch in alchemical_archives:
            arch_results = _extract_results_from_archive(arch)
            # Merge results from this archive into the combined dictionary
            for key, results_list in arch_results.items():
                raw_results[key].extend(results_list)

        results_source = "archive"
    elif network or network_key:
        # Load from network files or network keys + results files
        if not results_dir:
            raise click.UsageError(
                "Must specify --results_dir when using --network or --network-key"
            )

        networks_to_load = []

        if network:
            logger.info(f"Loading {len(network)} network file(s)")
            for net_file in network:
                logger.info(f"  - {net_file.name}")
                networks_to_load.append(
                    AlchemicalNetwork.from_json(net_file.as_posix())
                )

        if network_key:
            logger.info(f"Loading {len(network_key)} network(s) from alchemiscale")
            from alchemiscale import AlchemiscaleClient, ScopedKey

            client = AlchemiscaleClient(api_url="https://api.alchemiscale.org")
            for key in network_key:
                logger.info(f"  - {key}")
                scoped_key = ScopedKey.from_str(key)
                networks_to_load.append(client.get_network(scoped_key))

        # Verify all networks are identical
        network_obj = networks_to_load[0]
        for i, net in enumerate(networks_to_load[1:], start=1):
            if net != network_obj:
                raise ValueError(
                    f"Network {i} differs from the first network. "
                    "All networks must be identical."
                )

        logger.info(
            f"Loading results from {len(results_dir)} director{'ies' if len(results_dir) != 1 else 'y'}"
        )
        raw_results = _extract_results_from_files(results_dir)
        results_source = "files"

    # Build a set of expected transformations from the network
    transformations_to_run = set()
    for transformation in network_obj.edges:
        ligand_a_name = transformation.mapping.componentA.name
        ligand_b_name = transformation.mapping.componentB.name
        # get the phase
        if transformation.stateA.contains(ProteinComponent):
            phase = "complex"
        else:
            phase = "solvent"
        transformations_to_run.add((ligand_a_name, ligand_b_name, phase))

    logger.info(f"Found {len(transformations_to_run)} transformations in network")

    # Initialize output structure
    gathered_results = {"dg": [], "ddg": []}

    # Get system group and system name from the first edge
    if len(list(network_obj.edges)) > 0:
        transformation = list(network_obj.edges)[0]
        mapping_annotations = transformation.mapping.annotations
        extracted_system_group = mapping_annotations.get("system_group", "unknown")
        extracted_system_name = mapping_annotations.get("system_name", "unknown")
    else:
        extracted_system_group = "unknown"
        extracted_system_name = "unknown"

    # Use CLI options if provided, otherwise use extracted values
    final_system_group = system_group if system_group else extracted_system_group
    final_system_name = system_name if system_name else extracted_system_name

    logger.info(f"System: {final_system_group} / {final_system_name}")

    # Check that all simulations in the alchemical network have an associated result
    found_results = set()
    for key, results in raw_results.items():
        lig_a_name, lig_b_name = key
        for phase, result in results:
            result_key = (lig_a_name, lig_b_name, phase)
            if result_key not in transformations_to_run:
                raise ValueError(
                    f"Found results for transformation {result_key} which is not in the alchemical network"
                )
            found_results.add((lig_a_name, lig_b_name, phase))

    # Check for missing transformations
    missing_transformations = transformations_to_run - found_results
    if missing_transformations:
        raise ValueError(
            f"Missing results for transformations: {missing_transformations}"
        )

    # First pass: Build gathered_results with both ddg_uncertainty and mbar_err
    # We'll decide globally which uncertainty type to use after processing all edges
    for key, results in raw_results.items():
        lig_a_name, lig_b_name = key
        entry_data = {
            "ligand_a": lig_a_name,
            "ligand_b": lig_b_name,
            "system_group": final_system_group,
            "system_name": final_system_name,
        }

        # group the results by phase
        complex_results = [result for phase, result in results if phase == "complex"]
        solvent_results = [result for phase, result in results if phase == "solvent"]

        assert len(complex_results) == len(solvent_results), (
            f"Found different number of complex and solvent results for {key}"
        )

        # Extract estimates based on source type
        if results_source == "archive":
            complex_data = [
                result["estimate"]
                for result in complex_results
                if result["estimate"] is not None
            ]
            solvent_data = [
                result["estimate"]
                for result in solvent_results
                if result["estimate"] is not None
            ]
        else:  # files
            complex_data = [
                result["estimate"].m_as(unit.kilocalories_per_mole)
                for result in complex_results
            ]
            solvent_data = [
                result["estimate"].m_as(unit.kilocalories_per_mole)
                for result in solvent_results
            ]

        if complex_data and solvent_data:
            n_repeats = (len(complex_data), len(solvent_data))
            if n_repeats[0] < MIN_ALLOWED_REPEATS and n_repeats[0] > 1:
                raise ValueError(
                    f"Complex leg {key} is does not meet minimum number of repeats requirement. Must be 1 or at least {MIN_ALLOWED_REPEATS}."
                )
            if n_repeats[1] < MIN_ALLOWED_REPEATS and n_repeats[1] > 1:
                raise ValueError(
                    f"Solvent leg {key} is does not meet minimum number of repeats requirement. Must be 1 or at least {MIN_ALLOWED_REPEATS}."
                )

            complex_data *= unit.kilocalories_per_mole
            solvent_data *= unit.kilocalories_per_mole

            complex_dg = np.mean(complex_data)
            if n_repeats[0] == 1:
                complex_dg_uncertainty = np.nan * unit.kilocalories_per_mole
            else:
                complex_dg_uncertainty = np.std(complex_data)

            solvent_dg = np.mean(solvent_data)
            if n_repeats[1] == 1:
                solvent_dg_uncertainty = np.nan * unit.kilocalories_per_mole
            else:
                solvent_dg_uncertainty = np.std(solvent_data)

            # get the combined ddg and uncertainty
            entry_data["ddg"] = complex_dg - solvent_dg
            entry_data["ddg_uncertainty"] = np.sqrt(
                complex_dg_uncertainty**2 + solvent_dg_uncertainty**2
            )

            # add the raw values for debugging
            entry_data["dgs_complex"] = list(complex_data)
            entry_data["dgs_solvent"] = list(solvent_data)

            # extract overlap and mixing matrices, and calculate mbar uncertainty
            # mbar_errors will store the MBAR estimate errors for each repeat
            complex_mbar_errors = []
            solvent_mbar_errors = []

            for phase_results, label, mbar_errors in zip(
                [complex_results, solvent_results],
                ["complex", "solvent"],
                [complex_mbar_errors, solvent_mbar_errors],
            ):
                mbar_overlap_elements = []
                replica_mixing_elements = []

                for phase_result in phase_results:
                    if results_source == "archive":
                        # Archive results have estimate_error, mbar_overlap and replica_mixing already extracted
                        if phase_result.get("estimate_error") is not None:
                            mbar_errors.append(phase_result["estimate_error"])

                        if phase_result["mbar_overlap"] is not None:
                            mbar_overlap_elements.append(
                                np.diagonal(
                                    phase_result["mbar_overlap"], offset=1
                                ).min()
                            )
                        if phase_result["replica_mixing"] is not None:
                            replica_mixing_elements.append(
                                np.diagonal(
                                    phase_result["replica_mixing"], offset=1
                                ).min()
                            )
                    else:
                        # File results have nested structure
                        result_key = list(
                            phase_result["protocol_result"]["data"].keys()
                        )[0]

                        # Extract MBAR estimate error
                        estimate_error = phase_result["protocol_result"]["data"][
                            result_key
                        ][0]["outputs"].get("unit_estimate_error")
                        if estimate_error is not None:
                            mbar_errors.append(
                                estimate_error.m_as(unit.kilocalories_per_mole)
                            )

                        overlap_matrix = phase_result["protocol_result"]["data"][
                            result_key
                        ][0]["outputs"]["unit_mbar_overlap"]["matrix"]
                        mbar_overlap_elements.append(
                            np.diagonal(overlap_matrix, offset=1).min()
                        )

                        mixing_matrix = phase_result["protocol_result"]["data"][
                            result_key
                        ][0]["outputs"]["replica_exchange_statistics"]["matrix"]
                        replica_mixing_elements.append(
                            np.diagonal(mixing_matrix, offset=1).min()
                        )

                if mbar_overlap_elements:
                    entry_data[f"{label}_smallest_mbar_overlaps"] = (
                        mbar_overlap_elements
                    )

                if replica_mixing_elements:
                    entry_data[f"{label}_smallest_replica_mixing"] = (
                        replica_mixing_elements
                    )

                if mbar_errors:
                    mbar_errors_with_units = [
                        err * unit.kilocalories_per_mole for err in mbar_errors
                    ]
                    entry_data[f"{label}_mbar_errors"] = mbar_errors_with_units

            if np.isnan(entry_data["ddg_uncertainty"]):
                if (
                    "complex_mbar_errors" in entry_data
                    and "solvent_mbar_errors" in entry_data
                ):
                    if all([x == 1 for x in n_repeats]):
                        entry_data["mbar_error"] = np.sqrt(
                            entry_data["complex_mbar_errors"][0] ** 2
                            + entry_data["solvent_mbar_errors"][0] ** 2
                        )
                    else:
                        raise ValueError(
                            f"Uncertainty between complex/solvent {n_repeats} repeats is NaN "
                        )
                else:
                    raise ValueError(
                        f"No uncertainty can be derived from complex/solvent {n_repeats} repeats, and MBAR errors are unavailable."
                    )

            gathered_results["ddg"].append(entry_data)

    # Determine globally which uncertainty type to use for ALL edges
    # Check if all edges have valid ddg_uncertainty (from repeats)
    all_have_repeat_uncertainty = all(
        not np.isnan(result["ddg_uncertainty"].magnitude)
        for result in gathered_results["ddg"]
    )

    edges_without_repeat_uncertainty = [
        (result["ligand_a"], result["ligand_b"])
        for result in gathered_results["ddg"]
        if np.isnan(result["ddg_uncertainty"].magnitude)
    ]

    edges_without_min_repeats = [
        (result["ligand_a"], result["ligand_b"])
        for result in gathered_results["ddg"]
        if len(result["dgs_complex"]) < MIN_ALLOWED_REPEATS
    ]

    # Decide which uncertainty to use for the entire network
    if all_have_repeat_uncertainty:
        if edges_without_min_repeats:
            raise ValueError(
                f"Some edges have fewer than {MIN_ALLOWED_REPEATS}: {edges_without_min_repeats}"
            )
        else:
            use_mbar_err = False
            logger.info(
                "All edges have repeat-based uncertainties. Using ddg_uncertainty for all edges."
            )
    else:
        use_mbar_err = True
        logger.warning(
            f"Some edges lack repeat-based uncertainties. Using mbar_err for ALL edges. "
            f"Edges without repeat data: {edges_without_repeat_uncertainty}"
        )

    # Second pass: Build FEMap with consistent uncertainty type
    fe_map = FEMap()
    for result in gathered_results["ddg"]:
        if use_mbar_err:
            uncertainty = result["mbar_error"]
        else:
            uncertainty = result["ddg_uncertainty"]

        fe_map.add_relative_calculation(
            labelA=result["ligand_a"],
            labelB=result["ligand_b"],
            value=result["ddg"],
            uncertainty=uncertainty,
        )

    # check if the network is connected and we can calculate the DGs
    if fe_map.check_weakly_connected():
        # generate the absolute values for the map centered around zero
        # these should be shifted when comparing with experiment
        try:
            fe_map.generate_absolute_values()
            abs_df = fe_map.get_absolute_dataframe()
            for _, row in abs_df.iterrows():
                entry_data = {
                    "ligand": row["label"],
                    "dg": row["DG (kcal/mol)"] * unit.kilocalories_per_mole,
                    "dg_uncertainty": row["uncertainty (kcal/mol)"]
                    * unit.kilocalories_per_mole,
                    "system_group": system_group,
                    "system_name": system_name,
                    "source": row["source"],
                }
                gathered_results["dg"].append(entry_data)
        except Exception as e:
            logger.warning(
                f"Could not generate absolute values (DG) for the alchemical map: {e}. "
                "This may occur when uncertainties are NaN (single replicate data)."
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
    multiple=True,
)
@click.option(
    "--network",
    help="Path to alchemical network JSON file (can be specified multiple times or used with --archive)",
    type=click.Path(
        exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path
    ),
    multiple=True,
)
@click.option(
    "--network-key",
    help="Scope key of network from alchemiscale (can be specified multiple times)",
    type=str,
    multiple=True,
)
@click.option(
    "--results_dir",
    help="Directory containing transformation results (can be specified multiple times)",
    multiple=True,
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path
    ),
)
@click.option(
    "--output_dir",
    help="Directory to write the results JSON to",
    type=click.Path(
        exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path
    ),
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
def main(
    archive, network, network_key, results_dir, output_dir, system_group, system_name
):
    """CLI wrapper for run_generate_results.

    Gather transformation results and compute DDG/DG values.

    Can accept input from either:
    - One or more alchemical archives (--archive) with embedded results (can be specified multiple times)
    - One or more network files (--network) with separate results directories (--results_dir)
    - One or more network keys (--network-key) from alchemiscale with separate results directories (--results_dir)
    """
    run_generate_results(
        archive=archive,
        network=network,
        network_key=network_key,
        results_dir=results_dir,
        output_dir=output_dir,
        system_group=system_group,
        system_name=system_name,
    )


if __name__ == "__main__":
    main()
