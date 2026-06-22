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


@click.command()
@click.option(
    "--archive",
    help="Path to alchemical archive (.json.bz2 file)",
    type=click.Path(
        exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path
    ),
    default=None,
)
@click.option(
    "--network",
    help="Path to alchemical network JSON file (required if not using --archive)",
    type=click.Path(
        exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path
    ),
    default=None,
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
def main(archive, network, results_dir, output_dir, system_group, system_name):
    """
    Gather transformation results and compute DDG/DG values.

    Can accept input from either:
    - An alchemical archive (--archive) with embedded results
    - A network file (--network) with separate results directories (--results_dir)
    """
    # Validate input
    if archive and network:
        raise click.UsageError("Specify either --archive or --network, not both")

    if archive:
        # Load from archive
        print(f"Loading archive: {archive.name}")
        alchemical_archive = _load_archive(archive)
        network_obj = alchemical_archive.network

        print("Extracting results from archive")
        raw_results = _extract_results_from_archive(alchemical_archive)
        results_source = "archive"
    else:
        # Load from network + results files
        if not network:
            raise click.UsageError("Must specify either --archive or --network")
        if not results_dir:
            raise click.UsageError("Must specify --results_dir when using --network")

        print(f"Loading network: {network.name}")
        network_obj = AlchemicalNetwork.from_json(network.as_posix())

        print(
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

    print(f"Found {len(transformations_to_run)} transformations in network")

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

    print(f"System: {final_system_group} / {final_system_name}")

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

    # Build FEMap and extract DDG/DG values
    fe_map = FEMap()
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
            complex_dg = np.mean(complex_data) * unit.kilocalories_per_mole
            if len(complex_data) == 1:
                complex_dg_uncertainty = np.nan * unit.kilocalories_per_mole
            else:
                complex_dg_uncertainty = (
                    np.std(complex_data) * unit.kilocalories_per_mole
                )

            solvent_dg = np.mean(solvent_data) * unit.kilocalories_per_mole
            if len(solvent_data) == 1:
                solvent_dg_uncertainty = np.nan * unit.kilocalories_per_mole
            else:
                solvent_dg_uncertainty = (
                    np.std(solvent_data) * unit.kilocalories_per_mole
                )

            # get the combined ddg and uncertainty
            entry_data["ddg"] = complex_dg - solvent_dg
            entry_data["ddg_uncertainty"] = np.sqrt(
                complex_dg_uncertainty**2 + solvent_dg_uncertainty**2
            )

            # add the raw values for debugging
            entry_data["dgs_complex"] = complex_data
            entry_data["dgs_solvent"] = solvent_data

            # extract overlap and mixing matrices
            for phase_results, label in zip(
                [complex_results, solvent_results], ["complex", "solvent"]
            ):
                mbar_overlap_elements = []
                replica_mixing_elements = []

                for phase_result in phase_results:
                    if results_source == "archive":
                        # Archive results have mbar_overlap and replica_mixing already extracted
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

            gathered_results["ddg"].append(entry_data)

            # also add the DDG data to the FEMAP
            fe_map.add_relative_calculation(
                labelA=lig_a_name,
                labelB=lig_b_name,
                value=entry_data["ddg"],
                uncertainty=entry_data["ddg_uncertainty"],
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

    print(f"Writing results to: {output_file}")
    print(
        f"Done! Found {len(gathered_results['ddg'])} DDG entries and {len(gathered_results['dg'])} DG entries"
    )


if __name__ == "__main__":
    main()
