from collections import defaultdict
import json

import numpy as np
from gufe.tokenization import JSON_HANDLER
from openff.units import unit
from cinnabar import FEMap

from openfe_benchmarks.data._benchmark_systems import get_benchmark_data_system


def _asfe_result_key(result: dict) -> str:
    """Build a consistent ASFE label and experimental lookup key."""
    ligand = result.get("ligand", result.get("solute"))
    if ligand is None:
        raise ValueError("ASFE result entry must contain either 'ligand' or 'solute'")

    solvent = result.get("solvent")
    return f"{ligand},{solvent}" if solvent else ligand


def build_femap_from_relative_results(
    results: list[dict],
    source: str = "",
) -> dict[tuple[str, str], FEMap]:
    """
    Build FEMaps for each of the unique combinations of system_group and system_name in the DDG results and add experimental data
    for each of the ligands present in the DDG results.

    Parameters
    ----------
    results: list[dict]
        A list of relative binding free energy estimates which should include at least the following entries:
         - ligand_a: str
         - ligand_b: str
         - system_group: str
         - system_name: str
         - ddg: Quantity
         - ddg_uncertainty: Quantity
    source: str, default=''
        Source of computed data, used in ``cinnabar`` for results comparisons between different BenchmarkResults.

    Returns
    -------
    dict[tuple[str, str], FEMap]
        A dictionary mapping each unique combination of system_group and system_name to an FEMap with calculated and experimental reference data.
    """
    # get the unique combinations of system_group and system_name
    results_by_system_key = defaultdict(list)
    for result in results:
        key = (result["system_group"], result["system_name"])
        results_by_system_key[key].append(result)

    femaps_by_system_key = {}
    for system_key, system_results in results_by_system_key.items():
        system_group, system_name = system_key
        benchmark_data = get_benchmark_data_system(system_group, system_name)
        if benchmark_data.reference_data is None:
            raise ValueError(
                f"No reference data available for benchmark system {system_group}/{system_name}"
            )
        unique_ligands = set()

        # Check if all edges have valid ddg_uncertainty (not NaN)
        edges_no_uncertainty = [
            (result["ligand_a"], result["ligand_b"])
            for result in system_results
            if np.isnan(result["ddg_uncertainty"].magnitude)
        ]
        if edges_no_uncertainty:
            raise ValueError(
                f"Not all edges have ddg_uncertainty for {system_group} {system_name}: {edges_no_uncertainty}"
            )

        femap = FEMap()
        for result in system_results:
            ligand_a = result["ligand_a"]
            ligand_b = result["ligand_b"]
            # record the ligands added to the femap
            unique_ligands.update([ligand_a, ligand_b])
            femap.add_relative_calculation(
                labelA=ligand_a,
                labelB=ligand_b,
                value=result["ddg"],
                uncertainty=result["ddg_uncertainty"],
                source=source,
            )

        # add experimental data for each of the ligands in the results
        experimental_file = benchmark_data.reference_data["experimental_binding_data"]
        with open(experimental_file) as f:
            experimental_data = json.load(f, cls=JSON_HANDLER.decoder)

        for ligand in unique_ligands:
            exp_data = experimental_data.get(ligand, None)
            if exp_data is not None:
                femap.add_experimental_measurement(
                    label=ligand,
                    value=exp_data["dg"],
                    uncertainty=exp_data.get(
                        "uncertainty", 0 * unit.kilocalories_per_mole
                    ),
                    source="experimental",
                )

        femaps_by_system_key[system_key] = femap
    return femaps_by_system_key


def build_femap_from_absolute_results(
    results: list[dict],
    calculation_type: str = "asfe",
    source: str = "",
) -> dict[tuple[str, str], FEMap]:
    """
    Build FEMaps for each of the unique combinations of system_group and system_name in the absolute results
    and add experimental data where available.

    Parameters
    ----------
    results: list[dict]
        A list of absolute free energy estimates. ASFE and RBFE reference-value
        records are both supported. Each record should include at least:
         - ligand: str
         - system_group: str
         - system_name: str
         - dg: Quantity
         - dg_uncertainty: Quantity

    calculation_type: str, default='asfe'
        Type of calculation: 'asfe' for absolute solvation, 'rbfe' for relative binding
    source: str, default=''
        Source of computed data, used in ``cinnabar`` for results comparisons between different BenchmarkResults.

    Returns
    -------
    dict[tuple[str, str], FEMap]
        A dictionary mapping each unique combination of system_group and system_name to an FEMap with calculated
        and experimental data (where available).
    """
    # Detect format from data if not specified.
    # The presence of `solvent` distinguishes ASFE records from RBFE absolute references.
    if results and calculation_type == "asfe":
        first_result = results[0]
        if "solvent" not in first_result:
            calculation_type = "rbfe"

    # get the unique combinations of system_group and system_name
    results_by_system_key = defaultdict(list)
    for result in results:
        key = (result["system_group"], result["system_name"])
        results_by_system_key[key].append(result)

    femaps_by_system_key = {}
    for system_key, system_results in results_by_system_key.items():
        system_group, system_name = system_key
        benchmark_data = get_benchmark_data_system(system_group, system_name)
        if benchmark_data.reference_data is None:
            raise ValueError(
                f"No reference data available for benchmark system {system_group}/{system_name}"
            )

        # Check if all molecules have valid uncertainty (not NaN).
        missing_required_fields = [
            _asfe_result_key(result)
            for result in system_results
            if np.isnan(result["dg_uncertainty"].magnitude)
        ]
        if missing_required_fields:
            raise ValueError(
                "Absolute results must include both 'dg' and 'dg_uncertainty' "
                f"for {system_group} {system_name}: {missing_required_fields}"
            )

        missing_uncertainty = [
            _asfe_result_key(result)
            for result in system_results
            if np.isnan(result["dg_uncertainty"].magnitude)
        ]
        if missing_uncertainty:
            raise ValueError(
                "Not all molecules have uncertainty "
                f"for {system_group} {system_name}: {missing_uncertainty}"
            )

        femap = FEMap()
        for result in system_results:
            label = _asfe_result_key(result)
            femap.add_absolute_calculation(
                label=label,
                value=result["dg"],
                uncertainty=result["dg_uncertainty"],
                source=source,
            )

        # Add experimental data if available.
        # ASFE records map to solvation data; RBFE absolute references map to
        # binding data.
        experimental_key = (
            "experimental_solvation_free_energy_data"
            if calculation_type == "asfe"
            else "experimental_binding_data"
        )
        if experimental_key not in benchmark_data.reference_data:
            femaps_by_system_key[system_key] = femap
            continue

        experimental_file = benchmark_data.reference_data[experimental_key]
        with open(experimental_file) as f:
            experimental_data = json.load(f, cls=JSON_HANDLER.decoder)
        n_experimental_points = 0
        for result in system_results:
            label = _asfe_result_key(result)
            exp_data = experimental_data.get(label, None)
            if exp_data is not None:
                femap.add_experimental_measurement(
                    label=label,
                    value=exp_data["dg"],
                    uncertainty=exp_data.get(
                        "uncertainty", 0 * unit.kilocalories_per_mole
                    ),
                    source="experimental",
                )
                n_experimental_points += 1

        if n_experimental_points == 0 and calculation_type == "asfe":
            raise ValueError("No experimental data points were found")

        femaps_by_system_key[system_key] = femap

    return femaps_by_system_key
