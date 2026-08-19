from collections import defaultdict
import json

import numpy as np
from gufe.tokenization import JSON_HANDLER
from openff.units import unit
from cinnabar import FEMap

from openfe_benchmarks.data._benchmark_systems import get_benchmark_data_system


def build_femap_from_relative_results(
    results: list[dict],
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
    unique_ligands = set()
    for system_key, system_results in results_by_system_key.items():
        system_group, system_name = system_key
        benchmark_data = get_benchmark_data_system(system_group, system_name)

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
            )

        # add experimental data for each of the ligands in the results
        experimental_file = benchmark_data.reference_data["experimental_binding_data"]
        experimental_data = json.load(open(experimental_file), cls=JSON_HANDLER.decoder)

        for ligand in unique_ligands:
            exp_data = experimental_data.get(ligand, None)
            if exp_data is not None:
                femap.add_experimental_measurement(
                    label=ligand,
                    value=exp_data["dg"],
                    uncertainty=exp_data.get(
                        "uncertainty", 0 * unit.kilocalories_per_mole
                    ),
                )

        femaps_by_system_key[system_key] = femap
    return femaps_by_system_key


def build_femap_from_absolute_results(
    results: list[dict],
    calculation_type: str = "asfe",
) -> dict[tuple[str, str], FEMap]:
    """
    Build FEMaps for each of the unique combinations of system_group and system_name in the absolute results
    and add experimental data where available.

    Parameters
    ----------
    results: list[dict]
        A list of absolute free energy estimates. Format depends on calculation_type:

        For ASFE (absolute solvation):
         - solute: str
         - solvent: str
         - system_group: str
         - system_name: str
         - dg or estimate: Quantity
         - dg_uncertainty or estimate_error: Quantity

        For RBFE (absolute binding reference values):
         - ligand: str
         - system_group: str
         - system_name: str
         - dg: Quantity
         - dg_uncertainty: Quantity

    calculation_type: str, default='asfe'
        Type of calculation: 'asfe' for absolute solvation, 'rbfe' for relative binding

    Returns
    -------
    dict[tuple[str, str], FEMap]
        A dictionary mapping each unique combination of system_group and system_name to an FEMap with calculated
        and experimental data (where available).
    """
    # Detect format from data if not specified
    if results and calculation_type == "asfe":
        # Check if this is actually RBFE format
        first_result = results[0]
        if "ligand" in first_result and "solute" not in first_result:
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

        # Determine molecule key and value/error keys based on calculation type
        if calculation_type == "asfe":
            molecule_key = "solute"
            has_solvent = True
        else:  # rbfe
            molecule_key = "ligand"
            has_solvent = False

        # Check if all molecules have valid uncertainty (not NaN)
        molecules_no_uncertainty = [
            result[molecule_key]
            for result in system_results
            if np.isnan(
                result["dg_uncertainty"].magnitude
                if "dg" in result
                else result["estimate_error"].magnitude
            )
        ]
        if molecules_no_uncertainty:
            raise ValueError(
                f"Not all {molecule_key}s have uncertainty for {system_group} {system_name}: {molecules_no_uncertainty}"
            )

        femap = FEMap()
        for result in system_results:
            value_key = "dg" if "dg" in result else "estimate"
            err_key = "dg_uncertainty" if value_key == "dg" else "estimate_error"

            # Build label based on format
            if has_solvent:
                label = f"{result[molecule_key]},{result['solvent']}"
            else:
                label = result[molecule_key]

            femap.add_absolute_calculation(
                label=label,
                value=result[value_key],
                uncertainty=result[err_key],
                source="Computational",
            )

        # Add experimental data if available
        # For ASFE: experimental solvation free energy data
        # For RBFE: experimental binding free energy data (if available)
        experimental_key = (
            "experimental_solvation_free_energy_data"
            if calculation_type == "asfe"
            else "experimental_binding_free_energy_data"
        )

        if experimental_key in benchmark_data.reference_data:
            experimental_file = benchmark_data.reference_data[experimental_key]
            with open(experimental_file, "r") as f:
                experimental_data = json.load(f, cls=JSON_HANDLER.decoder)
            n_experimental_points = 0

            for result in system_results:
                # Build label to match experimental data format
                if has_solvent:
                    label = f"{result[molecule_key]},{result['solvent']}"
                else:
                    label = result[molecule_key]

                exp_data = experimental_data.get(label, None)
                if exp_data is not None:
                    femap.add_experimental_measurement(
                        label=label,
                        value=exp_data["dg"],
                        uncertainty=exp_data.get(
                            "uncertainty", 0 * unit.kilocalories_per_mole
                        ),
                    )
                    n_experimental_points += 1

            if n_experimental_points == 0 and calculation_type == "asfe":
                # Only raise error for ASFE where experimental data is expected
                raise ValueError(
                    "No experimental data points were found for ASFE submission"
                )

        femaps_by_system_key[system_key] = femap

    return femaps_by_system_key
