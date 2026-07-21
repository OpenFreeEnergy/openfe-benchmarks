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
        Optionally may include:
         - mbar_std: Quantity (used as fallback if ddg_uncertainty is NaN)

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

        # First, determine globally which uncertainty type to use for this system
        # Check if all edges have valid ddg_uncertainty (not NaN)
        all_have_repeat_uncertainty = all(
            not np.isnan(result["ddg_uncertainty"].magnitude)
            for result in system_results
        )

        # Decide which uncertainty to use for the entire network
        use_mbar_std = not all_have_repeat_uncertainty

        femap = FEMap()
        for result in system_results:
            ligand_a = result["ligand_a"]
            ligand_b = result["ligand_b"]
            # record the ligands added to the femap
            unique_ligands.update([ligand_a, ligand_b])
            ddg = result["ddg"]

            # Use consistent uncertainty type for all edges in this system
            if use_mbar_std and "mbar_std" in result:
                ddg_uncertainty = result["mbar_std"]
            else:
                ddg_uncertainty = result["ddg_uncertainty"]

            femap.add_relative_calculation(
                labelA=ligand_a,
                labelB=ligand_b,
                value=ddg,
                uncertainty=ddg_uncertainty,
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
                        "uncertainty", 0 * unit.kilocalorie_per_mole
                    ),
                )

        femaps_by_system_key[system_key] = femap
    return femaps_by_system_key
