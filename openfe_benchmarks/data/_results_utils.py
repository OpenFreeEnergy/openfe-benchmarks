from cinnabar import FEMap
from openfe_benchmarks.data._benchmark_systems import get_benchmark_data_system
from collections import defaultdict
import json
from gufe.tokenization import JSON_HANDLER
from openff.units import unit

def build_femap_from_relative_results(results: list[dict]) -> dict[tuple[str, str], FEMap]:
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
         - DDG: Quantity
         - DDG_uncertainty: Quantity

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
        femap = FEMap()
        for result in system_results:
            ligand_a = result["ligand_a"]
            ligand_b = result["ligand_b"]
            # record the ligands added to the femap
            unique_ligands.update([ligand_a, ligand_b])
            ddg = result["DDG"]
            ddg_uncertainty = result["DDG_uncertainty"]
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
                    uncertainty=exp_data.get("uncertainty", 0 * unit.kilocalorie_per_mole),
                )

        femaps_by_system_key[system_key] = femap
    return femaps_by_system_key