"""
This module provides an example script for planning relative binding free energy (RBFE) calculations.
It demonstrates how to process benchmark systems, compile alchemical transformations, and save the
resulting network to a JSON file.

The alchemical_network.json file is a consolidation of each of the transformation files. It is the
latter files that can be run on an HPC simply with:
`openfe quickrun path/to/transformation.json -o results.json -d working-directory`
See the [openfe tutorials](https://docs.openfree.energy/en/latest/tutorials/rbfe_cli_tutorial.html)
for more details.
"""

import os
import logging
import json

import openfe
from openfe import SolventComponent, ProteinComponent
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol

from openfe_benchmarks.data import get_benchmark_data_system
from openfe_benchmarks.scripts import utils as ofebu

logger = logging.getLogger(__name__)


def _configure_example_logging(level=logging.INFO):
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s")
    )
    # Attach to this module's logger so output appears when running the example
    logger.addHandler(handler)
    logger.setLevel(level)
    # Optionally enable package-wide logs:
    logging.getLogger("openfe_benchmarks").setLevel(level)


SOLVENT = SolventComponent(positive_ion="Na", negative_ion="Cl", neutralize=True)
BENCHMARK_SET = "mcs_docking_set"
BENCHMARK_SYS = "hne"
PARTIAL_CHARGE = "nagl_openff-gnn-am1bcc-1.0.0.pt"  # for the ligand and cofactors
FORCEFIELD = "openff-2.3.0"  # available [openmmforcefields SystemGenerator](https://github.com/openmm/openmmforcefields?tab=readme-ov-file#automating-force-field-management-with-systemgenerator)
LIG_NETWORK_FILE = "industry_benchmarks_network"
FILENAME_ALCHEMICALNETWORK = (
    f"alchemical_network_{BENCHMARK_SET}_{BENCHMARK_SYS}_nacl.json"
)
OUTPUT_DIR = "outputs"


def process_components(benchmark_sys):
    """
    Process the components of a benchmark system.

    Parameters
    ----------
    benchmark_sys : BenchmarkData
        The benchmark system to process.

    Returns
    -------
    tuple
        A tuple containing the ligand network, ligand dictionary, protein component, and cofactors.
    """

    if not benchmark_sys.ligand_networks:
        raise ValueError("Valid ligand network JSON file is required.")
    lig_network = openfe.LigandNetwork.from_json(
        file=str(benchmark_sys.ligand_networks[LIG_NETWORK_FILE])
    )
    ligand_dict = ofebu.process_sdf(
        benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=True
    )
    if benchmark_sys.protein is None:
        raise ValueError("Valid protein pdb is required.")
    protein = ProteinComponent.from_pdb_file(str(benchmark_sys.protein))

    cofactors = None
    if benchmark_sys.cofactors is not None:
        cofactors = ofebu.process_sdf(
            benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=False
        )

    # Apply custom partial charges to ligands and cofactors here

    return lig_network, ligand_dict, protein, cofactors


def compile_network_transformations(
    ligand_network, solvent, ligands_by_name, protein, cofactors
):
    """
    Compile alchemical transformations for a given network.

    Parameters
    ----------
    ligand_network : LigandNetwork
        The alchemical network containing edges to transform.
    solvent : SolventComponent
        The solvent component for the transformations.
    ligands_by_name : dict
        Dictionary mapping ligand names to SmallMoleculeComponent objects.
    protein : ProteinComponent
        The protein component for the transformations.
    cofactors : list or None
        List of cofactor components, if any.

    Returns
    -------
    list
        A list of alchemical transformations.
    """
    transformations = []
    for edge in ligand_network.edges:
        new_edge = openfe.LigandAtomMapping(
            componentA=ligands_by_name[edge.componentA.name],
            componentB=ligands_by_name[edge.componentB.name],
            componentA_to_componentB=edge.componentA_to_componentB,
            annotations=edge.annotations,
        )

        # create the transformations for the bound and solvent legs
        for leg in ["solvent", "complex"]:
            system_a_dict = {"ligand": new_edge.componentA, "solvent": solvent}
            system_b_dict = {"ligand": new_edge.componentB, "solvent": solvent}
            if leg == "complex":
                system_a_dict["protein"] = protein
                system_b_dict["protein"] = protein

                if cofactors is not None:
                    for i, cofactor in enumerate(cofactors):
                        cofactor_name = f"cofactor_{i}"
                        system_a_dict[cofactor_name] = cofactor
                        system_b_dict[cofactor_name] = cofactor

            system_a = openfe.ChemicalSystem(system_a_dict)
            system_b = openfe.ChemicalSystem(system_b_dict)

            name = f"{leg}_{new_edge.componentA.name}_{new_edge.componentB.name}"

            # Create protocol with adaptive settings
            # adaptive transformation settings are only supported for RelativeHybridTopologyProtocol currently
            protocol_settings = RelativeHybridTopologyProtocol.default_settings()
            protocol_settings.forcefield_settings.small_molecule_forcefield = FORCEFIELD
            transformation_protocol = RelativeHybridTopologyProtocol(
                settings=RelativeHybridTopologyProtocol._adaptive_settings(
                    stateA=system_a,
                    stateB=system_b,
                    mapping=new_edge,
                    initial_settings=protocol_settings,
                ),
            )

            # Create transformation
            transformation = openfe.Transformation(
                stateA=system_a,
                stateB=system_b,
                mapping=new_edge,
                protocol=transformation_protocol,
                name=name,
            )
            transformations.append(transformation)

    return transformations


def main():
    """
    Main function to process the benchmark system and save the alchemical network.

    This function retrieves the benchmark system, processes its components, compiles the transformations,
    and saves the resulting alchemical network to a JSON file.
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    benchmark_sys = get_benchmark_data_system(BENCHMARK_SET, BENCHMARK_SYS)
    lig_network, ligand_dict, protein, cofactors = process_components(benchmark_sys)

    transformations = compile_network_transformations(
        lig_network, SOLVENT, ligand_dict, protein, cofactors
    )

    # Can be used as input for Alchemiscale
    alchem_network = openfe.AlchemicalNetwork(edges=transformations)
    alchem_network.to_json(file=os.path.join(OUTPUT_DIR, FILENAME_ALCHEMICALNETWORK))


#    # Can be run HPC3
#    for transformation in alchem_network.edges:
#        transformation.to_json(os.path.join(OUTPUT_DIR, f"{transformation.name}.json"))


# __________________________ Validation Function ________________________________


def validate_rbfe_network(network_file):
    """Validate RBFE network against BenchmarkData expectations.

    Checks:
    - Exact number of edges (2 per ligand pair: complex + solvent)
    - Transformations match expected ligand network
    - Each transformation has protein, solvent, and cofactors (if present)
    - Ligand charges match expected values from BenchmarkData
    """

    errors = []

    # Load generated alchemical network
    try:
        with open(network_file) as f:
            network_json = json.load(f)
    except Exception as e:
        return [f"Failed to load network: {e}"]

    if not isinstance(network_json, list) or len(network_json) == 0:
        return ["Invalid network JSON structure"]

    logger.info(
        "Loaded alchemical network JSON with %d top-level items", len(network_json)
    )

    # Last element is the AlchemicalNetwork
    network_data = network_json[-1][1]
    edges = network_data.get("edges", [])
    logger.info("AlchemicalNetwork contains %d edges", len(edges))

    # Get benchmark data for validation
    benchmark_sys = get_benchmark_data_system(BENCHMARK_SET, BENCHMARK_SYS)
    expected_lig_network = openfe.LigandNetwork.from_json(
        file=str(benchmark_sys.ligand_networks[LIG_NETWORK_FILE])
    )

    # Check exact number of edges (2 per ligand edge: complex + solvent)
    expected_edge_count = len(expected_lig_network.edges) * 2
    actual_edge_count = len(edges)
    if actual_edge_count != expected_edge_count:
        errors.append(
            f"Expected exactly {expected_edge_count} edges "
            f"({len(expected_lig_network.edges)} ligand pairs × 2 legs), "
            f"got {actual_edge_count}"
        )
    else:
        logger.info(
            "Edge count correct: %d edges (expected %d)",
            actual_edge_count,
            expected_edge_count,
        )

    # Build ChemicalSystem lookup from network JSON to check components
    chem_systems = {}
    for item in network_json:
        if len(item) == 2 and "ChemicalSystem" in item[0]:
            chem_systems[item[0]] = item[1]

    # Validate each transformation
    transformation_names = set()
    complex_legs = []
    solvent_legs = []

    for edge_ref in edges:
        edge_key = edge_ref.get(":gufe-key:")
        # Find the transformation details
        transformation = None
        for item in network_json:
            if len(item) == 2 and item[0] == edge_key:
                transformation = item[1]
                break

        if not transformation:
            errors.append(f"Could not find transformation data for {edge_key}")
            continue

        name = transformation.get("name", "")
        transformation_names.add(name)

        # Check leg type (complex vs solvent)
        if name.startswith("complex_"):
            complex_legs.append(name)
        elif name.startswith("solvent_"):
            solvent_legs.append(name)

        # Get the ChemicalSystem components for stateA
        stateA_key = transformation.get("stateA", {}).get(":gufe-key:")
        if stateA_key and stateA_key in chem_systems:
            components = chem_systems[stateA_key].get("components", {})

            # Check required components based on leg type
            if name.startswith("complex_"):
                # Complex leg must have protein, solvent, ligand
                required = ["protein", "solvent", "ligand"]
                for comp in required:
                    if comp not in components:
                        errors.append(
                            f"Transformation '{name}' missing {comp} component"
                        )

                # Check for cofactors if benchmark system has them
                if benchmark_sys.cofactors is not None:
                    has_cofactor = any("cofactor" in k for k in components.keys())
                    if not has_cofactor:
                        errors.append(
                            f"Transformation '{name}' missing cofactor component"
                        )
                    else:
                        logger.info("Transformation '%s' includes cofactors", name)

                # Log success if none of the required components were missing
                missing = [c for c in required if c not in components]
                if not missing:
                    logger.info(
                        "Transformation '%s' complex leg has required components", name
                    )

            elif name.startswith("solvent_"):
                # Solvent leg must have solvent and ligand
                required = ["solvent", "ligand"]
                for comp in required:
                    if comp not in components:
                        errors.append(
                            f"Transformation '{name}' missing {comp} component"
                        )
                else:
                    logger.info(
                        "Transformation '%s' solvent leg has required components", name
                    )

    # Validate that we have both complex and solvent legs for each ligand pair
    expected_ligand_pairs = len(expected_lig_network.edges)
    if len(complex_legs) != expected_ligand_pairs:
        errors.append(
            f"Expected {expected_ligand_pairs} complex legs, got {len(complex_legs)}"
        )
    if len(solvent_legs) != expected_ligand_pairs:
        errors.append(
            f"Expected {expected_ligand_pairs} solvent legs, got {len(solvent_legs)}"
        )

    # Validate that ligands have partial charges from BenchmarkData
    # Load expected ligands to get count
    expected_ligands = ofebu.process_sdf(
        benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=True
    )

    # Map ligands in the network by their declared name (molprops['ofe-name'])
    lig_map = {}
    for item in network_json:
        if len(item) == 2 and "SmallMoleculeComponent" in item[0]:
            ligand_data = item[1]
            molprops = ligand_data.get("molprops", {}) or {}
            ofe_name = (
                str(molprops.get("ofe-name"))
                if molprops.get("ofe-name") is not None
                else None
            )

            # partial charges are stored under molprops e.g. 'atom.dprop.PartialCharge'
            has_charges = False
            for k in molprops.keys():
                if (
                    "partialcharge" in k.lower()
                    or "partial_charge" in k.lower()
                    or "partial_charge" in k
                ):
                    val = molprops.get(k)
                    if val is not None and (
                        not (isinstance(val, str) and val.strip() == "")
                    ):
                        has_charges = True
                        break

            if ofe_name:
                lig_map[ofe_name] = lig_map.get(ofe_name, False) or has_charges

    expected_keys = set(str(k) for k in expected_ligands.keys())
    found_keys = set(lig_map.keys())

    missing_expected = expected_keys - found_keys
    if missing_expected:
        errors.append(f"Ligands not found in network: {sorted(list(missing_expected))}")

    ligands_without_charges = [k for k, v in lig_map.items() if not v]
    if ligands_without_charges:
        errors.append(
            f"Ligands missing partial charges: {sorted(ligands_without_charges)}"
        )
    return errors


if __name__ == "__main__":
    _configure_example_logging(level=logging.INFO)
    main()
    validate_rbfe_network(os.path.join(OUTPUT_DIR, FILENAME_ALCHEMICALNETWORK))
