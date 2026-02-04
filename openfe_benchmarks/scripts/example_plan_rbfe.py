"""
This module provides an example script for planning relative binding free energy (RBFE) calculations.
It demonstrates how to process benchmark systems, compile alchemical transformations, and save the resulting network to a JSON file.
"""

import openfe
from openfe import SolventComponent, ProteinComponent
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol

from openfe_benchmarks.data import get_benchmark_data_system
from openfe_benchmarks.scripts import utils as ofebu

SOLVENT = SolventComponent(positive_ion='Na', negative_ion='Cl', neutralize=True)
BENCHMARK_SET = "mcs_docking_set"
BENCHMARK_SYS = "hne"
PARTIAL_CHARGE = "nagl_openff-gnn-am1bcc-1.0.0.pt"
FORCEFIELD = 'openff-2.3.0'
FILENAME_ALCHEMICALNETWORK = f"alchemical_network_{BENCHMARK_SET}_{BENCHMARK_SYS}_nacl.json"

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

    if benchmark_sys.network is None:
        raise ValueError("Valid protein network.json is required.")
    lig_network = openfe.LigandNetwork.from_json(file=str(benchmark_sys.network))
    ligand_dict = ofebu.process_sdf(benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=True)
    if benchmark_sys.protein is None:
        raise ValueError("Valid protein pdb is required.")
    protein = ProteinComponent.from_pdb_file(str(benchmark_sys.protein))

    cofactors = None
    if benchmark_sys.cofactors is not None:
        cofactors = ofebu.process_sdf(benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=False)

    return lig_network, ligand_dict, protein, cofactors

def compile_network_transformations(ligand_network, solvent, ligands_by_name, protein, cofactors):
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
            system_a_dict = {
                "ligand": new_edge.componentA, "solvent": solvent
            }
            system_b_dict = {
                "ligand": new_edge.componentB, "solvent": solvent
            }
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
            transformation_protocol = RelativeHybridTopologyProtocol
            protocol_settings = RelativeHybridTopologyProtocol.default_settings()

            if isinstance(transformation_protocol, RelativeHybridTopologyProtocol):
                # adaptive transformation settings are only supported for RelativeHybridTopologyProtocol currently
                protocol_settings = transformation_protocol._adaptive_settings(
                    stateA=system_a,
                    stateB=system_b,
                    mapping=new_edge,
                    initial_settings=protocol_settings,
                )

            transformation_protocol = RelativeHybridTopologyProtocol(settings=protocol_settings)

            # Create transformation
            transformation = openfe.Transformation(
                stateA=system_a,
                stateB=system_b,
                mapping=new_edge,
                protocol=transformation_protocol,
                name=name
            )
            transformations.append(transformation)

    return transformations


def main():
    """
    Main function to process the benchmark system and save the alchemical network.

    This function retrieves the benchmark system, processes its components, compiles the transformations,
    and saves the resulting alchemical network to a JSON file.
    """
    benchmark_sys = get_benchmark_data_system(BENCHMARK_SET, BENCHMARK_SYS)
    lig_network, ligand_dict, protein, cofactors = process_components(benchmark_sys)
    
    transformations = compile_network_transformations(
        lig_network, SOLVENT, ligand_dict, protein, cofactors
    )
    alchem_network = openfe.AlchemicalNetwork(edges=transformations)
    alchem_network.to_json(file=FILENAME_ALCHEMICALNETWORK)

if __name__ == "__main__":
    main()