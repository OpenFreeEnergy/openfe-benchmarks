import json

from rdkit import Chem

from openff.units import unit
import openfe
from openfe import (
    SmallMoleculeComponent, SolventComponent, ProteinComponent,
)
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol

from openfe_benchmarks.data import get_benchmark_data_system
from openfe_benchmarks.scripts import utils as ofebu

SOLVENT = SolventComponent(positive_ion='Na', negative_ion='Cl', neutralize=True)
BENCHMARK_SET = "mcs_docking_set"
BENCHMARK_SYS = "hne"
PARTIAL_CHARGE = "nagl_openff-gnn-am1bcc-1.0.0.pt"
FORCEFIELD = 'openff-2.3.0'
FILENAME = f"network_{BENCHMARK_SET}_{BENCHMARK_SYS}_nacl.json"

def process_components(benchmark_sys):

    if benchmark_sys.network is None:
        raise ValueError("Valid protein network.json is required.")
    network0 = openfe.LigandNetwork.from_json(file=str(benchmark_sys.network))
    ligand_dict = ofebu.process_sdf(benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=True)
    if benchmark_sys.protein is None:
        raise ValueError("Valid protein pdb is required.")
    protein = ProteinComponent.from_pdb_file(str(benchmark_sys.protein))

    cofactors = None
    if benchmark_sys.cofactors is not None:
        cofactors = ofebu.process_sdf(benchmark_sys.ligands[PARTIAL_CHARGE], return_dict=False)

    return network0, ligand_dict, protein, cofactors

def compile_network_transformations(network, solvent, ligands_by_name, protein, cofactors):
    transformations = []
    for edge in network.edges:
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
    
            transformation_name = "test" + "_" + system_a.name + "_" + system_b.name
            if "vacuum" in transformation_name: # usually detected with transformation_name, likely doesn't apply here
                protocol_settings.nonbonded_method = "nocutoff"

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
    benchmark_sys = get_benchmark_data_system(BENCHMARK_SET, BENCHMARK_SYS)
    network0, ligand_dict, protein, cofactors = process_components(benchmark_sys)
    
    transformations = compile_network_transformations(
        network0, SOLVENT, ligand_dict, protein, cofactors
    )
    network = openfe.AlchemicalNetwork(edges=transformations)
    network.to_json(file=FILENAME)

if __name__ == "__main__":
    main()